/* prototype.c --- 
 * 
 * Filename: prototype.c
 * Description: 
 * Author: O.TISSOT 
 * Created: ven. mai  9 14:48:36 2014 (+0200)
 * Version: 
 * Last-Updated: sam. mai 10 20:37:30 2014 (+0200)
 *           By: Olivier
 *     Update #: 159
 * Compatibility: 
 * ScaLAPACK and PBLAS (and their dependancies)
 */

/* Commentary: 
 * A simple prototype of the algorithm studied in the project.
 * It just do the low-rank factorization part.
 * I used examples  written by Kelly McQuighan for the
 * Kobe-Brown summer school on high performance simulations.
 * They are available on her personnal page :
 * http://www.dam.brown.edu/people/kellym
 */


/* Code: */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "scalapack_headers.h"
#include "blacs_headers.h"
#include "pblas_headers.h"
#include "scalapack_tools_headers.h"
#include "pspblasinfo.h"
#include "psmatgen.h"
#include "stringTools.h"
#include "pblasIOtools.h"

static int imax( int a, int b ){
  if (a>b) return(a); else return(b);
}

static int imin( int a, int b ){
  if (a<b) return(a); else return(b);
}

//------MAIN-----------------
int main(int argc, char *argv[])
{
  int i_am, k, m, n, num_proc_col, num_procs, num_proc_row, block_size, sys_param[6], status;
  int descA[9], descB[9], descC[9], info, maxmp, maxkp, maxmkn, minmpn;
  int kp, kq, mp, my_col, my_row, nq, contxt, izero = 0, ione = 1, len_mat_name, lwork, ltau;
  int iaseed = 100, ibseed = 200, outfile_name_len;
  FILE* outfile=NULL;
  char folder[] = "data";
  const char * infile_name = "data/params.dat";
  char outfile_name[100], usr_info[100];
  double *A=NULL, *B=NULL, *C=NULL, *A0=NULL, *C0=NULL, *tau=NULL, *work=NULL;
  double done = 1.0e0, dmone = -1.0e0, resid = 0.0e0, eps = 1e-12, wall_time = 0.0e0;

  //****************************************************
  //                INITIALIZE BLACS
  //****************************************************

  Cblacs_pinfo( &i_am, &num_procs );
  // read data from file BLAS.dat. Exit program if file not formatted properly
  pspblasinfo( &m, &n, &k, &block_size, &num_proc_row, &num_proc_col, sys_param, &i_am, &num_procs, folder, infile_name, outfile_name, usr_info, &outfile_name_len, &status );
  if ( status < 0 ) { 
    return EXIT_FAILURE;
  }

  if (i_am == 0){
    outfile = fopen( "data/results.out", "w" );

    // output data
    fprintf ( outfile, "This is a prototype of the randomized low-rank factorization described in :\n" );
    fprintf ( outfile, " ''Randomized algorithms for for low-rank matrix factorizations : sharp performance bounds'', E. Candes & R. Witten, 2013\n" );
    fprintf ( outfile, "It is available on E. Candes's web page.\n" );
    fprintf ( outfile, "It is written in parallel thanks to ScaLAPACK by Olivier Tissot. It is a part of a project for Laura Grigori's class about parallel computing given in 2nd semester of 2013-2014 of Master 2 AN&EDP of UPMC. A little report explains the methodology of the work. If you want it, please send me an e-mail : oli.tissot@gmail.com.\n" );
    fprintf ( outfile, "The code is inspired from examples written by Kelly McQuighan for a summer school on high performance simulations. They are availables on her webpage.\n" );
    fprintf ( outfile, "****************************************************************************************************\n");
    fprintf ( outfile, "An explanation of the input/output parameters is as follows :\n" );
    fprintf ( outfile, "m\t\t: The number of rows in the matrices A and C.\t\t\tm\t\t= %d\n", m );
    fprintf ( outfile, "n\t\t: The number of columns in the matrices B and C.\t\tn\t\t= %d\n", n );
    fprintf ( outfile, "k\t\t: The number of rows of B and the number of columns of A.\tk\t\t= %d\n", k );
    fprintf ( outfile, "block_size\t: The size of the square blocks the matrices A, B, and C are\n\t\t  split into.\t\t\t\t\t\t\tblock_size\t= %d\n", block_size );
    fprintf ( outfile, "num_proc_row\t: The number of process rows.\t\t\t\t\tnum_proc_row\t= %d\n", num_proc_row );
    fprintf ( outfile, "num_proc_col\t: The number of process columns.\t\t\t\tnum_proc_col\t= %d\n", num_proc_col );
    fprintf ( outfile, "****************************************************************************************************\n");
    fprintf ( outfile, "\nProgram progress:\n" );
  }

  if (num_proc_row*num_proc_col > num_procs){
    if (i_am==0) {
      fprintf ( stdout, "\n****************************************************************************************************\n\n");
      fprintf ( stdout, "ERROR : we do not have enough processes available to make a nprow x npcol process grid\n\n" );
      fprintf ( stdout, "****************************************************************************************************\n");
    }
    return EXIT_FAILURE;
  }

  Cblacs_get( -1, 0, &contxt );
  Cblacs_gridinit( &contxt, "Row", num_proc_row, num_proc_col );
  Cblacs_gridinfo( contxt, &num_proc_row, &num_proc_col, &my_row, &my_col );

  // only continue if I'm on the grid
  if ( ( my_row < num_proc_row ) & ( my_col < num_proc_col ) ) {

   // determine local matrix dimensions on each processor
   // usage: A_loc = mp x kq, B_loc = kp x nq, C_loc = mp x nq
    mp = numroc_( &m, &block_size, &my_row, &izero, &num_proc_row );
    kp = numroc_( &k, &block_size, &my_row, &izero, &num_proc_row );
    kq = numroc_( &k, &block_size, &my_col, &izero, &num_proc_col );
    nq = numroc_( &n, &block_size, &my_col, &izero, &num_proc_col );
    maxmp = imax(1,mp);
    maxkp = imax(1,kp);
    maxmkn = imax(imax(m,k),n);
    minmpn = imin(n, m);
    lwork = maxmkn*maxmkn;
    ltau = minmpn;

    // initialize descriptors for the matrices A, B, C and R
    descinit_( descA, &m, &k, &block_size, &block_size, &izero, &izero, &contxt, &maxmp, &info );
    descinit_( descB, &k, &n, &block_size, &block_size, &izero, &izero, &contxt, &maxkp, &info );
    descinit_( descC, &m, &n, &block_size, &block_size, &izero, &izero, &contxt, &maxmp, &info );

    if (i_am==0) fprintf ( outfile, "BLACS functions numroc and descinit completed successfully.\n");

    //****************************************************
    //                INITIALIZE MATRICES
    //****************************************************
    // allocate memory for the local part of matrices A, B, C and R
    A   = (double*) malloc(maxmp*kq*sizeof(double));
    A0  = (double*) malloc(maxmp*kq*sizeof(double));
    B   = (double*) malloc(maxkp*nq*sizeof(double));
    C   = (double*) malloc(maxmp*nq*sizeof(double));
    C0  = (double*) malloc(maxmp*nq*sizeof(double));
    tau = (double*) malloc(ltau*sizeof(double));
    work = (double*) malloc(lwork*sizeof(double));
    
    // generate random matrices A, B, and zero matrix C
    psmatgen("N", "N", A, 1, 1, descA, m, k, my_row, my_col, num_proc_row, num_proc_col, iaseed/*+time(NULL)*/, &status );
    if (status < 0) {
      Cblacs_exit( 0 );
      return EXIT_FAILURE;
    }
    psmatgen("N", "N", B, 1, 1, descB, k, n, my_row, my_col, num_proc_row, num_proc_col, ibseed/*+time(NULL)*/, &status );
    if (status < 0) {
      Cblacs_exit( 0 );
      return EXIT_FAILURE;
    }
    psmatgen("Z", "N", C, 1, 1, descC, m, n, my_row, my_col, num_proc_row, num_proc_col, 0, &status );
    if (status < 0) {
      Cblacs_exit( 0 );
      return EXIT_FAILURE;
    }
    pdlacpy_ ( "All", &m, &k, A, &ione, &ione, descA, A0, &ione, &ione, descA );
    if (i_am==0) fprintf ( outfile, "Matrices A, B, C and R successfully created.\n");  
  
    len_mat_name = 1;

    //****************************************************
    //        FIRST STEP : PRODUCT MATRIX-MATRIX
    //****************************************************

    wall_time -= Cdwalltime00( );

    if (i_am==0) {
      fprintf ( outfile, "\n****************************************************************************************************" );
      fprintf ( outfile, "\nFirst step : product matrix-matrix\n" );
      fprintf ( outfile, "****************************************************************************************************\n\n" );
      fprintf ( outfile, "Matrix A:\n" );
    }
    pdlaprnt2( &m, &k, A, &ione, &ione, descA, &izero, &izero, "A", outfile, work, len_mat_name );
    if (i_am==0) fprintf ( outfile, "\nMatrix G:\n" );
    pdlaprnt2( &k, &n, B, &ione, &ione, descB, &izero, &izero, "G", outfile, work, len_mat_name );

    pdgemm_ ( "N", "N", &m, &n, &k, &done, A0, &ione, &ione, descA, B, &ione, &ione, descB, &done, C, &ione, &ione, descC );
    if (i_am==0) fprintf ( outfile, "\nH := A*G\n" );
    pdlaprnt2( &m, &n, C, &ione, &ione, descC, &izero, &izero, "H", outfile, work, len_mat_name );

    //****************************************************
    //          SECOND STEP : QR DECOMPOSITION
    //****************************************************

    if (i_am==0) {
      fprintf ( outfile, "\n****************************************************************************************************" );
      fprintf ( outfile, "\nSecond step : QR decomposition\n" );
      fprintf ( outfile, "****************************************************************************************************\n\n" );
      fprintf ( outfile, "Matrix H:\n" );
    }
    pdlaprnt2( &m, &n, C, &ione, &ione, descC, &izero, &izero, "H", outfile, work, len_mat_name );
    
    /* QR decomposition */
    pdgeqrf_( &m, &n, C, &ione, &ione, descC, tau, work, &lwork, &info );
    /* Copy in buffer */
    pdlacpy_ ( "All", &m, &n, C, &ione, &ione, descC, C0, &ione, &ione, descC );
    /* To print Q we need to explicit it */
    pdorgqr_( &m, &n, &ltau, C, &ione, &ione, descC, tau, work, &lwork, &info );
    if(i_am==0) fprintf ( outfile, "\nMatrix B:\n" );
    pdlaprnt2( &m, &n, C, &ione, &ione, descC, &izero, &izero, "B", outfile, work, len_mat_name );


    //****************************************************
    //          THIRD STEP : C=Q^t*A 
    //****************************************************

    if (i_am==0) {
      fprintf ( outfile, "\n****************************************************************************************************" );
      fprintf ( outfile, "\nThird step : C=Q^t*A\n" );
      fprintf ( outfile, "****************************************************************************************************\n\n" );
      fprintf ( outfile, "Matrix C:\n" );
    }
    /* C = Q^t*A */
    pdormqr_( "L", "T", &m, &n, &ltau, C0, &ione, &ione, descC, tau, A0, &ione, &ione, descA, work, &lwork, &info );
    pdlaprnt2( &n, &k, A0, &ione, &ione, descA, &izero, &izero, "C", outfile, work, len_mat_name );


    //****************************************************
    //          CHECK : A=B*C? 
    //****************************************************

    if (i_am==0) {
      fprintf ( outfile, "\n****************************************************************************************************" );
      fprintf ( outfile, "\nCheck : A~B*C\n" );
      fprintf ( outfile, "****************************************************************************************************\n\n" );
      fprintf ( outfile, "A-B*C :\n" );
    }
    /* A = A - B*C */
    pdgemm_ ( "N", "N", &m, &k, &n, &dmone, C, &ione, &ione, descC, A0, &ione, &ione, descA, &done, A, &ione, &ione, descA );
    pdlaprnt2( &m, &k, A, &ione, &ione, descA, &izero, &izero, "A", outfile, work, len_mat_name );
    resid = pdlange_( "F", &m, &k, A, &ione, &ione, descA, work );

    wall_time += Cdwalltime00( );

    if(i_am==0) {
      fprintf ( outfile, "\n||A-B*C||_F = " );
      fprintf ( outfile, "%g\n", resid );
      if(n==k && resid < eps) fprintf ( outfile, "It seems to be alright : if l == n the approximation is eps-exact !\n" );
      else fprintf ( outfile, "I seems to have a problem we should get residual == eps machine if l == n... You need to debug :-)\n" );
    }

    //****************************************************
    //                FINALIZE ALL
    //****************************************************

    free(A);
    free(B);
    free(C);
    free(work);

    Cblacs_gridexit( contxt );
  }

  if (i_am==0) {
    fprintf ( outfile, "\n****************************************************************************************************" );
    fprintf ( outfile, "\nWall_time was approximately %f seconds.\n", wall_time );
    fprintf ( outfile, "****************************************************************************************************\n\n" );
    fprintf ( outfile, "End of prototype.");
    fclose ( outfile );
  }

  Cblacs_exit( 0 );
  return 0;
}


/* prototype.c ends here */
