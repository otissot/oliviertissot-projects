// PSMATGEN : Parallel Real Double precision MATrix GENerator.
//  Generate (or regenerate) a distributed matrix A (or sub-matrix of A).
// Original fortran code from http://acts.nersc.gov/scalapack/hands-on/
// Conversion into C by Kelly McQuighan, with modifications that don't
// use the psmatgeninc.f file. 
// last modified 8/1/13
// - modifs by Olivier Tissot 5/9/2014:
//   * adding zero matrix generator
// Notes from Fortran code
//  =====
//
//  The code is originally developed by David Walker, ORNL,
//  and modified by Jaeyoung Choi, ORNL.
//
//  Reference: G. Fox et al.
//  Section 12.3 of "Solving problems on concurrent processors Vol. I"
//
//  MATRICES ARE COLUMN-MAJOR!!!
//
//  Usage
//  ======
//
// psmatgen( 	double* A, 
//		int* DESC,
//		char* AFORM, 
//		char* DIAG, 		
//		int IROFF, 
//		int IRNUM, 
//		int ICOFF, 
//		int ICNUM, 
//		int MYROW, 
//		int MYCOL, 
//		int NPROW, 
//		int NPCOL)
//
//  Arguments
//  =========
//
//  AFORM   (global input) CHARACTER*1
//          if AFORM = "S" : A is returned is a symmetric matrix.
//          if AFORM = "H" : A is returned is a Hermitian matrix.
//                     "I" : ones along diagonal given by DIAG is returned
//          otherwise a random matrix is generated.
//          UNSUPPORTED OPTIONS if AFORM = "T" : A is overwritten with the transpose of
//                                               what would normally be generated.
//                              if AFORM = "C" : A is overwritten with the conjugate trans-
//                                               pose of what would normally be generated.
//
//  DIAG    (global input) CHARACTER*1
//          if DIAG = "D" : A is diagonally dominant.
//          if DIAG = "N" : A is not diagonally dominant.
//          if AFORM = "I" then DIAG can also take any number between (-m+1) < DIAG < (m-1)
//            note that DIAG must still be enclosed in quotes. Also, DIAG < 0 implies
//            ones will be on a subdiagonal, whereas DIAG > 0 implies ones will be on 
//            a superdiagonal.
//
//  DESC    (global input)
//	    Description of the local matrix. Contains 9 values as follows:
//	    DESC[0] = DTYPE
//	    DESC[1] = CTXT
//	    DESC[2] = M
//	    DESC[3] = N
//	    DESC[4] = MB
//	    DESC[5] = NB
//	    DESC[6] = RSRC
//	    DESC[7] = CSRC
//	    DESC[8] = LLD
//
//  M       (global input) INTEGER
//          The number of rows in the generated distributed matrix.
//
//  N       (global input) INTEGER
//          The number of columns in the generated distributed
//          matrix.
//
//  MB      (global input) INTEGER
//          The row blocking factor of the distributed matrix A.
//
//  NB      (global input) INTEGER
//          The column blocking factor of the distributed matrix A.
//
//  LLD     (local input) INTEGER
//          The leading dimension of the array containing the local
//          pieces of the distributed matrix A.
//
//  A       (local output) REAL            , pointer into the local
//          memory to an array of dimension ( LDA, * ) containing the
//          local pieces of the distributed matrix.
//
//  IA   (global input) INTEGER
//       The global initial row index of the submatrix to be generated
//       Note: as written, this function only supports if IA is a multiple
//       of mb
//
//  JA   (global input) INTEGER
//          The global initial column index of the submatrix to be generated
//          Note: as written, this function only supports if JA is a multiple
//          of nb
//
//  IAROW   (global input) INTEGER
//          The row processor coordinate which holds the first block
//          of the distributed matrix A.
//
//  IACOL   (global input) INTEGER
//          The column processor coordinate which holds the first
//          block of the distributed matrix A.
//
//  IRNUM   (local input) INTEGER
//          The number of local rows to be generated.
//
//  ICNUM   (local input) INTEGER
//          The number of local columns to be generated.
//
//  MYROW   (local input) INTEGER
//          The row process coordinate of the calling process.
//
//  MYCOL   (local input) INTEGER
//          The column process coordinate of the calling process.
//
//  NPROW   (global input) INTEGER
//          The number of process rows in the grid.
//
//  NPCOL   (global input) INTEGER
//          The number of process columns in the grid.
//
//  =====================================================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <blacs_headers.h>
#include <pblas_headers.h>
#include <scalapack_tools_headers.h>
#include "psmatgen.h"
#include "stringTools.h"
#include "pblasIOtools.h"
//#include "psmatgeninc.h"
//#include "mkl_scalapack.h"

static int imax( int a, int b ){
  if (a>b) return(a); else return(b);
}

static double dabs( double a) {
  if (a>0) return a; else return (-1*a);
}

void psmatgen( char* aform, char* diag,  double* A, int ia, int ja, int* desc, int msub, int nsub, int myrow, int mycol, int nprow, int npcol, int iseed, int* status){
  int m, n, mb, nb, lld, rsrc, csrc, i_am, moffglob, noffglob, irnum, icnum;
  int mp, SYM, HERM, DIAGDOM, NDIAG, NOTRAN, EYE, ZERO, info = 0, infoDesc=0;
  int mrrow, mrcol, moff, noff, mend, nend;
  int jk, ic, ioffc, i, ir, ik, ioffr, j, offset=0;
  int mendglob, nendglob;
  double one = 1.0e0, two = 2.0e0, zero=0.0e0, maxmn;
  double dnb, dmb, dia, dja, dmsub, dnsub, dnprow, dnpcol;
  double dmoffglob, dnoffglob, dmendglob, dnendglob;
  int extramb, extranb;

  i_am = myrow*nprow+mycol;
  // get necessary variables from desc
  //ctxt = desc[1];
  m = desc[2];
  n = desc[3];
  mb = desc[4];
  nb = desc[5];
  rsrc = desc[6];
  csrc = desc[7];
  lld = desc[8];
  ia--; ja--; // because indeces for C start at zero

  // test the input arguments
  // rsrc gives the process row over which the first row of the matrix is distributed
  // csrc gives the process column over which the first column of the matrix is distributed
  mp = numroc_( &m, &mb, &myrow, &rsrc, &nprow );
  //nq = numroc_(&n, &nb, &mycol, &csrc, &npcol );

  SYM = ( (*aform == 'S') || (*aform == 's') );
  HERM = ( (*aform == 'H') || (*aform == 'h') );
  EYE = ( (*aform == 'I') || (*aform == 'i') );
  ZERO = ( (*aform == 'Z') || (*aform == 'z') );
  NOTRAN = ( (*aform == 'N') || (*aform == 'n') );
  DIAGDOM = ( (*diag == 'D') || (*diag == 'd') );
  NDIAG = ( (*diag == 'N') || (*diag == 'n') );

  if ( (mb != nb )&&(SYM || HERM) ) {
    if ( i_am==0 ) fprintf ( stdout, "Symmetric matrices with rowNB not equal to colNB is not supported! Creating random matrix instead.\n" );
    SYM = 0; HERM = 0; EYE = 0; NOTRAN = 1;
  }
  if ( (m != n )&& EYE ) {
    if ( i_am==0 ) fprintf ( stdout, "Identity matrix with N not equal to M is not supported! Creating random matrix instead.\n" );
    SYM = 0; HERM = 0; EYE = 0; NOTRAN = 1;
  }
  if ( !DIAGDOM && !NDIAG ) offset=atoi(diag);

  if (!(SYM || HERM || NOTRAN || EYE || ZERO) ) info=1;
  else if (!DIAGDOM && !NDIAG && !EYE ) info=2;
  else if ( EYE && (abs(offset) >= m) ) info=2;
  else if ( SYM || HERM ) {
    if (m != n) { info=6; infoDesc=4; }
    else if (mb != nb) {info=6; infoDesc=6; }
  }
  else if (m < 0) { info = 6; infoDesc=3; }
  else if (n < 0) { info = 6; infoDesc=4; }
  else if (mb < 1) { info = 6; infoDesc=5; }
  else if (nb < 1) { info = 6; infoDesc=6; }
  else if (lld < 0 ) { info = 6; infoDesc=9; }
  else if ( (rsrc < 0) || (rsrc >= nprow) ) { info=6; infoDesc=7; }
  else if ( (csrc < 0) || (csrc >= npcol) ) { info=6; infoDesc=8; }
  else if ( (ia % mb) > 0) info=4;
  else if ( msub > (m-ia) ) info=7;
  else if ( (ja % mb) > 0) info=5;
  else if ( nsub > (n-ja) ) info=8;
  else if ( (myrow < 0) || (myrow >= nprow) ) info=9;
  else if ( (mycol < 0) || (mycol >= npcol) ) info=10;

  if ( info != 0) {
    if ( info == 6 ) {
      fprintf ( stdout, "Processor {%d, %d}: On entry to psmatgen.c, parameter number %d, which is an array, had an illegal value.\n", myrow, mycol, info ) ;
      fprintf ( stdout, "In particular, the %dth element of DESC is invalid.\n", infoDesc );
    } else if ( info == 1 ) { 
      if ( (myrow==rsrc) && (mycol==csrc) ){
      fprintf ( stdout, "Processor {%d, %d}: On entry to psmatgen.c, parameter number %d had an illegal value.\n", myrow, mycol, info ) ;
      fprintf ( stdout, "Currently supported options are: ''N'', ''S'', ''H'', ''I'', ''Z''.\n");
      }
    } else if ( info == 2 ) {
      if ( (myrow==rsrc) && (mycol==csrc) ){
      fprintf ( stdout, "Processor {%d, %d}: On entry to psmatgen.c, parameter number %d had an illegal value.\n", myrow, mycol, info ) ;
      fprintf ( stdout, "Currently supported options are: ''D'', ''N''.\n");
      fprintf ( stdout, "If AFORM=''I'' then DIAG may take a number between (-m+1) < DIAG < (m-1).\n");
      fprintf ( stdout, "Note that the values of DIAG must still be passed as a string: ''DIAG''.\n");
      }
    } else {
      fprintf ( stdout, "Processor {%d, %d}: On entry to psmatgen.c, parameter number %d had an illegal value.\n", myrow, mycol, info ) ;
    }
    *status = -1;
    return;
  }

  dnb = (double) nb; dmb = (double) mb;
  dia = (double) ia; dja = (double) ja; 
  dmsub = (double) msub; dnsub = (double) nsub;
  dnprow = (double) nprow; dnpcol = (double) npcol;

  // parameter values are valid so set some local constants
  mrrow = (nprow + myrow - rsrc) % nprow; // my process row within matrix distribution
  mrcol = (npcol + mycol - csrc) % npcol; // my process column within matrix distribution
  moffglob = floor ( dia / dmb );
  noffglob = floor ( dja / dnb );
  mendglob = ceil ( (dia+dmsub) / dmb );
  nendglob = ceil ( (dja+dnsub) / dnb );
 // double dmendglob, dnoffglob, dnendglob, dmoffglob;
  dmoffglob = (double) moffglob; dnoffglob = (double) noffglob;
  dmendglob = (double) mendglob; dnendglob = (double) nendglob;
  //  fprintf(stdout, "I am {%d,%d}. Moffglob=%d, noffglob=%d, mendglob=%d, nendglob=%d.\n", myrow, mycol, moffglob, noffglob, mendglob, nendglob);

  moff   = floor( dmoffglob / dnprow ); // my block number for row offset
  noff   = floor( dnoffglob / dnpcol ); // my block number for column offset
  if ( mrrow < ( moffglob % nprow) ) moff++;
  if ( mrcol < ( noffglob % npcol) ) noff++;
  mend   = ceil( dmendglob / dnprow ); // ending row block index, local # of blocks
  nend   = ceil( dnendglob / dnpcol ); // end column block index, local # of blocks
  if ( mrrow > ( (mendglob+nprow-1) % nprow ) ) mend--;
  if ( mrcol > ( (nendglob+npcol-1) % npcol ) ) nend--;

  irnum = (mend-moff)*mb; icnum = (nend-noff)*nb;
  if ( (msub % mb) != 0 ) {
    extramb = ceil(dmsub/dmb);
    if ( ( (extramb+nprow-1) % nprow ) == mrrow ) irnum+= (msub % mb);
  }
  if ( (nsub % nb) != 0 ) {
    extranb = ceil(dnsub/dnb);
    if ( ( (extranb+npcol-1) % npcol ) == mrcol ) icnum+= (nsub % nb);
  }

  //  fprintf(stdout, "I am {%d, %d}. Irnum=%d, icnum=%d. Moff=%d, mend=%d, noff=%d, nend=%d.\n", myrow, mycol, irnum, icnum, moff, mend, noff, nend);

  // symmetric of Hermitian matrix will be generated
  if ( SYM || HERM ) {

    // first generate low triangular part

    // first generate lower triangular part


    jk = 0; // keeps track of column
    for ( ic = noff; ic < nend; ic++ ){
      if (jk > icnum) break;
      ioffc = (ic*npcol+mrcol)*nb;
      for (i = 0; i < nb; i++ ) {
        if (jk > icnum) break;
        ik = 0; //keeps track of rows
        for ( ir = moff; ir < mend; ir++ ){
          if (ik > irnum ) break;
          ioffr = (ir*nprow+mrrow) * mb;
          for ( j = 0; j < mb; j++ ) {
            if (ik > irnum ) break;
            if ( (ioffr+j+offset) >= (ioffc+i) ) {
              srand( iseed+ioffr+i+m*(ioffc+j) );
              *(A+(jk+noff*nb)*mp+ik+moff*mb) = one - two*rand()/RAND_MAX;
            } else {
              srand( iseed+(ioffr+i)*n+ioffc+j );
              *(A+(jk+noff*nb)*mp+ik+moff*mb) = one - two*rand()/RAND_MAX;
            }
            ik=ik+1;
          }
        }
        jk = jk+1;
      }
    }

  } // end if symmetric matrix

  // a random matrix is generated
  else if ( EYE ) {

    jk = 0; // keeps track of column
    for ( ic = noff; ic < nend; ic++ ){
      if (jk > icnum) break;
      ioffc = (ic*npcol+mrcol)*nb;
      for (i = 0; i < nb; i++ ) {
        if (jk > icnum) break;
        ik = 0; //keeps track of rows
        for ( ir = moff; ir < mend; ir++ ){
          if (ik > irnum ) break;
          ioffr = (ir*nprow+mrrow) * mb;
          for ( j = 0; j < mb; j++ ) {
            if (ik > irnum ) break;
            if ( (ioffr+j+offset) == (ioffc+i) ) {
              *(A+(jk+noff*nb)*mp+ik+moff*mb) = one;
            } else {
              *(A+(jk+noff*nb)*mp+ik+moff*mb) = zero;
            }
            ik=ik+1;
          }
        }
        jk = jk+1;
      }
    }

  } // end identity matrix
  else if( ZERO ) {

    jk = 0; // keeps track of column
    for ( ic = noff; ic < nend; ic++ ){
      if (jk > icnum) break;
      ioffc = (ic*npcol+mrcol)*nb;
      for (i = 0; i < nb; i++ ) {
        if (jk > icnum) break;
        ik = 0; //keeps track of rows
        for ( ir = moff; ir < mend; ir++ ){
          if (ik > irnum ) break;
          ioffr = (ir*nprow+mrrow) * mb;
          for ( j = 0; j < mb; j++ ) {
            if (ik > irnum ) break;
            *(A+(jk+noff*nb)*mp+ik+moff*mb) = zero;
            ik=ik+1;
          }
        }
        jk = jk+1;
      }
    }
  } //end else: zero matrix
  else {

    // initialize random number generator
    srand(i_am+iseed);

    jk = 0; // keeps track of column
    for ( ic = noff; ic < nend; ic++ ){
      if (jk > icnum) break;
      ioffc = (ic*npcol+mrcol)*nb;
      for (i = 0; i < nb; i++ ) {
        if (jk > icnum) break;
        ik = 0; //keeps track of rows
        for ( ir = moff; ir < mend; ir++ ){
          if (ik > irnum ) break;
          ioffr = (ir*nprow+mrrow) * mb;
          for ( j = 0; j < mb; j++ ) {
            if (ik > irnum ) break;
            *(A+(jk+noff*nb)*mp+ik+moff*mb) = one - two*rand()/RAND_MAX;
            ik=ik+1;
          }
        }
        jk = jk+1;
      }
    }
  } //end else: random matrix
      
  // diagonally dominant matrix will be generated
  if ( DIAGDOM ) {

    if (mb != nb ) {
      if (i_am==0) fprintf ( stdout, "Diagonally dominant matrices with rowNB not equal to colNB is not supported! Command ignored.\n" );
      return;
    }

    if (EYE) {
      if (i_am==0) fprintf ( stdout, "Identity matrix is automatically diagonally dominant! Command ignored.\n" );
      return;
    }
  // since the magnitude of each element is less than 1, adding the max( n, m ) to the diagonal is sufficient
  // Def: diagonally dominant if the magnitude of the diagonal element in a row (column) is greater than 
  //      the sum of magnitude of all other elements in the row (column). In other words,
  //      ||a_{ii}|| >= \Sum_{j} a_{ij}, \Sum_{j} a_{ji}
    maxmn = (double) imax(m,n);
    jk = 0; // keeps track of column
    for ( ic = noff; ic < nend; ic++ ){
      if (jk > icnum) break;
      ioffc = (ic*npcol+mrcol)*nb;
      for (i = 0; i < nb; i++ ) {
        if (jk > icnum) break;
        ik = 0; //keeps track of rows
        for ( ir = moff; ir < mend; ir++ ){
          if (ik > irnum ) break;
          ioffr = (ir*nprow+mrrow) * mb;
          for ( j = 0; j < mb; j++ ) {
            if (ik > irnum ) break;
            if ( (ioffr+j+offset) == (ioffc+i) ) {
//              fprintf(stdout, "I am {%d,%d}. ic=%d, ir=%d, ioffr=%d, ioffc=%d.\n", myrow, mycol, ic, ir, ioffr, ioffc);
              *(A+(jk+noff*nb)*mp+ik+moff*mb) = dabs(*(A+(jk+noff*nb)*mp+ik+moff*mb)) + maxmn;
            }
            ik=ik+1;
          }
        }
        jk = jk+1;
      }
    }

  }

  return;
}
