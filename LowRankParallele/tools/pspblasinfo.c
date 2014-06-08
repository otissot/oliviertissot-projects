// -- PBLAS Example code --
// Original fortran code from http://acts.nersc.gov/scalapack/hands-on/
// Conversion into C by Kelly McQuighan,
// last modified 7/31/13
// - modifs by Olivier Tissot 5/9/2014 :
//   * fix bug : bad syntax for broadcast receive
// reads system parameters from the file BLAS.dat

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <blacs_headers.h>
#include <pblas_headers.h>
#include "pspblasinfo.h"
#include <string.h>
#include "stringTools.h"
#include "pblasIOtools.h"
#include "nameFile.h"

// returns -1 if unsuccessful, 0 if successful
void pspblasinfo( int* m, int* n, int* k, int* block_size, int* num_proc_row, int* num_proc_col, int* sys_param, int* i_am, int* num_procs, char folder[], const char * infile_name, char* outfile_name, char* usr_info, int* outfile_name_len, int* status ){

  char line[100], info[100], trash[100];
  char suffix_out[]="out", file_prefix[100];
  int contxt, i;
  FILE *infile;

  Cblacs_get( -1, 0, &contxt );
  Cblacs_gridinit( &contxt, "Row-major", 1, *num_procs );

  // processor 0 reads data from PBLAS.dat and broadcasts to other processors.
  if ( *i_am==0 ) {
    infile = fopen(infile_name, "r");
    if (infile == NULL) {
      fprintf(stdout, "***********************************************************\n");
      fprintf(stdout, "Error: input file %s does not exist !\n", infile_name);
      fprintf(stdout, "***********************************************************\n");
      *status = -1;
      Cigebs2d( contxt, "All", " ", 1, 1, status, 1 );
      return;
    }
    // header line
    if (fgets(line, 100, infile)==NULL ) {
      fprintf(stdout, "\n***********************************************************\n\n");
      fprintf(stdout, "Error: input file does not contain any data!\n");
      fprintf(stdout, "\n***********************************************************\n\n");
      *status = -1;
      Cigebs2d( contxt, "All", " ", 1, 1, status, 1 );
      return;
    } 
    // user info
    if (fgets(usr_info, 100, infile)==NULL ) {
      fprintf(stdout, "\n***********************************************************\n\n");
      fprintf(stdout, "Error: input file does not contain any data!\n");
      fprintf(stdout, "\n***********************************************************\n\n");
      *status = -1;
      //Cigebs2d( contxt, "All", " ", 1, 1, status, 1 );
      return;
    }
    // output file name
    if (fgets(line, 100, infile)==NULL )  {
      fprintf(stdout, "\n***********************************************************\n\n");
      fprintf(stdout, "Error: input file should contain the output file name on the 3rd line!\n");
      fprintf(stdout, "\n***********************************************************\n\n");
      *status = -1;
      Cigebs2d( contxt, "All", " ", 1, 1, status, 1 );
      return;
    }
    if ( sscanf(line,"%[^\t]%[^\n]", file_prefix, trash) < 1 ) {
      fprintf(stdout, "\n***********************************************************\n\n");
      fprintf(stdout, "Error: input file not formatted correctly! 3rd line should read (output file name) (tab) (additional text).\n");
      fprintf(stdout, "\n***********************************************************\n\n");
      *status = -1;
      Cigebs2d( contxt, "All", " ", 1, 1, status, 1 );
      return;
    }

    //  system parameters
    for (i=0; i<6; i++ ) {

      if (fgets(line, 100, infile) == NULL ) {
        fprintf(stdout, "\n***********************************************************\n\n");
        fprintf(stdout, "Error: input file not formatted correctly! Lines 4-13 should contain Scalapack matrix information in the format\n\n(value of m) (tab) (additional text)\n");
        fprintf(stdout, "in the following order: m, n, nrhs, mb, nb, nb_rhs, max_lldA, max_lldB, nprow, npcol\n");
        fprintf(stdout, "\n***********************************************************\n\n");
        *status = -1;
        Cigebs2d( contxt, "All", " ", 1, 1, status, 1 );
      return;
      }
      if ( sscanf(line,"%[^\t]%[^\n]", info, trash) < 2 ) {
        fprintf(stdout, "\n***********************************************************\n\n");
        fprintf(stdout, "Error: input file not formatted correctly on %dth line! It should read\n\n(value) (tab) (additional text)\n", i+4);
        fprintf(stdout, "\n***********************************************************\n\n");
        *status = -1;
        Cigebs2d( contxt, "All", " ", 1, 1, status, 1 );
        return;
      }
      sys_param[i] = atoi(info);
    }

    // completed successfully. broadcast parameters to other processes
    *status = 0;
    Cigebs2d( contxt, "All", " ", 1, 1, status, 1 );  
    Cigebs2d( contxt, "All", " ", 6, 1, sys_param, 6);
    
    fclose(infile);

  } else {
 
    // receive system variables and store in corresponding variables
    Cigebr2d( contxt, "All", " ", 1, 1, status, 1, 0, 0 );
    if ( *status < 0 ) return;
    Cigebr2d( contxt, "All", " ", 6, 1, sys_param, 6, 0, 0 );
  
  } // end if-then-else

  // store system parameters in their appropriate values
  *m            = sys_param[0];
  *n            = sys_param[1];
  *k            = sys_param[2];
  *block_size   = sys_param[3];
  *num_proc_row = sys_param[4];
  *num_proc_col = sys_param[5];
  
  if (*i_am == 0) {
    nameFile ( folder, file_prefix, "", suffix_out, outfile_name );
  }
  Cblacs_gridexit( contxt );

  return;
}
