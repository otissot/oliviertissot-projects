#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <blacs_headers.h>
#include <pblas_headers.h>
#include <scalapack_tools_headers.h>
#include <scalapack_headers.h>
#include "matinit.h"

//****************************************
// MATINIT generates and distributes matrices A and B
// depicted in Figure 2.1 and 2.2 from 
// http://www.netlib.org/scalapack/slug/node28.html#SECTION04230000000000000000
// to a 2x3 process grid
//
// Converted into C by Kelly McQuighan
// 8/2/13

//***************************************
void matinit( double* AA, int* descA, double* B, int* descB, int myrow, int mycol )
{
double S, C, A, L, P, K;
int lldA, lldB, i;

  S=19.0e0;
  C=3.0e0;
  A=1.0e0;
  L=12.0e0;
  P=16.0e0;
  K=11.0e0;

  lldA = descA[8];
  lldB = descB[8];
  for (i=0; i<lldB; i++) *(B+i) = 0.0e0;

  if ( ( myrow == 0 ) && ( mycol == 0 ) ) {
    *( AA+0        ) =  S;
    *( AA+1        ) = -S;
    *( AA+2        ) = -S;
    *( AA+3        ) = -S;
    *( AA+4        ) = -S;
    *( AA+0+  lldA ) =  C;
    *( AA+1+  lldA ) =  C;
    *( AA+2+  lldA ) = -C;
    *( AA+3+  lldA ) = -C;
    *( AA+4+  lldA ) = -C;
    *( AA+0+2*lldA ) =  A;
    *( AA+1+2*lldA ) =  A;
    *( AA+2+2*lldA ) =  A;
    *( AA+3+2*lldA ) =  A;
    *( AA+4+2*lldA ) = -A;
    *( AA+0+3*lldA ) =  C;
    *( AA+1+3*lldA ) =  C;
    *( AA+2+3*lldA ) =  C;
    *( AA+3+3*lldA ) =  C;
    *( AA+4+3*lldA ) = -C;
  } else if ( ( myrow == 0 ) && ( mycol == 1 ) ) {
    *( AA+0        ) =  A;
    *( AA+1        ) =  A;
    *( AA+2        ) = -A;
    *( AA+3        ) = -A;
    *( AA+4        ) = -A;
    *( AA+0+  lldA ) =  L;
    *( AA+1+  lldA ) =  L;
    *( AA+2+  lldA ) = -L;
    *( AA+3+  lldA ) = -L;
    *( AA+4+  lldA ) = -L;
    *( AA+0+2*lldA ) =  K;
    *( AA+1+2*lldA ) =  K;
    *( AA+2+2*lldA ) =  K;
    *( AA+3+2*lldA ) =  K;
    *( AA+4+2*lldA ) =  K;
  } else if ( ( myrow == 0 ) && ( mycol == 2 ) ) {
    *( AA+0        ) =  A;
    *( AA+1        ) =  A;
    *( AA+2        ) =  A;
    *( AA+3        ) = -A;
    *( AA+4        ) = -A;
    *( AA+0+  lldA ) =  P;
    *( AA+1+  lldA ) =  P;
    *( AA+2+  lldA ) =  P;
    *( AA+3+  lldA ) =  P;
    *( AA+4+  lldA ) = -P;
  } else if ( ( myrow == 1 ) && ( mycol == 0 ) ) {
    *( AA+0        ) = -S;
    *( AA+1        ) = -S;
    *( AA+2        ) = -S;
    *( AA+3        ) = -S;
    *( AA+0+  lldA ) = -C;
    *( AA+1+  lldA ) = -C;
    *( AA+2+  lldA ) = -C;
    *( AA+3+  lldA ) =  C;
    *( AA+0+2*lldA ) =  A;
    *( AA+1+2*lldA ) =  A;
    *( AA+2+2*lldA ) =  A;
    *( AA+3+2*lldA ) = -A;
    *( AA+0+3*lldA ) =  C;
    *( AA+1+3*lldA ) =  C;
    *( AA+2+3*lldA ) =  C;
    *( AA+3+3*lldA ) =  C;

    *B                   = 1.0e0;
  } else if ( ( myrow == 1 ) && ( mycol == 1 ) ) {
    *( AA+0        ) =  A;
    *( AA+1        ) = -A;
    *( AA+2        ) = -A;
    *( AA+3        ) = -A;
    *( AA+0+  lldA ) =  L;
    *( AA+1+  lldA ) =  L;
    *( AA+2+  lldA ) = -L;
    *( AA+3+  lldA ) = -L;
    *( AA+0+2*lldA ) =  K;
    *( AA+1+2*lldA ) =  K;
    *( AA+2+2*lldA ) =  K;
    *( AA+3+2*lldA ) =  K;
  } else if ( ( myrow == 1 ) && ( mycol == 2 ) ) {
    *( AA+0        ) =  A;
    *( AA+1        ) =  A;
    *( AA+2        ) = -A;
    *( AA+3        ) = -A;
    *( AA+0+  lldA ) =  P;
    *( AA+1+  lldA ) =  P;
    *( AA+2+  lldA ) = -P;
    *( AA+3+  lldA ) = -P;
  }

return;
}
