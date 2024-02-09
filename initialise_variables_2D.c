
#include "TP06_OpSplit_2D.h"

void initialise_variables_2D( double **u, int N )
{

  int n,m;

  /* Indices for u array */
  int V =      1;
  int M =      2;
  int H =      3;
  int J =      4;
  int R =      5;
  int S =      6;
  int D =      7;
  int F =      8;
  int F2 =     9;
  int FCass = 10;
  int Xr1 =   11;
  int Xr2 =   12;
  int Xs =    13;
  int RR =    14;
  int OO =    15;
  int CaSS =  16;
  int CaSR =  17;
  int Cai =   18;
  int Nai =   19;
  int Ki =    20;
  int fib_volt = 21;
  int rkv = 22; 
  int skv = 23;  

  /* Initial Gate Conditions */

  for (n = 1; n <= N; n++)
    {
    u[n][V] = -86.2;
    u[n][M] = 0.0;
    u[n][H] = 0.75;
    u[n][J] = 0.75;
    u[n][Xr1] = 0.0;
    u[n][Xr2] = 1.0;
    u[n][Xs] = 0.0;
    u[n][R] = 0.0;
    u[n][S] = 1.0;
    u[n][D] = 0.0;
    u[n][F] = 1.0;
    u[n][F2] = 1.0;
    u[n][FCass] = 1.0;
    u[n][RR] = 1.0;
    u[n][OO] = 0.0;
    u[n][Cai] = 0.00007;
    u[n][CaSR] = 3.0; //1.3;
    u[n][CaSS] = 0.00007;
    u[n][Nai] = 7.67;
    u[n][Ki] = 138.3;
    u[n][fib_volt] = -24.4;
    u[n][rkv] = 0;
    u[n][skv] = 1;
    }
    
}
