/***************************************************************

  checkpoint_read.c

  This module is a component of the FK4V 3D package

  Version 1.0

  (c) Richard Clayton, University of Sheffield 2005
      r.h.clayton@sheffield.ac.uk


  reads parameter space from binary file "checkpointXXXX.out"

***************************************************************/
#include "TP06_OpSplit_2D.h"

int checkpoint_read( double **u, double *time, int *t, int count, int N)
{
  int elements_to_read, i;
  int n, m, M;

  double u_n_m;
  char fname[80];

  FILE *chkpt_file;

  M = NUM_STATES;
  sprintf( fname,"%s%06d.out",CHKPTROOT,count );
  printf("opening checkpoint file %s\n",fname);
  chkpt_file = fopen( fname, "rb" );

  fread( time, sizeof(double),1,chkpt_file);
  fread( t, sizeof(int),1,chkpt_file);

  printf("read time = %g, t = %d\n",*time,*t);
  elements_to_read = N * M;
  printf("reading %d elements\n",elements_to_read);

  i = 0;
  u_n_m = 0.0;
  for (n = 1; n <= N; n++)
  {
	for (m = 1; m <= M; m++)
    {
      i += fread( &u_n_m, sizeof(double), 1, chkpt_file );
      u[n][m] = u_n_m;
    }
  }

  printf("read %d elements\n", i);
  if (ferror(chkpt_file)) perror("error reading data");

  fclose(chkpt_file);
  return (1);
}
