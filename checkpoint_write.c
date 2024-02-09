/***************************************************************

  checkpoint_write.c

  This module is a component of the FK4V 3D package

  Version 1.0

  (c) Richard Clayton, University of Sheffield 2005
      r.h.clayton@sheffield.ac.uk


  writes parameter space to binary file "checkpointXXXX.out"
  includes workaround to avoid problem associate with writing
  more than 32767 elements

***************************************************************/
#include "TP06_OpSplit_2D.h"

int checkpoint_write( double **u, double time, int t, int count, int N )
{
  int n, m, M;
  int elements_to_write, i;

  double u_n_m;

  char fname[80];
  FILE *chkpt_file;

  M = NUM_STATES;
  sprintf( fname,"%s%06d.out",CHKPTROOT,count );
  printf("opening checkpoint file %s\n",fname);
  chkpt_file = fopen( fname, "wb" );

  fwrite( &time, sizeof(double), 1, chkpt_file );
  fwrite( &t, sizeof(int), 1, chkpt_file );
  printf("write time = %g, t = %d\n",time,t);

  elements_to_write = N * M;
  printf("writing %d elements\n", elements_to_write);

  i = 0;
  for (n = 1; n <= N; n++)
  {
	for (m = 1; m <= M; m++)
    {
	   u_n_m = u[n][m];
	   i += fwrite(&u_n_m, sizeof(double), 1, chkpt_file );
    }
  }
  printf("written %d elements\n", i);
  if (ferror(chkpt_file)) perror("error writing data");

  fclose(chkpt_file);
  return (1);
}
