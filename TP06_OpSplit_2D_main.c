
/********************************************************************

 TP06_OpSplit_2D_main.c

 Author        S.Sridhar and R.H. Clayton (s.seshan@sheffield.ac.uk)

 This program evaluates the dynamics coupled myocyte-fibroblast motif-2 
 
The parameters for numerical solution are contained in the header file
*********************************************************************/

#include "TP06_OpSplit_2D.h"

void writeData(char *fname, double *dataToWrite, int **geom, int nrows, int ncols);

int main(int argc, char **argv)
{
  /* constants */
  const int tmax = BCL1*20*1000; //1000 for dt = 0.001 NUM_ITERATIONS;			// duration of simulation
  const int num_states = NUM_STATES;        // number of states stored at each grid point
  const int num_params = NUM_PARAMS;        // number of model parameters (inputs)
  const int voltage_steps = VOLTAGE_STEPS;  // voltage steps in the lookup table
  const int num_lookup = NUM_LOOKUP;		// number of variables in lookup table
  const int nrows = ROWS;
  const int ncols = COLUMNS;
  const int RC = ROWS * COLUMNS;            // total number of grid points
  const int V = 1;                          // index of membrane voltage
  const int fib_volt =21;
  const double dtlong = DT;					// long time step for diffusion
  const double half_dtlong = DT/2.0;		// half time step for diffusion
  //const double D = DIFFUSION;				// diffusion coefficient
  const double dx2 = DX*DX;					// dx squared
    
  /* variables */
  int N;
  int **geom, **nneighb;    				// arrays to store geometry and nearest neighbours
  int t, n, m, dummy;						// array indices
  int k, ko, kmax;					        // parameters for adaptive timestep
  int row, col;								// array indices
  int *celltype;							// specify myocyte (1) or fibroblast (0)
  int stfcount = 0;						    // index for stf output

  double dtshort;							// adaptive short time step for ODE solution
  double **u, **lookup;						// arrays for storing model state and lookup table
  double time;
  double timems;
  double *stimCurrent;// = 0.0;
  double *dVdt, *new_Vm, *old_Vm, *U, *U1;		// arrays for storing state during updates
  double dummy1, dummy2;
  double dV1, dVf1, fib_to_myo1, myo_to_fib1;
  double dV2, dVf2, fib_to_myo2, myo_to_fib2;
  double *params;
  double *D;                               // array to hold local diffusion coefficient
  double stimCurrent1, stimCurrent2;  
    
  const double bcl1 = BCL1; 
  const double bcl2 = BCL1;                    // basic cycle length for pacing
  //const int numBeats = NUMBEATS;           // number of S1 stimuli
  //const double s1s2 = S1S2;                // coupling interval to s2 stimulus
  int s1Beat = 0;
    
  double *threshold;                       // array to store threshold for APD90 detection
  double *maxV, *minV;                     // max and min voltage at each point
  double *upStrokeTimeS1;                  // array to store upstroke times for final S1 beat
  double *downStrokeTimeS1;                // array to store downstroke times for final S1 beat
  double *upStrokeTimeS2;                  // array to store upstroke times for S2 beat
  double *downStrokeTimeS2;                // array to store downstroke times for S2 beat
    
  FILE *egPtr; 
  char outputFile[80];					   // filename for outputs

 /*declare and initialise fibroblasts parameters */
 double myo_to_fib_long1, fib_to_myo_long1, myo_to_fib_long2, fib_to_myo_long2 ;
 int NUMFIB = 4;
 double FIB_CAP = 50; //50
 double MYO_CAP = 185; //185
 double G_str1GAP = 0;
 double G_str2GAP = G_STRONG1;
 double G_wk1GAP = G_WEAK1;
 double G_wk2GAP = 0; //G_WEAK1;
 int start1 = START; // Starting time of stimulation at cell 1
 int delay1 = DELAY*1000; // Delay in starting time of stimulation at cell 2
 //printf("done\n");
 
  /* Create geometry and nearest neighbour arrays */
  geom = imatrix( 1, nrows, 1, ncols );
  nneighb = imatrix( 1, RC, 1, 8 );
  //printf("initialising geometry arrays\n");
  N = 2 ; //initialise_geometry_2D( geom, nrows, ncols, nneighb );  
    
  /* Initialise arrays */
  u = fmatrix( 1, N, 1, num_states );
  lookup = fmatrix( 0, num_lookup, 0, voltage_steps );
  dVdt = fvector( 1, N );
  new_Vm = fvector( 1, N );
  old_Vm = fvector( 1, N );
  U = fvector(1, num_states);
  U1 = fvector(1, num_states);
  params = fvector(1, num_params);
  celltype = ivector(1, N);
  D = fvector(1, N);
  stimCurrent = fvector(1, N);  
  /* arrays for apd90 detection */
  maxV = fvector(1, N);
  minV = fvector(1, N);
  threshold = fvector(1, N);
  upStrokeTimeS1 = fvector(1, N);
  downStrokeTimeS1 = fvector(1, N);
  upStrokeTimeS2 = fvector(1, N);
  downStrokeTimeS2 = fvector(1, N);

  
  /* open files for output */
  sprintf(outputFile,"%sVm.txt",OUTPUTFILEROOT);
  egPtr = fopen(outputFile,"w");
  
  /* Initialise model state and paramaters */
  //printf("initialising ...\n");
  initialise_variables_2D( u, N );
  //initialise_parameters_2D( params );
 // printf("done\n");
    
  
  /* set celltype to be 1 throughout  and stimCurrent to 0 throughout*/
  m = 0;
  for (n = 1; n <= N; n++)
     {
      celltype[n] = 1;
      stimCurrent[n]=0.0;
      }
      
    
  /* initialise upstroke and downstroke arrays */
  for (n = 1; n <= N; n++)
    {
      maxV[n] = -90.0;
      minV[n] = 90.0;
      threshold[n] = -73.0;
      upStrokeTimeS1[n] = 0.0;
      downStrokeTimeS1[n] = 0.0;
      upStrokeTimeS2[n] = 0.0;
      downStrokeTimeS2[n] = 0.0;
    }
    
  /* Create lookup table */
  //dummy = create_TP06_lookup_OpSplit_2D( lookup );
//  printf("done\n");

  /* last bit of initialisation */
  t = 0;
  time = 0.0;
    
  if (CHKPT_READ)
    {
  	dummy = checkpoint_read( u, &time, &t, CHKPT_READ_TIME, N );
  	stfcount = ceil(time);
  	t = t + 1;
    }

/******************************************/
/*               Main loop                */
/******************************************/
  //printf("entering main loop\n");

  while (t < tmax)
	  {
      time = t * dtlong;
      t++;

            
            if (t%(int)(bcl1/dtlong) > start1 && t%(int)(bcl1/dtlong) <=(start1+1000)) //duration = 10 iterations
                 stimCurrent1 = -52.0;//-52.0;
            else stimCurrent1 = 0.0;  
                       
           if (t%(int)(bcl2/dtlong) > (start1+delay1) && t%(int)(bcl2/dtlong) <=(start1+delay1+1000))
                stimCurrent2 =  -52.0;
           else stimCurrent2  = 0.0;
          
               

          // uncomment these lines to implement adaptive time step
          // this implementation provides good agreement with standard scheme for dt=0.01 ms
          // apart from delay of ~0.1 ms in onset of AP upstroke
   	      //if (dVdt[n] > 0.01) ko = 5; else ko = 1;
             // kmax = ko + floor(fabs(dVdt[n]) * 20.0); // kmax varies between ko and 10
             //if (kmax > ceil(dtlong/0.00005)) //0.01 
            // kmax = 20;
             //kmax = dtlong/0.00005; //0.01

		  // comment this line to implement adaptive time step
		   kmax = 10;    

		  dtshort = dtlong / (double) kmax;

		  /* store state of current point in U temporarily*/
		  for (m = 1; m <= num_states; m++){
                     U[m] = u[1][m];U1[m]=u[2][m];
                     }

          /* integrate ODEs using Rush and Larsen scheme */
         
	    for (k = 1; k <= kmax; k++)
	    {
                    
                  dV1 = calculate_TP06_current_OpSplit( U, dtshort, lookup, celltype[1], stimCurrent1 );
                  dVf1 = calculate_MacFib_current( U, dtshort);  
                  
                  myo_to_fib1 = (G_str1GAP/MYO_CAP) * (U[fib_volt] - U[V]);
                  fib_to_myo1 = (G_str1GAP/FIB_CAP) * (U[V] - U[fib_volt]); 
                 
                  myo_to_fib_long1 = (G_wk1GAP/MYO_CAP) * (U1[fib_volt] - U[V]);
                  fib_to_myo_long1 = (G_wk1GAP/FIB_CAP) * (U[V]-U1[fib_volt]); 
                 
                  dV2 = calculate_TP06_current_OpSplit1( U1, dtshort, lookup, celltype[1], stimCurrent2 );
                  dVf2 = calculate_MacFib_current1( U1, dtshort);  
                  
                  myo_to_fib2 = (G_str2GAP/MYO_CAP) * (U1[fib_volt] - U1[V]);
                  fib_to_myo2 = (G_str2GAP/FIB_CAP) * (U1[V] - U1[fib_volt]); 
                 
                  myo_to_fib_long2 = (G_wk2GAP/MYO_CAP) * (U[fib_volt]-U1[V]); 
                  fib_to_myo_long2 = (G_wk2GAP/FIB_CAP) * (U1[V]-U[fib_volt]);
                    

                       
                   
 	    	 U[V] = U[V] + dtshort * (-dV1 + NUMFIB*myo_to_fib1  + NUMFIB*myo_to_fib_long1); //4 fibroblasts in each fibroblast unit 
                 U[fib_volt] = U[fib_volt] + dtshort * (-dVf1 + fib_to_myo1 + fib_to_myo_long2);
                 U1[V] = U1[V] + dtshort * (-dV2 + NUMFIB*myo_to_fib2  + NUMFIB*myo_to_fib_long2); //4 fibroblasts in each fibroblast unit 
                 U1[fib_volt] = U1[fib_volt] + dtshort * (-dVf2 + fib_to_myo2 + fib_to_myo_long1);  
	    }

         
             
		  /* update state u with new values stored in U */
		  for (m = 1; m <= num_states; m++)
		       {
	   	 	  u[1][m] = U[m];u[2][m]=U1[m];
	   	 	}
/* end of step 2 */
      
      
/* output stf file every 1 ms */
	 
	  if (t % 1000 == 0)
          {
              printf("%5.10f\t%5.10f\t%5.10f\t%5.10f\t%5.10f\t%5.10f\n",u[1][1],u[1][21],u[2][1],u[2][21],u[1][18],u[2][18]);
          }
 
         
/* write to checkpoint file if needed */
      if (CHKPT_WRITE)
        {
    	  if (t == CHKPT_WRITE_TIME)
    	    {
    		dummy = checkpoint_write( u, time, t, CHKPT_WRITE_TIME, N );
    		exit(1);
    	    }
        }

  }
 // printf("leaving main loop\n");

  /* save upstroke and downstroke data to files */
  /*sprintf( outputFile,"%supStrokeTimeS1.stf",OUTPUTFILEROOT);
  writeData(outputFile, upStrokeTimeS1, geom, nrows, ncols);  
  sprintf( outputFile,"%sdownStrokeTimeS1.stf",OUTPUTFILEROOT);
  writeData(outputFile, downStrokeTimeS1, geom, nrows, ncols);  
  sprintf( outputFile,"%supStrokeTimeS2.stf",OUTPUTFILEROOT);
  writeData(outputFile, upStrokeTimeS2, geom, nrows, ncols);  
  sprintf( outputFile,"%sdownStrokeTimeS2.stf",OUTPUTFILEROOT);
  writeData(outputFile, downStrokeTimeS2, geom, nrows, ncols);
  sprintf( outputFile,"%sthreshold.stf",OUTPUTFILEROOT);
  writeData(outputFile, threshold, geom, nrows, ncols);  
  sprintf( outputFile,"%sdiffusion.stf",OUTPUTFILEROOT);
  writeData(outputFile, D, geom, nrows, ncols);
  */              
  fclose(egPtr);
    
  /* end of main loop */

  /* Free memory */
  free_fmatrix(lookup,0,num_lookup,0,voltage_steps);
  free_fmatrix(u,1,N,1,num_states);
  free_imatrix(geom, 1, nrows, 1, ncols);
  free_imatrix(nneighb, 1, N, 1, 8);
  free_fvector(dVdt, 1, N );
  free_fvector(new_Vm, 1, N );
  free_fvector(old_Vm, 1, N );
  free_fvector(U, 1, N);
  free_fvector(U1, 1, N);
  free_fvector(params, 1, num_params);
  free_ivector(celltype, 1, N);
  free_fvector(stimCurrent,1,N);
              
  free_fvector(maxV, 1, N);
  free_fvector(minV, 1, N);
  free_fvector(threshold, 1, N);
  free_fvector(upStrokeTimeS1, 1, N);
  free_fvector(downStrokeTimeS1, 1, N);
  free_fvector(upStrokeTimeS2, 1, N);
  free_fvector(downStrokeTimeS2, 1, N);
    
  //fclose(egptr);

}

void writeData(char *fname, double *dataToWrite, int **geomarray, int nrows, int ncols)
{
    int lay, row, col;
    int index, outint;
    double outdouble;
    FILE *stf_file;
    
    stf_file = fopen( fname, "w" );
    
    //fprintf(stf_file, "NAME Vm\n");
   // fprintf(stf_file, "RANK 2\n");
   // fprintf(stf_file, "DIMENSIONS %d %d\n",ncols, nrows);
   // fprintf(stf_file, "BOUNDS %d %d %d %d\n", 0, ncols-1, 0, nrows-1);
    //fprintf(stf_file, "SCALAR\n");
    //fprintf(stf_file, "DATA\n");

    for (row = 1;row <= nrows; row+=1)
        {
        for (col = 1; col <= ncols; col+=1)
            {
            index = (geomarray[row][col] > 0)?geomarray[row][col]:0;
            if (index > 0)
                outdouble = dataToWrite[index];
            else
                outdouble = 0.0;
            fprintf(stf_file, "%4.2f ", outdouble);
            }
        fprintf( stf_file, "\n");
        }
    fclose(stf_file);
}

