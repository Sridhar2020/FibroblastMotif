/* Header file for TP06 in 2D with isotropic diffusion */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define ROWS            2      //340
#define COLUMNS         1
#define NUM_STATES      23      /* number of states in u array */
#define NUM_ITERATIONS  10000000  /* maximum duration of simulation (60000 is 6s) */
//#define NUM_POINTS      100   /* Length of cable */
#define DT              0.001    /* maximum timestep (ms) */
#define DX              0.25    /* 0.25 mm : Assume dx = dy */
#define DIFFUSION       0     /* mm2/ms :isotropic diffusion */
#define VMOFFSET        1001.0  /* offset for lookup table */
#define GAIN            10.0    /* gain for lookup table */
#define NUM_LOOKUP      40      /* number of variables in lookup table*/
#define NUM_PARAMS      50      /* number of cell model parameters */
#define VOLTAGE_STEPS   2001    /* number of voltage steps in lookup table */
#define CM              1.0     /* Membrane capacitance (uF/cm2) */

#define BCL1            300
//#define BCL2            350.0  
#define G_STRONG1       0.0
//#define G_STRONG2        0.0 
#define G_WEAK1         1.5
//#define G_WEAK2         0.0
#define START            0
#define DELAY            0

/* Macros for numerical recipes routines */
#define FREE_ARG       	char*
#define NR_END         	1

/* diffusion file */
#define DIFFUSIONFILE   "DiffusionCoefficient.txt"

/* Macros for output */
#define OUTPUTFILEROOT  "TP06_2D_"
#define STFFILEROOT     "STFfiles/TP06_2D_"

/* checkpointing */
/* read/write times should be multiples  f 5 ms */
#define CHKPTROOT 		    "checkpoint"
#define CHKPT_READ	        0
#define CHKPT_READ_TIME     200000 // 4000 ms
#define CHKPT_WRITE	        0
#define CHKPT_WRITE_TIME    200000 // 4000 ms

/* forward declaration of all functions used */

/* PDE solver */
int initialise_geometry_2D( int **geom, int nrows, int ncols, int **nneighb );
void initialise_variables_2D( double **u, int N );
void initialise_spiral_2D( double **u, int N, int ny, int nx );
void initialise_diffusion_2D( double *D, int nrows, int ncols );
    
int create_TP06_lookup_OpSplit_2D( double **lookup );
double calculate_TP06_current_OpSplit( double *U, double dt, double **lookup, int celltype, double stimCurrent );
double calculate_TP06_current_OpSplit1( double *U, double dt, double **lookup, int celltype, double stimCurrent );
/*fibroblast current*/
double calculate_MacFib_current( double *U, double dt); 
double calculate_MacFib_current1( double *U, double dt); 

/*fibroblast coupling */
double conjugate_coupling_fib2myo( double **u, int n, int **nneighb );
double conjugate_coupling_myo2fib( double **u, int n, int **nneighb );
  
double diffusion_2D( double **u, int **nneighb, int n, int N, double D, double dx2 );
int free_arrays( double **u, int N, int num_parameters, double **lookup, int voltage_steps, int num_lookup );

/* checkpointing */
int checkpoint_write( double **u, double time, int t, int count, int N );
int checkpoint_read( double **u, double *time, int *t, int count, int N);

/* Numerical recipes routines */
double *fvector( long nl, long nh );
int *ivector( long nl, long nh );
int **imatrix( long nrl, long nrh, long ncl, long nch );
double **fmatrix( long nrl, long nrh, long ncl, long nch );
int ***i3dmatrix( long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_ivector( int *m, long nl, long nh );
void free_imatrix( int **m, long nrl, long nrh, long ncl, long nch );
void free_fvector( double *m, long nl, long nh );
void free_fmatrix( double **m, long nrl, long nrh, long ncl, long nch );
void free_i3dmatrix( int ***m, long nrl, long nrh, long ncl, long nch, long ndl, long ndh );
void nrerror( char error_text[] );

int stfout_2D( double **u, int **geomarray, int stfcount, int nx, int ny );

/* RGB file output */
//int write_rgb( int nx, int ny, int N, double **u, int **geom, char *fname, int I, double iMax, double iMin );
//void putbyte( FILE *outfile, unsigned char value );
//void putshort( FILE *outfile, unsigned short value );
//static int putlong( FILE *outfile, unsigned long value );

/*** END ***/
