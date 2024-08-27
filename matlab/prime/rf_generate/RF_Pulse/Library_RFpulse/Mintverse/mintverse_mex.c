
#define MEX_COMPILE

/*
#define DEBUG
*/

#include <mex.h>
#include "mintverse.c"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])

	/*  function [bv,gv] = mintverse_mex(b,g,dt,bmax,gmax,smax,dtout,emax)*/

 
{
double *br;	/* 	real part of B1 pulse (arbitrary units)*/
double *bi;	/* 	real part of B1 pulse (arbitrary units)*/
double *g;	/*	Gradient pulse (arbitrary units) */
double *dt;	/*	Time increments (arbitrary units) */
long n;		/* 	Number of points. */
double bmax;	/*	Maximum b1-field amplitude (units of *b) */
double gmax;	/* 	maximum gradient amplitude, (units of *g) */
double smax;	/*	Maximum slew rate, (units of *g per unit of *dt) */
double tout;	/*	Desired Output sample rate */
double emax;	/* 	Max energy (optional). 	*/

double *brwork;	/* 	Output B1 pulse (arbitrary units)*/
double *biwork;	/* 	Output B1 pulse (arbitrary units)*/
double *gwork;	/*	Output Gradient pulse (arbitrary units) */

long nout;	/* 	Number of output points. 	*/
double *brout;	/*	Matlab Output B1 pulse (real)		*/
double *biout;	/*	Matlab Output B1 pulse (imag)		*/
double *gout;	/*	Matlab Output gradient pulse.		*/

long count;	/* 	For loops.	*/
int bcomplex;	/*	1 if b is complex */

/*	Check Arguments 	*/

bcomplex = 0;

if (nrhs < 3)
	mexErrMsgTxt("Incorrect Input Argument List");

if (nlhs < 2 )
	mexErrMsgTxt("Function should have at least 2 outputs.");


/*================	Get required input parameters================	*/

br = mxGetPr(prhs[0]);
n = (long) (mxGetM(prhs[0]) * mxGetN(prhs[0]));
if (mxIsComplex(prhs[0]))
	{
	bi = mxGetPi(prhs[0]);
	bcomplex = 1;
	}
else
	{
	bi = (double *) malloc(n*sizeof(double));
	for (count=0; count <n; count++)
		bi[count]=0.0;
	bcomplex = 0;
	}

g = mxGetPr(prhs[1]);
if ( (long) (mxGetM(prhs[1])* mxGetN(prhs[1]))  != n)
	{
	printf("Error:  B1 and Gradient should be the same length. \n");
	}

if ( (mxGetM(prhs[2])* mxGetN(prhs[2]))  > 1)	/* dt array passed */
	dt= mxGetPr(prhs[2]);
else							
		/* single dt value - allocate array and set to value */
	{
	dt = (double *) malloc(n*sizeof(double));
	*dt= *mxGetPr(prhs[2]);
	for (count=0; count<n; count++)
		dt[count]=dt[0];
	}



/*================ Get optional input parameters  ===============	*/

/* 	If a default is used, a message is displayed. 	*/

if (nrhs > 3)
	bmax = *mxGetPr(prhs[3]);
else
	{
	bmax = 0.14;	/* Default maximum in G */
	printf("No max B1 amplitude given - Default is %g \n",bmax);
	}

if (nrhs > 4)
	gmax = *mxGetPr(prhs[4]);
else
	{
	gmax = 3.9;	/* Default maximum in G/cm */
	printf("No max gradient amplitude given - Default is %g \n",gmax);
	}

if (nrhs > 5)
	smax = *mxGetPr(prhs[5]);
else
	{
	smax = 14500;	/* Default maximum in G/cm/s */
	printf("No max slew rate given - Default is %g \n",smax);
	}

if (nrhs > 6)
	tout = *mxGetPr(prhs[6]);
else
	{
	tout = *dt;	/* Default output rate. */
	printf("No output sample step given - Default is %g \n",tout);
	}

if (nrhs > 7)
	emax = *mxGetPr(prhs[7]);
else
	{
	emax = -1.0;	/* Maximum pulse energy. */
	printf("No max B1 energy constraint given.  \n",tout);
	}


/*	Call C-function here */

mintverse(br,bi,g,dt,n,bmax,gmax,smax,tout,emax,&nout,&brwork,&biwork, &gwork);


/*	Allocate Matlab outputs, and copy function outputs to them. */

if (bcomplex==0)
	plhs[0] = mxCreateDoubleMatrix(nout,1,mxREAL);
else
	plhs[0] = mxCreateDoubleMatrix(nout,1,mxCOMPLEX);

plhs[1] = mxCreateDoubleMatrix(nout,1,mxREAL);
brout = mxGetPr(plhs[0]);
if (bcomplex==1)
	biout = mxGetPi(plhs[0]);
gout = mxGetPr(plhs[1]);

#ifdef DEBUG
	printf("mex function freeing memory. \n");
#endif

bcopy(brwork,brout,nout*sizeof(double));
if (bcomplex==1)
	bcopy(biwork,biout,nout*sizeof(double));
bcopy(gwork,gout,nout*sizeof(double));

/*	Free up space that was allocated here and in mintverse(). */

free(brwork);
free(biwork);
free(gwork);

if (  (mxGetM(prhs[2])* mxGetN(prhs[2]))  == 1)
	free(dt);

}




