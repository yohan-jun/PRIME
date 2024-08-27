
#define MEX_COMPILE

/*
#define DEBUG
*/

#include <mex.h>
#include "mintverse.c"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])

	/*  function [gv] = verse_mex(b,g,gv) */
 
{
double *br;	/* 	B1 pulse (arbitrary units)*/
double *bi;	/* 	B1 pulse (arbitrary units)*/
double *g;	/*	Gradient pulse (arbitrary units) */
long n;		/* 	Number of input points. */
long nv;	/* 	Number of VERSE points. */
double *bvr;	/* 	VERSE B1 pulse (arbitrary units)*/
double *bvi;	/* 	VERSE B1 pulse (arbitrary units)*/
double *gv;	/*	VERSE Gradient pulse (arbitrary units) */

int galloc = 0;		/* 1 if g is allocated here. */
long count;	/* 	For loops.	*/
int bcomplex=0;	/* 	1 if B is complex. */

/*	Check Arguments 	*/

if (nrhs < 3)
	mexErrMsgTxt("Incorrect Input Argument List");

if (nlhs < 1 )
	mexErrMsgTxt("Function should have at least 1 output.");


/*================	Get required input parameters================	*/


br = mxGetPr(prhs[0]);
n = (long) (mxGetM(prhs[0])* mxGetN(prhs[0]));
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

printf("verse_mex - %d points \n",n);

if ( (long) (mxGetM(prhs[1])* mxGetN(prhs[1]))  == 1)	/* g level given. */
	{
	printf("verse_mex - gradient level %g. \n",*g);
	g = (double *) malloc(n * sizeof(double));
	galloc = 1;			/* indicate g is allocated. */
	g[0] = *mxGetPr(prhs[1]);			/* 1st point. */
	for (count=1; count < n; count++)
		g[count] = g[0];			/* Set to 1st value */
	}
else if ( (long) (mxGetM(prhs[1])* mxGetN(prhs[1]))  != n)
	{
	printf("Error:  B1 and Gradient should be the same length. \n");
	}

gv = mxGetPr(prhs[2]);
nv = (long) (mxGetM(prhs[2])* mxGetN(prhs[2]));
printf("verse_mex - %d output points. \n",nv);


/*	Allocate Matlab outputs, and copy function outputs to them. */

if (bcomplex==0)
	{
	plhs[0] = mxCreateDoubleMatrix(nv,1,mxREAL);
	bvr = mxGetPr(plhs[0]);
        bvi = (double *) malloc(nv*sizeof(double));
	}
else
	{
	plhs[0] = mxCreateDoubleMatrix(nv,1,mxCOMPLEX);
	bvr = mxGetPr(plhs[0]);
	bvi = mxGetPi(plhs[0]);
	}


/*	Call C-function here */

verse(br,bi,g,n,gv,bvr,bvi,nv);

	
if (galloc == 1)
	free(g);

}




