#ifdef MATLAB_MEX_FILE

#include "ML_passrng.h"

#ifdef ES_DEV_ENV
	#include "../SMat/smat.meal.h"
	#include "../SMat/smat_meal_passrng.h"
	#include "../SMat/smat_meal_passrng_hpp.h"
#else
	#include "smat.meal.h"
	#include "smat_meal_passrng.h"
	#include "smat_meal_passrng_hpp.h"
#endif	//	#ifdef ES_DEV_ENV

#include "matrix.h"

	void ML_pass_runif (int nlhs, mxArray* plhs[], int nrhs, mxArray *prhs[])
	{
		mxArray *aData = prhs [0] ;

		int n = mxGetM (aData) * mxGetN (aData) ;
		double *pData = mxGetPr (aData) ;

		pass_runif (pData, n) ;
	}

	void ML_pass_rnorm (int nlhs, mxArray* plhs[], int nrhs, mxArray *prhs[])
	{
		mxArray *aData = prhs [0] ;

		int n = mxGetM (aData) * mxGetN (aData) ;
		double *pData = mxGetPr (aData) ;

		pass_rnorm (pData, n) ;
	}

#endif	//	#ifdef MATLAB_MEX_FILE
