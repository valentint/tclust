#ifdef MATLAB_MEX_FILE

#include "ML_package.h"
#include "tclust.h"
#include "tkmeans.h"

R_MEAL_SETTINGS ("heinrich_fritz@hotmail.com") ;	//	settings for the R meal - implementation

	void mexFunction(int nlhs, mxArray* plhs[], int nrhs, mxArray *prhs[])
	{
		if (nrhs < 1)
		{
			meal_error ("At least one argument expected.") ;
			return ;
		}

		mxArray *pSwitch = prhs [0] ;

		if (mxGetM (pSwitch) != 1 || 
			mxGetN (pSwitch) != 1)
		{
			meal_error ("First argument has to be of length 1.") ;
			return ;
		}

		int nSwitch = (int) *mxGetPr (pSwitch) ;


		switch (nSwitch)
		{
		case 0:
			ML_tkmeans (nlhs, plhs, nrhs - 1, plhs + 1) ;
			break ;
		}
	}

/////////////////////////////////////
//	exporting functions to Matlab  //
/////////////////////////////////////

	void ML_tkmeans (int nlhs, mxArray* plhs[], int nrhs, mxArray *prhs[])
		//	(int *pnParIn, int *pnParOut, double *pdParIn, double *pdParOut, double *pdX, double *pdM, int *pnAssign, double *pdClustSize, double *pdWeights, double *pdObjER, int *pnConvER)
	{
		if (nrhs != 4)
			meal_error ("4 input arguments expected") ;
		if (nlhs != 4)
			meal_error ("4 output arguments expected") ;

			//	pnParIn
		mxArray *pmxnParIn = prhs[0] ;

		if (mxGetM (pmxnParIn) * mxGetN (pmxnParIn) != 5)
			meal_error ("Length of vector pnParIn is expected to be 5.") ;

		double *pnParIn = mxGetPr (pmxnParIn) ;

		int k = (int) pnParIn [0], nIter = (int) pnParIn [3] ;

			//	pdParIn
		mxArray *pmxdParIn = prhs[1] ;

		if (mxGetM (pmxdParIn) * mxGetN (pmxdParIn) != 1)
			meal_error ("Length of vector pdParIn is expected to be 1.") ;
		double *pdParIn = mxGetPr (pmxdParIn) ;

			//	pdX
		mxArray *pmxdX = prhs[2] ;

		int n = mxGetM (pmxdParIn), p = mxGetN (pmxdParIn) ;

		double *pdX = mxGetPr (pmxdX) ;

			//	pnParOut
		int pnParOut [4] ;

			//	pdParOut
		mxArray *pmxdParOut = plhs[1] = mxCreateDoubleMatrix (1, 1, mxREAL) ;
		double *pdParOut = mxGetPr (pmxdParOut) ;

			//	pdM
		mxArray *pmxdM = plhs[2] = mxCreateDoubleMatrix (p, k, mxREAL) ;
		double *pdM = mxGetPr (pmxdM) ;

			//	pnAssign
		SVecN vnAssign (n) ;
		int *pnAssign = vnAssign ;

			//	pdClustSize
		mxArray *pmxdClustSize = plhs[4] = mxCreateDoubleMatrix (k, 1, mxREAL) ;
		double *pdClustSize = mxGetPr (pmxdClustSize) ;

			//	pdWeights
		mxArray *pmxdWeights = plhs[5] = mxCreateDoubleMatrix (k, 1, mxREAL) ;
		double *pdWeights = mxGetPr (pmxdWeights) ;

			//	pdObjER
		mxArray *pmxdObjER = plhs[6] = mxCreateDoubleMatrix (nIter, 1, mxREAL) ;
		double *pdObjER = mxGetPr (pmxdObjER) ;

			//	pnConvER
		SVecN vnConvER (nIter) ;
		int *pnConvER = vnConvER ;

		mxArray *pmxdX = prhs[2] ;

		int n = mxGetM (pmxdParIn), p = mxGetN (pmxdParIn) ;

		double *pdX = mxGetPr (pmxdX) ;

		TRY (
				CTKMeans (	  n
							, p
							, k
							, pdParIn[0]	//	dAlpha
							, pdParIn[1]	//	dZeroTol
							, pdX
							, pnAssign
							, pdClustSize
							, pdWeights
							, pnParIn[1]	//	nEqualWeights
							, pnParIn[2]	//	nTrace
							, pdM
				).SetPtr (    pnParOut + 0	//	pnConvCount
				            , pnParOut + 1	//	pnIterSuccess
				            , pnParOut + 2	//	pnCode
				            , pnParOut + 3	//	pnErrExc
				            , pdParOut + 0	//	pdBestObj
				).calc (	  nIter			//	nIter
							, pnParIn[4]	//	nKSteps
							, pnConvER
							, pdObjER
				) ;
			)

			//	pnParOut
		mxArray *pmxnParOut = plhs[0] = mxCreateDoubleMatrix (4, 1, mxREAL) ;
		Int2Double (mxGetPr (pmxnParOut), pnParOut, 4) ;

			//	pnAssign
		mxArray *pmxnAssign = plhs[3] = mxCreateDoubleMatrix (n, 1, mxREAL) ;
		Int2Double (mxGetPr (pmxnAssign), pnAssign, n) ;

			//	pnConvER
		mxArray *pmxnConvER = plhs[7] = mxCreateDoubleMatrix (nIter, 1, mxREAL) ;
		Int2Double (mxGetPr (pmxnConvER), pnConvER, nIter) ;
	}

#endif	//	#ifdef MATLAB_MEX_FILE
