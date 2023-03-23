#include "R_package.h"
#include "tclust.h"
#include "tkmeans.h"

#ifdef ES_DEV_ENV
	#include "..\..\..\RDev\R_meal.h"
#else
	#include "R_meal.h"
#endif	//	#ifdef ES_DEV_ENV

R_MEAL_SETTINGS ("heinrich_fritz@hotmail.com") ;	//	settings for the R meal - implementation

////////////////////////////////
//	exporting functions to R  //
////////////////////////////////

	void C_tkmeans (int *pnParIn, int *pnParOut, double *pdParIn, double *pdParOut, double *pdX, double *pdM, int *pnAssign, double *pdClustSize, double *pdWeights, double *pdObjER, int *pnConvER)
	{
		TRY (
				CTKMeans (	  pnParIn[0]	//	n
							, pnParIn[1]	//	p
							, pnParIn[2]	//	k
							, pdParIn[0]	//	dAlpha
							, pdParIn[1]	//	dZeroTol
							, pdX
							, pnAssign
							, pdClustSize
							, pdWeights
							, pnParIn[3]	//	nEqualWeights
							, pnParIn[4]	//	nTrace
							, pdM
				).SetPtr (    pnParOut + 0	//	pnConvCount
				            , pnParOut + 1	//	pnIterSuccess
				            , pnParOut + 2	//	pnCode
				            , pnParOut + 3	//	pnErrExc
				            , pdParOut + 0	//	pdBestObj
				).calc (	  pnParIn[5]	//	nIter
							, pnParIn[6]	//	nKSteps
							, pnConvER
							, pdObjER
				) ;
			)
	}

	void C_tclust (int *pnParIn, int *pnParOut, double *pdParIn, double *pdParOut, double *pdX, double *pdM, double *pdS, int *pnAssign, double *pdClustSize, double *pdWeights, double *pdZ, double *pdObjER, int *pnConvER)
	{
		TRY (
                //VT::08.05.2018
                //VT::22.03.2023
                // meal_printf("\nMY-TRACE ... Entering C_tclust()/CTClust()...\n");
                // meal_printf("\n%d %d %d \n", pnParIn[0], pnParIn[1], pnParIn[2]);
				CTClust (pnParIn, pnParOut, pdParIn, pdParOut, pdX, pdM, pdS, pnAssign, pdClustSize, pdWeights, pdZ, pdObjER, pnConvER) ;
			)
	}

	void OptVectors_C (int *pdwParamIn, double *pdU, double *pdD, double *pdS, double *pdCSize)
	{
		double adParams [] = {0, 0, 1e-16} ;

		t_size	p = pdwParamIn[0],
				K = pdwParamIn[1] ;

		SMatD mU (pdU, p, p) ;
		SMatD mDd (pdD, p, K) ;

		TRY (
			CTClust obj (pdwParamIn, adParams, pdS, pdCSize) ;
			obj.OptVectors (mU, mDd) ;
		)

	}

	void restr_dir_C (int *pdwParamIn, int *pdwParamOut, double *pdParamIn, double *pdS, double *pdCSize, double *pdEVec, double *pdEval)
	{
		TRY (
			CTClust obj (pdwParamIn, pdParamIn, pdS, pdCSize, pdEVec, pdEval) ;
			*pdwParamOut = obj.restr_dir () ;
			)

	}
	
	void restr_prop_C (int *pdwParamIn, int *pdwParamOut, double *pdParamIn, double *pdS, double *pdCSize, double *pdEVec, double *pdEval)
	{
		TRY (
			CTClust obj (pdwParamIn, pdParamIn, pdS, pdCSize, pdEVec, pdEval) ;
			*pdwParamOut = obj.restr_prop () ;
		)
	}


	void sme_eigen_sym_2x2_norm (double * const pdEval, double *const pdEVec, const double *const pd, const double *pdZeroTol)
	{
		sme_eigen_sym_2x2_norm_raw_NC (pdEval, pdEVec, pd, *pdZeroTol) ;
	}


	void RestrictEigenValues (int *pnParIn, int *pnParOut, double *pdParIn, double *pdParOut, double *pdEv, double *pdClustSize)
	{
		t_size p = pnParIn [0], k = pnParIn [1] ;
		TRY (
			pnParOut [0] = RestrictEigenValues (!SMatD (pdEv, p, k), SVecD (pdClustSize, k), pdParIn[0], pdParIn[1], pdParOut[1]) ;
		)
	}

	void RestrictEigenValues_deter (int *pnParIn, int *pnParOut, double *pdParIn, double *pdParOut, double *pdEv, double *pdClustSize)
	{
		t_size p = pnParIn [0], k = pnParIn [1] ;
		TRY (
			pnParOut [0] = RestrictEigenValues_deter (!SMatD (pdEv, p, k), SVecD (pdClustSize, k), pdParIn[0], pdParIn[1], pdParOut[1]) ;
		)
	}
