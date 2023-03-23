#ifdef ES_DEV_ENV
	#include "../../../RDev/R.Inc.h"
	#include "../../../SMat/smat.def.h"
#else
	#include "R.Inc.h"
	#include "smat.def.h"
#endif

	EXPORT void C_tkmeans (int *pnParIn, int *pnParOut, double *pdParIn, double *pdParOut, double *pdX, double *pdM, int *pnAssign, double *pdClustSize, double *pdWeights, double *pdObjER, int *pnConvER) ;
	EXPORT void C_tclust (int *pnParIn, int *pnParOut, double *pdParIn, double *pdParOut, double *pdX, double *pdM, double *pdS, int *pnAssign, double *pdClustSize, double *pdWeights, double *pdZ, double *pdObjER, int *pnConvER) ;

	EXPORT void restr_dir_C (int *pdwParamIn, int *pdwParamOut, double *pdParamIn, double *pdS, double *pdCSize, double *pdEVec, double *pdEval) ;
	EXPORT void restr_prop_C (int *pdwParamIn, int *pdwParamOut, double *pdParamIn, double *pdS, double *pdCSize, double *pdEVec, double *pdEval) ;
	EXPORT void RestrictEigenValues_deter	(int *pnParIn, int *pnParOut, double *pdParIn, double *pdParOut, double *pdEv, double *pdClustSize) ;
	EXPORT void RestrictEigenValues			(int *pnParIn, int *pnParOut, double *pdParIn, double *pdParOut, double *pdEv, double *pdClustSize) ;

	EXPORT void OptVectors_C (int *pdwParamIn, double *pdU, double *pdD, double *pdS, double *pdCSize) ;

	EXPORT void sme_eigen_sym_2x2_norm (double * const pdEval, double *const pdEVec, const double *const pd, const double *pdZeroTol) ;
	
