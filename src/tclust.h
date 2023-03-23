

/*
#ifdef ES_DEV_ENV
	#include "../../../SMat/smat.h"
#else
	#include "smat.h"
#endif
//	or 
*/
#include "clust.h"

	class CTClust
	{
	public:

		CTClust (int *pdwParamIn, double *pdParamIn, double *pdS, double *pdCSize) ;
		CTClust (int *pdwParamIn, double *pdParamIn, double *pdS, double *pdCSize, double *pdEVec, double *pdEval) ;
		CTClust (int *pnParIn, int *pnParOut, double *pdParIn, double *pdParOut, double *pdX, double *pdM, double *pdS, int *pnAssign, double *pdClustSize, double *pdWeights, double *pdZ, double *pdObjER, int *pnConvER) ;
		~CTClust () ;
		void calc () ;

		BOOL restrEval () ;

		BOOL restr_diff_ax () ;
		BOOL restr_equal () ;
		BOOL restr_dir () ;
		BOOL restr_prop () ;
		BOOL restr_none () ;

		void OptVectors (SMatD &mU, const SCMatD &mIDiagD) ;

	protected:

//	Member Functions  //

		void SetAllCovmatsIdent () ;

		void FindInitClustAssignment () ;
		void FindInitClustSize () ;
		void FindInitClustSize_R () ;

		void CalcDensity (const SCMatD &mX, const SVecD &vDens, const SVecD &vCurM, const SCVecD &vEVal, const SCMatD &mEVec, const double dFact = 1) ;
		void FindNearestClust (const SVecD &vDisc, const SVecN &vInd) ;
		void FindNearestClust_old (const SVecD &vDisc, const SVecN &vInd) ;
		void FindNearestClust_new (const SVecD &vDisc, const SVecN &vInd) ;	
		void select_cluster (double &dDisc, int &nInd, const SCVecD &row) ;
		void select_cluster_old (double &dDisc, int &nInd, const SCVecD &row) ;
		void select_cluster_new (double &dDisc, int &nInd, const SCVecD &row) ;
		BOOL FindClustAssignment () ;
		BOOL FindClustAssignment_cat () ;
		BOOL FindClustAssignment_fuzzy () ;													//	temp5
		void SetCatZ (const SVecD &z, int nIdx) ;
		void CalcFuzzyRow (const SCVecD &ll, const SVecD &z, double &dDisc, int &nInd) ;	//	temp1

		BOOL CheckParams () ;
		void CheckRestrOk () ;

		void EstimClustParams () ;
		void EstimClustParams_cat () ;
		void EstimClustParams_fuzzy () ;
		double CalcObjFunc () ;
		double CalcObjFunc_cat () ;
		double CalcObjFunc_fuzzy () ;
		void LoadCluster (SMatD &c, t_size k) ;

		void CheckRestrictions () ;

		void CalcClusterSize_cat () ;

		void FindOutliers (const SVecD &vDisc, const SVecN &vInd) ;		//	temp2
		void FindOutliers_old (const SVecD &vDisc, const SVecN &vInd) ;	
		void FindOutliers_new (const SVecD &vDisc, const SVecN &vInd) ;	//	temp2

		void SaveCurResult (double dCurObj, int nCode = 0) ;

//	Member Variables  //

		t_size m_n, m_p, m_K ;																		//	int input parameters
		int	m_nFuzzy, m_nIter, m_nKSteps, m_nEqualWeights, m_nRestr, m_nDeter, m_nTrace ;
		t_size m_dwIterTune1, m_dwIterTune2, m_dwIterTune3 ;
		t_size m_dwOVV ;

		int &m_nConvCount, &m_nIterSuccess, &m_nCode, &m_nErrExc, &m_dwCountRestrOk ;				//	int output parameters
		int * const m_pnConvER ;
		const double m_dAlpha, m_dRestrFactor, m_dM, m_dZeroTol ; 									//	double input parameters

		double &m_dBestObj, &m_dUnRestrFactBest ;													//	double output parameters
		double * const m_pdObjER ;

																									//	some constants
		const double m_dDensFact, m_dPInv, m_dRestrFactp1p, m_dZeroTolSqrt, m_dMm1Inv ;
		double m_dUnRestrFact ;
		const t_size m_dwNoTrim, m_dwTrim ;


		SVecN m_vInd, m_vIndBest, m_vIndOld, m_vCurInd, m_vRank ;
		SVecD m_vWeights, m_vBestWeights, m_vClustSize, m_vClustSizeBest, m_vDisc ;//, m_vDiscSorted ;

		SMatD m_mCurM, m_mBestM, m_mX, m_mLL, m_mEVal, m_mZ, m_mZ_best, m_mZOld ;
		SCMatArrayD m_amEVec, m_amCurS, m_amBestEVec, m_amBestS ;

		SDataRef_Static m_aTemp [15] ;																//	2do: replace with tempRef ()
	} ;

	BOOL RestrictEigenValues_deter (const SVMatD &mEV, const SCVecD & vClustSize, double dFact, double dZeroTol, double &dUnRestrFact) ;
	BOOL RestrictEigenValues (const SVMatD &mEV, const SCVecD & vClustSize, double dFact, double dZeroTol, double &dUnRestrFact) ;
