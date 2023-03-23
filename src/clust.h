#ifndef INC_CLUST_H
#define INC_CLUST_H


#ifdef ES_DEV_ENV
	#include "../../../SMat/smat.h"
#else
	#include "smat.h"
#endif


	class CClust
	{
	protected:

		CClust () ;
	public:

		CClust (t_size n, t_size p, t_size k, double dAlpha, double dZeroTol, double *pdX, int *pnAssign, double *pdClustSize, double *pdClustWeights, int nEqualWeights = 0, int nTrace = 0) ;
		void calc (int nIter, int nKSteps, int * const pnConvER, double * const pdObjER) ;

		CClust &SetPtr (int *pnConvCount = NULL, int *pnIterSuccess = NULL, int *pnCode = NULL, int *pnErrExc = NULL, double *pdBestObj = NULL) ;

		~CClust () ;


	protected:

//	Member Functions  //

		virtual double CalcObjFunc () = 0 ;
		virtual void EstimClustParams () = 0 ;
		virtual void EstimInitClustParams (int k, const SCVecN &vNIdx) = 0 ;
		virtual void SaveCurResult () {} ;
		virtual BOOL RestrictClustParam () { return TRUE ; }
		virtual void SetSingularIniParams () { }
		virtual void CalcDensity (const SCMatD &mX, const SVecD &vDens, t_size k, const double dFact = 1) = 0 ;
		virtual BOOL FindClustAssignment () = 0 ;
		virtual void FindOutliers (const SVecD &vDisc, const SVecN &vInd) = 0 ;
		virtual void FindNearestClust (const SVecD &vDisc, const SVecN &vInd) = 0 ;
		virtual void select_cluster (double &dDisc, int &nInd, const SCVecD &row) = 0 ;


		void SaveCurResult0 (double dCurObj, int nCode) ;

		virtual void FindInitClustAssignment () ;
		void FindInitClustSize_R () ;
		void FindInitClustSize () ;


		void SetCatZ (const SVecD &z, int nIdx) ;

		BOOL CheckParams () ;
		void CheckRestrOk () ;

		void LoadCluster (SMatD &c, t_size k) ;

//	Member Variables  //

		t_size m_n, m_p, m_K ;																		//	int input parameters
		int	m_nEqualWeights, m_nTrace ;

		int *m_pnConvCount, *m_pnIterSuccess, *m_pnCode, *m_pnErrExc ;								//	int output parameters

		const double m_dAlpha, m_dZeroTol ;		 													//	double input parameters

		double m_dBestObj ;
		double *m_pdBestObj ;																		//	double output parameters
																									//	some constants

		const double m_dDensFact, m_dPInv, m_dZeroTolSqrt ;
		const t_size m_dwNoTrim, m_dwTrim ;

		SVecN m_vInd, m_vIndBest, m_vIndOld, m_vCurInd, m_vRank ;
		SVecD m_vWeights, m_vBestWeights, m_vClustSize, m_vClustSizeBest, m_vDisc ;

		SMatD m_mCurM, m_mBestM, m_mX, m_mLL/*, m_mEVal, m_mZ, m_mZ_best, m_mZOld*/ ;

		SDataRef_Static m_aTemp [15] ;																//	2do: replace with tempRef ()
	} ;

	class CClust_N : virtual public CClust
	{
	public:

	protected:

		virtual void FindOutliers (const SVecD &vDisc, const SVecN &vInd) ;
		virtual void FindNearestClust (const SVecD &vDisc, const SVecN &vInd) ;
		virtual void select_cluster (double &dDisc, int &nInd, const SCVecD &row) ;

	} ;

	class CClust_M : virtual public CClust
	{
	public:
		CClust_M (t_size p, t_size k, double *pdM) ;

	protected:

		virtual void EstimInitClustParams (int k, const SCVecN &vNIdx) ;
		virtual void SaveCurResult () ;

		virtual void CalcDensity (const SCMatD &mX, const SVecD &vDens, t_size k, const double dFact = 1) ;

//	Member Variables  //

		SMatD m_mCurM, m_mBestM ;

	} ;


	class CClust_C : virtual public CClust
	{
	public:

	protected:
		virtual double CalcObjFunc () ;
		virtual BOOL FindClustAssignment () ;
		void CalcClusterSize () ;

	} ;

	class CClust_CM : public CClust_C, public CClust_M
	{
	public:
		CClust_CM (t_size p, t_size k, double *pdM) ;

	protected:
		virtual void EstimClustParams () ;
		virtual double CalcObjFunc () { return CClust_C::CalcObjFunc () ; }
		virtual void EstimInitClustParams (int k, const SCVecN &vNIdx) { CClust_M::EstimInitClustParams (k, vNIdx) ;}
		virtual void SaveCurResult () ;
		virtual void CalcDensity (const SCMatD &mX, const SVecD &vDens, t_size k, const double dFact = 1) { CClust_M::CalcDensity (mX, vDens, k, dFact) ; }
		virtual BOOL FindClustAssignment () { return CClust_C::FindClustAssignment () ; } ;
	} ;

	class UOP	//	user defined operators for class EO
	{
	public:
		class AaC_BpaC				{ CALC_3_2(void) { a = c; b += c; } } ;
		class inc_a_if_b_equals_c	{ CALC_3_1(void) { if (b == (TB) c) a += 1 ; } } ;
		class inc_a_if_b_leq_c		{ CALC_3_1(void) { if (b <= (TB) c) a += 1 ; } } ;
		class Apa_logB				{ CALC_2_1(void) { a += log (b) ; } } ;				//	used in CalcObjFunc;
		class Apa_sqr_BmC			{ CALC_3_1(void) { a += sm_sqr (b * c) ; } } ;		//	used in CalcDensity;
		class Apa_sqr_BsC			{ CALC_3_1(void) { a += sm_sqr (b - c) ; } } ;		//	used in CalcDensity;
		class Aa_Bm_exp_Adm2		{ CALC_2_1(void) { a = b * exp (a / -2) ; } } ;		//	used in CalcDensity;
		class neg_log				{ CALC_2_1(void) { a = -log (b) ; } } ;
		class neg_log_0set0			{ CALC_2_1(void) { a = (b > 0) ? -log (b) : 0 ; } } ;
		class Apa_log_B				{ CALC_2_1(void) { a += log (b) ; } } ;
		class Apa_log_B_limit0		{ CALC_2_1(void) { a += log ((b >= 0) ? b : 0) ; } } ;
//		class Apa_sqrt_BmC			{ CALC_3_1(void) { a += sqrt (b * c) ; } } ;
//		class AmP_exp_Bd2			{ CALC_2_1(void) { a *= exp (b / 2) ; } } ;			//	used in CalcDensity;
	} ;

#endif	//	#ifndef INC_CLUST_H
