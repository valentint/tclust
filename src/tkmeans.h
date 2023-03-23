#include "clust.h"


	class CTKMeans : public CClust_CM, public CClust_N
	{
		public:

			CTKMeans (t_size n, t_size p, t_size k, double dAlpha, double dZeroTol, double *pdX, int *pnAssign, double *pdClustSize, double *pdClustWeights, int nEqualWeights, int nTrace, double *pdM) ;

		protected:

			virtual void EstimClustParams () { CClust_CM::EstimClustParams () ; }
			virtual void EstimInitClustParams (int k, const SCVecN &vNIdx) { CClust_CM::EstimInitClustParams (k, vNIdx) ;}
			virtual void SaveCurResult () {CClust_CM::SaveCurResult () ; }
			virtual BOOL FindClustAssignment () { return CClust_CM::FindClustAssignment () ; } ;

			virtual double CalcObjFunc () ; 
			virtual void CalcDensity (const SCMatD &mX, const SVecD &vDens, t_size k, const double dFact = 1) ;

            //VT::08.05.2018
    		virtual void FindInitClustAssignment () ;

			virtual void FindOutliers (const SVecD &vDisc, const SVecN &vInd) { CClust_N::FindOutliers (vDisc, vInd) ; }
			virtual void FindNearestClust (const SVecD &vDisc, const SVecN &vInd) {CClust_N::FindNearestClust (vDisc, vInd) ; }
			virtual void select_cluster (double &dDisc, int &nInd, const SCVecD &row) {CClust_N::select_cluster (dDisc, nInd, row) ; }
	} ;
