#include "clust.h"

//////////////
//  CClust  //
//////////////

	CClust::CClust ()
		: m_dAlpha (0), m_dZeroTol (0), m_dDensFact (0), m_dPInv (0), m_dZeroTolSqrt (0), m_dwNoTrim (0), m_dwTrim (0)
	{
		//	This constructor must never be called!!
		THROW (0) ;
	}

	CClust::CClust (t_size n, t_size p, t_size k, double dAlpha, double dZeroTol, double *pdX, int *pnAssign, double *pdClustSize, double *pdClustWeights, int nEqualWeights, int nTrace)
				//	parameters
		: m_n (n), m_p (p), m_K (k)
		, m_nEqualWeights (nEqualWeights), m_nTrace (nTrace)
		, m_dAlpha (dAlpha)
		, m_dZeroTol (dZeroTol)

			//	computing constants
		, m_dDensFact (pow (2 * meal_PI (), m_p / -2.0))
		, m_dPInv (1.0 / m_p)
		, m_dZeroTolSqrt (sqrt(m_dZeroTol))
		, m_dwNoTrim ((t_size) floor (m_n * (1-m_dAlpha)))
		, m_dwTrim (m_n-m_dwNoTrim )

			//	initializing vectors...
		, m_vInd (m_n), m_vIndBest (pnAssign, m_n), m_vIndOld (m_n), m_vCurInd (m_p + 1), m_vRank (m_n)
		, m_vWeights (m_K), m_vBestWeights (pdClustWeights, m_K), m_vClustSize (m_K), m_vClustSizeBest (pdClustSize, m_K), m_vDisc (m_n)

			//	initializing matrices...
		, m_mX (pdX, m_n, m_p), m_mLL (m_n, m_K)//, m_mEVal (m_p, m_K), m_mZ (m_n, m_K), m_mZOld (m_n, m_K)
	{
		meal_GetRNGstate () ;
	}

	CClust::~CClust ()
	{
		meal_PutRNGstate () ;
	}

#define SET_SAFE(PTR, VAL) {if (PTR) *PTR = VAL ; }
#define SET_SAFE_PRT(PTRS, PTR) {if (PTR) PTRS = PTR ; }
#define SET_SAFE_ARR(PTR, IDX, VAL) {if (PTR) PTR[IDX] = VAL ; }


	CClust &CClust::SetPtr (int *pnConvCount, int *pnIterSuccess, int *pnCode, int *pnErrExc, double *pdBestObj)
	{
		SET_SAFE_PRT (m_pnConvCount, pnConvCount)
		SET_SAFE_PRT (m_pnIterSuccess, pnIterSuccess)
		SET_SAFE_PRT (m_pnCode, pnCode)
		SET_SAFE_PRT (m_pnErrExc, pnErrExc)
		SET_SAFE_PRT (m_pdBestObj, pdBestObj)

		return *this ;
	}

	void CClust::calc (int nIter, int nKSteps, int * const pnConvER, double * const pdObjER)
	{	
		if (!CheckParams ())
		{
			SET_SAFE (m_pnErrExc, 1) ;
			return  ;
		}
		int i, j ;

		m_vIndOld.Reset (NAI) ;

		double dLastObj = 0 ;

		//CheckRestrictions () ;	//	2do: -> virtual

		for (i = 0; i < nIter; i++)
		{
            //VT::08.05.2018
            // meal_printf("MY-TRACE ... Entering FindInitClustAssignment() from CCClust::calc() ...\n");

			FindInitClustAssignment () ;	//	2do: -> virtual; call FindInitClustSize_R
			//FindInitClustSize_R () ;	
			//FindInitClustSize () ;

			for (j = 0; 1; j++)
			{
				if (!RestrictClustParam ())	//	2do: rename to restrClustParam
				{
					if (i)	//	2do: BUG: shouldn't this be "j" !?!?
					{
						SaveCurResult0 (CalcObjFunc (), 2) ;
						return ;
					}
					else
						SetSingularIniParams () ;	//	was SetAllCovMatsIdent
				}

				if (!FindClustAssignment () ||		//	returns false, if the cluster assignment has not changed -> break 
					j == nKSteps)					//	max number of iterations reached -> break
					break ;

				if (int(m_nTrace) >= 2)
				{
					double dCurObj = CalcObjFunc () ;

					if (j && dLastObj > dCurObj)
						meal_printf ("Objective function dropped from %.10f to %.10f in (%d/%d)\n", dLastObj, dCurObj, i, j) ;
					else
						meal_printf ("Objective function %.10f in (%d/%d)\n", dCurObj, i, j) ;

					dLastObj = dCurObj ;
				}

				EstimClustParams () ;
			}

			BOOL bConv = j < nKSteps ;	//	did this iteration converge?

			if (bConv && m_pnConvCount)
				++ *m_pnConvCount ;

			SET_SAFE_ARR (pnConvER, i, bConv) ;

			double dCurObj = CalcObjFunc () ;

			SET_SAFE_ARR (pdObjER, i, dCurObj) ;

			if (!i || dCurObj > m_dBestObj)			//	store if the current solution is the best solution so far..
				SaveCurResult0 (dCurObj, j >= nKSteps) ;
		}
	}

	BOOL CClust::CheckParams ()
	{		//	checks the input parameters for validity

		if (m_n < 1)
		{
			meal_printf ("Input parameter error: n must be >= 1\n") ;
			return FALSE ;
		}
		if (m_p < 1)
		{
			meal_printf ("Input parameter error: p must be >= 1\n") ;
			return FALSE ;
		}
		if (m_K < 1)
		{
			meal_printf ("Input parameter error: k must be >= 1\n") ;
			return FALSE ;
		}

		if (m_dAlpha < 0 || m_dAlpha > 1)
		{
			meal_printf ("Input parameter error: alpha must be between 0 and 1\n") ;
			return FALSE ;
		}

		if (!m_dwNoTrim || m_dAlpha >= 1.0)
		{
			meal_printf ("Input parameter error: alpha was chosen too large (all observations were trimmed)\n") ;
			return FALSE ;
		}
		return TRUE ;
	}

	void CClust::FindInitClustAssignment ()
	{		//	finds initial cluster assignment
		t_size k ;

		SVecN vTempN (m_aTemp [0], m_n) ;
		SMatD mDCurCluster (m_aTemp [0], m_p + 1, m_p) ;
		SVecN vNCurIdx (m_aTemp [1], m_p + 1) ;

        //VT::08.05.2018
        // meal_printf("MY-TRACE ... In CClust::FindInitClustAssignment() ...\n");

		for (k = 0; k < m_K; ++k)
		{								//	for all clusters
										//	finds p+1 observations for forming the initial cluster
			SampleNoReplace (m_p + 1, m_mX.nrow (), vNCurIdx, vTempN) ;

			EstimInitClustParams (k, vNCurIdx) ;


/*			SVecD vDCurMean = m_mCurM.GetColRef (k) ;	//	2do: use m_avCurM [k]	-> Implement	SCVecArray
			vDCurMean.Reset (0) ;
										//	stores the subsetted X matrix into MDCurCluster, and sums up it's columns to vcDCurMean.
			EO<UOP::AaC_BpaC>::MsVetMcdVcei (!mDCurCluster, *vDCurMean, m_mX, vNCurIdx) ;

										//	Divide the colSums of MDCurCluster by it's number of rows. -> vDCurMean is a mean vector.
			EO<SOP::a_divide>::VSc (*vDCurMean, mDCurCluster.nrow ()) ;

										//	Centers the mDCurCluster - matrix
			EO<SOP::a_minus>::MVcet (!mDCurCluster, vDCurMean) ;

										//	Calculating the covariance matrix
			cov_centered_NC (!m_amCurS[k], mDCurCluster, dCorrFact) ;
*/		}

		FindInitClustSize_R () ;
	}


	void CClust::LoadCluster (SMatD &c, t_size k)
	{
		int nCSize = 0 ;
		EO<UOP::inc_a_if_b_equals_c>::SScVc (nCSize, k, m_vInd) ;

		c.Require (nCSize, m_p) ;

		double *pdCurColX = m_mX ;
		double *const pdEndX = m_mX.GetDataEnd () ;
		int * const pnInd = m_vInd ;
		int * pnCurInd ;
		int * const pnEndInd = m_vInd.GetDataEnd () ;
		double *pdC = c ;

		while (pdCurColX < pdEndX)
		{
			pnCurInd = pnInd ;

			for (; pnCurInd < pnEndInd; ++pnCurInd)
				if (*pnCurInd == (int) k)
				{
					*pdC = pdCurColX [pnCurInd - pnInd] ;
					++pdC ;
				}

			pdCurColX += m_mX.GetColInc () ;
		}
			//	now we should have cluster k in matrix c
	}



	void CClust::FindInitClustSize_R ()
	{	//	as FindInitClustSize, but returns the same results as the current R - code 2 be deleted!!

		if (m_nEqualWeights)			//	for equal sized clusters
		{
			m_vClustSize.Reset (m_dwNoTrim / (double) m_K) ;
			m_vWeights.Reset (1.0 / (double) m_K) ;
		}
		else	
		{
			SVecD v (m_aTemp[0], m_K) ;
			SVecN vOrder (m_aTemp[1], m_K) ;

			t_size k, n ;

			runif_r (*v) ;		//	calc cumulated sum of cluster propability (proportional cluster size)
			meal_sort_order (v, vOrder, v.size ()) ;
			cumsum_r (*v) ;
			EO<SOP::a_divide>::VSc (*v, v.GetValue (0)) ;

//RtprintVec (v);										//	this is as in R (reverse order)
			m_vClustSize.Reset (0) ;

			double *pdV = v ;
			double *pdClustSize = m_vClustSize ;
			for (n = m_dwNoTrim; n > 0; n--)			//	
			{
				double dCur = meal_unif_rand () ;		//	get rand. cluster assignment
														//	1- runif for R compatibility (gives similar results - still differences - why?)
				for (k = m_K - 1; k != NAI; k--)		//	check in which group it can be found
					if (dCur <= pdV[k])
					{
						pdClustSize [k] += 1 ;
						break ;
					}
			}

			v.Copy_NC (m_vClustSize) ;
			int *pnOrder = vOrder ;

			for (k = m_K - 1; k != NAI; k--)
				pdClustSize[m_K - 1 - pnOrder[k]] = pdV[k] ;

			EO<SOP::divide_r>::VScVc (*m_vWeights, m_dwNoTrim, m_vClustSize) ;

//	if (int(m_dwTrace) >= 10)
//	{
//		RtprintVec (m_vClustSize, "ClustSize:\t");
//		RtprintVec (m_vWeights, "ClustWeights:\t");
//	}
		}
	}

	void CClust::FindInitClustSize ()
	{
		if (m_nEqualWeights)			//	for equal sized clusters
		{
			m_vClustSize.Reset (m_dwNoTrim / (double) m_K) ;
			m_vWeights.Reset (1.0 / (double) m_K) ;
		}
		else	
		{
			SVecD v (m_aTemp[0], m_K) ;

			t_size k, n ;

			runif_r (*v) ;
			sort (*v) ;							//	sorting is good when m_K is huge.
												//	when checking which cluster a guess belongs to, the big clusters are checked earlier
			cumsum_r (*v) ;
			EO<SOP::a_divide>::VSc (*v, v.GetValue (0)) ;	//	norm vector v (sum = 1 - i.e. first element is 1, as this these are cumulated sums)

			double *pdV = v ;
			double *pdClustSize = m_vClustSize ;

			m_vClustSize.Reset (0) ;
			for (n = m_dwNoTrim; n > 0; n--)			//	
			{
				double dCur = meal_unif_rand () ;		//	get rand. cluster assignment
														//		1- runif for R compatibility (gives similar results - still differences - why?)
				for (k = m_K - 1; k != NAI -1; k--)		//	check in which group it can be found
					if (dCur <= pdV[k])
					{
						pdClustSize [k] += 1 ;
						break ;
					}
			}
			EO<SOP::divide_r>::VScVc (*m_vWeights, m_dwNoTrim, m_vClustSize) ;
		}
	}

	void CClust::SaveCurResult0 (double dCurObj, int nCode)
	{
		SET_SAFE (m_pnCode, nCode) ;

		m_dBestObj = dCurObj ;
		SET_SAFE (m_pdBestObj, dCurObj) ;

		m_vBestWeights.Copy_NC (m_vWeights) ;
		m_vClustSizeBest.Copy_NC (m_vClustSize) ;
		m_vIndBest.Copy_NC (m_vInd) ;

		SaveCurResult () ;
	}

////////////////
//  CClust_N  //
////////////////

	void CClust_N::FindOutliers (const SVecD &vDisc, const SVecN &vInd)
	{
		if (!m_dwTrim)
			return ;
		int * const pnInd = vInd ;								//	setting the TrimIdx - array
		SVecN vOrder (m_aTemp[0], vDisc.size ()) ;
		int * const pnOrder = vOrder ;
		double *pdDisc = vDisc ;
		meal_sort_order (pdDisc, pnOrder, vDisc.size ()) ;		//	calculate the order of m_vDisc and store it in vOrder

		t_size i ;
		for (i = m_dwTrim - 1; i != NAI; i--)
			pnInd[pnOrder [i]] = -1 ;
	}


	void CClust_N::FindNearestClust (const SVecD &vDisc, const SVecN &vInd)	//	old version - without tie-handling!
	{
		t_size c ;

		vInd.Reset (0) ;
		vDisc.Copy_NC (m_mLL.GetColRef (0)) ;

		if (m_K == 1)		//	if we only have one cluster we already have found the solution by assigning the first cluster to each observation.
			return ;

		const double *pdCurLL = m_mLL.GetData (0, 1) ;
		double *pdCurDisc ;
		int *pnCurInd, *pnCurIndEnd = vInd.GetDataEnd () ;

		for (c = 1; c < m_mLL.ncol (); ++c) 
		{
			pnCurInd = vInd ;
			pdCurDisc = vDisc ;

			while (pnCurInd < pnCurIndEnd)
			{
				if (sm_setmax_b (*pdCurDisc, *pdCurLL))
					*pnCurInd = c ;
				++ pnCurInd ;
				++ pdCurDisc ;
				++ pdCurLL ;
			}
		}
	}

	void CClust_N::select_cluster (double &dDisc, int &nInd, const SCVecD &row)	//	no ties - handling
	{
		double const *pdRow = row, *pdMax = pdRow, *const pdEndRow = row.GetDataEnd ()  ;

		for (++pdRow; pdRow < pdEndRow; ++pdRow)
			if (*pdMax < *pdRow)
				pdMax = pdRow ;
		nInd = pdMax - row ;
		dDisc = *pdMax ;
	}

////////////////
//  CClust_T  //
////////////////

	void CClust_T::FindOutliers (const SVecD &vDisc, const SVecN &vInd)
	{
		if (!m_dwTrim)
			return ;

		SVecN vOrder (m_aTemp[0], vDisc.size ()) ;
		int * const pnOrder = vOrder, * const pnInd = vInd ;

		double *pdDisc = vDisc ;
		meal_sort_order (pdDisc, pnOrder, vDisc.size ()) ;		//	calculate the order of m_vDisc and store it in vOrder

		ASSERT (m_dwTrim < m_n) ;								//	just to be sure..

		double dThreshold = pdDisc[m_dwTrim - 1] + m_dZeroTol ;

		t_size i ;
		if (pdDisc[m_dwTrim] > dThreshold)						//	we don't have any ties around the trimming threshold.
		{
			for (i = m_dwTrim - 1; i!= NAI; i--)
				pnInd [pnOrder[i]] = -1 ;
			return ;
		}

		double dThresholdL = pdDisc[m_dwTrim - 1] - m_dZeroTol ;

		for (i = 0; i < m_n && pdDisc[i] < dThresholdL; ++i)
			pnInd [pnOrder[i]] = -1 ;

		t_size dwSmaller = i ;									//	1 + number of smaller values than dThresholdL

		for (; i < m_n; ++i)
			if (pdDisc[i] > dThreshold)
				break ;

		ASSERT (m_dwTrim > dwSmaller) ;

		t_size dwProblem = i - dwSmaller ;						//	the number of ties around the trimming threshold
		t_size dwChoose = m_dwTrim - dwSmaller ;				//	the number of ties which have to be trimmed 

		SVecN vIdx (m_aTemp[1], dwChoose) ;
		SVecN vTemp (m_aTemp [2], dwProblem) ;

												//	randomly draw these values
		int *pnIdx = vIdx, * const pnEndIdx = vIdx.GetDataEnd () ; 
		SampleNoReplace (dwChoose, dwProblem, pnIdx, vTemp) ;

		for (; pnIdx < pnEndIdx; ++pnIdx)
			pnInd [pnOrder [dwSmaller + *pnIdx]] = -1 ;
	}

	void CClust_T::FindNearestClust (const SVecD &vDisc, const SVecN &vInd)	//	new version - implements tie-handling
	{
		t_size r ;

		double *pdCurDisc = vDisc ;
		int *pnCurInd = vInd ;

		SVecD curRow (m_aTemp[1], m_K) ;

		for (r = 0; r < vDisc.size(); ++r)									//	for each row of m_mLL
		{
			CopyRow (*curRow, m_mLL, r) ;
			select_cluster (*pdCurDisc, *pnCurInd, curRow) ;				//	select the cluster corresponding to the largest value of the current row of m_mLL
			++pdCurDisc ;
			++pnCurInd ;
		}
	}

	void CClust_T::select_cluster (double &dDisc, int &nInd, const SCVecD &row)	//	temp0
	{
		double const *pdRow = row, *pdMax = pdRow, *const pdEndRow = row.GetDataEnd ()  ;

		int nTies = 0 ;
		++pdRow ;
		for (; pdRow < pdEndRow; ++pdRow)											//	find the maximum
		{
			if (*pdMax <= *pdRow)
			{
				nTies = *pdRow - *pdMax <= m_dZeroTol ;
				pdMax = pdRow ;
			}
		}

		if (!nTies)															//	we do have exactly one maximum (no ties at the maximum!)
		{
			nInd = pdMax - row ;
			dDisc = *pdMax ;
			return ;
		}
																			//	everything >= dThreshold is considered as possible maximum
		const double dThreshold = *pdMax - m_dZeroTol ;

		SVec<double const *> vIdx (m_aTemp [0], row.size ()) ;

		const double ** const pdIdx = vIdx, ** pdCurIdx = pdIdx ;

		pdMax = row ;														//	from now on pdMax refers to the beginning of vector "row".
		for (pdRow = pdMax; pdRow < pdEndRow; ++pdRow)
			if (dThreshold <= *pdRow)										//	if the current value is greater or equal than the the threshold
			{
				*pdCurIdx = pdRow ;											//	save the corresponding pointer
				++pdCurIdx ;
			}

		nInd = pdCurIdx - vIdx ;												//	debug line -> 2do: delete!
		nInd = int (meal_unif_rand () * (pdCurIdx - vIdx)) ;					//	get a random index 

		nInd = pdIdx [nInd] - pdMax ;										//	select a maximum according to the drawn random index.
		dDisc = pdMax[nInd] ;
	}

////////////////
//  CClust_M  //
////////////////

	CClust_M::CClust_M (double *pdM) 
		: m_mCurM (m_p, m_K), m_mBestM (pdM, m_p, m_K)
	{

	}

/*	CClust_M::CClust_M (t_size n, t_size p, t_size k, double dAlpha, double dZeroTol, double *pdX, int *pnAssign, double *pdClustSize, double *pdClustWeights, double *pdObjER, int *pnConvER, double *pdM)
		: CClust (n, p, k, dAlpha, dZeroTol, pdX, pnAssign, pdClustSize, pdClustWeights, pdObjER, pnConvER)
		, m_mCurM (m_p, m_K), m_mBestM (pdM, m_p, m_K)
	{

	}
*/
	void CClust_M::EstimInitClustParams (int k, const SCVecN &vNIdx)
	{
		SVecD vDCurMean = m_mCurM.GetColRef (k) ;	//	2do: use m_avCurM [k]	-> Implement	SCVecArray

		vDCurMean.Reset (0) ;

		EO<SOP::a_plus>::VetMcdVcei(*vDCurMean, m_mX, vNIdx) ;

		//	Divide the colSums of MDCurCluster by it's number of rows. -> vDCurMean is a mean vector.
		EO<SOP::a_divide>::VSc (*vDCurMean, vNIdx.size ()) ;
	}

	void CClust_M::SaveCurResult ()
	{
		m_mBestM.Copy_NC (m_mCurM) ;
	}

	void CClust_M::CalcDensity (const SCMatD &mX, const SVecD &vDens, t_size k, const double dFact)
	{
		const SVecD &vCurM = m_mCurM.GetColRef (k) ;

		vDens.Reset (0) ;

		EO<UOP::Apa_sqr_BsC>::VMcVct_NC (*vDens, mX, vCurM) ;				//	vDens <- rowSums (sqr (mX - vCurM))

		EO<UOP::Aa_Bm_exp_Adm2>::VSc (*vDens, dFact * m_dDensFact) ;
	}

////////////////
//  CClust_S  //	
////////////////

	CClust_S::CClust_S (double *pdM, double *pdS)
		: CClust_M (pdM)
		, m_dPd_Pp1 ((m_p) / (m_p  + 1.0))
		, m_mEVal (m_p, m_K)
		, m_amEVec (m_p, m_p, m_K), m_amCurS (m_p, m_p, m_K),m_amBestS (pdS, m_p, m_p, m_K)
	{
		
	}

/*	CClust_S::CClust_S (t_size n, t_size p, t_size k, double dAlpha, double dZeroTol, double *pdX, int *pnAssign, double *pdClustSize, double *pdClustWeights, double *pdObjER, int *pnConvER, double *pdM, double *pdS)
		: CClust (n, p, k, dAlpha, dZeroTol, pdX, pnAssign, pdClustSize, pdClustWeights, pdObjER, pnConvER)
		, CClust_M (n, p, k, dAlpha, dZeroTol, pdX, pnAssign, pdClustSize, pdClustWeights, pdObjER, pnConvER, pdM)
		, m_amEVec (m_p, m_p, m_K), m_amCurS (m_p, m_p, m_K),m_amBestS (pdS, m_p, m_p, m_K)
		, m_dPd_Pp1 ((m_p) / (m_p  + 1.0))
		, m_mEVal (m_p, m_K)
		{
			
		}
*/
	void CClust_S::EstimInitClustParams (int k, const SCVecN &vNIdx)
	{
		SMatD mDCurCluster (m_aTemp [0], m_p + 1, m_p) ;

		SVecD vDCurMean = m_mCurM.GetColRef (k) ;	//	2do: use m_avCurM [k]	-> Implement	SCVecArray
		vDCurMean.Reset (0) ;
									//	stores the subsetted X matrix into MDCurCluster, and sums up it's columns to vcDCurMean.
		EO<UOP::AaC_BpaC>::MsVetMcdVcei (!mDCurCluster, *vDCurMean, m_mX, vNIdx) ;

									//	Divide the colSums of MDCurCluster by it's number of rows. -> vDCurMean is a mean vector.
		EO<SOP::a_divide>::VSc (*vDCurMean, mDCurCluster.nrow ()) ;

									//	Centers the mDCurCluster - matrix
		EO<SOP::a_minus>::MVcet (!mDCurCluster, vDCurMean) ;

									//	Calculating the covariance matrix
		cov_centered_NC (!m_amCurS[k], mDCurCluster, m_dPd_Pp1) ;
	
	}

	void CClust_S::SaveCurResult ()
	{
		t_size k ;
		for (k = m_K - 1; k != NAI; k--)
			m_amBestS[k].Copy_NC (m_amCurS[k]) ;	//	2do: use tensor!

		CClust_M::SaveCurResult () ;
	}

	void CClust_S::SetSingularIniParams ()
	{
		//	setting all covariance matrices to ident matrix and correspondingly changes the eigenvalues & vectors

		t_size k ;
		for (k = m_K - 1; k != NAI; k--)
		{
			SetDiag_sq_NC (!m_amCurS[k]) ;
			SetAntiDiag_sq_NC (!m_amEVec[k]) ;
		}
		m_mEVal.Reset (1) ;
	}

	void CClust_S::CalcDensity (const SCMatD &mX, const SVecD &vDens, t_size k, const double dFact)
	{

		const SVecD &vCurM = m_mCurM.GetColRef (k) ;
		const SCVecD &vEVal = m_mEVal.GetColRef (k) ;
		const SCMatD &mEVec = m_amEVec [k] ;

		SMatD mZ (m_aTemp[0], mX.dim ()) ;
		SVecD vTemp (m_aTemp[1], m_p) ;
		SMatD mXc (m_aTemp[2], mX.dim ()) ;

		EO<SOP::subtract>::MMcVct_NC (!mXc, mX, vCurM) ;					//	centering X

//		m_mTempNP1.Reshape (mX.nrow (), mX.ncol ()) ;
//		FC_ElOp<FC::FC_minus, DWORD>::OpMV_row (mX, vCurM, m_mTempNP1) ;	// TempNP1 = X centered

		sme_matmult_NC (mXc, mEVec, !mZ) ;

//		m_mTempNP2.Reshape (mX.nrow (), mX.ncol ()) ;
//		matmultmat (m_mTempNP1, mEVec, m_mTempNP2) ;						//	X scaled (Z)

//		m_vTempN1.Reshape (m_p) ;

		EO<SOP::pow_r>::VScVc (*vTemp, -0.5, vEVal) ;						//	vTemp <- eval^(0.5)

//		FC_ElOp<FC::FC_pow, DWORD>::OpVE (vEVal, -0.5, m_vTempN1) ;			//	Gamma ^-0.5

		vDens.Reset (0) ;
		EO<UOP::Apa_sqr_BmC>::VMcVct_NC (*vDens, mZ, vTemp) ;				//	vDens <- rowSums (sqr (mZ %*% vTemp))

//		FC_ElOpAs<FC::FC_multiply>::OpMV_row (m_mTempNP2, m_vTempN1) ;		//	Z %*% Gamm ^-0.5
//		FC_ElOpAs<FC::FC_sqr>::OpM (m_mTempNP2) ;							//	sqr (Z %*% Gamm ^-0.5)
//		rowSums (m_mTempNP2, vDens) ;										//	 == mahalanobis

		const double dDet = prod (vTemp) ;
		//EO<SOP::a_multiply>::SVc (dDet, vTemp) ;							//	dDet <- prod (exp (eval^-(0.5)/2))

//		FC_ElOpAs<FC::FC_multiply>::OpVE (vDens, -0.5) ;					//	mahalanobis * 0.5
//		FC_ElOpAs<FC::FC_exp>::OpV (vDens) ;								//	exp (maha / 2)
//		double dDet = prod (m_vTempN1) ;									//	||Sigma||

		EO<UOP::Aa_Bm_exp_Adm2>::VSc (*vDens, dDet * dFact * m_dDensFact) ;
	}

////////////////
//  CClust_C  //
////////////////

	double CClust_C::CalcObjFunc ()
	{
		double dObj = 0 ;

		t_size k ;

		SVecD vDensity (m_aTemp[3], 0) ;
		SMatD mCurX (m_aTemp[4], m_n, m_p) ;

//		double *pdClustSize = m_vClustSize ;

		for (k = m_K - 1; k != NAI; --k)
		//for (k = 0; k < m_K; ++k)
		{
			LoadCluster (mCurX, k) ;

			t_size dwClustSize = mCurX.nrow () ;	//	2do: is this necessary? only summing up the according cells would be sufficient..

			if (!dwClustSize)
				continue ;

			//double &dClustSize = pdClustSize [k] = dwClustSize ;		//	2do: is this assignment necessary? 
//			double &dClustSize = pdClustSize [k];

//			if (dClustSize <= m_dZeroTol)
//				continue ;

			vDensity.Require (dwClustSize) ;

			CalcDensity (mCurX, vDensity, k) ;

			EO<UOP::Apa_logB>::SVc (dObj, vDensity) ;

			if (!m_nEqualWeights)	//	this can be avoided when equal-sized clusters are expected (only adding a fixed value to each object. function)
				dObj += dwClustSize * log (dwClustSize / (double) m_dwNoTrim) ;
		}

		return dObj ;
	}


	BOOL CClust_C::FindClustAssignment ()
	{		//	calculates the vector m_vIndTrim
		t_size k ;
		for (k = m_K - 1; k != NAI; k--)
			CalcDensity (m_mX, m_mLL.GetColRef (k), k, m_vWeights (k)) ;		//	2do: Implement VecArray

		SVecD vDisc (m_aTemp[3], m_n) ;

		FindNearestClust (vDisc, m_vInd) ;

		FindOutliers (vDisc, m_vInd) ;

		if (equal (m_vInd, m_vIndOld))				//	check the change of group assignment AFTER the calculation of outliers. 
			return FALSE ;							//	otherwise the alg will stop after the 1st iteration when k == 1, 
		m_vIndOld.Copy_NC (m_vInd) ;				//	since there wouldn't even be the change of the change of the cluster assignment

		CalcClusterSize () ;						//	Cluster sizes are always calculated. Just don't calc the weigths for (m_dwEqualWeights)

		return TRUE ;
	}

	void CClust_C::CalcClusterSize ()
	{
		m_vClustSize.Reset (0) ;

		double * const pdClustSize = m_vClustSize ;
		int *pnInd = m_vInd ;
		int * const pnEndInd = m_vInd.GetDataEnd () ;

		for (; pnInd < pnEndInd; ++pnInd)
			if (*pnInd != -1)		//	not an outlier
				pdClustSize[*pnInd] += 1 ;

//		t_size k ;
//		for (k = m_p - 1; k != NAI; k--)
//			NotifyClusterSize (k, pdClustSize[k]) ;

		if (!m_nEqualWeights)
			EO<SOP::divide_r>::VScVc (*m_vWeights, m_dwNoTrim, m_vClustSize) ;	//	m_vWeights <- m_vClustSize / m_dwNoTrim
	}



////////////////
//  CClust_F  //
////////////////

	CClust_F::CClust_F (double dM, double *pdZ)
		: m_dM (dM)
		, m_dMm1Inv (1 / (m_dM - 1))
		, m_mZ (m_n, m_K), m_mZ_best (pdZ, m_n, m_K), m_mZOld (m_n, m_K)
	{

	}

	void CClust_F::SaveCurResult ()
	{
		m_mZ_best.Copy_NC (m_mZ) ;
	}

	double CClust_F::CalcObjFunc ()
	{
		SVecD vDensity (m_aTemp[3], m_n) ;

		t_size k ;
		double dObj = 0, *pdDens, * const pdEndDens = vDensity.GetDataEnd () ;
		const double *pdClustSize = m_vClustSize ;

//		for (k = m_K - 1; k != NAI; --k)
		for (k = 0; k < m_K; ++k)
		{
			const double &dClustSize = pdClustSize [k] ;
			if (dClustSize <= m_dZeroTol)
				continue ;

			CalcDensity (m_mX, vDensity, k) ;	//	2do: introduce SCVecArray

			double *pdZ = m_mZ.GetData (0, k) ;
			double dZSum = 0 ;

			pdDens = vDensity ;
			while (pdDens < pdEndDens)
			{
				if (*pdDens <= 0.0)
				{
					dZSum += *pdZ ;
					if (dZSum > m_dZeroTol)
						return meal_NegInf () ;
				}
				else
					dObj += *pdZ * log (*pdDens) ;

				++pdZ ;
				++pdDens ;
			}

			if (!m_nEqualWeights)	//	this can be avoided when equal-sized clusters are expected (only adding a fixed value to each object. function)
				dObj += dClustSize * log (m_vWeights (k)) ;
		}

		return dObj ;
	}


	BOOL CClust_F::FindClustAssignment ()
	{
		t_size k ;
		for (k = m_K - 1; k != NAI; k--)
			CalcDensity (m_mX, m_mLL.GetColRef (k), k, m_vWeights (k)) ;

		SVecD vLLRow (m_aTemp[3], m_K) ;
		SVecD vZRow (m_aTemp[4], m_K) ;
		SVecD vDisc (m_aTemp[5], m_n) ;

		double *pdDisc = vDisc ;
		int *pnAssign = m_vInd ;

		for (k = 0; k < m_n; ++k)
		{
			CopyRow (*vLLRow, m_mLL, k) ;
			CalcFuzzyRow (vLLRow, vZRow, *pdDisc, *pnAssign) ;
			CopyRow (!m_mZ, vZRow, k) ;

			++pdDisc ;
			++pnAssign ;
		}

		FindOutliers (vDisc, m_vInd) ;

		pnAssign = m_vInd ;
		for (k = 0; k < m_n; ++k)														//	setting the rows of m_mZ which correspond to trimmed observations to zero.
		{
			if (*pnAssign == -1)
				ResetRow (!m_mZ, 0, k) ;
			++pnAssign ;
		}

		if (equal (m_mZ, m_mZOld))														//	the m_mZ matrix didn't change 
			return FALSE ;

		m_mZOld.Copy (m_mZ) ;															//	save current m_mZ - matrix

		colSums_NC (*m_vClustSize, m_mZ) ;												//	calculate cluster sizes

//		double * const pdClustSize = m_vClustSize ;
//		for (k = m_p - 1; k != NAI; k--)
//			NotifyClusterSize (k, pdClustSize[k]) ;

		if (!m_nEqualWeights)
			EO<SOP::divide_r>::VScVc (*m_vWeights, sum (m_vClustSize), m_vClustSize) ;	//	calculate cluster weights (if necessary)

		return TRUE ;
	}

	void CClust_F::CalcFuzzyRow (const SCVecD &ll, const SVecD &z, double &dDisc, int &nInd)	//	temp1
	{
		ASSERT (ll.size () == z.size()) ;
		select_cluster (dDisc, nInd, ll) ;

		if (dDisc >= 1
			//|| dDisc <= 0
		)
		{
			SetCatZ (z, nInd) ;
			if (dDisc > 0)
				dDisc = log (dDisc) ;
			return ;
		}

		SVecD ll_log (m_aTemp[1], ll.dim ()) ;
		EO<UOP::neg_log_0set0>::VVc_NC (*ll_log, ll) ;

		double *const pdL = ll_log,  *pdL1 = pdL, *pdL2, * const pdEndL = ll_log.GetDataEnd () ;
		double *pdZ = z ;
		double dSumZ = 0 ;

		dDisc = 0 ;

		while  (pdL1 < pdEndL)
		{
			*pdZ = 0 ;
			if (*pdL1 > 0)		//	only if we have anything to add..
			{

				for (pdL2 = pdL; pdL2 < pdEndL; ++pdL2)
					if (*pdL2 > 0)
						*pdZ += pow (*pdL1 / *pdL2, m_dMm1Inv) ;

				if (*pdZ > 0)
				{
					dSumZ += *pdZ ;
					*pdZ = pow (*pdZ, -m_dM) ;
					dDisc -= *pdL1 * *pdZ ;	//	*pdL1 is the negative log..
				}
			}
			++pdL1 ;
			++pdZ ;
		}

		if (dSumZ <= 0)
		{
			z.Reset (1.0 / m_K) ;
			dDisc = 0 ;
			EO<UOP::Apa_log_B_limit0>::SVc (dDisc, ll) ;	//	again calculating logarithms, because zeros have been removed earliear.
			dDisc /= m_K ;
//			SetCatZ (z, nInd) ;								//	old, wrong implementation
//			dDisc = -pdL[nInd] ;
		}
	}

	void CClust_F::SetCatZ (const SVecD &z, int nIdx)
	{
		double *pdZ = z, *const pdEndZ = z.GetDataEnd () ;
		
		while (pdZ < pdEndZ)
		{
			*pdZ = !nIdx ;
			--nIdx ;
			++pdZ ;
		}
	}

/////////////////
//  CClust_CS  //
/////////////////

	CClust_CS::CClust_CS (double *pdM, double *pdS)
		: CClust_S (pdM, pdS)
	{

	}

	void CClust_CS::EstimClustParams ()
	{		//	calculates the mean and the cov matrix of the clusters currently defined by m_vIndTrim
		t_size k ;

		SMatD mCurX (m_aTemp[0]) ;						//	2do: eval max clustersize and require this matrix.

		for (k = m_K - 1; k != NAI; k--)
//		for (k = 0; k < m_K; k++)
		{
			t_size dwCurClustSize = (t_size) m_vClustSize (k) ;

			if (!dwCurClustSize)
				continue ;

			SVecD vDCurMean = m_mCurM.GetColRef (k) ;	//	2do: use m_avCurM [k]	-> Implement	SCVecArray

			vDCurMean.Reset (0) ;

			mCurX.Require (dwCurClustSize, m_p) ;

											//	stores the subsetted X matrix into mCurX, and sums up it's columns to vcDCurMean.
			EO<UOP::AaC_BpaC>::MsVetMcdScgVceg_NC (!mCurX, *vDCurMean, m_mX, k, m_vInd) ;

											//	Divide the colSums of mCurX by it's number of rows. -> vDCurMean is a mean vector.
			EO<SOP::a_divide>::VSc (*vDCurMean, dwCurClustSize) ;

			if (dwCurClustSize > 1)
			{
											//	Centers the mCurX - matrix
				EO<SOP::a_minus>::MVcet (!mCurX, vDCurMean) ;

											//	Calculating the covariance matrix
				cov_centered_NC (!m_amCurS[k], mCurX, ((dwCurClustSize - 1.0) / dwCurClustSize)) ;
			}
			else
				m_amCurS[k].Reset (0) ;
		}
	}

	void CClust_CS::SaveCurResult ()
	{
		CClust_C::SaveCurResult () ;
		CClust_S::SaveCurResult () ;
	}


/////////////////
//  CClust_CM  //
/////////////////

	CClust_CM::CClust_CM (double *pdM)
		: CClust_M (pdM)
	{

	}

	void CClust_CM::EstimClustParams ()
	{		//	calculates the mean and the cov matrix of the clusters currently defined by m_vIndTrim
		t_size k ;

		for (k = m_K - 1; k != NAI; k--)
//		for (k = 0; k < m_K; k++)
		{
			t_size dwCurClustSize = (t_size) m_vClustSize (k) ;

			if (!dwCurClustSize)
				continue ;

			SVecD vDCurMean = m_mCurM.GetColRef (k) ;	//	2do: use m_avCurM [k]	-> Implement	SCVecArray

			vDCurMean.Reset (0) ;

											//	stores the subsetted X matrix into mCurX, and sums up it's columns to vcDCurMean.
			EO<SOP::a_plus>::VetMcdScgVceg_NC (*vDCurMean, m_mX, k, m_vInd) ;

											//	Divide the colSums of mCurX by it's number of rows. -> vDCurMean is a mean vector.
			EO<SOP::a_divide>::VSc (*vDCurMean, dwCurClustSize) ;
		}
	}

	void CClust_CM::SaveCurResult ()
	{
		CClust_C::SaveCurResult () ;
		CClust_M::SaveCurResult () ;
	}

/////////////////
//  CClust_FS  //
/////////////////

	CClust_FS::CClust_FS (double dM, double *pdZ, double *pdM, double *pdS)
		: CClust_F (dM, pdZ), CClust_S (pdM, pdS)
	{
		
	}

	void CClust_FS::EstimClustParams ()
	{
		t_size k ;

		SMatD mXc (m_aTemp[0], m_mX.dim ()) ;
		const double *pClustSize = m_vClustSize ;

		for (k = m_K - 1; k != NAI; k--)
		{
			const double &dCurClustSize = pClustSize [k] ;
			if (dCurClustSize > m_dZeroTol)
			{
				SVecD vDCurMean = m_mCurM.GetColRef (k) ;	//	2do: use m_avCurM [k]	-> Implement	SCVecArray

					//R	iter$center[k,] = (t(iter$z_ij[,k]) %*% X) / iter$csize[k]
				vDCurMean.Reset (0) ;
				EO<SOP::ApaBmC>::VtMcVc_NC (*vDCurMean, m_mX, m_mZ.GetColRef (k)) ;
				EO<SOP::a_divide>::VSc (*vDCurMean, dCurClustSize) ;

					//R	X.c <- (X - matrix (iter$center[k,], ncol = pa$p, nrow = pa$n, byrow = TRUE))
				EO<SOP::subtract>::MMcVct_NC (!mXc, m_mX, vDCurMean) ;

					//R	iter$sigma[,,k] <- (t(X.c * iter$z_ij[,k]) %*% X.c) / iter$csize[k]
				sme_matmult_at_diagb_a (mXc, m_mZ.GetColRef (k), !m_amCurS [k]) ;
				EO<SOP::a_divide>::VSc (*m_amCurS [k], dCurClustSize) ;
			}
			else
				m_amCurS[k].Reset (0) ;
		}
	}

	void CClust_FS::SaveCurResult ()
	{
		CClust_F::SaveCurResult () ;
		CClust_S::SaveCurResult () ;
	}


/////////////////
//  CClust_FM  //
/////////////////

	CClust_FM::CClust_FM (double dM, double *pdZ, double *pdM)
		: CClust_F (dM, pdZ), CClust_M (pdM)
	{

	}

	void CClust_FM::EstimClustParams ()
	{
		t_size k ;

		SMatD mXc (m_aTemp[0], m_mX.dim ()) ;
		const double *pClustSize = m_vClustSize ;

		for (k = m_K - 1; k != NAI; k--)
		{
			const double &dCurClustSize = pClustSize [k] ;

			if (dCurClustSize <= m_dZeroTol)
				continue ;

			SVecD vDCurMean = m_mCurM.GetColRef (k) ;	//	2do: use m_avCurM [k]	-> Implement	SCVecArray

				//R	iter$center[k,] = (t(iter$z_ij[,k]) %*% X) / iter$csize[k]
			vDCurMean.Reset (0) ;
			EO<SOP::ApaBmC>::VtMcVc_NC (*vDCurMean, m_mX, m_mZ.GetColRef (k)) ;
			EO<SOP::a_divide>::VSc (*vDCurMean, dCurClustSize) ;
		}
	}

	void CClust_FM::SaveCurResult ()
	{
		CClust_F::SaveCurResult () ;
		CClust_M::SaveCurResult () ;
	}
