#include "tclust.h"

	CTClust::CTClust (int *pdwParamIn, double *pdParamIn, double *pdS, double *pdCSize)
		: m_n (pdwParamIn[0]), m_p (pdwParamIn[1]), m_K (pdwParamIn[2]), m_nDeter (pdwParamIn[6]), m_dwIterTune1 (pdwParamIn[2]), m_dwIterTune2 (pdwParamIn[3]), m_dwIterTune3 (pdwParamIn[4]), m_dwOVV (pdwParamIn[5])
//		, m_dwTrace (2)
		, m_nConvCount (m_nFuzzy), m_nIterSuccess (m_nFuzzy), m_nCode (m_nFuzzy), m_nErrExc (m_nFuzzy), m_dwCountRestrOk(m_nFuzzy)
		, m_pnConvER (NULL)

		, m_dAlpha (pdParamIn[0]), m_dRestrFactor (pdParamIn[1])

		, m_dM (0)
		, m_dZeroTol (pdParamIn[2])
		, m_dBestObj (pdParamIn[0]), m_dUnRestrFactBest (pdParamIn[0])
		, m_pdObjER (NULL)

		, m_dDensFact (pow (2 * meal_PI (), m_p / -2.0))
		, m_dPInv (1.0 / m_p)
		, m_dRestrFactp1p (pow (m_dRestrFactor, m_dPInv))
		, m_dZeroTolSqrt (sqrt(m_dZeroTol))
		, m_dMm1Inv (1 / (m_dM - 1))
		, m_dwNoTrim ((t_size) floor (m_n * (1-m_dAlpha)))
		, m_dwTrim (m_n-m_dwNoTrim)

		, m_vClustSize (pdCSize, m_K)
		, m_amCurS (pdS, m_p, m_p, m_K)
	{
		meal_GetRNGstate () ;
                //VT::22.03.2023
                //  meal_printf("\nMY-TRACE ... CTClust() constructor 1 ...\n");
                //  meal_printf("\n%d %d %d \n", m_n, m_p, m_K);
	}

	CTClust::CTClust (int *pdwParamIn, double *pdParamIn, double *pdS, double *pdCSize, double *pdEVec, double *pdEval)
		: m_n (pdwParamIn[0]), m_p (pdwParamIn[1]), m_K (pdwParamIn[2]), m_nDeter (pdwParamIn[6]), m_dwIterTune1 (pdwParamIn[2]), m_dwIterTune2 (pdwParamIn[3]), m_dwIterTune3 (pdwParamIn[4]), m_dwOVV (pdwParamIn[5])
		, m_nConvCount (m_nFuzzy), m_nIterSuccess (m_nFuzzy), m_nCode (m_nFuzzy), m_nErrExc (m_nFuzzy), m_dwCountRestrOk (m_nFuzzy)
		, m_pnConvER (NULL)

		, m_dAlpha (pdParamIn[0]), m_dRestrFactor (pdParamIn[1])

		, m_dM (0)
		, m_dZeroTol (pdParamIn[2])
		, m_dBestObj (pdParamIn[0]), m_dUnRestrFactBest (pdParamIn[0])
		, m_pdObjER (NULL)

		, m_dDensFact (pow (2 * meal_PI (), m_p / -2.0))
		, m_dPInv (1.0 / m_p)
		, m_dRestrFactp1p (pow (m_dRestrFactor, m_dPInv))
		, m_dZeroTolSqrt (sqrt(m_dZeroTol))
		, m_dMm1Inv (1 / (m_dM - 1))
		, m_dwNoTrim ((t_size) floor (m_n * (1-m_dAlpha)))
		, m_dwTrim (m_n-m_dwNoTrim)

		, m_vClustSize (pdCSize, m_K)
		, m_mEVal (m_p, m_K)
		, m_amEVec (m_p, m_p, m_K), m_amCurS (pdS, m_p, m_p, m_K)
	{
		meal_GetRNGstate () ;
                //VT::22.03.2023
                //  meal_printf("\nMY-TRACE ... CTClust() constructor 2 ...\n");
                //  meal_printf("\n%d %d %d \n", m_n, m_p, m_K);
	}

	CTClust::CTClust (int *pnParIn, int *pnParOut, double *pdParIn, double *pdParOut, double *pdX, double *pdM, double *pdS, int *pnAssign, double *pdClustSize, double *pdWeights, double *pdZ, double *pdObjER, int *pnConvER)
				//	parameters
		: m_n (pnParIn[0]), m_p (pnParIn[1]), m_K (pnParIn[2])
		, m_nFuzzy (pnParIn[3]), m_nIter (pnParIn[4]), m_nKSteps (pnParIn[5])
		, m_nEqualWeights (pnParIn[6]), m_nRestr (pnParIn[7]), m_nDeter (pnParIn[8]), m_nTrace (pnParIn[9])
		, m_dwIterTune1 (pnParIn[10]), m_dwIterTune2 (pnParIn[11]), m_dwIterTune3 (pnParIn[12])
		, m_dwOVV (pnParIn[13])
		, m_nConvCount (pnParOut[0] = 0), m_nIterSuccess (pnParOut[1] = 0), m_nCode (pnParOut[2]), m_nErrExc (pnParOut[3]), m_dwCountRestrOk (pnParOut[4] = 0)
		, m_pnConvER (pnConvER)

		, m_dAlpha (pdParIn[0]), m_dRestrFactor (pdParIn[1]), m_dM (pdParIn[2]), m_dZeroTol (pdParIn[3])
		, m_dBestObj (pdParOut[0]), m_dUnRestrFactBest (pdParOut[1])
		, m_pdObjER (pdObjER)

			//	constants
		, m_dDensFact (pow (2 * meal_PI (), m_p / -2.0))
		, m_dPInv (1.0 / m_p)
		, m_dRestrFactp1p (pow (m_dRestrFactor, m_dPInv))
		, m_dZeroTolSqrt (sqrt(m_dZeroTol))
		, m_dMm1Inv (1 / (m_dM - 1))
		, m_dwNoTrim ((t_size) floor (m_n * (1-m_dAlpha)))
		, m_dwTrim (m_n-m_dwNoTrim )

			//	initializing vectors...
		, m_vInd (m_n), m_vIndBest (pnAssign, m_n), m_vIndOld (m_n), m_vCurInd (m_p + 1), m_vRank (m_n)
		, m_vWeights (m_K), m_vBestWeights (pdWeights, m_K), m_vClustSize (m_K), m_vClustSizeBest (pdClustSize, m_K), m_vDisc (m_n)//, m_vDiscSorted (m_n)

			//	initializing matrices...
		, m_mCurM (m_p, m_K), m_mBestM (pdM, m_p, m_K), m_mX (pdX, m_n, m_p), m_mLL (m_n, m_K), m_mEVal (m_p, m_K), m_mZ (m_n, m_K), m_mZOld (m_n, m_K)

			//	initializing matrix arrays...
		, m_amEVec (m_p, m_p, m_K), m_amCurS (m_p, m_p, m_K),m_amBestS (pdS, m_p, m_p, m_K) //, m_amBestEVec (m_p, m_p, m_K), 
	{
		if (m_nFuzzy)
			m_mZ_best.Set (pdZ, m_n, m_K) ;

//			if (m_dwTrace >= 1)
//			{
//				Rprintf ("nParams: %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n", m_n, m_p, m_K, m_dwIter, m_dwKSteps, m_dwEqualWeights, m_dwRestr, m_dwTrace, m_dwConvCount, m_dwIterSuccess, m_dwCode) ;
//				Rprintf ("dParams: %f, %f, %f, %f\n", m_dAlpha, m_dFactor, m_dZeroTol, m_dBestObj) ;
//
//			}

			//	splitting matrices and tensors into vector - and matrix arrays
		meal_GetRNGstate () ;

                //VT::22.03.2023
                //  meal_printf("\nMY-TRACE ... CTClust() constructor 3 ...\n");
                //  meal_printf("\n%d %d %d \n", m_n, m_p, m_K);

		calc () ;
	}

	CTClust::~CTClust ()
	{
		meal_PutRNGstate () ;
	}

	void CTClust::calc ()
	{	
		if (!CheckParams ())
		{
			m_nErrExc = 1 ;
			return ;
		}
		int &i = m_nIterSuccess = 0 , j ;

		m_vIndOld.Reset (NAI) ;

		double dLastObj = 0 ;

		CheckRestrictions () ;

		for (i = 0; i < m_nIter; i++)
		{
            //VT::08.05.2018
            // meal_printf("MY-TRACE ... Entering FindInitClustAssignment() from CTClust::calc() ...\n");

			FindInitClustAssignment () ;
			FindInitClustSize_R () ;
			//FindInitClustSize () ;

			for (j = 0; 1; j++)
			{
			     // meal_printf ("\n---- Converged vs. successful: %d in %d (%d/%d)\n", m_nConvCount, m_nIterSuccess, i, j) ;

				if (!restrEval ())
				{
					if (j)	//	2do: BUG: shouldn't this be "j" !?!? -> 20120515: changed to j
					{
						SaveCurResult (CalcObjFunc (), 2) ;
						return ;
					}
					else
						SetAllCovmatsIdent () ;
				}

				if (!FindClustAssignment () ||		//	returns false, if the cluster assignment has not changed -> break 
					j == m_nKSteps)					//	max number of iterations reached -> break
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

			m_pnConvER[i] = j < m_nKSteps ;
			if (m_pnConvER[i])						//	did this iteration converge?
				m_nConvCount ++ ;

			double dCurObj = CalcObjFunc () ;
			m_pdObjER[i] = dCurObj ;

			CheckRestrOk () ;

			if (!i || dCurObj > m_dBestObj)			//	store if the current solution is the best solution so far..
				SaveCurResult (dCurObj, j >= m_nKSteps) ;
		
			// meal_printf ("\n++++ Converged vs. successful: %d in %d (%d/%d)\n", m_nConvCount, m_nIterSuccess, i, j) ;
        }
	}

	void CTClust::CheckRestrOk ()
	{
		if (m_dUnRestrFact < m_dRestrFactor)
			++m_dwCountRestrOk ;
	}

	BOOL CTClust::CheckParams ()
	{		//	checks the input parameters for validity

		if (m_nRestr < 0 || m_nRestr > 4)
		{
			meal_printf ("Input parameter error: The restriction type must be between 0 and 4\n") ;
			return FALSE ;
		}

		if (m_dRestrFactor < 1)
		{
			meal_printf ("Input parameter error: The restriction factor must be >= 1\n") ;
			return FALSE ;
		}
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
		if (m_nFuzzy && m_dM < 1)
		{
			meal_printf ("Input parameter error: m must  be >= 1\n") ;
			return FALSE ;
		}

		if (!m_dwNoTrim || m_dAlpha >= 1.0)
		{
			meal_printf ("Input parameter error: alpha was chosen too small (all observations were trimmed)\n") ;
			return FALSE ;
		}
		return TRUE ;
	}


	void CTClust::FindInitClustAssignment ()
	{		//	finds initial cluster assignment
		t_size k ;

		SVecN vTempN (m_aTemp [0], m_n) ;
		SMatD mDCurCluster (m_aTemp [0], m_p + 1, m_p) ;
		SVecN vNCurIdx (m_aTemp [1], m_p + 1) ;

		const double dCorrFact = (m_p) / (m_p  + 1.0) ;
		for (k = 0; k < m_K; k++)
		{								//	for all clusters
										//	finds p+1 observations for forming the initial cluster
			SampleNoReplace (m_p + 1, m_mX.nrow (), vNCurIdx, vTempN) ;
										//	
			SVecD vDCurMean = m_mCurM.GetColRef (k) ;	//	2do: use m_avCurM [k]	-> Implement	SCVecArray
			vDCurMean.Reset (0) ;
										//	stores the subsetted X matrix into MDCurCluster, and sums up it's columns to vcDCurMean.
			EO<UOP::AaC_BpaC>::MsVetMcdVcei (!mDCurCluster, *vDCurMean, m_mX, vNCurIdx) ;

										//	Divide the colSums of MDCurCluster by it's number of rows. -> vDCurMean is a mean vector.
			EO<SOP::a_divide>::VSc (*vDCurMean, mDCurCluster.nrow ()) ;

										//	Centers the mDCurCluster - matrix
			EO<SOP::a_minus>::MVcet (!mDCurCluster, vDCurMean) ;

										//	Calculating the covariance matrix
			cov_centered_NC (!m_amCurS[k], mDCurCluster, dCorrFact) ;
		}
	}

	void CTClust::FindInitClustSize_R ()
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

	void CTClust::FindInitClustSize ()
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

	void CTClust::SaveCurResult (double dCurObj, int nCode)
	{
		m_nCode = nCode ;

		m_dBestObj = dCurObj ;

		m_vClustSizeBest.Copy_NC (m_vClustSize) ;
		m_vBestWeights.Copy_NC (m_vWeights) ;

		m_vIndBest.Copy_NC (m_vInd) ;

		m_mBestM.Copy_NC (m_mCurM) ;
		m_dUnRestrFactBest = m_dUnRestrFact ;

		if (m_nFuzzy)
			m_mZ_best.Copy_NC (m_mZ) ;

		t_size k ;
		for (k = m_K - 1; k != NAI; k--)
			m_amBestS[k].Copy_NC (m_amCurS[k]) ;	//	2do: use tensor!
	}

	void CTClust::SetAllCovmatsIdent ()
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

	void CTClust::LoadCluster (SMatD &c, t_size k)
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

	double CTClust::CalcObjFunc ()	//	temp3
	{
		if (m_nFuzzy)
			return CalcObjFunc_fuzzy () ;
		else
			return CalcObjFunc_cat () ;
	}

	double CTClust::CalcObjFunc_cat ()
	{
		double dObj = 0 ;

		t_size k ;

		SVecD vDensity (m_aTemp[3], 0) ;
		SMatD mCurX (m_aTemp[4], m_n, m_p) ;

		double *pdClustSize = m_vClustSize ;

		for (k = m_K - 1; k != NAI; --k)
		//for (k = 0; k < m_K; ++k)
		{
			LoadCluster (mCurX, k) ;

			t_size dwClustSize = mCurX.nrow () ;
			double &dClustSize = pdClustSize [k] = dwClustSize ;		//	2do: is this assignment necessary? 

			if (dClustSize <= m_dZeroTol)
				continue ;

			vDensity.Require (dwClustSize) ;

			CalcDensity (mCurX, vDensity, m_mCurM.GetColRef (k), m_mEVal.GetColRef (k), m_amEVec [k]) ;	//	2do: introduce SCVecArray

			EO<UOP::Apa_logB>::SVc (dObj, vDensity) ;

			if (!m_nEqualWeights)	//	this can be avoided when equal-sized clusters are expected (only adding a fixed value to each object. function)
				dObj += dwClustSize * log (dwClustSize / (double) m_dwNoTrim) ;
		}

		return dObj ;
	}

	double CTClust::CalcObjFunc_fuzzy ()
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

			CalcDensity (m_mX, vDensity, m_mCurM.GetColRef (k), m_mEVal.GetColRef (k), m_amEVec [k]) ;	//	2do: introduce SCVecArray

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

	void CTClust::CalcDensity (const SCMatD &mX, const SVecD &vDens, const SVecD &vCurM, const SCVecD &vEVal, const SCMatD &mEVec, const double dFact)
	{
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

	BOOL CTClust::FindClustAssignment ()
	{
		if (m_nFuzzy)
			return FindClustAssignment_fuzzy () ;
		return FindClustAssignment_cat () ;		
	}

	BOOL CTClust::FindClustAssignment_cat ()
	{		//	calculates the vector m_vIndTrim
		t_size k ;
		for (k = m_K - 1; k != NAI; k--)
			CalcDensity (m_mX, m_mLL.GetColRef (k), m_mCurM.GetColRef (k), m_mEVal.GetColRef (k), m_amEVec [k], m_vWeights (k)) ;		//	2do: Implement VecArray

		SVecD vDisc (m_aTemp[3], m_n) ;

		FindNearestClust (vDisc, m_vInd) ;

		FindOutliers (vDisc, m_vInd) ;

		if (equal (m_vInd, m_vIndOld))				//	check the change of group assignment AFTER the calculation of outliers. 
			return FALSE ;							//	otherwise the alg will stop after the 1st iteration when k == 1, 
		m_vIndOld.Copy_NC (m_vInd) ;				//	since there wouldn't even be the change of the change of the cluster assignment

		CalcClusterSize_cat () ;					//	Cluster sizes are always calculated. Just don't calc the weigths for (m_dwEqualWeights)

		return TRUE ;
	}

	void CTClust::CalcClusterSize_cat ()
	{
		m_vClustSize.Reset (0) ;

		double *pdClustSize = m_vClustSize ;
		int *pnInd = m_vInd ;
		int * const pnEndInd = m_vInd.GetDataEnd () ;

		for (; pnInd < pnEndInd; ++pnInd)
			if (*pnInd != -1)		//	not an outlier
				pdClustSize[*pnInd] += 1 ;

		if (!m_nEqualWeights)
			EO<SOP::divide_r>::VScVc (*m_vWeights, m_dwNoTrim, m_vClustSize) ;	//	m_vWeights <- m_vClustSize / m_dwNoTrim
	}

	BOOL CTClust::FindClustAssignment_fuzzy ()
	{
		t_size k ;
		for (k = m_K - 1; k != NAI; k--)
			CalcDensity (m_mX, m_mLL.GetColRef (k), m_mCurM.GetColRef (k), m_mEVal.GetColRef (k), m_amEVec [k], m_vWeights (k)) ;		//	2do: Implement VecArray

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

		if (!m_nEqualWeights)
			EO<SOP::divide_r>::VScVc (*m_vWeights, sum (m_vClustSize), m_vClustSize) ;	//	calculate cluster weights (if necessary)

		return TRUE ;
	}

	void CTClust::CalcFuzzyRow (const SCVecD &ll, const SVecD &z, double &dDisc, int &nInd)	//	temp1
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

	void CTClust::SetCatZ (const SVecD &z, int nIdx)
	{
		double *pdZ = z, *const pdEndZ = z.GetDataEnd () ;
		
		while (pdZ < pdEndZ)
		{
			*pdZ = !nIdx ;
			--nIdx ;
			++pdZ ;
		}
	}

	void CTClust::FindOutliers (const SVecD &vDisc, const SVecN &vInd)
	{
		FindOutliers_old (vDisc, vInd) ;
	}

	void CTClust::FindOutliers_old (const SVecD &vDisc, const SVecN &vInd)
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

	void CTClust::FindOutliers_new (const SVecD &vDisc, const SVecN &vInd)
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

	void CTClust::FindNearestClust (const SVecD &vDisc, const SVecN &vInd)
	{
		FindNearestClust_old (vDisc, vInd) ;
	}

	void CTClust::FindNearestClust_old (const SVecD &vDisc, const SVecN &vInd)	//	old version - without tie-handling!
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

	void CTClust::FindNearestClust_new (const SVecD &vDisc, const SVecN &vInd)	//	new version - implements tie-handling!
	{
		t_size r ;

		double *pdCurDisc = vDisc ;
		int *pnCurInd = vInd ;

		SVecD curRow (m_aTemp[1], m_K) ;

		for (r = 0; r < vDisc.size(); ++r)									//	for each row of m_mLL
		{
			CopyRow (*curRow, m_mLL, r) ;
//			select_cluster_old (*pdCurDisc, *pnCurInd, curRow) ;			//	select the cluster corresponding to the largest value of the current row of m_mLL
			select_cluster (*pdCurDisc, *pnCurInd, curRow) ;				//	select the cluster corresponding to the largest value of the current row of m_mLL
			++pdCurDisc ;
			++pnCurInd ;
		}
	}

	void CTClust::select_cluster (double &dDisc, int &nInd, const SCVecD &row)		//	temp0
	{
		select_cluster_old (dDisc, nInd, row) ;
	}

	void CTClust::select_cluster_old (double &dDisc, int &nInd, const SCVecD &row)	//	no ties - handling
	{
		double const *pdRow = row, *pdMax = pdRow, *const pdEndRow = row.GetDataEnd ()  ;

		for (++pdRow; pdRow < pdEndRow; ++pdRow)
			if (*pdMax < *pdRow)
				pdMax = pdRow ;
		nInd = pdMax - row ;
		dDisc = *pdMax ;
	}

	void CTClust::select_cluster_new (double &dDisc, int &nInd, const SCVecD &row)	//	temp0
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


	void CTClust::EstimClustParams ()
	{
		if (m_nFuzzy)
			EstimClustParams_fuzzy () ;
		else
			EstimClustParams_cat () ;
	}

	void CTClust::EstimClustParams_cat ()
	{		//	calculates the mean and the cov matrix of the clusters currently defined by m_vIndTrim
		t_size k ;

		SMatD mCurX (m_aTemp[0]) ;						//	2do: eval max clustersize and require this matrix.

		for (k = m_K - 1; k != NAI; k--)
//		for (k = 0; k < m_K; k++)
		{
			t_size dwCurClustSize = (t_size) m_vClustSize (k) ;

			SVecD vDCurMean = m_mCurM.GetColRef (k) ;	//	2do: use m_avCurM [k]	-> Implement	SCVecArray

			if (dwCurClustSize > m_dZeroTol)
			{
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
	}

	void CTClust::EstimClustParams_fuzzy ()
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
