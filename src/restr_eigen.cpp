
#include "math.h"
#include "tclust.h"

	BOOL CheckClusterSingularity (const SCMatD &mEV, const SCVecD & vClustSize, double dZeroTol) ;
	BOOL CheckClusterSingularity_NC (const SCMatD &mEV, const SCVecD & vClustSize, double dZeroTol) ;
	double CalcDiff_log (const SCMatD &ev, const SCVecD &ks, const double &dLower, const double &dUpper) ;
	double CalcDiff_log_NC (const SCMatD &ev, const SCVecD &ks, const double &dLower, const double &dUpper) ;
	void ZeroGroupsGetMeanEigenvalues (const SMatD &mEV, const SCVecD & vClustSize, const double &dZeroTol) ;
	void GetCheckArray (const SCMatD &mEV, const SCVecD & vClustSize, double dFact, SVecD &vdCheckEV, double dMax) ;

	double CalcDiff_log (const SCMatD &ev, const SCVecD &ks, const double &dLower, const double &dUpper)
	{
		THROW (ev.ncol () == ks.size ()) ;

//if (!(dLower <= dUpper))
//	meal_printf  ("dLower: %d\tdUpper: %d\n", dLower, dUpper) ;
//		THROW (dLower <= dUpper) ;
		return CalcDiff_log_NC (ev, ks, dLower, dUpper) ;
	}

	double CalcDiff_log_NC (const SCMatD &ev, const SCVecD &ks, const double &dLower, const double &dUpper)
	{
		ASSERT (ev.ncol () == ks.size ()) ;
		ASSERT (dLower <= dUpper) ;

		const double dLogLower = log (dLower) ;
		const double dLogUpper = log (dUpper) ;

		double dRetC = 0, dRet = 0 ;

		//const t_size rm1 = ev.nrow () - 1 ;
		t_size r ;

		const double * pdEv = ev ;
		const double * const pdEndEv = ev.GetDataEnd () ;
		const double * pdKs = ks ;

		while (pdEv < pdEndEv)									//	iterates the clusters (coloumns)
		{
			dRetC = 0 ;
			for (r = ev.nrow (); r > 0; r--)					//	iterates the EVs of the current cluster
			{
				const double &dCur = *pdEv ;

				if (dCur < dLower)
					dRetC += dLogLower + dCur / dLower ;
				else if (dCur > dUpper)
					dRetC += dLogUpper + dCur / dUpper ;
				else
					dRetC += log (dCur) + 1 ;
				++ pdEv ;
			}

			dRet += *pdKs * dRetC ;
			pdKs++ ;
		}

		return dRet ;
	}

	void CalcRST (const SCVecD &vec, double dL, double dU, t_size &dwR, double &dS, double &dT)
	{
		dwR = 0 ;
		dS = dT = 0 ;

		const double *pdVec = vec ;
		const double * const pdEndVec = vec.GetDataEnd () ;

		while  (pdVec < pdEndVec)
		{
			const double &dCur = *pdVec ;

			if (dCur < dL)
				dS += dCur ;
			if (dCur > dU)
				dT += dCur ;
			if (dCur < dL || dCur > dU)
				++dwR ;
			++pdVec ;
		}
	}

	void ZeroGroupsGetMeanEigenvalues (const SMatD &mEV, const SCVecD & vClustSize, const double &dZeroTol)
	{
			// --> just care about groups without observations. set their eigenvalues to the mean of all other EVs

		t_size k, dwCount = 0 ;
		double dSum = 0 ;
		const double *pdClustSize = vClustSize ;

		for (k = vClustSize.size () - 1; k != NAI; k--)
			if (pdClustSize [k] > dZeroTol)
			{
				dwCount ++ ;
				dSum += sum (mEV.GetColRef (k)) ;
			}

		dSum /= dwCount * mEV.nrow () ;

		for (k = vClustSize.size () - 1; k != NAI; k--)
			if (pdClustSize [k] <= dZeroTol)
				mEV.GetColRef (k).Reset (dSum) ;
	}

	
	void GetCheckArray (const SCMatD &mEV, const SCVecD & vClustSize, double dFact, SVecD &vdCheckEV, double dMax)
	{
		vdCheckEV.Require (mEV.size () * 2 + 2) ;

		double *pdCheckEv = vdCheckEV ;

		*pdCheckEv = 0 ;
		++pdCheckEv ;
		*pdCheckEv = dMax ;
		++pdCheckEv ;

		const double *pdEv = mEV ;
		const double *const pdEndEv = mEV.GetDataEnd () ;

		while (pdEv < pdEndEv)
		{
			*pdCheckEv = *pdEv ;
			++pdCheckEv ;
			*pdCheckEv = *pdEv / dFact;
			++pdCheckEv ;
			++pdEv ;
		}

		sort (*vdCheckEV) ;
			//	calcs the average of all the neigbours of vdCheckEV

		pdCheckEv = vdCheckEV ;

		t_size c ;
		for (c = 1; c < vdCheckEV.size (); c++)
			pdCheckEv[c-1] = (pdCheckEv[c-1] + pdCheckEv[c]) / 2 ;
		vdCheckEV.Reshape_NC (vdCheckEV.size () - 1) ;
	}

	BOOL CheckClusterSingularity (const SCMatD &mEV, const SCVecD & vClustSize, double dZeroTol)
	{
		THROW (mEV.ncol () == vClustSize.size ()) ;
		return CheckClusterSingularity_NC (mEV, vClustSize, dZeroTol) ;
	}

	BOOL CheckClusterSingularity_NC (const SCMatD &mEv, const SCVecD & vClustSize, double dZeroTol)
	{
		ASSERT (mEv.ncol () == vClustSize.size ()) ;

		t_size k ;
		const t_size dwColInc = mEv.GetColInc () ;
		const double *pdClustSize = vClustSize ;

		const double *const pdEv = mEv ;
		const double *pdCurEv, *pdEndCurEv ;

		for (k = mEv.ncol () - 1; k != NAI; k--)
			if (pdClustSize [k] > dZeroTol)
			{
				pdCurEv = pdEv + dwColInc * k ;
				pdEndCurEv = pdCurEv + dwColInc ;
				while (pdCurEv < pdEndCurEv)
				{
					if (*pdCurEv > dZeroTol)
						return TRUE ;
					++pdCurEv ;
				}
			}
		return FALSE ;
	}

	BOOL RestrictEigenValues (const SVMatD &mEV, const SCVecD & vClustSize, double dFact, double dZeroTol, double &dUnRestrFact)
	{
		//SVecD &vTempNPp2, 

		if (!CheckClusterSingularity_NC (mEV, vClustSize, dZeroTol))
			return FALSE ;	//	all eigenvalues of all clusters with at least one observation are 0

		THROW (dFact >= 1) ;

		double	dMin = 0, dMax = 0 ;
		t_size i ;

//		mEV.MinMax (dMin, dMax) ;	//	bug: only clusters with vClustSize > 0 shall be considered here...

		BOOL bFoundOne = FALSE ;	//	XXXC 20101108	calculating the min / max of Clustersizes only for clusters with a size > 0
		for (i = mEV.ncol () - 1; i != NAI; i--)
			if (vClustSize (i) > dZeroTol)
			{
				minmax (mEV.GetColRef (i), dMin, dMax, !bFoundOne) ;
				bFoundOne = TRUE ;
			}

		dUnRestrFact = dMax / dMin ;

		const t_size dwKm1 = vClustSize.size () - 1;

		if (//dMin <= dZeroTol ||
			dUnRestrFact > dFact)	//	min (ev) < max (ev) / dFact --> so we have to restrict the eigenvalues
		{
			//IVecD &vdCheckEV = vTempNPp2 ;
			//SVecD &vdCheckEV = SVecD::TempFree_NC (0, mEV.size () + 2) ;
			SVecD vdCheckEV (tempRef (0), mEV.size () + 2) ;

			GetCheckArray (mEV, vClustSize, dFact, vdCheckEV, dMax) ;

			double dMinVal = 0, dMinM = 0 ;

			double dSumU = 0, dSumL = 0 ;

			t_size dwR, c, r, dwCur = vdCheckEV.size () ;
			double dS, dT ;

			const double * const pdClustSize = vClustSize ;

			for (c = 0; c < dwCur; c++)
			//for (c = vdCheckEV.size () - 1; c != NAI; c--)
			{
				dSumU = dSumL = 0 ;
				double &dCurL = vdCheckEV (c), dCurU = dCurL * dFact ;

				for (r = dwKm1; r != NAI; r--)
				{
					CalcRST (mEV.GetColRef (r), dCurL, dCurU,dwR, dS, dT) ; 

					dSumU += pdClustSize [r] * (dS + dT / dFact) ;
					dSumL += pdClustSize [r] * dwR ;
				}

				dSumU /= dSumL ;

				double dCurVal = CalcDiff_log (mEV, vClustSize, dSumU, dSumU * dFact) ;

				if (!c || dMinVal > dCurVal)
				{
					dMinVal = dCurVal ;
					dMinM = dSumU ;
				}
			}

			limit_NC (*mEV, dMinM, dMinM * dFact) ;
		}
		else	//	we don't have to restrict eigenvalues.
			ZeroGroupsGetMeanEigenvalues (mEV, vClustSize, dZeroTol) ;

		return CheckClusterSingularity (mEV, vClustSize, dZeroTol) ;
	}

	void HandleEVsingularities (const SMatD &mEv, double dZeroTol)
	{
		double dColProd, dCurMax, dCurMin = 0 ;

		const double m_dpInv = 1.0 / mEv.nrow () ;			//	-> to CTClust

		double *pdCur, *pdCurCol, *pdCurColEnd ;
		const double * const pEndEv = mEv.GetDataEnd () ;

																		//	calculating the min of the determinant vector, setting all values < dZeroTol to dZeroTol
		//for (i = k - 1; i != NAI; i--)								//	for all coloumns (clusters)

		for (pdCurCol = mEv; pdCurCol < pEndEv; pdCurCol = pdCurColEnd)
		{
			dCurMin = *pdCurCol ;
			dCurMax = *pdCurCol ;

			pdCurColEnd = pdCurCol + mEv.GetColInc () ;

			//oldC	for (j = 0; j < p; j++)							//	for each eigenvalue of the current cluster
			for (pdCur = pdCurCol; pdCur < pdCurColEnd; ++pdCur)
			{
				sm_setmax (dCurMax, *pdCur) ;
				if (*pdCur <= dZeroTol)								//	cur value is smaller than zero tol
					*pdCur = dCurMin = dZeroTol ;
				else if (*pdCur < dCurMin)							//	or smaller than the current minimum (or the first value added)
					dCurMin = *pdCur ;
			}

			dColProd = 1 ;
			if (dCurMin / dCurMax <= dZeroTol)
			{
				dCurMin /= dZeroTol ;

				//oldC	curCol.Limit_U (dCurMin / dZeroTol) ;			//	was dCurMin * 1e15
				for (pdCur = pdCurCol; pdCur < pdCurColEnd; ++pdCur)
				{
					if (*pdCur > dCurMin)
						*pdCur = dCurMin ;
					dColProd *= *pdCur ;
				}
			}
			else
				for (pdCur = pdCurCol; pdCur < pdCurColEnd; ++pdCur)	//	2do: use prod - function..
					dColProd *= *pdCur ;

			//oldC	curCol /= pow (prod (curCol), dpInv) ;
			dColProd = pow (dColProd, -m_dpInv) ;

			for (pdCur = pdCurCol; pdCur < pdCurColEnd; ++pdCur)
				*pdCur *= dColProd ;
		}
	}

	void DeterMinMax (const SCVecD &vDeter, const SCVecD &vClustSize, double &dMin, double &dMax, double dZeroTol)
	{
		double const * const pdClustSize = vClustSize ;
		double const * const pdDeter = vDeter ;

#ifdef NOWARNINGS
		dMin = dMax = 0 ;
#endif

		t_size k ;

		BOOL bFoundOne = FALSE ;
		for (k = vDeter.size () - 1; k != NAI; k--)
			if (pdClustSize [k] > dZeroTol)
			{
				if (!bFoundOne)
				{
					dMin = dMax = pdDeter[k] ;
					bFoundOne = TRUE ;
				}
				else
					SOP::a_minmax::Calc (dMin, dMax, pdDeter[k]) ;
			}
	}

	BOOL RestrictEigenValues_deter (const SVMatD &mEv, const SCVecD & vClustSize, double dFact, double dZeroTol, double &dUnRestrFact)
	{
		t_size p = mEv.nrow (), k = mEv.ncol () ;

		const double m_dpInv = 1.0/p ;			//	-> to CTClust

		//static CDataREF_Static m_aTemp [1] ;	//	-> to CTClust
		//SMatD mDeter (m_aTemp[0], 1, k) ;

		//SMatD &mDeter = SMatD::TempFree_NC (1, 1, k) ;
		SMatD mDeter (tempRef (1), 1, k) ;
		SVecD vDeter (*mDeter, k) ;

		colProds (*vDeter, mEv) ;

		if (!CheckClusterSingularity_NC (mDeter, vClustSize, dZeroTol))
		{
			mEv.Reset (0) ;
			return FALSE ;
		}

		HandleEVsingularities (mEv, dZeroTol) ;
		//double dMin = 0, dMax = 0 ;
		double dMin, dMax ;
		//minmax (vDeter, dMin, dMax) ;
		DeterMinMax (vDeter, vClustSize, dMin, dMax, dZeroTol) ;

		dUnRestrFact = dMax / dMin ;

		//mDeter ^= m_dpInv ;
		if (dUnRestrFact <= dFact)
		{
			ZeroGroupsGetMeanEigenvalues (mDeter, vClustSize, dZeroTol) ;
			EO<SOP::a_pow>::VSc (*mDeter, m_dpInv) ;
		}
		else
		{
			EO<SOP::a_pow>::VSc (*mDeter, m_dpInv) ;
			RestrictEigenValues (!mDeter, vClustSize, pow (dFact, m_dpInv), dZeroTol, dUnRestrFact) ;
		}

		EO<SOP::a_multiply>::MVcet (!mEv, vDeter) ;
		//mEv.byrow () *= vDeter ;

		return CheckClusterSingularity (mEv, vClustSize, dZeroTol) ;
	}

