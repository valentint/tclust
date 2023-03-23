#include "tclust.h"

	BOOL CTClust::restrEval ()
	{
        //VT::09.08.2017
	// we have a problem in restr_equal, i.e. restriction on Sigma.
	// the last two restrictions are not used.
        if(m_nTrace > 0)
            meal_printf ("TRACE... Restriction: %d\n", m_nRestr) ;

		switch (m_nRestr)
		{
		case 0:	return restr_diff_ax() ;     // eigen
		case 1: return restr_dir() ;         // deter
		case 2:	return restr_equal() ;       // sigma
		case 3: return restr_prop() ;        // suppressed
		case 4: return restr_none();         // suppressed
		}
		THROW (0) ;
		return FALSE ;
	}

	BOOL CTClust::restr_none ()
	{
		t_size k ;
		for (k = m_K - 1; k != NAI; --k)
			//eigen_sqr (m_amCurS[k], m_avEVal[k], m_amEVec[k], FALSE) ;
			sme_eigen_sqr_NC (!m_amCurS[k], m_mEVal.GetColRef (k), !m_amEVec[k], FALSE) ;			//	2do: use SCVecArray instead!

		double dMin = 0, dMax = 0; 
		minmax (m_mEVal, dMin, dMax) ;
		return dMin / dMax > m_dZeroTol ;

	}

	BOOL CTClust::restr_diff_ax ()
	{
		t_size k ;
		for (k = m_K - 1; k != NAI; --k)
			//eigen_sqr (m_amCurS[k], m_avEVal[k], m_amEVec[k], FALSE) ;
			sme_eigen_sqr_NC (!m_amCurS[k], m_mEVal.GetColRef (k), !m_amEVec[k], FALSE) ;			//	2do: use SCVecArray instead!

		limit_l<double> (*m_mEVal, 0) ;
//		m_mEVal.Limit_L (0) ;

		if (m_nDeter)
		{
			if (!RestrictEigenValues_deter (!m_mEVal, m_vClustSize, m_dRestrFactor, m_dZeroTol, m_dUnRestrFact))
				return FALSE ;
		}	
		else if (!RestrictEigenValues (!m_mEVal, m_vClustSize, m_dRestrFactor, m_dZeroTol, m_dUnRestrFact))		//	trims the eigenvalues of the used clusters (given by m_vUsedK)
			return FALSE ;

			//	recalculates (recomposes) all covariance matrices

		for (k = 0; k < m_K; k++)
			sme_matmult_a_diagb_at_NC (m_amEVec[k], m_mEVal.GetColRef (k), !m_amCurS [k]) ;		//	2do: use SCVecArray instead!
///		{
//			m_mTempNP1.Reshape (m_p, m_p) ;
//			FC_ElOp<FC::FC_multiply, double>::OpMV_row (m_amEVec[k], m_avEVal[k], m_mTempNP1) ;	//	m_avEVal now is the squareroot of actually trimmed the EV
//			matmultmat (m_mTempNP1, t (m_amEVec[k]), m_amCurS [k]) ;
//		}
		return TRUE ;
	}

    // Restriction 'sigma'
	BOOL CTClust::restr_equal ()
	{
		t_size k ;

        //VT::09.08.2017
        if(m_nTrace > 0)
            meal_printf ("TRACE... Restriction=sigma, m_K: %d\n", m_K) ;

		EO<SOP::a_multiply>::VSc (*m_amCurS [0], m_vClustSize (0) / m_dwNoTrim) ;

		for (k = m_K - 1; k; k--)
			EO<SOP::ApaBmC>::VScVc (*m_amCurS [0], m_vClustSize (k) / m_dwNoTrim, m_amCurS [k]) ;

		sme_eigen_sqr_NC (m_amCurS[0], m_mEVal.GetColRef (0), !m_amEVec [0], FALSE) ;	//	2do: implement VecArray

		limit_l (*m_mEVal, 0) ;

		for (k = m_K - 1; k; k--)
		{
            //VT::09.08.2017
            if(m_nTrace > 0)
                meal_printf ("TRACE... Restriction=sigma, k: %d\n", k) ;

			m_mEVal.GetColRef (k).Copy_NC (m_mEVal.GetColRef (0)) ;

			m_amCurS [k].Copy_NC (m_amCurS [0]) ;
			m_amEVec [k].Copy_NC (m_amEVec [0]) ;
		}

		return max (m_mEVal.GetColRef (0)) > m_dZeroTol ;
		//return m_avEVal[0].max () > 0 ;
	}


	void CTClust::CheckRestrictions ()
	{
		switch (m_nRestr)
		{
		default: THROW (FALSE) ;
			
		case 0:	//	diffax
			if (m_p == 1)								//	1d
			{
				if(m_dRestrFactor == 1)					//		and all variances are supposed to be equal
					m_nRestr = 2 ;						//			use restr_equal
				if (m_nDeter)								//		using the deter-criterion does not make sense with 1d data
					m_nDeter = FALSE ;
			}
			break ;
		case 1:	//	dir
			if (m_p == 1)								//	1d
				if (m_dRestrFactor == 1)				//		and all variances are supposed to be equal
					m_nRestr = 2 ;						//			use restr_equal
				else
				{
					m_nRestr = 1 ;						//		use normal restr_diffax, as in 1d data the direction does not matter.
					m_nDeter = FALSE ;
				}
			else if (m_dRestrFactor == 1 && !m_nDeter)	//	searching for spherical clusters -> the direction does not matter.
				m_nRestr = 1 ;
			break ;
		case 2:	break ;									//	all cov matrices are identical
		case 3:	//	prop
			if (m_dRestrFactor == 1)					//	all cov structures are supposed to be equal
				m_nRestr = 2 ;							//		use restr_equal
			else if (m_p == 1)							//	1d
			{											//		use the simple eigenvalue restriction
				m_nRestr = 1 ;
				m_nDeter = FALSE ;
			}
			break ;
		case 4:
			break ;
		}		
	}

	BOOL CTClust::restr_prop ()
	{
		const SMatD &mDu_best = m_amEVec [0] ;
		SMatD mDu (m_aTemp [8], m_p, m_p) ;

		SVecD vDsf_best (m_aTemp [5], m_K, m_p) ;
		SVecD vDst_best (m_aTemp [5], m_K) ;


		SVecD vDsf (m_aTemp [6], m_K, m_p) ;
		SVecD vDst (m_aTemp [6], m_K) ;

		SMatD mDdd (m_aTemp [7], m_p, m_K) ;

		SVecD m_vTemp1 (m_aTemp [10], sm_max (m_p, m_K)) ;

		SMatD m_mTemp1 (m_aTemp [11], m_p + 1, m_p) ;
		SMatD m_mTemp2 (m_aTemp [12], m_p, m_p) ;

		double * const pdst = vDst ;
		const double * const pdCluStSize = m_vClustSize ;
		double dMinRes = 0 ;

		t_size i, j, k ;
		for (i = 0; i < m_dwIterTune1; i++)
		{
			m_mTemp1.Reshape_NC (m_p + 1, m_p) ;
			m_mTemp2.Reshape_NC (m_p, m_p) ;

				//R ee <- matrix(runif(pa$p*pa$p, min=0, max=1),pa$p,pa$p)
			//sme_runif (*m_mTemp1) ;
			rnorm (*m_mTemp1) ;

				//R w <- eigen(t(ee)%*%ee)
			sme_tmatmult_NC (m_mTemp1, m_mTemp1, !m_mTemp2, TRUE, FALSE) ;

			m_mTemp2.Reshape_NC (m_p, m_p) ;
			sme_eigen_sqr_NCV (!m_mTemp2, vDsf, !mDu, TRUE) ;

			EO<SOP::a_divide>::VSc (*vDsf, pow (prod (vDsf), 1.0 / m_p)) ;
			//EO<OPA::multiply>::vsc (vDsf, 1 / pow (prod (*vDsf), 1 / m_p)) ;	//	2check: faster?

			for (j = 0; j < m_dwIterTune2; j++)
			{
//				m_vTemp1.Reshape_NC (m_p) ;
//				EO<SOP::inv>::vevc_NC (*m_vTemp1, vDsf) ;

				m_mTemp2.Reshape (m_p, m_p) ;
				for (k = m_K - 1; k != NAI; k--)
				{
					//R	dd[,i] <- diag (t(u) %*%iter$sigma[,,i] %*% u)
//					sme_matmult_NC (mDu, m_amCurS [k], m_mTemp2) ;							//	2do: only calculate diagonal - impr. runtime!
//					sme_matmult_diag_NC (m_mTemp2, mDu, mDdd.GetColRef (k), FALSE, TRUE) ;

					sme_matmult_a_b_at_NC (mDu, m_amCurS [k], !m_mTemp2, TRUE, FALSE) ;
					SVecD vDiag (mDdd.GetColRef (k)) ;
					CopyDiag_NC (vDiag, m_mTemp2) ;

//					pdst[k] = sumprod (m_vTemp1, vDiag) / m_p ;
				}

				//R	st[i]=(sum(dd[,i]/sf)/pa$p)
				vDst.Reset (0) ;
				EO<SOP::ApaBdC>::VtMcVc_NC (*vDst, mDdd, vDsf) ;
				EO<SOP::a_divide>::VSc (*vDst, m_p) ;

				RestrictEigenValues (!SMatD (*vDst, 1, m_K), m_vClustSize, m_dRestrFactp1p, m_dZeroTol, m_dUnRestrFact) ;

				//R for (j in 1:pa$p)
				//R		sf2[j] <- sum (matrix (iter$csize,nrow=1) * ((st^(-1 / pa$p)) * dd[j,]))

				m_vTemp1.Reshape_NC (m_K) ;
				EO<SOP::divide>::VVcVc_NC (*m_vTemp1, m_vClustSize, vDst) ;
				//EO<SOP::divide>::vevcevc_NC (*m_vTemp1, m_vClustSize, vDst) ;
				vDsf.Reset (0) ;
				EO<SOP::ApaBmC>::VMcVct_NC (*vDsf, mDdd, m_vTemp1) ;

				//R sf3 <- prod(sf2)
				//R sf <- sf2 / sf3^(1 / pa$p)
				EO<SOP::a_divide>::VSc (*vDsf, pow (prod (vDsf), 1.0 / m_p)) ;

				//R	for (i in 1:pa$K)
				//R		dd[, i] <- t (sf) * st[i]^(1 / pa$p)
				EO<SOP::multiply>::MVcVct_NC (!mDdd, vDsf, vDst) ;

				OptVectors (mDu, mDdd) ;
			}
			double dRes = 0 ;

			m_vTemp1.Reshape_NC (m_p) ;
			//EO<SOP::inv>::vevc_NC (*m_vTemp1, vDsf) ;
			EO<SOP::inv>::VVc_NC(*m_vTemp1, vDsf) ;

			m_mTemp1.Reshape_NC (m_p, m_p) ;
			sme_matmult_a_diagb_at_NC (mDu, m_vTemp1, !m_mTemp1) ;
			for (k = m_K - 1; k != NAI; k--)
				dRes += pdCluStSize [k] * (log (pow (pdst[k], (double) m_p)) + sme_sum_matmult_diag_NC (m_mTemp1, m_amCurS[k]) / pdst[k]) ;

			if (!i || dRes < dMinRes)
			{
				mDu_best.Copy_NC (mDu) ;
				vDsf_best.Copy_NC (vDsf) ;
				vDst_best.Copy_NC (vDst) ;
				dMinRes = dRes ;
			}
		}

		if (max (vDsf_best) <= m_dZeroTol)
			return FALSE ;

		//R iter$sigma[,,i] <- st_best[i]^(1/pa$p)*u_best%*%diag(sf_best)%*%t(u_best)
			//	calculating the base covariance matrix.
		m_mTemp1.Reshape (m_p, m_p) ;
		sme_matmult_a_diagb_at_NC (mDu_best, vDsf_best, !m_mTemp1) ;

		double * const pdst_best = vDst_best ;
			//	multiply the base covariance matrix with the corresponding factor st[k]
		for (k = m_K - 1; k != NAI; k--)
		{
			EO<SOP::multiply>::VScVc_NC (*m_amCurS [k], pdst_best [k], m_mTemp1) ;
//			EO<SOP::multiply>::VeVcS_NC (*m_amCurS [k], m_mTemp1, pdst_best [k]) ;
			if (k)
				m_amEVec [k].Copy_NC (mDu_best) ;								//	reconstruct eigenvectors
		}
		EO<SOP::multiply>::MVcVct_NC (!m_mEVal, vDsf_best, vDst_best) ;		//	reconstruct eigenvalues

		return TRUE ;
	}


	BOOL CTClust::restr_dir ()
	{
		double dMinRes = 0 ;

		//SMatD mDBestU (m_aTemp [5], m_p, m_p) ;
		const SMatD &mDBestU = m_amEVec[0] ;//(m_aTemp [5], m_p, m_p) ;

		//SMatD mDBestD (m_aTemp [6], m_p, m_K) ;
		SMatD &mDBestD  = m_mEVal ; //(m_aTemp [6], m_p, m_K) ;

		SMatD mDdd (m_aTemp [7], m_p, m_K) ;

		SMatD mDu (m_aTemp [8], m_p, m_p) ;

		t_size i, j, k ;

		SMatD m_mTemp1 (m_aTemp [11], m_p + 1, m_p) ;
		SMatD m_mTemp2 (m_aTemp [12], m_p, m_p) ;
		SVecD m_vTemp1 (m_aTemp [10], m_p) ;

		for (i = 0; i < m_dwIterTune1; i++)
		{
			m_mTemp1.Reshape_NC (m_p + 1, m_p) ;

				//R ee <- matrix(runif(pa$p*pa$p, min=0, max=1),pa$p,pa$p)
			//sme_runif (*m_mTemp1) ;
			rnorm (*m_mTemp1) ;

				//R w <- eigen(t(ee)%*%ee)
			sme_tmatmult_NC (m_mTemp1, m_mTemp1, !m_mTemp2, TRUE, FALSE) ;

			m_mTemp2.Reshape_NC (m_p, m_p) ;
			sme_eigen_sqr_NCV (!m_mTemp2, m_vTemp1, !mDu, TRUE) ;	//	what 4 - are these eigenvectors actually ever used?

			m_mTemp2.Reshape_NC (m_p, m_p) ;
			for (j = 0; j < m_dwIterTune2; j++)
			{
				//R	for (i in 1:pa$K) 
				for (k = m_K - 1; k != NAI; k--)
				{
					//R	dd[,i] <- diag (t(u) %*%iter$sigma[,,i] %*% u)
//					sme_matmult_NC (mDu, m_amCurS [k], m_mTemp2) ;							//	2do: only calculate diagonal - this should save time!
//					sme_matmult_diag_NC (m_mTemp2, mDu, mDdd.GetColRef (k), FALSE, TRUE) ;
					sme_matmult_a_b_at_NC (mDu, m_amCurS [k], !m_mTemp2, TRUE, FALSE) ;
					CopyDiag_NC (mDdd.GetColRef (k), m_mTemp2) ;
				}

				//R	dd <- f.restr.eigen (dd, iter$csize, restr.fact, pa$zero.tol)
				if (m_nDeter)
					RestrictEigenValues_deter (!mDdd, m_vClustSize, m_dRestrFactor, m_dZeroTol, m_dUnRestrFact) ;
				else
					RestrictEigenValues (!mDdd, m_vClustSize, m_dRestrFactor, m_dZeroTol, m_dUnRestrFact) ;

				//R	u <- optvectors (iter, pa, u, dd, n.iter = n_iter2, bFast = bFast)
				OptVectors  (mDu, mDdd) ;
			}

			//R res <- 0
			double dRes = 0 ;

			m_mTemp2.Reshape_NC (m_p, m_p) ;
			m_vTemp1.Reshape_NC (m_p) ;

			//R for (i in 1:pa$K)
			for (k = m_K - 1; k != NAI; k--)
			{
				//R res <- res + iter$csize[i] * sum (diag (u %*% diag (1 / dd[,i]) %*% t (u) %*% iter$sigma[,,i]))
				m_vTemp1.Copy_NC (mDdd.GetColRef (k)) ;
				set_inv (*m_vTemp1) ;

				sme_matmult_a_diagb_at_NC (mDu, m_vTemp1, !m_mTemp2) ;
				dRes += m_vClustSize (k) * sme_sum_matmult_diag_NC (m_mTemp2, m_amCurS [k]) ;
			}

			//R	if (inic == 1 || res < resMin) { ... }
			if (!i || dRes < dMinRes)
			{
				mDBestU.Copy_NC (mDu) ;
				mDBestD.Copy_NC (mDdd) ;
				dMinRes = dRes ;
			}
		}

		//R for (i in 1:pa$K) 
		for (k = m_K - 1; k != NAI; k--)
		{
			//R iter$sigma[,,i] <- u_%*%diag(d_[,i])%*%t(u_)
			sme_matmult_a_diagb_at_NC (mDBestU, mDBestD.GetColRef (k), !m_amCurS [k]) ;
			if (k)
				m_amEVec [k].Copy_NC (mDBestU) ;
		}

		return max (mDBestD) > m_dZeroTol ;
	}

	void CTClust::OptVectors (SMatD &mU, const SCMatD &mIDiagD)
	{
		SMatD a_mj			(m_aTemp[0], m_p, m_p) ; 
		SMatD mUp1			(m_aTemp[0], m_p, m_p) ; 

		SMatD a_2mj			(m_aTemp[1], 2, 2) ;

		SMatD mUp2			(m_aTemp[2], m_p, m_p) ; 
		SMatD mUjm			(m_aTemp[2], m_p, 2) ; 

		SMatD mUU			(m_aTemp[3], m_p, 2) ; 

		SMatD mIDiagDInv	(m_aTemp[4], m_p, m_K) ;

		mIDiagDInv.Copy_NC	(mIDiagD) ;		//	2do	copy and invert can be done in one run!
		set_inv (*mIDiagDInv) ;

		SMatD m_mTemp1 (m_aTemp [11], m_p, m_p) ;
		SMatD m_mTemp2 (m_aTemp [12], m_p, 2) ;
		SVecD m_vTemp1 (m_aTemp [10], m_p) ;


		t_size n, j, m, k ;

		for (n = 0; n < m_dwIterTune3; n++)		//	the iteration-direction seems to influence the result - why is that??
			for (m = 0; m < m_p; m++)
				for (j = 0; j < m; j++)
				{
					a_mj.Reset (0) ;													//	V a_mj
					//R	for (i in 1:pa$K)
					for (k = m_K - 1; k != NAI; k--)
					{
						//R	a_mj <- a_mj + (d[j,j,i] - d[m,m,i]) / (d[j,j,i] * d[m,m,i]) * iter$csize[i] * iter$sigma[,,i]	##	weighted sum of covmats./
						const SCVecD &vD = mIDiagD.GetColRef (k) ;

						const double &dj = vD(j), &dm = vD(m) ;

						double dCurFact = (dj - dm) / (dj * dm) * m_vClustSize (k) ;	

						//EO<UOP::cPAbMa>::VSV_D (*m_amCurS[k], dCurFact, *a_mj) ;
						EO<SOP::ApaBmC>::VScVc_NC (*a_mj, dCurFact, m_amCurS[k]) ;
					}	//	m_mTemp_PxP_2 <==> a_mj

						//R a2_mj <- t (u[,c(j,m)]) %*% a_mj %*% u[,c(j,m)]	##	 ?

					mUjm.CopyCol_NC (0, j, mU) ;										//							V mUjm
					mUjm.CopyCol_NC (1, m, mU) ;

					sme_matmult_a_b_at_NC (mUjm, a_mj, !a_2mj, TRUE, FALSE) ;			//	A a_mj,			V a_2mj

						//	eigen_sqr (m_mTemp2x2_1, m_vTemp2, m_mTemp2x2_2, TRUE) ;
					m_mTemp2.Reshape_NC (2, 2) ;
					m_vTemp1.Reshape_NC (2) ;
														//	use REALLY fast method for 2x2 evec and eval
					sme_eigen_sym_2x2_norm_raw_NC (*m_vTemp1, *m_mTemp2, a_2mj, m_dZeroTolSqrt) ;//				A a_2mj	
//					sme_eigen_sqr_NC (a_2mj, m_vTemp1, m_mTemp2, TRUE) ;

						//R	uu <- u[,c(j,m)] %*% v$vectors
					sme_matmult_NC (mUjm, m_mTemp2, !mUU) ;								//							A mUjm					V mUU

					if (!m_dwOVV)
					{
						mU.CopyCol_NC (j, 0, mUU) ;	//	NEW
						mU.CopyCol_NC (m, 1, mUU) ;	//	NEW
					}
					else if (m_dwOVV == 1)
					{
							//R	up1 <- up2 <- u
						mUp1.Copy_NC (mU) ;													//			V mUp1
						mUp2.Copy_NC (mU) ;													//									V mUp2

							//R	up1[,c(j, m)] <- up2[,c(m, j)] <- uu
						mUp1.CopyCol_NC (j, 0, mUU) ;
						mUp1.CopyCol_NC (m, 1, mUU) ;
						mUp2.CopyCol_NC (m, 0, mUU) ;
						mUp2.CopyCol_NC (j, 1, mUU) ;

							//R resp1 <- resp2 <- 0
						double dResp1 = 0, dResp2 = 0 ;

						m_mTemp1.Reshape_NC (m_p, m_p) ;
						m_vTemp1.Reshape_NC (m_p) ;

	//					double dTemp1, dTemp2  ;

							//R	for (i in 1:pa$K)
						for (k = m_K - 1; k != NAI; k--)
						//for (k = 0; k < m_K; k++)
						{
							//R	resp1 <- resp1 + iter$csize[i] * sum (diag (up1 %*% diag (1 / diag (d[,,i])) %*% t (up1) %*% iter$sigma[,,i]))
							sme_matmult_a_diagb_at_NC (mUp1, mIDiagDInv.GetColRef (k), !m_mTemp1) ;
							//dResp1 += m_vClustSize(k) * sme_sum_matmult_diag_NC(m_mTemp1, m_amCurS[k]) ;
							dResp1 += m_vClustSize(k) * sme_sum_diag_At_matmult_B_NC(m_mTemp1, m_amCurS[k]) ;

	//						double dTemp1 = sme_sum_diag_At_matmult_B_NC(m_mTemp1, m_amCurS[k]) ;
	//						double dTemp2 = sme_sum_matmult_diag_NC(m_mTemp1, m_amCurS[k]) ;

							//R resp2 <- resp2 + iter$csize[i] * sum (diag (up2 %*% diag (1 / diag (d[,,i])) %*% t (up2) %*% iter$sigma[,,i]))
							sme_matmult_a_diagb_at_NC (mUp2, mIDiagDInv.GetColRef (k), !m_mTemp1) ;
							//dResp2 += m_vClustSize(k) * sme_sum_matmult_diag_NC(m_mTemp1, m_amCurS[k]) ;
							dResp2 += m_vClustSize(k) * sme_sum_diag_At_matmult_B_NC(m_mTemp1, m_amCurS[k]) ;
						}

	//if (m_dwTrace >= 2 && fabs (dResp1 - dResp2) < 1 && dResp1 != dResp2)			//	if dResp1 are almost similar dResp2, we might run into numerical instabilities
	//	Rprintf ("dResp1 - dResp2: %e\n", dResp1 - dResp2) ;

						if (fabs (dResp1 - dResp2) <= m_dZeroTolSqrt)				//	for numerical stability: if dResp1 and dResp2 are almost equal
							dResp1 = dResp2 ;										//	set them equal.

							//R idx.chk <- if(resp1 <= resp2) 1:2 else 2:1		## defines a "check index"
							//R neg1 <- t (uu[, idx.chk[1]]) %*% u[, j] < 0		## this "check index" is used here for checking which of the terms are negative
							//R neg2 <- t (uu[, idx.chk[2]]) %*% u[, m] < 0
							//R u <- if(resp1 <= resp2) up1 else up2

						BOOL bNeg1, bNeg2 ;
						if  (dResp1 <= dResp2) 
						{
							bNeg1 = sumprod (mUU.GetColRef (0), mU.GetColRef (j)) < 0 ;
							bNeg2 = sumprod (mUU.GetColRef (1), mU.GetColRef (m)) < 0 ;
							mU.Copy_NC (mUp1) ;												//			A mUp1
						}
						else
						{
							bNeg1 = sumprod (mUU.GetColRef (1), mU.GetColRef (j)) < 0 ;
							bNeg2 = sumprod (mUU.GetColRef (0), mU.GetColRef (m)) < 0 ;		//													A mUU
							mU.Copy_NC (mUp2) ;												//									A mUp2
						}

							//R if (neg1) u[, j] <- -u[, j]
						if (bNeg1)
							set_neg (*mU.GetColRef (j)) ;

							//R if (neg2) u[, m] <- -u[, m]
						if (bNeg2)
							set_neg (*mU.GetColRef (m)) ;
					}
//	Print_NC (mU) ;
//	eal_printf ("\n") ;
				}
	}
