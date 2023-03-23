/*
    SMat - Simple Matrix Classes v0.1beta
    Copyright (C) 2011 by Heinrich Fritz (heinrich_fritz@hotmail.com)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#define SMAT_FLAG_NO_INI
#include "smat.h"
#undef SMAT_FLAG_NO_INI

//#include <cstdlib>
#include "smat.meal.h"

	void *smat_malloc (t_size dwSize)	//	2do: move to meal (use e.g. R_malloc)
	{
		return new char [dwSize] ;
	}

	void smat_free (void *pPtr)
	{
		delete [] (char *) pPtr ;
	}

/*	void *smat_realloc (void *pOld, t_size dwOld, t_size dwNew)	//	2do: remove
	{
		void *pNew = smat_malloc (dwNew) ;

		if (dwOld && pOld)
			memcpy (pNew, pOld, sm_min (dwOld, dwNew)) ;
		delete [] (char *) pOld ;
		return pNew ;

	}
*/

////////////////
//	SMat_EXC  //
////////////////

	void SMat_EXC::OnException ()
	{
		meal_OnException (GetDate (), GetFile (), GetLine ()) ;
	}

	void SMat_EXC::OnUException ()
	{
		meal_OnUException () ;
	}

////////////////
//	SDataRef  //
////////////////

	SDataRef::SDataRef ()	//	2do: rename to SDataRef or CDataRef
		: m_pData (NULL), m_pDataEnd (NULL), m_dwRef (0), m_dwCount (0), m_bOwner  (TRUE), m_bStatic (FALSE)
	{

	}

	SDataRef::SDataRef (const t_size dwCount)
		: m_dwRef (0), m_bOwner (TRUE), m_bStatic (FALSE)
	{
		Alloc_NF (dwCount) ;
	}


	SDataRef::SDataRef (const t_size dwCount, void * const pData)
		: m_pData (pData), m_pDataEnd ((char *) pData + dwCount), m_dwRef (0), m_dwCount (dwCount), m_bOwner (FALSE), m_bStatic (FALSE)
	{
		ASSERT (pData) ;
	}

	SDataRef::~SDataRef ()
	{
		Free () ;
	}

	void SDataRef::Free ()
	{
		if (IsOwner ())
			//free (m_pData) ;
			smat_free (m_pData) ;
			//delete [] (char *) m_pData ;
		m_pData = NULL ;
		m_pDataEnd = NULL ;
		m_dwCount = 0 ;
	}

	void SDataRef::FreeIfIdle ()
	{
		if ((m_bStatic && m_dwRef <= 1) || 
			!m_dwRef)
			Free () ;
	}

	void Deref (SDataRef *&pRef)
	{
		SDataRef::sDeref (pRef) ;
	}

	void SDataRef::sDeref (SDataRef *&pRef)
	{
		if (pRef->Deref ())
			delete pRef ;
		pRef = NULL ;
	}

	SDataRef *SDataRef::Ref (SDataRef *&pRef)
	{
		if (pRef == this)
			return this ; 
		if (pRef)
			pRef->Deref () ;
		return Ref_NDR (pRef) ;
	}

	SDataRef *SDataRef::Ref_NDR (SDataRef *&pRef)	//	No DeRereference
	{
		IncRef () ;
		return pRef = this ;
	}

	BOOL SDataRef::Deref ()
	{
		--m_dwRef ;
		return !m_dwRef/* && !m_bStatic*/ ;
	}

	SDataRef *SDataRef::Recreate (t_size dwSize, SDataRef *&pRef)
	{
		THROW (IsOwner ()) ;

		if (GetRef () <= 1 || IsStatic ())
			Alloc (dwSize) ;
		else
			(new SDataRef (dwSize))->Ref (pRef) ;
		return pRef ;
	}

	void SDataRef::Alloc (t_size dwSize)
	{
		Free () ;
		Alloc_NF (dwSize) ;
	}

	void SDataRef::Alloc_NF (t_size dwSize)
	{
		if (dwSize)
		{
			GetDataRef () = smat_malloc (dwSize) ;
//			GetDataRef () = new char [dwSize] ;
			GetSizeRef () = dwSize ;
			GetDataEndRef () = (char *) GetData () + GetSize () ;
		}
		else
		{
			GetDataRef () = NULL ;
			GetSizeRef () = 0 ;
			GetDataEndRef () = NULL ;
		}
	}

	BOOL SDataRef::Require (t_size dwSize, SDataRef *&pRef)
	{
		if (dwSize <= GetSize ())
			return FALSE ;
		Recreate (dwSize, pRef) ;
		return TRUE ;
	}

	SDataRef &SDataRef::Empty ()
	{
		static SDataRef_Static emptyDR (0, FALSE) ;
		return emptyDR ;
	}

///////////////////////
//	SDataRef_Static  //
///////////////////////
	
	SDataRef_Static::SDataRef_Static (const t_size dwCount, const BOOL bStatic)
		: t_base (dwCount)
	{
		IncRef () ;
		if (bStatic)
			SetStatic () ;
	}

	SDataRef_Static::~SDataRef_Static ()
	{
		DecRef () ;
		ASSERT (!m_dwRef) ;
	}

	SDataRef_Static &SDataRef_Static::Require (t_size dwSize)
	{
		SDataRef *pFoo ;
		t_base::Require (dwSize, pFoo) ;
		return *this ;
	}

////////////////////
//	SDataRefCont  //
////////////////////

	SDataRefCont::SDataRefCont ()
		:m_ppData (NULL), m_dwSize (0)
	{
		
	}

	SDataRefCont::~SDataRefCont ()
	{
		Free () ;
	}

	void SDataRefCont::Require (t_size dwCount)
	{
		if (dwCount <= GetSize ())
			return ;

		t_pitem *pNew = new t_pitem [dwCount] ;			//	do a realloc for dataRef () with size dwCount
		if (GetSize ())
			memcpy (pNew, GetData (), GetMemSize ()) ;
		delete [] GetData () ;
		dataRef () = pNew ;

		//m_ppData = (t_pitem *) realloc (m_ppData, dwCount * sizeof (t_pitem)) ;
		//m_ppData = (t_pitem *) smat_realloc (m_ppData, m_dwSize * sizeof (t_pitem), dwCount * sizeof (t_pitem)) ;

		t_size i ;
		for (i = GetSize (); i < dwCount; ++i)
			GetData () [i] = new t_item ;
		sizeRef () = dwCount ;
	}

	SDataRef_Static &SDataRefCont::Item (t_size dwIdx)
	{
		Require (dwIdx + 1) ;
		return *GetData ()[dwIdx] ;
	}

	void SDataRefCont::Free ()
	{
		t_size i ;
		for (i = GetSize () - 1; i != NAI; i--)
			delete GetData () [i] ;
		//free (GetData ()) ;
		//smat_free (GetData ()) ;
		delete [] GetData () ;
		dataRef () = NULL ;
		sizeRef () = 0 ;
	}

	void SDataRefCont::FreeIfIdle ()
	{
		t_size i ;
		for (i = GetSize () - 1; i != NAI; i--)
			Item (i).FreeIfIdle () ;
	}

////////////////////
//	CDataCont_NT  //
////////////////////

	t_size &CDataCont_NT::GetInstanceCount ()
	{
		static t_size dwInstanceCount = 0 ;
		return dwInstanceCount ;
	}

////////////////////////
//	Temp - Functions  //
////////////////////////

	SDataRefCont &GetTempCont ()
	{
		static SDataRefCont cont;
		return cont ;
	}

	void RequireTemp (t_size dwCount)
	{
		GetTempCont().Require (dwCount) ;
	}

	SDataRef_Static &tempRef (t_size dwIdx)
	{
		return GetTempCont().Item (dwIdx) ;
	}


	SDataRefCont::CRefRange &GetPermTempRefRange ()
	{
		static SDataRefCont::CRefRange obj ;
		return obj ;
	}

	void FreeTempCont ()
	{
		GetTempCont ().FreeIfIdle () ;
	}

////////
////////
////////	//	smateasymath.h
////////		//	2do: sort this file!
////////

	void sme_matmult_R (const SCMatD &a, const SCMatD &b, SVMatD &c)
	{
		THROW (a.ncol () == b.nrow ()) ;
		c.Require (a.nrow (), b.ncol ()) ;

		sme_matmult_NC (a, b, c) ;
	}

	void sme_matmult (const SCMatD &a, const SCMatD &b, const SVMatD &c)
	{
		THROW (a.ncol () == b.nrow ()) ;
		THROW (a.nrow () == c.nrow () && c.ncol () == b.ncol ()) ;

		sme_matmult_NC (a, b, c) ;
	}

	void sme_matmult_NC (const SCMatD &a, const SCMatD &b, const SVMatD &c)
	{
		ASSERT (a.ncol () == b.nrow ()) ;
		ASSERT (a.nrow () == c.nrow () && c.ncol () == b.ncol ()) ;

		double one = 1.0, zero = 0.0;
		if (a.nrow () > 0 && a.ncol () > 0 && b.nrow () > 0 && b.ncol () > 0) {
		meal_gemm(
			"N", "N", 
			a.nrowPtrS (), b.ncolPtrS (), b.nrowPtrS (), 
			&one,
			a.GetData (), a.nrowPtrS (), 
			b.GetData (), b.nrowPtrS (), 
			&zero, c.GetData (), a.nrowPtrS ()
			);
		}
		else 
			c.Reset (0) ;
	}

	void sme_tmatmult_R (const SCMatD &a, const SCMatD &b, SVMatD &c, const BOOL bTransA, const BOOL bTransB)
	{
		c.Require (a.GetDim (bTransA), b.GetDim (!bTransB)) ;

		sme_tmatmult_NC (a, b, c, bTransA, bTransB) ;
	}

	void sme_tmatmult (const SCMatD &a, const SCMatD &b, const SVMatD &c, const BOOL bTransA, const BOOL bTransB)
	{
		THROW (a.GetDim (!bTransA) == b.GetDim (bTransB)) ;
		THROW (a.GetDim (bTransA) == c.nrow () && c.ncol () == b.GetDim (!bTransB)) ;

		sme_tmatmult_NC (a, b, c, bTransA, bTransB) ;
	}

	void sme_tmatmult_NC (const SCMatD &a, const SCMatD &b, const SVMatD &c, const BOOL bTransA, const BOOL bTransB)
	{
		//ASSERT (a.ncol () == b.nrow ()) ;
		ASSERT (a.GetDim (!bTransA) == b.GetDim (bTransB)) ;
		ASSERT (a.GetDim (bTransA) == c.nrow () && c.ncol () == b.GetDim (!bTransB)) ;
		//ASSERT (a.nrow () == c.nrow () && c.ncol () == b.ncol ()) ;

		double one = 1.0, zero = 0.0;
		if (a.nrow () > 0 && a.ncol () > 0 && b.nrow () > 0 && b.ncol () > 0) {
		meal_gemm(
			bTransA ? "T" : "N", bTransB ? "T" : "N", 
			a.GetDimPtrS_NC (bTransA), b.GetDimPtrS_NC (!bTransB), b.GetDimPtrS_NC (bTransB), 
			&one,
			a.GetData (), a.nrowPtrS (), 
			b.GetData (), b.nrowPtrS (), 
			&zero, c.GetData (), a.GetDimPtrS_NC (bTransA)
			);
		}
		else 
			c.Reset (0) ;
	}

	void sme_diag_R (const SVMatD &a, SVecD &c)
	{
		c.Require (sm_min<t_size> (a.nrow (), a.ncol ())) ;
		sme_diag_NC (a, c) ;
	}

	void sme_diag (const SVMatD &a, SVecD &c)
	{
		THROW (c.size () == ::sm_min (a.nrow (), a.ncol ())) ;
		sme_diag_NC (a, c) ;
	}

	void sme_diag_NC (const SVMatD &a, SVecD &c)
	{
		const t_size inc = a.GetColInc () + 1 ;

		double *pCur = a.GetData () ;
		t_size i ;
		for (i = 0; i < c.size (); i++)
		{
			c(i) = *pCur ;
			pCur += inc ;
		}
	}

/////////////////////
//	smat.random.h  //
/////////////////////

	double runif () { return meal_unif_rand () ; }
	double rnorm () { return meal_norm_rand () ; }
	double rexp () { return meal_exp_rand () ; }

	void runif (const SVData<double> &a)
	{
		runif_raw (a, a.GetDataEnd ()) ;
	}

	void runif_raw (double *p, double *pEnd)
	{
		for (; p < pEnd; ++p)
			*p = runif () ;
	}

	void runif_r (const SVData<double> &a)
	{
		runif_r_raw (a, a.GetDataEnd ()) ;
	}

	void runif_r_raw (double *p, double *pEnd)
	{
		for (--pEnd; p <= pEnd; --pEnd)
			*pEnd = runif () ;
	}

	void rnorm (const SVData<double> &a)
	{
		rnorm_raw (a, a.GetDataEnd ()) ;
	}

	void rnorm_raw (double *p, double *pEnd)
	{
		for (; p < pEnd; ++p)
			*p = rnorm () ;
	}

	void runif_raw (double *d, int l, double dL, double dU)
	{
		runif_raw (d, d + l, dL, dU) ;	
	}

	void runif_raw (double *d, double * const pEnd, double dL, double dU)
	{
		dU -= dL ;

		for ( ;d < pEnd; ++d)
			*d = runif () * dU + dL ;
	}

	void SampleNoReplace(int k, int n, int *y, int *x)
	{
		int i, j;
		for (i = n - 1; i != -1; i--)
			x[i] = i ;

		for (i = 0; i < k; i++)
		{
			j = int (n * runif()) ;
			y[i] = x[j] ;
			x[j] = x[--n];
		}
	}

///////////////////
//	Eigenvalues  //
///////////////////

	void sme_eigen_sqr_R	(const SCMatD &A,		SVecD &vVal,		SVMatD &mVec, BOOL bOrder)	{	ASSERT_TEMPRANGE (4, 4) ;	sme_eigen_sqr_RV	(!SMatD (tempRef (4), A), vVal, mVec, bOrder) ;	}
	void sme_eigen_sqr		(const SCMatD &A, const SVecD &vVal, const	SVMatD &mVec, BOOL bOrder)	{	ASSERT_TEMPRANGE (4, 4) ;	sme_eigen_sqr_V		(!SMatD (tempRef (4), A), vVal, mVec, bOrder) ;	}
	void sme_eigen_sqr_NC	(const SCMatD &A, const SVecD &vVal, const	SVMatD &mVec, BOOL bOrder)	{	ASSERT_TEMPRANGE (4, 4) ;	sme_eigen_sqr_NCV	(!SMatD (tempRef (4), A), vVal, mVec, bOrder) ;	}

	void sme_eigen_sqr_RV (const SVMatD &A, SVecD &vVal, SVMatD &mVec, BOOL bOrder)
	{
		t_size p = A.nrow () ;
		THROW (p == A.ncol ()) ;
		vVal.Redim (p) ;
		mVec.Redim (p, p) ;

		sme_eigen_sqr_NC (A, vVal, mVec, bOrder) ;
	}

	void sme_eigen_sqr_V (const SVMatD &A, const SVecD &vVal, const SVMatD &mVec, BOOL bOrder)	//	no redim
	{
		t_size  p = A.nrow () ;
		THROW (p == A.ncol ()) ;
		THROW (vVal.size () == p) ;
		THROW (mVec.nrow () == p && mVec.ncol () == p) ;

		sme_eigen_sqr_NC (A, vVal, mVec, bOrder) ;
	}

	
	void sme_eigen_sqr_NCV (const SVMatD &A, const SVecD &vVal, const SVMatD &mVec, BOOL bOrder)	//	no check
	{
		{
			ASSERT (A.nrow () == A.ncol ()) ;
			ASSERT (vVal.size () == A.ncol ()) ;
			ASSERT (mVec.nrow () == A.ncol () && mVec.ncol () == A.ncol ()) ;
		}
		int p = A.ncol () ;

//		SVecD &vWi = SVecD::TempFree_NC (1, p) ;
		SVecD vWi (tempRef (0), p) ;

		int nInfo ;
		int nWork = -1 ;

		double dWork ;

		meal_geev ("V", "N", &p, NULL, &p, NULL, NULL, NULL, &p, NULL, &p, &dWork, &nWork, &nInfo) ;

		nWork = int(dWork) ;
//		SVecD &vTemp = SVecD::TempFree_NC (2, nWork) ;
		SVecD vTemp (tempRef (1), nWork) ;

		if (bOrder)
		{
			ASSERT_TEMPRANGE (0, 3) ;
//			SMatD &mTempVec = SMatD::TempFree_NC (3, p, p) ;
			SMatD mTempVec (tempRef (2), p, p) ;

			meal_geev ("V", "N", &p, A, &p, vVal, vWi, mTempVec, &p, NULL, &p, vTemp, &nWork, &nInfo) ;

			// re-order
//			SVecN &vOrder = SVecN::TempFree_NC (0, p) ;
			SVecN vOrder (tempRef (3), p) ;
			meal_sort_order_rev (vVal, vOrder, p) ;

			mVec.CopyCol_Order_NC (mTempVec, *vOrder) ;
		}
		else
		{
			ASSERT_TEMPRANGE (0, 1) ;
			meal_geev ("V", "N", &p, A, &p, vVal, vWi, mVec, &p, NULL, &p, vTemp, &nWork, &nInfo) ;
		}

		THROW (!nInfo) ;
	}

	void sme_eigen_sym_2x2_norm_raw (double * const pdEval, double *const pdEVec, const double *const pd, const double &dZeroTol)
	{
		THROW (pd[2] == pd[1]) ;
		sme_eigen_sym_2x2_norm_raw_NC (pdEval, pdEVec, pd, dZeroTol) ;
	}

	void sme_eigen_sym_2x2_norm_raw_NC (double * const pdEval, double *const pdEVec, const double *const pd, const double &dZeroTol)
	{
		const double &a = pd[0], &b = pd[2], &d = pd[3] ;
//		ASSERT (b == pd[1]) ;

		const double	&dDet = pdEVec[0] = a * d - b * b,	//	"abusing" pdEVec as temporary memory
						&dTrace = pdEVec[1] = a + d ;
		double			&dTemp = pdEVec[2] ;

		double &dL1 = pdEval[0],
				&dL2 = pdEval[1] ;

		dTemp = sqrt (sm_sqr (dTrace) / 4.0 - dDet) ;		//	calculating eigenvalues
		dL2 = dTrace / 2 ;
		dL1 = dL2 + dTemp ;
		dL2 -= dTemp ;

		if (fabs (b) / (fabs (a) + fabs (d)) > dZeroTol)	//	calculate some kind of condition number, which is checked against a zero tolerance
		{
			pdEVec[0] = dL1 - d ;							//	calculating and norming the eigenvalues
			pdEVec[1] = sqrt (sm_sqr (pdEVec[0]) + sm_sqr (b)) ;
			pdEVec[0] /= pdEVec[1] ;
			pdEVec[1] = b / pdEVec[1] ;

			pdEVec[2] = dL2 - d ;
			pdEVec[3] = sqrt (sm_sqr (pdEVec[2]) + sm_sqr (b)) ;
			pdEVec[2] /= pdEVec[3] ;
			pdEVec[3] = b / pdEVec[3] ;
		}
		else
		{
			pdEVec[0] = pdEVec[3] = 1 ;
			pdEVec[1] = pdEVec[2] = 0 ;
		}
	}

///////////////
//	MatMult  //
///////////////

	void sme_matmult_a_diagb_at_R (const SCMatD &a, const SCVecD &b, SVMatD &c)
	{
		THROW (a.ncol () == b.size ()) ;
		c.Require (a.nrow (), a.nrow ()) ;
		sme_matmult_a_diagb_at_NC (a, b, c) ;
	}

	void sme_matmult_a_diagb_at	 (const SCMatD &a, const SCVecD &b, const SVMatD &c)
	{
		THROW (a.ncol () == b.size ()) ;
		THROW (c.nrow () == a.nrow () && c.ncol () == a.nrow ()) ;
		sme_matmult_a_diagb_at_NC (a, b, c) ;
	}

	void sme_matmult_a_diagb_at_NC (const SCMatD &a, const SCVecD &b, const SVMatD &c)
	{
		ASSERT (a.ncol () == b.size ()) ;
		ASSERT (c.nrow () == a.nrow () && c.ncol () == a.nrow ()) ;

//		SMatD &temp = SMatD::TempFree (0, a.nrow (), a.ncol ()) ;
		SMatD temp (tempRef (0), a.nrow (), a.ncol ()) ;
//		EO<OP::multiply>::McVcMd_byrow_NC (a, b, temp) ;
		EO<SOP::multiply>::MMcVct_NC (!temp, a, b) ;
		sme_tmatmult_NC (temp, a, c, FALSE, TRUE) ;
	}

	void sme_matmult_at_diagb_a_R (const SCMatD &a, const SCVecD &b, SVMatD &c)
	{
		THROW (a.nrow () == b.size ()) ;
		c.Require (a.ncol (), a.ncol ()) ;
		sme_matmult_at_diagb_a_NC (a, b, c) ;
	}

	void sme_matmult_at_diagb_a (const SCMatD &a, const SCVecD &b, const SVMatD &c)
	{
		THROW (a.nrow () == b.size ()) ;
		THROW (c.nrow () == a.ncol () && c.ncol () == a.ncol ()) ;
		sme_matmult_at_diagb_a_NC (a, b, c) ;
	}

	void sme_matmult_at_diagb_a_NC (const SCMatD &a, const SCVecD &b, const SVMatD &c)
	{
		ASSERT (a.nrow () == b.size ()) ;
		ASSERT (c.nrow () == a.ncol () && c.ncol () == a.ncol ()) ;

//		SMatD &temp = SMatD::TempFree (0, a.nrow (), a.ncol ()) ;
		SMatD temp (tempRef (0), a.nrow (), a.ncol ()) ;
		//EO<OP::multiply>::McVcMd_byrow_NC (a, b, temp) ;
		EO<SOP::multiply>::MMcVc (!temp, a, b) ;
		sme_tmatmult_NC (temp, a, c, TRUE, FALSE) ;
	}


////////////////////
//	matmult_a_at  //
////////////////////

	void sme_matmult_a_at_R		(const SCMatD &a, SVMatD &b, BOOL bTransA)
	{
		b.Require (a.GetDim (bTransA), a.GetDim (bTransA)) ;
		sme_matmult_a_at_NC (a, b, bTransA) ;
	}

	void sme_matmult_a_at		(const SCMatD &a, const SVMatD &b, BOOL bTransA)
	{
		THROW (b.nrow () == b.GetDim (bTransA) && b.ncol () == a.GetDim (bTransA)) ;
		sme_matmult_a_at_NC (a, b, bTransA) ;
	}

	void sme_matmult_a_at_NC	(const SCMatD &a, const SVMatD &b, BOOL bTransA)
	{
		ASSERT (b.nrow () == a.GetDim (bTransA) && b.ncol () == a.GetDim (bTransA)) ;

		sme_tmatmult_NC (a, a, !b, bTransA, !bTransA) ;
	}


//////////////////////
//	matmult_a_b_at  //
//////////////////////

	void sme_matmult_a_b_at_R		(const SCMatD &a, const SCMatD &b, SVMatD &c, BOOL bTransA, BOOL bTransB)
	{
		THROW (b.nrow () == b.ncol ()) ;
		THROW (a.GetDim (!bTransA) == b.nrow ()) ;
		c.Require (a.GetDim (bTransA), a.GetDim (bTransA)) ;
		sme_matmult_a_b_at_NC (a, b, c, bTransA, bTransB) ;
	}

	void sme_matmult_a_b_at		(const SCMatD &a, const SCMatD &b, SVMatD &c, BOOL bTransA, BOOL bTransB)
	{
		THROW (b.nrow () == b.ncol ()) ;
		THROW (a.GetDim (!bTransA) == b.nrow ()) ;
		THROW (c.nrow () == a.GetDim (bTransA) && c.ncol () == a.GetDim (bTransA)) ;
		sme_matmult_a_b_at_NC (a, b, c, bTransA, bTransB) ;
	}

	void sme_matmult_a_b_at_NC	(const SCMatD &a, const SCMatD &b, SVMatD &c, BOOL bTransA, BOOL bTransB)
	{
		ASSERT (b.nrow () == b.ncol ()) ;
		ASSERT (a.GetDim (!bTransA) == b.nrow ()) ;
		ASSERT (c.nrow () == a.GetDim (bTransA) && c.ncol () == a.GetDim (bTransA)) ;

//		SMatD &mTemp = SMatD::TempFree_NC (0, a.GetDim (bTransA), b.GetDim (!bTransB)) ;
		SMatD mTemp (tempRef (0), a.GetDim (bTransA), b.GetDim (!bTransB)) ;

		sme_tmatmult_NC (a, b, !mTemp, bTransA, bTransB) ;
		sme_tmatmult_NC (mTemp, a, c, FALSE, !bTransA) ;
	}

////////////////////
//	matmult_diag  //
////////////////////
	
	void sme_matmult_diag_R (const SCMatD &a, const SCMatD &b, SVecD &c)
	{
		THROW (a.ncol () == b.nrow ()) ;
		c.Require (sm_min (a.nrow (), b.ncol ())) ; 
		sme_matmult_diag_NC (a, b, c) ;
	}

	void sme_matmult_diag (const SCMatD &a, const SCMatD &b, const SVecD &c)
	{
		THROW (a.ncol () == b.nrow ()) ;
		THROW (sm_min (a.nrow (), b.ncol ()) == c.size ()) ;
		sme_matmult_diag_NC (a, b, c) ;
	}

	void sme_matmult_diag_NC (const SCMatD &a, const SCMatD &b, const SVecD &c)
	{
		ASSERT (a.ncol () == b.nrow ()) ;
		const t_size dwMin = sm_min (a.nrow (), b.ncol ()) ; 
		ASSERT (dwMin == c.size ()) ;

		c.Reset (0) ;

		const double *pdB = b.GetDataEnd () ;

		t_size i, j ;
		double *pdC = c.GetDataEnd () ;
		for (j = dwMin - 1; j != NAI; j--)
		{
			pdC-- ;
			for (i = a.ncol () - 1; i != NAI; i--)
				*pdC += a (j, i) * *--pdB ;
		}

/*
		const double *pdB = b ;

		t_size i, j ;

		double *pdC = c ;
		for (j = 0; j < dwMin; j++)
		{
			for (i = 0; i < a.ncol (); i++)
			{
				*pdC += a (j, i) * *pdB ;
				pdB++ ;
			}
			pdC++ ;
		}

*/
	}

	double sme_sum_matmult_diag (const SCMatD &a, const SCMatD &b)
	{
		double c = 0 ;
		sme_sum_matmult_diag (a, b, c) ;
		return c ;
	}

	double sme_sum_matmult_diag_NC (const SCMatD &a, const SCMatD &b)
	{
		double c = 0 ;
		sme_sum_matmult_diag_NC (a, b, c) ;
		return c ;
	}


	void sme_sum_matmult_diag (const SCMatD &a, const SCMatD &b, double &c)
	{
		THROW (a.ncol () == b.nrow ()) ;
		sme_sum_matmult_diag_NC (a, b, c) ;
	}

	void sme_sum_matmult_diag_NC (const SCMatD &a, const SCMatD &b, double &c)
	{
		ASSERT (a.ncol () == b.nrow ()) ;
		const t_size dwMin = sm_min (a.nrow (), b.ncol ()) ; 

//		if (c < 0)
		{
			c = 0 ;

			const double *pdB = b ;
			t_size i, j ;
			double dTemp ;
			for (j = 0; j < dwMin; j++)
			{
				dTemp = 0 ;
				for (i = 0; i < a.ncol (); i++)
					dTemp += a (j, i) * *pdB++ ;
				c += dTemp ;
			}
		}
/*		else
		{
			c = 0 ;
			const double *pdB = b.GetDataEnd () ;
			t_size i, j ;
			double dTemp ;
			for (j = dwMin - 1; j != NAI; j--)
			{
				dTemp = 0 ;
				for (i = a.ncol () - 1; i != NAI; i--)
					dTemp += a (j, i) * *--pdB ;
				c += dTemp ;
			}
		}
*/	}

	double sme_sum_diag_At_matmult_B (const SCMatD &a, const SCMatD &b)
	{
		double dSum ;
		sme_sum_diag_Bt_matmult_C (dSum, a, b) ;
		return dSum ;
	}

	double sme_sum_diag_At_matmult_B_NC (const SCMatD &a, const SCMatD &b)
	{
		double dSum ;
		sme_sum_diag_Bt_matmult_C_NC (dSum, a, b) ;
		return dSum ;
	}

	void sme_sum_diag_Bt_matmult_C (double &a, const SCMatD &b, const SCMatD &c)
	{	//	calculates  sum (diag (t (A) %*% B))
		THROW (b.nrow () == c.nrow ()) ;
		sme_sum_diag_Bt_matmult_C_NC (a, b, c) ;
	}

	void sme_sum_diag_Bt_matmult_C_NC (double &a, const SCMatD &b, const SCMatD &c)
	{	//	calculates  sum (diag (t (A) %*% B))

		ASSERT (b.nrow () == c.nrow ()) ;
		const t_size dwDim = sm_min (b.ncol (), b.ncol ()) ;

		double const * pdB = b, * const pdEndB = pdB + dwDim * b.nrow ();

		a = 0 ;
		EO<SOP::ApaBmC>::SVcVc_raw (a, pdB, pdEndB, c.GetData ()) ;
	}

////////////////////////
//	Calc Covariances  //
////////////////////////

	void cov_centered_R (SVMatD &a, const SCMatD &b, const double &dFact)
	{
		a.Require (b.ncol (), b.ncol ()) ;
		cov_centered_NC (a, b) ;
	}

	void cov_centered (const SVMatD &a, const SCMatD &b, const double &dFact)
	{
		THROW (a.nrow () == a.ncol ()) ;
		THROW (a.nrow () == b.ncol ()) ;
		cov_centered_NC (a, b) ;
	}

	void cov_centered_NC (const SVMatD &a, const SCMatD &b, const double &dFact)
	{		//	2do: think about only calculating half of it!
		ASSERT (a.nrow () == a.ncol ()) ;
		ASSERT (a.nrow () == b.ncol ()) ;

		sme_tmatmult_NC (b, b, a, TRUE, FALSE) ;
		EO<SOP::a_multiply>::VSc (*a, dFact / (b.nrow () - 1.0)) ;
	}

////////////
//	misc  //
////////////

	double median (const SCData<double> &a)
	{
		ASSERT_TEMPRANGE (0, 0) ;
		SVecD temp (tempRef (0), a.size ()) ;	//	2do: should be copied by constructor!
		temp.Copy_NC (a) ;
		return median_V (*temp) ;
	}

	double median_V (const SVData<double> &a)	//	2do: make this a template function!
	{
		int n = a.size () ;

		double *pA = a ;
		if (n < 3)
		{
			if (!n) 
				return meal_NaN () ;
			if (n == 1)
				return pA[0] ;
			return (pA[0] + pA[1]) / 2 ;
		}

		int nHalf = (n + 1) >> 1 ;

			//	odd length
		if (n & 1)
			return psort_V(a, nHalf-1);

			//	even length
		const double dTemp = psort_V(a, nHalf - 1) ;
		return (dTemp + min (pA + nHalf, pA + n)) / 2 ;
	}

	double mad0 (const SVData<double> &a) 
	{
		const double dCenter = median_V (a) ;
		EO<SOP::Aa_abs_AsB>::VSc (a, dCenter) ;

		return median_V (a) ;
	}

	double medianabs_V (const SVData <double> &a)
	{
		EO<SOP::a_abs>::V (a) ;
		return median_V (a) ;
	}

	double mad_V (const SVData<double> &a)
	{
		return mad0 (a) * 1.482602218505602 ;
	}

	double mad (const SCData<double> &a)
	{
		ASSERT_TEMPRANGE (0, 0) ;
		SVecD temp (tempRef (0), a.size ()) ;	//	2do: should be copied by the constructor!
		temp.Copy_NC (a) ;
		return mad_V (*temp) ;
	}

	void sort_R (const SCData<double> &a, SVecD &b)
	{
		b.Require (a.size ()) ;
		sort_NC (a, b) ;
	}

	void sort (const SCData<double> &a, const SVecD &b)
	{
		THROW (a.size () == b.size ()) ;
		sort_NC (a, b) ;
	}

	void sort_NC (const SCData<double> &a, const SVecD &b)
	{
		b.Copy_NC (a) ;
		meal_sort (b, b.size ()) ;
	}

	void sort (const SVData<double> &a)
	{
		meal_sort (a, a.size ()) ;
	}

	void sort_order (const SVData<double> &a, const SVData<int> &b)
	{
		THROW (a.size () == b.size ()) ;
		sort_order_NC (a, b) ;
	}

	void sort_order_NC (const SVData<double> &a, const SVData<int> &b)
	{
		ASSERT (a.size () == b.size ()) ;
		meal_sort_order (a, b, b.size ()) ;
	}

	void norm2 (double &dNorm, const SCData<double> &a)
	{
		dNorm = 0 ;
		EO<SOP::Apa_sqr_B>::SVc (dNorm, a) ;
		dNorm = sqrt (dNorm) ;
	}

	double norm2 (const SCData<double> &a)
	{
		double dNorm ;
		norm2 (dNorm, a) ;
		return dNorm ;
	}

//////////////////////////
//	Printing Functions  //
//////////////////////////

	void Print (const double &v)
	{
		meal_printf ("%f", v) ;
	}

	void Print (const float &v)
	{
		meal_printf ("%f", v) ;
	}

	void Print (const int &v)
	{
		meal_printf ("%d", v) ;
	}

	void Print (const t_size &v)
	{
		meal_printf ("%d", v) ;
	}


/*		//	2do: move to smat.test.cpp - file
	EXPORT void ex_mad (int *pnParIn, double *pdParOut, double *pdData)
	{
		*pdParOut = mad_V (*SVecD (pdData, *pnParIn)) ;
	}
*/
/*		//	2do: move to smat.test.cpp - file
	EXPORT void sme_tmatmult (int *pnParIn, double *pdA, double *pdB, double *pdC)	//	2do: put in different file!
	{		
		const t_size 
			na = pnParIn [0], pa = pnParIn [1],
			nb = pnParIn [2], pb = pnParIn [3],
			nc = pnParIn [4], pc = pnParIn [5] ;
		BOOL bTransA = pnParIn[6], bTransB = pnParIn[7] ;

		SCMatD a (pdA, na, pa) ;
		SCMatD b (pdB, nb, pb) ;
		SMatD c (pdC, nc, pc) ;

		sme_tmatmult_NC (a, b, !c, bTransA, bTransB) ;
	}
*/
/*		//	2do: move to smat.test.cpp - file
	EXPORT void ex_median (int *pnParIn, double *pdParOut, double *pdData)
	{
		*pdParOut = median_V (*SVecD (pdData, *pnParIn)) ;
	}
*/
/*		//	2do: move to smat.test.cpp - file
 	EXPORT void sme_matmult_a_diagb_at (int *pnParIn, double *pdA, double  *pdB, double *pdC)
	{
		int n = pnParIn [0], p = pnParIn [1] ;

		SCMatD a (pdA, n, p) ;
		SCVecD b (pdB, p) ;
		SMatD c (pdC, n, p) ;

		sme_matmult_a_diagb_at_NC (a, b, !c) ;
	}
*/
/*		//	2do: move to smat.test.cpp - file
	EXPORT void sme_eigen_sqr (int *pnParIn, double *pdA, double *pdVal, double *pdVec)	//	2do: put into different File
	{
		int p = pnParIn[0] ;
		BOOL bOrder = pnParIn[1] ;
		sme_eigen_sqr (!SMatD (pdA, p, p), SVecD (pdVal, p), !SMatD (pdVec, p, p), bOrder) ;
	}
*/
