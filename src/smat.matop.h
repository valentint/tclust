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

#ifndef SMAT_MATOP_H
#define SMAT_MATOP_H



	void sme_matmult_R				(const SCMatD &a, const SCMatD &b, SVMatD &c) ;
	void sme_matmult				(const SCMatD &a, const SCMatD &b, const SVMatD &c) ;
	void sme_matmult_NC				(const SCMatD &a, const SCMatD &b, const SVMatD &c) ;

	void sme_tmatmult_R				(const SCMatD &a, const SCMatD &b, SVMatD &c, const BOOL bTransA, const BOOL bTransB) ;
	void sme_tmatmult				(const SCMatD &a, const SCMatD &b, const SVMatD &c, const BOOL bTransA, const BOOL bTransB) ;
	void sme_tmatmult_NC			(const SCMatD &a, const SCMatD &b, const SVMatD &c, const BOOL bTransA, const BOOL bTransB) ;

	void sme_matmult_a_at_R			(const SCMatD &a, SVMatD &b, BOOL bTransA) ;
	void sme_matmult_a_at			(const SCMatD &a, const SVMatD &b, BOOL bTransA) ;
	void sme_matmult_a_at_NC		(const SCMatD &a, const SVMatD &b, BOOL bTransA) ;

	void sme_matmult_a_b_at_R		(const SCMatD &a, const SCMatD &b, SVMatD &c, BOOL bTransA = FALSE, BOOL bTransB = FALSE) ;
	void sme_matmult_a_b_at		(const SCMatD &a, const SCMatD &b, SVMatD &c, BOOL bTransA = FALSE, BOOL bTransB = FALSE) ;
	void sme_matmult_a_b_at_NC	(const SCMatD &a, const SCMatD &b, SVMatD &c, BOOL bTransA = FALSE, BOOL bTransB = FALSE) ;


	void sme_matmult_a_diagb_at_R	(const SCMatD &a, const SCVecD &b, SVMatD &c) ;
	void sme_matmult_a_diagb_at		(const SCMatD &a, const SCVecD &b, const SVMatD &c) ;
	void sme_matmult_a_diagb_at_NC	(const SCMatD &a, const SCVecD &b, const SVMatD &c) ;

	void sme_matmult_at_diagb_a_R	(const SCMatD &a, const SCVecD &b, SVMatD &c) ;
	void sme_matmult_at_diagb_a		(const SCMatD &a, const SCVecD &b, const SVMatD &c) ;
	void sme_matmult_at_diagb_a_NC	(const SCMatD &a, const SCVecD &b, const SVMatD &c) ;

	void sme_matmult_diag_R		(const SCMatD &a, const SCMatD &b, SVecD &c) ;			//	C <- diag (A %*% B)
	void sme_matmult_diag		(const SCMatD &a, const SCMatD &b, const SVecD &c) ;	//	C <- diag (A %*% B)
	void sme_matmult_diag_NC	(const SCMatD &a, const SCMatD &b, const SVecD &c) ;	//	C <- diag (A %*% B) # NC

	void sme_sum_matmult_diag (const SCMatD &a, const SCMatD &b, double &c) ;		//	C <- sum (C <- diag (A %*% B))
	void sme_sum_matmult_diag_NC (const SCMatD &a, const SCMatD &b, double &c) ;	//	C <- sum (C <- diag (A %*% B)) # NC
	double sme_sum_matmult_diag (const SCMatD &a, const SCMatD &b) ;
	double sme_sum_matmult_diag_NC (const SCMatD &a, const SCMatD &b) ;

	double	sme_sum_diag_At_matmult_B		(			const SCMatD &a, const SCMatD &b) ;
	double	sme_sum_diag_At_matmult_B_NC	(			const SCMatD &a, const SCMatD &b) ;

	void	sme_sum_diag_Bt_matmult_C		(double &a, const SCMatD &b, const SCMatD &c) ;
	void	sme_sum_diag_Bt_matmult_C_NC	(double &a, const SCMatD &b, const SCMatD &c) ;

	void sme_diag_R		(const SVMatD &a, SVecD &c) ;										//	2do: make this a template function
	void sme_diag		(const SVMatD &a, SVecD &c) ;
	void sme_diag_NC	(const SVMatD &a, SVecD &c) ;

//////////////////////////
//	Matrix Diagonals	//
//////////////////////////

	template <class TA, class TB>
	void SetDiag_R (SVMat <TA> &a, const SCData<TB> &b)
	{
		a.Require (b.size (), b.size ()) ;
		SetDiag_NC (a, b) ;
	}

	template <class TA, class TB>
	void SetDiag (const SVMat <TA> &a, const SCData<TB> &b)
	{
		THROW (a.ncol() == b.size ()) ;
		THROW (a.nrow() == b.size ()) ;
		SetDiag_NC (a, b) ;
	}

	template <class TA, class TB>
	void SetDiag_NC (const SVMat <TA> &a, const SCData<TB> &b)
	{										//2do: implement for non square matrices!
		ASSERT (a.ncol() == b.size ()) ;
		ASSERT (a.nrow() == b.size ()) ;
		
		t_size i, j ;

		double *pA = a.GetDataEnd () ;
		const double *pB = b.GetDataEnd () ;

		for (i = a.ncol () - 2; i != NAI; i--)
		{
			*--pA =  *--pB ;
			for (j = a.ncol () - 1; j != NAI; j--)
				*--pA = 0 ;
		}

		*--pA = *--pB ;
	}

	template <class TA>
	void SetDiag (const SVMat <TA> &a)
	{
		t_size dwR, dwC ;
		t_size const dwREnd = a.nrow () ;
		TA * pA = a ;
		TA const * const pAEnd = a.GetDataEnd () ;

		for (dwC = 0; pA < pAEnd; ++dwC)
			for (dwR = 0; dwR < dwREnd; ++dwR)
			{
				*pA = (dwC == dwR) ? 1.0 : 0.0 ;
				++pA ;
			}
	}


	template <class TA>
	void SetDiag_sq (const SVMat <TA> &a)
	{
		THROW (a.ncol () == a.nrow ()) ;
		SetDiag_sq_NC (a) ;
	}

	template <class TA>
	void SetDiag_sq_NC (const SVMat <TA> &a)
	{
		ASSERT (a.ncol () == a.nrow ()) ;
		const t_size inc = a.GetColInc () ;
		TA *pA = a ;
		TA *pEndA = a.GetDataEnd () ; //pA + a.nrow () * a.ncol () ;

		*pA = 1 ;
		++pA ;

		for (; pA < pEndA;)
		{
			pA = Reset (pA, pA + inc) ;
			*pA = 1 ;
			++pA ;
		}

//		Reset (pA, a.GetDataEnd ()) ;
	}

	template <class TA>
	void SetAntiDiag_sq (const SVMat <TA> &a)
	{
		THROW (a.ncol () == a.nrow ()) ;
		SetAntiDiag_sq_NC (a) ;
	}

		template <class TA>
	void SetAntiDiag_sq_NC (const SVMat <TA> &a)
	{
		ASSERT (a.ncol () == a.nrow ()) ;
		const t_size inc = a.GetColInc () - 2 ;
		TA *pA = a ;
		TA *pEndA = pA + a.nrow () * (a.ncol () - 1) + 1 ;

		*pA = 0 ;
		++pA ;

		for (; pA < pEndA; ++pA)
		{
			pA = Reset (pA, pA + inc) ;
			*pA = 1 ;
		}

		Reset (pA, a.GetDataEnd ()) ;
	}



#endif	//	#ifndef SMAT_MATOP_H
