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

#ifndef SMAT_MISC_H
#define SMAT_MISC_H

#include "smat.base.h"
#include "smat.elop.h"

/////////////////////
//	Min Max Limit  //
/////////////////////

	template <class TA>
	void limit (const SVData<TA> &a, const TA &l, const TA &u)
	{
		THROW (l <= u) ;
		limit_nc (a, l, u) ;
	}

	template <class TA, class TB, class TC>
	void limit_NC (const SVData<TA> &a, const TB &l, const TC &u)
	{
		ASSERT (l <= u) ;
		EO<SOP::a_limit>::VScSc (a, l, u) ;
	}

	template <class TA, class TB>
	void limit_l (const SVData<TA> &a, const TB &l)
	{
		EO<SOP::a_limit_l>::VSc (a, l) ;
	}

	template <class TA, class TB>
	void limit_u (const SVData<TA> &a, const TB &u)
	{
		EO<SOP::a_limit_u>::VSc (a, u) ;
	}

	template <class TA>
	const TA &max (const SCData<TA> &a)
	{
		return EO<SOP::is_greater>::Vc_transitive (a) ;
	}

	template <class TA>
	const TA &min (const SCData<TA> &a)
	{
		return EO<SOP::is_less>::Vc_transitive (a) ;
	}

	template <class TA>
	void minmax (const SCData<TA> &a, TA &min, TA &max, BOOL bInit= TRUE)
	{
		if (!a.size ())
			return ;
		if (bInit)
		{
			const double *pA = a ;
			min = max = *pA ;
			EO<SOP::a_minmax>::SSVc_raw (min, max, pA + 1, a.GetDataEnd ()) ;
		}
		else
			EO<SOP::a_minmax>::SSVc (min, max, a) ;
	}

	template <class TA>
	t_size getMaxIdx (const SCData <TA> &a)
	{
		t_size idx = 0 ;
		EO<SOP::a_max_idx>::SVc (idx, a) ;
		return idx ;
	}

	template <class TA>
	TA min (const TA *pA, t_size n)
	{
		return min (pA, pA + n) ;
	}

	template <class TA>
	TA min (const TA *pA, TA const * const pEndA)
	{
		TA min = *pA ;

		for (++pA; pA < pEndA; ++pA)
			if (*pA < min)
				min = * pA ;
		return min ;
	}

	template <class TA>
	const TA *minP (const TA *pA, t_size n)
	{
		return minP (pA, pA + n) ;
	}

	template <class TA>
	const TA *minP (const TA *pA, TA const * const pEndA)
	{

		const TA *pMin = pA ;

		for (++pA; pA < pEndA; ++pA)
			if (*pA < *pMin)
				pMin = pA ;
		return pMin ;
	}


	template <class TA>
	const TA max (const TA *pA, t_size n)
	{
		return max (pA, pA + n) ;
	}

	template <class TA>
	const TA max (const TA *pA, TA const * const pEndA)
	{
		TA max = *pA ;

		for (++pA; pA < pEndA; ++pA)
			if (*pA > max)
				max = *pA ;
		return max ;
	}

	template <class TA>
	const TA *maxP (const TA *pA, t_size n)
	{
		return maxP (pA, pA + n) ;
	}

	template <class TA>
	const TA *maxP (const TA *pA, TA const * const pEndA)
	{

		const TA *pMax = pA ;

		for (++pA; pA < pEndA; ++pA)
			if (*pA > *pMax)
				pMax = pA ;
		return pMax ;
	}

////////////
//	Misc  //
////////////

	template <class TA>
	void set_neg (const SVData<TA> &a)
	{
//		EO<OPA::neg>::vd (a) ;
		EO<SOP::a_neg>::V (a) ;
	}

	template <class TA>
	void set_inv (const SVData<TA> &a)
	{
//		EO<OPA::inv>::vd (a) ;
		EO<SOP::a_inv>::V (a) ;
	}

	template <class TA, class TB>
	BOOL equal (const SCData <TA> &a, const SCData <TB> &b)
	{
		if (a.size () != b.size ())
			return FALSE ;

		TA const *pA = a ;
		TA const * const pEndA = a.GetDataEnd () ;
		TB const *pB = b ;

		while (pA < pEndA)
		{
			if (*pA != *pB)
				return FALSE ;
			++pA ;
			++pB ;
		}

		return TRUE ;
	}

	template <class TA, class TB>
	t_size CountMatches (const SCData <TA> &a, const TB &b)
	{
		t_size n = 0 ;
		EO<SOP::inc_a_if_b_equals_c>::SScVc (n, b, a) ;
		return n ;
	}

	template <class TA>
	t_size CountTrue (const SCData <TA> &a)
	{
		t_size n = 0 ;
		EO<SOP::inc_a_if_b>::SVc(n, a) ;
		return n ;
	}

//////////////////////////
//	Printing Functions  //
//////////////////////////

	void Print (const double &v) ;
	void Print (const float &v) ;
	void Print (const int &v) ;
	void Print (const t_size &v) ;

	template <class TA>
	void Print_NC (const t_size dwRe, const SCMat <TA> &a, const t_size dwCe)
	{
		Print_NC (0, dwRe, a, 0, dwCe) ;
	}

	template <class TA>
	void Print_NC (const SCMat <TA> &a, const t_size dwCs, const t_size dwCe)
	{
		Print_NC (0, a.nrow (), a, dwCs, dwCe) ;
	}

	template <class TA>
	void Print_NC (const t_size dwRs, const t_size dwRe, const SCMat <TA> &a)
	{
		Print_NC (dwRs, dwRe, a) ;
	}
	
	template <class TA>
	void Print_NC (const SCMat <TA> &a, const t_size dwCe)
	{
		Print_NC (0, a.nrow (), a, 0, dwCe) ;
	}

	template <class TA>
	void Print_NC (const t_size dwRe, const SCMat <TA> &a)
	{
		Print_NC (0, dwRe, a, 0, a.ncol ()) ;
	}

	template <class TA>
	void Print_NC (const SCMat <TA> &a)
	{
		Print_NC (0, a.nrow (), a, 0, a.ncol ()) ;
	}

	template <class TA>
	void Print_NC (const t_size dwRs, const t_size dwRe, const SCMat <TA> &a, const t_size dwCs, const t_size dwCe)
	{
		ASSERT (dwRs <= dwRe) ;
		ASSERT (dwCs <= dwCe) ;
		ASSERT (dwRe <= a.nrow ()) ;
		ASSERT (dwCe <= a.ncol ()) ;

		t_size r, c ;
		const TA *pA ;
		for (r = dwRs; r < dwRe; ++r)
		{
			pA = a.GetData (r, dwCs) ;
			for (c = dwCs; c < dwCe; ++c)
			{
				Print (*pA) ;
				meal_printf ("\t") ; 
				pA += a.GetColInc () ;
			}
			meal_printf ("\n") ; 
		}
		meal_printf ("\n") ; 
	}

	template <class TA>
	void Print_NC (const SCData <TA> &a)
	{
		const TA *pA = a, * const pEndA = a.GetDataEnd () ;

		while (pA < pEndA)
		{
				Print (*pA) ;
				meal_printf ("\t") ; 
				++pA ;
		}
		meal_printf ("\n") ;
	}
#endif	//	#ifndef SMAT_MISC_H
