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

#ifndef SMAT_STAT_H
#define SMAT_STAT_H
#include "smat.base.h"


	void cov_centered_R		(		SVMatD &a, const SCMatD &b, const double &dFact = 1) ;	//	2do: make this a template function
	void cov_centered		(const	SVMatD &a, const SCMatD &b, const double &dFact = 1) ;
	void cov_centered_NC	(const	SVMatD &a, const SCMatD &b, const double &dFact = 1) ;

	double median (const SCData<double> &a) ;
	double median_V (const SVData<double> &a) ;

	double mad0 (const SVData<double> &a) ;
	double mad (const SCData<double> &a) ;
	double mad_V (const SVData<double> &a) ;

	double medianabs_V (const SVData <double> &a) ;


	void norm2 (double &dNorm, const SCData<double> &a) ;
	double norm2 (const SCData<double> &a) ;

	
	
	template <class TA, class TB>
	inline void mean (TA & a, const SCData<TB> &b)
	{
		TA s = 0 ;
		EO<SOP::a_plus>::SVc (s, b) ;
		a = s / b.size () ;
	}

	template <class TA>
	inline TA mean (const SCData<TA> &a)
	{
		TA ret ;
		mean (ret, a) ;
		return ret ;
	}

	template <class TA, class TB>
	inline void sd (TA & a, const SCData<TB> &b)
	{
		var (a, b) ;
		a = (TA) sqrt ((double) a) ;
	}

	template <class TA, class TB>
	inline void sd_st (TA & a, const SCData<TB> &b)
	{
		var_st (a, b) ;
		a = (TA) sqrt ((double) a) ;
	}

	template <class TA, class TB>
	inline void var (TA & a, const SCData<TB> &b)
	{
		var_raw (a, b) ;
		a /= (b.size () - 1) ;
	}

	template <class TA, class TB>
	inline void var_raw (TA & a, const SCData<TB> &b)
	{
		TA m, d = 0;
		mean (m, b) ;

		EO<SOP::Apa_sq_BsC>::SScVc (d, m, b) ;
		a = d ;
	}

	template <class TA, class TB>
	inline void var_st (TA & a, const SCData<TB> &b)
	{
		var_st_raw (a, b) ;
		a /= (b.size () - 1) ;
	}

	template <class TA, class TB>
	inline void var_st_raw (TA & a, const SCData<TB> &b)
	{
		TA ss = 0, s = 0 ;
		EO<SOP::ApaC_Bpa_sq_C>::SSVc (s, ss, b) ;
		a = ss  - b.size () * sm_sqr (s / b.size ()) ;
	}


	template <class TA, class TB>
	TA sumprod (const SCData <TA> &a, const SCData <TB> &b)
	{
		TA ret = 0 ;
//		EO<OP::CpaAmB>::vcdvcds (a, b, ret) ;
		EO<SOP::ApaBmC>::SVcVc (ret, a, b) ;
		return ret ;
	}

	template <class TA, class TB, class TC>
	void sumprod (const SCData <TA> &a, const SCData <TB> &b, TC &c)
	{ //EO<OP::CpaAmB>::vcdvcds (a, b, c) ;
		EO<SOP::ApaBmC>::SVcVc (c, a, b) ;	}

	template <class TA>
	TA sum (const SCData<TA> &a)
	{
		TA ret = 0 ;
		sum (a, ret) ;
		return ret ;
	}

	template <class TA, class TB>
	void sum (const SCData<TA> &a, TB &sum)
	{
		sum = 0 ;
		EO<SOP::a_plus>::SVc (sum, a) ;
	}

	template <class TA, class TB>
	void sum (TA *pA, t_size n, TB &sum)
	{
		sum = 0 ;
		EO<SOP::a_plus>::SVc_raw (sum, pA, pA + n) ;
	}

	template <class TA>
	TA sum (TA *pA, t_size n)
	{
		TA ret = 0 ;
		sum (pA, n, ret) ;
		return ret ;
	}

	template <class TA>
	void cumsum (const SVData <TA> &a)
	{
		EO<SOP::a_plus>::V_pairs (a) ;
	}

	template <class TA>
	void cumsum_r (const SVData <TA> &a)
	{
		EO<SOP::a_plus>::V_pairs_r (a) ;
	}

	template <class TA, class TB>
	void colSums_R (SVec<TA> &a, const SCMat<TB> &b)
	{
		a.Requires (b.ncol ()) ;
		colSums_NC (a, b) ;
	}

	template <class TA, class TB>
	void colSums (const SVData<TA> &a, const SCMat<TB> &b)
	{
		THROW (a.size () == b.ncol ()) ;
		colSums_NC (a, b) ;
	}

	template <class TA, class TB>
	void colSums_NC (const SVData<TA> &a, const SCMat<TB> &b)
	{
		a.Reset (0) ;
		EO<SOP::a_plus>::VetMcd_NC (a, b) ;
	}

////////////////
//	Products  //
////////////////

	template <class TA>
	TA prod (const SCData<TA> &a)
	{
		TA ret = 1 ;
		prod (a, ret) ;
		return ret ;
	}

	template <class TA>
	void prod (const SCData<TA> &a, TA &prod)
	{
		prod = 1 ;
		EO<SOP::a_multiply>::SVc (prod, a) ;
	}

	template <class TA, class TB>
	void colProds_R (SVec<TA> &a, const SCMat<TB> &b)
	{
		a.Requires (b.ncol ()) ;
		colProds_NC (a, b) ;
	}

	template <class TA, class TB>
	void colProds (const SVData<TA> &a, const SCMat<TB> &b)
	{
		THROW (a.size () == b.ncol ()) ;
		colProds_NC (a, b) ;
	}

	template <class TA, class TB>
	void colProds_NC (const SVData<TA> &a, const SCMat<TB> &b)
	{
		a.Reset (1) ;
		EO<SOP::a_multiply>::VetMcd_NC (a, b) ;
	}




#endif	//	#ifndef SMAT_STAT_H
