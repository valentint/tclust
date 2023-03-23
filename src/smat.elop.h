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

#ifndef SMAT_ELOP_H
#define SMAT_ELOP_H


#include "smat.base.h"
#include <math.h>		//	move to smat.def.h
#include "smat.meal.h"	//	?? move to smat.def.h ??


#define CALC_1_0(TYPE)	public: template <class TA>						 inline static TYPE Calc (const TA &a)
#define CALC_1_1(TYPE)	public: template <class TA>						 inline static TYPE Calc (		TA &a)

#define CALC_2_0(TYPE)	public: template <class TA, class TB>			 inline static TYPE Calc (const TA &a, const	TB &b)
#define CALC_2_1(TYPE)	public: template <class TA, class TB>			 inline static TYPE Calc (		TA &a, const	TB &b)
#define CALC_2_2(TYPE)	public: template <class TA, class TB>			 inline static TYPE Calc (		TA &a,			TB &b)

#define CALC_3_0(TYPE)	public: template <class TA, class TB, class TC>	 inline static TYPE Calc (const TA &a, const	TB &b, const	TC &c)
#define CALC_3_1(TYPE)	public: template <class TA, class TB, class TC>	 inline static TYPE Calc (		TA &a, const	TB &b, const	TC &c)
#define CALC_3_2(TYPE)	public: template <class TA, class TB, class TC>	 inline static TYPE Calc (		TA &a,			TB &b, const	TC &c)
#define CALC_3_3(TYPE)	public: template <class TA, class TB, class TC>	 inline static TYPE Calc (		TA &a,			TB &b,			TC &c)

#define CALC_4_0(TYPE)	public: template <class TA, class TB, class TC, class TD>	 inline static TYPE Calc (const TA &a, const	TB &b, const	TC &c, const	TD &d)
#define CALC_4_1(TYPE)	public: template <class TA, class TB, class TC, class TD>	 inline static TYPE Calc (		TA &a, const	TB &b, const	TC &c, const	TD &d)
#define CALC_4_2(TYPE)	public: template <class TA, class TB, class TC, class TD>	 inline static TYPE Calc (		TA &a,			TB &b, const	TC &c, const	TD &d)
#define CALC_4_3(TYPE)	public: template <class TA, class TB, class TC, class TD>	 inline static TYPE Calc (		TA &a,			TB &b,			TC &c, const	TD &d)
#define CALC_4_4(TYPE)	public: template <class TA, class TB, class TC, class TD>	 inline static TYPE Calc (		TA &a,			TB &b,			TC &c,			TD &d)

#define CALC_5_1(TYPE)	public: template <class TA, class TB, class TC, class TD, class TE>	 inline static TYPE Calc (		TA &a, const	TB &b, const	TC &c, const	TD &d, const	TE &e)

	//	some functions: Renamed 4 Clang (clang complains when calling ::exp from a method of a class named "exp" )
	static	inline double r4c_exp (double a0) { return ::exp (a0) ; }
	static	inline double r4c_log (double a0) { return ::log (a0) ; }
	static	inline double r4c_pow (double a0, double a1) { return ::pow (a0, a1) ; }

	class SOP		//	standard operators
	{				//	2do: renamo to OP


	public:
		class assign		{ CALC_2_1(void)	{ a = (TA) b ; }	} ;
		class add			{ CALC_3_1(void)	{ a = (TA)		(b + c) ; }	} ;
		class subtract		{ CALC_3_1(void)	{ a = (TA)		(b - c) ; }	} ;
		typedef subtract sub ;

		class divide		{ CALC_3_1(void)	{ a = (TA)		(b / c) ; }	} ;
		class divide_r		{ CALC_3_1(void)	{ a = (TA)		(c / b) ; }	} ;
		class multiply		{ CALC_3_1(void)	{ a = (TA)		(b * c) ; }	} ;
		class pow			{ CALC_3_1(void)	{ a = (TA) ::r4c_pow(b, c)	; } } ;
		class pow_r			{ CALC_3_1(void)	{ a = (TA) ::r4c_pow(c, b)	; } } ;

		class mod			{ CALC_3_1(void)	{ a = (TA)		(b % c)	; }	} ;
		class neg			{ CALC_2_1(void)	{ a = -b ; }	} ;
		class inv			{ CALC_2_1(void)	{ a = (TA) (1 / b) ; }	} ;

		class exp			{ CALC_2_1(void)	{ a = (TA) ::r4c_exp ((double) b) ; }	} ;
		class log			{ CALC_2_1(void)	{ a = (TA) ::r4c_log ((double) b) ; }	} ;
		class pow2			{ CALC_2_1(void)	{ a = (TA) (b * b) ; }	} ;
		class sign			{ CALC_2_1(void)	{ a = (b > 0) ? (TA) 1 : (b < 0) ? (TA) -1: (TA) 0 ; }	} ;
		class abs			{ CALC_2_1(void)	{ a = (b < 0) ? (TA) -b : (TA) b  ; }	} ;

/*		class or			{ CALC_3_1(void)	{ c = (TC) (a | b) ; }	} ;
		class and			{ CALC_3_1(void)	{ c = (TC) (a & b) ; }	} ;
		class OR			{ CALC_3_1(void)	{ c = (TC) (a || b) ; }	} ;
		class AND			{ CALC_3_1(void)	{ c = (TC) (a && b) ; }	} ;
		class xor			{ CALC_3_1(void)	{ c = (TC) (a ^ b) ; }	} ;
		class NOT			{ CALC_2_1(void)	{ b = !a ; }	} ;
		class not			{ CALC_2_1(void)	{ b = ~a ; }	} ;
*/
		class gr			{ CALC_3_1(void)	{ a = (TA) (b > c) ; }	} ;
		class greq			{ CALC_3_1(void)	{ a = (TA) (b >= c) ; }	} ;
		class eq			{ CALC_3_1(void)	{ a = (TA) (b == c) ; }	} ;
		class le			{ CALC_3_1(void)	{ a = (TA) (b < c) ; }	} ;
		class leeq			{ CALC_3_1(void)	{ a = (TA) (b <= c) ; }	} ;

		//	2 arguments
		class a_plus		{ CALC_2_1(void)	{ a += (TA) b ; } } ;
		class a_minus		{ CALC_2_1(void)	{ a -= (TA) b ; } } ;
		class a_multiply	{ CALC_2_1(void)	{ a *= (TA) b ; } } ;
		class a_divide		{ CALC_2_1(void)	{ a /= (TA) b ; } } ;
		class a_divide_r	{ CALC_2_1(void)	{ a =  ((TA) b / a) ; } } ;
		class a_mod			{ CALC_2_1(void)	{ a %= (TA) b ; } } ;

		class a_or			{ CALC_2_1(void)	{ a |= (TA) b ; } } ;
		class a_and			{ CALC_2_1(void)	{ a &= (TA) b ; } } ;
		class a_xor			{ CALC_2_1(void)	{ a ^= (TA) b ; } } ;
		class a_OR			{ CALC_2_1(void)	{ a = a || (TA) b ; } } ;
		class a_AND			{ CALC_2_1(void)	{ a = a &&(TA) b ; } } ;

		//	1 argument
		class a_neg			{ CALC_1_1(void)	{ a = -a ; } } ;
		class a_exp			{ CALC_1_1(void)	{ a = (TA) ::r4c_exp ((double) a) ; }	} ;
		class a_log			{ CALC_1_1(void)	{ a = (TA) ::r4c_log ((double) a) ; }	} ;
		class a_sqrt		{ CALC_1_1(void)	{ a = (TA) ::sqrt ((double) a) ; }	} ;
		class a_pow2		{ CALC_1_1(void)	{ a = sm_sqr(a) ; }	} ;
		class a_pow3		{ CALC_1_1(void)	{ a = sm_sqr(a) * a ; }	} ;
		class a_pow4		{ CALC_1_1(void)	{ a = sm_sqr (sm_sqr (a)) ; }	} ;
		class a_abs			{ CALC_1_1(void)	{ if (a < 0) a = -a ; }	} ;
		class a_pow			{ CALC_2_1(void)	{ a = ::r4c_pow ((double) a, (double)b) ; } } ;


		class a_not			{ CALC_1_1(void)	{ a = ~a ; }	} ;
		class a_NOT			{ CALC_1_1(void)	{ a = !a ; }	} ;
		class a_inv			{ CALC_1_1(void)	{ a = 1 / a ; }	} ;

		class a_limit		{ CALC_3_1(void)	{ ASSERT (b <= c) ; if (a < b) a = b ; else if (a > c) a = c ; }	} ;
		class a_limit_l		{ CALC_2_1(void)	{ if (a < b) a = b ; }	} ;
		class a_limit_u		{ CALC_2_1(void)	{ if (a < b) a = b ; }	} ;
		class a_minmax		{ CALC_3_2(void)	{ ASSERT (a <= b) ;	if (a > c) a = c ; else if (b < c) b = c ; }	} ;

		class is_greater	{ CALC_2_0(BOOL)	{ return a >  b ; } } ;
		class is_greatereq	{ CALC_2_0(BOOL)	{ return a >= b ; } } ;
		class is_less		{ CALC_2_0(BOOL)	{ return a <  b ; } } ;
		class is_lesseq		{ CALC_2_0(BOOL)	{ return a <= b ; } } ;

// class 2 operators
		class ApaBmC		{ CALC_3_1(void)	{ a += (TA) (b * c) ; }	} ;
		class ApaBdC		{ CALC_3_1(void)	{ a += (TA) (b / c) ; }	} ;
		class ApaCdB		{ CALC_3_1(void)	{ a += (TA) (c / b) ; }	} ;
		class Apa1dB		{ CALC_2_1(void)	{ a += (TA) (1 / b) ; }	} ;
		class AsaBmC		{ CALC_3_1(void)	{ a -= (TA) (b * c) ; }	} ;

		class Aa_abs_AsB	{ CALC_2_1(void)	{ a = fabs (a - b) ; } } ;

		class Apa_sqr_B		{ CALC_2_1(void)	{ a += sm_sqr (b) ; } } ;
		class Apa_sqr_BsC		{ CALC_3_1(void)	{ a += sm_sqr (b-c) ; } } ;
		class Apa_sqrt_B	{ CALC_2_1(void)	{ a += ::sqrt (b) ; } } ;
		class Apa_sqr_CdB	{ CALC_3_1(void)	{ a += (TA) sm_sqr (c / b) ; }	} ;
		class BdaC_Apa_sqr_B{ CALC_3_2(void)	{ b /= c ; a += (TA) sm_sqr (b) ; }	} ;
		class Apa_sq_BsC	{ CALC_3_1(void)	{ a += sm_sqr (b - c) ; }	} ;		//	calculation of var
		class ApaC_Bpa_sq_C	{ CALC_3_2(void)	{ a += c; b += sm_sqr (c) ; }	} ;	//	calculation of var using steiner

//		class ApaBmB		{ CALC_2_1(void)	{ a += sm_sqr (b) ; } } ;
//		class ApaBmB_BmaC	{ CALC_3_2(void)	{ a += sm_sqr (b) ; b *= c ; } } ;
//		class ApaBmBpCmC_BmaC	{ CALC_3_2(void)	{ a += sm_sqr (b) + sm_sqr (c) ; b *= c ; } } ;

		class ApaBmB		{ CALC_2_1(void)	{ a += ::r4c_pow (b, 2.0) ; } } ;
		class ApaBmB_BmaC	{ CALC_3_2(void)	{ a += ::r4c_pow (b, 2.0) ; b *= c ; } } ;
		class ApaBmBpCmC_BmaC	{ CALC_3_2(void)	{ a += ::r4c_pow (b, 2.0) + ::r4c_pow (c, 2.0) ; b *= c ; } } ;

		class AmaB_BmaC		{ CALC_3_2(void)	{ a *= b ; b *= c; } } ;
		class AmaBmD_BmaC	{ CALC_4_2(void )	{ a *= b * d ; b *= c; } } ;

		class ApaBmB_BpaCmD { CALC_4_2(void )	{ a += sm_sqr (b); b += c * d ; } } ;
	

		class inc_a_if_b				{ CALC_2_1(void) { if (b) a += 1 ; } } ;
		class inc_a_if_b_equals_c		{ CALC_3_1(void) { if (b == (TB) c) a += 1 ; } } ;
		class inc_a_if_b_less_c			{ CALC_3_1(void) { if (b < (TB) c) a += 1 ; } } ;
		class inc_a_if_b_lesseq_c		{ CALC_3_1(void) { if (b <= (TB) c) a += 1 ; } } ;
		class inc_a_if_b_greater_c		{ CALC_3_1(void) { if (b > (TB) c) a += 1 ; } } ;
		class inc_a_if_b_greatereq_c	{ CALC_3_1(void) { if (b >= (TB) c) a += 1 ; } } ;
		class a_max_idx					{ CALC_2_1(void) { if (a < (t_size) b) a = (t_size) b ; } } ;
	} ;

	template<class F>
	class EO
	{
	public:

		template <class TA, class TB, class TC>
		static void MMcVc_R (const SVMat<TA> &a, const SCMat<TB> &b, const SCData<TC> &c)
		{
			THROW (a.nrow () == c.size ()) ;
			a.Require (b) ;
			MMcVc_NC (a, b, c) ;
		}

		template <class TA, class TB, class TC>
		static void MMcVc (const SVMat<TA> &a, const SCMat<TB> &b, const SCData<TC> &c)
		{
			THROW (a.EqualDims (b)) ;
			THROW (a.nrow () == c.size ()) ;
			MMcVc_NC (a, b, c) ;
		}

		template <class TA, class TB, class TC>
		static void MMcVc_NC (const SVMat<TA> &a, const SCMat<TB> &b, const SCData<TC> &c)
		{
			ASSERT (a.EqualDims (b)) ;
			ASSERT (a.nrow () == c.size ()) ;

			TA * pA = a ;
			TA * const pEndA = a.GetDataEnd () ;
			TB const * pB = b ;
			TC const * pC ;
			TC const * const pEndC = c.GetDataEnd () ;

			while (pA < pEndA)
			{
				pC = c ;
				while (pC < pEndC)
				{
					F::Calc (*pA, *pB, *pC) ;
					++pA ;
					++pB ;
					++pC ;
				}
			}
		}

		template <class TA, class TB, class TC>
		static void MMcVct_R (const SVMat<TA> &a, const SCMat<TB> &b, const SCData<TC> &c)
		{
			THROW (a.ncol () == c.size ()) ;
			a.Require (b) ;
			MMcVct_NC (a, b, c) ;
		}

		template <class TA, class TB, class TC>
		static void MMcVct (const SVMat<TA> &a, const SCMat<TB> &b, const SCData<TC> &c)
		{
			THROW (a.EqualDims (b)) ;
			THROW (a.ncol () == c.size ()) ;
			MMcVct_NC (a, b, c) ;
		}

		template <class TA, class TB, class TC>
		static void MMcVct_NC (const SVMat<TA> &a, const SCMat<TB> &b, const SCData<TC> &c)
		{
			ASSERT (a.EqualDims (b)) ;
			ASSERT (a.ncol () == c.size ()) ;

			TA * pA = a, * pColEndA ;
			TA * const pEndA = a.GetDataEnd () ;
			TB const * pB = b ;
			TC const * pC = c ;

			while (pA < pEndA)
			{
				pColEndA = pA + a.GetColInc () ;
				while (pA < pColEndA)
				{
					F::Calc (*pA, *pB, *pC) ;
					++pA ;
					++pB ;
				}
				++pC ;
			}
		}

		template <class TA, class TB, class TC, class TD>
		static void MVMcVct (const SVMat<TA> &a, const SVData<TB> &b, const SCMat<TC> &c, const SCData<TD> &d)
		{
			THROW (a.EqualDims (c)) ;
			THROW (b.size () == c.nrow ()) ;
			THROW (d.size () == c.ncol ()) ;
			MVMcVct_NC (a, b, c, d) ;
		}

		template <class TA, class TB, class TC, class TD>
		static void MVMcVct_NC (const SVMat<TA> &a, const SVData<TB> &b, const SCMat<TC> &c, const SCData<TD> &d)
		{
			ASSERT (a.EqualDims (c)) ;
			ASSERT (b.size () == c.nrow ()) ;
			ASSERT (d.size () == c.ncol ()) ;

			TA * pA = a, * const pEndA = a.GetDataEnd () ;
			TB * pB, * const pStartB = b, *const pEndB = b.GetDataEnd () ;
			const TC * pC = c ;
			const TD * pD = d ;

			while (pA < pEndA)
			{
				pB = pStartB ;
				while (pB < pEndB)
				{
					F::Calc (*pA, *pB, *pC, *pD) ;
					++pA ;
					++pB ;
					++pC ;
				}

				++pD ;
			}

			
		}
		
		template <class TA, class TB, class TC>
		static void MVcVct_R (const SVMat<TA> &a, const SCData<TB> &b, const SCData<TC> &c)
		{
			a.Require (b.size (), c.size ()) ;
			MVcVct_NC (a, b, c) ;
		}

		template <class TA, class TB, class TC>
		static void MVcVct (const SVMat<TA> &a, const SCData<TB> &b, const SCData<TC> &c)
		{
			THROW (a.nrow () == b.size ()) ;
			THROW (a.ncol () == c.size ()) ;
			MVcVct_NC (a, b, c) ;
		}
		
		template <class TA, class TB, class TC>
		static void MVcVct_NC (const SVMat<TA> &a, const SCData<TB> &b, const SCData<TC> &c)
		{
			ASSERT (a.nrow () == b.size ()) ;
			ASSERT (a.ncol () == c.size ()) ;

					TA *		pA = a ;
			const	TA * const	pEndA = a.GetDataEnd () ;
			const	TB *		pB ;
			const	TB * const	pEndB = b.GetDataEnd () ;
			const	TC *		pC = c ;

			while (pA < pEndA)
			{
				pB = b ;
				while (pB < pEndB)
				{
					F::Calc (*pA, *pB, *pC) ;
					++pA ;
					++pB ;
				}
				++pC ;
			}
		}

		

		template <class TA, class TB, class TC, class TD>
		static void MsMcVcVbc_R (const SVMat<TA> &a, const SCMat<TB> &b, const SCVec<TC> &c, const SCVec<TD> &d)
		{
			THROW (a.ncol () == b.ncol ()) ;
			THROW (b.nrow () == c.size ()) ;
			THROW (b.nrow () == d.size ()) ;
			d.Require (CountTrue (d), b.ncol ()) ;

			MsMcVcVbc_NC (a, b, c, d) ;
		}


		template <class TA, class TB, class TC, class TD>
		static void MsMcVcVbc (const SVMat<TA> &a, const SCMat<TB> &b, const SCVec<TC> &c, const SCVec<TD> &d)
		{
			THROW (a.ncol () == b.ncol ()) ;
			THROW (b.nrow () == c.size ()) ;
			THROW (b.nrow () == d.size ()) ;
			THROW (CountTrue (d) == a.nrow ()) ;

			MsMcVcVbc_NC (a, b, c, d) ;
		}

		template <class TA, class TB, class TC, class TD>
		static void MsMcVcVbc_NC (const SVMat<TA> &a, const SCMat<TB> &b, const SCVec<TC> &c, const SCVec<TD> &d)
		{	//	b for binary index!
			ASSERT (a.ncol () == b.ncol ()) ;
			ASSERT (b.nrow () == c.size ()) ;
			ASSERT (b.nrow () == d.size ()) ;
			ASSERT (CountTrue (d) == a.nrow ()) ;

			TA *pA = a ;
			const TB *pB = b, * const pEndB = b.GetDataEnd () ;
			const TC * const pStartC = c ;
			const TD * const pStartD = d, * const pEndD = d.GetDataEnd () ;

			while (pB < pEndB)
			{
				const TD *pD = pStartD ;
				const TC *pC = pStartC ;

				while (pD < pEndD)
				{
					if (*pD)
					{
						F::Calc (*pA, *pB, *pC) ;
						++pA ;
					}
					++pB ;
					++pC ;
					++pD ;
				}
			}
		}



		template <class TA, class TB, class TC, class TD>
		static void MsMcVctVbc_R (const SVMat<TA> &a, const SCMat<TB> &b, const SCVec<TC> &c, const SCVec<TD> &d)
		{
			THROW (a.ncol () == b.ncol ()) ;
			THROW (b.ncol () == c.size ()) ;
			THROW (b.nrow () == d.size ()) ;
			d.Require (CountTrue (d), b.ncol ()) ;

			MsMcVctVbc_NC (a, b, c, d) ;
		}


		template <class TA, class TB, class TC, class TD>
		static void MsMcVctVbc (const SVMat<TA> &a, const SCMat<TB> &b, const SCVec<TC> &c, const SCVec<TD> &d)
		{
			THROW (a.ncol () == b.ncol ()) ;
			THROW (b.ncol () == c.size ()) ;
			THROW (b.nrow () == d.size ()) ;
			THROW (CountTrue (d) == a.nrow ()) ;

			MsMcVctVbc_NC (a, b, c, d) ;
		}

		template <class TA, class TB, class TC, class TD>
		static void MsMcVctVbc_NC (const SVMat<TA> &a, const SCMat<TB> &b, const SCVec<TC> &c, const SCVec<TD> &d)
		{	//	b for binary index!
			ASSERT (a.ncol () == b.ncol ()) ;
			ASSERT (b.ncol () == c.size ()) ;
			ASSERT (b.nrow () == d.size ()) ;
			ASSERT (CountTrue (d) == a.nrow ()) ;

			TA *pA = a ;
			const TB *pB = b, * const pEndB = b.GetDataEnd () ;
			const TC *pC = c ;
			const TD * const pStartD = d, * const pEndD = d.GetDataEnd () ;

			while (pB < pEndB)
			{
				const TD *pD = pStartD ;

				while (pD < pEndD)
				{
					if (*pD)
					{
						F::Calc (*pA, *pB, *pC) ;
						++pA ;
					}
					++pB ;
					++pD ;
				}
				++pC ;
			}
		}

		template <class TA, class TB, class TC, class TD, class TE>	
		static void MsVetMcdScgVceg_R (SVMat<TA> &a, const SVData<TB> &b, const SCMat<TC> &c, const TD &d, const SCData<TE> &e)
		{	//	g for group. 
			THROW (c.ncol () == b.nsize ()) ;
			THROW (c.nrow () == e.size ()) ;
			a.Require (CountMatches (e, d), c.ncol ()) ;

			MsVetMcdScgVceg_NC (a, b, c, d, e) ;
		}

		template <class TA, class TB, class TC, class TD, class TE>	
		static void MsVetMcdScgVceg (const SVMat<TA> &a, const SVData<TB> &b, const SCMat<TC> &c, const TD &d, const SCData<TE> &e)
		{	//	g for group. 

			THROW (a.ncol () == c.ncol ()) ;
			THROW (a.ncol () == b.size ()) ;
			THROW (c.nrow () == e.size ()) ;
			THROW (CountMatches (e, d) == a.nrow ()) ;

			MsVetMcdScgVceg_NC (a, b, c, d, e) ;
		}

		template <class TA, class TB, class TC, class TD, class TE>	
		static void MsVetMcdScgVceg_NC (const SVMat<TA> &a, const SVData<TB> &b, const SCMat<TC> &c, const TD &d, const SCData<TE> &e)
		{	//	g for group. 

			ASSERT (a.ncol () == c.ncol ()) ;
			ASSERT (a.ncol () == b.size ()) ;
			ASSERT (c.nrow () == e.size ()) ;
			ASSERT (CountMatches (e, d) == a.nrow ()) ;

			//t_size dwColIncA = a.GetColInc () ;
			TA * pA = a ;
			TA const * const pEndA = a.GetDataEnd () ;
			TB * pB = b ;
			TC const * pC = c ;
			TE const * pE ;
			TE const * pEndE = e.GetDataEnd () ;

			for (; pA < pEndA; )	//	for each column of A
			{
				pE = e ;
				while (pE < pEndE)
				{
					if (d == (TD) *pE)
					{
						F::Calc (*pA, *pB, *pC) ;
						++pA ;
					}
					++pC ;
					++pE ;
				}
				++pB ;

			}
		}

		template <class TA, class TB, class TC, class TD>	
		static void MsVetMcdVcei_R (const SVMat<TA> &a, const SVData<TB> &b, const SCMat<TC> &c, const SCData<TD> &d)
		{	//	Main dimension: matrix c
			//	Constant index vector d of size c.nrow () 
			//	Matrix a refers to the columns and subsetted rows of matrix c
			//	Vector d is of size c.ncol () 
			//	F.Calc is executed for each element of c[,d]

			THROW (b.size () == c.ncol ()) ;
			THROW (getMaxIdx (d) < c.nrow ()) ;
			a.Require (d.size (), c.ncol ()) ;

			MsVetMcdVcei_NC (a, b, c, d) ;
		}

		template <class TA, class TB, class TC, class TD>	
		static void MsVetMcdVcei (const SVMat<TA> &a, const SVData<TB> &b, const SCMat<TC> &c, const SCData<TD> &d)
		{	//	Main dimension: matrix c
			//	Constant index vector d of size c.nrow () 
			//	Matrix a refers to the columns and subsetted rows of matrix c
			//	Vector d is of size c.ncol () 
			//	F.Calc is executed for each element of c[,d]

			THROW (a.ncol () == c.ncol ()) ;
			THROW (a.nrow () == d.size ()) ;
			THROW (b.size () == c.ncol ()) ;
			THROW (getMaxIdx (d) < c.nrow ()) ;
			//2do: check array d for max index!

			MsVetMcdVcei_NC (a, b, c, d) ;
		}

		template <class TA, class TB, class TC, class TD>	
		static void MsVetMcdVcei_NC (const SVMat<TA> &a, const SVData<TB> &b, const SCMat<TC> &c, const SCData<TD> &d)
		{
			ASSERT (a.ncol () == c.ncol ()) ;
			ASSERT (a.nrow () == d.size ()) ;
			ASSERT (b.size () == c.ncol ()) ;
			ASSERT (getMaxIdx (d) < c.nrow ()) ;

			//t_size dwColIncA = a.GetColInc () ;
			t_size dwColIncC = c.GetColInc () ;
			TA * pA = a ;
			TA const * const pEndA = a.GetDataEnd () ;
			TB * pB = b ;
			TC const * pC = c ;
			TD const * pD ;
			TD const * pEndD = d.GetDataEnd () ;

			for (; pA < pEndA; )
			{
				for (pD = d; pD < pEndD; ++pD)
				{
					F::Calc (*pA, *pB, pC [(t_size) *pD]) ;
					++pA ;
				}
				++pB ;
				pC += dwColIncC ;
			}
		}

		template <class TA, class TB>
		static void MVcet (const SVMat<TA> &a, const SCData<TB> &b)
		{
			THROW (a.ncol () == b.size ()) ;
			MVcet_NC (a, b) ;
		}

		template <class TA, class TB>
		static void MVcet_NC (const SVMat<TA> &a, const SCData<TB> &b)
		{
			ASSERT (a.ncol () == b.size ()) ;

			TA *pA = a ;
			TA const * const pEndA = a.GetDataEnd ();
			TA const * pColEndA ;
			const t_size dwColIncA = a.GetColInc () ;
			TB const * pB = b ;

			for (; pA < pEndA; ++pB)
				for (pColEndA = pA + dwColIncA; pA < pColEndA; ++pA)
					F::Calc (*pA, *pB) ;
		}
		

		template <class TA, class TB, class TC>							//	matrix impl
		static inline void SSVc (TA &a, TB &b, const SCData<TC> &c)
		{
			SSVc_raw (a, b, c.GetData (), c.GetDataEnd ()) ;
		}

		template <class TA, class TB, class TC>							//	raw impl.
		static void SSVc_raw (TA &a, TB &b, const TC * pC, TC const * const pEndC)
		{
			while (pC < pEndC)
			{
				F::Calc (a, b, *pC) ;
				++ pC ;
			}
		}

		template <class TA, class TB, class TC, class TD>				//	matrix impl
		static inline void SSVcVc (TA &a, TB &b, const SCData<TC> &c, const SCData<TD> &d)
		{
			THROW (c.size () == d.size ()) ;
			SSVcVc_NC (a, b, c, d) ;
		}
		
		template <class TA, class TB, class TC, class TD>				//	matrix impl
		static inline void SSVcVc_NC (TA &a, TB &b, const SCData<TC> &c, const SCData<TD> &d)
		{
			ASSERT (c.size () == d.size ()) ;
			SSVcVc_raw (a, b, c.GetData (), c.GetDataEnd (), d.GetData ()) ;
		}

		template <class TA, class TB, class TC, class TD>				//	raw impl.
		static void SSVcVc_raw (TA &a, TB &b, const TC * pC, TC const * const pEndC, const TD * pD)
		{
			while (pC < pEndC)
			{
				F::Calc (a, b, *pC, *pD) ;
				++ pC ;
				++ pD ;
			}
		}

		template <class TA, class TB, class TC, class TD, class TE>
		static inline void SScScVcVc_NC (TA &a, const TB &b, const TC &c, const SCData<TD> &d, const SCData<TE> &e)
		{
			ASSERT (d.size () == e.size ()) ;

			const TD *pD = d, * const pEndD = d.GetDataEnd () ;
			const TE *pE = e ;
			while (pD < pEndD)
			{
				F::Calc (a, b, c, *pD, *pE) ;
				++ pD ;
				++ pE ;
			}
		}

		template <class TA, class TB, class TC>							//	matrix impl
		static inline void SScVc (TA &a, const TB &b, const SCData<TC> &c)
		{
			SScVc_raw (a, b, c.GetData (), c.GetDataEnd ()) ;
		}

		template <class TA, class TB, class TC>							//	raw impl.
		static void SScVc_raw (TA &a, const TB &b, TC const * pC, TC const * const pEndC)
		{
			while (pC < pEndC)
			{
				F::Calc (a, b, *pC) ;
				++ pC ;
			}
		}

		template <class TA, class TB>
		static void SVc (TA &a, const SCData<TB> &b)
		{
			SVc_raw (a, b.GetData (), b.GetDataEnd ()) ;
		}

		template <class TA, class TB>
		static void SVc_raw (TA &a, TB const * pB, TB const * const pEndB)
		{
			for  (; pB < pEndB; ++pB)
				F::Calc (a, *pB) ;
		}

		template <class TA, class TB>
		static void SV (TA &a, SCData<TB> &b)
		{
			TB *pB ;
			TB * const pEnd = b.GetDataEnd () ;

			for  (pB = b; pB < pEnd; ++pB)
				F::Calc (a, *pB) ;
		}

		template <class TA, class TB, class TC>							//	matrix impl
		static inline void SVcVc (TA &a, const SCData<TB> &b, const SCData<TC> &c)
		{
			THROW (b.size () == c.size ()) ;
			SVcVc_NC (a, b, c) ;
		}

		template <class TA, class TB, class TC>							//	matrix impl
		static inline void SVcVc_NC (TA &a, const SCData<TB> &b, const SCData<TC> &c)
		{
			ASSERT (b.size () == c.size ()) ;
			SVcVc_raw (a, b.GetData (), b.GetDataEnd (), c.GetData ()) ;
		}

		template <class TA, class TB, class TC>							//	raw impl.
		static void SVcVc_raw (TA &a, TB const * pB, TB const * const pEndB, TC const * pC)
		{
			while (pB < pEndB)
			{
				F::Calc (a, *pB, *pC) ;
				++ pB ;
				++ pC ;
			}
		}

		template <class TA, class TB, class TC>							//	raw impl.
		static void SVSc (TA &a, const SVData<TB> &b, const TC &c)
		{
			TB *pB = b, * const pEndB = b.GetDataEnd () ;

			for (; pB < pEndB; ++pB)
				F::Calc (a, *pB, c) ;
		}


		template <class TA>
		static void V (const SVData<TA> &a)
		{
			TA *pa = a.GetData () ;
			const TA *const paEnd = a.GetDataEnd () ;

			for (; pa < paEnd; pa++)
				F::Calc (*pa) ;
		}

		template <class TA, class TB, class TC>							//	matrix impl
		static inline void SVVc (TA &a, const SVData<TB> &b, const SCData<TC> &c)
		{
			THROW (b.size () == c.size ()) ;
			SVVc_NC (a, b, c) ;
		}

		template <class TA, class TB, class TC>							//	matrix impl
		static inline void SVVc_NC (TA &a, const SVData<TB> &b, const SCData<TC> &c)
		{
			ASSERT (b.size () == c.size ()) ;
			SVVc_raw (a, b.GetData (), b.GetDataEnd (), c.GetData ()) ;
		}

		template <class TA, class TB, class TC>							//	raw impl.
		static void SVVc_raw (TA &a, TB * pB, TB * const pEndB, TC const * pC)
		{
			while (pB < pEndB)
			{
				F::Calc (a, *pB, *pC) ;
				++ pB ;
				++ pC ;
			}
		}

		template <class TB, class TC, class TD, class TE>	
		static void VetMcdScgVceg (const SVData<TB> &b, const SCMat<TC> &c, const TD &d, const SCData<TE> &e)
		{	//	g for group. 

			THROW (c.nrow () == e.size ()) ;

			VetMcdScgVceg_NC (b, c, d, e) ;
		}

		template <class TB, class TC, class TD, class TE>	
		static void VetMcdScgVceg_NC (const SVData<TB> &b, const SCMat<TC> &c, const TD &d, const SCData<TE> &e)
		{	//	g for group. 

			ASSERT (c.nrow () == e.size ()) ;

			//t_size dwColIncA = a.GetColInc () ;
			TB * const pEndB = b.GetDataEnd () ;
			TB * pB = b ;
			TC const * pC = c ;
			TE const * pE ;
			TE const * pEndE = e.GetDataEnd () ;

			for (; pB < pEndB; )	//	for each column of A
			{
				pE = e ;
				while (pE < pEndE)
				{
					if (d == (TD) *pE)
						F::Calc (*pB, *pC) ;
					++pC ;
					++pE ;
				}
				++pB ;

			}
		}

		template <class TB, class TC, class TD>	
		static void VetMcdVcei (const SVData<TB> &b, const SCMat<TC> &c, const SCData<TD> &d)
		{	//	Main dimension: matrix c
			//	Constant index vector d of size c.nrow () 
			//	Vector d is of size c.ncol () 
			//	F.Calc is executed for each element of c[,d]

			THROW (b.size () == c.ncol ()) ;
			THROW (getMaxIdx (d) < c.nrow ()) ;
			//2do: check array d for max index!

			VetMcdVcei_NC (b, c, d) ;
		}

		template <class TB, class TC, class TD>	
		static void VetMcdVcei_NC (const SVData<TB> &b, const SCMat<TC> &c, const SCData<TD> &d)
		{
			ASSERT (b.size () == c.ncol ()) ;
			ASSERT (getMaxIdx (d) < c.nrow ()) ;

			//t_size dwColIncA = a.GetColInc () ;
			t_size dwColIncC = c.GetColInc () ;

			TB * pB = b ;
			TB const * pEndB = b.GetDataEnd () ;

			TC const * pC = c ;
			TD const * pD ;
			TD const * pEndD = d.GetDataEnd () ;

			for (; pB < pEndB; )
			{
				for (pD = d; pD < pEndD; ++pD)
					F::Calc (*pB, pC [(t_size) *pD]) ;
				++pB ;
				pC += dwColIncC ;
			}
		}



		template <class TA, class TB, class TC>
		static void VVcVc (const SVData<TA> &a, const SCData<TB> &b, const SCData<TC> &c)
		{
			THROW (a.size  () == b.size ()) ;
			THROW (a.size  () == c.size ()) ;
			VVcVc_NC (a, b, c) ;
		}

		template <class TA, class TB, class TC>
		static void VVcVc_NC (const SVData<TA> &a, const SCData<TB> &b, const SCData<TC> &c)
		{
			ASSERT (a.size  () == b.size ()) ;
			ASSERT (a.size  () == c.size ()) ;

			const TA *pEndA = a.GetDataEnd () ;
			TA *pA = a.GetData () ;
			const TB * pB = b.GetData () ;
			const TC * pC = c.GetData () ;

			while (pA < pEndA)
			{
				F::Calc (*pA, *pB, *pC) ;
				++pA ;
				++pB ;
				++pC ;
			}
		}

		template <class TA, class TB>
		static void VMc (const SVData<TA> &a, const SCMat<TB> &b)
		{
			THROW (a.size () == b.nrow ()) ;

			VMc_NC (a, b) ;
		}

		template <class TA, class TB>
		static void VMc_NC (const SVData<TA> &a, const SCMat<TB> &b)
		{
			ASSERT (a.size () == b.nrow ()) ;

			TA * pA, * const pStartA = a, * const pEndA = a.GetDataEnd () ;
			const TB * pB = b, * const pEndB = b.GetDataEnd () ;

			while (pB < pEndB)
			{
				pA  = pStartA ;
				while (pA < pEndA)
				{
					F::Calc (*pA, *pB) ;
					++pA ;
					++pB ;
				}
			}
		}



		template <class TA, class TB, class TC>
		static void VMcVct (const SVData<TA> &a, const SCMat<TB> &b, const SCData<TC> &c)
		{
			THROW (a.size () == b.nrow ()) ;
			THROW (c.size () == b.ncol ()) ;

			VMcVct_NC (a, b, c) ;
		}

		template <class TA, class TB, class TC>
		static void VMcVct_NC (const SVData<TA> &a, const SCMat<TB> &b, const SCData<TC> &c)
		{
			ASSERT (a.size () == b.nrow ()) ;
			ASSERT (c.size () == b.ncol ()) ;

			TA * const  pStartA = a, * const pEndA = a.GetDataEnd () ;
			const TB * pB = b, * const pEndB = b.GetDataEnd () ;
			const TC * pC = c ;

			while (pB < pEndB)
			{
				TA *pA = pStartA ;
				while (pA < pEndA)
				{
					F::Calc (*pA, *pB, *pC) ;
					++pA ;
					++pB ;
				}
				++pC ;
			}
		}

		template <class TA, class TB>
		static inline void VetMcd (const SVData<TA> &a, const SCMat<TB> &b)
		{
			THROW (a.size () == b.ncols ()) ;
			VetMcd_NC (a, b) ;
		}

		template <class TA, class TB>
		static inline void VetMcd_NC (const SVData<TA> &a, const SCMat<TB> &b)
		{
			ASSERT (a.size () == b.ncol ()) ;

			TA *pA = a ;
			TA *const pEndA = a.GetDataEnd () ;
			const TB *pbColEnd, *pB = b ;

			for ( ; pA < pEndA; ++pA)
				for (pbColEnd = pB + b.nrow (); pB < pbColEnd; ++pB)
					F::Calc (*pA, *pB) ;
		}


		template <class TA, class TB, class TC>
		static void VtMcVc (const SVData<TA> &a, const SCMat<TB> &b, const SCData<TC> &c)
		{
			THROW (a.size () == b.ncol ()) ;
			THROW (c.size () == b.nrow ()) ;

			VtMcVc_NC (a, b, c) ;
		}

		template <class TA, class TB, class TC>
		static void VtMcVc_NC (const SVData<TA> &a, const SCMat<TB> &b, const SCData<TC> &c)
		{
			ASSERT (a.size () == b.ncol ()) ;
			ASSERT (c.size () == b.nrow ()) ;

					TA	*		pA = a ;
			const	TB	*		pB = b ;
			const	TC	*		pC = c ;
			const	TB 	* const pEndB = b.GetDataEnd () ;
			const	TC  * const pEndC = c.GetDataEnd () ;

			while (pB < pEndB)
			{
				pC = c ;
				while (pC < pEndC)
				{
					F::Calc (*pA, *pB, *pC) ;
					++pB ;
					++pC ;
				}
				++pA ;
			}
		}


		template <class TA, class TB>
		static void VSc (const SVData<TA> &a, const TB &b)
		{
			TA *		pA = a ;
			TA * const	pEnd = a.GetDataEnd () ;

			while (pA < pEnd)
			{
				F::Calc (*pA, b) ;
				++pA ;
			}
		}

		template <class TA, class TB, class TC>
		static void VScSc (const SVData<TA> &a, const TB &b, const TC &c)
		{
			TA *		pA = a ;
			TA * const	pEndA = a.GetDataEnd () ;
			while (pA < pEndA)
			{
				F::Calc (*pA, b, c) ;
				++ pA ;
			}
		}

		template <class TA, class TB, class TC, class TD>
		static void VScScVc (const SVData<TA> &a, const TB &b, const TC &c, const SCData<TD> &d)
		{
			THROW (a.size () == d.size ()) ;
			VScScVc_NC (a, b, c, d) ;
		}

		template <class TA, class TB, class TC, class TD>
		static void VScScVc_NC (const SVData<TA> &a, const TB &b, const TC &c, const SCData<TD> &d)
		{
			ASSERT (a.size () == d.size ()) ;

			TA *  pA = a, * const pEndA = a.GetDataEnd () ;
			TD const * pD = d ;

			while (pA < pEndA)
			{
				F::Calc (*pA, b, c, *pD) ;
				++ pA ;
				++ pD ;
			}
		}

		
		template <class TA, class TB, class TC, class TD>
		static void SVScVc (TA &a, const SVData<TB> &b, const TC &c, const SCData<TD> &d)
		{
			THROW (b.size () == d.size ()) ;
			SVScVc_NC (a, b, c, d) ;
		}

		template <class TA, class TB, class TC, class TD>
		static void SVScVc_NC (TA &a, const SVData<TB> &b, const TC &c, const SCData<TD> &d)
		{
			ASSERT (b.size () == d.size ()) ;
			TB *		pB = b, * const	pEndB = b.GetDataEnd () ;
			TD const *	pD = d ;

			while (pB < pEndB)
			{
				F::Calc (a, *pB, c, *pD) ;
				++pB ;
				++pD ;
			}
		}

		template <class TA, class TB, class TC>
		static void VScVc (const SVData<TA> &a, const TB &b, const SCData<TC> &c)
		{
			THROW (a.size () == c.size ()) ;
			VScVc_NC (a, b, c) ;
		}

		template <class TA, class TB, class TC>
		static void VScVc_NC (const SVData<TA> &a, const TB &b, const SCData<TC> &c)
		{
			ASSERT (a.size () == c.size ()) ;
			TA *		pA = a ;
			TA * const	pEndA = a.GetDataEnd () ;
			TC const *	pC = c ;

			while (pA < pEndA)
			{
				F::Calc (*pA, b, *pC) ;
				++pA ;
				++pC ;
			}
		}

		template <class TA, class TB>
		static void VtMc (const SVData<TA> &a, const SCMat<TB> &b)
		{
			THROW (a.size () == b.ncol ()) ;
			VtMc_NC (a, b) ;
		}

		template <class TA, class TB>
		static void VtMc_NC (const SVData<TA> &a, const SCMat<TB> &b)
		{
			ASSERT (a.size () == b.ncol ()) ;
			
			TA *pA = a, * const pEndA = a.GetDataEnd () ;
			const TB *pB = b ;

			while (pA < pEndA)
			{
				TB const * const pColEndB = pB + b.GetColInc () ;
				while (pB < pColEndB)
				{
					F::Calc (*pA, *pB) ;
					++pB ;
				}
				++pA ;
			}
		}



		template <class TA, class TB, class TC, class TD>
		static void VtMcVcVc (const SVData<TA> &a, const SCMat<TB> &b, const SCData<TC> &c, const SCData<TD> &d)
		{
			THROW (a.size () == b.ncol ()) ;
			THROW (c.size () == b.nrow ()) ;
			THROW (d.size () == b.nrow ()) ;
			VtMcVcVc_NC (a, b, c, d) ;
		}

		template <class TA, class TB, class TC, class TD>
		static void VtMcVcVc_NC (const SVData<TA> &a, const SCMat<TB> &b, const SCData<TC> &c, const SCData<TD> &d)
		{
			ASSERT (a.size () == b.ncol ()) ;
			ASSERT (c.size () == b.nrow ()) ;
			ASSERT (d.size () == b.nrow ()) ;
			
			TA *pA = a, * const pEndA = a.GetDataEnd () ;
			const TB *pB = b ;
			const TC *pC, * const pStartC = c, * const pEndC = c.GetDataEnd () ;
			const TD *pD, * const pStartD = d ;

			while (pA < pEndA)
			{
				pC = pStartC ;
				pD = pStartD ;
				while (pC < pEndC)
				{
					F::Calc (*pA, *pB, *pC, *pD) ;
					++pB ;
					++pC ;
					++pD ;
				}
				++pA ;
			}

		}


		template <class TA, class TB, class TC>
		static void VsVcVbc (const SVData <TA> &a, const SCData<TB> &b, const SCData<TC> &c)
		{
			THROW (b.size () == c.size ()) ;
			THROW (CountTrue (c) == a.size ()) ;
			VsVcVbc_NC (a, b, c) ;
		}

		template <class TA, class TB, class TC>
		static void VsVcVbc_NC (const SVData <TA> &a, const SCData<TB> &b, const SCData<TC> &c)
		{
			ASSERT (b.size () == c.size ()) ;
			ASSERT (CountTrue (c) == a.size ()) ;
			
			TA *pA = a ;
			const TB *pB  = b, * const pEndB = b.GetDataEnd () ;
			const TC *pC = c ;

			while (pB < pEndB)
			{
				if (*pC)
				{
					F::Calc (*pA, *pB) ;
					++ pA ;
				}
				++ pB ;
				++ pC ;
			}

		}

		template <class TA, class TB>
		static void VVc (const SVData <TA> &a, const SCData<TB> &b)
		{
			THROW (a.size () == b.size ()) ;
			VVc_NC (a, b) ;
		}

		template <class TA, class TB>
		static void VVc_NC (const SVData <TA> &a, const SCData<TB> &b)
		{
			ASSERT (a.size () == b.size ()) ;

			VVc_raw (a.GetData (), a.GetDataEnd (), b.GetData ()) ;
		}

		template <class TA, class TB>
		static void VVc_raw (TA *pdA, TA * const pdEndA, TB const *pdB)
		{
			while (pdA < pdEndA)
			{
				F::Calc (*pdA, *pdB) ;				
				++pdA ;
				++pdB ;
			}
		}

		template <class TA, class TB, class TC>
		static void VVSc (const SVData <TA> &a, const SVData<TB> &b, const TC &c)
		{
			THROW (a.size () == b.size ()) ;
			VVSc_NC (a, b, c) ;
		}
		
		template <class TA, class TB, class TC>
		static void VVSc_NC (const SVData <TA> &a, const SVData<TB> &b, const TC &c)
		{
			ASSERT (a.size () == b.size ()) ;

			TA *pA = a ;
			TB *pB  = b, * const pEndB = b.GetDataEnd () ;

			while (pB < pEndB)
			{
				F::Calc (*pA, *pB, c) ;
				++ pA ;
				++ pB ;
			}

		}


		template <class TA, class TB, class TC, class TD>
		static void VVScSc (const SVData <TA> &a, const SVData<TB> &b, const TC &c, const TD &d)
		{
			THROW (a.size () == b.size ()) ;
			VVScSc_NC (a, b, c, d) ;
		}
		
		template <class TA, class TB, class TC, class TD>
		static void VVScSc_NC (const SVData <TA> &a, const SVData<TB> &b, const TC &c, const TD &d)
		{
			ASSERT (a.size () == b.size ()) ;

			TA *pA = a ;
			TB *pB  = b, * const pEndB = b.GetDataEnd () ;

			while (pB < pEndB)
			{
				F::Calc (*pA, *pB, c, d) ;
				++ pA ;
				++ pB ;
			}

		}

		
		
		
		


		template <class TA>
		static const TA &Vc_transitive (const SCData <TA> &a)
		{
			const TA *			pA = a ;
			const TA * const	pEndA = a.GetDataEnd () ;
			const TA *			pRet = pA ;

			while (++pA < pEndA)
				if (F::Calc (*pA, *pRet))
					pRet = pA ;

			return *pRet ;
		}

		template <class TA>
		static TA &V_transitive (const SVData <TA> &a)
		{
					TA	*		pA		= a ;
			const	TA 	* const	pEndA	= a.GetDataEnd () ;
					TA	*		pRet	= pA ;

			while (++pA < pEndA)
				if (F::Calc (*pA, *pRet))
					pRet = pA ;

			return *pRet ;
		}

		template <class TA>
		static void V_pairs (const SVData <TA> &a)
		{
					TA	*		pA		= a ;
			const	TA 	* const	pEndA	= a.GetDataEnd () ;

			while (++pA < pEndA)
				F::Calc (pA[0], pA[-1]) ;
		}

		template <class TA>
		static void V_pairs_r (const SVData <TA> &a)
		{	//	computes each pair in reverse order. size = 10: (8, 9); (7, 8); ...; (0, 1)
					TA	*		pA		= a.GetDataEnd () ;
			const	TA 	* const	pEndA	= a ;

			while (--pA > pEndA)
				F::Calc (pA[-1], pA[0]) ;
		}


	} ;

#endif	//	#ifndef SMAT_ELOP_H
