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

//	R.meal.cpp
//	R Mathematical Environment Abstraction Layer

#define USE_FC_LEN_T

#ifdef R_PACKAGE_FILE


#define R_USE_C99_IN_CXX

#include "R_meal.h"




/////////////////////	
//	CRmealSettings  //
/////////////////////	

	CRmealSettings &GetRealSettings ()
	{
		static CRmealSettings settings ;
		return settings ;
	}

	CRmealSettings::CRmealSettings ()
		: m_szEmail ("<NA>")
	{

	}

	CRmealSettings::CRmealSettings (const char *szEmail)
	{
		if (szEmail)	GetRealSettings ().m_szEmail = szEmail ;
	}


#include <R.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

//	LAPACK
	void meal_geev (const char* jobvl, const char* jobvr, const int* n, double* a, const int* lda, double* wr, double* wi, double* vl, const int* ldvl, double* vr, const int* ldvr, double* work, const int* lwork, int* info)
	{ F77_CALL(dgeev)(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info FCONE FCONE) ; }

	void meal_gesv (const int* n, const int* nrhs, double* a, const int* lda, int* ipiv, double* b, const int* ldb, int* info)
	{ F77_CALL(dgesv)(n, nrhs, a, lda, ipiv, b, ldb, info) ; }

	void meal_gesvd (const char* jobu, const char* jobvt, const int* m, const int* n, double* a, const int* lda, double* s, double* u, const int* ldu, double* vt, const int* ldvt, double* work, const int* lwork, int* info)
	{ F77_CALL(dgesvd)(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info FCONE FCONE) ; }


//	SORT
	void meal_sort (double *d, int l)
	{	R_qsort (d, 1, l) ;	}

	void meal_sort_order (double *d, int *o, int l)
	{
		int i ;
		for (i = l - 1; i != -1; i--)
			o[i] = i ;
		rsort_with_index(d, o, l);
	}

	void meal_sort_order_rev (double *d, int *o, int l)
	{
		int i ;
		for (i = l - 1; i != -1; i--)
			o[i] = i ;
		rsort_with_index(d, o, l);				//	2do: use r_qsort_I instead!

		double dTemp ;
		int nTemp ;

		for (i = 0, --l; i < l; ++i, --l)		//	2do: implement an reverse order - function. this should be as fast! (could be called twice)
		{
//			sm_swap (d[i], d[l], dTemp) ;		//	2do: check if this works..
//			sm_swap (o[i], o[l], nTemp) ;

			dTemp = d[i] ;
			d[i] = d[l] ;
			d[l] = dTemp ;

			nTemp = o[i]; ;
			o[i] = o[l]; ;
			o[l] = nTemp; ;
		}
	}

////////////////////////
//	Random Generator  //
////////////////////////

	void meal_PutRNGstate ()	{ PutRNGstate () ;	}
	void meal_GetRNGstate ()	{ GetRNGstate () ;	}

	double meal_unif_rand ()	{ return unif_rand () ;	}
	double meal_norm_rand ()	{ return norm_rand () ;	}
	double meal_exp_rand  ()	{ return exp_rand () ;	}

////////////////////////////////////
//	special values amd constants  //
////////////////////////////////////

	double	meal_NaN		() { return R_NaN ; }
	double	meal_PosInf	() { return R_PosInf ; }
	double	meal_NegInf	() { return R_NegInf ; }
	double	meal_NaReal	() { return R_NaReal ; }
	int		meal_NaInt	() { return R_NaInt ; }

	double  meal_PI		() { return M_PI ; }

//////////////////////////
//	printing functions  //
//////////////////////////

	void meal_printf (const char *sz, ...)
	{
		va_list va_l ;
		va_start (va_l, sz) ;
		Rvprintf (sz, va_l) ;
	}

	// VT::07.12.2023 - fix warning format string is not a string literal (potentially insecure)
    /*
	void meal_warning (const char *sz)
	{
		Rf_warning (sz) ;
	}

	void meal_error (const char *sz)
	{
		Rf_error (sz) ;
	}
    */
    
	void *meal_alloc (size_t n, int s)
	{
		return calloc (n, s) ;
	}

	void meal_free (void *p)
	{
		Free (p) ;
	}

//////////////////
//	Exceptions  //
//////////////////

	void meal_OnException (const char * szDate, const char * szFile, int nLine)
	{
		meal_printf (
			"\n"
			"  An exception occurred.\n"
			"  Please contact the author (%s), providing\n"
			"  the following information:\n"
			"\n"
			"    - The R-code which caused the problem\n"
			"    - Eventually used data sets and the state of the random generator (seed)\n"
			"    - R version number\n"
			"    - Package version number\n"
			"    - File:    %s\n"
			"    - Line:    %d\n"
			"\n"
			"  Your contribution is appreciated!\n\n",
			GetRealSettings ().GetEmail (), szFile, nLine) ;
		
        // VT::07.12.2023 - fix warning format string is not a string literal (potentially insecure)
		// meal_error ("An exception has occurred.") ;
        Rf_error("An exception has occurred.");
	}

	void meal_OnUException ()
	{
		meal_printf (
			"\n"
			"  An unknown exception occurred.\n"
			"  Please contact the author (%s), providing\n"
			"  the following details:\n"
			"\n"
			"    - The R-code which caused the problem\n"
			"    - Eventually used data sets and the state of the random generator (seed)\n"
			"    - R version number\n"
			"    - Package version number\n"
			"\n"
			"  Your contribution is appreciated!\n\n",
			GetRealSettings ().GetEmail ()) ;

		// VT::07.12.2023 - fix warning format string is not a string literal (potentially insecure)
		// meal_error ("An unknown exception has occurred.") ;
        Rf_error("An unknown exception has occurred.");
	}
#endif	//	#ifdef R_PACKAGE_FILE
