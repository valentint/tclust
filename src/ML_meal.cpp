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

//	ML_meal.cpp
//	MATLAB Mathematical Environment Abstraction Layer

#ifdef MATLAB_MEX_FILE
#include "ML_meal.h"

#ifdef ES_DEV_ENV
	#include "../SMat/smat.h"
#else
	#include "smat.h"
#endif

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


	void meal_geev (const char* jobvl, const char* jobvr, const int* n, double* a, const int* lda, double* wr, double* wi, double* vl, const int* ldvl, double* vr, const int* ldvr, double* work, const int* lwork, int* info)
	{ FORTRAN_WRAPPER(dgeev)((char *) jobvl, (char *) jobvr, (int *) n, a, (int *) lda, wr, wi, vl, (int *) ldvl, vr, (int *) ldvr, work, (int *) lwork, info) ; }

	void meal_gemm (const char *transa, const char *transb, const int *m, const int *n, const int *k, const double *alpha, const double *a, const int *lda, const double *b, const int *ldb, const double *beta, double *c, const int *ldc)
	{ 
#ifdef OLD_ML
		THROW (0) ;
#else
		FORTRAN_WRAPPER(dgemm)((char *) transa, (char *)transb, (int *)m, (int *) n, (int *) k, (double *) alpha, (double *) a, (int *) lda, (double *) b, (int *) ldb, (double *) beta, (double *) c, (int *) ldc) ; 
#endif
	}

	void meal_dgesv (const int* n, const int* nrhs, double* a, const int* lda, int* ipiv, double* b, const int* ldb, int* info)
	{ FORTRAN_WRAPPER(dgesv)((int *) n, (int *) nrhs, a, (int *) lda, ipiv, b, (int *) ldb, info) ; }

	void meal_dgesvd (const char* jobu, const char* jobvt, const int* m, const int* n, double* a, const int* lda, double* s, double* u, const int* ldu, double* vt, const int* ldvt, double* work, const int* lwork, int* info)
	{ FORTRAN_WRAPPER(dgesvd)((char *) jobu, (char *) jobvt, (int *)m, (int *)n, a, (int *)lda, s, u, (int *)ldu, vt, (int *)ldvt, work, (int *)lwork, info) ; }

	void meal_sort (double *d, int l)
	{
		sme_qsort (d, l) ;
	}

	void meal_sort_order (double *d, int *o, int l)
	{
		sme_qsortI (d, o, l) ;
	}

	void meal_sort_order_rev (double *d, int *o, int l)
	{
		sme_qsortI (d, o, l, TRUE) ;
	}
/*
////////////////////////
//	Random Generator  //
////////////////////////

	void meal_PutRNGstate ()	{ THROW(0) ; }
	void meal_GetRNGstate ()	{ THROW(0) ; 	}

	double meal_unif_rand ()	{ THROW(0) ; return 0 ;	}
	double meal_norm_rand ()	{ THROW(0) ; return 0 ;	}
	double meal_exp_rand  ()	{ THROW(0) ; return 0 ;	}
*/
////////////////////////////////////
//	special values amd constants  //
////////////////////////////////////

	double	meal_NaN	() { return mxGetNaN () ; }
	double	meal_PosInf	() { return mxGetInf () ; }
	double	meal_NegInf	() { return -mxGetInf (); }
	double	meal_NaReal	() { THROW (0); return 0 ; }
	int		meal_NaInt	() { THROW (0); return 0 ; }

	double  meal_PI		() { return 3.141592653589793238462643383279502884197169399375 ; } //	{ return utGetPI () ; }	//	where did utGetPI go?

//////////////////////////
//	printing functions  //
//////////////////////////

	void meal_printf (const char *sz, ...)
	{
//		va_list va_l ;				//	further arguments not supported yet...
//		va_start (va_l, sz) ;
		mexPrintf (sz/*, va_l*/) ;
	}

	void meal_warning (const char *sz)
	{
		mexWarnMsgTxt (sz) ;
	}

	void meal_error (const char *sz)
	{
		mexErrMsgTxt (sz) ;
	}

	void *meal_alloc (size_t n, int s)
	{
		return new char [n * s] ;
	}

	void meal_free (void *p)
	{
		delete [] (char *) p ;
	}


//////////////////
//	Exceptions  //
//////////////////

	void meal_OnException (const char * szDate, const char * szFile, int nLine)
	{

		mexPrintf (
			"\n"
			"  An exception occurred.\n"
			"  Please contact the author (%s), providing\n"
			"  the following details:\n"
			"\n"
			"\tR version number\n"
			"\tPackage version number\n"
			"\tBuild date:\t%s\n"
			"\tFile:\t\t%s\n"
			"\tLine:\t\t%d\n"
			"\n"
			"  If possible please include the code which caused this error, including\n"
			"  eventual source data and the state of the random generator (seed) before\n"
			"  experiencing this issue.\n"
			"\n"
			"\tYour contribution is appreciated!\n\n",
			GetRealSettings ().GetEmail (), szDate, szFile, nLine) ;
		meal_error ("An exception has occurred.") ;
	}

	void meal_OnUException ()
	{
		mexPrintf (
			"\n"
			"  An unknown exception occurred.\n"
			"  Please contact the author (%s), providing\n"
			"  the following details:\n"
			"\n"
			"\tR version number\n"
			"\tPackage version number\n"
			"\n"
			"  If possible please include the code which caused this error, including\n"
			"  eventual source data and the state of the random generator (seed) before\n"
			"  experiencing this issue.\n\n"
			"\n"
			"\tYour contribution is appreciated!\n\n",
			GetRealSettings ().GetEmail ()) ;
		meal_error ("An unknown exception has occurred.") ;
	}

#endif	//	#ifdef MATLAB_MEX_FILE
