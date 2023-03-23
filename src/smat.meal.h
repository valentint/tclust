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

//	smat.meal.h
//	Mathematical Environment Abstraction Layer

#ifndef SMAT_MEAL_H
#define SMAT_MEAL_H


////////////
//  BLAS  //
////////////

//	Level 1
	double meal_dot (const int *n, const double *dx, const int *incx, const double *dy, const int *incy);
	double meal_nrm2 (const int *n, const double *dx, const int *incx);
	void meal_scal (const int *n, const double *alpha, double *dx, const int *incx) ;
	void meal_axpy(const int *n, const double *alpha, const double *dx, const int *incx, double *dy, const int *incy) ;

//	Level 2

#ifdef IMPL_BLAS_R_CONST_ERROR
	void meal_ger (const int *m, const int *n, const double *alpha, double *x, const int *incx, double *y, const int *incy, double *a, const int *lda) ;
#else
	void meal_ger (const int *m, const int *n, const double *alpha, const double *x, const int *incx, const double *y, const int *incy, double *a, const int *lda) ;
#endif

//	Level 3
	void meal_gemm (const char *transa, const char *transb, const int *m, const int *n, const int *k, const double *alpha, const double *a, const int *lda, const double *b, const int *ldb, const double *beta, double *c, const int *ldc) ;

//////////////
//  LAPACK  //
//////////////

//svd
	void meal_gesv (const int* n, const int* nrhs, double* a, const int* lda, int* ipiv, double* b, const int* ldb, int* info) ;
	void meal_gesvd (const char* jobu, const char* jobvt, const int* m, const int* n, double* a, const int* lda, double* s, double* u, const int* ldu, double* vt, const int* ldvt, double* work, const int* lwork, int* info) ;

//invert
	void meal_geev(const char* jobvl, const char* jobvr, const int* n, double* a, const int* lda, double* wr, double* wi, double* vl, const int* ldvl, double* vr, const int* ldvr, double* work, const int* lwork, int* info) ;

/////////////////////
//	Sort Routines  //
/////////////////////

	void meal_sort (double *d, int l) ;
	void meal_sort_order (double *, int *, int) ;
	void meal_sort_order_rev (double *d, int *o, int l) ;

///////////////////////////////
//  Random Number Generator  //
///////////////////////////////

	void meal_PutRNGstate () ;
	void meal_GetRNGstate () ;

	double meal_unif_rand () ;
	double meal_norm_rand ();
	double meal_exp_rand  ();

//	void meal_runif (double *d, int l) ;
//	void meal_runif (double *d, int l, double dL, double dU) ;
//	void meal_runif_r (double *d, int l) ;
//	void meal_SampleNoReplace(int k, int n, int *y, int *x) ;

////////////////////////////////////
//	special values and constants  //
////////////////////////////////////

	double	meal_NaN	() ;
	double	meal_PosInf	() ;
	double	meal_NegInf	() ;
	double	meal_NaReal	() ;
	int		meal_NaInt	() ;

	double  meal_PI () ;

//////////////////////////
//	printing functions  //
//////////////////////////

	void meal_printf (const char *, ...) ;
	void meal_warning (const char *) ;
	void meal_error (const char *) ;

//////////////////
//	Exceptions  //
//////////////////

	void meal_OnException (const char * szDate, const char * szFile, int nLine) ;
	void meal_OnUException () ;

#endif	//	#ifndef SMAT_MEAL_H
