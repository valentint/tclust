#define USE_FC_LEN_T

#ifdef R_PACKAGE_FILE

#include <R_ext/RS.h>
#include <R_ext/BLAS.h>

#include "Rversion.h"


#ifndef _MB_CONST
	#define _MB_CONST	const 
#endif

#ifndef _MB_TYPE_D
	#define _MB_TYPE_D double
#endif


#ifndef _MB_INT
	#define _MB_INT	int
#endif

#ifndef _MB_CHAR
	#define _MB_CHAR	char
#endif


//	Level 1 BLAS


	_MB_TYPE_D meal_asum (const _MB_INT *n, const _MB_TYPE_D *dx, const _MB_INT *incx)
	{
		return F77_CALL (dasum)
		(
			(_MB_CONST _MB_INT *)		n,
			(_MB_CONST _MB_TYPE_D *)	dx,
			(_MB_CONST _MB_INT *)		incx
		) ;
	}

	void meal_axpy(const _MB_INT *n, const _MB_TYPE_D *alpha, const _MB_TYPE_D *dx, const _MB_INT *incx, _MB_TYPE_D *dy, const _MB_INT *incy)
	{
		F77_CALL (daxpy)
		(
			(_MB_CONST _MB_INT *)		n,
			(_MB_CONST _MB_TYPE_D *)	alpha,
			(_MB_CONST _MB_TYPE_D *)	dx,
			(_MB_CONST _MB_INT *)		incx,
			(_MB_TYPE_D *)				dy,
			(_MB_CONST _MB_INT *)		incy
		) ; 
	}

	void meal_copy (const _MB_INT *n, const _MB_TYPE_D *dx, const _MB_INT *incx, _MB_TYPE_D *dy, const _MB_INT *incy)
	{
		F77_CALL (dcopy)
		(
			(_MB_CONST _MB_INT *)		n,
			(_MB_CONST _MB_TYPE_D *)	dx,
			(_MB_CONST _MB_INT *)		incx,
			(_MB_TYPE_D *)				dy,
			(_MB_CONST _MB_INT *)		incy
		) ;
	}

	_MB_TYPE_D meal_dot (const _MB_INT *n, const _MB_TYPE_D *dx, const _MB_INT *incx, const _MB_TYPE_D *dy, const _MB_INT *incy)
	{
		return F77_CALL (ddot)
		(
			(_MB_CONST _MB_INT *)		n,
			(_MB_CONST _MB_TYPE_D *)	dx,
			(_MB_CONST _MB_INT *)		incx,
			(_MB_CONST _MB_TYPE_D *)	dy,
			(_MB_CONST _MB_INT *)		incy
		) ; 
	}

	_MB_TYPE_D meal_nrm2 (const _MB_INT *n, const _MB_TYPE_D *dx, const _MB_INT *incx)
	{
		return F77_CALL (dnrm2)
		(
			(_MB_CONST _MB_INT *)		n,
			(_MB_CONST _MB_TYPE_D *)	dx,
			(_MB_CONST _MB_INT *)		incx
		) ;
	}

	void meal_rot (const _MB_INT *n, _MB_TYPE_D *dx, const _MB_INT *incx, _MB_TYPE_D *dy, const _MB_INT *incy, const _MB_TYPE_D *c, const _MB_TYPE_D *s)
	{
		F77_CALL (drot)
		(
			(_MB_CONST _MB_INT *)		n,
			(_MB_TYPE_D *)				dx,
			(_MB_CONST _MB_INT *)		incx,
			(_MB_TYPE_D *)				dy,
			(_MB_CONST _MB_INT *)		incy,
			(_MB_CONST _MB_TYPE_D *)	c,
			(_MB_CONST _MB_TYPE_D *)	s
		) ;
	}

	void meal_rotg (const _MB_TYPE_D *a, const _MB_TYPE_D *b, _MB_TYPE_D *c, _MB_TYPE_D *s)
	{
		F77_CALL (drotg)
		(
			(_MB_CONST _MB_TYPE_D *)	a,
			(_MB_CONST _MB_TYPE_D *)	b,
			(_MB_TYPE_D *)				c,
			(_MB_TYPE_D *)				s
		) ;
	}

	void meal_rotm (const _MB_INT *n, _MB_TYPE_D *dx, const _MB_INT *incx, _MB_TYPE_D *dy, const _MB_INT *incy, const _MB_TYPE_D *dparam)
	{
		F77_CALL (drotm)
		(
			(_MB_CONST _MB_INT *)		n,
			(_MB_TYPE_D *)				dx,
			(_MB_CONST _MB_INT *)		incx,
			(_MB_TYPE_D *)				dy,
			(_MB_CONST _MB_INT *)		incy,
			(_MB_CONST _MB_TYPE_D *)	dparam
		) ;
	}

	void meal_rotmg (const _MB_TYPE_D *dd1, const _MB_TYPE_D *dd2, const _MB_TYPE_D *dx1, const _MB_TYPE_D *dy1, _MB_TYPE_D *param)
	{
		F77_CALL (drotmg)
		(
			(_MB_CONST _MB_TYPE_D *)	dd1,
			(_MB_CONST _MB_TYPE_D *)	dd2,
			(_MB_CONST _MB_TYPE_D *)	dx1,
			(_MB_CONST _MB_TYPE_D *)	dy1,
			(_MB_TYPE_D *)				param
		) ;
	}

	void meal_scal (const _MB_INT *n, const _MB_TYPE_D *alpha, _MB_TYPE_D *dx, const _MB_INT *incx)
	{
		F77_CALL(dscal)
		(
			(_MB_CONST _MB_INT *)	n,
			(_MB_CONST _MB_TYPE_D *)	alpha,
			(_MB_TYPE_D *)			dx,
			(_MB_CONST _MB_INT *)	incx
		) ; 
	}

	void meal_swap (const _MB_INT *n, _MB_TYPE_D *dx, const _MB_INT *incx, _MB_TYPE_D *dy, const _MB_INT *incy)
	{
		F77_CALL (dswap)
		(
			(_MB_CONST _MB_INT *)		n,
			(_MB_TYPE_D *)				dx,
			(_MB_CONST _MB_INT *)		incx,
			(_MB_TYPE_D *)				dy,
			(_MB_CONST _MB_INT *)		incy
		) ;
	}

	_MB_INT meal_iamax (const _MB_INT *n, const _MB_TYPE_D *dx, const _MB_INT *incx)
	{
		return F77_CALL (idamax)
		(
			(_MB_CONST _MB_INT *)		n,
			(_MB_CONST _MB_TYPE_D *)	dx,
			(_MB_CONST _MB_INT *)		incx
		) ;
	}


//	Level 2 BLAS

	void meal_symv (const _MB_CHAR *uplo, const _MB_INT *n, const _MB_TYPE_D *alpha, const _MB_TYPE_D *a, const _MB_INT *lda, const _MB_TYPE_D *x, const _MB_INT *incx, const _MB_TYPE_D *beta, _MB_TYPE_D *y, const _MB_INT *incy)
	{
		F77_CALL (dsymv)
		(
			(_MB_CONST _MB_CHAR *)		uplo,
			(_MB_CONST _MB_INT *)		n,
			(_MB_CONST _MB_TYPE_D *)	alpha,
			(_MB_CONST _MB_TYPE_D *)	a,
			(_MB_CONST _MB_INT *)		lda,
			(_MB_CONST _MB_TYPE_D *)	x,
			(_MB_CONST _MB_INT *)		incx,
			(_MB_CONST _MB_TYPE_D *)	beta,
			(_MB_TYPE_D *)				y,
			(_MB_CONST _MB_INT *)		incy
            FCONE
		);
	}

	void meal_tbmv (const _MB_CHAR *uplo, const _MB_CHAR *trans, const _MB_CHAR *diag, const _MB_INT *n, const _MB_INT *k, const _MB_TYPE_D *a, const _MB_INT *lda, _MB_TYPE_D *x, const _MB_INT *incx)
	{
		F77_CALL (dtbmv)
		(
			(_MB_CONST _MB_CHAR *)		uplo,
			(_MB_CONST _MB_CHAR *)		trans,
			(_MB_CONST _MB_CHAR *)		diag,
			(_MB_CONST _MB_INT *)		n,
			(_MB_CONST _MB_INT *)		k,
			(_MB_CONST _MB_TYPE_D *)	a,
			(_MB_CONST _MB_INT *)		lda,
			(_MB_TYPE_D *)				x,
			(_MB_CONST _MB_INT *)		incx
            FCONE
            FCONE
            FCONE
		);
	}

	void meal_tpmv (const _MB_CHAR *uplo, const _MB_CHAR *trans, const _MB_CHAR *diag, const _MB_INT *n, const _MB_TYPE_D *ap, _MB_TYPE_D *x, const _MB_INT *incx)
	{
		F77_CALL (dtpmv)
		(
			(_MB_CONST _MB_CHAR *)		uplo,
			(_MB_CONST _MB_CHAR *)		trans,
			(_MB_CONST _MB_CHAR *)		diag,
			(_MB_CONST _MB_INT *)		n,
			(_MB_CONST _MB_TYPE_D *)	ap,
			(_MB_TYPE_D *)				x,
			(_MB_CONST _MB_INT *)		incx
            FCONE
            FCONE
            FCONE
		); 
	}

	void meal_trmv (const _MB_CHAR *uplo, const _MB_CHAR *trans, const _MB_CHAR *diag, const _MB_INT *n, const _MB_TYPE_D *a, const _MB_INT *lda, _MB_TYPE_D *x, const _MB_INT *incx)
	{
		F77_CALL (dtrmv)
		(
			(_MB_CONST _MB_CHAR *)		uplo,
			(_MB_CONST _MB_CHAR *)		trans,
			(_MB_CONST _MB_CHAR *)		diag,
			(_MB_CONST _MB_INT *)		n,
			(_MB_CONST _MB_TYPE_D *)	a,
			(_MB_CONST _MB_INT *)		lda,
			(_MB_TYPE_D *)				x,
			(_MB_CONST _MB_INT *)		incx
            FCONE
            FCONE
            FCONE
		);
	}

	void meal_tbsv (const _MB_CHAR *uplo, const _MB_CHAR *trans, const _MB_CHAR *diag, const _MB_INT *n, const _MB_INT *k, const _MB_TYPE_D *a, const _MB_INT *lda, _MB_TYPE_D *x, const _MB_INT *incx)
	{
		F77_CALL (dtbsv)
		(
			(_MB_CONST _MB_CHAR *)		uplo,
			(_MB_CONST _MB_CHAR *)		trans,
			(_MB_CONST _MB_CHAR *)		diag,
			(_MB_CONST _MB_INT *)		n,
			(_MB_CONST _MB_INT *)		k,
			(_MB_CONST _MB_TYPE_D *)	a,
			(_MB_CONST _MB_INT *)		lda,
			(_MB_TYPE_D *)	x,
			(_MB_CONST _MB_INT *)		incx
            FCONE
            FCONE
            FCONE
		);
	}

	void meal_tpsv (const _MB_CHAR *uplo, const _MB_CHAR *trans, const _MB_CHAR *diag, const _MB_INT *n, const _MB_TYPE_D *ap, _MB_TYPE_D *x, const _MB_INT *incx)
	{
		F77_CALL (dtpsv)
		(
			(_MB_CONST _MB_CHAR *)		uplo,
			(_MB_CONST _MB_CHAR *)		trans,
			(_MB_CONST _MB_CHAR *)		diag,
			(_MB_CONST _MB_INT *)		n,
			(_MB_CONST _MB_TYPE_D *	)	ap,
			(_MB_TYPE_D *)				x,
			(_MB_CONST _MB_INT *)		incx
            FCONE
            FCONE
            FCONE
		);
	}

	void meal_trsv (const _MB_CHAR *uplo, const _MB_CHAR *trans, const _MB_CHAR *diag, const _MB_INT *n, const _MB_TYPE_D *a, const _MB_INT *lda, _MB_TYPE_D *x, const _MB_INT *incx)
	{
		F77_CALL (dtrsv)
		(
			(_MB_CONST _MB_CHAR *)		uplo,
			(_MB_CONST _MB_CHAR *)		trans,
			(_MB_CONST _MB_CHAR *)		diag,
			(_MB_CONST _MB_INT *)		n,
			(_MB_CONST _MB_TYPE_D *)	a,
			(_MB_CONST _MB_INT *)		lda,
			(_MB_TYPE_D *)				x,
			(_MB_CONST _MB_INT *)		incx
            FCONE
            FCONE
            FCONE
		);
	}

#if defined(R_VERSION) && R_VERSION >= R_Version(2, 12, 0)

	void meal_ger (const _MB_INT *m, const _MB_INT *n, const _MB_TYPE_D *alpha, const _MB_TYPE_D *x, const _MB_INT *incx, const _MB_TYPE_D *y, const _MB_INT *incy, _MB_TYPE_D *a, const _MB_INT *lda)
	{
		F77_CALL (dger)
		(
			(_MB_CONST _MB_INT *)		m,
			(_MB_CONST _MB_INT *)		n,
			(_MB_CONST _MB_TYPE_D *)	alpha,
			(_MB_CONST _MB_TYPE_D *)	x,
			(_MB_CONST _MB_INT *)		incx,
			(_MB_CONST _MB_TYPE_D *)	y,
			(_MB_CONST _MB_INT *)		incy,
			(_MB_TYPE_D *)				a,
			(_MB_CONST _MB_INT *)		lda
		);
	}

#else

	void meal_ger (const _MB_INT *m, const _MB_INT *n, const _MB_TYPE_D *alpha, _MB_TYPE_D *x, const _MB_INT *incx, _MB_TYPE_D *y, const _MB_INT *incy, _MB_TYPE_D *a, const _MB_INT *lda)
	{
		F77_CALL (dger)
		(
			(_MB_CONST _MB_INT *)		m,
			(_MB_CONST _MB_INT *)		n,
			(_MB_CONST _MB_TYPE_D *)	alpha,
			(_MB_TYPE_D *)	x,
			(_MB_CONST _MB_INT *)		incx,
			(_MB_TYPE_D *)	y,
			(_MB_CONST _MB_INT *)		incy,
			(_MB_TYPE_D *)				a,
			(_MB_CONST _MB_INT *)		lda
		);
	}


#endif

	void meal_syr (const _MB_CHAR *uplo, const _MB_INT *n, const _MB_TYPE_D *alpha, const _MB_TYPE_D *x, const _MB_INT *incx, _MB_TYPE_D *a, const _MB_INT *lda)
	{
		F77_CALL (dsyr)
		(
			(_MB_CONST _MB_CHAR *)		uplo,
			(_MB_CONST _MB_INT *)		n,
			(_MB_CONST _MB_TYPE_D *)	alpha,
			(_MB_CONST _MB_TYPE_D *)	x,
			(_MB_CONST _MB_INT *)		incx,
			(_MB_TYPE_D *)				a,
			(_MB_CONST _MB_INT *)		lda
            FCONE
		);
	}

	void meal_spr (const _MB_CHAR *uplo, const _MB_INT *n, const _MB_TYPE_D *alpha, const _MB_TYPE_D *x, const _MB_INT *incx, _MB_TYPE_D *ap)
	{
		F77_CALL (dspr)
		(
			(_MB_CONST _MB_CHAR *)		uplo,
			(_MB_CONST _MB_INT *)		n,
			(_MB_CONST _MB_TYPE_D *)	alpha,
			(_MB_CONST _MB_TYPE_D *)	x,
			(_MB_CONST _MB_INT *)		incx,
			(_MB_TYPE_D *)				ap
            FCONE
		);
	}

	void meal_syr2 (const _MB_CHAR *uplo, const _MB_INT *n, const _MB_TYPE_D *alpha, const _MB_TYPE_D *x, const _MB_INT *incx, const _MB_TYPE_D *y, const _MB_INT *incy, _MB_TYPE_D *a, const _MB_INT *lda)
	{
		F77_CALL (dsyr2)
		(
			(_MB_CONST _MB_CHAR *)		uplo,
			(_MB_CONST _MB_INT *)		n,
			(_MB_CONST _MB_TYPE_D *)	alpha,
			(_MB_CONST _MB_TYPE_D *)	x,
			(_MB_CONST _MB_INT *)		incx,
			(_MB_CONST _MB_TYPE_D *)	y,
			(_MB_CONST _MB_INT *)		incy,
			(_MB_TYPE_D *)				a,
			(_MB_CONST _MB_INT *)		lda
            FCONE
		);
	}

	void meal_spr2 (const _MB_CHAR *uplo, const _MB_INT *n, const _MB_TYPE_D *alpha, const _MB_TYPE_D *x, const _MB_INT *incx, const _MB_TYPE_D *y, const _MB_INT *incy, _MB_TYPE_D *ap)
	{
		F77_CALL (dspr2)
		(
			(_MB_CONST _MB_CHAR *)		uplo,
			(_MB_CONST _MB_INT *)		n,
			(_MB_CONST _MB_TYPE_D *)	alpha,
			(_MB_CONST _MB_TYPE_D *)	x,
			(_MB_CONST _MB_INT *)		incx,
			(_MB_CONST _MB_TYPE_D *)	y,
			(_MB_CONST _MB_INT *)		incy,
			(_MB_TYPE_D *)				ap
            FCONE
		);
	}



//	Level 3 BLAS

	void meal_gemm (const _MB_CHAR *transa, const _MB_CHAR *transb, const _MB_INT *m, const _MB_INT *n, const _MB_INT *k, const _MB_TYPE_D *alpha, const _MB_TYPE_D *a, const _MB_INT *lda, const _MB_TYPE_D *b, const _MB_INT *ldb, const _MB_TYPE_D *beta, _MB_TYPE_D *c, const _MB_INT *ldc)
	{
		F77_CALL (dgemm)
		(
			(_MB_CONST _MB_CHAR *)		transa,
			(_MB_CONST _MB_CHAR *)		transb,
			(_MB_CONST _MB_INT *)		m,
			(_MB_CONST _MB_INT *)		n,
			(_MB_CONST _MB_INT *)		k,
			(_MB_CONST _MB_TYPE_D *)	alpha,
			(_MB_CONST _MB_TYPE_D *)	a,
			(_MB_CONST _MB_INT *)		lda,
			(_MB_CONST _MB_TYPE_D *)	b,
			(_MB_CONST _MB_INT *)		ldb,
			(_MB_CONST _MB_TYPE_D *)	beta,
			(_MB_TYPE_D *)				c,
			(_MB_CONST _MB_INT *)		ldc
            FCONE
            FCONE
		) ;
	}

	void meal_trsm (const _MB_CHAR *side, const _MB_CHAR *uplo, const _MB_CHAR *transa, const _MB_CHAR *diag, const _MB_INT *m, const _MB_INT *n, const _MB_TYPE_D *alpha, const _MB_TYPE_D *a, const _MB_INT *lda, _MB_TYPE_D *b, const _MB_INT *ldb)
	{
		F77_CALL (dtrsm)
		(
			(_MB_CONST _MB_CHAR *)		side,
			(_MB_CONST _MB_CHAR *)		uplo,
			(_MB_CONST _MB_CHAR *)		transa,
			(_MB_CONST _MB_CHAR *)		diag,
			(_MB_CONST _MB_INT *)		m,
			(_MB_CONST _MB_INT *)		n,
			(_MB_CONST _MB_TYPE_D *)	alpha,
			(_MB_CONST _MB_TYPE_D *)	a,
			(_MB_CONST _MB_INT *)		lda,
			(_MB_TYPE_D *)				b,
			(_MB_CONST _MB_INT *)		ldb
            FCONE
            FCONE
            FCONE
            FCONE
		);
	}

	void meal_trmm (const _MB_CHAR *side, const _MB_CHAR *uplo, const _MB_CHAR *transa, const _MB_CHAR *diag, const _MB_INT *m, const _MB_INT *n, const _MB_TYPE_D *alpha, const _MB_TYPE_D *a, const _MB_INT *lda, _MB_TYPE_D *b, const _MB_INT *ldb)
	{
		F77_CALL (dtrmm)
		(
			(_MB_CONST _MB_CHAR *)		side,
			(_MB_CONST _MB_CHAR *)		uplo,
			(_MB_CONST _MB_CHAR *)		transa,
			(_MB_CONST _MB_CHAR *)		diag,
			(_MB_CONST _MB_INT *)		m,
			(_MB_CONST _MB_INT *)		n,
			(_MB_CONST _MB_TYPE_D *)	alpha,
			(_MB_CONST _MB_TYPE_D *)	a,
			(_MB_CONST _MB_INT *)		lda,
			(_MB_TYPE_D *)				b,
			(_MB_CONST _MB_INT *)		ldb
            FCONE
            FCONE
            FCONE
            FCONE
		);
	}

	void meal_symm (const _MB_CHAR *side, const _MB_CHAR *uplo, const _MB_INT *m, const _MB_INT *n, const _MB_TYPE_D *alpha, const _MB_TYPE_D *a, const _MB_INT *lda, const _MB_TYPE_D *b, const _MB_INT *ldb, const _MB_TYPE_D *beta, _MB_TYPE_D *c, const _MB_INT *ldc)
	{
		F77_CALL (dsymm)
		(
			(_MB_CONST _MB_CHAR *)		side,
			(_MB_CONST _MB_CHAR *)		uplo,
			(_MB_CONST _MB_INT *)		m,
			(_MB_CONST _MB_INT *)		n,
			(_MB_CONST _MB_TYPE_D *)	alpha,
			(_MB_CONST _MB_TYPE_D *)	a,
			(_MB_CONST _MB_INT *)		lda,
			(_MB_CONST _MB_TYPE_D *)	b,
			(_MB_CONST _MB_INT *)		ldb,
			(_MB_CONST _MB_TYPE_D *)	beta,
			(_MB_TYPE_D *)				c,
			(_MB_CONST _MB_INT *)		ldc
            FCONE
            FCONE
		);
	}

	void meal_syrk (const _MB_CHAR *uplo, const _MB_CHAR *trans, const _MB_INT *n, const _MB_INT *k, const _MB_TYPE_D *alpha, const _MB_TYPE_D *a, const _MB_INT *lda, const _MB_TYPE_D *beta, _MB_TYPE_D *c, const _MB_INT *ldc)
	{
		F77_CALL (dsyrk)
		(
			(_MB_CONST _MB_CHAR *)		uplo,
			(_MB_CONST _MB_CHAR *)		trans,
			(_MB_CONST _MB_INT *)		n,
			(_MB_CONST _MB_INT *)		k,
			(_MB_CONST _MB_TYPE_D *)	alpha,
			(_MB_CONST _MB_TYPE_D *)	a,
			(_MB_CONST _MB_INT *)		lda,
			(_MB_CONST _MB_TYPE_D *)	beta,
			(_MB_TYPE_D *)				c,
			(_MB_CONST _MB_INT *)		ldc
            FCONE
            FCONE
		);
	}

	void meal_syr2k (const _MB_CHAR *uplo, const _MB_CHAR *trans, const _MB_INT *n, const _MB_INT *k, const _MB_TYPE_D *alpha, const _MB_TYPE_D *a, const _MB_INT *lda, const _MB_TYPE_D *b, const _MB_INT *ldb, const _MB_TYPE_D *beta, _MB_TYPE_D *c, const _MB_INT *ldc)
	{
		F77_CALL (dsyr2k)
		(
			(_MB_CONST _MB_CHAR *)		uplo,
			(_MB_CONST _MB_CHAR *)		trans,
			(_MB_CONST _MB_INT *)		n,
			(_MB_CONST _MB_INT *)		k,
			(_MB_CONST _MB_TYPE_D *)	alpha,
			(_MB_CONST _MB_TYPE_D *)	a,
			(_MB_CONST _MB_INT *)		lda,
			(_MB_CONST _MB_TYPE_D *)	b,
			(_MB_CONST _MB_INT *)		ldb,
			(_MB_CONST _MB_TYPE_D *)	beta,
			(_MB_TYPE_D *)				c,
			(_MB_CONST _MB_INT *)		ldc
            FCONE
            FCONE
		);
	}

#undef _MB_TYPE_D
#undef _MB_SYM
#undef _MB_INT
#undef _MB_CHAR

#endif	//	#ifdef R_PACKAGE_FILE
