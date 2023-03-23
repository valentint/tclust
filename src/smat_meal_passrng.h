#ifndef SMAT_MEAL_PASSRNG_H
#define SMAT_MEAL_PASSRNG_H

#ifdef ES_DEV_ENV
	#include "../SMat/smat.h"
#else
	#include "smat.h"
#endif

	class CPassRng
	{
	public:

		CPassRng () ;

		double Get () ;
		void Set (double *pData, t_size n) ;

	protected:
		SVecD m_data ;
		double *m_pCurData, *m_pDataEnd ;
	} ;

	void pass_runif (double *pd, int n) ;
	void pass_rnorm (double *pd, int n) ;
	void pass_rexp  (double *pd, int n) ;


#endif	//	#ifndef SMAT_MEAL_PASSRNG_H
