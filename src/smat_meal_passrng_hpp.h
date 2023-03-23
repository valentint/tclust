#ifndef SMAT_MEAL_PASSRNG_HPP
#define SMAT_MEAL_PASSRNG_HPP

#ifndef SMAT_MEAL_H
	#error "File 'smat.meal.h' not included. Please include before including 'smat_meal_passrng.hpp'."
#endif

#ifndef SMAT_MEAL_PASSRNG_H
	#error "File 'smat_meal_passrng.h' not included. Please include before including 'smat_meal_passrng.hpp'."
#endif

//	when including this files make sure, that "smat_meal_passrng.h" and "smat.meal.h" has been included first!

	CPassRng::CPassRng ()
		: m_pCurData (NULL)
		, m_pDataEnd (NULL)
	{

	}

	double CPassRng::Get ()
	{
		THROW (m_pCurData < m_pDataEnd) ;
		return *(m_pCurData++) ;
	}

	void CPassRng::Set (double *pData, t_size n)
	{
		m_data.Require (n) ;
		m_data.Copy (pData, n) ;

		m_pCurData = m_data ;
		m_pDataEnd = m_data.GetDataEnd () ;
	}

	CPassRng &GetPassRng_runif ()
	{
		static CPassRng obj ;
		return obj ;
	}

	CPassRng &GetPassRng_rnorm ()
	{
		static CPassRng obj ;
		return obj ;
	}

	CPassRng &GetPassRng_rexp ()
	{
		static CPassRng obj ;
		return obj ;
	}

	void pass_runif (double *pd, int n) { GetPassRng_runif ().Set (pd, n) ; }
	void pass_rnorm (double *pd, int n) { GetPassRng_rnorm ().Set (pd, n) ;	}
	void pass_rexp  (double *pd, int n) { GetPassRng_rexp  ().Set (pd, n) ;	}

	void meal_PutRNGstate () {}
	void meal_GetRNGstate () {}

	double meal_unif_rand () { return GetPassRng_runif ().Get () ; }
	double meal_norm_rand () { return GetPassRng_rnorm ().Get () ; }
	double meal_exp_rand  () { return GetPassRng_rexp  ().Get () ; }

#endif	//	#ifndef SMAT_MEAL_PASSRNG_HPP
