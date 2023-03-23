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

//	ML_meal.h
//	MATLAB Mathematical Environment Abstraction Layer

#ifndef ML_MEAL_H
#define ML_MEAL_H
#ifdef MATLAB_MEX_FILE

#include "matrix.h"
#include "mex.h"

#include "lapack.h"
#include "blas.h"

//#include "mwutil.h"

#ifdef SMAT_H

	void Int2Double (double *pd, const int *pn, int n)
	{
		EO<SOP::assign>::VVc_raw (pd, pd + n, pn) ;
	}

	void Double2Int (int *pn, const double *pd, int n)
	{
		EO<SOP::assign>::VVc_raw (pn, pn + n, pd) ;
	}

#endif	//	#ifdef SMAT_H

	class CRmealSettings 
	{
	public:
		CRmealSettings () ;
		CRmealSettings (const char *szEmail) ;
		const char *GetEmail () { return m_szEmail ; }
	protected:
		const char *m_szEmail ;
	} ;

#endif	//	#ifdef MATLAB_MEX_FILE

#endif	//	#ifndef ML_MEAL_H
