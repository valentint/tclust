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


	void sme_eigen_sqr_R	(const SCMatD &A,		SVecD &vVal,		SVMatD &mVec, BOOL bOrder) ;	//	temp4
	void sme_eigen_sqr		(const SCMatD &A, const SVecD &vVal, const	SVMatD &mVec, BOOL bOrder) ;	//	temp4
	void sme_eigen_sqr_NC	(const SCMatD &A, const SVecD &vVal, const	SVMatD &mVec, BOOL bOrder) ;	//	temp4

	void sme_eigen_sqr_RV	(const SVMatD &A,		SVecD &vVal,		SVMatD &mVec, BOOL bOrder) ;	//	temp3
	void sme_eigen_sqr_V		(const SVMatD &A, const	SVecD &vVal, const	SVMatD &mVec, BOOL bOrder) ;	//	temp3
	void sme_eigen_sqr_NCV	(const SVMatD &A, const	SVecD &vVal, const	SVMatD &mVec, BOOL bOrder) ;	//	temp3

	void sme_eigen_sym_2x2_norm_raw		(double * const pdEval, double *const pdEVec, const double *const pd, const double &dZeroTol) ;
	void sme_eigen_sym_2x2_norm_raw_NC	(double * const pdEval, double *const pdEVec, const double *const pd, const double &dZeroTol) ;

//2do: implement solver
