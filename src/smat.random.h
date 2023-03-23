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

#ifndef SMAT_RANDOM_H
#define SMAT_RANDOM_H

	double runif () ;
	double rnorm () ;
	double rexp () ;

	void runif_raw (double *p, double *pEnd) ;
	void runif (const SVData<double> &a) ;

	void runif_r_raw (double *p, double *pEnd) ;
	void runif_r (const SVData<double> &a) ;

	
	void rnorm (const SVData<double> &a) ;
	void rnorm_raw (double *p, double *pEnd) ;

	void runif_raw (double *d, int l, double dL, double dU) ;
	void runif_raw (double *d, double * const pEnd, double dL, double dU) ;

	void SampleNoReplace(int k, int n, int *y, int *x) ;


#endif	//	#ifndef SMAT_RANDOM_H
