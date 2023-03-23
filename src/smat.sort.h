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

#ifndef SMAT_SORT_H
#define SMAT_SORT_H
#include "smat.base.h"

	void sort_R (const SCData<double> &a, SVecD &b) ;
	void sort (const SCData<double> &a, const SVecD &b) ;
	void sort_NC (const SCData<double> &a, const SVecD &b) ;

	void sort (const SVData<double> &a) ;
	void sort_order (const SVData<double> &a, const SVData<int> &b) ;
	void sort_order_NC (const SVData<double> &a, const SVData<int> &b) ;

	template <class T>
	T psort_V (const SVData <T> &a, t_size k)
	{
		T *pA = a ;

		double pivot, swapBuf ;

		t_size l = 0, lr = a.size () - 1, jnc, j ;

		while (l < lr)
		{
			pivot = pA[k] ;

			jnc = l ;
			j = lr ;

			while (jnc <= j)
			{
				while (pA[jnc] < pivot)
					++jnc ;

				while (pA[j] > pivot)
					--j ;

				if(jnc <= j)
				{
					sm_swap (pA [jnc], pA[j], swapBuf) ;
					++jnc ;
					--j ;
				}
			}

			if (j < k)
				l = jnc ;
			if (k < jnc)
				lr = j ;
		}

		return pA [k] ;
	}

#endif	//	#ifndef SMAT_SORT_H
