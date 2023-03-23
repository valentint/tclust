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

#ifndef SMAT_MEM_H
#define SMAT_MEM_H

	template <class TA, class TB>
	void CopyRow_R (SVData<TA> &a, const SCMat<TB> &b, const t_size &c)	//	2do: rename to CopyRow_R
	{
		THROW (c < b.nrow  ()) ;
		a.Require (b.ncol()) ;

		CopyRow_NC (a, b, c) ;
	}

	template <class TA, class TB>
	void CopyRow (const SVData<TA> &a, const SCMat<TB> &b, const t_size &c)
	{
		THROW (c < b.nrow  ()) ;
		THROW (a.size () == b.ncol()) ;
		CopyRow_NC (a, b, c) ;
	}

	template <class TA, class TB>
	void CopyRow_NC (const SVData<TA> &a, const SCMat<TB> &b, const t_size &c)
	{
		ASSERT (c < b.nrow  ()) ;
		ASSERT (a.size () == b.ncol()) ;

		TA *pA = a ;//, * const pEndA = a.GetDataEnd () ;
		TB const *pB = b.GetData (c, 0), * const pEndB = b.GetDataEnd () ;

		while (pB < pEndB)
		{
			*pA = *pB ;
			++pA ;
			pB += b.GetColInc () ;
		}
	}

	template <class TA, class TB>
	void CopyRow (const SVMat<TA> &a, const SCData<TB> &b, const t_size &c)
	{
		THROW (c < a.nrow  ()) ;
		THROW (a.ncol () == b.size ()) ;
		CopyRow_NC (a, b, c) ;
	}

	template <class TA, class TB>
	void CopyRow_NC (const SVMat<TA> &a, const SCData<TB> &b, const t_size &c)
	{
		ASSERT (c < a.nrow  ()) ;
		ASSERT (a.ncol () == b.size ()) ;
		
		TA *pA = a.GetData (c, 0), * const pEndA = a.GetDataEnd () ;
		TB const * pB = b ;

		while (pA < pEndA)
		{
			*pA = *pB ;
			pA += a.GetColInc () ;
			++pB ;
		}
	}

	template <class TA, class TB>
	void ResetRow (const SVMat<TA> &a, const TB &b, const t_size &c)
	{
		THROW (c < a.nrow  ()) ;
		ResetRow_NC (a, b, c) ;
	}

	template <class TA, class TB>
	void ResetRow_NC (const SVMat<TA> &a, const TB &b, const t_size &c)
	{
		ASSERT (c < a.nrow  ()) ;

		TA *pA = a.GetData (c, 0), * const pEndA = a.GetDataEnd () ;

		while (pA < pEndA)
		{
			*pA = b ;
			pA += a.GetColInc () ;
		}
	}



	template <class TA>
	TA *Reset (TA *pA, const TA *pEndA, const TA &val)
	{
		for (; pA < pEndA; pA++)
			*pA = val ;
		return pA ;
	}

	template <class TA>
	TA *Reset (TA *pA, const TA *pEndA)
	{
		for (; pA < pEndA; pA++)
			*pA = 0 ;
		return pA ;
	}

	template <class TA, class TB>
	void Copy (TA * pA, TB const * pB, t_size n)	//	the "master" - copy function
	{
		//memcpy (pA, pB, sizeof (T) * n) ;	//	2do: implement with a for - loop using "=". -> Copy constructor would be used!

		TA * const pEndA = pA + n ;
		while (pA < pEndA)
		{
			*pA = *pB ;
			++pA ;
			++pB ;
		}
	}

	template <class TA, class TB>
	void Copy (const SVData <TA> &a, const TB * const pB)
	{
		Copy (a.GetData (), pB, a.size ()) ;
	}

	template <class TA, class TB>
	void Copy (TA * const pA, const SCData <TB> &b)
	{
		Copy (pA, b.GetData (), b.size ()) ;
	}

	template <class TA, class TB>
	void Copy_R (SVVec <TA> &a, const SCData<TB> &b)
	{
		a.Require (b) ;
		Copy_NC (a, b) ;
	}

	template <class TA, class TB>
	void Copy (const SVData <TA> &a, const SCData<TB> &b)
	{
		THROW (a.size () == b.size ()) ;
		Copy_NC (a, b) ;
	}

	template <class TA, class TB>
	void Copy_NC (const SVData <TA> &a, const SCData<TB> &b)
	{
		ASSERT (a.size () == b.size ()) ;
		Copy (a.GetData (), b.GetData (), b.size ()) ; 
	}
	
	template <class TA, class TB>
	void CopyDiag_R (SVec <TA> &v, const SCMat <TB> &m)
	{
		v.Require (m.GetMinDim ()) ;
		CopyDiag (v, m) ;
	}

	template <class A, class B>
	void CopyDiag (const SVec <A> &v, const SCMat <B> &m)
	{
		THROW (v.nsize () == m.GetMinDim ()) ;
		CopyDiag_NC (v, m) ;
	}

	template <class A, class B>
	void CopyDiag_NC (const SVec <A> &v, const SCMat <B> &m)
	{
		ASSERT (v.size () == m.GetMinDim ()) ;
		t_size i ;
		for (i = v.size () - 1; i != NAI; i--)
			v (i) = m (i, i) ;
	}

	template <class TA, class TB>
	void CopyCol (const SVData<TA> &a, const SCMat<TB> &b, t_size dwCol)
	{
		THROW (dwCol < b.ncol ()) ;
		CopyCol_NC (a, b, dwCol) ;
	}

	template <class TA, class TB>
	void CopyCol_NC (const SVData<TA> &a, const SCMat<TB> &b, t_size dwCol)
	{
		ASSERT (dwCol < b.ncol ()) ;
		Copy (a.GetData (), b.GetData_NC (0, dwCol), b.nrow ()) ;
	}

	template <class TA, class TB>
	void CopyCol (const SVMat<TA> &a, const SCMat<TB> &b, t_size dwStart, t_size dwEnd)
	{
		THROW (dwStart <= dwEnd) ;
		THROW (dwEnd <= b.ncol ()) ;
		THROW (a.nrow () == b.nrow ()) ;
		THROW (a.ncol () == dwEnd - dwStart) ;
		CopyCol_NC (a, b, dwStart, dwEnd) ;
	}

	template <class TA, class TB>
	void CopyCol_NC (const SVMat<TA> &a, const SCMat<TB> &b, t_size dwStart, t_size dwEnd)
	{
		ASSERT (dwStart <= dwEnd) ;
		ASSERT (dwEnd <= b.ncol ()) ;
		ASSERT (a.nrow () == b.nrow ()) ;
		ASSERT (a.ncol () == dwEnd - dwStart) ;
		Copy (a.GetData (), b.GetData_NC (0, dwStart), b.nrow () * (dwEnd - dwStart)) ;
	}



#endif	//	#ifndef SMAT_MEM_H
