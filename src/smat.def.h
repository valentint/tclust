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

#ifndef ES_SMAT_DEF_H
#define ES_SMAT_DEF_H

typedef unsigned int UINT ;
typedef UINT BOOL ;

#define FALSE 0
#define TRUE 1

#ifndef NULL
	#define NULL 0
#endif

typedef  UINT t_size ;	//	2do: use size_t instead
typedef  UINT t_count ;

#define NAI	((t_size) -1)	//	not an index

#ifdef _DEBUG
	#define ASSERT(a)	THROW(a)
#else
	#define ASSERT(a)
#endif

class SMat_EXC
{
	public:
		SMat_EXC (char const * const szDate, char const * const szFile, const int nLine)
			: m_szDate (szDate), m_szFile (szFile), m_nLine (nLine) {}
		const char * const GetDate () const { return m_szDate ; }
		const char * const GetFile () const { return m_szFile ; }
		const int GetLine () const { return m_nLine ; }

		void OnException () ;
		static void OnUException () ;

	private:
		const char * const m_szDate, * const m_szFile ;
		const int m_nLine ;
} ;

#ifdef _DEBUG
	#define TRY(C)	C
#else
	#define TRY(C)	try {C} catch (SMat_EXC exc)  { exc.OnException () ; } catch (...) { SMat_EXC::OnUException () ; }
#endif

#define THROW(a) {if (!(a)){THROW_BASE}}
#define THROW_BASE { { throw (SMat_EXC (__DATE__, __FILE__, __LINE__)) ; } }

template <class T> inline const T sm_sqr (const T &a) { return a * a ; }

template <class T> void sm_swap (T &a, T &b)	//	-> move to COP?
{
	T c (a) ;
	a = b ;
	b = c ;
}

template <class T> void sm_swap (T &a, T &b, T&temp)
{
	temp = a;
	a = b;
	b = temp ;
}

///////////////
//	min max  //	//	-> move to COP?
///////////////

template <class T> const T &sm_min (const T &a, const T &b)
{
	if (a < b)
		return a ;
	return b ;
}

template <class T> const T &sm_max (const T &a, const T &b)
{
	if (a > b)
		return a ;
	return b ;
}

template <class TA, class TB>
TA &sm_setmax_t (TA &a, const TB &b)
{
	if (a < b)
		a = b ;
	return a ;
}

template <class TA, class TB>
TA &sm_setmin_t (TA &a, const TB &b)
{
	if (a > b)
		a = b ;
	return a ;
}		

template <class TA, class TB>
BOOL sm_setmax_b (TA &a, const TB &b)
{
	if (a >= b)
		return FALSE ;
	a = b ;
	return TRUE ;
}

template <class TA, class TB>
BOOL sm_setmin_b (TA &a, const TB &b)
{
	if (a <= b)
		return FALSE; 
	a = b ;
	return TRUE ;
}

template <class TA, class TB>
void sm_setmax (TA &a, const TB &b)
{
	if (a < b)
		a = b ;
}

template <class TA, class TB>
void sm_setmin (TA &a, const TB &b)
{
	if (a > b)
		a = b ;
}


#endif	//#ifndef ES_SMAT_DEF_H
