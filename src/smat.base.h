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

#ifndef ES_SMAT_BASE_H
#define ES_SMAT_BASE_H
#include "smat.def.h"
#include <memory.h>

#ifdef _DEBUG
	#define ASSERT_TEMPRANGE(L, U) SDataRefCont::CRefRange temp_assert_temprange(L, U, &GetPermTempRefRange ())
#else
	#define ASSERT_TEMPRANGE(L, U)
#endif


//////////////////////////
//	Forward References  //
//////////////////////////////

	template <class T> class SMat ;
	template <class T> class SCMat ;
	template <class T> class SVMat ;
	template <class T> class SVec ;
	template <class T> class SCVec ;
	template <class T> class SVVec ;
	template <class T> class SCData ;
	template <class T> class SVData ;

///////////////////////////////
//	Basic type declarations  //	
///////////////////////////////
#define SMAT_TYPE (TYPE, TOKEN)				\
											\
	typedef SMat	<TYPE> SMat##TOKEN ;	\
	typedef SCMat	<TYPE> SCMat##TOKEN ;	\
	typedef SVMat	<TYPE> SVMat##TOKEN ;	\
	typedef SVec	<TYPE> SVec##TOKEN ;	\
	typedef SCVec	<TYPE> SCVec##TOKEN ;	\
	typedef SVVec	<TYPE> SVVec##TOKEN ;	\
	typedef SCData	<TYPE> SCData##TOKEN ;	\
	typedef SVData	<TYPE> SVData##TOKEN ;

		//2do: use SMAT_TYPE - macro

	//	double
	typedef SMat<double> SMatD ;
	typedef SCMat<double> SCMatD ;
	typedef SVMat<double> SVMatD ;
	typedef SVec<double> SVecD ;
	typedef SCVec<double> SCVecD ;
	typedef SVVec<double> SVVecD ;
	typedef SCData<double> SCDataD ;
	typedef SVData<double> SVDataD ;

	//	int
	typedef SMat<int> SMatN ;
	typedef SCMat<int> SCMatN ;
	typedef SVMat<int> SVMatN ;
	typedef SVec<int> SVecN ;
	typedef SCVec<int> SCVecN ;
	typedef SVVec<int> SVVecN ;
	typedef SCData<int> SCDataN ;
	typedef SVData<int> SVDataN ;

////////////////
//	SDataRef  //
////////////////

	class SDataRef
	{
	protected:
		SDataRef () ;
//		SDataRef (const t_size dwRef) ;
	public:
		SDataRef (const t_size dwCount) ;
//		SDataRef (const t_size dwCount, const t_size dwRef) ;
		SDataRef (const t_size dwCount, void * const pData) ;
		~SDataRef () ;

		SDataRef *Reference (SDataRef *&pRef) ;
		SDataRef *Dereference (SDataRef *&pRef) ;

		static void sDeref (SDataRef *&pRef) ;
		SDataRef *Ref (SDataRef *&pRef) ;
		SDataRef *Ref_NDR (SDataRef *&pRef)	;//	No DeRereference

		inline t_count GetRef () const { return m_dwRef ; }
		inline t_size GetSize () const { return GetSizeRef () ; }
		inline BOOL IsOwner () const { return m_bOwner ; }
		inline void *GetEnd () const { return GetDataEndRef () ; }

		BOOL Require (t_size s, SDataRef *&pRef) ;
		SDataRef *Recreate (t_size dwSize, SDataRef *&pRef) ;
		inline void *GetData () const { return GetDataRef () ; }
		inline BOOL IsStatic () const { return m_bStatic ; }
		
		static SDataRef &Empty () ;

		void FreeIfIdle () ;

		//inline t_size GetItemSize () { return m_dwItemSize ; }
	protected:
		inline void SetStatic () { m_bStatic = TRUE ; }

		inline void IncRef () { ++m_dwRef ; }
		inline void DecRef () { --m_dwRef ; }

		void Alloc (t_size dwSize) ;
		void Alloc_NF (t_size dwSize) ;

//		void CalcDataEnd ();

		void Free () ;
		void CheckFree () ;
		BOOL Deref () ;

		void *m_pData, *m_pDataEnd ;
		t_size m_dwRef, m_dwCount ;
		//const t_size m_dwItemSize ;
		const BOOL m_bOwner ;
		BOOL m_bStatic ;

	private:
		void *& GetDataRef () { return m_pData ; }
		void *& GetDataEndRef () { return m_pDataEnd ; }
		t_size & GetSizeRef () { return m_dwCount ; } 

		void * const& GetDataRef () const { return m_pData ; }
		void * const& GetDataEndRef () const { return m_pDataEnd ; }
		const t_size & GetSizeRef () const { return m_dwCount ; } 
	} ;

	void Deref (SDataRef *&pRef) ;

	class SDataRef_Static: public SDataRef
	{
		typedef SDataRef t_base ;
	public:
		SDataRef_Static (const t_size dwCount = 0, const BOOL bStatic = TRUE) ;
		~SDataRef_Static () ;

		SDataRef_Static &Require (t_size dwSize) ;
	protected:
		
	} ;

	class SDataRefCont
	{
		typedef SDataRef_Static t_item ;
		typedef t_item* t_pitem ;
	public:

		SDataRefCont () ;
		~SDataRefCont () ;

		class CRefRange
		{
		public:
			CRefRange (const t_size dwL = NAI, const t_size dwU = NAI, CRefRange *pPerm = NULL)
				: m_dwL (dwL), m_dwU (dwU), m_pPerm (pPerm)
			{
				ASSERT (dwL <= dwU) ;
				if (pPerm)
				{
					ASSERT (m_pPerm->GetL () > dwU) ;
					pPerm->Get (m_dwLOld, m_dwUOld) ;
					pPerm->Set (dwL, dwU) ;
				}
			}
			
			~CRefRange ()
			{
				if (m_pPerm)
					m_pPerm->Set (m_dwLOld, m_dwUOld) ;
			}

			t_size GetL () const { return m_dwL ; } 
			t_size GetU () const { return m_dwU ; } 

			void Get (t_size &dwL, t_size &dwU)	const { dwL = m_dwL, dwU = m_dwU ; }
			void Set (t_size dwL, t_size dwU)		{ m_dwL = dwL, m_dwU = dwU ; }


		protected:
			t_size m_dwL, m_dwU ;
			t_size m_dwLOld, m_dwUOld ;
			CRefRange *m_pPerm ;
		} ;

		inline t_size GetSize () { return sizeRef () ; }
		inline t_size GetMemSize () { return GetSize () * sizeof (t_pitem) ; }
		inline t_pitem * GetData () { return dataRef (); }

		void Require (t_size dwCount) ;
		t_item &Item (t_size dwIdx) ;

		void FreeIfIdle ();

	protected:
		void Free () ;

		inline t_size &sizeRef () { return m_dwSize ; }
		inline t_pitem* &dataRef () { return m_ppData; }

		t_pitem *m_ppData ;
		t_size m_dwSize ;
	} ;

	SDataRefCont &GetTempCont () ;
	SDataRefCont::CRefRange &GetPermTempRefRange () ;
	void RequireTemp (t_size dwCount) ;
	SDataRef_Static &tempRef (t_size dwIdx) ;
	void FreeTempCont () ;

	template <class T>
	T *tempRef (t_size dwIdx, T * & p, t_size dwSize)
	{
		SDataRef_Static &ref = tempRef (dwIdx) ;
		ref.Require (dwSize * sizeof (T)) ;
		return p = (T *) ref.GetData () ;
	}

	template <class T>
	T *tempArray (t_size dwIdx, t_size dwSize)
	{
		SDataRef_Static &ref = tempRef (dwIdx) ;
		ref.Require (dwSize * sizeof (T)) ;
		return (T *) ref.GetData () ;
	}

////////////////////
//	CDataCont_NT  //
////////////////////

	class CDataCont_NT
	{
	public:
		CDataCont_NT ()
		{
			++GetInstanceCount () ;
		}

		~CDataCont_NT ()
		{
			if (!--GetInstanceCount ())
				FreeTempCont () ;
		}

	private:
		static t_size &GetInstanceCount () ;
	} ;

#define TEMP_GUARD CDataCont_NT __TEMP_GUARD


/////////////////
//	SVData  //
/////////////////

	template <class T>
		class SVData	:public CDataCont_NT//	No Dimension
	{
		typedef SVData <T> t_this  ; 

		typedef SDataRef_Static t_tempRef ;

	public:

//	Destuctor  //

		~SVData ()
		{
			SDataRef::sDeref (m_pDataRef) ;
		}

//	Constructor  //

		SVData ()
		{
			Ref_NDR (SDataRef::Empty ()) ;
			ResetOffsetSize () ;
		}

		SVData (SDataRef &ref)
		{
			Ref_NDR (ref) ;
			t_this::ResetOffset (Bytes2Size (ref.GetSize ())) ;
		}

		SVData (SDataRef_Static &ref)
		{
			Ref_NDR (ref) ;
			t_this::ResetOffset (0) ;
		}

		SVData (SDataRef &ref, t_size dwSize)
		{
			Ref_NDR (ref) ;
			t_this::ResetOffset (dwSize) ;
		}

		SVData (SDataRef_Static &ref, t_size dwSize)
			: m_dwOffset (0)
		{
			Ref_NDR (ref) ;
			t_this::Require (dwSize) ;
		}

		SVData (SDataRef &ref, t_size dwOffset, t_size dwSize)
			: m_dwOffset (dwOffset)
		{
			Ref_NDR (ref) ;
			t_this::SetSize (dwSize) ;
		}

		SVData (SDataRef_Static &ref, t_size dwOffset, t_size dwSize)
		{
			Ref_NDR (ref) ;
			t_this::Require (dwOffset, dwSize) ;
		}

		SVData (const t_this &p)
		{
			Ref_NDR (p.GetDataRef ()) ;
			t_this::ResetOffsetSize () ;
		}

		SVData (const t_this &p, t_size dwSize)
		{
			Ref_NDR (p.GetDataRef ()) ;
			t_this::ResetOffset (dwSize) ;
		}

		SVData (const t_this &p, t_size dwOffset, t_size dwSize)
			: m_dwOffset (dwOffset)
		{
			Ref_NDR (p.GetDataRef ()) ;
			t_this::SetSize (dwSize) ;
		}

		SVData (const t_size dwSize)
		{
			Ref_NDR (*new SDataRef (Size2Bytes (dwSize))) ;
			t_this::ResetOffset_NC (dwSize) ;
		}

		SVData (T * const pData, const t_size dwSize)
		{
			Ref_NDR (*new SDataRef (Size2Bytes (dwSize), pData)) ;
			t_this::ResetOffset_NC (dwSize) ;
		}

		SVData (SDataRef_Static &ref, const t_this &p)		//	this is supposed to copy the data from p to ref.
		{
			Ref_NDR (ref) ;
			if (&ref == &p.GetDataRef ())						//	they act on the same reference. thus no need to copy the data
				t_this::SetOffset_NC (p.GetOffset (), p.GetSize ()) ;
			else
			{
				m_dwOffset = 0 ;
				Require (p.GetSize ()) ;
				memcpy (ref.GetData (), p.GetData (), Size2Bytes (p.GetSize ())) ;
			}
		}	

		t_this & operator = (const t_this &p)
		{
			Ref (p.GetDataRef ()) ;
			GetSizeRef		() = p.GetSize () ;
			GetOffsetRef	() = p.GetOffset () ;
			GetEndRef		() = p.GetEnd () ;
			return *this ;
		}


		inline t_size size () const { return GetSizeRef () ; }

		inline T &operator () (const t_size dwIdx) const { return GetData (dwIdx) ; }
		inline operator T * () const { return GetData () ; }
		inline T *GetDataEnd	() const { return GetRawData () + GetEnd () ; }
		inline T *GetData		() const { return GetRawData () + GetOffset	() ; }


		inline T *GetData (t_size dwIdx) const
		{
			THROW (dwIdx < GetSize ()) ;
			return GetData_NC (dwIdx) ;
		}

		inline T *GetData_NC (t_size dwIdx) const
		{
			ASSERT (dwIdx < GetSize ()) ;
			return GetData () + dwIdx ;
		}

		void Reset (const T& v) const
		{
			T *pCur ;
			T *const pEnd = t_this::GetDataEnd () ;
			for (pCur = t_this::GetData (); pCur < pEnd; pCur++)
				*pCur = v ;
		}


		void Set (T * const pData, const t_size dwSize)
		{
			Ref (*new SDataRef (Size2Bytes (dwSize), pData)) ;
			t_this::ResetOffset_NC (dwSize) ;
		}

	protected:
		void Redim				()												{ t_this::ResetOffsetSize 	()					; }
		void Redim				(const t_size dwSize)							{ t_this::ResetOffset		(dwSize)			; }
		void Redim_NC			(const t_size dwSize)							{ t_this::ResetOffset_NC	(dwSize)			; }

		void Reshape			(t_size dwSize)									{ t_this::SetSize			(dwSize)			; }
		void Reshape			(t_size dwOffset, t_size dwSize)				{ t_this::SetOffset			(dwOffset, dwSize)	; }
		void Reshape_NC			(t_size dwSize)									{ t_this::SetSize_NC		(dwSize)			; }
		void Reshape_NC			(t_size dwOffset, t_size dwSize)				{ t_this::SetOffset_NC		(dwOffset, dwSize)	; }

		void Recreate (t_size dwSize)
		{
			t_this::GetDataRef ()->Recreate (Size2Bytes (dwSize), t_this::GetDataRef ()) ;
			t_this::ResetOffset (dwSize) ;
		}

		void Require (t_size dwSize)
		{
			if (t_this::GetDataRef ().Require (Size2Bytes (dwSize), m_pDataRef))
				t_this::ResetOffset (dwSize) ;
			else
				t_this::SetSize (dwSize) ;
		}

		void Require (t_size dwOffset, t_size dwSize)
		{
			t_this::GetDataRef ().Require (Size2Bytes (dwOffset + dwSize), m_pDataRef) ;
			t_this::SetOffset_NC (dwOffset, dwSize) ;
		}

		inline T &Item (const t_size dwIdx) const
		{
			THROW (dwIdx < GetSize ()) ;
			return GetData () [dwIdx] ;
		}

		inline T &Item_NC (const t_size dwIdx) const
		{
			ASSERT (dwIdx < GetSize ()) ;
			return GetData () [dwIdx] ;
		}



		void SetDataRef (SDataRef &ref) { ref.Ref (m_pDataRef) ; }
		void SetDataRef_NDR	(SDataRef &ref) { ref.Ref_NDR (m_pDataRef) ; }

		inline SDataRef &GetDataRef () const { return *m_pDataRef; }
		inline T *GetRawData () const { return (T *) GetDataRef ().GetData () ; }
		inline t_size GetRawSize () const { return GetDataRef ().GetSize () ; }
		inline t_size GetOffset () const { return m_dwOffset ; }
		inline t_size GetSize () const { return size () ; }
		inline t_size GetEnd () const { return GetEndRef () ; }

		void ResetOffsetSize ()	{ ResetOffset_NC (Bytes2Size (GetRawSize ())) ; }

		inline BOOL GetDataIntegrity () { return Size2Bytes (GetEnd ()) <= GetRawSize () ; }

		void Ref (SDataRef &ref)		{ ref.Ref (m_pDataRef) ; }
		void Ref_NDR (SDataRef &ref)	{ ref.Ref_NDR (m_pDataRef) ; }

		BOOL HasRawCap (t_size dwSize) { return Size2Bytes (dwSize) <= GetRawSize () ; }

		inline t_size Size2Bytes (t_size dwSize)	{ return dwSize * sizeof (T) ; }
		inline t_size Bytes2Size (t_size dwBytes)	{ return dwBytes / sizeof (T) ; }


		void SetOffset (t_size dwOffset, t_size dwSize)
		{
#ifdef PEDANTIC
			THROW (dwOffset + dwSize <= GetSize () + GetOffset ()) ;
#else
			THROW (HasRawCap (dwOffset + dwSize)) ;
#endif
			m_dwOffset = dwOffset ;
			SetSize_NC (dwSize) ;
		}

		void SetOffset_NC (t_size dwOffset, t_size dwSize)
		{
#ifdef PEDANTIC
			ASSERT (dwOffset + dwSize <= GetSize () + GetOffset ()) ;
#else
			ASSERT (HasRawCap (dwOffset + dwSize)) ;
#endif
			m_dwOffset = dwOffset ;
			SetSize_NC (dwSize) ;
		}

		void Offset (t_size dwOffset, t_size dwSize)
		{
#ifdef PEDANTIC
			THROW (dwOffset + dwSize <= GetSize ()) ;
#else
			THROW (HasRawCap (m_dwOffset + dwOffset + dwSize)) ;
#endif
			m_dwOffset += dwOffset ;
			SetSize_NC (dwSize) ;
		}

		void Offset_NC (t_size dwOffset, t_size dwSize)
		{
			
#ifdef PEDANTIC
			ASSERT (dwOffset + dwSize <= GetSize ()) ;
#else
			ASSERT (HasRawCap (m_dwOffset + dwOffset + dwSize)) ;
#endif
			m_dwOffset += dwOffset ;
			SetSize_NC (dwSize) ;
		}

		void SetSize (t_size dwSize)
		{
			THROW (HasRawCap (GetOffset () + dwSize)) ;
			SetSize_NC (dwSize) ;
		}

		void SetSize_NC (t_size dwSize)
		{
			ASSERT (HasRawCap (GetOffset () + dwSize)) ;
			GetSizeRef () = dwSize ;
			GetEndRef () = m_dwOffset + GetSize () ;
		}

		void ResetOffset (t_size dwSize)
		{
			THROW (HasRawCap (dwSize)) ;
			ResetOffset_NC (dwSize) ;
		}

		void ResetOffset_NC (t_size dwSize)
		{
			ASSERT (HasRawCap (dwSize)) ;
			GetOffsetRef () = 0  ;
			GetEndRef () = GetSizeRef () = dwSize ;
		}

	private:
		inline t_size &GetSizeRef	() { return m_dwSize	; }
		inline t_size &GetOffsetRef	() { return m_dwOffset	; }
		inline t_size &GetEndRef	() { return m_dwEnd		; }

		inline const t_size &GetSizeRef		() const { return m_dwSize		; }
		inline const t_size &GetOffsetRef	() const { return m_dwOffset	; }
		inline const t_size &GetEndRef		() const { return m_dwEnd		; }

		SDataRef *m_pDataRef ;

		t_size m_dwSize, m_dwOffset, m_dwEnd ;

	} ;

	template <class T>
	class SCData : protected SVData<T>
	{
		typedef SCData<T> t_this ;
		typedef SVData<T> t_base ;

	protected:
//	Conststructors  //

		SCData (const t_base &dat)													: t_base (dat)						{ }
		SCData (const t_base &dat, const t_size dwSize)								: t_base (dat, dwSize)				{ }
		SCData (const t_base &dat, const t_size dwOffset, const t_size dwSize)		: t_base (dat, dwOffset, dwSize)	{ }

	public:
		SCData ()																										{ }
		SCData (const t_this &p)													: t_base (p)						{ }
		SCData (const t_size dwSize)												: t_base (dwSize)					{ }
		SCData (T * const pData, const t_size dwSize)								: t_base (pData, dwSize)			{ }
		SCData (SDataRef_Static &ref)												: t_base (ref)						{ }
		SCData (SDataRef &ref, const t_size dwSize)									: t_base (ref, dwSize)				{ }
		SCData (SDataRef_Static &ref, const t_size dwSize)							: t_base (ref, dwSize)				{ }
		SCData (SDataRef &ref, const t_size dwOffset, const t_size dwSize)			: t_base (ref, dwOffset, dwSize)	{ }
		SCData (SDataRef_Static &ref, const t_size dwOffset, const t_size dwSize)	: t_base (ref, dwOffset, dwSize)	{ }

		SCData (SDataRef		&ref, const t_this &p)								: t_base (ref, p)					{ }	//	this is supposed to copy the data from p.		
		SCData (SDataRef_Static &ref, const t_this &p)								: t_base (ref, p)					{ }	//	this is supposed to copy the data from p.

		inline const T *GetDataEnd	()				const { return t_base::GetDataEnd	() ; }
		inline const T *GetData		()				const { return t_base::GetData		() ; }
		inline const T *GetData		(t_size dwIdx)	const { return t_base::GetData		(dwIdx) ; }
		inline const T *GetData_NC	(t_size dwIdx)	const { return t_base::GetData_NC	(dwIdx) ; }
		inline t_size size			()				const { return t_base::size			() ; }
		inline t_size GetSize		()				const { return t_this::size			() ; }

		inline const T &operator () (const t_size dwIdx) const { return *GetData (dwIdx) ; }
		inline operator const T * () const { return GetData () ; }
	};

////////////////
//	CDimCont  //
////////////////

	template <t_size P>
	class CDimCont
	{
		typedef CDimCont <P> t_this ;
	public:

		CDimCont (const t_this &dc)
		{
			memcpy (GetDims (), dc.GetDims (), sizeof (m_adwDim)) ;
		}

		CDimCont () {}

		BOOL EqualDims (const t_this &dc) const 
		{
			t_size i ;
			for (i = 0; i < P; i++)
				if (GetDim(i) != dc.GetDim (i))
					return FALSE ;
			return TRUE ;
		}

		const t_size size () const  { return DimProd () ; }

		const t_size DimProd () const
		{
			t_size ret = 1 ;
			t_size i ;
			for (i = 0; i < P; i++)
				ret *= GetDim_NC (i) ;
			return ret ;
		}

		void SetDim (const t_this &dc)
		{
			memcpy (GetDims (), dc.GetDims (), sizeof (m_adwDim)) ;
		}

		const t_this &dim () const { return *this ; }

		inline		t_size		GetDim			()			const	{					return P			; }
		inline const t_size		GetDim			(t_size p)	const	{ THROW  (p < P) ;	return m_adwDim [p] ; }
		inline const t_size		GetDim_NC		(t_size p)	const	{ ASSERT (p < P) ;	return m_adwDim [p] ; }

		inline const t_size *	GetDimPtr		(t_size p)	const	{ THROW  (p < P) ;	return m_adwDim + p ; }
		inline const t_size *	GetDimPtr_NC	(t_size p)	const	{ ASSERT (p < P) ;	return m_adwDim + p ; }
		inline const int *		GetDimPtrS_NC	(t_size p)	const	{ ASSERT (p < P) ;	return (int *) (m_adwDim + p) ; }
	protected:
		inline		 t_size &	GetDimRef		(t_size p)			{ THROW  (p < P) ;	return m_adwDim [p] ; }
		inline		 t_size &	GetDimRef_NC	(t_size p)			{ ASSERT (p < P) ;	return m_adwDim [p] ; }

		inline t_size *GetDims () { return m_adwDim ; }
		inline const t_size *GetDims () const { return m_adwDim ; }

		t_size m_adwDim [P] ;
	} ;


	typedef SVData <int> SSVecN ;
///////////////
//  Vectors  //
///////////////

	template <class T>
		class SCVec : public SCData<T>, public CDimCont<1>
	{
		typedef SCData<T> t_base ;
		typedef SCVec<T> t_this ;
		typedef CDimCont <1> tb_DimCont ;

		typedef SVData<T> tb_DataCont ;

	public:

//	Constructors  //

		SCVec (const tb_DataCont &dat)												: t_base (dat)								{ SetDim_NC (0) ; }
		SCVec (const tb_DataCont &dat, const t_size dwSize)							: t_base (dat, dwSize)						{ SetDim_NC (dwSize) ; }
		SCVec (const tb_DataCont &dat, const t_size dwOffset, const t_size dwSize)	: t_base (dat, dwOffset, dwSize)			{ SetDim_NC (dwSize) ; }

		SCVec ()																												{ SetDim_NC (0) ; }
		SCVec (const t_this &p)														: t_base (p), tb_DimCont (p)				{ /*SetDim_NC (p.size ()) ;*/ }
		SCVec (const t_size dwSize)													: t_base (dwSize)							{ SetDim_NC (dwSize) ; }
		SCVec (T * const pData, const t_size dwSize)								: t_base (pData, dwSize)					{ SetDim_NC (dwSize) ; }
		SCVec (SDataRef_Static &ref)												: t_base (ref)								{ SetDim_NC (0) ; }
		SCVec (SDataRef &ref, const t_size dwSize)									: t_base (ref, dwSize)						{ SetDim_NC (dwSize) ; }
		SCVec (SDataRef_Static &ref, const t_size dwSize)							: t_base (ref, dwSize)						{ SetDim_NC (dwSize) ; }
		SCVec (SDataRef &ref, const t_size dwOffset, const t_size dwSize)			: t_base (ref, dwOffset, dwSize)			{ SetDim_NC (dwSize) ; }
		SCVec (SDataRef_Static &ref, const t_size dwOffset, const t_size dwSize)	: t_base (ref, dwOffset, dwSize)			{ SetDim_NC (dwSize) ; }

		SCVec (SDataRef &ref,			const t_this &p)							: t_base (ref, p), tb_DimCont (p)	{ }
		SCVec (SDataRef_Static &ref,	const t_this &p)							: t_base (ref, p), tb_DimCont (p)	{ }

		SCVec (SDataRef &ref,			const tb_DimCont &p)						: t_base (ref, p.size ()), tb_DimCont (p)	{ }
		SCVec (SDataRef_Static &ref,	const tb_DimCont &p)						: t_base (ref, p.size ()), tb_DimCont (p)	{ }

		inline t_size size () const { return t_this::GetDim_NC (0) ; }

//	Data Access  //
		inline const T &operator ()	(const t_size &dwIdx)	const { return t_this::Item			(dwIdx)						; }
		inline const T &Item		(const t_size &dwIdx)	const { return t_base::Item			(dwIdx)						; }
		inline const T &Item_NC		(const t_size &dwIdx)	const { return t_base::Item_NC		(dwIdx)						; }

		inline const T *GetDataEnd	()						const { return t_base::GetDataEnd	()							; }
		inline const T *GetData		()						const { return t_base::GetData		()							; }
		inline const T *GetData		(const t_size &dwIdx)	const { return t_base::GetData 		(dwIdx)	; }
		inline const T *GetData_NC	(const t_size &dwIdx)	const { return t_base::GetData_NC	(dwIdx)	; }

		inline const T GetValue		(const t_size &dwIdx)	const { return *tb_DataCont::GetData 	(dwIdx) ; }
		inline const T GetValue_NC	(const t_size &dwIdx)	const { return *tb_DataCont::GetData_NC	(dwIdx) ; }


	protected:

		const t_this &operator = (const t_this &p ) const { THROW (0) ;}	//	this MUST never be called, as you can't change a constant matrix! Thus it's protected!

        // VT::25.06.2017
        //
		// Fix to compile on gcc-7: fix proposed by Prof. Ripley:
		// change 0 to 0U (unsigned) and comment out the const function (second line)
		// 
		inline		    t_size &nsizeRef ()		{ return t_this::GetDimRef_NC (0U) ; }
//		inline const	t_size &nsizeRef () const	{ return t_this::GetDimRef_NC (0U) ; }

		inline void SetDim (const t_size dwSize)
		{
			THROW (dwSize <= t_base::GetSize ()) ;
			t_this::nsizeRef () = dwSize ;
		}

		inline void SetDim_NC (const t_size dwSize)
		{
			ASSERT (dwSize <= t_base::GetSize ()) ;
			t_this::nsizeRef () = dwSize ;
		}

		inline void SetDim (const tb_DimCont & m)
		{
			THROW (m.DimProd () <= t_base::GetSize ()) ;
			tb_DimCont::SetDim (m) ;
		}

		inline void SetDim_NC (const tb_DimCont & m)
		{
			ASSERT (m.DimProd () <= t_base::GetSize ()) ;
			tb_DimCont::SetDim (m) ;
		}
	} ;

	template <class T> class SVVec ;
	template <class T>
		class SVec : public SCVec<T>
	{
		typedef SVec <T> t_this ;
		typedef SCVec <T> tc_this ;
		typedef SVVec <T> tv_this ;
		typedef SCVec<T> t_base ;
		typedef SVData<T> tb_DataCont ;
		typedef SCData<T> tb_DataContC ;

		typedef CDimCont <1> tb_DimCont ;

//		typedef CTempContainer <t_this> t_TempContainer ;
//		friend class CTempContainer <t_this> ;

	public:
		typedef T t_data ;

//	Constructors  //

		SVec ()																											{ }
		SVec (const t_this &p)														: t_base (p)						{ }
		SVec (const t_size dwSize)													: t_base (dwSize)					{ }
		SVec (T * const pData, const t_size dwSize)									: t_base (pData, dwSize)			{ }
		SVec (SDataRef_Static &ref)													: t_base (ref)						{ }	
		SVec (SDataRef &ref,		const t_size dwSize)							: t_base (ref, dwSize)				{ }
		SVec (SDataRef_Static &ref, const t_size dwSize)							: t_base (ref, dwSize)				{ }
		SVec (SDataRef &ref,		const t_size dwOffset, const t_size dwSize)		: t_base (ref, dwOffset, dwSize)	{ }
		SVec (SDataRef_Static &ref, const t_size dwOffset, const t_size dwSize)		: t_base (ref, dwOffset, dwSize)	{ }

		SVec (SDataRef &ref,		const tc_this &p)								: t_base (ref, p)					{ }
		SVec (SDataRef_Static &ref, const tc_this &p)								: t_base (ref, p)					{ }

		SVec (SDataRef &ref,		const tb_DimCont &p)							: t_base (ref, p)					{ }
		SVec (SDataRef_Static &ref, const tb_DimCont &p)							: t_base (ref, p)					{ }

		SVec (const tb_DataCont &dat)												: t_base (dat)						{ }
		SVec (const tb_DataCont &dat, const t_size dwSize)							: t_base (dat, dwSize)				{ }
		SVec (const tb_DataCont &dat, const t_size dwOffset, const t_size dwSize)	: t_base (dat, dwOffset, dwSize)	{ }

		t_this &operator = (const t_this &p )
		{
			tb_DimCont::operator = (p) ;
			tb_DataCont::operator = (p) ;
			return *this ;
		}

//  Re-Creation / Re-structuring  //

		void Reshape (const t_size dwSize)
		{
			t_base::Reshape (dwSize) ;
			t_this::SetDim_NC (dwSize) ;
		}

		void Reshape_NC (const t_size dwSize)
		{
			t_base::Reshape (dwSize) ;
			t_this::SetDim_NC (dwSize) ;
		}

		void Reshape (const t_size dwOffset, const t_size dwSize)
		{
			t_base::Reshape (dwOffset, dwSize) ;
			t_this::SetDim_NC (dwSize) ;
		}

		void Redim (const t_size dwSize)
		{
			t_base::Redim (dwSize) ;
			t_this::SetDim_NC (dwSize) ;
		}

		void Redim_NC (const t_size dwSize)
		{
			t_base::Redim_NC (dwSize) ;
			t_this::SetDim_NC (dwSize) ;
		}

		void Recreate (const t_size dwSize)
		{
			t_base::Recreate (dwSize) ;
			t_this::SetDim_NC (dwSize) ;
		}

		void Require (const tc_this &v) {Require (v.size ()) ; }

		void Require (const t_size dwSize)
		{
			t_base::Require (dwSize) ;
			t_this::SetDim_NC (dwSize) ;
		}

//	Cast Operators  //

		const	tb_DataCont &operator *() const	{ return *(const	tb_DataCont *) this ; }
				tb_DataCont &operator *()		{ return *(			tb_DataCont *) this ; }

		const	tv_this		&operator !() const {return * (const	tv_this*)		this ; }
				tv_this		&operator !()		{return * (			tv_this*)		this ; }

		operator T * () const { return GetData () ; }

				tc_this &Const ()		{return *(		tc_this) this ; }
		const	tc_this &Const () const {return *(const	tc_this) this ; }


//	Copy Operations  //		//	2do: move to global functions!
		void Copy_R (const tc_this &vec)
		{
			Require (vec) ;
			Copy_NC (vec) ;
		}

		void Copy (const tc_this &vec) const
		{
			THROW (this->EqualDims (vec)) ;
			Copy_NC (vec) ;
		}

		void Copy (const tb_DataContC &vec) const
		{
			ASSERT (t_this::size () == vec.size ()) ;
			Copy_NC (vec) ;
		}

		void Copy_NC (const tc_this &vec) const
		{
			ASSERT (this->EqualDims (vec)) ;
			memcpy (t_this::GetData (), vec.GetData (), t_this::GetSize () * sizeof (T)) ;
		}

		void Copy_NC (const tb_DataContC &vec) const
		{
			ASSERT (t_this::size () == vec.size ()) ;
			memcpy (t_this::GetData (), vec.GetData (), t_this::GetSize () * sizeof (T)) ;
		}

		void Copy (T const * const p, t_size n) const
		{
			THROW (n < t_this::size ()) ;
			Copy_NC (p, n) ;
		}

		void Copy_NC (T const * const p, t_size n) const
		{
			ASSERT (n < t_this::size ()) ;
			memcpy (t_this::GetData (), p, n * sizeof (T)) ;
		}

		void Copy_NC (T const * const p) const
		{
			memcpy (t_this::GetData (), p, t_this::GetSize () * sizeof (T)) ;
		}

//	Data Access  //
		inline T &operator ()	(const t_size &dwIdx)	const { return t_this::Item				(dwIdx) ; }
		inline T &Item			(const t_size &dwIdx)	const { return tb_DataCont::Item		(dwIdx) ; }
		inline T &Item_NC		(const t_size &dwIdx)	const { return tb_DataCont::Item_NC		(dwIdx) ; }

		inline T *GetDataEnd	()						const { return tb_DataCont::GetDataEnd	()							; }
		inline T *GetData		()						const { return tb_DataCont::GetData		()							; }
		inline T *GetData		(const t_size &dwIdx)	const { return tb_DataCont::GetData		(t_this::GetIdx (dwIdx))	; }
		inline T *GetData_NC	(const t_size &dwIdx)	const { return tb_DataCont::GetData_NC	(t_this::GetIdx (dwIdx))	; }

		inline void Reset (const T &v) const { t_base::Reset (v) ; }

		void Set (T *pData, t_size dwSize)
		{
			t_base::Set (pData, dwSize) ;
			t_base::SetDim_NC (dwSize) ;
		}

		const t_this GetDataRef (t_size dwStart, t_size dwEnd) const
		{
			ASSERT (dwStart <= dwEnd) ;
			ASSERT (dwEnd <= t_this::size ()) ;
			return t_this (**this, dwStart, dwEnd - dwStart) ;
		}



		protected:
//		const t_this &operator = (t_this &p ) const { THROW (FALSE) ; return NULL ; }	//	this MUST never be called, as you can't change a constant matrix/vector! Thus it's protected!
	} ;

	template <class T>
		class SVVec : public SVec <T>
	{
		typedef SVec <T> t_base ;
		typedef SVVec <T> t_this ;
		SVVec () {}	//	the constructor is private. this type cannot be created.
	public:

		t_this &operator = (const t_base &p )
		{
			t_base::operator = (p) ;
			return *this ;
		}
	} ;

//////////////
//	Matrix  //
//////////////

	template <class T>
		class SCMat : public SCData<T>, public CDimCont<2>
	{
		typedef SCData<T> t_base ;
		typedef SCMat <T> t_this ;

		typedef CDimCont <2> tb_DimCont ;
		typedef SVData<T> tb_DataCont ;

	public:

//	Constructors  //

		SCMat (const tb_DataCont &dat)																	: t_base (dat)								{ SetDim_NC (0, 0) ; }
		SCMat (const tb_DataCont &dat, const t_size dwRows, const t_size dwCols)						: t_base (dat, dwRows * dwCols)				{ SetDim_NC (dwRows, dwCols) ; }
		SCMat (const tb_DataCont &dat, const t_size dwOffset, const t_size dwRows, const t_size dwCols)	: t_base (dat, dwOffset, dwRows * dwCols)	{ SetDim_NC (dwRows, dwCols) ; }

		SCMat ()																																	{ SetDim_NC (0, 0) ; }
		SCMat (const t_this &p)																			: t_base (p), tb_DimCont (p)				{ }
		SCMat (const t_size dwRows, const t_size dwCols)												: t_base (dwRows * dwCols)					{ SetDim_NC (dwRows, dwCols) ; }
		SCMat (T * const pData, const t_size dwRows, const t_size dwCols)								: t_base (pData, dwRows * dwCols)			{ SetDim_NC (dwRows, dwCols) ; }
		SCMat (SDataRef_Static &ref)																	: t_base (ref)								{ SetDim_NC (0, 0) ; }
		SCMat (SDataRef &ref, const t_size dwRows, const t_size dwCols)									: t_base (ref, dwRows * dwCols)				{ SetDim_NC (dwRows, dwCols) ; }
		SCMat (SDataRef_Static &ref, const t_size dwRows, const t_size dwCols)							: t_base (ref, dwRows * dwCols)				{ SetDim_NC (dwRows, dwCols) ; }
		SCMat (SDataRef &ref, const t_size dwOffset, const t_size dwRows, const t_size dwCols)			: t_base (ref, dwOffset, dwRows * dwCols)	{ SetDim_NC (dwRows, dwCols) ; }
		SCMat (SDataRef_Static &ref, const t_size dwOffset, const t_size dwRows, const t_size dwCols)	: t_base (ref, dwOffset, dwRows * dwCols)	{ SetDim_NC (dwRows, dwCols) ; }

		SCMat (SDataRef &ref,			const t_this &p)												: t_base (ref, p), tb_DimCont (p)	{ }
		SCMat (SDataRef_Static &ref,	const t_this &p)												: t_base (ref, p), tb_DimCont (p)	{ }

		SCMat (SDataRef &ref,			const tb_DimCont &p)											: t_base (ref, p.size ()), tb_DimCont (p)	{ }
		SCMat (SDataRef_Static &ref,	const tb_DimCont &p)											: t_base (ref, p.size ()), tb_DimCont (p)	{ }

//	Data Access  //

		inline const T &operator ()	(const t_size &dwRow, const t_size &dwCol) const			{ return Item					(dwRow, dwCol) ; }
		inline const T &Item		(const t_size &dwRow, const t_size &dwCol) const			{ return t_base::Item	(GetIdx (dwRow, dwCol)) ; }
		inline const T &Item_NC		(const t_size &dwRow, const t_size &dwCol) const			{ return t_base::Item_NC (GetIdx (dwRow, dwCol)) ; }

		inline const T *GetDataEnd	()											const { return t_base::GetDataEnd	()						; }
		inline const T *GetData	()												const { return t_base::GetData		()						; }
		inline const T *GetData	(const t_size &dwRow, const t_size &dwCol)		const { return t_base::GetData		(GetIdx (dwRow, dwCol)) ; }
		inline const T *GetData_NC (const t_size &dwRow, const t_size &dwCol)	const { return t_base::GetData_NC	(GetIdx (dwRow, dwCol)) ; }

		inline const T GetValue		(const t_size &dwRow, const t_size &dwCol)	const { return *tb_DataCont::GetData 	(t_this::GetIdx (dwRow, dwCol)) ; }
		inline const T GetValue_NC	(const t_size &dwRow, const t_size &dwCol)	const { return *tb_DataCont::GetData_NC	(t_this::GetIdx (dwRow, dwCol)) ; }

		const SCVec<T> GetColRef (const t_size dwCol) const
		{
			ASSERT (dwCol < ncol ()) ;
			return SCVec <T> (*this, GetIdx (0, dwCol), t_this::nrow ()) ;
		}

/*		const SCMat<T> GetColsRef (t_size dwStart, t_size dwEnd) const
		{
			ASSERT (dwStart <= dwEnd) ;
			ASSERT (dwEnd < ncol ()) ;
			return SCMat <T> (*this, GetIdx (0, dwStart), t_this::nrows () * (dwEnd - dwStart + 1)) ;
		}
*/

//	Index Operations  //
		T *IncCol_NC (T*& p) const { return p += t_this::nrows () ; }
		t_size &IncCol_NC (t_size & n) const { return n += t_this::nrows () ; }
		T *DecCol_NC (T*& p) const { return p -= t_this::nrows () ; }
		t_size &DecCol_NC (t_size & n) const { return n -= t_this::nrows () ; }

//	Dim Operations  //
		inline const t_size GetIdx (const t_size &nRow, const t_size &nCol) const 
		{ return nRow + GetColInc () * nCol ; }

		inline const t_size size () const { return t_base::size () ; }

		inline const t_size nrow () const { return t_this::GetDim_NC (0) ; }
		inline const t_size ncol () const { return t_this::GetDim_NC (1) ; }
		inline const t_size *nrowPtr () const { return t_this::GetDimPtr_NC (0) ; }
		inline const t_size *ncolPtr () const { return t_this::GetDimPtr_NC (1) ; }

		inline const int *nrowPtrS () const { return (const int *) t_this::GetDimPtr_NC (0) ; }
		inline const int *ncolPtrS () const { return (const int *) t_this::GetDimPtr_NC (1) ; }
		const t_size GetColInc () const { return t_this::nrow () ; } 

		const t_size GetMinDim () const { return (nrow () < ncol ()) ? nrow () : ncol () ; }
		const t_size GetMaxDim () const { return (nrow () > ncol ()) ? nrow () : ncol () ; }
		const t_size GetDimDiff () const { return nrow () - ncol () ; }
	protected:

		const t_this &operator = (const t_this &p ) const { THROW (FALSE) ;}	//	this MUST never be called, as you can't change a constant matrix! Thus it's protected!

		void SetDim (const tb_DimCont & m)
		{
			THROW (m.DimProd () <= t_base::GetSize ()) ;
			tb_DimCont::SetDim (m) ;
		}

		void SetDim_NC (const tb_DimCont & m)
		{
			ASSERT (m.DimProd () <= t_base::GetSize ()) ;
			tb_DimCont::SetDim (m) ;
		}

		void SetDim (const t_size dwRows, const t_size dwCols)
		{
			THROW (dwRows * dwCols <= t_base::GetSize ()) ;
			t_this::nrowRef () = dwRows ;
			t_this::ncolRef () = dwCols ;
		}

		void SetDim_NC (const t_size dwRows, const t_size dwCols)
		{
			ASSERT (dwRows * dwCols <= t_base::GetSize ()) ;
			t_this::nrowRef () = dwRows ;
			t_this::ncolRef () = dwCols ;
		}

		inline t_size &nrowRef () { return t_this::GetDimRef_NC (0) ; }
		inline t_size &ncolRef () { return t_this::GetDimRef_NC (1) ; }

	} ;

////////////
//	SMat  //
////////////

template <class T> class SVMat ;

	template <class T>
		class SMat : public SCMat<T> //SMatRef<T>//, public CDataOwner<T, 2>
	{
		typedef SMat<T> t_this ;
		typedef SCMat<T> tc_this ;

		typedef SCMat<T> t_base ;
		typedef SVData<T> tb_DataCont ;

		typedef CDimCont<2> tb_DimCont ;

		typedef SVMat<T> tv_this ;

	public:
		typedef T t_data ;

//	Constructors  //

		SMat ()																																	{ }
		SMat (const t_this &p)																			: t_base (p)							{ }
		SMat (const t_size dwRows, const t_size dwCols)													: t_base (dwRows, dwCols)				{ }
		SMat (T * const pData, const t_size dwRows, const t_size dwCols)								: t_base (pData, dwRows, dwCols)		{ }
		SMat (SDataRef_Static &ref)																		: t_base (ref)							{ }
		SMat (SDataRef &ref, const t_size dwRows, const t_size dwCols)									: t_base (ref, dwRows, dwCols)			{ }
		SMat (SDataRef_Static &ref, const t_size dwRows, const t_size dwCols)							: t_base (ref, dwRows, dwCols)			{ }
		SMat (SDataRef &ref, const t_size dwOffset, const t_size dwRows, const t_size dwCols)			: t_base (ref, dwOffset, dwRows, dwCols){ }
		SMat (SDataRef_Static &ref, const t_size dwOffset, const t_size dwRows, const t_size dwCols)	: t_base (ref, dwOffset, dwRows, dwCols){ }

		SMat (SDataRef &ref,		const tc_this &p)													: t_base (ref, p)						{ }
		SMat (SDataRef_Static &ref, const tc_this &p)													: t_base (ref, p)						{ }

		SMat (SDataRef &ref,		const tb_DimCont &p)												: t_base (ref, p)						{ }
		SMat (SDataRef_Static &ref, const tb_DimCont &p)												: t_base (ref, p)						{ }

		SMat (const tb_DataCont &dat)																	: t_base (dat)							{ }
		SMat (const tb_DataCont &dat, const t_size dwRows, const t_size dwCols)							: t_base (dat, dwRows, dwCols)			{ }
		SMat (const tb_DataCont &dat, const t_size dwOffset, const t_size dwRows, const t_size dwCols)	: t_base (dat, dwOffset, dwRows, dwCols){ }

		t_this &operator = (const t_this &p )
		{
			tb_DimCont::operator = (p) ;
			tb_DataCont::operator = (p) ;
			return *this ;
		}

//  Re-Creation / Re-structuring  //

		void Reshape (const t_size dwRows, const t_size dwCols)
		{
			t_base::Reshape (dwRows * dwCols) ;
			t_this::SetDim_NC (dwRows, dwCols) ;
		}

		void Reshape_NC (const t_size dwRows, const t_size dwCols)
		{
			t_base::Reshape_NC (dwRows * dwCols) ;
			t_this::SetDim_NC (dwRows, dwCols) ;
		}

		void Reshape (const t_size dwOffset, const t_size dwRows, const t_size dwCols)
		{
			t_base::Reshape (dwOffset, dwRows * dwCols) ;
			t_this::SetDim_NC (dwRows, dwCols) ;
		}

		void Redim (const t_size dwRows, const t_size dwCols)
		{
			t_base::Redim (dwRows * dwCols) ;
			t_this::SetDim_NC (dwRows, dwCols) ;
		}

		void Redim_NC (const t_size dwRows, const t_size dwCols)
		{
			t_base::Redim_NC (dwRows * dwCols) ;
			t_this::SetDim_NC (dwRows, dwCols) ;
		}

		void Recreate (const t_size dwRows, const t_size dwCols)
		{
			t_base::Recreate (dwRows * dwCols) ;
			t_this::SetDim_NC (dwRows, dwCols) ;
		}

		void Require (const t_size dwRows, const t_size dwCols)
		{
			t_base::Require (dwRows * dwCols) ;
			t_this::SetDim_NC (dwRows, dwCols) ;
		}

		void Require (const tb_DimCont &mat)
		{
			t_base::Require (mat.DimProd ()) ;
			t_this::SetDim (mat) ;
		}

//	Data Access  //

		inline T &operator ()	(const t_size &dwRow, const t_size &dwCol)	const { return Item						(dwRow, dwCol)					; }
		inline T &operator ()	(const t_size &dwIdx)						const { return Item						(dwIdx)					; }
		inline T &Item			(const t_size &dwRow, const t_size &dwCol)	const { return tb_DataCont::Item		(t_this::GetIdx (dwRow, dwCol))	; }
		inline T &Item_NC		(const t_size &dwRow, const t_size &dwCol)	const { return tb_DataCont::Item_NC		(t_this::GetIdx (dwRow, dwCol))	; }

		inline T *GetDataEnd	()											const { return tb_DataCont::GetDataEnd	()								; }
		inline T *GetData		()											const { return tb_DataCont::GetData		()								; }
		inline T *GetData		(const t_size &dwRow, const t_size &dwCol)	const { return tb_DataCont::GetData 	(t_this::GetIdx (dwRow, dwCol)) ; }
		inline T *GetData_NC	(const t_size &dwRow, const t_size &dwCol)	const { return tb_DataCont::GetData_NC	(t_this::GetIdx (dwRow, dwCol)) ; }



		inline void Reset (const T &v) const { t_base::Reset (v) ; }

		const SVec<T> GetColRef (t_size dwCol) const
		{
			ASSERT (dwCol < t_this::ncol ()) ;
			return SVec<T> (**this, t_this::GetIdx (0, dwCol), t_this::nrow ()) ;
		}

		const t_this GetColRef (t_size dwStart, t_size dwEnd) const
		{
			ASSERT (dwStart <= dwEnd) ;
			ASSERT (dwEnd <= t_this::ncol ()) ;
			return t_this (**this, t_this::GetIdx (0, dwStart), t_this::nrow (), dwEnd - dwStart) ;
		}


//	Cast Operators  //
		const	tb_DataCont &operator *() const	{ return *(const	tb_DataCont *)	this ; }
				tb_DataCont &operator *()		{ return *(			tb_DataCont *)	this ; }

		const	tv_this		&operator !() const {return * (const	tv_this*)		this ; }
				tv_this		&operator !()		{return * (			tv_this*)		this ; }

		operator T * () const { return GetData () ; }

				tc_this &Const ()		{return *(		tc_this *) this ; }
		const	tc_this &Const () const {return *(const	tc_this *) this ; }

//	Copy Operations  //

		void CopyCol (const t_size dwDest, const t_size dwSource, const tc_this &m) const
		{
			THROW (t_this::nrow () == t_this::m.nrow ()) ;
			THROW (dwSource < t_this::m.ncol ()) ;
			THROW (dwDest < t_this::ncol ()) ;

			CopyCol_NC (dwDest, dwSource, m) ;
		}

		void CopyCol_NC (const t_size dwDest, const t_size dwSource, const tc_this &m) const
		{
			ASSERT (t_this::nrow () == m.nrow ()) ;
			ASSERT (dwSource < m.ncol ()) ;
			ASSERT (dwDest < t_this::ncol ()) ;

			//::Copy (t_this::GetData (0, dwDest), m.GetData (0, dwSource), t_this::nrow ()) ;
			memcpy (t_this::GetData () + t_this::GetIdx (0, dwDest), m.GetData () + m.GetIdx (0, dwSource), t_this::nrow () * sizeof (T)) ;
		}

		void Copy_R (const tc_this &mat)
		{
			Require (mat) ;
			Copy_NC (mat) ;
		}

		void Copy (const tc_this &mat) const
		{
			THROW (this->EqualDims (mat)) ;
			Copy_NC (mat) ;
		}

		void Copy_NC (const tc_this &mat) const
		{
			ASSERT (this->EqualDims (mat)) ;
			memcpy (t_this::GetData (), mat.GetData (), t_this::GetSize () * sizeof (T)) ;
		}

		void CopyCol_Order_R (const tc_this &mat, const SSVecN &order)	//2do: move to public scope
		{
			Require (mat.ncol (), order.GetSize ()) ;
			CopyCol_Order_NC (mat, order) ;
		}

		void CopyCol_Order (const tc_this &mat, const SSVecN &order) const
		{
			THROW (t_this::ncol () == order.GetSize ()) ;
			THROW (t_this::nrow () == mat.nrow ()) ;

			CopyCol_Order_NC (mat, order) ;
		}

		void CopyCol_Order_NC (const tc_this &mat, const SSVecN &order) const
		{
			ASSERT (t_this::ncol () == order.size ()) ;
			ASSERT (t_this::nrow () == mat.nrow ()) ;

			t_size i ;
			const int *pnData = order.GetData () ;
			for (i = order.size () - 1; i != NAI; i--)
			{
				ASSERT ((unsigned) pnData [i] < t_this::ncol ()) ;
				CopyCol_NC (i, pnData [i], mat) ;
			}
		}

		void Set (T *pData, t_size dwrow, t_size dwcol)
		{
			t_base::Set (pData, dwrow * dwcol) ;
			t_base::SetDim_NC (dwrow, dwcol) ;
		}

	protected:

	} ;


	template <class T>
		class SVMat : public SMat<T>
	{
		typedef SMat<T> t_base ;
		typedef SVMat <T> t_this ;
		SVMat () {}	//	the constructor is private. this type cannot be created.
		public:

		t_this &operator = (const t_base &p )
		{
			t_base::operator = (p) ;
			return *this ;
		}
	} ;


//	SCMatArray  //

	template <class T>
	class SCMatArray
	{
		typedef SMat <T> t_item ;
		typedef const SMat <T> tc_item ;
		typedef SCMatArray<T> t_this ;
	public:

		SCMatArray (const t_size dwRows, const t_size dwCols, const t_size dwSize)
		{
			SDataRef *pDataRef = new SDataRef (dwRows * dwCols * dwSize * sizeof (T)) ;
			FillMats_ND (pDataRef, dwRows, dwCols, dwSize) ;
		}

		SCMatArray (T *pData, const t_size dwRows, const t_size dwCols, const t_size dwSize)
		{
			SDataRef *pDataRef = new SDataRef (dwRows * dwCols * dwSize * sizeof (T), pData) ;
			FillMats_ND (pDataRef, dwRows, dwCols, dwSize) ;
		}

		SCMatArray (SDataRef_Static &ref , const t_size dwRows, const t_size dwCols, const t_size dwSize)
		{
			ref.Require (dwRows * dwCols * dwSize * sizeof (T)) ;
			FillMats_ND (&ref, dwRows, dwCols, dwSize) ;
		}

		SCMatArray (SDataRef &ref, const t_size dwRows, const t_size dwCols, const t_size dwSize)
		{
			THROW (ref.GetSize () <= dwRows * dwCols * dwSize * sizeof (T)) ;
			FillMats_ND (&ref, dwRows, dwCols, dwSize) ;
		}

		void Create (const t_size dwRows, const t_size dwCols, const t_size dwSize)
		{
			SDataRef *pDataRef = new SDataRef (dwRows * dwCols * dwSize * sizeof (T)) ;
			FillMats (pDataRef, dwRows, dwCols, dwSize) ;
		}

		void Create (T *pData, const t_size dwRows, const t_size dwCols, const t_size dwSize)
		{
			SDataRef *pDataRef = new SDataRef (dwRows * dwCols * dwSize * sizeof (T), pData) ;
			FillMats (pDataRef, dwRows, dwCols, dwSize) ;
		}		

		SCMatArray () : m_apData (NULL), m_dwSize (0) { }

		~SCMatArray ()
		{	
			t_this::Free () ;
		}

		inline tc_item &operator [] (const t_size idx) const { return Item (idx) ; }

		tc_item &Item (t_size idx)	 const { THROW (idx < t_this::GetSize ()) ; return *m_apData [idx] ; }
		tc_item &Item_NC (t_size idx) const { ASSERT	(idx < t_this::GetSize ()) ; return *m_apData [idx] ; }

		const t_size &GetSize () const { return m_dwSize ; }

	protected:
		void Free ()
		{
			int i ;
			for (i = GetSize () - 1; i != -1; i--)
				delete m_apData [i] ;
			delete [] m_apData ;
			m_apData = NULL ;
			m_dwSize = 0 ;
		}

		void FillMats (SDataRef *pDataRef, const t_size dwRows, const t_size dwCols, const t_size dwSize)
		{
			Free () ;
			FillMats_ND (pDataRef, dwRows, dwCols, dwSize) ;
		}

		void FillMats_ND (SDataRef *pDataRef, const t_size dwRows, const t_size dwCols, const t_size dwSize)
		{
			const t_size dwMatExt = dwRows * dwCols ; 
			t_size dwOffset = 0 ;
			t_size i ;
			m_dwSize = dwSize ;
			m_apData = new t_item *[dwSize] ;

			for (i = 0; i < dwSize; i++)
			{
				m_apData[i] = new t_item (*pDataRef, dwOffset, dwRows, dwCols) ;
				dwOffset += dwMatExt ;
			}			
		}

		t_item **m_apData ;
		t_size m_dwSize ;
	} ;

	typedef SCMatArray<double> SCMatArrayD ;
	typedef SCMatArray<int> SCMatArrayN ;

/////////////////
//  SMat Sort  //	basic sort routines..
/////////////////

	template <class T> 
	class CQSortComp
	{
	public:
		static int compare (const void * elem1, const void * elem2)
		{
			if (*(T *) elem1 < *(T *) elem2)
				return -1 ;
			if (*(T *) elem1 > *(T *) elem2)
				return 1 ;
			return 0 ;
		}

		static int compare_rev (const void * elem1, const void * elem2)
		{
			return -compare_rev (elem1, elem2) ;
		}

		static int compare_p (const void * elem1, const void * elem2)
		{
			return compare (* (T **)elem1, * (T **)elem2) ;
		}

		static int compare_p_rev (const void * elem1, const void * elem2)
		{
			return -compare (* (T **)elem1, * (T **)elem2) ;
		}
	} ;

	template <class T>
	void sme_qsort (T *p, t_size dwLen, BOOL bDecr = FALSE)
	{
		if (bDecr)
			qsort (p, dwLen, sizeof (T), CQSortComp<T>::compare_rev) ;	
		else
			qsort (p, dwLen, sizeof (T), CQSortComp<T>::compare) ;
	}

	template <class T>
	void sme_qsortI (T *p, int *pnIdx, t_size dwLen, BOOL bDecr = FALSE)
	{
		ASSERT_TEMPRANGE (0, 0) ;
		SDataRef_Static &tr = tempRef (0) ;

		tr.Require (sm_max (sizeof (T*), sizeof (T)) * dwLen) ;
		T **pIdx = (T **) tr.GetData () ;

		t_size i = 0 ;

		for (i = dwLen - 1; i != NAI; i--)
			pIdx[i] = p + i ;

		if (bDecr)
			qsort (pIdx, dwLen, sizeof (T *), CQSortComp<T>::compare_p_rev) ;
		else
			qsort (pIdx, dwLen, sizeof (T *), CQSortComp<T>::compare_p) ;


		for (i = dwLen - 1; i != NAI; i--)
			pnIdx[i] = pIdx[i] - p ;

		//tr.Require (sizeof (T) * dwLen) ,
		T *pBuf = (T *) tr.GetData () ;
		memcpy (pBuf, p, sizeof (T) * dwLen) ;
		
		for (i = dwLen - 1; i != NAI; i--)
			p[i] = pBuf[pnIdx[i]] ;

	}


	template <class T>
		t_size which_max1 (T const *p, t_size n)
	{
		T const * const pEnd = p + n ;
		T max = *p ;
		T const * pMax = p, *pCur = p + 1;
		for (; pCur < pEnd; ++pCur)
			if (sm_setmax_b (max, *pCur))
				pMax = p ;
		return pMax - p ;

	}

	template <class T>
		t_size which_max1 (const SCData <T> &a)
	{
		return which_max1 (a.GetData (), a.size ()) ;
	}

#endif	//#ifndef ES_SMAT_BASE_H


