#pragma once

template<class T> class CEnumerator;
template<class T> class CEnumInfo;

template<class T> class CThreadEnumerator
{
public:
	CK CThreadEnumerator()							{ reset(); }
	CK ~CThreadEnumerator()							{ release(); }
	CK void setupThreadForBIBD(const CEnumerator<T> *pMaster, size_t nRow, int threadIdx);
	void EnumerateBIBD(designRaram *pParam, const CEnumerator<T> *pMaster);
	CK inline t_threadCode code() const				{ return m_code; }
	CK inline CEnumerator<T> *enumerator() const	{ return m_pEnum; }
	CK inline CEnumInfo<T> *enumInfo() const		{ return enumerator()? enumerator()->enumInfo() : NULL; }
	CK inline void reInit()							{ release(); reset(); }
#if WRITE_MULTITHREAD_LOG
	inline void setThreadID(int id)					{ m_threadID = id; }
	inline int threadID() const					    { return m_threadID; }
#endif
	CK inline void setCode(t_threadCode code)		{ m_code = code; }
#if USE_BOOST && USE_POOL
	inline void setThread(boost::thread *pThread)	{ m_pTread = pThread; }
	inline boost::thread *getThread() const			{ return m_pTread; }
#endif

#if !CONSTR_ON_GPU
	void LaunchCanonicityTesting(const CEnumerator<T> *pEnum)
		{ enumerator()->GPU_CanonChecker()->LaunchCanonicityTesting(enumInfo(), pEnum); }
#endif
private:
	CK inline void release()						{ delete enumerator(); }
	CK inline void reset()							{ m_pEnum = NULL; setCode(t_threadNotUsed); }

	CEnumerator<T> *m_pEnum;
	t_threadCode m_code;
#if USE_BOOST && USE_POOL
	boost::thread *m_pTread;
#endif
#if WRITE_MULTITHREAD_LOG
	int m_threadID;
#endif
};

#if USE_THREADS
#include "C_tDesignEnumerator.h"
#include "PBIBD_Enumerator.h"
#include "EnumInfo.h"
template<class T> class C_InSysEnumerator;
template<class T> class C_tDesignEnumerator;
template<class T> class CBIBD_Enumerator;

template<class T>
void CThreadEnumerator<T>::setupThreadForBIBD(const CEnumerator<T> *pMaster, size_t nRow, int threadIdx)
{
	if (pMaster->IS_enumerator()) {
		auto *pInsSysEnum = (const C_InSysEnumerator<T> *)(pMaster);
		auto *pInSys = static_cast<const C_InSys<T> *>(pInsSysEnum->matrix());
		if (pInsSysEnum->isTDesign_enumerator(3)) {
			auto pSlaveDesign = new C_tDesign<T>((const C_tDesign<T> *)(pInSys), nRow);
			m_pEnum = new C_tDesignEnumerator<T>(pSlaveDesign, true, pMaster->noReplicatedBlocks(), threadIdx, NUM_GPU_WORKERS);
		}
		else {
			if (pInsSysEnum->isTDesign_enumerator(2)) {
				if (pInsSysEnum->isPBIB_enumerator()) {
					auto pSlaveBIBD = new C_PBIBD<T>((const C_PBIBD<T> *)(pInSys), nRow);
					m_pEnum = new CPBIBD_Enumerator<T>(pSlaveBIBD, true, pMaster->noReplicatedBlocks(), threadIdx, NUM_GPU_WORKERS);
				}
				else {
					auto pSlaveBIBD = new C_BIBD<T>((const C_BIBD<T> *)(pInSys), nRow);
					m_pEnum = new CBIBD_Enumerator<T>(pSlaveBIBD, true, pMaster->noReplicatedBlocks(), threadIdx, NUM_GPU_WORKERS);
				}
			}
		}

		m_pEnum->setEnumInfo(new CInsSysEnumInfo<T>());
	}

	setCode(t_threadUndefined);
}

template<class T>
void CThreadEnumerator<T>::EnumerateBIBD(designRaram *pParam, const CEnumerator<T> *pMaster)
{
	thread_message(threadID(), "threadEnumerate START", code(), m_pEnum);
	m_pEnum->Enumerate(pParam, false, m_pEnum->enumInfo(), pMaster, &m_code);
	thread_message(threadID(), "threadEnumerate DONE", code());
}
#endif