#pragma once

#include "IG_Enumerator.h"
#define NEW_CODE    0
Class2Def(CEnumerator);
Class2Def(CEnumInfo);

Class2Def(CThreadEnumerator)
{
public:
	CK CThreadEnumerator()							{ m_pEnum = NULL; reset();
#if WRITE_MULTITHREAD_LOG
		extern int threadCntr;
		m_threadID = threadCntr++;
#endif
	}
	CK ~CThreadEnumerator()							{ release(); }
	CK void setupThreadForBIBD(const EnumeratorPntr pMaster, T nRows, int threadIdx);
	void EnumerateBIBD(designParam *pParam, EnumeratorPntr pMaster);
	CK inline t_threadCode code() const				{ return m_code; }
	CK inline EnumeratorPntr enumerator() const		{ return m_pEnum; }
	CK inline EnumInfoPntr enumInfo() const			{ return enumerator()? enumerator()->enumInfo() : NULL; }
#if NEW_CODE
	CK inline void reInit()							{ reset(); }
#else
	CK inline void reInit()							{ release(); reset(); }
#endif
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
	void LaunchCanonicityTesting(const EnumeratorPntr pEnum)
		{ enumerator()->CanonCheckerGPU()->LaunchCanonicityTesting(enumInfo(), pEnum); }
#endif
private:
	CK inline void release()						{ delete enumerator(); }
#if NEW_CODE
	CK inline void reset()							{ setCode(t_threadNotUsed); }
#else
	CK inline void reset()							{ m_pEnum = NULL; setCode(t_threadNotUsed); }
#endif

	EnumeratorPntr m_pEnum;
	t_threadCode m_code;
#if USE_BOOST && USE_POOL
	boost::thread *m_pTread;
#endif
#if WRITE_MULTITHREAD_LOG
	int m_threadID;
#endif
};

Class2Def(CThreadEnumPool)
{
 public:
	CK CThreadEnumPool(Class2(CThreadEnumerator) * *pntr, size_t numb) : m_ppThreadEnums(pntr), m_numbMax(numb) {}
	CK inline void pushToPool(ThreadEnumeratorPntr pEnum)	{ assert(m_numbUsed < m_numbMax);
															  m_ppThreadEnums[m_numbUsed++] = pEnum;
															}
	CK inline size_t poolSize() const						{ return m_numbUsed; }
	CK ThreadEnumeratorPntr popFromPool()					{ return m_ppThreadEnums[--m_numbUsed]; }
private:
	Class2(CThreadEnumerator)** m_ppThreadEnums = NULL;
	const size_t m_numbMax = 0;
	size_t m_numbUsed = 0;

};

#if USE_THREADS
#include "C_tDesignEnumerator.h"
#include "PBIBD_Enumerator.h"
#include "CombBIBD_Enumerator.h"
#include "EnumInfo.h"

Class2Def(C_InSysEnumerator);
Class2Def(C_tDesignEnumerator);
Class2Def(CBIBD_Enumerator);
Class2Def(CCombBIBD_Enumerator);

FClass2(CThreadEnumerator, void)::setupThreadForBIBD(const EnumeratorPntr pMaster, T nRow, int threadIdx) {
	if (pMaster->IS_enumerator()) {
		if (!m_pEnum) {
		auto *pInSys = pMaster->matrix();
		InSysPntr pSlaveDesign = NULL;
		const uint enumFlags = pMaster->enumFlags() | t_matrixOwner;
		switch (pInSys->objectType()) {
		case t_objectType::t_BIBD:
			    pSlaveDesign = new C_BIBD<T, S>((const C_BIBD<T, S> *)(pInSys), nRow);
				m_pEnum = new CBIBD_Enumerator<T, S>(pSlaveDesign, enumFlags, threadIdx, NUM_GPU_WORKERS);
				break;
		case t_objectType::t_Kirkman_Triple:
		case t_objectType::t_CombinedBIBD:
				pSlaveDesign = new CCombinedBIBD<T, S>((const Class2(CCombinedBIBD) *)(pInSys), nRow);
				m_pEnum = new CCombBIBD_Enumerator<T, S>(pSlaveDesign, enumFlags, threadIdx, NUM_GPU_WORKERS);
				break;
		case t_objectType::t_tDesign: {
				pSlaveDesign = new C_tDesign<T, S>((const Class2(C_tDesign) *)(pInSys), nRow);
				m_pEnum = new C_tDesignEnumerator<T, S>(pSlaveDesign, enumFlags, threadIdx, NUM_GPU_WORKERS);
				break;
			}
		case t_objectType::t_PBIBD:
				pSlaveDesign = new C_PBIBD<T, S>((const Class2(C_PBIBD) *)(pInSys), nRow);
				m_pEnum = new CPBIBD_Enumerator<T, S>(pSlaveDesign, enumFlags, threadIdx, NUM_GPU_WORKERS);
				break;
		case t_objectType::t_SemiSymmetricGraph:
			    pSlaveDesign = new CSemiSymmetricGraph<T, S>((const Class2(CSemiSymmetricGraph) *)(pInSys), nRow);
				m_pEnum = new CIG_Enumerator<T, S>(pSlaveDesign, pMaster->designParams(), enumFlags, threadIdx, NUM_GPU_WORKERS);
				m_pEnum->CloneMasterInfo(pMaster, nRow);
				break;
		}

		pSlaveDesign->setObjectType(pInSys->objectType());
		m_pEnum->setEnumInfo(new Class2(CInsSysEnumInfo)());
		} else {
			// Releasing some memory if needed
			m_pEnum->releaseRowStuff(nRow);
			((InSysPntr)(m_pEnum->matrix()))->DuplicateMasterMatrix(pMaster->matrix(), nRow);
			nRow = 0;
		}
	}

	setCode(t_threadUndefined);
}

FClass2(CThreadEnumerator, void)::EnumerateBIBD(designParam *pParam, EnumeratorPntr pMaster) {
	THREAD_MESSAGE(threadID(), -1, "threadEnumerate START", code(), m_pEnum);
	m_pEnum->Enumerate(pParam, false, m_pEnum->enumInfo(), pMaster, &m_code);
	THREAD_MESSAGE(threadID(), -1, "threadEnumerate DONE", code());
}
#endif
