#pragma once

#include "IG_Enumerator.h"

IClass2Def(Enumerator);
IClass2Def(EnumInfo);

IClass2Def(ThreadEnumerator)
{
public:
	CK CThreadEnumerator()							{ reset(); }
	CK ~CThreadEnumerator()							{ release(); }
	CK void setupThreadForBIBD(const EnumeratorPntr pMaster, size_t nRow, int threadIdx);
	void EnumerateBIBD(designParam *pParam, const EnumeratorPntr pMaster);
	CK inline t_threadCode code() const				{ return m_code; }
	CK inline EnumeratorPntr enumerator() const		{ return m_pEnum; }
	CK inline EnumInfoPntr enumInfo() const			{ return enumerator()? enumerator()->enumInfo() : NULL; }
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
	void LaunchCanonicityTesting(const EnumeratorPntr pEnum)
		{ enumerator()->GPU_CanonChecker()->LaunchCanonicityTesting(enumInfo(), pEnum); }
#endif
private:
	CK inline void release()						{ delete enumerator(); }
	CK inline void reset()							{ m_pEnum = NULL; setCode(t_threadNotUsed); }

	EnumeratorPntr m_pEnum;
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
#include "CombBIBD_Enumerator.h"
#include "EnumInfo.h"

IClass2Def(_InSysEnumerator);
IClass2Def(_tDesignEnumerator);
IClass2Def(BIBD_Enumerator);
IClass2Def(CombBIBD_Enumerator);

TClass2(ThreadEnumerator, void)::setupThreadForBIBD(const EnumeratorPntr pMaster, size_t nRow, int threadIdx) {
	if (pMaster->IS_enumerator()) {
		auto *pInsSysEnum = (const IClass2(_InSysEnumerator) *)(pMaster);
		auto *pInSys = static_cast<const InSysPntr>(pInsSysEnum->matrix());
		const InSysPntr pSlaveDesign;
		const uint enumFlags = pMaster->enumFlags() | t_matrixOwner;
		switch (pInSys->objectType()) {
		case t_BIBD:
			    pSlaveDesign = new C_BIBD<T, S>((const C_BIBD<T, S> *)(pInSys), nRow);
				m_pEnum = new CBIBD_Enumerator<T, S>(pSlaveDesign, enumFlags, threadIdx, NUM_GPU_WORKERS);
				break;
		case t_CombinedBIBD:
				pSlaveDesign = new CCombinedBIBD<T, S>((const IClass2(CombinedBIBD) *)(pInSys), nRow);
				m_pEnum = new CCombBIBD_Enumerator<T, S>(pSlaveDesign, enumFlags, threadIdx, NUM_GPU_WORKERS);
				break;
		case t_tDesign: {
				const auto pSlaveTDesign = new C_tDesign<T, S>((const IClass2(_tDesign) *)(pInSys), nRow);
				m_pEnum = new C_tDesignEnumerator<T, S>(pSlaveTDesign, enumFlags, threadIdx, NUM_GPU_WORKERS);
				break;
			}
		case t_PBIBD:
				pSlaveDesign = new C_PBIBD<T, S>((const IClass2(_PBIBD) *)(pInSys), nRow);
				m_pEnum = new CPBIBD_Enumerator<T, S>(pSlaveDesign, enumFlags, threadIdx, NUM_GPU_WORKERS);
				break;
		case t_SemiSymmetricGraph:
			    pSlaveDesign = new CSemiSymmetricGraph<T, S>((const IClass2(SemiSymmetricGraph) *)(pInSys), nRow);
				m_pEnum = new CIG_Enumerator<T, S>(pSlaveDesign, pMaster->designParams(), enumFlags, threadIdx, NUM_GPU_WORKERS);
				m_pEnum->CloneMasterInfo(pMaster, nRow);
				break;
		}

		m_pEnum->setEnumInfo(new IClass2(InsSysEnumInfo)());
	}

	setCode(t_threadUndefined);
}

TClass2(ThreadEnumerator, void)::EnumerateBIBD(designParam *pParam, const EnumeratorPntr pMaster) {
	thread_message(threadID(), "threadEnumerate START", code(), m_pEnum);
	m_pEnum->Enumerate(pParam, false, m_pEnum->enumInfo(), pMaster, &m_code);
	thread_message(threadID(), "threadEnumerate DONE", code());
}
#endif