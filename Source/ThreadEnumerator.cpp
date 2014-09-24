#include "C_tDesignEnumerator.h"

#if USE_THREADS
#include "EnumInfo.h"

void CThreadEnumerator::release()
{
	delete enumerator();
	delete enumInfo();
}

void CThreadEnumerator::reset()
{
	m_pEnum = NULL;
	m_pInfo = NULL;
	setCode(t_threadNotUsed);
}

void CThreadEnumerator::setupThreadForBIBD(const CEnumerator *pMaster, size_t nRow)
{
	if (pMaster->isIS_enumerator()) {
		const C_InSysEnumerator *pInsSysEnum = dynamic_cast<const C_InSysEnumerator *>(pMaster);
		const C_InSys *pInSys = dynamic_cast<const C_InSys *>(pInsSysEnum->matrix());
		if (pInsSysEnum->isTDesign_enumerator(3)) {
			C_tDesign *pSlaveDesign = new C_tDesign(dynamic_cast<const C_tDesign *>(pInSys), nRow);
			m_pEnum = new C_tDesignEnumerator(pSlaveDesign, true);
		}
		else
		if (pInsSysEnum->isTDesign_enumerator(2)) {
			C_BIBD *pSlaveBIBD = new C_BIBD(dynamic_cast<const C_BIBD *>(pInSys), nRow);
			m_pEnum = new CBIBD_Enumerator(pSlaveBIBD, true);
		}
		m_pInfo = new CInsSysEnumInfo();
	}

	setCode(t_threadUndefined);
}

void CThreadEnumerator::EnumerateBIBD(const CEnumerator *pMaster)
{
	thread_message(threadID(), "threadEnumerate START", code(), m_pEnum);
	m_pEnum->Enumerate(false, m_pInfo, pMaster, &m_code);
	thread_message(threadID(), "threadEnumerate DONE", code());
}
#endif