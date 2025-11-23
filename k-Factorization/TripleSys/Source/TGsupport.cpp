#include "TopGun.h"

void RunThread(int threadId, eThreadStartMode iMode,
	TopGun *pMaster, CStorageSet<tchar>* secondRowsDB, ctchar* mstart0, ctchar* mfirst, sLongLong* pcnt, int iThread, CRowStorage *pRowStorage)
{
	alldata sys(*pMaster, pMaster->paramPtr(), 0, pRowStorage);
	sys.Run(threadId, iMode, secondRowsDB, mstart0, mfirst, pMaster->nRowsStart(), pcnt, 0, iThread);
	pMaster->updateMatrixDB(sys.matrixDB());
}
void TopGun::waitAllThreadFinished()
{
	int i = 0;
	for (auto& t : threads) {
		if (threadActive[i])
		{
			t.join();
			if (m_cnt[i * 2] < 0)
				printfRed("\n*** Internal error in thread %d (%zd) ***\n", i + 1, m_cnt[i * 2]);
		}
		i++;
	}
}
void TopGun::threadStopped(int iTask)
{
	const auto idx = iTask * 2;
	m_cntTotal[idx] += m_cnt[idx];
	m_cntTotal[idx + 1] += m_cnt[idx + 1];
	m_cnt[idx] = m_cnt[idx + 1] = 0;
	threadActive[iTask] = false;
	//printf("thread %d stopped\n", iTask + 1);
}
void TopGun::startThread(int iTask, int iTaskId, eThreadStartMode iMode, CRowStorage* pRowStorage)
{
	assert(this != nullptr);
	assert(mstart != nullptr);
	assert(mfirst != nullptr);
	assert(m_cnt != nullptr);
	assert(iMode == eCalculateMatrices && pRowStorage != nullptr || iMode != eCalculateMatrices && pRowStorage == nullptr);

	m_cnt[iTask * 2] = -1;
	m_cnt[iTask * 2 + 1] = 0;
	threads[iTask] = std::thread{ RunThread, iTaskId, iMode,
		this, m_pSecondRowsDB, mstart, mfirst, m_cnt + iTask * 2, iTask, pRowStorage};
	threadActive[iTask] = true;
#if 0
	printfRed("*** Thread %d ", iTask + 1);
	printTable("Start matrix", mstart, nRowsStart, m_numPlayers, 0, m_groupSize, true);
#endif
}

