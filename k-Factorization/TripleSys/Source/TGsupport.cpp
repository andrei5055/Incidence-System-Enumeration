#include "TopGun.h"

void RunThread(int threadNumber, int iMode,
	TopGun *pMaster, tchar* mstart0, tchar* mstart, sLongLong* pcnt, bool firstThread)
{
	alldata sys(*pMaster, pMaster->paramPtr());
	sys.Run(threadNumber, iMode, mstart0, mstart, pMaster->nRowsStart(), pMaster->nRowsOut(), pcnt, 0, firstThread);
	if (pMaster->descrStorage())
		pMaster->addMatrixDB(sys.matrixDB());
	else
		pMaster->transferMatrixDB(sys.matrixDB());
}
void TopGun:: waitAllThreadFinished()
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
void TopGun::startThread(int iTask)
{
	m_cnt[iTask * 2] = -1;
	m_cnt[iTask * 2 + 1] = 0;
	threads[iTask] = std::thread{ RunThread, ++m_iTaskSeq, eCalcResult,
		this, mstart, mstart, m_cnt + iTask * 2, iTask == 0};
	threadActive[iTask] = true;
#if 0
	printfRed("*** Thread %d ", iTask + 1);
	printTable("Start matrix", mstart, nRowsStart, m_numPlayers, 0, m_groupSize, true);
#endif
	mstart += mStartMatrixSize;
	m_iMatrix += 1;
}

void TopGun::myTemporaryCheck()
{
	tchar* mStartLinks; // needed only for this test function
	mStartLinks = new tchar[nMatricesMax() * mLinksSize];
	if (!mStartLinks) {
		printfRed("*** Test: Not enough memory for initial %d-rows %d links. Exit\n", m_numPlayers, nMatricesMax());
		myExit(1);
	}
	auto* lnki = mStartLinks;
	auto* stri = startMatrix;

	int nm = nMatrices < 20000 ? nMatrices : 20000;
	int is = 1;//nMatrices / nm;
	int mstep = mStartMatrixSize* is;
	int lstep = mLinksSize * is;
	char* lnkt = new char[mLinksSize];

	for (int i = 0; i < nm; i++, stri += mstep, lnki += lstep)
	{
		linksFromMatrix(lnki, stri, nRowsStart());
	}
	for (int iRow = 3; iRow <= nRowsStart() && iRow < m_numDays; iRow++)
	{
		lnki = mStartLinks;
		stri = startMatrix;
		for (int i = 0; i < nm - 1; i++, stri += mstep, lnki += lstep)
		{
			if ((i % 100) == 0)
				printf(" %d%%\r", (i + 2) * 100 / nm);
			if (*lnki == 0)
				continue;

			for (int k = 0; k < mLinksSize; k++)
			{
				char a = lnki[k];
				lnkt[k] = (a == unset) ? a : (a >= iRow ? unset : 0);
			}
			auto* strj = startMatrix;
			auto* lnkj = mStartLinks;
			for (int j = i + 1; j < nm; j++, strj += mstep, lnkj += lstep)
			{
				if (*lnkj == 0)
					continue;
#if 1
				bool bSame = (memcmp(stri + m_numPlayers * iRow, strj + m_numPlayers * iRow,
					m_numPlayers * (m_numDays - iRow)) == 0);
#else
				bool bSame = true;
				for (int k = 0; k < mLinksSize; k++)
				{
					char b = lnkj[k];
					b = (b == unset) ? b : (b >= iRow ? unset : 0);
					if (lnkt[k] != b)
					{
						bSame = false;
						break;
					}
				}
#endif
#if 0
				if (bSame && iRow == nRowsStart)
				{
					printTableColor("Links a", lnki, m_numPlayers, m_numPlayers);
					printTableColor("Links b", lnkj, m_numPlayers, m_numPlayers);
					printTableColor("r i", stri, nRowsStart, m_numPlayers);
					printTableColor("r j", strj, nRowsStart, m_numPlayers);
				}
#endif
				if (bSame && memcmp(stri, strj, m_numPlayers * iRow) != 0)
				{
#if 0
					if (iRow == nRowsStart)
					{
						printTableColor("Links a", lnki, m_numPlayers, m_numPlayers);
						printTableColor("Links b", lnkj, m_numPlayers, m_numPlayers);
						printTableColor("r i", stri, nRowsStart, m_numPlayers);
						printTableColor("r j", strj, nRowsStart, m_numPlayers);
						myExit(0);
					}
#endif
					*(lnkj) = 0;
				}
			}
		}
		int cnti = 0;
		lnki = mStartLinks;
		for (int i = 0; i < nm; i++, lnki += lstep)
		{
			if (*lnki == 0)
			{
				cnti++;
				*lnki = unset;
			}
		}

		printf("%.1f%% improvement (%dx%dx%d)\n", cnti * 100.0 / nm, m_numPlayers, iRow, m_groupSize);
	}
	delete[] mStartLinks;
	mStartLinks = NULL;
	delete[] lnkt;
	//myExit(0);
}