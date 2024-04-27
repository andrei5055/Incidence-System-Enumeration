#include "TopGun.h"

void RunThread(int threadNumber, int iMode,
	char* mstart0, char* mstart, int nRowsStart, int nRowsOut, sLongLong* pcnt, bool bPrint)
{
	alldata sys(nPlayers);
	sys.Run(threadNumber, iMode, mstart0, mstart, nRowsStart, nRowsOut, pcnt, bPrint);
}
void TopGun:: waitAllThreadFinished()
{
	int i = 0;
	for (auto& t : threads) {
		if (threadActive[i])
		{
			t.join();
			if (cnt[i * 2] < 0)
				printf("\n*** Internal error in thread %d (%zd) ***\n", i + 1, cnt[i * 2]);
		}
		i++;
	}
}
void TopGun::threadStopped(int iTask)
{
	cntTotal[iTask * 2] += cnt[iTask * 2];
	cntTotal[iTask * 2 + 1] += cnt[iTask * 2 + 1];
	cnt[iTask * 2] = cnt[iTask * 2 + 1] = 0;
	threadActive[iTask] = false;
	//printf("thread %d stopped\n", iTask + 1);
}
void TopGun::startThread(int iTask)
{
	cnt[iTask * 2] = -1;
	cnt[iTask * 2 + 1] = 0;
	threads[iTask] = std::thread{ RunThread, ++iTaskSeq, eCalcResult,
		mstart, mstart, nRowsStart, nRowsOut, &cnt[iTask * 2], false };
	threadActive[iTask] = true;
#if 0
	printf("*** Thread %d ", iTask + 1);
	printTable("Start matrix", mstart, nRowsStart, m_numPlayers, 0, m_groupSize, true);
#endif
	mstart += mStartMatrixSize;
	iMatrix += 1;
}
#if 1
int TopGun::getStartMatrices()
{
	int nmsAll = 0;
	int nf = 1;
	int nfr = 0;
	char* s = startMatrix;
	int n = nMatricesMax;
	int nerr = 0;
	int nms = 0;
	while (1)
	{
		char logfn[256];
		createFolderAndFileName(logfn, sizeof(logfn), StartFolder, ResultNameFormat,
			m_numPlayers, nRowsStart, m_groupSize, nf++);
		if ((nms = readStartData(logfn, s, n, m_numPlayers, nRowsStart, m_groupSize)) == 0)
		{
			if (nerr == 0)
				printf("\n");
			printf("\rCant open file with 'Start' matrices: %s", logfn);
			if (++nerr > 100)
			    break;
			continue;
		}
		nerr = 0;
		nfr++;
		printf("\n%d %d-rows 'Start' matrices loaded from file %s", nms, nRowsStart, logfn);
		nmsAll += nms;
		s += nms * mStartMatrixSize;
		n -= nms;
	}
	nMatrices = nmsAll;
	printf("\n%d %d-rows 'Start' matrices loaded from %d files\n", nMatrices, nRowsStart, nfr);
	//exit(0);
	/**
	createStartFolderAndFileName(startMatricesFullFileName, sizeof(startMatricesFullFileName), StartMatricesFolder, startMatricesFileName,
		m_numPlayers, nRowsStart, m_groupSize);
	
	RunThread(0, eCalcStart, ImproveResults, (char*)0, startMatrix, nRowsStart, 0, dNumMatrices, true);
	nMatrices = (int)dNumMatrices[0];
	if (nMatrices > nMatricesMax)
	{
		printf("Number of initial %d-rows matrices > limit (%d). Exit\n", nRowsStart, nMatricesMax);
		exit(1);
	}

	saveStartData(startMatricesFullFileName, startMatrix, nMatrices, m_numPlayers, nRowsStart, m_groupSize);
	printf("%d %d-rows 'Start' matrices saved to file %s\n", nMatrices, nRowsStart, startMatricesFullFileName);
	exit(0);**/
	return nMatrices;
}
#else
int TopGun::getStartMatrices()
{
	createStartFolderAndFileName(startMatricesFullFileName, sizeof(startMatricesFullFileName), StartMatricesFolder, startMatricesFileName,
		m_numPlayers, nRowsStart, m_groupSize);

	if ((nMatrices = readStartData(startMatricesFullFileName, startMatrix, nMatricesMax, m_numPlayers, nRowsStart, m_groupSize)) > 0)
	{
		printf("%d %d-rows 'Start' matrices loaded from file %s\n", nMatrices, nRowsStart, startMatricesFullFileName);
	}
	else
	{
		RunThread(0, eCalcStart, ImproveResults, (char*)0, startMatrix, nRowsStart, 0, dNumMatrices, true);
		nMatrices = (int)dNumMatrices[0];
		if (nMatrices > nMatricesMax)
		{
			printf("Number of initial %d-rows matrices > limit (%d). Exit\n", nRowsStart, nMatricesMax);
			exit(1);
		}
		saveStartData(startMatricesFullFileName, startMatrix, nMatrices, m_numPlayers, nRowsStart, m_groupSize);
		printf("%d %d-rows 'Start' matrices saved to file %s\n", nMatrices, nRowsStart, startMatricesFullFileName);
	}
}
#endif
void TopGun::myTemporaryCheck()
{
	char* lnki = mStartLinks;
	char* stri = startMatrix;

	int nm = nMatrices < 20000 ? nMatrices : 20000;
	int is = 1;//nMatrices / nm;
	int mstep = mStartMatrixSize* is;
	int lstep = mLinksSize * is;

	for (int i = 0; i < nm; i++, stri += mstep, lnki += lstep)
	{
		linksFromMatrix(lnki, stri, nRowsStart, m_numPlayers);
	}
	for (int iRow = 3; iRow <= nRowsStart && iRow < m_numDays; iRow++)
	{
		lnki = mStartLinks;
		stri = startMatrix;
		for (int i = 0; i < nm - 1; i++, stri += mstep, lnki += lstep)
		{
			if ((i % 100) == 0)
				printf(" %d%%\r", (i + 2) * 100 / nm);
			if (*lnki == 0)
				continue;
			char lnkt[nPlayers * nPlayers];
			for (int k = 0; k < mLinksSize; k++)
			{
				char a = lnki[k];
				lnkt[k] = (a == unset) ? a : (a >= iRow ? unset : 0);
			}
			char* strj = startMatrix;
			char* lnkj = mStartLinks;
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
						exit(0);
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
	//exit(0);
}