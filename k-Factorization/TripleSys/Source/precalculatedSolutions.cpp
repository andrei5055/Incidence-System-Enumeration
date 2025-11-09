#include <iostream>
#include "TripleSys.h"
#include "SRGToolkit.h"

#include <mutex>
extern std::mutex mtxLinks; // The mutex to protect the shared resource
extern CStorageIdx<tchar>** mpLinks;
extern int SemiPhase;
extern int NumMatricesProcessed;
#define CheckMissingMatrix 0
#if CheckMissingMatrix
tchar missingMatrix[] = {
  0,  1,   2,  3,   4,  5,   6,  7,   8,  9,  10, 11,  12, 13, 
  0,  2,   1,  4,   3,  6,   5,  8,   7, 10,   9, 12,  11, 13, 
  0,  3,   1,  5,   2,  7,   4,  9,   6, 12,   8, 11,  10, 13, 
  0,  4,   1,  7,   2, 11,   3, 12,   5,  9,   6, 10,   8, 13, 
  0,  5,   1, 12,   2, 13,   3, 10,   4,  6,   7,  8,   9, 11, 
  0,  6,   1, 11,   2,  9,   3,  4,   5, 12,   7, 13,   8, 10, 
  0,  7,   1, 10,   2,  5,   3, 11,   4, 13,   6,  9,   8, 12, 
  0,  8,   1,  2,   3, 13,   4, 10,   5,  6,   7,  9,  11, 12, 
  0,  9,   1, 13,   2,  8,   3,  5,   4,  7,   6, 11,  10, 12, 
  0, 10,   1,  6,   2,  4,   3,  8,   5, 11,   7, 12,   9, 13, 
  0, 11,   1,  3,   2, 12,   4,  8,   5,  7,   6, 13,   9, 10, 
  0, 12,   1,  9,   2, 10,   3,  7,   4, 11,   5, 13,   6,  8, 
  0, 13,   1,  8,   2,  6,   3,  9,   4, 12,   5, 10,   7, 11
};
tchar canonicalLinks[] = {
  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  0,  0, 
  1,  0,  1,  1,  1,  1,  1,  1,  1,  1,  0,  1,  0,  0, 
  1,  1,  0,  1,  1,  1,  1,  1,  1,  0,  0,  0,  1,  1, 
  1,  1,  1,  0,  1,  1,  1,  0,  0,  1,  1,  1,  1,  0, 
  1,  1,  1,  1,  0,  1,  0,  1,  0,  1,  1,  0,  1,  1, 
  1,  1,  1,  1,  1,  0,  0,  0,  1,  0,  1,  1,  1,  1, 
  1,  1,  1,  1,  0,  0,  0,  1,  0,  1,  1,  1,  1,  1, 
  1,  1,  1,  0,  1,  0,  1,  0,  1,  1,  0,  1,  1,  1, 
  1,  1,  1,  0,  0,  1,  0,  1,  0,  1,  1,  1,  1,  1, 
  1,  1,  0,  1,  1,  0,  1,  1,  1,  0,  1,  1,  0,  1, 
  1,  0,  0,  1,  1,  1,  1,  0,  1,  1,  0,  1,  1,  1, 
  0,  1,  0,  1,  0,  1,  1,  1,  1,  1,  1,  0,  1,  1, 
  0,  0,  1,  1,  1,  1,  1,  1,  1,  0,  1,  1,  0,  1, 
  0,  0,  1,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0, 
};
#endif
int alldata::addMaskToDB(ll* msk, tchar* r, int nr, int nc, int np) {
	msk[0] = msk[1] = 0;
	memset(m_pGraph, 0, nc * nc);
	
	for (tchar i = 0; i < (nr * nc); i += 2) {
		auto const a = r[i], b = r[i + 1];
		m_pGraph[b * nc + a] = m_pGraph[a * nc + b] = 1;
	}

#if CheckMissingMatrix
	if (memcmp(missingMatrix, r, nr * nc) == 0) {
		printfRed("%d rows 'missing' submatrix\n", nr);
		printTable("'missing' submatrix", missingMatrix, nr, nc, m_groupSize);
		printTableColor("'missing submatrix' links", m_pGraph, nc, nc, 0);
	}
#endif
	if (np == 1)
		printTableColor("\nInput links", m_pGraph, nc, nc, 0);
	/*leo*/
	tchar tmp[16 * 16];
	memcpy(tmp, m_pGraph, nc * nc);
	auto* pCanonGraph = m_pGraphCanonizer->canonize_graph();
	int ic = memcmp(tmp, pCanonGraph, nc * nc);
	if (ic > 0) {
		printTableColor("\nInput links", tmp, nc, nc, 0);
		printTableColor("Canonized", pCanonGraph, nc, nc, 0);
		exit(1);
	}
	return ic < 0 ? 0 : 1;
	if (np == 1)
		printTableColor("Canonical links", pCanonGraph, nc, nc, 0);
#if CheckMissingMatrix
	if (memcmp(canonicalLinks, pCanonGraph, nc * nc) == 0 && 
		memcmp(missingMatrix, r, nr * nc) != 0) {
		printfYellow("Submatrix below has the same canonical Links table as the 'missing' submatrix\n");
		printTable("Submatrix with same as 'missing' submatrix canonical Links table", r, nr, nc, m_groupSize);
		printTable("'Missing' submatrix", missingMatrix, nr, nc, m_groupSize);
		printTableColor("Canonical Links", pCanonGraph, nc, nc, 0);
	}
#endif
	int m_lenMask = 16; // needed for assert in SetMask below
	for (tchar i = 0; i < nc; i++) {
		for (tchar j = i + 1; j < nc; j++) {
			if (pCanonGraph[i * nc + j]) {
				SetMask(msk, i, j);
			}
		}
	}
	std::lock_guard<std::mutex> lock(mtxLinks);
	int id = nr - 1;
	if (mpLinks == NULL) {
		NumMatricesProcessed = 1;
		mpLinks = new CStorageIdx<tchar>*[m_numDaysResult];
		memset(mpLinks, 0, m_numDaysResult * sizeof(mpLinks[0]));
	}
	if (!mpLinks[id]) {
		mpLinks[id] = new CStorageIdx<tchar>(50000000, 16);
	}
	bool bSame = false;
	NumMatricesProcessed++;
	if (mpLinks[id]->isProcessed((ctchar*)msk)) {
	//if (mpLinks[id]->findObject((ctchar*)msk, 0, mpLinks[id]->numObjects()) != UINT_MAX) {
#if CheckMissingMatrix
		if (memcmp(missingMatrix, r, nr * nc) == 0) {
			printfYellow("Canonical links table for 'missing' submatrix is already present in DB, Skipped\n");
		}
#endif
		bSame = true;
		SemiPhase++;
	}

	if (!(NumMatricesProcessed % np)) {
		printfGreen(" %d Total checked(%d are the same). Row:Links ", NumMatricesProcessed, SemiPhase);
		for (int i = 3; i < m_numDaysResult; i++)
			if (mpLinks[i])
				printfGreen("%d:%d ", i + 1, mpLinks[i]->numObjects());
		printf("\n");
	}
	if (bSame)
		return 0;
	return id;
}
CC ePrecalculateReturn alldata::precalculatedSolutions(eThreadStartMode iCalcMode)
{
	ll msk[2];
	int iDaySaved = 0;
	int iLastDay = 0;
ProcessPrecalculatedRow:
	if (iDaySaved) {
		if (m_pRowUsage->getMatrix2(result(), neighbors(), numDaysResult(), iDaySaved)) {
#if !USE_CUDA
			m_cTime = clock();
			for (int i = m_nPrecalcRows + 3; i < numDaysResult(); i++)
				m_rowTime[i] = m_cTime - m_iTime;
#endif
			iDay = numDaysResult() - 1;
			return eCheckCurrentMatrix;
		}
		iDay = iDaySaved;
		iDaySaved = 0;
		bPrevResult = false;
		m_playerIndex = 0;
	}
	if (bPrevResult) {
		iDay--;
		bPrevResult = false;
	}
	if (iDay < m_nPrecalcRows)
		m_secondPlayerInRow4First = 0;
	else {
		int ipx = 0;
		if (m_playerIndex)
		{
			iDay = m_playerIndex / m_numPlayers;
			ipx = m_playerIndex % m_numPlayers;
			m_playerIndex = 0;
			if (iDay < m_nPrecalcRows)
			{
				if (param(t_MultiThreading))
				{
					return eNoResult;
				}
				m_pRowStorage->init();
				m_secondPlayerInRow4 = m_secondPlayerInRow4First = 0;
				m_playerIndex = 0;
				m_nRows4 = m_nRows4Day = 0;
				return eContinue;
			}
			if (iDay >= numDaysResult())
			{
				//ASSERT_IF(1);
				return eNoResult;
			}
		}
		/**leo
		if (iLastDay >= iDay) {
			std::lock_guard<std::mutex> lock(mtxLinks);
			mpLinks[iLastDay]->updateRepo((ctchar*)msk); 
			iLastDay = 0;
		}*/
		if (m_lastRowWithTestedTrs >= iDay)
			m_lastRowWithTestedTrs = iDay - 1;
		const auto retVal = m_pRowUsage->getRow(iDay, ipx);
		switch (retVal) {
		case 2:
			iDaySaved = iDay;
			break;
		case 1:
			///m_pRowUsage->getMatrix(result(), neighbors(), iDay + 1);
			//printTable("tbl", result(), iDay + 1, m_numPlayers, groupSize());

			m_playerIndex = 0;
#if !USE_CUDA
			if (m_bPrint && iDay < m_nPrecalcRows + 3) {
				m_cTime = clock();
				m_rowTime[iDay] = m_cTime - m_iTime;
			}
#endif
			if (++iDay < numDaysResult() && !checkSubmatrix()) {

#if 0   // Temporary
				if (!p1f_counter || ((++m_p1f_counter) % p1f_counter))
#endif
#if !USE_CUDA
					if (!m_bPrint || m_cTime - m_rTime < ReportInterval)
#endif
						goto ProcessPrecalculatedRow;
			}

			m_pRowUsage->getMatrix(result(), neighbors(), iDay);
#if CheckMissingMatrix
			//if (memcmp(missingMatrix, result(), iDay * m_numPlayers) == 0) 
			//	iDay = iDay;
#endif
#if 1
			if ((m_test & 32) && checkSubmatrix() && m_groupSize == 2 && m_numPlayers <= 16) {
				/*leo
				if (iLastDay) {
					printfRed("iLastDay not 0(%d)\n", iLastDay);
					exit(1);
				}*/
				iLastDay = addMaskToDB(msk, result(0), iDay, m_numPlayers, 100000);
				if (!iLastDay)
					iDay--;
				goto ProcessPrecalculatedRow;
			}	
#endif
			iDay--;
#if !USE_CUDA
			m_cTime = clock();
			for (int i = m_nPrecalcRows + 3; i <= iDay; i++)
				m_rowTime[i] = m_cTime - m_iTime;
#endif
			return eCheckCurrentMatrix;
		case -1: // reported if requested row not in solution
		{
			// continue to case 0:
		}
		case 0:
#if !USE_CUDA
			if (iDay == 10 && (m_test & 64) && m_bPrint) {
				m_pRowUsage->getMatrix(result(), neighbors(), iDay);
				linksFromMatrix(links(), result(), iDay);
				printTable("Precalc matrix", result(), iDay, numPlayers(), m_groupSize);
				printTableColor("Links for matrix above", links(0), numPlayers(), numPlayers(), m_groupSize);
				static int cc;
				if (++cc > 5)
					exit(1);
			}
#endif
			iDay--;
			break;
		}
		goto ProcessPrecalculatedRow;
	}
	if (iCalcMode == eCalcResult) {
		m_precalcMode = eCalculateRows;
		m_pRowStorage->init();
		iDay++;
		bPrevResult = true;
		return eContinue;
	}
	return eNoResult;
}
