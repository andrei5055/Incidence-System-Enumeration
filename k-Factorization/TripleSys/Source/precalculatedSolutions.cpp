#include <iostream>
#include "TripleSys.h"

#include <mutex>
extern std::mutex mtxLinks; // The mutex to protect the shared resource
extern CStorageIdx<tchar>** mpLinks;
extern int SemiPhase;

CC ePrecalculateReturn alldata::precalculatedSolutions(eThreadStartMode iCalcMode)
{
	int iDaySaved = 0;
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
				//ASSERT(1);
				return eNoResult;
			}
		}

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
			if (++iDay < numDaysResult() && !checkCanonicity()) {

#if 0   // Temporary
				if (!p1f_counter || ((++m_p1f_counter) % p1f_counter))
#endif
#if !USE_CUDA
					if (!m_bPrint || m_cTime - m_rTime < ReportInterval)
#endif
						goto ProcessPrecalculatedRow;
			}

			m_pRowUsage->getMatrix(result(), neighbors(), iDay);
#if 1
			if (iDay == m_numDaysResult && (m_test & 32) && m_groupSize == 2 && m_numPlayers == 16) {
				ll msk[2];
				msk[0] = msk[1] = 0;
				tchar* r = result();
				int m_lenMask = 16;
				for (tchar i = 0; i < (m_numDaysResult * 16); i += 2) {
					SetMask(msk, r[i], r[i + 1]);
				}
				std::lock_guard<std::mutex> lock(mtxLinks);
				if (mpLinks == NULL) {
					mpLinks = new CStorageIdx<tchar>*[1];
					mpLinks[0] = new CStorageIdx<tchar>(5000000, 16);
				}
				if (mpLinks[0]->isProcessed((ctchar*)msk))
					SemiPhase++;
				else {
					if (!(mpLinks[0]->numObjects() % 100000000))
						printfGreen(" nMasks(%d rows)=%d, same=%d\n", iDay, mpLinks[0]->numObjects(), SemiPhase);
				}
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
#if 0 //!USE_CUDA
			if (iDay == 10 && m_bPrint) {
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