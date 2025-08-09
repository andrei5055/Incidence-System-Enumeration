#include <iostream>
#include "TripleSys.h"

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
			iDay--;
			break;
		}
		goto ProcessPrecalculatedRow;
	}
	if (iCalcMode == eCalcResult) {
		m_useRowsPrecalculation = eCalculateRows;
		m_pRowStorage->init();
		iDay++;
		bPrevResult = true;
		return eContinue;
	}
	return eNoResult;
}