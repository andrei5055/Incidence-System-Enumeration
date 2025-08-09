#include <iostream>
#include "TripleSys.h"

CC ePrecalculateReturn alldata::endOfRowPrecalculation(eThreadStartMode iCalcMode)
{
	if (m_nPrecalcRows > iDay) {
		m_useRowsPrecalculation = eCalculateRows;
		m_pRowStorage->init();
	}
	else if (m_nPrecalcRows == iDay) {
		ASSERT(m_secondPlayerInRow4 == 0);
		const bool bCBMP = !completeGraph();
		m_secondPlayerInRow4++;
		if (bCBMP && !m_groupSizeRemainder[m_secondPlayerInRow4])
			m_secondPlayerInRow4++;
		if (!m_nRows4Day && m_groupSize == 2)
		{
			m_secondPlayerInRow4 = m_secondPlayerInRow4First;
			bPrevResult = true;
			if (m_nRows4) {
				m_nRows4 = 0;
				if (iCalcMode == eCalculateRows) {
					nLoops = m_nRows4;
					return eNoResult;
				}
				m_pRowUsage->init();
				m_pRowStorage->reset();
			}
			return eContinue;
		}
		m_nRows4Day = 0;
		if (m_secondPlayerInRow4 <= m_secondPlayerInRow4Last)
			return eContinue;
		m_secondPlayerInRow4 = m_secondPlayerInRow4First = 0;
		if (m_nRows4) {
			iDay = m_nPrecalcRows;
			if (m_bPrint) {
				printf("Total number of precalculated row solutions = %5d\n", m_nRows4);
				m_lastRowWithTestedTrs = 0;
			}
			m_useRowsPrecalculation = eCalculateMatrices;
			m_playerIndex = 0;

			if (!m_pRowStorage->initCompatibilityMasks()) {
#if !USE_CUDA
				//printfRed("*** Unexpected error returned by initCompatibilityMask()\n");
				//ASSERT(1);
#endif
				if (param(t_MultiThreading)) {
					nLoops = 0;
					return eNoResult;
				}
				iDay = m_nPrecalcRows;
				m_useRowsPrecalculation = eCalculateRows;
				m_pRowStorage->init();
				m_nRows4 = 0;	iDay++;
				//linksFromMatrix(links(), result(), iDay);
				//iPlayer = m_numPlayers;
				bPrevResult = true;
				return eContinue;
			}

			if (iCalcMode == eCalculateRows) {
				nLoops = m_nRows4;
				return eNoResult;
			}
			m_pRowUsage->init();
			//						m_pRowStorage->reset();
			m_nRows4 = 0;
			return eContinue;
		}
	}
	return eOk;
}
CC void alldata::addPrecalculatedRow()
{
	if (!m_nRows4++ && iDay >= 2) {
		TrCycles trc;
		cyclesFor2Rows(m_TrCyclesFirst2Rows, &trc, neighbors(0), neighbors(1), result(0), result(1));
	}

	m_nRows4Day++;
#if 0
	printf("%6d:", m_nRows4);
	printTable("", neighbors(3), 1, m_numPlayers, 2);
	printTable("r", result(3), 1, m_numPlayers, 2);
#endif
#if 0
	bool bP1F = p1fCheck3(result(0), result(m_nPrecalcRows), neighbors(0), neighbors(m_nPrecalcRows));
	if (!bP1F)
		printf("not p1f\n");
	bP1F = p1fCheck3(result(2), result(m_nPrecalcRows), neighbors(2), neighbors(m_nPrecalcRows));
	if (!bP1F)
		printf("not p1f\n");
#endif
	ASSERT(!m_secondPlayerInRow4First);
#if 0
	if (iTest) {
		printTable("r3", result(m_nPrecalcRows), 1, m_numPlayers, m_groupSize);
		if (result(m_nPrecalcRows)[1] == 8)
			iTest = iTest;
	}
#endif
	tchar nb[MAX_PLAYER_NUMBER];
	u1fSetTableRow(nb, result(m_nPrecalcRows), true);
	const auto retVal = m_pRowStorage->addRow(result(m_nPrecalcRows), neighbors(m_nPrecalcRows), nb);
	if (!retVal) {
		m_playerIndex = m_nPrecalcRows * m_numPlayers;
		goBack();
		m_pRowStorage->init();
		m_secondPlayerInRow4 = m_secondPlayerInRow4First = 0;
		m_playerIndex = 0;
		m_nRows4 = m_nRows4Day = 0;
	}
}
CC void alldata::initPrecalculationData(eThreadStartMode iCalcMode, int nRowsStart)
{
	const bool bCBMP = !completeGraph();
	m_nPrecalcRows = param(t_useRowsPrecalculation);
	if (iCalcMode == eCalcSecondRow) {
		iCalcMode = eCalcResult;
		m_numDaysResult = 2;
		nRowsStart = m_nPrecalcRows = 0;
	}
	if (iCalcMode == eCalcResult)
		m_useRowsPrecalculation = (m_nPrecalcRows && m_nPrecalcRows <= 3 && m_groupSize <= 3 && nRowsStart <= m_nPrecalcRows) ? eCalculateRows : eDisabled;
	else
		m_useRowsPrecalculation = iCalcMode;

	tchar secondPlayerMax = m_numPlayers - (m_groupSize == 2 ? 1 : m_groupSize + 1 + m_numDays - m_numDaysResult);
	tchar m_secondPlayerInRow4First = 0;
	if (bCBMP) {
		secondPlayerMax = m_numPlayers - 1 - (m_groupSize - 2) * m_groupSize;
	}
	m_secondPlayerInRow4Last = MIN2(param(t_lastRowSecondPlayer), secondPlayerMax);
	if (!m_secondPlayerInRow4Last)
		m_secondPlayerInRow4Last = m_groupSize == 2 ? (bCBMP ? m_numDaysResult * 2 - 1 : m_numDaysResult) : secondPlayerMax;
	m_secondPlayerInRow4 = 0;

	m_nRows4 = 0;
	m_nRows4Day = 0;
#if 1 // preset automorphism groups
	const auto iCalc = m_useRowsPrecalculation;
	m_useRowsPrecalculation = eCalcResult;
	if (!param(t_autGroupNumb)) {
		if (param(t_useAutForPrecRows) > 1 && iDay >= param(t_useAutForPrecRows)) {
			auto* pGroupInfo = groupInfo(param(t_useAutForPrecRows));
			if (pGroupInfo) {
				cnvCheckNew(0, param(t_useAutForPrecRows), false); // create initial set of tr for first i rows
				pGroupInfo->copyIndex(*this);
				resetGroupOrder();
			}
		}
	}
	else {
		for (int i = firstGroupIdx(); i <= lastGroupIdx(); i++) {
			auto* pGroupInfo = groupInfo(i);
			if (!pGroupInfo)
				break;
			if (iDay < i)
				break;
			cnvCheckNew(0, i, false); // create initial set of tr for first i rows
			pGroupInfo->copyIndex(*this);
			resetGroupOrder();
		}
	}
	m_useRowsPrecalculation = iCalc;
#endif
}