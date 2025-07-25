#include "TripleSys.h"
#include "CheckCanon.h"

typedef enum {
	t_reasonUnknown,
	t_ordering,
	t_invertOrdering,
	t_changing_day_0,
	t_changing_day_0_group,
	t_playerPosition_1_4,
	t_NotThatPlayerInPosition_1_4,
	t_Statement_7,
	t_Statement_18,
	t_Statement_19,
	t_Statement_19_G,
} t_RejectionRreason;

static const char* reason[] = {
		"Reason unknown",
		"Ordering problem",
		"Inverted matrix has smaller code",
		"Improving code by changing day 0",
		"Improving code by changing day 0 and using group",
		"Only players 4 or 9 can be in position [1,4]",
		"Player #%d cannot be at position [1, 4]",
		"Player# in [1, 4] should be less than [1, 7]",
		"Incorrect order of players in a group with a leading player #%d",
		"Rejected by generalization of Statement 19 for the group #%d of day 1",
		"Rejected by generalization of Statement 19 for the group of degree #%d of day 1",
};

template class CCheckerCanon<SIZE_TYPE>;

template<typename T>
void renumberPlayers(T* pntr, size_t i, size_t iLast) {
	for (; i < iLast; i++)
		pntr[i] = pntr[pntr[i]];
}

template<typename T>
void groupOrdering(T* pElems, size_t numGroup, T* buffer, size_t groupSize, T* pDays = NULL) {
	// adj - adjustment of pointers used when comparing values.
	//    0 for comparing groups within a day 
	//    1 for comparing days
	T adj = pDays ? 1 : 0;
	const auto len = groupSize * sizeof(*pElems);
	const auto iMax = numGroup * groupSize;
	for (size_t i = 0; i < iMax; i += groupSize) {
		auto bestIdx = i;
		auto bestVal = *(pElems + i + adj);
		for (size_t j = i + groupSize; j < iMax; j += groupSize) {
			const auto curVal = *(pElems + j + adj);
			if (bestVal > curVal) {
				bestVal = curVal;
				bestIdx = j;
			}
		}

		if (bestIdx != i) {
			memcpy(buffer, pElems + i, len);
			memcpy(pElems + i, pElems + bestIdx, len);
			memcpy(pElems + bestIdx, buffer, len);
			if (pDays) {
				// Rearranging day indices
				const auto i1 = i / groupSize;
				const auto i2 = bestIdx / groupSize;
				const auto tmp = pDays[i1];
				pDays[i1] = pDays[i2];
				pDays[i2] = tmp;
			}
		}
	}
}

CheckerCanon(void)::sortTuples(T *players) const {
	elemOrdering(players, m_numElem, groupSize());
	groupOrdering(players, numGroups(), tmpBuffer(), groupSize());
}

#if DEBUG_NextPermut
int perm_cntr, matr_cntr;
bool flg = false;
#define M_CNTR 2172 //167 //496
#endif

CheckerCanon(bool)::CheckCanonicity(const T* result, int nDays, int* pGrpNumb, T* bResult) {
	// Input parameters:
	//    result - pointer to a sequence of lists, each containing "m_numElem" players
	//             for each day, players are divided into groups, each of which contains "n"m_groupSise" players
	//    nDays  - number of days (or lists, mentioned above)
	// Output parameter:
	//    pGrpNumb - pointer to the group number of the last row which needed to be change 
	//               if the return value of that method is false;  
	//    bResult (optional) - pointer to the array of (m_numElem + nDays) elements,
	//             when it's not NULL and the "result" is not a canonical one, then
	//             this array will contain the permutations of
	//             (a) players (as its first "m_numElem" elements)
	//             (b) days (starting with (bResult+m_numElem)'s element)
	//    

	setNumReason(t_RejectionRreason::t_reasonUnknown);
	*pGrpNumb = numGroups() * nDays - 2;
	setResultOut(bResult);
	setStudiedMatrix(result, nDays);
	setTrivialPerm(result);
	resetImprovedResultFlag();

	m_pDestMemory = bResult ? bResult : resultMemory();
	if (bResult)
		createDaySequence();

	if (!checkCanonicity())
		return false;

	if (!checkWithGroup(numGroups(), &CCheckerCanon<T>::checkReorderedGroups, result)) {
		*pGrpNumb = groupIndex();
		return false;
	}

	if (m_numDays == m_numDaysMax) {
#if DEBUG_NextPermut
		if (++matr_cntr == M_CNTR)
			matr_cntr += 0;
		perm_cntr = 0;
#endif
#if CHECK_WITH_GROUP
		T ttr[24];
		const auto retVal = checkWithGroup(m_numElem, &CCheckerCanon<T>::orderingMatrix, ttr, false);
		if (!retVal)
			*pGrpNumb = groupIndex();
#if DEBUG_NextPermut
		if (false && 187 >= matr_cntr && matr_cntr >= 162) {
			FOPEN_F(f, "../CCC.txt", "a");
			fprintf(f, "matr_cntr = %3d: retVal = %d  *pGrpNumb = %d\n", matr_cntr, retVal, *pGrpNumb);
			FCLOSE_F(f);
		}
#endif
		return retVal;
#endif	
	}

	return true;
}

CheckerCanon(int)::checkReorderedGroups(const T* permut, T nElem, const T* pMatr) {
	// Construct a permutation of players corresponding to the permutations of the groups of the first day. 
	auto pPlayers = playersPerm(1);
	const auto len = groupSize() * sizeof(*pPlayers);
	auto i = nElem;
	while (i--) {
		memcpy(pPlayers + i * groupSize(), pMatr + permut[i] * groupSize(), len);
	}

	auto pntrTo = destMemory();
	auto pntrFrom = pMatr;
	for (T iDay = 1; iDay < numDays(); iDay++)
		recodePlayers(pPlayers, pntrFrom += numElem(), pntrTo += numElem());

	// Copying trivial permutation of days
	memcpy(pntrTo + numElem(), pMatr, numDays() * sizeof(T));
	orderigRemainingDays(1, 0, destMemory());
	const int diff = memcmp(destMemory() + numElem(), pMatr + numElem(), lenRow() * (numDays() - 1));
	if (diff < 0) {
		auto pRes = destMemory();
		auto pInput = pMatr;

		int dayIdx = 1;
		auto len = lenRow();
		while (!memcmp(pRes += numElem(), pInput += numElem(), len))
			dayIdx++;

		// If the set of days was reordered, we need to adjust grpIdx 
		const auto pDayIdx = destMemory() + numElem() * numDays();
		for (T i = 1; i <= dayIdx; i++) {
			if (dayIdx < pDayIdx[i])
				dayIdx = pDayIdx[i];
		}

		T grpIdx = 0;
		len /= numGroups();
		while (!memcmp(pRes, pInput, len)) {
			pRes += groupSize();
			pInput += groupSize();
			grpIdx++;
		}

		for (T i = 0; i < grpIdx; i++) {
			if (grpIdx < permut[i])
				grpIdx = permut[i];
		}

		setGroupIndex(dayIdx * numGroups() + grpIdx);
		if (resultOut()) {
			memcpy(destMemory(), pMatr, lenRow());
			addImproveResultFlags(t_bResultFlags::t_readyToExplainMatr);
			const auto len = commentBufferLength();
			auto pBuffer = comment();
			SPRINTFS(pBuffer, comment(), len, "Reordering of players:\n");
			for (T i = 0; i < nElem; i++) {
				if (i == permut[i])
					continue;

				int idx = permut[i];
				SPRINTFS(pBuffer, comment(), len, " (");
				for (int k = 0; k < 2; k++) {
					idx *= groupSize();
					for (T j = 0; j < groupSize(); j++)
						SPRINTFS(pBuffer, comment(), len, "%2d,", idx + j);

					pBuffer--;
					SPRINTFS(pBuffer, comment(), len, k ? ")" : ") ==> (");
					idx = i;
				}
			}

			addImproveResultFlags(t_bResultFlags::t_readyToExplainTxt);
		}
	}

	return diff;
}

CheckerCanon(bool)::checkCanonicity() {
	T leadingPlayers[5];
	for (T iDay = 0; iDay < numDays(); iDay++) {
		if (!checkOrderingForDay(iDay)) {
			if (!resultOut())
				return false;

			// Make a trivial permutations for set of days from preconstructed (1, 0, 2, 3, 4, ...)
			auto pDest = destMemory() + lenResult();
			*pDest = 0;
			*++pDest = 1;
			reportTxtError(resultOut(), reason[t_RejectionRreason::t_ordering], NULL, iDay);
			return checkRemainingDays(iDay);
		}

		// Check canonicity of the codes for the other days
		if (iDay) {
			// Do this only when day 0 changed its place.
			setDayNumb(iDay);
			auto pRow = getMatrixRow(iDay);
			T maxVal;
			memcpy(leadingPlayers, trivialPerm(), (maxVal = groupSize()) * sizeof(leadingPlayers[0]));
			auto pPlayerPerm = playersPerm(3);
			memcpy(pPlayerPerm, pRow, lenRow());

			// Loop over different groups of the first day
			while (true) {
				sortTuples(playersPerm());
				if (!checkDay_1(iDay, pPlayerPerm))
					return false;

				if (maxVal == numElem())
					break;

				// Starting working with the next group
				maxVal = initNextSetOfGroups(maxVal, pRow, pPlayerPerm, leadingPlayers);
			}
		}
		else {
			// Check all remaining days for canonicity.
			char buffer[256];
			setReasonParam(-1);
			for (int j = 1; j < numDays(); j++) {
				if (checkDay(j))
					continue;

				if (!resultOut())
					return false;

				auto pReason = reason[numReason()];
				if (reasonParam() != -1) {
					sprintf_s(buffer, pReason, reasonParam());
					pReason = buffer;
				}

				const auto dayToBlame = numReason() == t_RejectionRreason::t_invertOrdering ? numDays() - 1 : j;
				return reportTxtError(resultOut(), pReason, NULL, dayToBlame);
			}
		}
	}

	return true;
}

CheckerCanon(T)::switchLeadingPlayersOfGroups(T placeIdx, T * playerPerm, const T* pLeaders) const {
	static T leaderPlace[] = { 0, 2, 1, 1, 0, 2, 2, 1, 0, 1, 2, 0, 2, 0, 1 }; // for now S(3) hardcoded
	auto pFrom = leaderPlace + placeIdx;
	for (T i = 0; i < groupSize(); i++)
		playerPerm[i * groupSize()] = pLeaders[pFrom[i]];

	recordTuples(playerPerm, playersPerm());
	return (placeIdx += 3) < countof(leaderPlace) ? placeIdx : ELEMENT_MAX;
}


CheckerCanon(T)::initNextSetOfGroups(T maxVal, const T* pRow, T *playerPerm, T *pLeaders) const {
	const auto minVal = maxVal;
	maxVal += groupSize();
	auto pRowLast = pRow + numElem();
	// Find groups with the elements in the interval [minVal, maxVal)
	auto tmpTo = playerPerm;
	auto tmpTo1 = playerPerm + groupSize() * groupSize();
	auto idx = 0;
	while (true) {
		T elem = 0;
		auto* pTo = &tmpTo1;
		auto k = 0;
		for (; k < groupSize(); k++) {
			elem = pRow[k];
			if (elem >= maxVal)
				break;
			if (elem >= minVal) {
				pTo = &tmpTo;
				break;
			}
		}

		// Copying this group to the beginning OR to the end of the array
		memcpy(*pTo, pRow, groupSize() * sizeof(*tmpTo));
		if (*pTo == tmpTo) {
			// Group with one of the elements we are looking for is found
			// Let's put it first
			if (k) {
				tmpTo[k] = tmpTo[0];
				tmpTo[0] = elem;
			}
							
			*pLeaders++ = elem; // Saving leading player

			if (++idx >= groupSize()) {
				// Copying remaining groups
				pRow += groupSize();
				memcpy(tmpTo1, pRow, (pRowLast - pRow) * sizeof(*tmpTo1));
				break;
			}
		}
		*pTo += groupSize();
		pRow += groupSize();
	}

	recordTuples(playerPerm, playersPerm());
	return maxVal;
}

CheckerCanon(bool)::checkWithGroup(T numElem, int (CCheckerCanon<T>::*func)(const T*, T, const T*), const T* pCurrentRow, bool symmetrical) {
	T* permut = permutation();
	T lenStab = 0;

	size_t nAuto = 0;
	// Copying trivial permutation
	const auto len = numElem * sizeof(T);
	memcpy(permut, trivialPerm(), len);
	memcpy(orbits(), trivialPerm(), len);

	const auto calcGroupOrder = numDays() == m_numDaysMax && numElem == this->numElem() || !symmetrical;
	const auto rowPermut = calcGroupOrder;

	T (CCheckerCanon<T>:: *nextPerm)(T*, const T*, T, T, T) = symmetrical? &CCheckerCanon<T>::next_permutation : &CCheckerCanon<T>::nextPermutationA;
	if (!symmetrical) {
		memset(m_pPermIndex, 0, 2 * numGroups() * sizeof(m_pPermIndex[0]));
		memcpy(m_pGroupPerm, trivialPerm(), numGroups() * sizeof(m_pGroupPerm[0]));
		if (USE_TRANSLATE_BY_LEO)
			m_pAD->initCheckByGroup(numDays(), 0);
		else
			memcpy(m_dayIdx, trivialPerm(), (m_nDaysToTest = numDays()) * sizeof(m_dayIdx[0]));

		if (!m_pSubGroup) {
			m_pSubGroup = new T[m_GroupOrder * groupSize()];
			const auto len = groupSize() * sizeof(m_pSubGroup[0]);
			memcpy(m_pSubGroup, trivialPerm(), len);
			auto pElemNext = m_pSubGroup;
			for (auto i = 1; i < m_GroupOrder; i++) {
				auto pElemPrev = pElemNext;
				memcpy(pElemNext += groupSize(), pElemPrev, len);
				CGroupOrder<T>::setStabilizerLength(ELEMENT_MAX);
				CGroupOrder<T>::next_permutation(pElemNext, NULL, groupSize());
			}
		}
	}
	
	auto pPerm = (symmetrical || !USE_EQUAL) ? permut : pCurrentRow;

	CGroupOrder<T>::setStabilizerLength(numElem - 1);
	CGroupOrder<T>::setStabilizerLengthAut(ELEMENT_MAX);
	CGroupOrder<T>::setGroupOrder(1);
	T nElem = ELEMENT_MAX;
#if DEBUG_NextPermut
	if (flg = (matr_cntr == M_CNTR))
		matr_cntr += 0;
#endif
#define PRINT_PERMUT  0
#define PRINT_PERMUT_ 0
#if PRINT_PERMUT || PRINT_PERMUT_
	char buffer[256], *ptr;
	static int ctr_canon = 0;
	size_t ctr = 0;
	size_t counter = 1;	
	const auto fName = "../auto_orb.txt";
	if (calcGroupOrder && matr_cntr == 1) {
		FOPEN_F(f, fName, "w");
		FCLOSE_F(f);
	}
#endif
#if PRINT_PERMUT_
	size_t counter = 1;

	FOPEN_F(f, fName, "w");
#endif
	bool autFlag = false;
	T first = 0;
	if (!symmetrical)
		goto check_func;  // Trivial permutation need to be checked

	while (true) {
		nElem = (this->*nextPerm)(permut, orbits(), numElem, nElem, lenStab);
		if (nElem == ELEMENT_MAX)
			break;

#if PRINT_PERMUT_
		ptr = buffer;
		SPRINTFD(ptr, buffer, "%5zd:", ++counter);
		for (T i = 0; i < numElem; i++)
			SPRINTFD(ptr, buffer, " %3d", permut[i]);

		_printf(f, false, "%s\n", buffer);
#endif
#if NEED_TO_DEBUG
		if (!symmetrical) {
			const auto frst = permut[0];
			if (first < frst) {
				if (frst != pOrb[frst]) {
					ggg++;
					nElem = IDX_MAX;
					continue;
				}
			}
		}
#endif
	check_func:
		const auto diff = (this->*func)(permut, numElem, pCurrentRow);
		if (diff < 0)
			return false;

		if (!diff) {
			// Automorphism found
			nElem = IDX_MAX;
			if (!symmetrical)
				continue;

			autFlag = true;
			CGroupOrder<T>::UpdateOrbits(pPerm, numElem, orbits(), rowPermut, calcGroupOrder);
#if PRINT_PERMUT
			if (calcGroupOrder) {
				FOPEN_F(f, fName, "a");
				ptr = buffer;
				SPRINTFD(ptr, buffer, "  permut = %4zd     :", ++counter);
				for (T i = 0; i < numElem; i++)
					SPRINTFD(ptr, buffer, " %3d", permut[i]);

				_printf(f, false, "%s\n", buffer);

				ptr = buffer;
				SPRINTFD(ptr, buffer, "matr_cntr = %4d     :", matr_cntr);
				for (T i = 0; i < numElem; i++)
					SPRINTFD(ptr, buffer, " %3d", pPerm[i]);

				_printf(f, false, "%s\n", buffer);

				ptr = buffer;
				SPRINTFD(ptr, buffer, "perm_cntr = %5d (%1zd):", perm_cntr, ++ctr);
				for (T i = 0; i < numElem; i++)
					SPRINTFD(ptr, buffer, " %3d", orbits()[i]);

				SPRINTFD(ptr, buffer, ":  |Aut(G)| = %zd\n", CGroupOrder<T>::groupOrder());
				_printf(f, false, buffer);
				FCLOSE_F(f);
			}
#endif
		} else
			nElem = numElem;
	}

#if 0
	if (!symmetrical) {
		FOPEN_F(f, "../auto_star.txt", matr_cntr != 1 ? "a" : "w");
		if (f) {
			fprintf(f, "%6d: nAuto = %zd\n", matr_cntr, nAuto);
			FCLOSE_F(f);
		}
	}
#endif

	if (autFlag && rowPermut && calcGroupOrder)
		CGroupOrder<T>::updateGroupOrder(numElem, orbits());

#if PRINT_PERMUT_
	FCLOSE_F(f);
#endif
#if PRINT_PERMUT
	if (calcGroupOrder) {
		FOPEN_F(f, fName, "a");
		sprintf_s(buffer, "%4d: |Aut(G)| = %zd\n", ++ctr_canon, CGroupOrder<T>::groupOrder());
		_printf(f, false, buffer);
		FCLOSE_F(f);
	}
#endif
	return true;
}

CheckerCanon(bool)::copyTuple(const T* res, T inc, bool doCheck) const {
	T prev;
	auto players = playersPerm();
	players[prev = res[0]] = inc;
	for (T j = 1; j < groupSize(); j++) {
		const auto next = res[j];
		ASSERT(prev == next);
		if (doCheck && prev > next)
			return false;

		players[prev = next] = j + inc;
	}
	return true;
}

CheckerCanon(bool)::rollBack(T* p_dayRes, T* p_dayIsUsed, int& j, int nDays) const
{
	while (j > 1) { // Do we need to go to previous day?
		p_dayIsUsed[p_dayRes[--j]] = 0;
		if (++p_dayRes[j] < nDays) {
			j--;
			return true;
		}

		p_dayRes[j] = 0;
	}

	return false;
}

CheckerCanon(bool)::checkOrderingForDay(T nDay) const {
	auto res = getMatrixRow(nDay);
	if (res[0] || !copyTuple(res))
		return false;

	T inc = 0;
	for (auto j = numGroups(); --j;) {
		if (!copyTuple(res += groupSize(), inc += groupSize()) ||
			m_players[*res] < m_players[*(res - groupSize())]) // Comparing first elements of the groups
			return false;
	}

	return true;
}

CheckerCanon(bool)::explainRejection(const T* players, T playerPrevID, T playerNewID, T firstDayID, bool doOutput, const T* pNewOrder) {
	auto* pDest = resultOut();
	if (!pDest)
		return false;              // we don't need an explanation for rejection

	// Create a conversion table.
	const T* result = studiedMatrix();
	memcpy(pDest, result, lenRow());  
	if (pNewOrder) {
		ASSERT(result == pNewOrder);
		T i = 1;
		while (playerPrevID < playerNewID) {
			pDest[pNewOrder[i]] = playerPrevID++;
			i += 2;
		}
	}
	else
		pDest[pDest[playerPrevID] = playerNewID] = playerPrevID;

	// Make a conversion with just created table.
	const auto numElem_2 = 2 * m_numElem;
	auto pDest_1 = pDest + m_numElem;
	memcpy(pDest_1, players, lenRow());
	renumberPlayers(pDest, m_numElem, numElem_2);

	groupOrdering(pDest_1, numGroups(), tmpBuffer(), groupSize());
	const auto diff = memcmp(pDest_1, result + m_numElem, lenRow());
	ASSERT(diff >= 0);

	if (result + m_numElem == players || doOutput || numDays() == 2) {
		addImproveResultFlags(t_bResultFlags::t_readyToExplainMatr);
		if (!doOutput) // When doOutput is true, the correct last rows have already been copied.
			memcpy(pDest + numElem_2, result + numElem_2, (lenResult() - numElem_2) * sizeof(*pDest));

		renumberPlayers(pDest, numElem_2, lenResult());

		orderigRemainingDays(2, 0, pDest);
		memcpy(pDest, result, lenRow());
		*(pDest + 1) = 1 - (*(pDest += lenResult()) = firstDayID);
	}

	return false;
}

CheckerCanon(bool)::checkPermutationOfFirstDayGroups(int numGroup, const T* pCurrentRow, bool useCurrentRow)
{
	if (useCurrentRow)
		memcpy(playersPerm(1), pCurrentRow, lenRow());

	if (checkWithGroup(numGroup, &CCheckerCanon<T>::checkPermutationOnGroups, pCurrentRow))
		return true;

	if (resultOut()) {
		setNumReason(t_RejectionRreason::t_Statement_19_G);
		setReasonParam(numGroup);	// the number of the groups used
	}
	return false;
}

CheckerCanon(int)::checkPermutationOnGroups(const T* permGroup, T numElem, const T* pCurrentRow)
{
	auto pTmp = playersPerm(1);
	const auto lenGroup = groupSize() * sizeof(*pTmp);
	for (T i = 0; i < numElem; i++)
		memcpy(pTmp + i * groupSize(), pCurrentRow + permGroup[i] * groupSize(), lenGroup);

	recordTuples(pTmp, m_tmpBuffer1);
	memcpy(pTmp = playersPerm(2), m_tmpBuffer1, lenRow());
	sortTuples(m_tmpBuffer1);
	const auto retVal = checkDayCode(0, 1, m_tmpBuffer1, false);
	if (retVal > 0)
		return retVal;

	if (retVal < 0 && !resultOut())
		return retVal;

	const auto rc = checkRemainingDays(dayNumb(), 0, NULL, pTmp);
	if (retVal < 0 || !rc)
		return -1;

	return 0;
}

CheckerCanon(bool)::checkPosition1_4(const T *players) {
	// Statement 7: In canonical matrix z1 < z2
	//    0  1  2    3  4  5    6  7  8 ....
	//    0  3  6    1 z1  *    2 z2 *
#if USE_STATEMENT_7
	if (groupSize() == 3 && players[4] > players[7]) {
		setNumReason(t_RejectionRreason::t_Statement_7);
		return explainRejection(players, 1, 2);
	}
#endif
#if USE_STATEMENT_18
	// Statement 18: For any group size s > 1, any n > 0 and any pairs (0<=i < j< s) in canonical matrix the position
	// of player P(s, n, i) in the second day is less than the position of player P(s, n, j) in the same day.
	//
	// Create a pReorder table of with groupSize() following pairs: (player_position, player_ID)
	T idx;
	auto pNewOrder = resultMemory();
	pNewOrder[idx = 0] = 0;  // Player number 0 is always on position [1,0]
	bool flag = false;
	T i = 1;
	for (; i < numElem(); i++) {
		T j = 0;
		while (players[++j] != i);

		if (i % groupSize()) {
			if (j < pNewOrder[idx])
				flag = true; // problem found
			idx += 2;
		}
		else {
			if (flag)
				break;		// problem on prervious group was found
			idx = 0;
		}

		pNewOrder[idx] = j;
		pNewOrder[idx + 1] = i;
	}

	if (flag) {
		groupOrdering(pNewOrder, groupSize(), tmpBuffer(), 2);
		setNumReason(t_RejectionRreason::t_Statement_18);
		const auto playerPrevID = i - groupSize();
		setReasonParam(playerPrevID);
		return explainRejection(players, playerPrevID, i, 0, false, pNewOrder);
	}
#endif

#if	(UsePos_1_4_condition & 2)
	if (players != studiedMatrix() + m_numElem)
		return true;
#endif
#if USE_STATEMENT_17
	// Statement 17: Only the players 4 or 9 could be in position[1, 4].
	// 
#if USE_STATEMENT_18 == 0
	// When USE_STATEMENT_18 != 0, these cases are covered
	// List of simple player substitutions (subst[i] <---> subst[i+1], for i%2 == 0) 
	// which will improve the matrix code, if player subst[i] is at position [1, 2]
	static T subst[] = { 8, 7,10, 9,11, 9 };
	for (int i = 0; i < sizeof(subst) / sizeof(subst[0]); i += 2) {
		const auto playerID = subst[i];
		if (players[4] == playerID) {  // player is on position #4 of day #1
			setNumReason(t_RejectionRreason::t_NotThatPlayerInPosition_1_4);
			setReasonParam(subst[i]);
			return explainRejection(players, playerID, subst[i + 1]);
		}
	}
#endif
	if (players[4] == 7) {
		if (!resultOut())
			return false;              // we don't need an explanation for rejection

		setNumReason(t_RejectionRreason::t_NotThatPlayerInPosition_1_4);
		setReasonParam(7);

		if (players == studiedMatrix() + m_numElem) {
			// Do this only when day 0 did not changed its place.
			// As described in Statement 17, swap days 0 and 1
			recordTuples(players, playersPerm());
			checkRemainingDays(1, -1);
			return explainRejection(playersPerm(), 1, 2, 1, true);
		}
		return false;
	}
#endif
#if USE_STATEMENT_19
	// Swaping group #0 of day 1 with all other groups
	auto pntr = getMatrixRow(1);
	const auto lenGroup = groupSize() * sizeof(T);
#if 1   // Both of these options work and result in similar rejection of matrices.
	auto pTmp = playersPerm(1);
	auto pntrFrom = pntr;
	auto pntrTo = pTmp;
	// NOTE: There is no point in replacing group #0 of the first day with the group # > groupSize() 
	// and leaving at least one of the groups #1, #2, ...groupSize()-1 on their places. 
	// In this case, we will not get the right leading group (which is (0, 3, 6) for triples) on the first day.
	auto pSecondRow = playersPerm(2);
	for (T i = 1; i < groupSize(); i++) {
		memcpy(pTmp, pntr, lenRow());
		memcpy(pTmp, pntrFrom += groupSize(), lenGroup);
		memcpy(pntrTo += groupSize(), pntr, lenGroup);
		recordTuples(pTmp, pSecondRow);
		sortTuples(pSecondRow);

		const int diff = checkDayCode(0, 1, pSecondRow, false);
		if (diff < 0) {
			if (!resultOut())
				return false;       // we don't need an explanation for rejection

			setNumReason(t_RejectionRreason::t_Statement_19);
			setReasonParam(i);		// the number of the group used
			return checkRemainingDays(1, diff, pTmp);
		}
	}
	return true;
#else
	// For some reason, using a symmetrical group acting on players #0 - #2
	// does not add any new rejections for numPlayers 15 or 21 cases.
	// But just in case, we will keep the following fragment...
	return checkPermutationOfFirstDayGroups(groupSize(), getMatrixRow(1), false);
#endif
#endif
}

CheckerCanon(void)::createDaySequence(T iDay) const {
	// Adding all unused days to the array.
	auto pDest = destMemory() + lenResult();
	*pDest++ = iDay;
	auto j = *pDest = 0;
	while (++j < numDays()) {
		if (j != iDay)
			*++pDest = j;
	}
}

CheckerCanon(int)::checkDayCode(int diff, T iDay, const T *secontRow, bool createDaySeq) {
	if (!diff)
		diff = memcmp(secontRow, studiedMatrix() + m_numElem, lenRow());

	if (!diff && numDays() > 2 || diff < 0 && resultOut()) {
		T* pDest = destMemory();
		if (diff != -9999) {
			// Saving two first days:
			memcpy(pDest, studiedMatrix(), lenRow());
			memcpy(pDest + m_numElem, secontRow, lenRow());
		}

		if (diff < 0) {
			if (diff == -9999) {
				// Adding 2 day indices
				*(pDest += lenResult()) = iDay;
				*++pDest = 0;
				reportTxtError(pDest - 1, reason[numReason()], pDest - 1, 2);
			}
			else
				addImproveResultFlags(t_bResultFlags::t_readyToExplainMatr);
		}
	}

	// Adding all unused days to the array
	if (createDaySeq && resultOut())
		createDaySequence(iDay);

	return diff;
}

CheckerCanon(bool)::checkRemainingDays(T iDay, int retVal, const T *pPerm, T *pPermPlayer) {
	// A function to record the remaining days and check them when the code for 
	// the first two days was built exactly the same as in the input matrix.
	const auto numElem_2 = 2 * m_numElem;
	auto* pDest = destMemory();
	if (iDay) {
		// Renumbering of players according to permutation of days: (0, iDay).
		if (!pPermPlayer) {
			if (!pPerm)
				pPerm = getMatrixRow(iDay);

			pPermPlayer = playersPerm(2);
			recordTuples(pPerm, pPermPlayer);
		}

		// Renumbering the set of players in the groups for all days except 0 and iDay.
		auto* pOut = pDest + 2 * m_numElem;
		auto* pIn = studiedMatrix();
		for (int j = 1; j < numDays(); j++) {
			pIn += m_numElem;
			if (j == iDay)
				continue;   // Skip day which is already used as 0's one. 

			recodePlayers(pPermPlayer, pIn, pOut);
			pOut += m_numElem;
		}

		if (retVal == -9999) {
			renumberPlayers(pDest, numElem_2, lenResult() - numElem_2);
			memcpy(pDest, studiedMatrix(), lenRow());
		}
	}

	orderigRemainingDays(2, 0, pDest);

	// Comparing all remaining days
	if (retVal == 0) {
		const auto lenCmp = (lenResult() - m_numElem2) * sizeof(T);
		retVal = memcmp(pDest + numElem_2, getMatrixRow(2), lenCmp);
	}

	if (retVal >= 0)
		return true;

	addImproveResultFlags(t_bResultFlags::t_readyToExplainMatr);
	return false;
}

CheckerCanon(bool)::checkDay_1(T iDay, const T* pPlayerPerm) {
	// iDay index if the day which replaced the day 0
	int diff = 0;
#if	(UsePos_1_4_condition & 2)
	const auto players = playersPerm();
#if UsePos_1_4_condition && ImproveResults
	if (!checkPosition1_4(players)){
		diff = -9999;
		if (numDays() == 2)
			return diff;
	}
#else
	if (players[4] != 4 && players[4] != 9) {
		diff = -9999;
		return diff;
	}
	if (players[4] > players[7])
		return -1;
#endif
#endif
	diff = checkDayCode(diff, iDay, playersPerm());
	if (diff < 0 && !resultOut())
		return false;

	if (diff <= 0 && !checkRemainingDays(iDay, diff, pPlayerPerm))
		return reportTxtError(resultOut(), reason[t_changing_day_0], NULL, iDay);

#if USE_CHANGING_DAY_0_GROUPS
	if (numDays() == m_numDaysMax && !checkPermutationOfFirstDayGroups(numGroups(), pPlayerPerm)) {
		if (!resultOut())
			return false;

		char* pAddExpl;
		reportTxtError(resultOut(), reason[t_changing_day_0_group], NULL, iDay, &pAddExpl);

		auto pBuffer = comment();
		const auto len = commentBufferLength();
		const auto* pPlayerPerm = playersPerm(1);
		SPRINTFS(pAddExpl, pBuffer, len, "\nPlayer permutation: (%d", *pPlayerPerm);

		for (T i = 1; i < numElem(); i++)
			SPRINTFS(pAddExpl, pBuffer, len, ", %d", *(pPlayerPerm + i));

		SPRINTFS(pAddExpl, pBuffer, len, ").");
		return false;
	}
#endif

	return true;
}

CheckerCanon(int)::orderingMatrix(T nDays, T numGroups, bool expected, bool invert, const T* permPlayer) {
	T* pDest = resultOut();
	if (!pDest)
		pDest = resultMemory();

	memcpy(pDest, studiedMatrix(), lenResult() * sizeof(*pDest));
	if (invert) {
		// Transformation: player #i new number is: numElem() - 1 - i
		auto pTmp = playersPerm(1);
		const auto j = numElem() - 1;
		for (auto i = numElem(); i--;)
			pTmp[i] = j - i;

		permPlayer = pTmp;
	}

	if (permPlayer)
		recodePlayers(permPlayer, pDest, pDest, lenResult());
	
	// Copying trivial permutations for days 
	memcpy(pDest + lenResult(), trivialPerm(), numDays() * sizeof(*pDest));
	orderigRemainingDays(nDays, numGroups, pDest);
	const auto diff = memcmp(pDest, studiedMatrix(), lenResult() * sizeof(*pDest));
	if (diff < 0) {
		setNumReason(invert ? t_RejectionRreason::t_invertOrdering : t_RejectionRreason::t_ordering);
		addImproveResultFlags(t_bResultFlags::t_readyToExplainMatr);
		return diff;
	}

	ASSERT((!expected) && (diff < 0));
	return diff;
}

CheckerCanon(bool)::checkDay(T iDay) {
	auto pMatrixRow = getMatrixRow(iDay);
	T* pDest = resultOut();
#if	(UsePos_1_4_condition & 1)
	if (iDay == 1) {
		if (!preordered()) {
			ASSERT(pDest==NULL);  // When matrices are not preordered, 
								  // we expect to have an external buffer
			setPreordered(true);
			if (orderingMatrix(0, 0, false) < 0)
				return false;
		}

		if (!checkPosition1_4(pMatrixRow))
			return false;

		if (orderingMatrix(0, 0, false, true) < 0)
			return false;
	}
#endif

	auto res(pMatrixRow);
	auto resEnd = res;
	T j = 0;
	for (; j < numGroups(); j++) {
		resEnd += groupSize();
		while (++res < resEnd && *(res - 1) < *res);

		if (res < resEnd)
			break;

		// Check ordering of the first elements of the groups
		if (j && *(res - groupSize()) < *(res - 2 * groupSize()))
			break;
	}

	if (j == numGroups())
		return true;

	if (!pDest)
		return false;

	return orderingMatrix(iDay, j) < 0;
}

CheckerCanon(void)::orderigRemainingDays(T daysOK, T groupsOK, T *pDest) const {
	// Function will reorder all elements in the groups except 
	//     - all groups of the first daysOK days
	//     - first groupsOK of the (daysOK+1)-th day 
	// and after that will reorder
	//     - the new groups inside coressponding day
	//     - all days in accordance with the second element of the first new group of the day. 

	const auto lenOK = daysOK * m_numElem;
	const auto nElemOK = lenOK + groupsOK * groupSize();
	// Ordering all remaining groups of the elements
	elemOrdering(pDest + nElemOK, lenResult() - nElemOK, groupSize());

	// Ordering groups by their first elements
	T i = daysOK;
	auto pntr = pDest + lenOK;
	while (i++ < numDays()) {
		groupOrdering(pntr, numGroups(), tmpBuffer(), groupSize());
		pntr += m_numElem;
	}

	const auto daysToOrder = numDays() - daysOK;
	if (daysToOrder > 1) {
		// Ordering days by the first two elements of their first group
		groupOrdering(pDest + lenOK, daysToOrder, tmpBuffer(), m_numElem, pDest + lenResult() + daysOK);
	}
}

CheckerCanon(bool)::reportTxtError(T *bBuffer, const char *pReason, const T *pDays, T nDay, char **ppAddExpl) {
	// For pDays 1= NULL, nDay is the number of days defined in pDays
	// Otherwise nDay is a day numbed OR indefined (-1)
	if (!bBuffer)
		return false;

	// Caller expects to get explanations or better result
	addImproveResultFlags(t_bResultFlags::t_readyToExplainTxt);
	const auto len = commentBufferLength();
	const auto copyLen = strlen(pReason);
	auto pBuffer = comment();
	memcpy_s(pBuffer, len, pReason, copyLen);
	if (!pDays && nDay == -1) // Do we need to add some day related information to the comment? 
		return false;         // No, we don't   

	// Adding day related information 
	auto* pBuff = pBuffer + copyLen;

	SPRINTFS(pBuff, pBuffer, len, ". The problem was detected for day");
	if (pDays && nDay > 1) {
		SPRINTFS(pBuff, pBuffer, len, "s: (%d", pDays[0]);
		for (int i = 1; i < nDay; i++)
			SPRINTFS(pBuff, pBuffer, len, ", %d", pDays[i]);
			
		SPRINTFS(pBuff, pBuffer, len, ").");
	} else
		SPRINTFS(pBuff, pBuffer, len, " %d.", nDay);

	if (ppAddExpl)
		*ppAddExpl = pBuff;

	return false;
}

void outInfo(FILE* f, const unsigned char* pInfo, int len, const char* pName) {
	char buf[256], * p = buf;
	SPRINTFD(p, buf, "%8s:", pName);
	for (int i = 0; i < len; i++)
		SPRINTFD(p, buf, " %d", pInfo[i]);
	fprintf(f, "%s\n", buf);
}

#if DEBUG_NextPermut
CheckerCanon(void)::printInfo(FILE* f, const T* perm, int idx) const {
	fprintf(f, "Returned from %d (%d)\n", idx, perm_cntr);
	outInfo(f, perm, 10, "Perm");
	outInfo(f, m_pGroupPerm, 5, "GrPerm");
	outInfo(f, m_pPermIndex, 5, "IdxPerm");
	outInfo(f, m_pNumPermutUsed, 5, "NPermUsd");
	fclose(f);
}
#endif

CheckerCanon(int)::orderingMatrix(const T* permut, T numElem, const T* pDummy) {
#if 0
	extern int matr_cntr;
	extern int cntr;
	static char file_name[32];
	if (cntr == 1)
		sprintf_s(file_name, "../permuts_%03d.txt", matr_cntr);

	char buffer[256], * pBuf = buffer;
	for (T j = 0; j < numElem; j++)
		SPRINTFD(pBuf, buffer, "%3d", permut[j]);

	FOPEN_F(f, file_name, cntr != 1 ? "a" : "w");
	const auto retVal = orderingMatrix(0, 0, false, false, permut);
	if (f) {
		fprintf(f, "%6d: %s  retVal = %d\n", cntr, buffer, retVal);
		FCLOSE_F(f);
	}

	if (retVal < 0 || cntr == 3839)
		cntr += 0;

	return retVal;
#else
#if USE_TRANSLATE_BY_LEO
	int dayMax = 0;  // Andrei: Previously we defined dayMax in cnvCheckKm1
	// For some reason, we don't do this. 
	// We can't use USE_TRANSLATE_BY_LEO = 1 until we fix this. 
	const auto ret = m_pAD->cnvCheckKm1(permut, numDays(), (T*)pDummy);
	setGroupIndex((dayMax + 1) * numGroups() - 2);
	return ret;
#else

	int ret = 1;
	auto* ttr = (T*)(pDummy);
	for (T n = 0; n < m_nDaysToTest; n++) {
		auto* resn = getMatrixRow(m_dayIdx[n]);
		for (T i = 0; i < numElem; i++)
			ttr[resn[i]] = permut[i];

		const auto retVal = orderingMatrix(0, 0, false, false, ttr);

		if (retVal > 0)
			continue;

		if (!retVal) {
			ret = 0;
			if (n)
				m_dayIdx[n--] = m_dayIdx[--m_nDaysToTest];

			continue;
		}

		static int cntr;
		T* pDest = resultOut();
		if (!pDest)
			pDest = resultMemory();

		const auto* pDays = pDest + lenResult();
		const auto* pInputMatr = studiedMatrix();
		T dayMax = *pDays;
		do {
			pDest += numElem;
			pInputMatr += numElem;
			if (dayMax < *++pDays)
				dayMax = *pDays;

		} while (!memcmp(pDest, pInputMatr, numElem));
		setGroupIndex((dayMax + 1) * numGroups() - 2);
#if PRINT_TRANSFORMED
		extern bool flg;
		if (flg)
			printTransformed(m_nDaysToTest, numElem, groupSize(), (const tchar*)permut, (const tchar*)ttr, (const tchar*)studiedMatrix(), (const char*)resultMemory(), n);
#endif
		return -1;
	}
	return ret;
#endif
#endif
}

CheckerCanon(T)::nextPermutationA(T* perm, const T* pOrbits, T nElem, T idx, T lenStab) {
#if DEBUG_NextPermut
	++perm_cntr;
	static int cntrMax = 3840;
#if USE_ORBTS == 0
	static int ctrIdx[] = ////, 173, 125, 317, 195
#if M_CNTR == 9
	{ 16, 13, 45, 29, 173, 125, 317, 195 };
#elif M_CNTR == 11
	{ 16, 13, 45, 28, 76, 49, 337, 97, 1249, 865};
#elif M_CNTR == 1955
	{1536, 1370, 2906, 1754};
#else
	{3840};
#endif
	static int jStart;
	if (matr_cntr == M_CNTR && jStart <= countof(ctrIdx) - 2 && perm_cntr >= ctrIdx[jStart]) {
		cntrMax -= ctrIdx[jStart] - (perm_cntr = ctrIdx[jStart + 1]);
		jStart += 2;
	}
#endif

	FILE* f = NULL;
	flg = false;
	if (false && matr_cntr == M_CNTR) {
		flg = true;
		if (/*perm_cntr == 2471 ||*/ perm_cntr == 1369 /*2310*/)
			flg = true;

		FOPEN_F(ff, "../ddd.txt", "a");
		f = ff;
	}
#endif
	const auto len = groupSize() * sizeof(T);
	auto iMinVar = numGroups() - 1;
	if (USE_ORBTS && idx == IDX_MAX) {
		// Permutation perm defines an automorphism.
		T j = 0;
		while (perm[j] == j) 
			j++;

		// We are preparing to reset all remaining groups.
		auto i = iMinVar = j / groupSize() + 1;
		memset(m_pNumPermutUsed + i, 0, (numGroups() - i) * sizeof(m_pNumPermutUsed[0]));

		// Utilizing either a decreasing or increasing ordering of m_pGroupPerm:
		const auto less = true;// m_pPermIndex[i - 1] < m_GroupOrder - 1;
		for (; i < numGroups(); i++) {
			// To do that we need for all remaining groups
			// (a) to set m_pPermIndex to their highest values
			m_pPermIndex[i] = m_GroupOrder;
			// (b) reorder all group indices  
			auto idx = m_pGroupPerm[i];
			for (auto k = i; ++k < numGroups();) {
				if ((idx < m_pGroupPerm[k]) == less) {
					m_pGroupPerm[i] = m_pGroupPerm[k];
					m_pGroupPerm[k] = idx;
					idx = m_pGroupPerm[i];
				}
			}
		}
	}

	bool orderTail = false;
	auto pPerm = perm + nElem;
	for (auto i = numGroups(); i--;) {
		pPerm -= groupSize();
		auto val = m_pGroupPerm[i] * groupSize();
	start_loop:
		auto const* pSubgrPerm = m_pSubGroup;
		if (++m_pPermIndex[i] < m_GroupOrder)
			pSubgrPerm += m_pPermIndex[i] * groupSize();
		else
			m_pPermIndex[i] = 0;

		for (T j = 0; j < groupSize(); j++)
			pPerm[j] = val + pSubgrPerm[j];

		if (m_pPermIndex[i]) {
			// Skip stabilizer elements
			if (USE_ORBTS && (!i || !memcmp(perm, trivialPerm(), i * len))) {
				auto k = i * groupSize();
				while (perm[k] == k) 
					k++;

			check_orb:
				if (perm[k] != pOrbits[perm[k]]) {
					if (m_pPermIndex[i] < groupSize() - 1 || idx == IDX_MAX) {
						idx = 0;
						goto start_loop;
					}

					if (orderTail) {
						// The tail of the permutation needs to be re-established.
						const auto j = i * groupSize();
						memcpy(m_pGroupPerm + i, trivialPerm() + i, (numGroups() - i) * sizeof(T));
						memcpy(perm + j , trivialPerm() + j, (nElem - j) * sizeof(T));
						m_pPermIndex[i] = 0;
						continue;
					}

					// Try to find a group to swap with the current one.
					auto j = i;
					auto pGr = pPerm;
					while (++j < numGroups() && val > *(pGr += groupSize()));

					const auto tmp = m_pGroupPerm[i];
					if (j < numGroups()) {
						// Swapping groups and their indices
						memcpy(pPerm, pGr, len);
						// Set the initial order of the elements 
						// of the current group on its new location.
						for (T j = 0; j < groupSize(); j++)
							pGr[j] = val + j;

						val = (m_pGroupPerm[i] = m_pGroupPerm[j]) * groupSize();
						m_pGroupPerm[j] = tmp;
						m_pPermIndex[i] = 0;
						m_pNumPermutUsed[i + 1] = 0;
						//memset(m_pNumPermutUsed + i, 0, (numGroups() - i) * sizeof(m_pNumPermutUsed[0]));
						goto check_orb;
					}
					
					// Impossible to swap
					if (!i)
						break;  // No more permutations 

					// Moving all groups and element of permutation to the left
					const auto nElem = numGroups() - i - 1;
					memcpy(pPerm, pPerm + groupSize(), nElem * len);
					memcpy(m_pGroupPerm + i, m_pGroupPerm + i + 1, nElem * sizeof(m_pGroupPerm[0]));
					m_pGroupPerm[nElem + i] = tmp;
					m_pPermIndex[i] = 0;
					m_pNumPermutUsed[i + 1] = 0;
					m_pNumPermutUsed[i] = 1;
					//memset(m_pNumPermutUsed + i, 0, (numGroups() - i) * sizeof(m_pNumPermutUsed[0]));
					// and the current group to last position.
					for (T j = 0; j < groupSize(); j++)
						pGr[j] = val + j;

					continue;
				}
			}
#if DEBUG_NextPermut
			if (f)
				printInfo(f, perm, 1);
#endif
			return 0;
		}

		if (i >= iMinVar)
			continue;

		if (m_pGroupPerm[i] < m_pGroupPerm[i + 1]) {
			CGroupOrder<T>::setStabilizerLength(ELEMENT_MAX);
			CGroupOrder<T>::next_permutation(m_pGroupPerm, NULL, numGroups());
			auto j = i;
			while (j < numGroups()) {
				memcpy(perm + j * groupSize(), trivialPerm() + m_pGroupPerm[j] * groupSize(), len);
				j++;
			}

			if (USE_ORBTS && (!i || !memcmp(perm, trivialPerm(), i * len))) {
				const auto j = i * groupSize();
				if (perm[j] != pOrbits[perm[j]]) {
					++m_pNumPermutUsed[i];
					orderTail = true;
					val = m_pGroupPerm[i] * groupSize();
					goto start_loop;
				}
			}

#if DEBUG_NextPermut
			if (f)
				printInfo(f, perm, 2);
#endif
			return 0;
		}

		if (++m_pNumPermutUsed[i] < m_GroupOrder) {
			auto j = numGroups();
			auto i1 = i - 1;
			revert(m_pGroupPerm, j, (T)i1);
			auto shift = i * groupSize();
			memcpy(tmpBuffer() + shift, perm + shift, (j - i) * len);
			while (++i1 < --j) {
				memcpy(perm + shift, tmpBuffer() + j * groupSize(), len);
				memcpy(perm + j * groupSize(), tmpBuffer() + shift, len);
				shift += groupSize();
			}
		} else
			m_pNumPermutUsed[i] = 0;
	}

#if DEBUG_NextPermut
	FCLOSE_F(f);
#endif
	return ELEMENT_MAX;
}

CheckerCanon(bool)::CheckPermutations(const T* result, const T* pMatrix, int nRows, int numFactors) {
	T* permPlayers = playersPerm();
	memcpy(permPlayers, result, lenRow());
	auto pOrbits = new T [numElem()];
	memcpy(permPlayers, result, lenRow());
	memcpy(pOrbits, result, lenRow());
	bool checkPermut = false;
	int errLine, errGroup, dubLine;
	size_t cntr = 1;
	size_t counter[2] = { 0 };
	const auto *pDayPerm = pMatrix + nRows * numElem();
	T nRow = 0;
	T lenStab = 0;
	T idx = ELEMENT_MAX;
	tchar lnks[MAX_PLAYER_NUMBER * MAX_PLAYER_NUMBER];
	printf(" I am going to the loop:\n");

	while (true) {
		if (_CheckMatrix((const tchar *)pMatrix, nRows, numElem(), groupSize(), lnks, true, &errLine, &errGroup, &dubLine, numFactors)) {
			counter[0]++;
			char buffer[256];
			auto pntr = buffer;
			for (T j = 0; j < numElem(); j++)
				pntr += sprintf_s(pntr, 256 - (pntr - buffer), "%3d", permPlayers[j]);

			pntr += sprintf_s(pntr, 256 - (pntr - buffer), "\n");
			printf(buffer);
		}
		else
			counter[1]++;

		const auto nElem = CGroupOrder<T>::next_permutation(permPlayers, pOrbits, numElem(), idx, lenStab);
		if (nElem == ELEMENT_MAX)
			break;

		cntr++;
		break;
		auto pntrTo = (T *)pMatrix + numElem();
		for (T i = 2; i < nRows; i++) {
			memcpy(pntrTo += numElem(), result + numElem() * pDayPerm[i], lenRow());
			recodePlayers(permPlayers, pntrTo, pntrTo);
		}
	}

	printf("Good and Bad counters are %zd - %zd\n", counter[0], counter[1]);
	delete[] pOrbits;
	return counter[1] == 0;
}

