
#include <assert.h>
#include "CheckCanon.h"

typedef enum {
	t_reasonUnknown,
	t_ordering,
	t_invertOrdering,
	t_changing_day_0,
	t_playerPosition_1_4,
	t_NotThatPlayerInPosition_1_4,
	t_Statement_7,
	t_Statement_18,
	t_Statement_19,
} t_RejectionRreason;

static const char* reason[] = {
		"Reason unknown",
		"Ordering problem",
		"Inverted matrix has smaller code",
		"Increasing code by changing day 0",
		"Only players 4 or 9 can be in position [1,4]",
		"Player #%d cannot be at position [1, 4]",
		"Player# in [1, 4] should be less than [1, 7]",
		"Incorrect order of players in a group with a leading player #%d",
		"Rejected by generalization of Statement 19 for the group #%d of day 1",
};

template class CCheckerCanon<SIZE_TYPE, SIZE_TYPE>;

template<typename T>
void renumberPlayers(T* pntr, size_t i, size_t iLast) {
	for (; i < iLast; i++)
		pntr[i] = pntr[pntr[i]];
}

template<typename T>
void elemOrdering(T* pElems, size_t numElem, size_t groupSize) {
	// Ordering elements in the groups od size groupSize
	auto j = numElem + groupSize;
	switch (groupSize) {
	case 2:
		// Ordering groups of pairs.
		for (; j -= 2; pElems += 2) {
			if (pElems[0] > pElems[1]) {
				const auto tmp = pElems[0];
				pElems[0] = pElems[1];
				pElems[1] = tmp;
			}
		}
		return;
	case 3:
		// Ordering groups of triples.
		for (; j -= 3; pElems += 3) {
			const auto tmp0 = pElems[0];
			const auto tmp1 = pElems[1];
			const auto tmp2 = pElems[2];
			if (tmp2 > tmp1) {
				if (tmp0 > tmp1) {
					pElems[0] = tmp1;
					if (tmp2 < tmp0) {
						pElems[1] = tmp2;
						pElems[2] = tmp0;
					}
					else
						pElems[1] = tmp0;
				}
			}
			else {
				if (tmp2 > tmp0) {
					pElems[1] = tmp2;
					pElems[2] = tmp1;
				}
				else {
					pElems[0] = tmp2;
					if (tmp0 < tmp1) {
						pElems[1] = tmp0;
						pElems[2] = tmp1;
					}
					else
						pElems[2] = tmp0;
				}
			}
		}
		return;
	}

	assert(false); // Not implemented for given groupSize
}

template<typename T>
void groupOrdering(T* pElems, size_t numGroup, T* buffer, size_t groupSize = GroupSize, T* pDays = NULL) {
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

CheckerCanon(void)::sortTuples() const {
	elemOrdering(m_players, m_numElem, groupSize());
	groupOrdering(m_players, numGroups(), tmpBuffer(), groupSize());
}

CheckerCanon(bool)::CheckCanonicity(const T *result, int nDays, T *bResult) {
	// Input parameters:
	//    result - pointer to a sequence of lists, each containing "m_numElem" players
	//             for each day, players are divided into groups, each of which contains "n"m_groupSise" players
	//    nDays  - number of days (or lists, mentioned above)
	// Output parameter:
	//    bResult (optional) - pointer to the array of (m_numElem + nDays) elements,
	//             when it's not NULL and the "result" is not a canonical one, then
	//             this array will contain the permutations of
	//             (a) players (as its first "m_numElem" elements)
	//             (b) days (starting with (bResult+m_numElem)'s element)
	//    

	T numReason = t_RejectionRreason::t_reasonUnknown;

	setResultOut(bResult);
	setStudiedMatrix(result, nDays);
	resetImprovedResultFlag();

	m_pDestMemory = bResult ? bResult : resultMemory();
	if (bResult)
		createDaySequence();

	for (int iDay = 0; iDay < nDays; iDay++) {
		if (!checkOrderingForDay(iDay)) {
			if (!bResult)
				return false;

			// Make a trivial permutations for set of days from preconstructed (1, 0, 2, 3, 4, ...)
			auto pDest = destMemory() + lenResult();
			*pDest = 0;
			*++pDest = 1;
			reportTxtError(bResult, reason[t_RejectionRreason::t_ordering], NULL, iDay);
			return checkRemainingDays(iDay);
		}

		// Check canonicity of the codes for the other days
		if (iDay) {
			// Do this only when day 0 changed its place.
			sortTuples();
			const auto retVal = checkDay_1(iDay, &numReason);
			if (nDays == 2)
				return retVal >= 0; // In this case there is nothing more to do than has already been done.

			if (retVal > 0)
				continue;

			if (retVal < 0 && !bResult)
				return false;  // The result has improved, but we don't need to know how.

			if (!checkRemainingDays(iDay, retVal)) 
				return reportTxtError(bResult, reason[t_RejectionRreason::t_changing_day_0], NULL, iDay);
		}
		else {
			// Check all remaining days for canonicity.
			T playerNumb = -1;

			for (int j = 1; j < nDays; j++) {
				if (!checkDay(j, &numReason, &playerNumb)) {
					if (!bResult)
						return false;

					char buffer[256];
					auto pReason = reason[numReason];
					if (playerNumb != -1) {
						sprintf_s(buffer, pReason, playerNumb);
						pReason = buffer;
					}

					const auto dayToBlame = numReason == t_RejectionRreason::t_invertOrdering? nDays - 1 : j;
					return reportTxtError(bResult, pReason, NULL, dayToBlame);
				}
			}
		}
	}

#if 0
	T* pOrbits = m_players + m_numElem;

		int j = 0;
		T k;
		while (true) {
			while (++j < nDays) {
				// Looking for the first unused day
				k = p_dayRes[j];
				while (k < nDays && p_dayIsUsed[k])
					k++;

				if (k == nDays) // Unable to find unused day.
					break;

				p_dayIsUsed[p_dayRes[j] = k] = 1;
				if (!checkDay_1(p_players, k, j, p_dayRes, p_dayIsUsed)) {
					return false;
				} else {
					break;
				}
			}

			if (!rollBack(p_dayRes, p_dayIsUsed, j, nDays))
				break;
		}

		p_dayIsUsed[iDay] = 0;
		continue;   // temporary
#endif
#if 0
		if (m_numDays != m_numDaysMax)
			return true;

		static int ccc = 0;
		if (!ccc++)
			return true;

		return checkWithGroup(result, m_numElem);
#else
	return true;
#endif
}

CheckerCanon(bool)::checkWithGroup(const T* result, T numElem/*, int (*func)(CCheckerCanon<T>, const T*)*/) {
	T* permut = permutation();
	T lenStab = 0;

	CGroupOrder<T>::setStabilizerLength(numElem - 1);
	CGroupOrder<T>::setStabilizerLengthAut(ELEMENT_MAX);

	// Copying trivial permutation
	const auto len = numElem * sizeof(T);
	memcpy(permut, result, len);
	memcpy(oprbits(), result, len);

	const auto calcGroupOrder = false;
	const auto rowPermut = false;
	size_t counter = 1;
	size_t ctr = 1;
	T nElem = numElem;
	T idx = ELEMENT_MAX;
	while (true) {
		nElem = nextPermutation(permut, oprbits(), numElem, idx, lenStab);
		if (nElem == ELEMENT_MAX)
			break;

		counter++;
#if 0
		char buffer[256], * ptr = buffer;
		SPRINTFD(ptr, buffer, "%5zd:", counter);
		for (T i = 0; i < numElem; i++)
			SPRINTFD(ptr, buffer, " %3d", permut[i]);

		printf("%s\n", buffer);
#endif
		const auto diff = orderingMatrix(0, 0, NULL, false, false, permut);
		if (diff < 0)
			return false;

		if (!diff) {
			// Automorphism found
			ctr++;
			CGroupOrder<T>::UpdateOrbits(permut, numElem, oprbits(), rowPermut, calcGroupOrder);
		}
	}

	return true;
}

CheckerCanon(bool)::copyTuple(const T* res, T inc, bool doCheck) const {
	T prev;
	m_players[prev = res[0]] = inc;
	for (T j = 1; j < groupSize(); j++) {
		const auto next = res[j];
		assert(prev != next);
		if (doCheck && prev > next)
			return false;

		m_players[prev = next] = j + inc;
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
	auto res = studiedMatrix() + nDay * m_numElem;
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
		return false;              // we don't need explanation for rejection

	// Create a conversion table.
	const T* result = studiedMatrix();
	memcpy(pDest, result, lenRow());  
	if (pNewOrder) {
		assert(result != pNewOrder);
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
	assert(diff < 0);

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

#if 0
CheckerCanon(bool)::checkPermutationOfFirstDaygroups()
{
	T permPlayers[3], orbits[3];
	memcpy(permPlayers, getMatrixRow(0), lenGroup);
	memcpy(orbits, permPlayers, lenGroup);
	memcpy(pTmp, pntr, lenRow());
	while (true) {
		const auto nElem = nextPermutation(permPlayers, orbits, groupSize(), idx, 0);
		if (nElem == ELEMENT_MAX)
			break;

		for (T i = 0; i < 3; i++)
			memcpy(pTmp + i * groupSize(), pntr + permPlayers[i] * groupSize(), lenGroup);

		recordTuples(pTmp);
		sortTuples();

		const int diff = checkDayCode(0, NULL, 1);
		if (diff < 0) {
			if (!resultOut())
				return false;       // we don't need explanation for rejection

			*pNumReason = t_RejectionRreason::t_Statement_19;
			*pNumPlayer = i;		// the number of the group used
			memcpy(m_players, pTmp, lenRow());
			return checkRemainingDays(1, diff, m_players);
		}
	}
}
#endif

CheckerCanon(bool)::checkPosition1_4(const T *players, T *pNumReason, T* pNumPlayer) {
	// Statement 7: In canonical matrix z1 < z2
	//    0  1  2    3  4  5    6  7  8 ....
	//    0  3  6    1 z1  *    2 z2 *
#if USE_STATEMENT_7
	if (players[4] > players[7]) {
		*pNumReason = t_RejectionRreason::t_Statement_7;
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
		*pNumReason = t_RejectionRreason::t_Statement_18;
		const auto playerPrevID = i - groupSize();
		if (pNumPlayer)
			*pNumPlayer = playerPrevID;

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
	for (int i = 0; i < countof(subst); i += 2) {
		const auto playerID = subst[i];
		if (players[4] == playerID) {  // player is on position #4 of day #1
			*pNumReason = t_RejectionRreason::t_NotThatPlayerInPosition_1_4
				* pNumPlayer = subst[i];
			return explainRejection(players, playerID, subst[i + 1]);
		}
	}
#endif
	if (players[4] == 7) {
		if (!resultOut())
			return false;              // we don't need explanation for rejection

		*pNumReason = t_RejectionRreason::t_NotThatPlayerInPosition_1_4;
		*pNumPlayer = 7;

		if (players == studiedMatrix() + m_numElem) {
			// Do this only when day 0 did not changed its place.
			// As described in Statement 17, swap
			checkOrderingForDay(1); // days 0 and 1
			checkRemainingDays(1, -1);
			return explainRejection(m_players, 1, 2, 1, true);
		}
		return false;
	}
#endif
#if USE_STATEMENT_19
	// Swaping group #0 of day 1 with all other groups
	auto pntr = getMatrixRow(1);
	const auto lenGroup = groupSize() * sizeof(T);
	auto pTmp = m_players + m_numElem;
#if 0
	auto pntrFrom = pntr;
	auto pntrTo = pTmp;
	// NOTE: There is no point in replacing group #0 of the first day with the group # > groupSize() 
	// and leaving at least one of the groups #1, #2, ...groupSize()-1 on their places. 
	// In this case, we will not get the right leading group (which is (0, 3, 6) for triples) on the first day.
	for (T i = 1; i < groupSize(); i++) {
		memcpy(pTmp, pntr, lenRow());
		memcpy(pTmp, pntrFrom += groupSize(), lenGroup);
		memcpy(pntrTo += groupSize(), pntr, lenGroup);
		recordTuples(pTmp);
		sortTuples();

		const int diff = checkDayCode(0, NULL, 1);
		if (diff < 0) {
			if (!resultOut())
				return false;       // we don't need explanation for rejection

			*pNumReason = t_RejectionRreason::t_Statement_19;
			*pNumPlayer = i;		// the number of the group used
			memcpy(m_players, pTmp, lenRow());
			return checkRemainingDays(1, diff, m_players);
		}
	}
#else
	// For some reason, using a symmetrical group acting on players #0 - #2
	// does not add any new rejections for numPlayers 15 or 21 cases.
	// But just in case, we will keep the following fragment...
	T permPlayers[3], orbits[3];
	memcpy(permPlayers, getMatrixRow(0), lenGroup);
	memcpy(orbits, permPlayers, lenGroup);
	memcpy(pTmp, pntr, lenRow());
	while (true) {
		const auto nElem = nextPermutation(permPlayers, orbits, groupSize(), idx, 0);
		if (nElem == ELEMENT_MAX)
			break;

		for (T i = 0; i < 3; i++)
			memcpy(pTmp + i * groupSize(), pntr + permPlayers[i] * groupSize(), lenGroup);

		recordTuples(pTmp);
		sortTuples();

		const int diff = checkDayCode(0, NULL, 1);
		if (diff < 0) {
			if (!resultOut())
				return false;       // we don't need explanation for rejection

			*pNumReason = t_RejectionRreason::t_Statement_19;
			*pNumPlayer = i;		// the number of the group used
			memcpy(m_players, pTmp, lenRow());
			return checkRemainingDays(1, diff, m_players);
		}
	}
#endif
#endif
	return true;
}

/*
{
	CheckerCanon(bool)::checkWithGroup(const T* result, T numElem)
}
*/

CheckerCanon(void)::createDaySequence(T iDay) const {
	// Adding all unused days to the array.
	auto pDest = destMemory() + lenResult();
	*pDest++ = iDay;
	*pDest = 0;
	int j = 0;
	while (++j < numDays()) {
		if (j != iDay)
			*++pDest = j;
	}
}

CheckerCanon(int)::checkDayCode(int diff, T* pNumReason, T iDay) {
	if (!diff)
		diff = memcmp(m_players, studiedMatrix() + m_numElem, lenRow());

	if (!diff && numDays() > 2 || diff < 0 && resultOut()) {
		T* pDest = destMemory();
		if (diff != -9999) {
			// Saving two first days:
			memcpy(pDest, studiedMatrix(), lenRow());
			memcpy(pDest + m_numElem, m_players, lenRow());
		}

		if (diff < 0) {
			if (diff == -9999) {
				// Adding 2 day indices
				*(pDest += lenResult()) = iDay;
				*++pDest = 0;
				reportTxtError(pDest - 1, reason[*pNumReason], pDest - 1, 2);
			}
			// Andrei: When these 2 lines are commented out the AND numDays()==2
			// improved matrix is not printed 
			else
				addImproveResultFlags(t_bResultFlags::t_readyToExplainMatr);
		}

		// Adding all unused days to the array
		createDaySequence(iDay);
	}

	return diff;
}

CheckerCanon(bool)::checkRemainingDays(T iDay, int retVal, const T *pPerm) {
	// A function to record the remaining days and check them when the code for 
	// the first two days was built exactly the same as in the input matrix.
	// 
	const auto numElem_2 = 2 * m_numElem;
	auto* pDest = destMemory();
	if (iDay) {
		// Renumbering of players according to permutaion of days: (0, iDay).
		if (!pPerm)
			pPerm = getMatrixRow(iDay);

		T* pOrbits = m_players + m_numElem;
		for (auto j = m_numElem; j--;)
			pOrbits[pPerm[j]] = j;

		// Renumbering the set of players in the groups for all days except 0 and iDay.
		auto* pOut = pDest + 2 * m_numElem;
		auto* pIn = studiedMatrix();
		for (int j = 1; j < numDays(); j++) {
			pIn += m_numElem;
			if (j == iDay)
				continue;   // Skip day which is already used as 0's one. 

			for (T i = 0; i < m_numElem; i++)
				pOut[i] = pOrbits[pIn[i]];

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

CheckerCanon(int)::checkDay_1(int iDay, T* pNumReason) {
	// iDay index if the day which replaced the day 0
	int diff = 0;
#if	(UsePos_1_4_condition & 2)
#if UsePos_1_4_condition && ImproveResults
	if (!checkPosition1_4(m_players, pNumReason)){
		diff = -9999;
		if (numDays() == 2)
			return diff;
	}
#else
	if (m_players[4] != 4 && m_players[4] != 9) {
		diff = -9999;
		return diff;
	}
	if (m_players[4] > m_players[7])
		return -1;
#endif
#endif
	return checkDayCode(diff, pNumReason, iDay);
}

CheckerCanon(int)::orderingMatrix(T nDays, T numGroups, T* pNumReason, bool expected, bool invert, const T* permPlayer) {
	T* pDest = resultOut();
	if (!pDest)
		pDest = resultMemory();

	memcpy(pDest, studiedMatrix(), lenResult() * sizeof(*pDest));
	if (invert) {
		// Transformation: player #i new number is: numElem() - 1 - i
		auto pTmp = m_players + m_numElem;
		const auto j = numElem() - 1;
		for (auto i = numElem(); i--;)
			pTmp[i] = j - i;

		permPlayer = pTmp;
	}

	if (permPlayer) {
		for (auto i = lenResult(); i--;)
			pDest[i] = permPlayer[pDest[i]];
	}

	auto pDays = pDest + lenResult();
	for (auto j = numDays(); j--; )
		pDays[j] = j;

	orderigRemainingDays(nDays, numGroups, pDest);
	const auto diff = memcmp(pDest, studiedMatrix(), lenResult() * sizeof(*pDest));
	if (diff < 0) {
		if (pNumReason) {
			*pNumReason = invert ? t_RejectionRreason::t_invertOrdering : t_RejectionRreason::t_ordering;
			addImproveResultFlags(t_bResultFlags::t_readyToExplainMatr);
		}
		return diff;
	}

	assert(expected || diff >= 0);
	return diff;
}

CheckerCanon(bool)::checkDay(T iDay, T *pNumReason, T* pNumPlayer) {
	auto pMatrixRow = getMatrixRow(iDay);
	T* pDest = resultOut();
#if	(UsePos_1_4_condition & 1)
	if (iDay == 1) {
		if (!preordered()) {
			assert(pDest);  // When matrices are not preordered, 
				            // we expect to have an external buffer
			setPreordered(true);
			if (orderingMatrix(0, 0, pNumReason, false) < 0)
				return false;
		}

		if (!checkPosition1_4(pMatrixRow, pNumReason, pNumPlayer))
			return false;

		if (orderingMatrix(0, 0, pNumReason, false, true) < 0)
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

	return orderingMatrix(iDay, j, pNumReason) < 0;
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

CheckerCanon(bool)::permutPlayers4Day(const T* p_players, const T* res, T numGroup, T* resDayOut) const {
	// Function re-arranges players on a given day
	const auto* resEnd = res;
	for (T j = 0; j < numGroup; j++) {
		resEnd += groupSize();
		while (++res < resEnd) {
			const auto diff = (int)*(res - 1) - *res;
			if (diff > 0)
				return false;
		}

		if (j) {
			if (*(res - 1) < *(res - groupSize() - 1))
				return false;
		}
		//	void groupOrdering(T * resPerm, T numElem, T groupSize)

/*
		for (T j = 0; j < groupSize(); j++, res++) {
			const auto diff = (int)*(res - 1) - *res;
			assert(diff);
			if (diff > 0)
				return false;

			p_players[res[j]] = j + inc;
		}
		*/
	}
	return true;
}

CheckerCanon(bool)::reportTxtError(T *bBuffer, const char *pReason, T *pDays, T nDay) {
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
	pBuff += sprintf_s(pBuff, len - (pBuff - pBuffer), ". The problem was detected for day");
	if (pDays && nDay > 1) {
		pBuff += sprintf_s(pBuff, len - (pBuff - pBuffer), "s: (%d", pDays[0]);
		for (int i = 1; i < nDay; i++)
			pBuff += sprintf_s(pBuff, len - (pBuff - pBuffer), ", %d", pDays[i]);
			
		pBuff += sprintf_s(pBuff, len - (pBuff - pBuffer), ").");
	} else
		pBuff += sprintf_s(pBuff, len - (pBuff - pBuffer), " %d", nDay);

	return false;
}

CheckerCanon(T)::nextPermutation(T* perm, const T* pOrbits, T nElem, T idx, T lenStab) {
	// Function generates next permutation for the k-system 
	// Find non-increasing suffix
	//const auto nRow = numRow();
	const auto stabLenght = CGroupOrder<T>::stabilizerLength();
	T temp, j, i = stabLenght;

	// Check if the algorithm, used immediately after 
	// some non-trivial automorphism was found
	if (idx == IDX_MAX && perm[i] == nElem - 1)
		idx = ELEMENT_MAX;

	if (idx == IDX_MAX) {
		// Firts call after some automorphism was found
		temp = perm[idx = i];
		for (j = nElem; --j > temp;)
			perm[j] = j;

		for (auto k = j++; k-- > i;)
			perm[k + 1] = k;
	}
	else {
		if (idx >= IDX_MAX) {
			j = i = nElem;
			while (--i > 0 && perm[i - 1] >= perm[i]);

			if (i == lenStab)
				return ELEMENT_MAX;  // no more permutations

			// Find successor to pivot
			temp = perm[--i];
			while (perm[--j] <= temp);
		}
		else {
			temp = perm[j = i = idx];
			while (++j < nElem && perm[j] <= temp);
			if (j >= nElem) {
				if (nElem > i - 2)
					revert(perm, nElem, i);

				return nextPermutation(perm, pOrbits, nElem);
			}
		}
	}

	if (stabLenght == i) {
		bool flag = false;
		auto k = j, tmp = perm[j];
		if (idx >= IDX_MAX) {
			while (k > i && *(pOrbits + perm[k]) != perm[k])
				k--;

			if (k != j) {
				if (!k)
					return ELEMENT_MAX;

				flag = k == i;
				tmp = perm[k--];
				while (++k < j)
					perm[k] = perm[k + 1];
			}
		}
		else {
			while (k < nElem && *(pOrbits + perm[k]) != perm[k])
				k++;

			if (k != j) {
				flag = k == nElem;
				if (flag) {
					if (!i)
						return ELEMENT_MAX;

					// Re-establish trivial permutation
					k = idx - 1;
					while (++k < j)
						perm[k] = k;
				}
				else {
					tmp = perm[k++];
					while (--k > j)
						perm[k] = perm[k - 1];
				}
			}
		}

		perm[j] = tmp;
		if (flag) {
			j = idx >= IDX_MAX ? nElem - 1 : i;
			temp = perm[--i];
			CGroupOrder<T>::setStabilizerLength(i);
		}
	}

	if (!(perm[j] % 3) && i < j) {
		i = i;
	}
	perm[i] = perm[j];
	perm[j] = temp;
	if (idx >= IDX_MAX - 1) {
		if (CGroupOrder<T>::stabilizerLength() > i)
			CGroupOrder<T>::setStabilizerLength(i);

		if (nElem > i - 2)
			revert(perm, nElem, i);
	}

	return i;
}

CheckerCanon(bool)::CheckPermutations(const T* result, const T* pMatrix, int nRows) {
	T* permPlayers = m_players;
	memcpy(permPlayers, result, lenRow());
	auto pOrbits = new T [numElem()];
	memcpy(permPlayers, result, lenRow());
	memcpy(pOrbits, result, lenRow());
	bool checkPermut = false;
	int errLine, errGroup, dubLine;
	size_t cntr = 1;
	size_t counter[2] = { 0 };
	char* pMatrixOut = (char*)(pMatrix);
	char *pDayPerm = pMatrixOut + nRows * numElem();
	T nRow = 0;
	T lenStab = 0;
	T idx = ELEMENT_MAX;
	printf(" I am going to the loop:\n");

	while (true) {
		if (_CheckMatrix(pMatrixOut, nRows, numElem(), true, &errLine, &errGroup, &dubLine)) {
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

		const auto nElem = nextPermutation(permPlayers, pOrbits, numElem(), idx, lenStab);
		if (nElem == ELEMENT_MAX)
			break;

		cntr++;
		break;
		auto pntrTo = (char *)pMatrix + numElem();
		for (T i = 2; i < nRows; i++) {
			memcpy(pntrTo += numElem(), result + numElem() * pDayPerm[i], lenRow());
			for (T j = 0; j < numElem(); j++)
				pntrTo[j] = permPlayers[pntrTo[j]];
		}

		//			continue;
				//	if (checkPermut = PermutResults(result, permPlayers, nDays))
				//		break;

					// Under the action of the current permutation of players, the set of
					// days is the same as before - an automorphism has been discovered
		//addAutomorphism(nElem = m_numElem, permPlayers, pOrbits);
	}

	printf("Good and Bad counters are %zd - %zd\n", counter[0], counter[1]);
	delete[] pOrbits;
	return counter[1] == 0;
}

