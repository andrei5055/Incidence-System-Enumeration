
#include <assert.h>
#include "CheckCanon.h"
typedef enum {
	t_reasonUnknown,
	t_ordering,
	t_playerPosition_1_4,
	t_NotThatPlayerInPosition_1_4,
	t_Statement_7,
	t_Statement_18,
	t_Statement_19,
} t_RejectionRreason;

static const char* reason[] = {
		"Reason unknown",
		"Ordering problem",
		"Only players 4 or 9 can be in position [1,4]",
		"Player #%d cannot be at position [1, 4]",
		"Player# in [1, 4] should be less than [1, 7]",
		"Incorrect order of players in a group with a leading player #%d",
		"Rejected by Statement 19"
};

template class CCheckerCanon<unsigned char, unsigned char>;

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

	T numReason;

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
				return false;
		}
		else {
			// Check all remaining days for canonicity.
			T playerNumb = -1;
			for (int j = 1; j < nDays; j++) {
				if (!checkDay(j, &numReason, &playerNumb)) {
					char buffer[256];
					auto pReason = reason[numReason];
					if (playerNumb != -1) {
						sprintf_s(buffer, pReason, playerNumb);
						pReason = buffer;
					}

					return reportTxtError(bResult, pReason, NULL, j);
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

		// Not ready yet
		T* permPlayers = NULL;
		permColumn = init(m_numElem, 0, false, pOrbits, &permPlayers, false, permColumn);

		T nElem = ELEMENT_MAX;
		/*
				if (check_trivial_row_perm) {
					nRow = startingRowNumb;
					goto try_permut;
				}
				*/
		while (true) {

			//		next_permut:
			nElem = nextPermutation(permPlayers, pOrbits, nElem, lenStab);
			if (nElem == ELEMENT_MAX || nElem < lenStabilizer())
				break;

			for (T iDay = 0; iDay < nDays; iDay++) {
				const auto* pDayRes = result + iDay * m_numElem;
				for (; nElem < m_numElem; nElem++) {
					// const auto* pRow = pMatr->GetRow(nElem);
				}
			}
		}
	}
#endif

#if 0
	// Not ready yet
	for (auto j = nDays; --j;) {
		// Check the possibility of modifying the last triple of the j-th day by the conversion 
		// i ==> (m_numPlayers - i - 1) to get the first triple of 1-st day: (0, 3, 6).  
		const auto* pntr = result + m_numElem * (j + 1) - 1;
		if (*pntr != (m_numElem - 1) || *(pntr - 1) != (m_numElem - 4) || *(pntr - 2) != (m_numElem - 7))
			continue; // It didn't happen on day j.

		// Try the same conversion on a full set of players and reorder new triples of the day j.
		pntr -= 3;
		const auto val = m_numElem - 1;
		const auto iLast = m_numElem - 3;
		for (T i = 0; i < iLast; i++)
			p_players[i] = val - *pntr--;

		const auto retVal = memcmp(p_players, result + m_numElem + 3, iLast * sizeof(*p_players)) >= 1;
		if (!retVal && bResult) {
			auto i = m_numElem;
			while (i--)
				bResult[i] = i;//  val - i;

			bResult[m_numElem] = 0;// j;
			bResult[m_numElem + 1] = j;//0;
			T k = 0;
			for (T i = 2; i < nDays; i++) {
				if (++k == j)
					k++;
				bResult[m_numElem + i] = k;
			}	
		}
		/*  val - i (0, j)
		 "    0  1  2    3  4  5    6  7  8    9 10 11   12 13 14 \n"
		 "    0  3  6    1  7 12    2  8  9    4 10 14    5 11 13 \n"
		 "    0  4  7    1  3  8    2 10 13    5  9 14    6 11 12 \n"
		 "    0  5  8    1  4 11    2  7 14    3 10 12    6  9 13 \n"
		 "    0  9 12    1  5 10    2  4  6    3  7 13    8 11 14 \n"
Improved Result #1: val - i (0, j)
		"   14 13 12   11 10  9    8  7  6    5  4  3    2  1  0 \n"
        "   14  5  2   13  9  4   12 10  8   11  7  1    6  3  0 \n"
Improved Result #1: i (0, j)
 "    0  1  2    3  4  5    6  7  8    9 10 11   12 13 14 \n"
 "    0  9 12    1  5 10    2  4  6    3  7 13    8 11 14 \n"
 "    0  3  6    1  7 12    2  8  9    4 10 14    5 11 13 \n"
 "    0  4  7    1  3  8    2 10 13    5  9 14    6 11 12 \n"
 "    0  5  8    1  4 11    2  7 14    3 10 12    6  9 13 \n"

 */
		return retVal;
	}
#endif
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
	bool copyTuplesOK = true;
	if (res[0] || !copyTuple(res)) {
		/*			T* pDest = resultOut();
					if (pDest) { //???
						memcpy(pDest, result, (iDay + 1) * m_numElem * sizeof(*pDest));
						elemOrdering(pDest += iDay * m_numElem, m_numElem, groupSize());
						memcpy(pDest += m_numElem, result + m_numElem * (iDay + 1), (numDays() - iDay) * m_numElem * sizeof(*pDest));
					}
					*/
		copyTuplesOK = false;
	}

	if (copyTuplesOK) {
		T inc = 0;
		for (auto j = numGroups(); --j;) {
			if (!copyTuple(res += groupSize(), inc += groupSize()) ||
				m_players[*res] < m_players[*(res - groupSize())]) // Comparing first elements of the groups
				copyTuplesOK = false;
		}
	}

	return copyTuplesOK;
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
	// List of simple player substitutions (subst[i] <---> subst[i+1], for i%2 == 0) 
	// which will improve the matrix code, if player subst[i] is at position [1, 2]
	static T subst[] = { 8, 7,10, 9,11, 9 };
	// aaa_3118.txt     subst = { 8, 7 }   file for UsePos_1_4_condition = 1
	// aaa_3118++.txt	subst = { 8, 7,10, 9,11, 9 };
	// aaa_2469.txt     using everything from this method, but only for UsePos_1_4_condition = 1

	for (int i = 0; i < countof(subst); i += 2) {
		const auto playerID = subst[i];
		if (players[4] == playerID) {  // player is on position #4 of day #1
			*pNumReason = t_RejectionRreason::t_NotThatPlayerInPosition_1_4
			if (pNumPlayer)
				*pNumPlayer = subst[i];
			return explainRejection(players, playerID, subst[i + 1]);
		}
	}
#endif
	if (players[4] == 7) {
		if (!resultOut())
			return false;              // we don't need explanation for rejection

		*pNumReason = t_RejectionRreason::t_NotThatPlayerInPosition_1_4;
		if (pNumPlayer)
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
	// Swaping group 0 of day 1 with the other groups
	auto pntr = getMatrixRow(1);
	auto pntrFrom = pntr;
	auto pntrTo = tmpBuffer();
	for (T i = 1; i < numGroups(); i++) {
		memcpy(tmpBuffer(), pntr, lenRow());
		memcpy(tmpBuffer(), pntrFrom += groupSize(), groupSize() * sizeof(T));
		memcpy(pntrTo += groupSize(), pntr, groupSize() * sizeof(T));
		recordTuples(tmpBuffer());
		sortTuples();

		const int diff = checkDayCode(0, NULL, 1);
		if (diff < 0) {
			*pNumReason = t_RejectionRreason::t_Statement_19;
			checkOrderingForDay(1); // days 0 and 1
			return checkRemainingDays(1, diff);
		}
	}
#endif
	return true;
}

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

CheckerCanon(bool)::checkRemainingDays(T iDay, int retVal, bool recordByMatrixRow) {
	// A function to record the remaining days and check them when the code for 
	// the first two days was built exactly the same as in the input matrix.
	// 
	const auto numElem_2 = 2 * m_numElem;
	auto* pDest = destMemory();
	if (iDay) {
		// Renumbering of players according to permutaion of days: (0, iDay).
		T* pOrbits = m_players + m_numElem;
		auto* pPerm = recordByMatrixRow? getMatrixRow(iDay) : m_players;
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

CheckerCanon(bool)::checkDay(T iDay, T *pNumReason, T* pNumPlayer) {
	auto pMatrixRow = getMatrixRow(iDay);
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


#if	(UsePos_1_4_condition & 1)
	if (iDay == 1 && !checkPosition1_4(pMatrixRow, pNumReason, pNumPlayer))
		return false;
#endif

	if (j == numGroups())
		return true;

	T* pDest = resultOut();
	if (pDest) {
		// Copying all days to the output
		memcpy(pDest, studiedMatrix(), lenResult() * sizeof(*pDest));
		pDest += lenResult();
		for (auto i = m_numElem; i--;)
			*(pDest + i) = i;
		// Unfinished code....
		orderigRemainingDays(iDay, j, pDest);
	}

	addImproveResultFlags(t_bResultFlags::t_readyToExplainMatr);
	*pNumReason = t_RejectionRreason::t_ordering;
	return false;
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

