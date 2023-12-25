
#include <assert.h>
#include "CheckCanon.h"
typedef enum {
	t_reasonUnknown,
	t_ordering,
	t_playrPosition_1_4,
	t_NotThatPlayerInPosition_1_4,
} t_RejectionRreason;

static const char* reason[] = {
		"Reason unknown",
		"Ordering problem",
		"Only players 4 or 9 can be in position [1,4]",
		"Current player cannot be at position [1, 4]",
};

template class CCheckerCanon<unsigned char, unsigned char>;

template<typename T>
void renumberPlayers(T* pntr, size_t i, size_t iLast) {
	for (; i < iLast; i++)
		pntr[i] = pntr[pntr[i]];
}

template<typename T>
void elemOrdering(T* pElems, size_t numElem, size_t groupSize)
{
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

	static char specialReason[256];

	T* pOrbits = m_players + m_numElem;
	T numReason;

	setResultOut(bResult);
	setNumDays(nDays);
	resetImprovedResultFlag();

	auto result2 = result + m_numElem2; 
	const auto lenCmp = (lenResult() - m_numElem2) * sizeof(*result);

	const auto lenGroup = groupSize();
	const auto numGroup = m_numElem / lenGroup;
	auto* pDest = bResult ? bResult : resultMemory();
	auto* res = result;
	bool copyTuplesOK = true;
	for (int iDay = 0; iDay < nDays; iDay++, res += lenGroup) {
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
			for (auto j = numGroup; --j;) {
				if (!copyTuple(res += lenGroup, inc += lenGroup) ||
					m_players[*res] < m_players[*(res - lenGroup)]) // Comparing first elements of the groups
					copyTuplesOK = false;
			}
		}

		if (!copyTuplesOK) {
			return false;
			// Some problems found when we had tr
			res = result + iDay * m_numElem;
	//		for (auto )
		}

		// Check canonicity of the codes for the other days
		if (iDay) {
			// Do this only when day 0 changed its place.
			elemOrdering(m_players, m_numElem, lenGroup);
			groupOrdering(m_players, numGroup, getTmpBuffer(), lenGroup);

			const auto retVal = checkDay_1(result, iDay, pDest, &numReason);
			if (nDays == 2)
				return retVal >= 0; // In this case there is nothing more to do than has already been done.

			if (retVal > 0)
				continue;

			if (retVal < 0 && !bResult)
				return false;  // The result has improved, but we don't need to know how.

			// Renumbering of players according to permutaion of days: (0, iDay).
			auto* pPerm = result + iDay * m_numElem;
			for (auto j = m_numElem; j--;)
				pOrbits[pPerm[j]] = j;

			// Renumbering the set of players in the groups for all days except 0 and iDay.
			auto* pOut = pDest + 2 * m_numElem;
			auto* pIn = result;
			for (int j = 1; j < nDays; j++) {
				pIn += m_numElem;
				if (j == iDay)
					continue;   // Skip day which is already used as 0's one. 

				for (T i = 0; i < m_numElem; i++)
					pOut[i] = pOrbits[pIn[i]];

				pOut += m_numElem;
			}

			orderigRemainingDays(2, 0, pDest);

			// Comparing all remaining days		
			if (retVal < 0 || USE_2_ROW_CANON == 0 && memcmp(pDest + m_numElem2, result2, lenCmp) < 0) {
				if (retVal == -9999) {
					const auto numElem_2 = 2 * m_numElem;
					renumberPlayers(pDest, numElem_2, lenResult() - numElem_2);
					memcpy(pDest, result, m_numElem * sizeof(*pDest));
				}

				addImproveResultFlags(t_bResultFlags::t_readyToExplainMatr);
				return false;
			}
		}
		else {
			// Check all remaining days for canonicity.
			for (int j = 1; j < nDays; j++) {
				if (!checkDay(result, j, &numReason)) {
					return reportTxtError(bResult, reason[numReason], NULL, j);
				}
			}
		}
	}

#if 0
		if (USE_2_ROW_CANON)
			continue;

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
				if (!checkDay_1(result, p_players, k, j, p_dayRes, p_dayIsUsed)) {
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

CheckerCanon(bool)::checkPosition1_4(const T* result, const T *players, T playerID, T *pNumReason) {
	if (players[4] == playerID) {
		auto* pDest = resultOut();
		if (pDest) {
			memcpy(pDest, result, m_numElem * sizeof(*pDest));
			memcpy(pDest + m_numElem, players, m_numElem * sizeof(*pDest));
			pDest[7] = playerID;
			pDest[playerID] = 7;

			const auto numElem_2 = 2 * m_numElem;
			renumberPlayers(pDest, m_numElem, numElem_2);
			groupOrdering(pDest + m_numElem, numGroups(), getTmpBuffer(), groupSize());
			const auto diff = memcmp(pDest + m_numElem, players, m_numElem);
			if (diff < 0) {
				*pNumReason = t_RejectionRreason::t_NotThatPlayerInPosition_1_4;
				if (result + m_numElem == players || numDays() == 2) {
					addImproveResultFlags(t_bResultFlags::t_readyToExplainMatr);
					memcpy(pDest + numElem_2, result + numElem_2, (lenResult() - numElem_2) * sizeof(*pDest));
					renumberPlayers(pDest, numElem_2, lenResult() - numElem_2);

					*(pDest + lenResult()) = 0;
					*(pDest + lenResult() + 1) = 1;
					orderigRemainingDays(2, 0, pDest);
					memcpy(pDest, result, m_numElem * sizeof(*pDest));
				}
			}
			else {
				assert(false);  // We should not be here
			}
		}

		return false;
	}

	return true;
}

CheckerCanon(int)::checkDay_1(const T *result, int iDay, T* pDest, T* pNumReason) {
	// iDay index if the day which replaced the day 0
	int diff = 0;
	const auto* resDay = result + m_numElem;
#if	(UsePos_1_4_condition & 2)
	if (false && !checkPosition1_4(result, m_players, 8, pNumReason)){
		diff = -9999;
		if (numDays() == 2)
			return diff;
	}
	//if (!checkPosition1_4(result, result + m_numElem, 8, pNumReason))
	/*
	if (m_players[4] != 4 && m_players[4] != 9) {
		diff = -9999;
	}
	*/
#endif

	if (!diff) {
		T t = -1;
		while (++t < m_numElem && !(diff = (int)m_players[t] - resDay[t]));
	}

	if (!diff && numDays() > 2 || diff < 0 && resultOut()) {
		if (diff != -9999) {
			// Saving two first days:
			memcpy(pDest, result, m_numElem * sizeof(*pDest));
			memcpy(pDest + m_numElem, m_players, m_numElem * sizeof(*pDest));
		}

		// Adding 2 day indices
		*(pDest += lenResult()) = iDay;
		*++pDest = 0;
		if (diff < 0) {
			if (diff == -9999)
				reportTxtError(pDest - 1, reason[*pNumReason], pDest - 1, 2);
			// Andrei: When these 2 lines are commented out the AND numDays()==2
			// improved matrix is not printed 
			else  
				addImproveResultFlags(t_bResultFlags::t_readyToExplainMatr);
		}

		// Adding all unused days to the array.
		int j = 0;
		while (++j < numDays()) {
			if (j != iDay)
				*++pDest = j;
		}
	}

	return diff;
}


CheckerCanon(bool)::checkDay(const T* result, T iDay, T *pNumReason) {
	const auto* res = result + iDay * m_numElem;
	const auto* resEnd = res;
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
	if (iDay == 1) {
		res = result + m_numElem;
		if (!checkPosition1_4(result, result + m_numElem, 8, pNumReason))
			return false;
			/*
		if (res[4] != 4 && res[4] != 9) {
			*pNumReason = t_RejectionRreason::t_playrPosition_1_4;
			return false;
		} */
	}
#endif

	if (j == numGroups())
		return true;

	T* pDest = resultOut();
	if (pDest) {
		// Copying all days to the output
		memcpy(pDest, result, lenResult() * sizeof(*pDest));
		pDest += lenResult();
		for (auto i = m_numElem; i--;)
			*(pDest + i) = i;
		// Unfinished code....
		orderigRemainingDays(iDay, j, pDest);
	}

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
		groupOrdering(pntr, numGroups(), getTmpBuffer(), groupSize());
		pntr += m_numElem;
	}

	const auto daysToOrder = numDays() - daysOK;
	if (daysToOrder) {
		// Ordering days by the first two elements of their first group
		groupOrdering(pDest + lenOK, daysToOrder, getTmpBuffer(), m_numElem, pDest + lenResult() + daysOK);
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

