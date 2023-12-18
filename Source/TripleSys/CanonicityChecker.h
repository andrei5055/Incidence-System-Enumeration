#pragma once
#include <assert.h>

#define TFunc2(x, ...)          template<typename T, typename S> __VA_ARGS__ x
#define Class2(x)               x<T,S>
#define Class2Def(x)            TFunc2(x, class)
#define FClass2(x, ...)			TFunc2(Class2(x), __VA_ARGS__)
#define CanonicityChecker(...)	FClass2(CCanonicityChecker, __VA_ARGS__)

#define countof(x)     sizeof(x)/sizeof(x[0])


Class2Def(CCanonicityChecker) {
public:
	CCanonicityChecker(T nRow, T nCol, T groupSize = GroupSize)
		: m_numElem(nCol), m_numElem2(2 * nCol), m_numDaysMax(nRow), m_groupSise(groupSize) {
		m_players = new T[2 * m_numElem];
		m_tmpBuffer = new T[m_numElem + nRow];
		m_pResutMemory = new T[(m_numElem + 1) * nRow];
	}
	~CCanonicityChecker()					{ delete[] m_players;
											  delete [] getTmpBuffer();
											  delete[] resultMemory();
	}
	bool CheckCanonicity(const T* result, int nLines, T *bResult=NULL);
private:
	inline auto groupSize() const			{ return m_groupSise; }
	auto stabiliserLengthExt() const		{ return m_nStabExtern; }
	void setStabiliserLengthExt(T len)		{ m_nStabExtern = len; }
	bool copyTuple(const T* res, T inc = 0) const;
	bool rollBack(T* p_dayRes, T* p_dayIsUsed, int& j, int nDays) const;
	inline void setNumDays(T nDays)			{ m_lenResult = (m_numDays = nDays) * m_numElem; }
	inline auto numDays() const				{ return m_numDays; }
	inline void setResultOut(T* pntr)		{ m_pResultOut = pntr; }
	inline auto resultOut() const			{ return m_pResultOut; }
	inline auto getTmpBuffer() const		{ return m_tmpBuffer; }
	inline auto resultMemory() const		{ return m_pResutMemory;  }
	inline auto lenResult()	const			{ return m_lenResult; }
	int checkDay_1(const T* result, int iDay, T *pDest) const;
	bool checkDay(const T* res, T iDay, T numGroup) const;
	void orderigRemainingDays(T daysOK, T groupsOK, T numGroup, T *pDest) const;
	bool permutPlayers4Day(const T* p_players, const T* resDayIn, T numGroup, T* resDayOut) const;

	T m_nStabExtern = 0;		// number of first elements of permutation which Canonicity Checker will not move
	T* m_players = NULL;
	const T m_numElem;			// The number of elements that will be the same for all partially constructed objects
								// (it is equal nCol for combinatorial designs or number of players for k-system)
	const T m_numElem2;			// This is twice the number of elements. 
	const T m_numDaysMax;
	const T m_groupSise;
	size_t m_lenResult;
	T m_numDays;
	T* m_pResultOut;
	T* m_tmpBuffer = NULL;		// Buffer uswd for groups and days ordering
	T* m_pResutMemory = NULL;	// Memory allocated to improve results
};

CanonicityChecker(bool)::CheckCanonicity(const T *result, int nDays, T *bResult) {
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

	T* pOrbits = m_players + m_numElem;

	setResultOut(bResult);
	setNumDays(nDays);
	auto result2 = result + m_numElem2; 
	const auto lenCmp = (lenResult() - m_numElem2) * sizeof(*result);

	const auto lenGroup = groupSize();
	const auto numGroup = m_numElem / lenGroup;
	auto* pDest = bResult ? bResult : resultMemory();
	auto* res = result;
	for (int iDay = 0; iDay < nDays; iDay++, res += lenGroup) {
		if (res[0] || !copyTuple(res)) {
			/*			T* pDest = resultOut();
						if (pDest) { //???
							memcpy(pDest, result, (iDay + 1) * m_numElem * sizeof(*pDest));
							elemOrdering(pDest += iDay * m_numElem, m_numElem, groupSize());
							memcpy(pDest += m_numElem, result + m_numElem * (iDay + 1), (numDays() - iDay) * m_numElem * sizeof(*pDest));
						}
						*/
			return false;
		}

		T inc = 0;
		for (auto j = numGroup; --j;) {
			if (!copyTuple(res += lenGroup, inc += lenGroup) ||
				m_players[*res] < m_players[*(res - lenGroup)]) // Comparing first elements of the groups
				return false;
		}

		// Check canonicity of the codes for the other days
		if (iDay) {
			// Do this only when day 0 changed its place.
			elemOrdering(m_players, m_numElem, lenGroup);
			groupOrdering(m_players, numGroup, getTmpBuffer(), lenGroup);

			const auto retVal = checkDay_1(result, iDay, pDest);
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

			orderigRemainingDays(2, 0, numGroup, pDest);

			// Comparing all remaining days		
			if (retVal < 0 || USE_2_ROW_CANON == 0 && memcmp(pDest + m_numElem2, result2, lenCmp) < 0)
				return false;
		}
		else {
			// Check all remaining days for canonicity.
			for (int j = 1; j < nDays; j++) {
				if (!checkDay(result, j, numGroup))
					return false;
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

CanonicityChecker(bool)::copyTuple(const T* res, T inc) const {
	m_players[res[0]] = inc;
	for (T j = 1; j < groupSize(); j++) {
		const auto diff = (int)res[j - 1] - res[j];
		assert(diff);
		if (diff > 0)
			return false;

		m_players[res[j]] = j + inc;
	}
	return true;
}

CanonicityChecker(bool)::rollBack(T* p_dayRes, T* p_dayIsUsed, int& j, int nDays) const
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

CanonicityChecker(int)::checkDay_1(const T *result, int iDay, T* pDest) const {
	// iDay index if the day which replaced the day 0
	const auto* resDay = result + m_numElem;
	int diff = 0;
	T t = -1;
	while (++t < m_numElem && !(diff = (int)m_players[t] - resDay[t]));

	if (!diff && numDays() > 2 || diff < 0 && resultOut()) {
		// Saving two first days:
		memcpy(pDest, result, m_numElem * sizeof(*pDest));
		memcpy(pDest + m_numElem, m_players, m_numElem * sizeof(*pDest));

		*(pDest += lenResult()) = iDay;
		*++pDest = 0;

		// Adding all unused days to the array.
		int j = 0;
		while (++j < numDays()) {
			if (j != iDay)
				*++pDest = j;
		}
	}

	return diff;
}

CanonicityChecker(bool)::checkDay(const T* result, T iDay, T numGroup) const {
	const auto* res = result + iDay * m_numElem;
	const auto* resEnd = res;
	T j = 0;
	for (; j < numGroup; j++) {
		resEnd += groupSize();
		while (++res < resEnd && *(res - 1) < *res);

		if (res < resEnd)
			break;

		// Check ordering of the first elements of the groups
		if (j && *(res - groupSize()) < *(res - 2 * groupSize()))
			break;
	}

	if (j == numGroup)
		return true;

	T* pDest = resultOut();
	if (pDest) {
		// Copying all days to the output
		memcpy(pDest, result, lenResult() * sizeof(*pDest));
		pDest += lenResult();
		for (auto i = m_numElem; i--;)
			*(pDest + i) = i;

		orderigRemainingDays(iDay, j, numGroup, pDest);
	}

	return false;
}

CanonicityChecker(void)::orderigRemainingDays(T daysOK, T groupsOK, T numGroup, T *pDest) const {
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
		groupOrdering(pntr, numGroup, getTmpBuffer(), groupSize());
		pntr += m_numElem;
	}

	const auto daysToOrder = numDays() - daysOK;
	if (daysToOrder) {
		// Ordering days by the first two elements of their first group
		groupOrdering<T>(pDest + lenOK, daysToOrder, getTmpBuffer(), m_numElem, pDest + lenResult() + daysOK);
	}
}

CanonicityChecker(bool)::permutPlayers4Day(const T* p_players, const T* res, T numGroup, T* resDayOut) const {
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

