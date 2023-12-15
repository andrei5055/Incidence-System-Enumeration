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
		: m_numElem(nCol), m_numDaysMax(nRow), m_lenResult(nCol * nRow), m_groupSise(groupSize) {
		m_pDayRes = new T[2 * (m_numElem + nRow)];
	}
	~CCanonicityChecker()					{ delete[] m_pDayRes; }
	bool CheckCanonicity(const T* result, int nLines, T *bResult=NULL);
private:
	inline T groupSize() const				{ return m_groupSise; }
	auto stabiliserLengthExt() const		{ return m_nStabExtern; }
	void setStabiliserLengthExt(T len)		{ m_nStabExtern = len; }
	bool copyTuple(const T* res, T* p_players, T inc = 0) const;
	bool rollBack(T* p_dayRes, T* p_dayIsUsed, int& j, int nDays) const;
	inline void setNumDays(T nDays)			{ m_numDays = nDays; }
	inline auto numDays() const				{ return m_numDays;}
	bool checkDay_1(const T* result, const T* p_players, T k, T j, T* pDest, T* p_dayRes, const T* p_dayIsUsed) const;

	T m_nStabExtern = 0;		// number of first elements of permutation which Canonicity Checker will not move
	T* m_pDayRes = NULL;
	const T m_numElem;			// The number of elements that will be the same for all partially constructed objects
								// (it is equal nCol for combinatorial designs or number of players for k-system)
	const T m_numDaysMax;
	const T m_groupSise;
	unsigned int m_lenResult;
	T m_numDays;
	int a = 0;
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
	//return true;
	/*
	"     0   1   2    3   4   5    6   7   8 \n"
	"     0   3   6    1   4   8    2   5   7 \n"
	*/
#if 0
		// return true  - continue, 
		//        false - stop (calculate new matrix)
		static bool ones = false;
		if (nLines > 5 && !ones)
		{
			ones = true;
			//return false;
		}
		return true;
#else
	const auto lenStab = stabiliserLengthExt();
	size_t startIndex = 0;
	T* permColumn = NULL;
	const auto* res = result;
	const auto lenMemory = m_numElem + 2 * nDays;
	T* p_dayRes = m_pDayRes;
	T* p_dayIsUsed = p_dayRes + nDays;
	T* p_players = p_dayIsUsed + nDays;
	T* pOrbits = p_players + m_numElem;

	if (bResult)
		setNumDays(nDays);
	
	if (result[19] == 9 && result[20] == 12)
		nDays = nDays;
	const auto lenGroup = groupSize();
	const auto numGroup = m_numElem / lenGroup;
	a++;
	if (a >= 1031348)
		a += 0;

	for (int iDay = 0; iDay < nDays; iDay++, res += lenGroup) {
		auto *pDest = bResult;
		if (res[0] || !copyTuple(res, p_players))
			return false;

		T inc = 0;
		for (auto j = numGroup; --j;) {
			if (!copyTuple(res += lenGroup, p_players, inc += lenGroup) ||
				p_players[*res] < p_players[*(res - lenGroup)]) // Comparing first elements of the groups
				return false;
		}

		// Check canonicity of the codes for the other days
		memset(p_dayRes, 0, 2 * nDays * sizeof(*p_dayRes));
		p_dayIsUsed[p_dayRes[0] = iDay] = 1;
		if (iDay) {
			// Do this only when day 0 changed its place.
			groupOrdering(p_players, m_numElem, lenGroup);

			p_dayIsUsed[p_dayRes[1] = 0] = 1;
			if (!checkDay_1(result, p_players, 0, 1, pDest, p_dayRes, p_dayIsUsed))
				return false;
		}

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
				if (!checkDay_1(result, p_players, k, j, pDest, p_dayRes, p_dayIsUsed)) {
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
#if 0
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
#endif
	}

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
#endif
}

CanonicityChecker(bool)::copyTuple(const T* res, T* p_players, T inc) const {
	p_players[res[0]] = inc;
	for (T j = 1; j < groupSize(); j++) {
		const auto diff = (int)res[j - 1] - res[j];
		assert(diff);
		if (diff > 0)
			return false;

		p_players[res[j]] = j + inc;
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

CanonicityChecker(bool)::checkDay_1(const T *result, const T *p_players, T k, T j, T* pDest, T *p_dayRes, const T *p_dayIsUsed) const
{
	const auto* resDayPerm = result + k * m_numElem;
	const auto* resDay = result + j * m_numElem;
	int diff = 0;
	T t = -1;
	while (++t < m_numElem && !(diff = (int)resDayPerm[p_players[t]] - resDay[t]));
	if (t >= m_numElem || diff >= 0)
		return true;

	if (pDest) {
		// Saving two first days:
		memcpy(pDest, result, m_numElem * sizeof(*pDest));
		memcpy(pDest += m_numElem, p_players, m_numElem * sizeof(*pDest));

		// Adding all unused days to the array.
		k = -1;
		while (++j < numDays()) {
			while (p_dayIsUsed[++k]);
			p_dayRes[j] = k;
		}
		memcpy(pDest += m_numElem, p_dayRes, numDays() * sizeof(*pDest));
	}

	return false;
}