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
	CCanonicityChecker(T nRow, T nCol, T groupSize=GroupSize)
		: m_numElem(nCol), m_groupSise(groupSize) {}
	bool CheckCanonicity(const T* result, int nLines, T *bResult=NULL);
private:
	inline T groupSize() const				{ return m_groupSise; }
	auto stabiliserLengthExt() const		{ return m_nStabExtern; }
	void setStabiliserLengthExt(T len)		{ m_nStabExtern = len; }
	const auto *orbits() const				{ return m_pObits; }
	bool copyTuple(const T* res, T* p_players, T inc = 0) const;
	bool rollBack(T* p_dayRes, T* p_dayIsUsed, int& j, int nDays) const;

	T m_nStabExtern = 0;		// number of first elements of permutation which Canonicity Checker will not move
	T* m_pObits[2][2];
	const T m_numElem;				   // The number of elements that will be the same for all partially constructed objects
									   // (it is equal nCol for combinatorial designs or number of players for k-system)
	const T m_groupSise;
};

CanonicityChecker(bool)::CheckCanonicity(const T *result, int nDays, T *bResult) {
	// Input parameters:
	//    result - pointer to a sequence of lists, each containing "m_numElem" players
	//             for each day, players are divided into groups, each of which contains "n"m_groupSise" players
	//    nDays  - number of days (or lists, mentioned above)
	//    bResult (optional) - pointer to the array of (m_numElem + nDays) elements,
	//             when it's not NULL and the "result" is not a canonical one, then
	//             this array will contain the permutations of
	//             (a) players (as its first "m_numElem" elements)
	//             (b) days (starting with (bResult+m_numElem)'s element)    
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
	static int cntr = 0;
	size_t startIndex = 0;
	T* permColumn = NULL;
	auto *pOrbits = orbits();
	const auto* res = result;
	const auto lenMemory = m_numElem + 2 * nDays;
	T buff[100];
	T* p_players = lenMemory <= countof(buff) ? buff : new T[lenMemory];
	T* p_dayRes = p_players + m_numElem;
	T* p_dayIsUsed = p_dayRes + nDays;

	const auto lenGroup = groupSize();
	const auto numGroup = m_numElem / lenGroup;
	for (int iDay = 0; iDay < nDays; iDay++, res += lenGroup) {
		if (res[0] || !copyTuple(res, p_players))
			return false;

		T inc = 0;
		for (auto j = numGroup; --j;) {
			if (!copyTuple(res += lenGroup, p_players, inc += lenGroup) ||
				p_players[*res] < p_players[*(res - lenGroup)]) // Comparing first elements of the groups
				return false;
		}

		if (lenGroup >= 3) {
			// Ordering last two elements of each 
			auto* res = p_players;
			for (auto j = numGroup; j--; res += lenGroup) {
				if (res[1] > res[2]) {  // TODO: Write code for lenGroup > 3
					const auto tmp = res[1];
					res[1] = res[2];
					res[2] = tmp;
				}
			}
		}

		// Check canonicity of the codes for the other days
		memset(p_dayRes, 0, 2 * nDays * sizeof(*p_dayRes));
		p_dayIsUsed[p_dayRes[0] = iDay] = 1;
		int j = 0;
		T k;
		while (true) {
			while (++j < nDays) {
				cntr++;
				// Looking for the first unused day
				k = p_dayRes[j];
				while (k < nDays && p_dayIsUsed[k])
					k++;

				if (k == nDays) // Unable to find unused day.
					break;

				p_dayIsUsed[p_dayRes[j] = k] = 1;
				const auto* resDayPerm = result + k * m_numElem;
				const auto* resDay = result + j * m_numElem;
				int diff = 0;
				T t = -1;
				while (++t < m_numElem && !(diff = (int)p_players[resDayPerm[t]] - resDay[t]));
				if (t < m_numElem) {
					if (diff < 0)
						return false;
					else
						break;
				}
			}

			if (rollBack(p_dayRes, p_dayIsUsed, j, nDays))
				continue;
			//			else

			/*
						if (j > 1) { // Do we need to go to previous day?
							p_dayIsUsed[p_dayRes[--j]] = 0;
							p_dayRes[j--]++;
							continue;
						}
			*/
			if (true || k == nDays) {// Unable to find unused day.
				break;   // there is no previous day for which a different choice can be made 
			}
			// automorphism found
		}

		p_dayIsUsed[iDay] = 0;
		continue;   // temporary
#if 0
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

	if (p_players != buff)
		delete[] p_players;

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