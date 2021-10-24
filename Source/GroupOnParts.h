#pragma once

Class1Def(CGroupHandle)
{
public:
	CGroupHandle()									{}
	~CGroupHandle()									{ delete[] m_pPermut; }
	void InitGroupHandle(S partIdx, S len, size_t groupOrder) {
		// Construct symmetrical group on len elements
		m_nPartIdx = partIdx;
		m_pPermut = new S[(m_nPermLength = len) * (m_nGroupOrder = groupOrder)];
		for (S i = 0; i < len; i++)
			m_pPermut[i] = i;

		bool buffer[32];
		bool* pIndex = len < countof(buffer) ? buffer : new bool[len];
		memset(pIndex, 0, len * sizeof(*m_pPermut));

		auto* currPerm = m_pPermut;
		for (S nPerm = 1; nPerm < groupOrder; nPerm++) {
			memcpy(currPerm + len, currPerm, len * sizeof(*m_pPermut));
			auto currPntr = (currPerm += len) + len;
			S minIdx = len;
			while (true) {
				pIndex[*--currPntr] = true; // mark as "non-used"
				if (minIdx > *currPntr)
					minIdx = *currPntr;
				if (*(currPntr - 1) > *currPntr)
					continue;

				// Find smallest "non-used" which is bigger, than current
				auto elemToReplace = len - (currPntr - currPerm);
				auto i = *--currPntr;
				while (!pIndex[++i]);
				pIndex[*currPntr] = true;
				if (minIdx > *currPntr)
					minIdx = *currPntr;

				pIndex[*currPntr = i] = false;
				minIdx--;
				while (true) {
					if (pIndex[++minIdx]) {
						pIndex[*++currPntr = minIdx] = false;
						if (!--elemToReplace)
							break;
					}
				}
				break;
			}
		}

		if (pIndex != buffer)
			delete[] pIndex;
	}
	CK inline const auto partIdx() const					{ return m_nPartIdx; }
	CK inline const auto groupOrder() const					{ return m_nGroupOrder; }
	CK inline const auto permLength() const					{ return m_nPermLength; }
	CK inline const auto* getPermutation(size_t idx) const	{ return m_pPermut + idx * permLength(); }
private:
	S m_nPartIdx;
	S m_nPermLength;
	size_t m_nGroupOrder;
	S* m_pPermut = NULL;
};

Class1Def(CGroupOnParts)
{
public:
	CGroupOnParts(const void* pOwner, const CVector<S>&lenghts, S minRow) :
		m_pOwner(pOwner), m_nNumGroups(lenghts.GetSize() / 3), m_nMinRowNumb(minRow) {
		m_pGroupHandles = new CGroupHandle<S>[numGroups()];
		const auto lenMax = lenghts.GetSize();
		for (size_t i = 0; i < lenMax; i += 3)
			m_pGroupHandles[i / 3].InitGroupHandle(lenghts.GetAt(i), lenghts.GetAt(i + 1), lenghts.GetAt(i + 2));
	}
	~CGroupOnParts()									{ delete[] m_pGroupHandles; }
	inline const void* owner() const					{ return m_pOwner; }
	CK inline const auto numGroups() const				{ return m_nNumGroups; }
	CK inline const auto* groupHandle(size_t idx) const { return m_pGroupHandles + idx; }
	CK inline auto useGroupOnParts(S nRow) const		{ return nRow > m_nMinRowNumb; }
	CK inline auto getStartingRowNumb() const			{ return m_nMinRowNumb; }
private:
	const void* m_pOwner;
	CGroupHandle<S>* m_pGroupHandles;
	const size_t m_nNumGroups;
	const S m_nMinRowNumb;     // minimal number of rows on which group will be used
};

