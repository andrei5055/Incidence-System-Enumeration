#pragma once

template<typename T>
class CGroupOrder {
public:
	CC inline auto groupOrder() const { return this ? m_nGroupOrder : 1; }
	CC inline void setGroupOrder(size_t val) { m_nGroupOrder = val; }
protected:
	CC void updateGroupOrder(const T numRow, const T * pOrb) {
		size_t len = 1;
		auto idx = stabilizerLengthAut();
		while (++idx < numRow) {
			if (*(pOrb + idx) == stabilizerLengthAut())
				len++;
		}

		setGroupOrder(len * groupOrder());
	}
	CC T udpdateStabLength(const T * permut, T lenPerm, const T * pOrb, bool calcOrder, bool rowPermut) {
		T idx = 0;
		while (idx == permut[idx])
			idx++;

		if (calcOrder) {
			if (rowPermut && stabilizerLengthAut() > idx)
				updateGroupOrder(lenPerm, pOrb);

			setStabilizerLengthAut(idx);
		}

		setStabilizerLength(idx);
		return idx;
	}
	CC inline void setStabilizerLength(T len) { m_nStabLength = len; }
	CC inline auto stabilizerLength() const { return m_nStabLength; }
	CC inline void setStabilizerLengthAut(T l) { m_nStabLengthAut = l; }
	CC inline auto stabilizerLengthAut() const { return m_nStabLengthAut; }
private:
	T m_nStabLength;
	T m_nStabLengthAut;
	size_t m_nGroupOrder;
};