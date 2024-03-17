#pragma once

#define IDX_MAX				(ELEMENT_MAX - 1)

template<typename T>
inline void revert(T* perm, T j, T i) {
	while (++i < --j) perm[i] ^= (perm[j] ^= (perm[i] ^= perm[j]));
}

template<typename T>
void Update_Orbits(const T* permut, T lenPerm, T* pOrb, T idx) {
	// Update orbits of elements
	do {
		auto i = *(pOrb + idx);
		auto j = *(pOrb + permut[idx]);
		if (j == i)
			continue;

		if (j < i) {
			i ^= j;
			i ^= j ^= i;
		}

		for (auto k = lenPerm; k--;) {
			if (*(pOrb + k) == j)
				*(pOrb + k) = i;
		}
	} while (++idx < lenPerm);
}

template<typename T>
class CGroupOrder {
public:
	CC inline auto groupOrder() const { return this ? m_nGroupOrder : 1; }
	CC inline void setGroupOrder(size_t val) { m_nGroupOrder = val; }
protected:
	void UpdateOrbits(const T* permut, const T lenPerm, T* pOrb, bool rowPermut, bool calcGroupOrder=false)
	{
		const T idx = udpdateStabLength(permut, lenPerm, pOrb, calcGroupOrder, rowPermut);
		Update_Orbits(permut, lenPerm, pOrb, idx);
	}

	CC void updateGroupOrder(const T numRow, const T * pOrb) {
		size_t len = 1;
		auto idx = stabilizerLengthAut();
		while (++idx < numRow) {
			if (*(pOrb + idx) == stabilizerLengthAut())
				len++;
		}

		setGroupOrder(len * groupOrder());
	}

	CC inline void setStabilizerLength(T len) { m_nStabLength = len; }
	CC inline auto stabilizerLength() const { return m_nStabLength; }
	CC inline void setStabilizerLengthAut(T l) { m_nStabLengthAut = l; }
	CC inline auto stabilizerLengthAut() const { return m_nStabLengthAut; }
	CC T next_permutation(T* perm, const T* pOrbits, T nRow, T idx = ELEMENT_MAX, T lenStab = 0) {
            // Function generates next permutation among those which stabilize first lenStab elements
            // We are using the algorithm from http://nayuki.eigenstate.org/res/next-lexicographical-permutation-algorithm/nextperm.java
            // taking into account that we don't need the the permutations which are equivalent with respect to already found orbits of the 
            // automorphis group acting on the matrix's rows.
            //
            // For instance, we found the automorphism
            //  (0, 1, 2, ..., i-1, pi, ...) and pi != i
            // it means that the elements pi and i are in the same orbits with respect to any stabilizer of the length < i.
            // It also means that on any place j 
            //  (0, 1, 2, ... , j-1,  j, ...) we don't need to try 
            //  (0, 1, 2, ... , j-1, pi, ...) when we already checked the permutation
            //  (0, 1, 2, ... , j-1,  i, ...)

            // Find non-increasing suffix
            T temp, i, j;

            // Check if the algorithm, used immediately after 
            // some non-trivial automorphism was found
            if (idx == IDX_MAX && perm[stabilizerLength()] == nRow - 1)
                idx = ELEMENT_MAX;

            if (idx == IDX_MAX) {
                // Firts call after some automorphism was found
                temp = perm[idx = i = stabilizerLength()];
                for (j = nRow; --j > temp;)
                    perm[j] = j;

                for (auto k = j++; k-- > i;)
                    perm[k + 1] = k;
            }
            else {
                if (idx >= IDX_MAX) {
                    j = i = nRow;
                    while (--i > 0 && perm[i - 1] >= perm[i]);

                    if (i == lenStab)
                        return ELEMENT_MAX;

                    // Find successor to pivot
                    temp = perm[--i];
                    while (perm[--j] <= temp);
                }
                else {
                    temp = perm[j = i = idx];
                    while (++j < nRow && perm[j] <= temp);
                    if (j >= nRow) {
                        revert(perm, nRow, i);
                        return next_permutation(perm, pOrbits, nRow);
                    }
                }
            }

            if (stabilizerLength() == i) {
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
                    while (k < nRow && *(pOrbits + perm[k]) != perm[k])
                        k++;

                    if (k != j) {
                        flag = k == nRow;
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
                    j = idx >= IDX_MAX ? nRow - 1 : i;
                    temp = perm[--i];
                    setStabilizerLength(i);
                }
            }

            perm[i] = perm[j];
            perm[j] = temp;
            if (idx >= IDX_MAX) {
                if (stabilizerLength() > i)
                    setStabilizerLength(i);

                revert(perm, nRow, i);
            }

            return i;
	}

private:
	CC T udpdateStabLength(const T* permut, T lenPerm, const T* pOrb, bool calcOrder, bool rowPermut) {
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

	T m_nStabLength;
	T m_nStabLengthAut;
	size_t m_nGroupOrder;
};
