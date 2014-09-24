#include <stdio.h>

class CPermutStorage
{
public:
	CPermutStorage();
	~CPermutStorage();
	void outputPermutations(FILE *file, size_t len) const;
	void adjustGenerators(int *pIdx, size_t lenIdx);
	size_t constructGroup();
	size_t findSolutionIndex(const VECTOR_ELEMENT_TYPE *pFirst, size_t idx, VECTOR_ELEMENT_TYPE *pMem, size_t *pCanonIdx, int &nCanon);
	inline size_t lenPerm() const					{ return m_nLenPerm; }
	inline void initPermutStorage()					{ setLenMemUsed(0); setLenPerm(SIZE_MAX); }
	void savePermut(const size_t lenPerm, const size_t *perm = NULL);
	void outputPerm(FILE *file, const size_t *perm, size_t lenPerm) const;
	inline size_t *getPermutByIndex(size_t i) const	{ return permutMemory() + i * lenPerm(); }
	inline size_t numPerm() const					{ return lenMemUsed() / lenPerm(); }
protected:
	inline void savePermut(const size_t *pPerm)		{ savePermut(lenPerm(), pPerm); }
private:
	inline size_t lenPermByte() const				{ return m_nLenPermByte;  }
	inline void setPermutMemory(size_t *pntr)       { m_pPermutMem = pntr; }
	inline size_t *permutMemory() const             { return m_pPermutMem; }
	inline size_t lenMemMax() const                 { return m_nLenMax; }
	inline size_t lenMemUsed() const                { return m_nLenUsed; }
	inline void setLenMemMax(size_t val)            { m_nLenMax = val; }
	inline void setLenMemUsed(size_t val)           { m_nLenUsed = val; }
	inline void setLenPerm(size_t val)				{ m_nLenPermByte = (m_nLenPerm = val) * sizeof(m_pPermutMem[0]); }
	inline void deallocateMemoryForPermut()			{ m_nLenUsed -= lenPerm(); }
	size_t *allocateMemoryForPermut(size_t lenPermut);
	void orderPermutations(size_t *pPermPerm);
	size_t *multiplyPermutations(size_t firstPermIdx, size_t secondPermIdx, size_t *pMultRes = NULL, size_t *pToIdx = NULL);
	void multiplyPermutations(size_t currPermIdx, size_t fromIdx, size_t toIdx, size_t permOrder, size_t lastIdx, size_t *pPermPerm);

	size_t *m_pPermutMem;
	size_t m_nLenMax;
	size_t m_nLenUsed;
	size_t m_nLenPerm;
	size_t m_nLenPermByte;
};