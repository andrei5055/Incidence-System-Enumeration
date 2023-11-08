#include "PermutStorage.h"
#include "matrix.h"

template class CPermutStorage<TDATA_TYPES>;

PermutStorage(void)::orderPermutations(size_t *pPermPerm)
{
	const auto nPerm = numPerm();
	assert(nPerm > 0);
	pPermPerm[0] = 0;
	pPermPerm[1] = 1;
	if (nPerm <= 2)
		return;

	for (size_t i = nPerm; i-- > 2;)
		pPermPerm[i] = i;

	for (size_t i = 0; i < nPerm; i++) {
		auto jBest = UINT64_MAX;
        size_t idx = pPermPerm[i];
		const auto *pFirst = permutMemory() + lenPerm() * idx;
		for (auto j = i + 1; j < nPerm; j++) {
			auto jdx = pPermPerm[j];
			const auto *pSecnd = permutMemory() + lenPerm() * jdx;
			for (int k = 0; true; k++) {
				const int diff = (int)*(pFirst + k) - (int)*(pSecnd + k);
				if (!diff)
					continue;

				if (diff > 0) {
					jBest = j;
					idx = jdx;
					pFirst = pSecnd;
				}

				break;
			}
		}

		if (idx == pPermPerm[i])
			continue;

		pPermPerm[jBest] = pPermPerm[i];
		pPermPerm[i] = idx;
	}
}

#if OUT_PERMUTATION
FILE* CPermutStorage<TDATA_TYPES>::m_pFile;

PermutStorage(void)::printPerm(const T* pPerm, bool savePerm, int add, T permLen) const {
	if (ccc < CHECK_CCC)
		return;

	if (!pPerm)
		return;

#if OUT_PERMUTATION == 1		// output of the permutations only for completely constructed matrix
	if (!savePerm)
		return;
#endif
	if (!permLen)
		permLen = lenPerm();

	char buf[256];
	const auto lenBuffer = 3 * permLen + 9;
	char* pBuffer = lenBuffer < sizeof(buf) ? buf : new char[lenBuffer];
	char* pBuf = pBuffer;
	if (this)
		pBuf += SNPRINTF(pBuf, lenBuffer - (pBuf - pBuffer), "%4d: ", (m_cntr += add));
	else
		pBuf += SNPRINTF(pBuf, lenBuffer - (pBuf - pBuffer), "====> ");

	for (T i = 0; i < permLen; i++)
		pBuf += SNPRINTF(pBuf, lenBuffer - (pBuf - pBuffer), "%3d", pPerm[i]);

	strcpy_s(pBuf, lenBuffer - (pBuf - pBuffer), "\n");
	outString(pBuffer, outFile());
	if (pBuffer != buf)
		delete[] pBuf;
}
#endif


