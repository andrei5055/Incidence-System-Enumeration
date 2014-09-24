#include "DataTypes.h"
#include "PermutStorage.h"

CPermutStorage::CPermutStorage()
{
	setPermutMemory(NULL);
	setLenMemUsed(0);
	setLenMemMax(0);
}

CPermutStorage::~CPermutStorage()
{
	delete [] permutMemory();
}

size_t *CPermutStorage::allocateMemoryForPermut(size_t lenPermut)
{
	const size_t lenUsed = lenMemUsed();
	const size_t newLength = lenUsed + lenPermut;

	if (lenMemMax() < newLength) {
		const size_t len = 2 * newLength;
		setLenMemMax(len);
		auto *permutMem = new size_t[len];
		if (permutMemory()) {
			memcpy(permutMem, permutMemory(), lenUsed * sizeof(*permutMemory()));
			delete[] permutMemory();
		}

		setPermutMemory(permutMem);
	}

	setLenMemUsed(newLength);
	return permutMemory() + lenUsed;
}

void CPermutStorage::savePermut(const size_t lenPermut, const size_t *perm)
{
	size_t *pMem = allocateMemoryForPermut(lenPermut);

	if (lenPerm() == SIZE_MAX)
		setLenPerm(lenPermut);

	if (!perm) {
		for (size_t i = lenPermut; i--;)
			*(pMem + i) = i;
	} else
		memcpy(pMem, perm, lenPermut * sizeof(*permutMemory()));
}

void CPermutStorage::outputPermutations(FILE *file, size_t lenPerm) const
{
	for (size_t len = 0; len < lenMemUsed(); len += lenPerm)
		outputPerm(file, permutMemory() + len, lenPerm);
}

void CPermutStorage::outputPerm(FILE *file, const size_t *perm, size_t lenPerm) const
{
	size_t len = lenPerm < 10 ? 2 : lenPerm < 100 ? 3 : lenPerm < 1000 ? 4 : 5;
	char format[8];
	SPRINTF(format, "%%%lud", len);
	len *= lenPerm;
	char *pBuffer = new char[len += 2];
	char *pBuf = pBuffer;
	for (size_t i = 0; i < lenPerm; i++)
		pBuf += sprintf_s(pBuf, len - (pBuf - pBuffer), format, *(perm + i));

	strcpy_s(pBuf, len - (pBuf - pBuffer), "\n");
	outString(pBuffer, file);
	delete[] pBuffer;
}

void CPermutStorage::orderPermutations(size_t *pPermPerm)
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
		size_t jBest = -1, idx = pPermPerm[i];
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

void CPermutStorage::adjustGenerators(int *pIdx, size_t lenIdx)
{
	// Adjustement of the generators of the automorphism group on columns according to 
	// non-unforcible group of columns defined by indeces in pIdx 
	size_t i, iMax = numPerm();
	int buff[256];
	int *pFrom = iMax <= countof(buff) ? buff : new int[iMax];
	int j = 0;
	for (i = 0; i < lenIdx; i++)
		pFrom[pIdx[i]] = j++;

	auto *pTo = permutMemory();
	for (i = 0; i < iMax; i++, pTo += lenIdx) {
		auto *pColPerm = getPermutByIndex(i);
		for (size_t j = 0; j < lenIdx; j++)
			pTo[j] = pFrom[pColPerm[pIdx[j]]];
	}

	if (pFrom != buff)
		delete pFrom;

	setLenPerm(lenIdx);
	setLenMemUsed(lenIdx * iMax);
}

size_t CPermutStorage::constructGroup()
{
	// Algorithm from http://www.sipria.ru/pdf/dt24114.pdf is implemented here.
	size_t buffer[16];
	const auto nPerm = numPerm();
	size_t *nextStart = nPerm < countof(buffer) ? buffer : new size_t[nPerm + 1];

	// Loop for all permutations except the trivial one
	auto lastIdx = nPerm;
	size_t *pTmpPerm = NULL;
	for (size_t i = 1; i < nPerm; i++) {
		// Construct all degrees of current permutations
		auto tmpPermIdx = i;
		auto lastIdxTmp = lastIdx;
		while (!(pTmpPerm = multiplyPermutations(tmpPermIdx, i, pTmpPerm, &lastIdxTmp)))
			tmpPermIdx = lastIdxTmp - 1;

		// Save index for next round of multiplication
		const auto jMax = lastIdx;
		nextStart[i] = lastIdx = lastIdxTmp;

		// Find left and right multiplication of current permutation and  all permultations from previous round
		for (size_t j = i+1; j < jMax; j++) {
			pTmpPerm = multiplyPermutations(j, i, pTmpPerm, &lastIdx);
			pTmpPerm = multiplyPermutations(i, j, pTmpPerm, &lastIdx);
		}
	}

	size_t lastPrev = SIZE_MAX;
	while (lastPrev != lastIdx) {
		lastPrev = lastIdx;
		for (size_t i = 1; i < nPerm; i++) {
			for (auto j = nextStart[i]; j < lastIdx; j++) {
				pTmpPerm = multiplyPermutations(j, i, pTmpPerm, &lastIdx);
				pTmpPerm = multiplyPermutations(i, j, pTmpPerm, &lastIdx);
			}
		}
	}

	if (pTmpPerm)
		deallocateMemoryForPermut();

	if (nextStart != buffer)
		delete [] nextStart;

	return numPerm();
}

size_t *CPermutStorage::multiplyPermutations(size_t firstPermIdx, size_t secondPermIdx, size_t *pMultRes, size_t *pToIdx)
{
	if (!pMultRes)
		pMultRes = allocateMemoryForPermut(lenPerm());

	// Since the permutations could be moved in allocateMemoryForPermut,
	// we need to access the permuts by their indices here
	const auto *pFirst = getPermutByIndex(firstPermIdx);
	const auto *pSecond = getPermutByIndex(secondPermIdx);

	// Do multiplication of two permutations here:
	for (auto k = lenPerm(); k--;)
		*(pMultRes + k) = *(pFirst + *(pSecond + k));

	// Compare with the previously constructed permutations
	for (auto i = *pToIdx; i--;) {
		if (!memcmp(getPermutByIndex(i), pMultRes, lenPermByte()))
			return pMultRes;   // this memory will be reused
	}

	++*pToIdx;
	return NULL;
}

void CPermutStorage::multiplyPermutations(size_t currPermIdx, size_t fromIdx, size_t toIdx, size_t permOrder, size_t lastIdx, size_t *pPermPerm)
{
	for (size_t j = fromIdx; j < toIdx; j++) {
		const auto permIdx = pPermPerm? pPermPerm[j] : j;
		multiplyPermutations(permIdx, currPermIdx);

		for (size_t p = 0; p < permOrder; p++)
			multiplyPermutations(permIdx, lastIdx + p);
	}
}

size_t CPermutStorage::findSolutionIndex(const VECTOR_ELEMENT_TYPE *pFirst, size_t idx, VECTOR_ELEMENT_TYPE *pMem, size_t *pCanonIdx, int &nCanon)
{
	const VECTOR_ELEMENT_TYPE *pCurrSolution = pFirst + lenPerm() * idx;
	VECTOR_ELEMENT_TYPE *pCanonical = (VECTOR_ELEMENT_TYPE *)pCurrSolution;
	const auto len = lenPerm() * sizeof(*pMem);

	// Loop for all permutations, except the trivial one
	const auto *pPermut = permutMemory();
	const auto *pPermutEnd = pPermut + lenMemUsed();
	while ((pPermut += lenPerm()) < pPermutEnd) {
		for (auto i = lenPerm(); i--;)
			pMem[i] = pCurrSolution[*(pPermut + i)];

		if (memcmp(pMem, pCanonical, len) <= 0)
			continue;

		// Copy current solution as canonical
		if (pCanonical != pCurrSolution) {
			VECTOR_ELEMENT_TYPE *pTmp = pMem;
			pMem = pCanonical;
			pCanonical = pTmp;
		} else {
			pCanonical = pMem;
			pMem += lenPerm();
		}
	}

	if (pCanonical == pCurrSolution)
		return pCanonIdx[nCanon++] = idx; // Current solution is the canonical one

	// Find corresponding canonical solution
	int i = nCanon;
	while (i-- && memcmp(pFirst + lenPerm() * pCanonIdx[i], pCanonical, len));

	if (i < 0) {
		// It could happen when we are using reordering of solutions by the automorphism group
		// BIBD(12, 44, 11, 3, 2) provides a lot's of such examples
		return SIZE_MAX;
	}

	return pCanonIdx[i];
}
