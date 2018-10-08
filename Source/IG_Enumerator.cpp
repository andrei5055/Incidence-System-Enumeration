
#include "IG_Enumerator.h"

template class CIG_Enumerator<MATRIX_ELEMENT_TYPE>;

template<class T>
bool CIG_Enumerator<T>::fileExists(const char *path, bool file) const {
	const bool retVal = CEnumerator<T>::fileExists(path, file);
	if (!file || !retVal)
		return retVal;

	// The answer is fake. When following statement is false, the file exist,
	// but we don't need the caller know that, becase for Inconsistent graphs
	// all outputs for same order graphs will be in the same file
	return designParams()->logFile != std::string(path);
}

template<class T>
bool CIG_Enumerator<T>::createNewFile(const char *fName) const {
	if (!fName)
		return firstPath();

	if (designParams()->logFile == fName)
		return false;

	designParams()->logFile = fName;
	return true;
}

template<class T>
bool CIG_Enumerator<T>::TestFeatures(CEnumInfo<T> *pEnumInfo, const CMatrixData<T> *pMatrix, int *pMatrFlags, CEnumerator<T> *pEnum) const
{
	if (!CPBIBD_Enumerator<T>::TestFeatures(pEnumInfo, pMatrix, pMatrFlags))
		return false;

	if (!this->groupIsTransitive())
		return false;

	// Check the transitivity of the stabilizer of the  
	// first element on the blocks to which it belongs
	permStorage()->CreateOrbits(permColStorage(), pMatrix, getRowOrbits(0));
	const T *pColOrbit = getColOrbits(0);
	const auto k = this->getInSys()->GetK();
	for (auto j = k; j--;) {
		if (*(pColOrbit + j))
			return false;
	}

	*pMatrFlags |= t_trahsitiveGroup;

	const auto nCols = pMatrix->colNumb();
	const auto nRows = pMatrix->rowNumb();
	const auto mult = nCols / nRows;

	const designParam *pDesignParam = designParams();
	if (pDesignParam->lambda().size() > 2) {
		auto degrees(pDesignParam->lambdaA());
		int r = pDesignParam->r;
		for (int i = 1; i < nRows; ++i) {
			const auto pRow = pMatrix->GetRow(i);
			int lambdaCur = 0;
			for (int j = 0; j < r; ++j)
				lambdaCur += *(pRow + j);

			const auto idx = findLambda(lambdaCur);
			if (idx < 0 || --degrees[idx] < 0)
				return false;
		}
	}

	// Need to check that transposed matrix is not isomorphic to constructed one
	C_InSys<T> transpMatr;
	transpMatr.InitTransposed(pMatrix, mult);
	CCanonicityChecker canonChecker(transpMatr.rowNumb(), nCols, 2);
	CanonizeByColumns(&transpMatr, NULL, &canonChecker);

	// Transposed matrix also should be transitive the the rows
	if (!canonChecker.groupIsTransitive())
		return false;

	auto pntr = pMatrix->GetDataPntr();
	bool retVal = false;
	if (mult > 1) {
		auto pntrTr = transpMatr.GetDataPntr() - nCols;
		for (int i = 0; !retVal && i < nRows; i++, pntr += nCols) {
			for (int j = 0; j < mult; ++j) {
				if (memcmp(pntr, pntrTr += nCols, nCols * sizeof(T))) {
					retVal = true;
					break;
				}
			}
		}
	}
	else
		retVal = memcmp(pntr, transpMatr.GetDataPntr(), transpMatr.lenData()) != 0;

	if (retVal && checkProperty(t_printTransposedMatrix)) {
		// Printing of the transposed matrix
		transpMatr.printOut(this->outFile(), nCols, pEnumInfo->constrCanonical() + 1, &canonChecker);
	}

	return retVal;
}

template<class T>
bool CIG_Enumerator<T>::makeFileName(char *buffer, size_t lenBuffer, const char *ext) const
{
	const auto dirLength = this->getDirectory(buffer, lenBuffer);
	const auto pParam = designParams();
	const auto nVertex = this->rowNumb() * pParam->r / pParam->k;
	SNPRINTF(buffer + dirLength, lenBuffer - dirLength, "%ss_of_order_%d%s", getObjName(), nVertex, ext ? ext : FILE_NAME(""));
	return true;
}

template<class T>
bool CIG_Enumerator<T>::prepareToFindRowSolution() {
	// In that function we
	// (a) calculate the intersection of last constructed row with all previously constructed rows
	auto pIntersection = rowIntersections();
	if (!pIntersection)
		return true;

	auto nRow = this->currentRowNumb() - 1;
	if (!nRow)
		return true;

	const auto v = rowNumb();
	const auto nCol = colNumb();
	const auto len = v * (v - 1) / 2;
	const auto r = this->getInSys()->GetR();
	T bufIdx[32];
	T *pIdx = r > countof(bufIdx) ? new T[r] : bufIdx;

	// Define the set of indices of the blocks, containing last element
	auto i = r;
	const auto *pCurrLastRow = matrix()->GetRow(nRow);
	auto j = colNumb();
	while (true) {
		if (!*(pCurrLastRow + --j))
			continue;

		*(pIdx + --i) = j;
		if (!i)
			break;
	}

	const auto &lambdaSet = designParams()->lambda();
	const auto *pCurrRow = matrix()->GetDataPntr();
	auto step = v;
	pIntersection += nRow - step;

	for (int i = 0; i < nRow; ++i) {
		// Define the intersection of the current and last rows of the matrix
		size_t lambda = 0;
		for (auto j = r; j--;) {
			if (*(pCurrRow + *(pIdx + j)))
				++lambda;
		}

		// Define and save the index of the lambda we just found
		int j = 0;
		while (lambdaSet[j++] != lambda);

		*(pIntersection += --step) = j - 1;
		pCurrRow += nCol;
	}

	const auto k = this->getInSys()->GetK();
	const bool retBal = (++nRow < k) || CheckConstructedBlocks(nRow, k);
	if (pIdx != bufIdx)
		delete[] pIdx;

	return retBal;
}

size_t calcSum(const std::vector<int> &lambdaA, const std::vector<int> &lambda) {
	size_t sum = 0;
	for (auto i = lambdaA.size(); i--;)
		sum += lambdaA[i] * lambda[i] * lambda[i];

	return sum;
}

template<class T>
size_t calcSum(const T *lambdaA, int len, int mult) {
	size_t sum = 0;
	for (auto i = mult; i <= len; i += mult)
		sum += lambdaA[i/mult] * i * i;

	return sum;
}

template<class T>
static int lexCompare(const std::vector<int> &lamA, const std::vector<int> &lam, const T *lambdaA, int idx, int mult)
{
	idx /= mult;
	for (auto j = lam.size(); j--;) {
		const auto val = lam[j];
		while (idx * mult > val && !lambdaA[idx])
			idx--;

		if (idx * mult > val)
			return -1;  //second is greater

		if (idx * mult < val || lambdaA[idx] < lamA[j])
			return 1;	// first is greater

		if (lambdaA[idx] > lamA[j])
			return -1;  //second is greater

		idx--;
	}

	return 0;			// both are the same
}

template<class T>
bool CIG_Enumerator<T>::CheckConstructedBlocks(T nRow, T k) const
{
	// We will check the parameters b(i) for all elements, which belong to the constructed block
	// They should be as defined in designParam::lambdaB
	const auto *pColOrbitIni = *(colOrbitsIni() + nRow);
	const CColOrbit<T> *pColOrbit = this->colOrbit(nRow);

	const auto nLambd = designParams()->lambdaB().size();
	const auto len = k + nLambd + 2 * nLambdas();

	T elementNumb[64] ;
	auto pElementNumb = len > countof(elementNumb) ? new T[len] : elementNumb;
	auto lambdaBCurrRow = pElementNumb + k;

	while (pColOrbit) {
		if (pColOrbit->columnWeight() == k) {
			// Define the current column number
			const size_t nColCurr = ((char *)pColOrbit - (char *)pColOrbitIni) / colOrbitLen();

			// Making the array of indices of all elements which belong to current block
			T i, j = k;
			// NOTE: the i-th element alwais belong to the block # nColCurr
			pElementNumb[--j] = i = nRow - 1;
			// Looking for remaining elements which also belong to the same block
			while (true) {
				const auto *pRow = matrix()->GetRow(--i) + nColCurr;
				if (*pRow) {
					*(pElementNumb + --j) = i;
					if (!j)
						break;
				}
			}

			// Ordered pair of element numbers 
			T sIdx, bIdx;
			for (T i = 0; i < k; ++i) {
				T *pIndex = &sIdx;
				bIdx = *(pElementNumb + i);

				memcpy(lambdaBCurrRow, lambdaBSrc(), nLambd * sizeof(*lambdaBCurrRow));
				for (T j = 0; j < k; ++j) {
					if (j == i) {
						sIdx = *(pElementNumb+i);
						pIndex = &bIdx;
						continue;
					}

					*pIndex = *(pElementNumb + j);

					// Both bIdx-th and sIdx-th elements belong to the current block and sIdx
					// The index of number of thier common blocks is the n-th element of rowIntersections()
					//        n = v*(v-1)/2 - (k-sIdx)*(k-sIdx-1) / 2 + bIdx - sIdx - 1
					// After simplification   
					const int n = rowNumb() * sIdx - sIdx * (sIdx + 3) / 2 + bIdx - 1;
					const auto idx = *(rowIntersections() + n);
					if (!lambdaBCurrRow[idx]--)
						return false;  // We have exceeded the number of possible intersections for the current idx
				}

				if (!DefineInterstructForBlocks(nColCurr, k, pElementNumb, i, lambdaBCurrRow + nLambd))
					return false;
			}
		}

		pColOrbit = pColOrbit->next();
	}

	return true;
}

template<class T>
bool CIG_Enumerator<T>::DefineInterstructForBlocks(size_t nColCurr, T k, const T *pElementNumb, T i, T *lambdaACurrCol) const
{
	// The structure of intersection of constructed block with the other (even not yet completely
	// constructed) blocks will be checked. It should be exactly the same as for block #0
	const auto nBlocks = colNumb();
	const auto *pRow = matrix()->GetRow(0);
	if (nColCurr || *(pElementNumb + i)) {
		auto lambdaBCurrCol = lambdaACurrCol + nLambdas();
		memset(lambdaACurrCol, 0, lenLambdas());
		T n = -1;
		while (++n < nBlocks) {
			if (n == nColCurr)
				continue;  //skip current block

			bool flag = false;
			const auto *pBlock = pRow + n;
			int lambdaCurr = 0;
			for (T j = 0; j < k; ++j) {
				if (*(pBlock + *(pElementNumb + j) * nBlocks)) {
					++lambdaCurr;
					if (i == j)
						flag = true;
				}
			}

			if (flag)
				++*(lambdaBCurrCol + lambdaCurr);

			++*(lambdaACurrCol + lambdaCurr);
		}

		return !memcmp(lambdaACurrCol, lambdaA(), lenLambdas());
	}

	// The block #0 was just constructed. Define the structure of intersections
	// of this block with the other blocks, containing element #0
	// In canonical matrix these blocks have #'s from 1 <= j < k
	memset(lambdaA(), 0, lenLambdas());
	const auto r = this->getInSys()->GetR();
	T n = 0;
	while (++n < nBlocks) {
		const auto *pBlock = pRow + n;
		// We know everything about first matrix row. Let's start with the second one
		int lambdaCurr = 0;
		for (T j = 1; j < k; ++j) {
			if (*(pBlock += nBlocks))
				++lambdaCurr;
		}

		if (n < r)
			++*(lambdaB() + lambdaCurr + 1);
		else
			++*(lambdaA() + lambdaCurr);
	}

	// Real value lambdaA(i) should be defined AND all values 
	// lambdaA and lambdaB should correspond each other
	for (T j = 1; j <= k; ++j) {
		if (!(*(lambdaA() + j) += *(lambdaB() + j)))
			continue;  // for zeros the following relation always holds

		if (*(lambdaA() + j) * j != k * *(lambdaB() + j))
			return false;
	}

	// According to Theorem 5.4 from Andrei Ivanov's thesis S(i^2 * a(i)) 
	// is a constant which is the same for transposed matrix. Let's check if  
	// it is true for just constructed lambdaA() and designParams()->lambdaA()
	// NOTE: In fact, since we simplified the design of IG with vertex replicas, 
	// we have to use NOT designParams()->lambdaA(), but corresponding initial  
	// array of these parameters which is still kept in
	auto pRowInterStruct = designParams()->InterStruct();
	const auto mult = pRowInterStruct->mult();
#define NEW  1
#if NEW
	const auto pCounterparts = pRowInterStruct->getNext()->Counterparts();
	for (const auto element : *pCounterparts) {
		const auto &lam = element->lambda();
		const auto &lamA = element->lambdaA();
		const auto comp = lexCompare(lamA, lam, lambdaA(), r, mult);
		if (comp < 0)
			return false;
		if (!comp)
			return true;
	}

	return false;
#else
	const auto sumForCols = calcSum(lambdaA(), r, mult);
	while (pRowInterStruct = pRowInterStruct->getNext()) {
		const auto &lam = pRowInterStruct->lambda();
		const auto &lamA = pRowInterStruct->lambdaA();
		if (sumForCols != calcSum(lamA, lam))
			continue;

		// if (lamA, lam) lexicografically greater than lambdaA() 
		if (lexCompare(lamA, lam, lambdaA(), r, mult) < 0)
			return false;

		break;
	}
	return pRowInterStruct != NULL;
#endif
}