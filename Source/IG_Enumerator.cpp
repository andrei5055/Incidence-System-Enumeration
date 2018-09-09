
#include "IG_Enumerator.h"

template class CIG_Enumerator<MATRIX_ELEMENT_TYPE>;

#define PRINT_DEBUG  0
#if PRINT_DEBUG
static bool fff(int nCols, int rowNumb, uchar *pMatr, FILE *file) {
	for (int j = 0; j < nCols; ++j) {
		int sum = 0;
		for (int row = 0; row < rowNumb; ++row)
			sum += *(pMatr + row * nCols + j);

		if (sum != 6) {
			fprintf(file, "Problem with column # %d\n", j);
			return true;
		}
	}
	return false;
}
#endif

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

	*pMatrFlags |= t_trahsitiveGroup;

	const auto nCols = pMatrix->colNumb();
	const auto nRows = pMatrix->rowNumb();
	const auto mult = nCols / nRows;

	const designRaram *pDesignParam = designParams();
	if (pDesignParam->lambda.size() > 2) {
		auto degrees(pDesignParam->lambdaA);
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
	// CMatrixData<T>  transpMatr;
	CMatrix<T> transpMatr(nCols, nCols, 1);
	transpMatr.InitTransposed(pMatrix, mult);
	const auto rowNumb = transpMatr.rowNumb();
	CMatrixCol<T> matrCol(&transpMatr);

	// For now we will consider only binary incidence systems
	// This part of the program will be a bit more complicated for general case
	assert(matrCol.rankMatr() <= 2);

	matrCol.initiateColOrbits(rowNumb, this->IS_enumerator());
	auto pColIdxMem = new T[nCols];

	const auto colOrbLen = matrCol.colOrbitLen();
	bool flag = false;

	CCanonicityChecker canonChecker(rowNumb, nCols, 2);
#if PRINT_DEBUG
	FOPEN(file, "C:\\Users\\andreii\\Calc\\fff.txt", "w");
	transpMatr.printOut(file, rowNumb, 99, NULL);
	int nMatr = 1;
#endif
	auto pMatr = matrCol.matrix()->GetDataPntr();
	T *pTmp = NULL;		// Memory for reordering the rows and columns of the matrix
						// If needed, it will be allocated. 
	T i = 0;
	while (true) {
		for (auto j = nCols; j--;)
			pColIdxMem[j] = j;

		CColOrbit<T> *pColOrbitNext = matrCol.colOrbits()[i];
		while (i < rowNumb) {
			T *pBeg = matrCol.matrix()->GetRow(i);
			CColOrbit<T> *pColOrbit = pColOrbitNext;
			auto pColOrbNext = pColOrbitNext = matrCol.colOrbits()[++i];
			auto pColIdx = pColIdxMem;

			// Loop over the orbits of columns
			while (pColOrbit) {
				// Each column orbits could be splited into sub-orbits
				// For binary incidence systems the algorythm is simple (and it is implemented here)
				// For general case we need to re-order columns of current orbits according to their
				// incidence with current element (corresponding to the current row of the matrix)

				const auto len = pColOrbit->length();
				auto idx = 0;
				auto idxLast = len;
				while (true) {
					// Find first zero from left
					while (idx < idxLast && pBeg[pColIdx[idx]])
						idx++;

					if (idx == idxLast)
						break;

					// Find first non-zero form right
					while (idx < --idxLast && !pBeg[pColIdx[idxLast]]);

					if (idx == idxLast)
						break;

					const auto i1 = pColIdx[idx++];
					const auto i2 = pColIdx[idxLast] - i1;

					auto pTmp = pBeg + i1;
					*pTmp = 1;
					*(pTmp + i2) = 0;

					int j = i;
					while (++j <= rowNumb) {
						auto tmp = *(pTmp += nCols);
						*pTmp = *(pTmp + i2);
						*(pTmp + i2) = tmp;
					}
				}

				pColOrbit = pColOrbit->next();
				if (i < rowNumb) {
					CColOrbit<T> *pNext = pColOrbit ? (CColOrbit<T> *)((char *)pColOrbNext + colOrbLen * len) : NULL;

					// Save the column's orbit information
					if (idx == 0 || idx == len) // Orbit was not splitted
						pColOrbNext->Init(len, pNext);
					else {
						CColOrbit<T> *pNxt = (CColOrbit<T> *)((char *)pColOrbNext + colOrbLen * idx);
						pColOrbNext->Init(idx, pNxt);
						pNxt->Init(len - idx, pNext);
					}

					pColOrbNext = pNext;
					pColIdx += len;
				}
#if PRINT_DEBUG
				if (fff(nCols, rowNumb, pMatr, file))
					break;
#endif
			}
		}

#if PRINT_DEBUG
		transpMatr.printOut(file, rowNumb, -nMatr, NULL);
		if (fff(nCols, rowNumb, pMatr, file))
			break;
#endif
		if (canonChecker.TestCanonicity(rowNumb, &matrCol, t_saveRowPermutations))
			break;  // Matrix is canonized

					// Reorder the rows and columns of the matrix according to  
					// permutations which were found during canonicity testing

					// Define the row number, where non-canonicity was noticed
		i = 0;
		const auto pPermRow = canonChecker.permRow();
		while (pPermRow[i] == i)
			i++;

		if (!pTmp)
			pTmp = new T[nCols * rowNumb];

		// Copy last (rowNum - i) rows of matrix into temporary buffer
		memcpy(pTmp, matrCol.matrix()->GetRow(i), nCols * (rowNumb - i));

#if PRINT_DEBUG
		fprintf(file, "i = %d\n", i);

		for (int k = 0; k < rowNumb; k++)
			fprintf(file, "%2d ", k);

		fprintf(file, "\n");

		auto ppp = pPermRow;
		for (int j = 0; j < 2; ++j) {
			for (int k = 0; k < rowNumb; k++)
				fprintf(file, "%2d ", ppp[k]);

			fprintf(file, "\n");
			ppp = canonChecker.permCol();
		}
#endif
		// Reorder the row and the columns:
		const auto pPermCol = canonChecker.permCol();
		for (T row = i; row < rowNumb; ++row) {
			auto pMatrTo = pMatr + row * nCols;
			const auto pMatrFrom = pTmp + (pPermRow[row] - i) * nCols;
			for (T j = 0; j < nCols; ++j)
				pMatrTo[j] = pMatrFrom[pPermCol[j]];
		}

#if PRINT_DEBUG
		transpMatr.printOut(file, rowNumb, nMatr++, &canonChecker);
		if (fff(nCols, rowNumb, pMatr, file))
			break;
#endif
	}
#if PRINT_DEBUG
	if (file)
		fclose(file);
#endif
	delete[] pColIdxMem;
	delete[] pTmp;

	// Transposed matrix also shoule be transitive the the rows
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

	if (retVal) {
		// Printing of the transposed matrix
		transpMatr.printOut(this->outFile(), nCols, 99, &canonChecker);
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
//	return true;
	// In that function we
	// (a) calculate the intersection of last constructed row with all previusly constructed rows
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

	const auto &lambdaSet = designParams()->lambda;
	const auto *pCurrRow = matrix()->GetDataPntr();
	auto step = v;
	pIntersection += nRow - step;

	for (int i = 0; i < nRow; ++i) {
		// Define the intersection of the current and last rows of the matrix;
		size_t lambda = 0;
		for (auto j = r; j--;) {
			if (*(pCurrRow + *(pIdx + j)))
				++lambda;
		}

		// Define and save the index of the lambda we just found
		int j = 0;
		while (lambdaSet[j] != lambda)
			j++;

		*(pIntersection += --step) = j;

		pCurrRow += nCol;
	}

	const auto k = this->getInSys()->GetK();
	const bool retBal = (++nRow < k) || CheckConstructedBlocks(nRow, k);
	if (pIdx != bufIdx)
		delete[] pIdx;

	return retBal;
}

template<class T>
bool CIG_Enumerator<T>::CheckConstructedBlocks(T nRow, T k)
{
	// The structure of intersection of constructed block with the other (even not yet completely
	// constructed) blocks will be checked. It should be exactly the same as for block #0
	// We will also check the parameters b(i) for all elements, which belong to the constructed block
	// They chould be as defined in designRaram::lambdaB

	const auto *pColOrbitIni = *(colOrbitsIni() + nRow);
	CColOrbit<T> *pColOrbit = this->colOrbit(nRow);
	T elementNumb[64];
	assert(k <= countof(elementNumb));
	while (pColOrbit) {
		if (pColOrbit->columnWeight() == k) {
			// Define the current column number
			const size_t nColCurr = ((char *)pColOrbit - (char *)pColOrbitIni) / colOrbitLen();

			// Making the array of indices of all elements which belong to current block
			T i, j = k;
			// NOTE: the i-th element alwais belong to the block # nColCurr
			elementNumb[--j] = i = nRow - 1;
			// Looking for remaining elements which also belong to the same block
			while (true) {
				const auto *pRow = matrix()->GetRow(--i) + nColCurr;
				if (*pRow) {
					elementNumb[--j] = i;
					if (!j)
						break;
				}
			}

			// Ordered pair of element numbers 
			T sIdx, bIdx;
			for (T i = 0; i < k; ++i) {
				T *pIndex = &sIdx;
				bIdx = elementNumb[i];
				auto lambdaB = designParams()->lambdaB;
				for (T j = 0; j < k; ++j) {
					if (j == i) {
						sIdx = elementNumb[i];
						pIndex = &bIdx;
						continue;
					}

					*pIndex = elementNumb[j];

					// Both bIdx-th and sIdx-th elements belong to the current block
					// The index of number of thier common blocks is the n-th element of rowIntersections()
					//        n = v*(v-1)/2 - (k-sIdx)*(k-sIdx-1) / 2 + bIdx - sIdx - 1
					// After simplification   
					const int n = rowNumb() * sIdx - sIdx * (sIdx + 3) / 2 + bIdx - 1;
					const auto idx = *(rowIntersections() + n);
					if (!lambdaB[idx]--)
						return false;  // We have exceeded the number of possible intersections for the current idx
				}
			}
		}

		pColOrbit = pColOrbit->next();
	}

	return true;
}
