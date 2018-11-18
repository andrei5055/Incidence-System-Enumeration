
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
	if (!CheckOrbits(permStorage(), getRowOrbits(0)))
		return false;

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
	CCanonicityChecker canonChecker(transpMatr.rowNumb(), nCols);
	CanonizeByColumns(&transpMatr, NULL, &canonChecker);

	// Transposed matrix also should be transitive on the rows
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
	auto nRow = this->currentRowNumb() - 1;
	if (!nRow)
		return true;

	const auto k = this->getInSys()->GetK();
	const auto r = this->getInSys()->GetR();
	const auto nLambd = designParams()->lambdaB().size();
	const auto len = r + k + nLambd + 2 * nLambdas();
	T bufIdx[64];
	T *pIdx = len > countof(bufIdx) ? new T[len] : bufIdx;

	uchar blockFlg[64];
	uchar *pBlockFlags = matrix()->colNumb() > countof(blockFlg)? new uchar [matrix()->colNumb()] : blockFlg;

	auto pIntersection = rowIntersections();
	if (!pIntersection) {
		const auto retVal = CheckTransitivityOnConstructedBlocks(nRow + 1, k, r, pIdx + r, pBlockFlags);
		if (pIdx != bufIdx)
			delete[] pIdx;

		return retVal;
	}

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

	const auto nCol = colNumb();
	const auto &lambdaSet = designParams()->lambda();
	const auto *pCurrRow = matrix()->GetDataPntr();
	auto step = rowNumb();
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

	const bool retBal = (++nRow < k) || CheckConstructedBlocks(nRow, k, pIdx + r);
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
bool CIG_Enumerator<T>::CheckOrbits(const CPermutStorage<T> *permRowStorage, T *pRowOrbits) const
{
	// Last block containing first element was just constructed
	// We can check the Aut(M) to see if it is transitive on first k blocks
	if (permRowStorage->isEmpty())
		return false;

	const int from = pRowOrbits ? 1 : 0;
	T *pColOrbit = getColOrbits(0);
	permRowStorage->CreateOrbits(permColStorage(), this->matrix(), pRowOrbits, pColOrbit, from);
	for (auto j = this->getInSys()->GetR(); j--;) {
		if (*(pColOrbit + j))
			return false;
	}

	return true;
}

template<class T>
void CIG_Enumerator<T>::FindAllElementsOfBlock(T nRow, size_t nColCurr, int j, T *pElementNumb) const
{
	// Making the array of indices of all elements which belong to current block
	// NOTE: the i-th element alwais belong to the block # nColCurr
	T i;
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
}

template<class T>
bool CIG_Enumerator<T>::CheckConstructedBlocks(T nRow, T k, T *pElementNumb)
{
	// We will check the parameters b(i) for all elements, which belong to the constructed block
	// They should be as defined in designParam::lambdaB
	const auto *pColOrbitIni = *(colOrbitsIni() + nRow);
	const CColOrbit<T> *pColOrbit = this->colOrbit(nRow);

	auto lastBlockIdx = this->getInSys()->GetR() - 1;
	const auto nLambd = designParams()->lambdaB().size();
	auto lambdaBCurrRow = pElementNumb + k;

	while (pColOrbit) {
		if (pColOrbit->columnWeight() == k) {
			// Define the current column number
			const size_t nColCurr = ((char *)pColOrbit - (char *)pColOrbitIni) / colOrbitLen();
			FindAllElementsOfBlock(nRow, nColCurr, k, pElementNumb);

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

			if (nColCurr == lastBlockIdx) {
				if (permRowStorage()) {
					if (!CheckOrbits(permRowStorage()))
						return false;

					setNumRows(nRow);
				}
			}
		}
		else
			lastBlockIdx = 0;  // To do not call of CheckOrbits

		pColOrbit = pColOrbit->next();
	}

	return true;
}

template<class T>
bool CIG_Enumerator<T>::CheckTransitivityOnConstructedBlocks(T nRow, T k, T r, T *pElementNumb, uchar *pBlockFlags)
{
	// We cannot make this check only when
	//   (a) we checked canonicity before
	//   (b) two elements may not have common blocks
	if (!permRowStorage() || designParams()->lambda()[0])
		return true;

	// We will check the parameters b(i) for all elements, which belong to the constructed block
	// They should be as defined in designParam::lambdaB
	const auto *pColOrbitIni = *(colOrbitsIni() + nRow);
	const auto *pColOrbit = this->colOrbit(nRow);
	const int colNumb = matrix()->colNumb();

	C_InSys<T> matr;  // Incidence System, which will be canonize
	T *pMatr = NULL;
	uchar *pElemFlags;

	bool exitFlag = false;
	pBlockFlags[0] = 0;
	while (pColOrbit) {
		// Define the current column number
		const auto nColCurr = ((char *)pColOrbit - (char *)pColOrbitIni) / colOrbitLen();
		const auto lenOrb = pColOrbit->length();
		const auto n = nColCurr + lenOrb;
		if (n < r)
			return true;	// Not all blocks containing first element are constructed yet

		if (pColOrbit->columnWeight() == k) {
			if (n == r) {
				if (!CheckOrbits(permRowStorage()))
					return false;

				setNumRows(nRow);
				return true;
			}

			// All blocks containing first element are constructed
			// If the same is TRUE for some other element, 
			// then corresponding submatrices should be identical.

			if (!pBlockFlags[0]) { // The flags for the block were not initialized yet				
				memset(pBlockFlags, 1, r);                  // First r blocks are constructed
				memset(pBlockFlags + r, 0, colNumb - r);	// the others were not checked yet
			}

			if (!pBlockFlags[nColCurr]) {
				for (auto i = lenOrb; i--;)
					pBlockFlags[nColCurr + i] = 1;	// the blocks of current orbits are constructed
			}

			// Construct the array of all elements containing current block
			FindAllElementsOfBlock(nRow, nColCurr, k, pElementNumb);
			bool flag = true;

			// Loop over all these elements
			int i = 0;
			for (; i < k; ++i) {
				// Check if all blocks of current element are constructed
				const auto currElem = pElementNumb[i];
				const auto *pCurrElem = matrix()->GetRow(currElem);

				// When we are here, all r blocks associated with the first elements 
				// are constructed. Therefore, there is no need to test them.  
				int j = colNumb;
				while (--j >= r) {
					if (!pCurrElem[j] || pBlockFlags[j] == 1)
						continue;

					if (pBlockFlags[j] > 1)
						break;   // The block was previously checked and it is NOT yet constructed

					// Check if block # j is constructed
					int kTmp = k;
					auto idx = nRow;
					const auto *pRow = matrix()->GetRow(idx-1) + j + colNumb;
					while (idx--) {
						if (!*(pRow -= colNumb)) {
							if (kTmp > idx) // Is it still possible to reach 0?
								break;		// No - this block is not constructed yet

							continue;
						}

						if (!--kTmp)
							break;
					}

					if (kTmp) {
						pBlockFlags[j] = 2;
						break;
					}

					pBlockFlags[j] = 1;
				}

				if (j >= r)
					continue;   

				// All blocks associated with the current element are constructed
				// Construct matrix with all elements associated with these blocks
				if (!pMatr) {
					assert(k < numRows() && numRows() < matrix()->rowNumb());

					matr.Init(numRows(), colNumb);
					pMatr = new T[colNumb * numRows()];
					if (!elementFlags()) {
						setElementFlags(new uchar[matrix()->rowNumb()]); // allocate maximal amount of flags
						setPermCols(new T[colNumb]);
					}

					pElemFlags = elementFlags();
				}

				T *pCurrRow = pMatr - colNumb;
				if (flag) {
					// Copying the rows, corresponding to k elements associated with the current block (we can do it once)
				    flag = false;   // in order not to do it next time
					
					// Create permutation of columns which corresponds to the maximal representation of the row of currElem 
					int idx1 = 0;
					int idx0 = k;
					T *pPermCols = permCols();
					auto * const pntr = matrix()->GetRow(currElem);
					for (j = 0; j < colNumb; ++j) {
						if (*(pntr + j)) {
							*(pPermCols + idx1++) = j;
							if (idx1 == k) {
								while (++j < colNumb)
									*(pPermCols + idx0++) = j;

								break;
							}
						}
						else
							*(pPermCols + idx0++) = j;
					}

					// Mark all elements as 'unused'
					memset(pElemFlags, 0, nRow * sizeof(*pElemFlags));
					for (j = 0; j < k; ++j) {
						pElemFlags[pElementNumb[j]] = 1;	// element is used
						memcpy(pCurrRow += colNumb, matrix()->GetRow(pElementNumb[j]), colNumb * sizeof(*pMatr));
					}
				}
				else
					pCurrRow += k * colNumb;

				// Copying all remaining rows
				auto numCopyed = k;
				for (j = nRow; j--;) {
					if (pElemFlags[j])
						continue;   // corresponding row is already there

					// To be chosen, the element should be in one of the r block, containing currElem
					// Therefore when rowIntersections() is NOT defined
					const auto *pRow = matrix()->GetRow(j);
					if (!rowIntersections()) {
						// Loop over the all columns
						auto col = colNumb; 
						while (col-- && (!pCurrElem[col] || !pRow[col]));

						if (col < 0)
							continue;
					}
					else {
						// Otherwise, corresponding value of rowIntersections() array should not be zero
						int sIdx, bIdx;
						if (j < currElem) {
							sIdx = j;
							bIdx = currElem;
						}
						else {
							sIdx = currElem;
							bIdx = j;
						}

						const int n = rowNumb() * sIdx - sIdx * (sIdx + 3) / 2 + bIdx - 1;
						if (!*(rowIntersections() + n)) // On zero's index we do have intersection by 0 blocks
							continue;                   // there is no common blocks with the currElem
					}

					if (numCopyed++ >= numRows())
						break; // the matrix has more row than expected

					memcpy(pCurrRow += colNumb, pRow, colNumb * sizeof(pMatr[0]));
				}

				// Check that constructed matrix has the expected number of rows
				if (j >= 0 || numCopyed < numRows())
					break;    

				// Matrix constructed, let's find its canonical representation
				matr.AssignData(pMatr);
				CCanonicityChecker canonChecker(matr.rowNumb(), colNumb);
				CanonizeByColumns(&matr, NULL, &canonChecker, true);

				if (memcmp(matrix()->GetDataPntr(), matr.GetDataPntr(), matr.lenData()))
					break;        // There is no isomorphism we expected to see 
			}

			if (i < k) {
				exitFlag = true;
				break;
			}
		} else {
			if (n == r)
			   return true;
        }

		if (exitFlag)
			break;

		pColOrbit = pColOrbit->next();
	}

	delete[] pMatr;
	return !exitFlag;
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
}