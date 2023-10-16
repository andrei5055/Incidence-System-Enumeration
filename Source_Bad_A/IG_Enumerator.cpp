
#include "IG_Enumerator.h"

template class CIG_Enumerator<TDATA_TYPES>;
bool isPrime(int k);

FClass2(CIG_Enumerator)::CIG_Enumerator(const InSysPntr pBIBD, const designParam *pParam, unsigned int enumFlags, bool firstPath, int treadIdx, uint nCanonChecker) :
	Class2(CPBIBD_Enumerator)(pBIBD, enumFlags, treadIdx, nCanonChecker), m_firstPath(firstPath) {
	const auto inSys = this->getInSys();
	const auto lambdaSet = inSys->GetNumSet(t_lSet);
	auto nLambd = lambdaSet->GetSize();

	m_pRowIntersections = nLambd > 2 || lambdaSet->GetAt(0) || lambdaSet->GetAt(1) != 1 ?
		new S[lenRowIntersection(this->rowNumb())] : NULL;

	// Allocate memory to store current lambdaA, lambdaB for blocks
	// This set of lambda's is not known, but it cannot have more than k + 1 elements
	const auto k = inSys->GetK();
	const auto r = getR();
	const auto len = k + 1;
	m_pLambda[0] = new S[2 * len + nLambd];
	m_pLambda[2] = (m_pLambda[1] = m_pLambda[0] + nLambd) + len;

	// Converting the values lamdaB into T format
	const auto &coeffB = pParam->lambdaB();
	for (auto i = nLambd; i--;)
		*(m_pLambda[0] + i) = coeffB[i];

	m_pElements = new CIncidenceStorage<T, S>(inSys->rowNumb(), r);
	m_pBlocks = new CIncidenceStorage<T, S>(inSys->colNumb(), k);
	setElementFlags(NULL);
	setPermCols(NULL);
	setNumRows(0);


	// Allocate memory to store the numbers of units in three areas, if needed
	if (--nLambd <= 2 && pParam->lambda()[nLambd] == 2 && k == r) {
		setAreaWeight(new uint[3]);
		memset(areaWeight(), 0, 3 * sizeof(*areaWeight()));
		const bool flag = /*false;*/ !pParam->InterStruct()->Counterparts();
		setStrongCheckUnitsInThreeAreas(flag && k % 3);
	} else
		setAreaWeight(NULL);
}

FClass2(CIG_Enumerator)::~CIG_Enumerator() {
	delete[] rowIntersections();
	delete[] lambdaBSrc();
	delete[] elementFlags();
	delete[] permCols();
	delete[] areaWeight();
}

FClass2(CIG_Enumerator, void)::CloneMasterInfo(const EnumeratorPntr p, size_t nRow) {
	const auto pMaster = static_cast<const Class2(CIG_Enumerator) *>(p);
	copyRowIntersection(pMaster->rowIntersections());

	// Copying the numbers of units which are in three areas, if we use this check
	if (pMaster->areaWeight())
		memcpy(areaWeight(), pMaster->areaWeight(), 3 * sizeof(*areaWeight()));

	// No need to clone remaining information from master if it is not yet constructed
	if (nRow < this->getInSys()->GetK())
		return;

	memcpy(lambdaA(), pMaster->lambdaA(), lenLambdas());
	setNumRows(pMaster->numRows());
}

FClass2(CIG_Enumerator, void)::copyRowIntersection(const S *pntr) {
	if (rowIntersections())
		memcpy(rowIntersections(), pntr, lenRowIntersection(this->rowNumb()) * sizeof(*pntr));
}

FClass2(CIG_Enumerator, bool)::fileExists(const char *path, bool file) const {
	const bool retVal = Class2(CEnumerator)::fileExists(path, file);
	if (!file || !retVal)
		return retVal;

	// The answer is fake. When following statement is false, the file exists,
	// but we don't need the caller knows that, because for semi-symmetric graphs
	// all outputs for same order graphs will be in the same file
	return this->designParams()->logFile != std::string(path);
}

FClass2(CIG_Enumerator, bool)::createNewFile(const char *fName) const {
	if (!fName)
		return firstPath();

	if (this->designParams()->logFile == fName)
		return false;

	this->designParams()->logFile = fName;
	return true;
}

FClass2(CIG_Enumerator, bool)::TestFeatures(EnumInfoPntr pEnumInfo, const MatrixDataPntr pMatrix, int *pMatrFlags, const EnumeratorPntr pEnum) const
{
	if (!Class2(CPBIBD_Enumerator)::TestFeatures(pEnumInfo, pMatrix, pMatrFlags))
		return false;

	if (!this->groupIsTransitive())
		return false;

	// Check the transitivity of the stabilizer of the  
	// first element on the blocks to which it belongs
	if (!CheckOrbits(this->permStorage(), this->getRowOrbits(0)))
		return false;

	*pMatrFlags |= t_transitiveGroup;

	const auto nCols = pMatrix->colNumb();
	const auto nRows = pMatrix->rowNumb();
	const auto mult = nCols / nRows;

	const designParam *pDesignParam = this->designParams();
	if (pDesignParam->lambda().size() > 2) {
		auto degrees(pDesignParam->lambdaA());
		const int r = pDesignParam->r;
		for (auto i = nRows; --i;) {
			const auto pRow = pMatrix->GetRow(i);
			int lambdaCur = 0;
			for (int j = 0; j < r; ++j)
				lambdaCur += *(pRow + j);

			const auto idx = this->findLambda(lambdaCur);
			if (idx < 0 || --degrees[idx] < 0)
				return false;
		}
	}

	// Need to check that transposed matrix is not isomorphic to constructed one
	Class2(C_InSys) transpMatr;
	transpMatr.InitTransposed(pMatrix, mult);
	Class2(CCanonicityChecker) canonChecker(transpMatr.rowNumb(), nCols);
	this->CanonizeByColumns(&transpMatr, NULL, &canonChecker);

	// Transposed matrix also should be transitive on the rows
	if (!canonChecker.groupIsTransitive())
		return false;

	auto pntr = pMatrix->GetDataPntr();
	bool retVal = false;
	if (mult > 1) {
		const auto *pntrTr = transpMatr.GetDataPntr() - nCols;
		for (auto i = nRows; i--; pntr += nCols) {
			for (auto j = mult; j--;) {
				if (memcmp(pntr, pntrTr += nCols, nCols * sizeof(T))) {
					retVal = true;
					i = 0;
					break;
				}
			}
		}
	}
	else
		retVal = memcmp(pntr, transpMatr.GetDataPntr(), transpMatr.lenData()) != 0;

	if (retVal && this->checkProperty(t_printTransposedMatrix)) {
		// Printing of the transposed matrix
		transpMatr.printOut(this->outFile(), nCols, pEnumInfo->numMatrOfType(t_design_type::t_canonical) + 1, &canonChecker);
	}

	return retVal;
}

FClass2(CIG_Enumerator, bool)::makeFileName(char *buffer, size_t lenBuffer, const char *ext) const
{
	const auto dirLength = this->getDirectory(buffer, lenBuffer);
	const auto pParam = this->designParams();
	const auto nVertex = this->rowNumb() * pParam->r / pParam->k;
	SNPRINTF(buffer + dirLength, lenBuffer - dirLength, "%ss_of_order_%d%s", getObjName(), nVertex, ext ? ext : FILE_NAME(""));
	return true;
}

FClass2(CIG_Enumerator, bool)::prepareToFindRowSolution() {
	// In that function we
	// (a) check the numbers of units in the areas
	// (b) calculate the intersection of last constructed row with all previously constructed rows
	S nRow = this->currentRowNumb() - 1;
	if (!nRow)
		return true;

	const auto r = this->getR();
	const auto k = this->getInSys()->GetK();
	if (checkUnitsInThreeAreas() && nRow > 1 && nRow < 2 * (k - 1)) {
		if (!checkThreeAreasUnits(r, k, nRow))
			return false;
	}

	const auto pMatrix = this->matrix();
	const auto nLambd = this->designParams()->lambdaB().size();
	const auto len = r + k + nLambd + 2 * nLambdas();
	S bufIdx[64] = {};
	S *pIdx = len > countof(bufIdx) ? new S[len] : bufIdx;

	uchar blockFlg[64] = {};
	uchar *pBlockFlags = pMatrix->colNumb() > countof(blockFlg)? new uchar [pMatrix->colNumb()] : blockFlg;

	auto pIntersection = rowIntersections();
	const auto retVal = CheckTransitivityOnConstructedBlocks(nRow + 1, k, r, pIdx + r, pBlockFlags);
	if (!pIntersection || !retVal) {
		if (pIdx != bufIdx)
			delete[] pIdx;
#if 0
		if (!retVal) {
			static int cntr = 0;
			const auto pMatrix = static_cast<const CMatrix<T> *>(this->matrix());
			OUT_MATRIX(pMatrix, outFile(), nRow + 1, cntr);
			if (++cntr >= 5)
				exit(0);
		}
#endif
		return retVal;
	}

	// Define the set of indices of the blocks, containing last element
	auto i = r;
	const auto *pCurrLastRow = pMatrix->GetRow(nRow);
	auto j = this->colNumb();
	while (true) {
		if (!*(pCurrLastRow + --j))
			continue;

		*(pIdx + --i) = j;
		if (!i)
			break;
	}

	const auto nCol = this->colNumb();
	const auto &lambdaSet = this->designParams()->lambda();
	const auto *pCurrRow = pMatrix->GetDataPntr();
	auto step = this->rowNumb();
	pIntersection += nRow - step;

	for (uint i = 0; i < nRow; ++i) {
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

FClass2(CIG_Enumerator, bool)::checkThreeAreasUnits(S r, S k, S nRow) const
{
	const auto *pRow = this->matrix()->GetRow(nRow);
	if (nRow == 2 && *(pRow + 1))
		return false;

	const auto lambdaB = this->designParams()->lambdaB();
	const auto b2 = lambdaB[lambdaB.size() - 1];
	const S valMax = b2 - 1;
	S fromIdx, nCol, nColLast, lim = 2 * (k - 1);
	bool flag;
	if (!strongCheckUnitsInThreeAreas()) {
		// Theorem 1: When b(2) != 0, b(3) == ... == b(k) == 0, following three (k-2)x(k-2) areas contain exactly b(2)-1 units each, 
		// no more than one unit in each row or column.
		//
		//    1 1 1 1 . . . . . 1| 0 0 . . . 0|0 0 ...
		//    1 1 0 0 . . . . . 0| 1 1 . . . 1|0 0 ...
		//    -------------------|------------|---------
		//    1 0 1 0 . . . . . 0|            |
		//    . 0 0 . . . . . . 0|            |
		//    . 0 . . . . . . . 0|            |
		//    1 0 . . . 1 . . . 0|   area 1   |
		//    1 0 0 . . 0 0 . . 0|            |
		//    . . . . . . . . . .|            |
		//    1 0 0 . . 0 0 . . 0|            |
		//    -------------------|------------|---------
		//    0 1|               |            |
		//    0 1|               |            |
		//    . .|    area 2     |   area 3   |
		//    . .|               |            |
		//    0 1|               |            |
		//    ---|---------------|------------|----------
		//    0 0 1 . . . . . . . . . . . . . . 
		//    
		//    Prof: Trivial.
		//
		const int step = r - 2;
		if (nRow < k) {
			fromIdx = 0;
			nColLast = r;
			flag = nRow == k - 1;
		}
		else {
			fromIdx = 1;
			nColLast = 2;
			flag = nRow == lim - 1;
		}

		bool retVal = true;
		const auto toIdx = 2 * fromIdx + 1;
		do {
			nColLast = (nCol = nColLast) + step;
			while (nCol < nColLast) {
				if (*(pRow + nCol))
					break;
				else
					nCol++;
			}

			if (nCol < nColLast && ++*(areaWeight() + fromIdx) > valMax ||
				flag && *(areaWeight() + fromIdx) != valMax)
				retVal = false;
		} while (++fromIdx < toIdx);

		return retVal;
	} 

	// Theorem 2: When 
	//    (a) b(2) != 0, b(3) == ... == b(k) == 0,
	//    (b) k % 3 != 0
	// Then:
	//    (a) following three (k-2)x(k-2) areas contain exactly b(2)-1 units each, 
	// no more than one unit in each row or column.
	//    (b) in area1 & area3 the units are in upper left part of diagonal.
	//
	//    1 1 1 1 . . . . . 1| 0 0 . . . 0|0 0 ...
	//    1 1 0 0 . . . . . 0| 1 1 . . . 1|0 0 ...
	//    -----------|-------|------------|---------
	//    1 0 1 0 . .|. . . 0|0 . . . .  0|
	//    . 0 0 . . .|. . . 0|            |
	//    . 0 . . . .|. . . 0|   area A   |
	//    1 0 . . . 1|. . . 0|0 . . . .  0|
	//    -----------|-------|------------|---------
	//    1 0 0 . . 0|0 . . 0|1 0         |
	//    . . . . . .|. . . .|0 1 area 1  |
	//    1 0 0 . . 0|0 . . 0|            |
	//    ---|-------|-------|------------|---------
	//    0 1|0 . . 0|1 0 . 0|            |
	//    0 1|       |0 1    |            |
	//    . .|area B |area 2 |   area 3   |
	//    . .|       |       |            |
	//    0 1|0 . . 0|       |            |
	//    ---|-------|-------|------------|----------
	//    0 0 1 . . . . . . . . . . . . . . 
	//    
	//    Prof: Trivial.
	//

	int step;
	if (nRow < k) {
		fromIdx = 0;
		nColLast = b2 + 1;
		flag = nRow == k - 1;
		step = 2 * r - 2 - nColLast;
	}
	else {
		fromIdx = 1;
		nColLast = 2;
		flag = nRow == lim - 1;
		step = r - 2;
	}

	bool retVal = true;
	const auto toIdx = 2 * fromIdx + 1;
	const auto colNumber = this->matrix()->colNumb();
	do {
		nColLast = (nCol = nColLast) + step;
		while (nCol < nColLast) {
			if (*(pRow + nCol))
				break;
			else
				nCol++;
		}

		if (nCol < nColLast) {
			if (nRow <= b2 || nCol <= b2 || nRow < k && nCol < r)
				return false;  // unit in area A or B

			if (2 * b2 <= nRow && nRow < k || nRow >= k && 2 * b2 <= nCol && nCol < r)
				return false;  // more units in area 1 or 2, than it's allowed

			if (nRow >= k && nCol >= r) {
				++*(areaWeight() + 2);
				if (flag && *(areaWeight() + 2) != valMax)
					return false;
			}

			if (*(pRow + nCol - colNumber))
				return false;  // two units in two consecutive rows
		}
		else {
			// all zeros in current row of the area
			if (b2 < nRow && nRow < 2 * b2 ||
				k <= nRow && nRow < k + b2 - 1 && nCol < r)
				return false;
		}

	} while (++fromIdx < toIdx);

	return retVal;
}

FClass2(CIG_Enumerator, void)::reset(T nRow, bool resetSolutions) {
	Class2(CEnumerator)::reset(nRow);
	if (nRow-- == numRows()) 
		setNumRows(0);

	const auto k = this->getInSys()->GetK();
	if (!checkUnitsInThreeAreas() || nRow <= 1 || nRow >= 2 * (k - 1))
		return;

	int idxLast, idx;
	const auto r = this->getR();
	const auto *pRow = this->matrix()->GetRow(nRow);
	if (!strongCheckUnitsInThreeAreas()) {
		int fromIdx, toIdx;
		if (nRow < k) {
			if (!*areaWeight())
				return;  // nothing to subtract

			fromIdx = 0;
			toIdx = 1;
			idxLast = r;
		}
		else {
			if (*(areaWeight() + 1)) {
				fromIdx = 1;
				idxLast = 2;
				toIdx = *(areaWeight() + 2) ? 3 : 2;
			}
			else if (*(areaWeight() + 2)) {
				fromIdx = 2;
				toIdx = 3;
				idxLast = r;
			}
			else
				return;  // nothing to subtract
		}

		do {
			idxLast = (idx = idxLast) + r - 2;
			while (idx < idxLast) {
				if (*(pRow + idx++)) {
					--*(areaWeight() + fromIdx);
					break;
				}
			}

		} while (++fromIdx < toIdx);
	}
	else {
		// When we are here, only the units from area 3 should be extracted
		auto *pArea3 = areaWeight() + 2;
		if (!*pArea3)
			return; // nothing to extract from area 3

		idx = r;
		idxLast = idx + r - 2;
		while (idx < idxLast) {
			if (*(pRow + idx++)) {
				--*pArea3;
				return;
			}
		}
	}
}

template<class T>
static int lexCompare(const std::vector<uint> &lamA, const std::vector<uint> &lam, const T *lambdaA, int idx, uint mult)
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

FClass2(CIG_Enumerator, bool)::CheckOrbits(const PermutStoragePntr permRowStorage, S *pRowOrbits) const
{
	// Last block containing first element was just constructed
	// We can check the Aut(M) to see if it is transitive on first k blocks
	if (permRowStorage->isEmpty())
		return false;

	const int from = pRowOrbits ? 1 : 0;
	auto *pColOrbit = this->getColOrbits(0);
	permRowStorage->CreateOrbits(this->permColStorage(), this->matrix(), pRowOrbits, pColOrbit, from);
	for (auto j = this->getR(); j--;) {
		if (*(pColOrbit + j))
			return false;
	}

	return true;
}

FClass2(CIG_Enumerator, void)::FindAllElementsOfBlock(S nRow, size_t nColCurr, int j, S *pElementNumb) const
{
	// Making the array of indices of all elements which belong to current block
	// NOTE: the i-th element alwais belong to the block # nColCurr
	auto i = pElementNumb[--j] = nRow - 1;
	// Looking for remaining elements which also belong to the same block
	const auto pMatrix = this->matrix();
	while (true) {
		const auto *pRow = pMatrix->GetRow(--i) + nColCurr;
		if (*pRow) {
			*(pElementNumb + --j) = i;
			if (!j)
				break;
		}
	}
}

FClass2(CIG_Enumerator, bool)::CheckConstructedBlocks(S nRow, S k, S *pElementNumb)
{
	// We will check the parameters b(i) for all elements, which belong to the constructed block
	// They should be as defined in designParam::lambdaB
	const auto *pColOrbitIni = *(this->colOrbitsIni() + nRow);
	const auto *pColOrbit = this->colOrbit(nRow);

	auto lastBlockIdx = this->getR() - 1;
	const auto nLambd = this->designParams()->lambdaB().size();
	auto lambdaBCurrRow = pElementNumb + k;

	const auto colOrbitLen = this->colOrbitLen();
	const auto rowNumb = this->rowNumb();
	const auto permRowStorage = this->permRowStorage();
	while (pColOrbit) {
		if (pColOrbit->columnWeight() == k) {
			// Define the current column number
			const size_t nColCurr = ((char *)pColOrbit - (char *)pColOrbitIni) / colOrbitLen;
			FindAllElementsOfBlock(nRow, nColCurr, k, pElementNumb);

			// Ordered pair of element numbers 
			S sIdx, bIdx;
			for (S i = 0; i < k; ++i) {
				auto *pIndex = &sIdx;
				bIdx = *(pElementNumb + i);

				memcpy(lambdaBCurrRow, lambdaBSrc(), nLambd * sizeof(*lambdaBCurrRow));
				for (S j = 0; j < k; ++j) {
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
					const int n = rowNumb * sIdx - sIdx * (sIdx + 3) / 2 + bIdx - 1;
					const auto idx = *(rowIntersections() + n);
					if (!lambdaBCurrRow[idx]--)
						return false;  // We have exceeded the number of possible intersections for the current idx
				}

				if (!DefineInterstructForBlocks(nColCurr, k, pElementNumb, i, lambdaBCurrRow + nLambd))
					return false;
			}

			if (nColCurr == lastBlockIdx) {
				if (permRowStorage) {
					if (!CheckOrbits(permRowStorage))
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

FClass2(CIG_Enumerator, bool)::CheckTransitivityOnConstructedBlocks(S nRow, S k, S r, S *pElementNumb, uchar *pBlockFlags)
{
	// We cannot make this check only when
	//   (a) we checked canonicity before
	//   (b) two elements may not have common blocks
	const auto permRowStorage = this->permRowStorage();
	if (!permRowStorage || this->designParams()->lambda()[0])
		return true;

	// We will check the parameters b(i) for all elements, which belong to the constructed block
	// They should be as defined in designParam::lambdaB
	const auto *pColOrbitIni = *(this->colOrbitsIni() + nRow);
	const auto *pColOrbit = this->colOrbit(nRow);
	const auto pMatrix = this->matrix();
	const auto colNumb = pMatrix->colNumb();
	const auto rowNumb = pMatrix->rowNumb();

	Class2(C_InSys) matr;  // Incidence System, which will be canonize
	T *pMatr = NULL;
	uchar *pElemFlags;

	bool exitFlag = false;
	pBlockFlags[0] = 0;
	while (pColOrbit) {
		// Define the current column number
		const auto nColCurr = ((char *)pColOrbit - (char *)pColOrbitIni) / this->colOrbitLen();
		const auto lenOrb = pColOrbit->length();
		const auto n = nColCurr + lenOrb;
		if (n < r)
			return true;	// Not all blocks containing first element are constructed yet

		if (pColOrbit->columnWeight() == k) {
			if (n == r) {
				if (!CheckOrbits(permRowStorage))
					return false;

				setNumRows(nRow);
				return true;
			}
			else {
				// Some block with # n > r is constructed
				if (!numRows())	  // Check if all r first blocks are constructed
					return true;  // They are not
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
			S i = 0;
			for (; i < k; ++i) {
				// Check if all blocks of current element are constructed
				const auto currElem = pElementNumb[i];
				const auto *pCurrElem = pMatrix->GetRow(currElem);

				// When we are here, all r blocks associated with the first elements 
				// are constructed. Therefore, there is no need to test them.  
				S j = colNumb;
				while (--j >= r) {
					if (!pCurrElem[j] || pBlockFlags[j] == 1)
						continue;

					if (pBlockFlags[j] > 1)
						break;   // The block was previously checked and it is NOT yet constructed

					// Check if block # j is constructed
					auto kTmp = k;
					auto idx = nRow;
					const auto *pRow = pMatrix->GetRow(idx-1) + j + colNumb;
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
					assert(k < numRows() && numRows() < rowNumb);

					matr.Init(numRows(), colNumb);
					pMatr = new T[colNumb * numRows()];
					if (!elementFlags()) {
						setElementFlags(new uchar[rowNumb]); // allocate maximal amount of flags
						setPermCols(new S[colNumb]);
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
					auto *pPermCols = permCols();
					auto * const pntr = pMatrix->GetRow(currElem);
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
						memcpy(pCurrRow += colNumb, pMatrix->GetRow(pElementNumb[j]), colNumb * sizeof(*pMatr));
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
					const auto *pRow = pMatrix->GetRow(j);
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

						const int n = rowNumb * sIdx - sIdx * (sIdx + 3) / 2 + bIdx - 1;
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
				Class2(CCanonicityChecker) canonChecker(matr.rowNumb(), colNumb);
				this->CanonizeByColumns(&matr, NULL, &canonChecker, true);

				if (memcmp(pMatrix->GetDataPntr(), matr.GetDataPntr(), matr.lenData()))
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

FClass2(CIG_Enumerator, bool)::DefineInterstructForBlocks(size_t nColCurr, S k, const S *pElementNumb, S i, S *lambdaACurrCol) const
{
	// The structure of intersection of constructed block with the other (even not yet completely
	// constructed) blocks will be checked. It should be exactly the same as for block #0
	const auto nBlocks = this->colNumb();
	const auto *pRow = this->matrix()->GetRow(0);
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
	const auto r = this->getR();
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
	auto pRowInterStruct = this->designParams()->InterStruct();
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