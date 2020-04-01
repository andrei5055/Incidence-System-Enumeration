#include "InsSysEnumerator.h"
#include "RightPartFilter.h"
#include "EnumInfo.h"

#if USE_EXRA_EQUATIONS
#include "EquSystem.h"
#endif

template class C_InSysEnumerator<TDATA_TYPES>;

FClass2(C_InSysEnumerator, void)::CanonizeByColumns(MatrixDataPntr pMatrix, S *pColIdxStorage, CanonicityCheckerPntr pCanonChecker, bool permCol) const
{
	const auto rowNumb = pMatrix->rowNumb();
	const auto nCols = pMatrix->colNumb();
	Class2(CMatrixCol) matrCol(pMatrix);
	// For now we will consider only binary incidence systems
	// This part of the program will be a bit more complicated for general case
	assert(matrCol.rankMatr() <= 2);

	matrCol.initiateColOrbits(rowNumb, 0, this->matrix()->partsInfo(), this->IS_enumerator());
	auto pColIdxMem = pColIdxStorage? pColIdxStorage : new S[nCols];

	const auto colOrbLen = matrCol.colOrbitLen();
	bool flag = false;

	auto pMatr = matrCol.matrix()->GetDataPntr();
	T *pTmp = NULL;		// Memory for reordering the rows and columns of the matrix
						// If needed, it will be allocated. 
	S i = 0;
	while (true) {
		bool colPermFound = false;
		for (auto j = nCols; j--;)
			pColIdxMem[j] = j;

		auto *pColOrbitNext = matrCol.colOrbits()[i];
		while (i < rowNumb) {
			T *pBeg = matrCol.matrix()->GetRow(i);
			auto *pColOrbit = pColOrbitNext;
			auto pColOrbNext = pColOrbitNext = matrCol.colOrbits()[++i];
			auto pColIdx = pColIdxMem;

			// Loop over the orbits of columns
			while (pColOrbit) {
				// Each column orbits could be splited into sub-orbits
				// For binary incidence systems the algorithm is simple (and it is implemented here)
				// For general case we need to re-order columns of current orbits according to their
				// incidence with current element (corresponding to the current row of the matrix)

				const auto len = pColOrbit->length();
				uint32_t idx = 0;
				auto idxLast = len;
				while (true) {
					// Find first zero from left
					while (idx < idxLast && pBeg[pColIdx[idx]])
						idx++;

					if (idx == idxLast)
						break;

					// Find first non-zero from right
					while (idx < --idxLast && !pBeg[pColIdx[idxLast]]);

					if (idx == idxLast)
						break;

					const auto i1 = pColIdx[idx];
					pColIdx[idx++] = pColIdx[idxLast];
					pColIdx[idxLast] = i1;
					colPermFound = true;
				}

				pColOrbit = pColOrbit->next();
				if (i < rowNumb) {
					CColOrbit<S> *pNext = pColOrbit ? (CColOrbit<S> *)((char *)pColOrbNext + colOrbLen * len) : NULL;

					// Save the column's orbit information
					if (idx == 0 || idx == len) // Orbit was not split
						pColOrbNext->Init(len, pNext);
					else {
						CColOrbit<S> *pNxt = (CColOrbit<S> *)((char *)pColOrbNext + colOrbLen * idx);
						pColOrbNext->Init(idx, pNxt);
						pNxt->Init(len - idx, pNext);
					}

					pColOrbNext = pNext;
				}

				pColIdx += len;
			}
		}

		if (!pCanonChecker)
			break;

		if (permCol && colPermFound) {
			// Permut of columns is needed
			S jMin = 0;
			while (jMin < nCols && pColIdxMem[jMin] == jMin)
				jMin++;

			if (jMin < nCols) {
				if (!pTmp)
					pTmp = new T[nCols * rowNumb];

				const auto len = (nCols - jMin) * sizeof(T);
				for (S i = 0; i < rowNumb; ++i) {
					auto *pRow = matrCol.matrix()->GetRow(i);
					memcpy(pTmp + jMin, pRow + jMin, len);
					for (S j = jMin; j < nCols; ++j)
						pRow[j] = pTmp[pColIdxMem[j]];
				}
			}
		}

		if (pCanonChecker->TestCanonicity(rowNumb, &matrCol, t_saveRowPermutations))
			break;  // Matrix is canonized

					// Reorder the rows and columns of the matrix according to  
					// permutations which were found during canonicity testing

					// Define the row number, where non-canonicity was noticed
		i = 0;
		const auto pPermRow = pCanonChecker->permRow();
		while (pPermRow[i] == i)
			i++;

		if (!pTmp)
			pTmp = new T[nCols * rowNumb];

		// Copy last (rowNum - i) rows of matrix into temporary buffer
		memcpy(pTmp, matrCol.matrix()->GetRow(i), nCols * (rowNumb - i));

		// Reorder the row and the columns:
		const auto pPermCol = pCanonChecker->permCol();
		for (S row = i; row < rowNumb; ++row) {
			auto pMatrTo = pMatr + row * nCols;
			const auto pMatrFrom = pTmp + (pPermRow[row] - i) * nCols;
			for (S j = 0; j < nCols; ++j)
				pMatrTo[j] = pMatrFrom[pPermCol[j]];
		}
	}

	if (!pColIdxStorage)
		delete[] pColIdxMem;

	delete[] pTmp;
}

FClass2(C_InSysEnumerator, void)::ConstructColumnPermutation(const MatrixDataPntr pMatrix)
{
	Class2(C_InSys) transformedMatr;
	transformedMatr.InitWithPermutedRows(pMatrix, this->permRow(), this->numRow());
	const auto colNumb = pMatrix->colNumb();
	CanonizeByColumns(&transformedMatr, this->permColStorage()->allocateMemoryForPermut(colNumb));
	this->permColStorage()->setLenPerm(colNumb);
}
