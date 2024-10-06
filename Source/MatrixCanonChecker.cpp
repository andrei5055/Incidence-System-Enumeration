#include "MatrixCanonChecker.h"
#include "RowSolution.h"
#include "designParam.h"

template class CMatrixCanonChecker<TDATA_TYPES>;

FClass2(CMatrixCanonChecker)::~CMatrixCanonChecker() {
	delete enumInfo();
	delete[] commonElemNumber();
	delete[] blockIdx();
	delete[] partIdx();
}

FClass2(CMatrixCanonChecker, ColOrbPntr)::MakeRow(T nRow, const T *pRowSolution, uint clean_flags, T partIdx) const
{
	auto* pRow = this->matrix()->ResetRowPart(nRow, partIdx, clean_flags);

	const auto* pColOrbit = this->colOrbit(nRow, partIdx);
	const auto* pNextRowColOrbit = this->colOrbitIni(nRow + 1, partIdx);
	const auto colOrbLen = this->colOrbitLen();

	const auto rank = CCanonicityChecker<T, S>::rank();
	const auto maxElement = rank - 1;
	const auto* pColOrbitIni = this->colOrbitIni(nRow, partIdx);

	const bool nextColOrbNeeded = clean_flags & t_getNextColOrb;
	ColOrbPntr pNextRowColOrbitNew = NULL;
	ColOrbPntr pColOrbitLast = NULL;
	while (pColOrbit) {
		ColOrbPntr pNewColOrbit = NULL;
		// Define the number of columns to start with
		const auto nColCurr = ((char*)pColOrbit - (char*)pColOrbitIni) / colOrbLen;
		auto lenRemaining = pColOrbit->length();
		auto* pRowCurr = pRow + nColCurr;
		for (auto i = rank; i--;) {
			const auto lenFragm = i ? pRowSolution[maxElement - i] : lenRemaining;
			if (!lenFragm)
				continue;

			if (nextColOrbNeeded) {
				if (!pNewColOrbit) {
					pNewColOrbit = (ColOrbPntr)((char*)pNextRowColOrbit + nColCurr * colOrbLen);
					if (pColOrbitLast)
						pColOrbitLast->setNext(pNewColOrbit);
					else
						pNextRowColOrbitNew = pNewColOrbit;
				}

				pNewColOrbit = (pColOrbitLast = pNewColOrbit)->InitOrbit(lenFragm, colOrbLen, pColOrbit, i);
			}

			if (!i)
				break;

			// Construct corresponding part of matrix's nRow's row
			rowSetFragm(pRowCurr, i, lenFragm);
			if (!(lenRemaining -= lenFragm))
				break;    // remaining parts of the current solution are 0's
		}

		pRowSolution += maxElement;
		pColOrbit = pColOrbit->next();
	}

	auto* ppUnforced = getUnforcedColOrbPntr(partIdx);
	if (ppUnforced) {
		// Set unforced (by non-zero values) elements:
		auto row = firstUnforcedRow();
		ppUnforced += this->shiftToUnforcedOrbit(row);
		for (; row <= nRow; row++) {
			const auto* pColOrbitIni = this->colOrbitIni(row, partIdx);
			for (int i = 1; i < rank; i++) {
				pColOrbit = *(ppUnforced + i);
				while (pColOrbit) {
					const auto nColCurr = ((char*)pColOrbit - (char*)pColOrbitIni) / colOrbLen;
					const auto lenFragm = pColOrbit->length();
					if (lenFragm > 1)
						rowSetFragm(pRow + nColCurr, i, lenFragm);
					else
						pRow[nColCurr] = i;

					pColOrbit = pColOrbit->next();
				}
			}
			ppUnforced += rank;
		}
	}

	if (pColOrbitLast)    // Not the last row of the matrix
		pColOrbitLast->setNext(NULL);

	return pNextRowColOrbitNew;
}

#define PRINT_INTERSECTIONS PRINT_SOLUTIONS
#if PRINT_INTERSECTIONS
static size_t cntr = 0;
#endif

FClass2(CMatrixCanonChecker, void)::ResetBlockIntersections(T nRow, T partIndex) {
	// Reset the intersections of the current block with blocks
	// containing the current element
	const auto b = matrix()->colNumb();
	const auto lenPart = b / classSize();
	auto* pBlockIdx = blockIdx() + lenPart * (--nRow);
	const auto nColAbs = *(pBlockIdx + partIndex);
	auto* pCommonElemNumber = commonElemNumber() + nColAbs * b;
	if (*(partIdx() + nRow) > partIndex)
		*(partIdx() + nRow) = partIndex;

#if PRINT_INTERSECTIONS
#define START_PRINTING 0 //42000000
	char buf[256];
	if (++cntr >= START_PRINTING) {
		sprintf_s(buf, "\nRow %d, part %d: Reset intersections of block %2d with (", nRow + 1, partIndex, nColAbs);
		outString(buf, outFile());
		if (cntr >= 345439) //345411)
			cntr += 0;
	}
#endif
	for (auto i = partIndex; i--;) {
		auto* pntr = pCommonElemNumber + *(pBlockIdx + i);
#if PRINT_INTERSECTIONS
		if (cntr >= START_PRINTING) {
			sprintf_s(buf, "%2d, ", *(pBlockIdx + i));
			outString(buf, outFile());
		}
		if (!*pntr) {
			sprintf_s(buf, "\n\n!!! PROBLEM:  nRow = %d  nColAbs = %d  nJ = %d\n\n", nRow, nColAbs, *(pBlockIdx + i));
			outString(buf, outFile());
			fclose(outFile());
		}
#endif
		assert(*pntr != 0);
		*pntr = 0;
	}
#if PRINT_INTERSECTIONS
	if (cntr >= START_PRINTING) {
		sprintf_s(buf, ")  cntr = %zd  ccc = %zd\n", cntr, ccc);
		outString(buf, outFile());
	}
#endif
}

FClass2(CMatrixCanonChecker, bool)::CheckBlockIntersections(
	T nRow, T b, const T* pRowSolution, T* pBlockIdx, T partIdx)
{
	// Function to check the block intersections for Kirkman Triple Systems
	const auto* pColOrbit = this->colOrbit(nRow, partIdx);
	const auto colOrbLen = this->colOrbitLen();
	const auto nColBase = partIdx * classSize();

	const auto* pColOrbitIni = this->colOrbitIni(nRow, partIdx);
	ColOrbPntr pColOrbitLast = NULL;
	while (pColOrbit) {
		if (*pRowSolution++) {
			// Define the number of columns to start with
			const auto nColCurr = ((char*)pColOrbit - (char*)pColOrbitIni) / colOrbLen;
			const auto nColAbs = static_cast<T>(nColBase + nColCurr);
			// Check the intersections of the current block with blocks
			// containing the current element
			auto* pCommonElemNumber = commonElemNumber() + nColAbs * b;
			for (auto i = partIdx; i--;) {
				auto* pntr = pCommonElemNumber + *(pBlockIdx + i);
				if (*pntr) {   // corresponding blocks already have one common elements.
					// Rollback of newly marked intersections.
					for (auto j = partIdx; --j > i;)
						*(pCommonElemNumber + *(pBlockIdx + j)) = 0;

					return false;
				}

				*pntr = 1; // the first intersection of corresponding blocks
			}

			*(pBlockIdx + partIdx) = nColAbs;
			break;
		}

		pColOrbit = pColOrbit->next();
	}

	return true;
}

FClass2(CMatrixCanonChecker, void)::GenerateBinaryColumnOrbits(T nRow, S *pRow, T *pColPermut) const {
	// Create column orbits for binary row number `nRow` 
	// of the matrix() and make it canonical, if possible
	if (!pRow)
		pRow = matrix()->GetRow(nRow);

	const auto* pColOrbit = colOrbit(nRow++);
	CColOrbit<T>* pColOrbitLast;
	const auto v = matrix()->rowNumb();
	const auto b = !pColPermut ? matrix()->colNumb() : 0;
	auto* pNewColOrbit = (nRow < v) ? colOrbitIni(nRow) : NULL;
	T j, jMax, lenOrb, type;
	j = jMax = 0;
	while (pColOrbit) {
		T n = 0;
		type = pRow[j];
		jMax += (lenOrb = pColOrbit->length());
		if (lenOrb > 1) {
			n = type;
			while (++j < jMax)
				n += pRow[j];

			if (type && n < lenOrb || n) {
				// If `type` is 0 and `n` is not, the orbit started with 0, but at least one 1 was found within it.
				type = 1;   // we need to change type
				// Orbit is divided into two parts
				T j2, j1 = (j2 = j) - lenOrb;
				const T jLast = j1 + n;
				while (j1 < jLast) {
					// Find first 0 in current fragment
					while (pRow[j1]) j1++;
					if (j1 == jLast)
						break;

					// Find last 1 in current fragment
					while (!pRow[--j2]);

					if (pNewColOrbit) {
						// Rearranging corresponding elements of the column permutation
						if (!pColPermut) {
							// Permutation of two columns in the remaining rows of the matrix
							auto *pNextRow = pRow;
							for (auto i1 = nRow; i1 < v; i1++) {
								pNextRow += b;
								if (pNextRow[j1] != pNextRow[j2])
									pNextRow[j2] = 1 - (pNextRow[j1] = pNextRow[j2]);
							}
						} else {
							auto tmp = pColPermut[j1];
							pColPermut[j1] = pColPermut[j2];
							pColPermut[j2] = tmp;
						}
					}

					pRow[j1++] = 1;   // When we going back, first we change index.
					pRow[j2] = 0;     // This is why we have no symmetry here.
				}

				if (pNewColOrbit) {
					pNewColOrbit = pNewColOrbit->InitOrbit(n, colOrbitLen(), pColOrbit, 1);
					type = 0;
				}
			}
			else
				n = 0;
		}
		else
			j = jMax;

		if (pNewColOrbit)
			pNewColOrbit = (pColOrbitLast = pNewColOrbit)->InitOrbit(lenOrb - n, colOrbitLen(), pColOrbit, type);

		pColOrbit = pColOrbit->next();
	}
	
	if (pNewColOrbit)
		pColOrbitLast->setNext(NULL);
}

FClass2(CMatrixCanonChecker, void)::GenerateColumnOrbits(T nRow, S* pRow, T* pColPermut) const {
	// Create column orbits for row number `nRow`  
	// of the matrix() and make it canonical, if possible
	// NOTE: This function is not yet written as it should be.
	if (!pRow)
		pRow = matrix()->GetRow(nRow);

	const auto* pColOrbit = colOrbit(nRow++);
	CColOrbit<T>* pColOrbitLast;
	const auto v = matrix()->rowNumb();
	const auto b = !pColPermut ? matrix()->colNumb() : 0;
	auto* pNewColOrbit = (nRow < v) ? colOrbitIni(nRow) : NULL;
	T j, jMax, lenOrb, type;
	j = jMax = 0;
	while (pColOrbit) {
		T n = 0;
		type = pRow[j];
		jMax += (lenOrb = pColOrbit->length());
		if (lenOrb > 1) {
			n = 1;
			while (++j < jMax)
				if (pRow[j] == type)
					n++;

			if (type && n < lenOrb) {
				// Orbit is divided into multiple parts

				T j2, j1 = (j2 = j) - lenOrb;
				const T jLast = j1 + n;
				while (j1 < jLast) {
					// Find first element which is not equal to type
					while (pRow[j1] == type) j1++;
					if (j1 == jLast)
						break;

					// Find last 1 in current fragment
					while (!pRow[--j2]);

					if (pNewColOrbit) {
						// Rearranging corresponding elements of the column permutation
						if (!pColPermut) {
							// Permutation of two columns in the remaining rows of the matrix
							auto* pNextRow = pRow;
							for (auto i1 = nRow; i1 < v; i1++) {
								pNextRow += b;
								if (pNextRow[j1] != pNextRow[j2])
									pNextRow[j2] = 1 - (pNextRow[j1] = pNextRow[j2]);
							}
						}
						else {
							auto tmp = pColPermut[j1];
							pColPermut[j1] = pColPermut[j2];
							pColPermut[j2] = tmp;
						}
					}

					pRow[j1++] = 1;   // When we going back, first we change index.
					pRow[j2] = 0;     // This is why we have no symmetry here.
				}

				if (pNewColOrbit) {
					// NOTE: We need to:
					// a) add `type` as the input parameter of the following method;
					// b) use the count of each type's occurrences in the previous rows, instead of the column's weight.
					pNewColOrbit = pNewColOrbit->InitOrbit(n, colOrbitLen(), pColOrbit, 1);
					type = 0;
				}
			}
			else
				n = 0;
		}
		else
			j = jMax;

		if (pNewColOrbit)
			pNewColOrbit = (pColOrbitLast = pNewColOrbit)->InitOrbit(lenOrb - n, colOrbitLen(), pColOrbit, type);

		pColOrbit = pColOrbit->next();
	}

	if (pNewColOrbit)
		pColOrbitLast->setNext(NULL);
}

#if 1
#define checkMatr(x, y, z)
#else
void checkMatr(unsigned char* pMatr, int v, int b)
{
	for (int i = 0; i < b; i++) {
		auto* pTo = pMatr + i;
		int l = 0;
		for (int j = 1; j < v; j++, pTo += b) {
			if (*pTo) {
				if (l++ == 3)
					break;
			}
		}
	}

	auto* pTo = pMatr;
	for (int i = 1; i < v; i++, pTo += b) {
		auto* pFrom = pTo;
		for (int j = i + 1; j < v; j++) {
			pFrom += b;
			int l = 0;
			for (auto k = b; k--;) {
				if (pFrom[k] && pTo[k]) {
					if (l)
						break;

					l++;
				}
			}
		}
	}
}
#endif

size_t lenToCompare = 0;
int compareRows(const void* arg1, const void* arg2) {
	return memcmp(arg2, arg1, lenToCompare);
}

FClass2(CMatrixCanonChecker, void)::sortRowsUpdateColumnOrbits(T v, T b, T nRowStart, bool initFlag)
{
	checkMatr(matrix()->GetRow(1), v, b);
	qsort(matrix()->GetRow(nRowStart), v - nRowStart, lenToCompare = b, compareRows);
	checkMatr(matrix()->GetRow(1), v, b);
	if (initFlag) {
		initiateColOrbits(v, 0, NULL, true);
		const auto iMax = CCanonicityChecker::rank();
		const auto lenOrb = b / iMax;
		const auto* pColOrbit = colOrbit(0);
		auto* pNewColOrbit = colOrbitIni(1);
		CColOrbit<T>* pPrev;
		for (T i = 0; i < iMax; i++) {
			pPrev = pNewColOrbit;
			pNewColOrbit = pPrev->InitOrbit(lenOrb, colOrbitLen(), pColOrbit, 0);
		}
		pPrev->setNext(NULL);
	}

	for (auto i = nRowStart; i < v; i++) {
		GenerateBinaryColumnOrbits(i, matrix()->GetRow(i));
		checkMatr(matrix()->GetRow(1), i, b);
	}
	checkMatr(matrix()->GetRow(1), v, b);
}

FClass2(CMatrixCanonChecker, S *)::CanonizeMatrix(int k, CanonicityCheckerPntr *ppClassGroup, T numClasses) {
	const auto b = matrix()->colNumb();
	const auto v = matrix()->rowNumb();
	T colBuffer[512], *pColumnBuf = b <= countof(colBuffer)? colBuffer : new T[b];
	const auto len = b * sizeof(pColumnBuf[0]);
	TestCanonParams<T, S> canonParam = { this, matrix(), 1 };

	S* pCompMatrix = NULL;
	S* pMatr = NULL;
	T* permClasses = NULL;
	T* pOrbits = NULL;
	T lenPerm, lenStab = 0;
	T nRowStart = 0;
	T numGroups = 0;
	const auto lenMatrix = (v - 1) * b;
	const auto lenData = lenMatrix * sizeof(S);
	CanonicityCheckerPntr pClassGroupHandle;

	if (k) {
		// As of 09/09/2024 we are here only when canonizing the K-SYSTEM's
		lenToCompare = b;
		sortRowsUpdateColumnOrbits(v, b, nRowStart = 1, true);
		// Preparing data for processing automorphisms associated with day permutations
		numGroups = (v - 1) / k;
		//numClasses = lenStab = lenPerm = (v - 2) / (k - 1);
		lenStab = lenPerm = numClasses;// = (v - 2) / (k - 1);
		pMatr = matrix()->GetRow(1);
	}

	int cmp = 1;
	T idx = IDX_MAX;
	while (true) {
		while (!TestCanonicity(v, &canonParam, t_saveRowPermutations)) {
			// Re-arrange the matrix rows according to the permutation found by permRow() 
			auto* pPermRow = permRow();
			bool adjustColumnOrbits = false;
			int firstAdjustedRow = -1;
			for (T tmp, from, i = nRowStart; i < v; i++) {
				T* pRow;
				if (i != (from = pPermRow[i])) {
					T* pFrom;
					memcpy(pColumnBuf, pFrom = pRow = matrix()->GetRow(tmp = i), len);
					do {
						pPermRow[tmp] = tmp;
						T* pTo = pFrom;
						memcpy(pTo, pFrom = matrix()->GetRow(tmp = from), len);
					} while ((from = pPermRow[from]) != i);

					pPermRow[tmp] = tmp;
					memcpy(pFrom, pColumnBuf, len);
					adjustColumnOrbits = true;
				}
				else {
					// Column orbits need to be adjusted only at least one previous row was moved
					if (!adjustColumnOrbits)
						continue;

					pRow = matrix()->GetRow(i);
				}

				if (firstAdjustedRow < 0)
					firstAdjustedRow = i;

				GenerateBinaryColumnOrbits(i, pRow);
			}

			if (k && adjustColumnOrbits) {
				checkMatr(pMatr, v, b);
				sortRowsUpdateColumnOrbits(v, b, firstAdjustedRow);
				checkMatr(pMatr, v, b);
			}
		}

		if (!k)
			break;

		if (pCompMatrix) {
			cmp = memcmp(pMatr, pCompMatrix, lenData);
			if (!cmp) {
				// Non-trivial automorphism on classes found
				pClassGroupHandle->addAutomorphism(permClasses);
				idx = IDX_MAX;
			}
			else {
				if (cmp < 0)
					memcpy(pMatr, pCompMatrix, lenData);
			}
		}
		else {
			// Matrix for comparison was not yet created
			pClassGroupHandle = new CCanonicityChecker<TDATA_TYPES>(lenPerm, 0, 0, t_outRowOrbits | t_outRowPermute);
			permClasses = new T[numClasses];
			pCompMatrix = new S[lenMatrix];
		}
		if (cmp == 1) {
			pClassGroupHandle->initOrbits(permClasses);
			memcpy(pCompMatrix, pMatr, lenData);
			checkMatr(pCompMatrix, v, b);
		}

		idx = pClassGroupHandle->next_permutation(permClasses, idx);
		if (idx == ELEMENT_MAX)
			break;		// there are no more permutations of the classes

		idx = lenPerm;

		T i = 0;
		while (permClasses[i] == i) i++;

		// Permute the column groups based on the class permutations just constructed.
		while (i < lenPerm) {
			auto* pTo = pMatr + numGroups * i;
			const auto* pFrom = pCompMatrix + numGroups * permClasses[i++];
			for (auto i = v; --i;) {
				memcpy(pTo, pFrom, numGroups * sizeof(pTo[0]));
				pTo += b;
				pFrom += b;
			}
		}

		sortRowsUpdateColumnOrbits(v, b, nRowStart = 1);
		checkMatr(pMatr, v, b);
	}

	pClassGroupHandle->updateOrderOfGroup();
	extraGroupOrder()->setGroupOrder(pClassGroupHandle->groupOrder());

	delete[] permClasses;
	delete[] pCompMatrix;
	if (ppClassGroup)
		*ppClassGroup = pClassGroupHandle;
	else
		delete pClassGroupHandle;

	if (pColumnBuf != colBuffer)
		delete[] pColumnBuf;

	return pMatr;
}
