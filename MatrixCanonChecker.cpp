#include "MatrixCanonChecker.h"
#include "RowSolution.h"

//template class CEnumInfo<TDATA_TYPES>;
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
			/*
			static int css;
			if (cntr > START_PRINTING) {
				css++;
				pCommonElemNumber = commonElemNumber() + 29 * b;
				fprintf(outFile(), "\ncss = %d: BlInt[29, 20] = %d", css, *(pCommonElemNumber + 20));
				fprintf(outFile(), "\ncss = %d:  (%d, %d)", css, *(pBlockIdx + 4), *(pBlockIdx + 5));
				if (css == 561037)
					css += 0;
			} */
			break;
		}

		pColOrbit = pColOrbit->next();
	}

	return true;
}

FClass2(CMatrixCanonChecker, void)::CreateColumnOrbits(T nRow, S *pRow, T *pColPermut) const {
	// Create column orbits for row number nRow of the matrix()
	// and make it canonical, if it is possible
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

			if (n && n < lenOrb) {
				// Orbit is split in two parts

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

size_t lenToCompare = 0;
int compareRows(const void* arg1, const void* arg2) {
	return memcmp(arg2, arg1, lenToCompare);
}

FClass2(CMatrixCanonChecker, void)::CanonizeMatrix(bool columnOrbitsConstructed) {
	const auto b = matrix()->colNumb();
	const auto v = matrix()->rowNumb();
	T colBuffer[512], *pColumnBuf = b <= countof(colBuffer)? colBuffer : new T[b];
	const auto len = b * sizeof(pColumnBuf[0]);
	TestCanonParams<T, S> canonParam = { this, matrix(), 1 };

	if (!columnOrbitsConstructed) {
		initiateColOrbits(v, 0, NULL, true);

		qsort(matrix()->GetRow(0), v, lenToCompare = b, compareRows);
		for (T i = 0; i < v; i++) {
			CreateColumnOrbits(i, matrix()->GetRow(i));
		}
	}

	while (!TestCanonicity(v, &canonParam, t_saveRowPermutations)) {
		// Rearrange the rows of matrix according to the permutation found by permRow() 
		auto* pPermRow = permRow();
		bool adjustColumnOrbits = false;
		for (T tmp, from, i = 0; i < v; i++) {
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
			} else {
				// Column orbits need to be adjusted only at least one previous row was moved
				if (!adjustColumnOrbits)
					continue;

				pRow = matrix()->GetRow(i);
			}

			CreateColumnOrbits(i, pRow);
		}
	}

	if (pColumnBuf != colBuffer)
		delete[] pColumnBuf;
}
