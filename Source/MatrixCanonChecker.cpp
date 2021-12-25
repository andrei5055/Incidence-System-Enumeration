#include "MatrixCanonChecker.h"

template class CMatrixCanonChecker<TDATA_TYPES>;

FClass2(CMatrixCanonChecker, ColOrbPntr)::MakeRow(T nRow, const T *pRowSolution, bool nextColOrbNeeded, T partIdx) const
{
	auto* pRow = this->matrix()->ResetRowPart(nRow, partIdx);
	if (nextColOrbNeeded)
		nextColOrbNeeded &= nRow + 1 < matrix()->rowNumb();

	const auto* pColOrbit = this->colOrbit(nRow, partIdx);
	const auto* pNextRowColOrbit = this->colOrbitIni(nRow + 1, partIdx);
	const auto colOrbLen = this->colOrbitLen();

	const auto rank = CCanonicityChecker<T, S>::rank();
	const auto maxElement = rank - 1;
	const auto* pColOrbitIni = this->colOrbitIni(nRow, partIdx);

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

