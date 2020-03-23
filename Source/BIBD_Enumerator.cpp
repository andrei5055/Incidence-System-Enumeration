#include "BIBD_Enumerator.h"

template class CBIBD_Enumerator<TDATA_TYPES>;

FClass2(CBIBD_Enumerator, int)::unforcedElement(const CColOrbit<S> *pOrb, int nRow) const
{
	const size_t diffWeight = this->getInSys()->GetK() - pOrb->columnWeight();
	if (diffWeight)
		return diffWeight == this->rowNumb() - nRow ? 1 : -1;

	// All units are there
	return 0;
}

FClass2(CBIBD_Enumerator, bool)::isValidSolution(const VECTOR_ELEMENT_TYPE* pSol) const
{
	// Check if solution is valid (for elimination of invalid solutions)
	auto rowNumb = this->currentRowNumb();
	if (rowNumb <= firtstNonfixedRowNumber())
		return true;

	// For canonical BIBD the number of blocks containing any of any three elements cannot be
	// bigger than the number of blocks containing first, second and third element.
	// Let's check it
	const auto lambda = this->lambda();
	const auto x0_3 = this->getX0_3();
	if (lambda == x0_3)     // Intersection of first three rows is maximal
		return true;		// Nothing to test

	const auto k = this->getInSys()->GetK();
	auto limit = rowNumb + k + 1;
	if (limit >= this->rowNumb() && rowNumb + 3 < this->rowNumb()) {
		// Theorem: The number of columns of canonical matrix which are forcible constructed by units cannot be bigger than 
		// number of blocks containing first, second and third element.
		// We should start to check this condition on the first row which tested solutions could create
		// first column forcible constructed by units: rowNumb >= this->rowNumb() - k - 1 AND
		// at least three rows need to be constructed: rowNumb < this->rowNumb() - 3
		const auto* pColOrbit = this->colOrbit(rowNumb);
		const auto* pRowSolution = pSol;
		limit -= this->rowNumb();
		auto nForcible = forcibleLambda(rowNumb, 0);
		while (pColOrbit) {
			if (pColOrbit->columnWeight() == limit) {
				// Define the number of new columns that will be enforceable completed by units
				const auto newEnforsed = pColOrbit->length() - *pRowSolution;
				if (newEnforsed && (nForcible += newEnforsed) > x0_3)
					return false;
			}

			pRowSolution++;
			pColOrbit = pColOrbit->next();
		}
	}

	this->MakeRow(pSol);

	// Define intersection of current row with previous one:
	const auto* pMathix = this->matrix();
	const auto* pCurrRow = pMathix->GetRow(rowNumb--);
	const auto* pPrevRow = pMathix->GetRow(rowNumb--);
	const auto colNumb = this->colNumb();

	int columns[32];
	assert(lambda <= countof(columns));
	int idx = 0;
	T j = 0;
	while (true) {
		if (pCurrRow[j] && pPrevRow[j]) {
			columns[idx++] = j;
			if (idx == lambda)
				break;
		}

		++j;
	}

	// Let's check the necessary and sufficient condition
	//     for first matrix row:
	if (columns[x0_3] < this->getR())
		return false;

	//     for remaining rows:
	do {
		pCurrRow = pMathix->GetRow(rowNumb);
		auto j = x0_3;
		for (auto i = lambda; i--;) {
			if (pCurrRow[columns[i]]) {
				if (!j--)
					return false;
			}
		}
	} while (--rowNumb);

	return true;
}

FClass2(CBIBD_Enumerator, void)::getEnumerationObjectKey(char *pInfo, int len) const {
	SNPRINTF(pInfo, len, "(%3" _FRMT", %2" _FRMT", %2" _FRMT")",
		this->rowNumb(), this->getInSys()->GetK(), this->getInSys()->lambda());
}

#if !CONSTR_ON_GPU
FClass2(CBIBD_Enumerator, bool)::makeFileName(char* buffer, size_t lenBuffer, const char* ext) const
{
	const auto inSys = this->getInSys();
	auto len = this->getDirectory(buffer, lenBuffer);
	len += SNPRINTF(buffer + len, lenBuffer - len, ME_FRMT"_" ME_FRMT"_", inSys->rowNumbExt(), inSys->GetK());
	len += this->addLambdaInfo(buffer + len, lenBuffer - len, ME_FRMT);
	SNPRINTF(buffer + len, lenBuffer - len, "%s", ext ? ext : FILE_NAME(""));
	return true;
}

FClass2(CBIBD_Enumerator, bool)::makeJobTitle(const designParam *pParam, char *buffer, int lenBuffer, const char *comment) const
{
	int lambdaSetSize = 0;
	auto len = getJobTitleInfo(buffer, lenBuffer);
	len += addLambdaInfo(buffer + len, lenBuffer - len, "%2" _FRMT, &lambdaSetSize);

	if (pParam->lambdaSizeMax() > lambdaSetSize) {
		auto maxSize = pParam->lambdaSizeMax() - lambdaSetSize;
		auto pBuf = buffer + len;
		pBuf += SNPRINTF(pBuf, lenBuffer - (pBuf - buffer), ")");
		while (maxSize-- > 0)
			pBuf += SNPRINTF(pBuf, lenBuffer - (pBuf - buffer), "   ");
	}
	else
		SNPRINTF(buffer + len, lenBuffer - len, ")%s", comment);

	return true;
}

FClass2(CBIBD_Enumerator, int)::getJobTitleInfo(char *buffer, int lenBuffer) const
{
	const auto v = this->rowNumb();
	const auto b = this->matrix()->colNumb();
	const auto k = this->getInSys()->GetK();
	return SNPRINTF(buffer, lenBuffer, "%s(%3" _FRMT", %3" _FRMT", %2" _FRMT", %2" _FRMT", ", getObjName(), v, b, b * k / v, k);
}

FClass2(CBIBD_Enumerator, int)::addLambdaInform(const Class1(CVector) *lambdaSet, char* buf, size_t lenBuffer, int* pLambdaSetSize) const
{
	const auto lambdaNumb = lambdaSet->GetSize();
	if (pLambdaSetSize)
		*pLambdaSetSize = static_cast<int>(lambdaNumb);

	const auto* pFrmt = "{%2d";
	int len = 0;
	for (size_t i = 0; i < lambdaNumb; i++) {
		len += SNPRINTF(buf + len, lenBuffer - len, pFrmt, lambdaSet->GetAt(i));
		pFrmt = ",%2d";
	}

	return len + SNPRINTF(buf + len, lenBuffer - len, "}");
}
#endif
