#include "BIBD_Enumerator.h"

template class CBIBD_Enumerator<MATRIX_ELEMENT_TYPE>;

template<class T>
int CBIBD_Enumerator<T>::unforcedElement(const CColOrbit<T> *pOrb, int nRow) const
{
	const size_t diffWeight = this->getInSys()->GetK() - pOrb->columnWeight();
	if (diffWeight)
		return diffWeight == this->rowNumb() - nRow ? 1 : -1;

	// All units are there
	return 0;
}

template<class T>
bool CBIBD_Enumerator<T>::isValidSolution(const VECTOR_ELEMENT_TYPE* pSol) const
{
//	return true;
	// Check if solution is valid (for elimination of invalid solutions)
	auto rowNumb = this->currentRowNumb();
	if (rowNumb <= 2)
		return true;

	const auto lambda = this->lambda();
	const auto x0_3 = this->getX0_3();
	if (lambda == x0_3)
		return true;		// Nothing to test

	this->MakeRow(pSol);

	// Define intersection of current row with previous one
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

	if (columns[x0_3] < this->getR())
		return false;

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

#if !CONSTR_ON_GPU
template<class T>
bool CBIBD_Enumerator<T>::makeFileName(char *buffer, size_t lenBuffer, const char *ext) const
{
	const auto dirLength = this->getDirectory(buffer, lenBuffer) ;
	SNPRINTF(buffer + dirLength, lenBuffer - dirLength, ME_FRMT"_" ME_FRMT"_" ME_FRMT"%s", this->rowNumb(),
						this->getInSys()->GetK(), this->getInSys()->lambda(), ext ? ext : FILE_NAME(""));
	return true;
}

template<class T>
bool CBIBD_Enumerator<T>::makeJobTitle(const designParam *pParam, char *buffer, int lenBuffer, const char *comment) const
{
	const auto v = this->rowNumb();
	const auto b = this->matrix()->colNumb();
	const auto k = this->getInSys()->GetK();
	auto len = SNPRINTF(buffer, lenBuffer, "%s(%3" _FRMT", %3" _FRMT", %2" _FRMT", %2" _FRMT", ", getObjName(), v, b, b * k / v, k);
	int lambdaSetSize = 0;
	len += addLambdaInfo(buffer + len, lenBuffer - len, &lambdaSetSize);

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
#endif
