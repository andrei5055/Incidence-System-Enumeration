#include "CombBIBD_Enumerator.h"

template class CCombBIBD_Enumerator<TDATA_TYPES>;

#if !CONSTR_ON_GPU
FClass2(CCombBIBD_Enumerator, int)::addLambdaInfo(char *buf, size_t lenBuffer, const char* pFormat, size_t *pLambdaSetSize) const {
	return addLambdaInform(paramSet(t_lSet), buf, lenBuffer, pLambdaSetSize);
}

FClass2(CCombBIBD_Enumerator, int)::getJobTitleInfo(char* buffer, int lenBuffer) const {
	const auto inSys = this->getInSys();
	return SNPRINTF(buffer, lenBuffer, "%s(%3" _FRMT", %2" _FRMT", ", getObjName(), inSys->rowNumbExt(), inSys->GetK());
}
#endif

FClass2(CCombBIBD_Enumerator, void)::getEnumerationObjectKey(char* pKey, int len) const {
	char buffer[64];
	this->makeJobTitle(this->designParams(), buffer, countof(buffer));
	SNPRINTF(pKey, len, "%s", buffer + strlen(this->getObjName()));
}

FClass2(CCombBIBD_Enumerator, char*)::getEnumerationObjectKeyA(char* pKey, int len) const {
	getEnumerationObjectKey(pKey, len);
	auto* pntr = strstr(pKey, "})");
	*pntr = '\0';
	return pKey;
}

FClass2(CCombBIBD_Enumerator, RowSolutionPntr)::setFirstRowSolutions() {
	auto pSolutions = this->rowStuff(1);
	const auto* pR_set = this->paramSet(t_rSet);
	for (auto i = this->numParts(); i--;) {
		auto pPartSolutions = pSolutions + i;
		pPartSolutions->InitSolutions(1, 1);
		pPartSolutions->AddElement(pR_set->GetAt(i));
	}

	return pSolutions;
}

FClass2(CCombBIBD_Enumerator, void)::CreateForcedRows() {
	// Combined BIBDs having n parts will be constructed with first "artificial" row:
	//  n, n, ...,n, n-1, ..., n-1, n-2,..., n-2,... 2,...2,1...,1
	// It's why we need to start enumeration with the row number 1.
	this->setCurrentRowNumb(1);
	CreateFirstRow();
}

FClass2(CCombBIBD_Enumerator, void)::CreateFirstRow(S* pFirstRow)
{
	const auto pInSys = this->getInSys();
	const auto k = pInSys->GetNumSet(t_kSet)->GetAt(0);
	const auto v = pInSys->rowNumb() - 1;
	const auto pR_set = this->paramSet(t_rSet);
	const auto iMax = static_cast<T>(this->numParts());
	if (!pFirstRow)
		pFirstRow = pInSys->GetRow(0);

	for (T i = 0; i < iMax; i++) {
		const auto b = pR_set->GetAt(i) * v / k;
		rowSetFragm(pFirstRow, iMax - i, b);
		pFirstRow += b;
	}
}

FClass2(CCombBIBD_Enumerator, CMatrixData<T, S> *)::CreateSpareMatrix(const EnumeratorPntr pMaster)
{
	const auto pMatr = pMaster ? pMaster->matrix() : this->matrix();
	MatrixDataPntr pSpareMatrix = new CMatrixData<T, S>();
	pSpareMatrix->Init(pMatr->rowNumb(), pMatr->colNumb());
	CreateFirstRow(pSpareMatrix->GetRow(0));
	auto *pPartsInformation = pSpareMatrix->InitPartsInfo(this->numParts());
	pPartsInformation->CopyPartInfo(pMatr->partsInfo());
	return pSpareMatrix;
}

FClass2(CCombBIBD_Enumerator, CGroupOnParts<T> *)::makeGroupOnParts(const CanonicityCheckerPntr owner) {
	auto lanbdaSet = paramSet(t_lSet);
	const int jMax = lanbdaSet->GetSize() - 1;
	CVector<T> lengths;
	T prevLambda = 0;
	uint count;
	uint factorial;
	for (int j = 0; j <= jMax; j++) {
		const auto lambda = lanbdaSet->GetAt(j);
		if (prevLambda == lambda) {
			factorial *= (++count);
			if (j < jMax)
				continue;
			j++;
		}

		if (prevLambda) {
			if (factorial > 1) {
				lengths.AddElement(j - count); 	// index of the first BIBD with the same lambda
				lengths.AddElement(count);		// number of BIBDs with the same lambda
				lengths.AddElement(factorial);  // group order
			}
		}
		prevLambda = lambda;
		factorial = count = 1;
	}
	if (!lengths.GetSize())
		return NULL;


	updateCanonicityChecker(rowNumb(), colNumb());
	auto pGroupOnParts = new CGroupOnParts<T>(owner, lengths, 3);
	InitGroupOderStorage(pGroupOnParts);
	return pGroupOnParts;
}
