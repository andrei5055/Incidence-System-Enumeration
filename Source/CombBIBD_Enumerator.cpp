#include "CombBIBD_Enumerator.h"

template class CCombBIBD_Enumerator<TDATA_TYPES>;

#if !CONSTR_ON_GPU
FClass2(CCombBIBD_Enumerator, int)::addLambdaInfo(char *buf, size_t lenBuffer, const char* pFormat, int *pLambdaSetSize) const {
	return addLambdaInform(paramSet(t_lSet), buf, lenBuffer, pLambdaSetSize);
}

FClass2(CCombBIBD_Enumerator, int)::getJobTitleInfo(char* buffer, int lenBuffer) const {
	const auto inSys = this->getInSys();
	return SNPRINTF(buffer, lenBuffer, "%s(%3" _FRMT", %2" _FRMT", ", getObjName(), inSys->rowNumbExt(), inSys->GetK());
}
#endif

FClass2(CCombBIBD_Enumerator, void)::getEnumerationObjectKey(char* pInfo, int len) const {
	char buffer[64];
	this->makeJobTitle(this->designParams(), buffer, countof(buffer));
	SNPRINTF(pInfo, len, "%s", buffer + strlen(this->getObjName()));
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
	// Combined BIBDs having n parts will be constructed with first "atrtificial" row:
	//  n, n, ...,n, n-1, ..., n-1, n-2,..., n-2,... 2,...2,1...,1
	// It's why we need to start enumeration with the row number 1.
	this->setCurrentRowNumb(1);
	const auto pInSys = this->getInSys();
	const auto k = pInSys->GetNumSet(t_kSet)->GetAt(0);
	const auto v = pInSys->rowNumb() - 1;
	const auto pR_set = this->paramSet(t_rSet);
	const auto iMax = static_cast<T>(this->numParts());
	auto* pRow = pInSys->GetRow(0);
	for (T i = 0; i < iMax; i++) {
		const auto b = pR_set->GetAt(i) * v / k;
		rowSetFragm(pRow, iMax - i, b);
		pRow += b;
	}
}
