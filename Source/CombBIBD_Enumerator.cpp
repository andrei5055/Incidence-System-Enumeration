#include "CombBIBD_Enumerator.h"

template class CCombBIBD_Enumerator<TDATA_TYPES>;

#if !CONSTR_ON_GPU
FClass2(CCombBIBD_Enumerator, int)::addLambdaInfo(char *buf, size_t lenBuffer, int *pLambdaSetSize) const {
	const auto lambdaSet = (static_cast<Class2(CCombinedBIBD) *>(this->getInSys()))->paramSet(t_lSet);
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