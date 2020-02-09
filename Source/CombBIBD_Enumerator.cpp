#include "CombBIBD_Enumerator.h"

template class CCombBIBD_Enumerator<MATRIX_ELEMENT_TYPE>;

#if !CONSTR_ON_GPU
template<class T>
int CCombBIBD_Enumerator<T>::addLambdaInfo(char *buf, size_t lenBuffer, int *pLambdaSetSize) const {
	const auto lambdaSet = (static_cast<CCombinedBIBD<T> *>(this->getInSys()))->paramSet(t_lSet);
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

template<class T>
int CCombBIBD_Enumerator<T>::getJobTitleInfo(char* buffer, int lenBuffer) const {
	const auto v = this->rowNumb() - 1;
	const auto k = this->getInSys()->GetK();
	return SNPRINTF(buffer, lenBuffer, "%s(%3" _FRMT", %2" _FRMT", ", getObjName(), v, k);
}
#endif

template<class T>
void CCombBIBD_Enumerator<T>::getEnumerationObjectKey(char* pInfo, int len) const {
	char buffer[64];
	makeJobTitle(this->designParams(), buffer, countof(buffer));
	SNPRINTF(pInfo, len, "%s", buffer + strlen(this->getObjName()));
}