#pragma once
#include "BIBD_Enumerator.h"

template<class T>
class CPBIBD_Enumerator : public CBIBD_Enumerator<T>
{
public:
	CK CPBIBD_Enumerator(const C_InSys<T> *pBIBD, unsigned int enumFlags = t_enumDefault, int treadIdx = -1, uint nCanonChecker = 0) :
		CBIBD_Enumerator<T>(pBIBD, enumFlags, treadIdx, nCanonChecker) {}
	CK virtual bool isPBIB_enumerator() const					{ return true; }
protected:
	CK virtual bool checkLambda(size_t lambdaCur) const			{ return findLambda(lambdaCur) >= 0; }
	CK int findLambda(size_t lambdaCur) const;
	CK virtual void ReportLamdaProblem(T i, T j, size_t lambda) const;
	CK const char *getObjName() const override					{ return "PBIBD"; }
	CK virtual int addLambdaInfo(char *buffer, size_t lenBuffer) const;
};

template<class T>
int CPBIBD_Enumerator<T>::findLambda(size_t lambdaCur) const {
	const auto lambdaSet = this->getInSys()->GetNumSet(t_lSet);
	for (int i = 0; i < lambdaSet->GetSize(); i++) {
		if (lambdaSet->GetAt(i) == lambdaCur)
			return i;
	}

	return -1;
}

template<class T>
int CPBIBD_Enumerator<T>::addLambdaInfo(char *buf, size_t lenBuffer) const {
	const auto lambdaSet = this->getInSys()->GetNumSet(t_lSet);
	auto pBuf = buf;
	for (size_t i = 0; i < lambdaSet->GetSize(); i++) {
		if (pBuf == buf)
			pBuf += SNPRINTF(buf, lenBuffer, "{%d", lambdaSet->GetAt(i));
		else
			pBuf += SNPRINTF(pBuf, lenBuffer - (pBuf - buf), ", %d", lambdaSet->GetAt(i));
	}
	pBuf += SNPRINTF(pBuf, lenBuffer - (pBuf - buf), "}");
	return (int)(pBuf - buf);
}

template<class T>
void CPBIBD_Enumerator<T>::ReportLamdaProblem(T i, T j, size_t lambda) const {
	char buf[128];
	addLambdaInfo(buf, sizeof(buf));
	OUT_STRING(buff, 256, "Wrong number of common units in the rows (" ME_FRMT ", " ME_FRMT "): %zu is not in %s\n",
		i, j, lambda, buf);
}
