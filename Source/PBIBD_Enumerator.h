#pragma once
#include "BIBD_Enumerator.h"

IClass2Def(PBIBD_Enumerator) : public IClass2(BIBD_Enumerator)
{
public:
	CK CPBIBD_Enumerator(const InSysPntr pBIBD, unsigned int enumFlags = t_enumDefault, int treadIdx = -1, uint nCanonChecker = 0) :
		IClass2(BIBD_Enumerator)(pBIBD, enumFlags, treadIdx, nCanonChecker) {}
	CK virtual bool isPBIB_enumerator() const					{ return true; }
protected:
	CK virtual bool checkLambda(size_t lambdaCur) const			{ return findLambda(lambdaCur) != -1; }
	CK size_t findLambda(size_t lambdaCur) const;
	CK virtual void ReportLamdaProblem(S i, S j, size_t lambda) const;
	CK const char *getObjName() const override					{ return "PBIBD"; }
	CK virtual int addLambdaInfo(char *buffer, size_t lenBuffer, int *pLambdaSetSize = NULL) const;
};

TClass2(PBIBD_Enumerator, size_t)::findLambda(size_t lambdaCur) const {
	const auto lambdaSet = this->getInSys()->GetNumSet(t_lSet);
	for (size_t i = 0; i < lambdaSet->GetSize(); i++) {
		if (lambdaSet->GetAt(i) == lambdaCur)
			return i;
	}

	return -1;
}

TClass2(PBIBD_Enumerator, int)::addLambdaInfo(char *buf, size_t lenBuffer, int *pLambdaSetSize) const {
	const auto lambdaSet = this->getInSys()->GetNumSet(t_lSet);
	if (pLambdaSetSize)
		*pLambdaSetSize = static_cast<int>(lambdaSet->GetSize());

	auto pBuf = buf;
	for (size_t i = 0; i < lambdaSet->GetSize(); i++) {
		if (pBuf == buf)
			pBuf += SNPRINTF(buf, lenBuffer, "{%2d", lambdaSet->GetAt(i));
		else
			pBuf += SNPRINTF(pBuf, lenBuffer - (pBuf - buf), ",%2d", lambdaSet->GetAt(i));
	}

	pBuf += SNPRINTF(pBuf, lenBuffer - (pBuf - buf), "}");
	return (int)(pBuf - buf);
}

TClass2(PBIBD_Enumerator, void)::ReportLamdaProblem(S i, S j, size_t lambda) const {
	char buf[128];
	addLambdaInfo(buf, sizeof(buf));
	OUT_STRING(buff, 256, "Wrong number of common units in the rows (" ME_FRMT ", " ME_FRMT "): %zu is not in %s\n",
		i, j, lambda, buf);
}
