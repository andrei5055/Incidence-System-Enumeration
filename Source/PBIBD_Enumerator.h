#pragma once
#include "BIBD_Enumerator.h"

Class2Def(CPBIBD_Enumerator) : public Class2(CBIBD_Enumerator)
{
public:
	CK CPBIBD_Enumerator(const InSysPntr pBIBD, unsigned int enumFlags = t_enumDefault, int treadIdx = -1, uint nCanonChecker = 0) :
		Class2(CBIBD_Enumerator)(pBIBD, enumFlags, treadIdx, nCanonChecker) {}
	CK virtual bool isPBIB_enumerator() const					{ return true; }
protected:
	CK virtual bool checkLambda(size_t lambdaCur) const			{ return findLambda(lambdaCur) != -1; }
	CK size_t findLambda(size_t lambdaCur) const;
	CK virtual void ReportLamdaProblem(S i, S j, S lambda) const;
	CK const char *getObjName() const override					{ return "PBIBD"; }
	CK virtual int addLambdaInfo(char *buffer, size_t lenBuffer, const char* pFrmt = NULL, size_t *pLambdaSetSize = NULL) const {
		return addLambdaInform(this->getInSys()->GetNumSet(t_lSet), buffer, lenBuffer, pLambdaSetSize);
	}
};

FClass2(CPBIBD_Enumerator, size_t)::findLambda(size_t lambdaCur) const {
	const auto lambdaSet = this->getInSys()->GetNumSet(t_lSet);
	for (size_t i = 0; i < lambdaSet->GetSize(); i++) {
		if (lambdaSet->GetAt(i) == lambdaCur)
			return i;
	}

	return -1;
}

FClass2(CPBIBD_Enumerator, void)::ReportLamdaProblem(S i, S j, S lambda) const {
	char buf[128];
	addLambdaInfo(buf, sizeof(buf));
	OUT_STRING(buff, 256, "Wrong number of common units in the rows (" ME_FRMT ", " ME_FRMT "): " ME_FRMT " is not in %s\n",
		i, j, lambda, buf);
}
