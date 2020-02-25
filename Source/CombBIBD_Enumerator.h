#pragma once
#include "BIBD_Enumerator.h"

Class2Def(CCombBIBD_Enumerator) : public Class2(CBIBD_Enumerator)
{
public:
	CCombBIBD_Enumerator(const InSysPntr pBIBD, uint enumFlags = t_enumDefault, int treadIdx = -1, uint nCanonChecker = 0) :
		Class2(CBIBD_Enumerator)(pBIBD, enumFlags, treadIdx, nCanonChecker) {}
protected:
	CK virtual const char* getObjName() const		{ return "CBIBD"; }
	virtual const char* getTopLevelDirName() const	{ return "Combined_BIBDs"; }
#if !CONSTR_ON_GPU
	CK virtual int addLambdaInfo(char *buffer, size_t lenBuffer, const char* pFrmt = NULL, int* pLambdaSetSize = NULL) const;
	CK virtual int getJobTitleInfo(char *buffer, int lenBuffer) const;
	virtual void getEnumerationObjectKey(char* pInfo, int len) const;
	virtual const char* getObjNameFormat() const	{ return "  %14s:      "; }
#endif
};

