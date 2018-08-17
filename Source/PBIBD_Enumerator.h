#pragma once
#include "BIBD_Enumerator.h"

template<class T>
class CPBIBD_Enumerator : public CBIBD_Enumerator<T>
{
public:
	CK CPBIBD_Enumerator(const C_InSys<T> *pBIBD, bool matrOwner = false, bool noReplicatedBlocks = false, int treadIdx = -1, uint nCanonChecker = 0) :
		CBIBD_Enumerator<T>(pBIBD, matrOwner, noReplicatedBlocks, treadIdx, nCanonChecker) {}
	CK virtual bool isPBIB_enumerator() const			{ return true; }
protected:
	CK virtual bool checkLambda(T lambdaCur) const;
	CK virtual void ReportLamdaProblem(T i, T j, T lambda) const;
  CK const char *getObjName() const             { return "PBIBD"; }
	CK virtual size_t addLambdaInfo(char *buffer, size_t lenBuffer) const;
};

template<class T>
bool CPBIBD_Enumerator<T>::checkLambda(T lambdaCur) const {
	const auto lambdaSet = this->getInSys()->GetNumSet(t_lSet);
	for (size_t i = 0; i < lambdaSet->GetSize(); i++) {
		if (lambdaSet->GetAt(i) == lambdaCur)
			return true;
	}

	return false;
}

template<class T>
size_t CPBIBD_Enumerator<T>::addLambdaInfo(char *buf, size_t lenBuffer) const {
	const auto lambdaSet = this->getInSys()->GetNumSet(t_lSet);
	auto pBuf = buf;
	for (size_t i = 0; i < lambdaSet->GetSize(); i++) {
		if (pBuf == buf)
			pBuf += SNPRINTF(buf, lenBuffer, "{%d", lambdaSet->GetAt(i));
		else
			pBuf += SNPRINTF(pBuf, lenBuffer - (pBuf - buf), ", %d", lambdaSet->GetAt(i));
	}
	pBuf += SNPRINTF(pBuf, lenBuffer - (pBuf - buf), "}");
	return pBuf - buf;
}

template<class T>
void CPBIBD_Enumerator<T>::ReportLamdaProblem(T i, T j, T lambda) const {
	char buf[128];
	addLambdaInfo(buf, sizeof(buf));
	OUT_STRING(buff, 256, "Wrong number of common units in the rows (" ME_FRMT ", " ME_FRMT "): " ME_FRMT " is not in %s\n",
		i, j, lambda, buf);
}

template<class T>
class CInconsistentGraph_Enumerator : public CPBIBD_Enumerator<T>
{
public:
	CK CInconsistentGraph_Enumerator(const C_InSys<T> *pBIBD, bool matrOwner = false, bool noReplicatedBlocks = false, int treadIdx = -1, uint nCanonChecker = 0) :
		CPBIBD_Enumerator<T>(pBIBD, matrOwner, noReplicatedBlocks, treadIdx, nCanonChecker) {}
protected:
	CK virtual bool TestFeatures(CEnumInfo<T> *pEnumInfo, const CMatrixData<T> *pMatrix, int *pMatrFlags = NULL, CEnumerator<T> *pEnum = NULL) const;
	CK const char *getObjName() const             { return "I-Graph"; }
};

template<class T>
bool CInconsistentGraph_Enumerator<T>::TestFeatures(CEnumInfo<T> *pEnumInfo, const CMatrixData<T> *pMatrix, int *pMatrFlags, CEnumerator<T> *pEnum) const
{
	if (!CPBIBD_Enumerator<T>::TestFeatures(pEnumInfo, pMatrix, pMatrFlags))
		return false;

	if (!this->groupIsTransitive())
		return false;

	*pMatrFlags |= t_trahsitiveGroup;

	// Need to check that transposed matrix is not isomorphic to constructed one
	CMatrixData<T> transpMatr;
	transpMatr.InitTransposed(pMatrix);
	transpMatr.GetDataPntr()[0] = 0;
	CMatrixCol<T> matrCol(&transpMatr);
	T level;
	bool flg = pEnum->TestCanonicity(transpMatr.rowNumb(), &matrCol, t_saveRowToChange + t_saveRowPermutations, &level);

	return true;
}