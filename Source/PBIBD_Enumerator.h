#pragma once
#include "BIBD_Enumerator.h"

template<class T>
class CPBIBD_Enumerator : public CBIBD_Enumerator<T>
{
public:
	CK CPBIBD_Enumerator(const C_InSys<T> *pBIBD, bool matrOwner = false, bool noReplicatedBlocks = false, int treadIdx = -1, uint nCanonChecker = 0) :
		CBIBD_Enumerator(pBIBD, matrOwner, noReplicatedBlocks, treadIdx, nCanonChecker) {}
	CK virtual bool isPBIB_enumerator() const			{ return true; }
protected:
	CK virtual bool checkLambda(T lambdaCur) const;
	CK virtual void ReportLamdaProblem(T i, T j, T lambda) const;

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
void CPBIBD_Enumerator<T>::ReportLamdaProblem(T i, T j, T lambda) const {
	char buf[128], *pBuf = buf;
	const auto lambdaSet = this->getInSys()->GetNumSet(t_lSet);
	bool first = true;
	for (size_t i = 0; i < lambdaSet->GetSize(); i++) {
		if (first)
			pBuf += SNPRINTF(buf, pBuf - buf, "{%d", lambdaSet->GetAt(i));
		else
			pBuf += SNPRINTF(buf, pBuf - buf, ", %d", lambdaSet->GetAt(i));
	}

	SNPRINTF(pBuf, pBuf - buf, "}");
	OUT_STRING(buff, 256, "Wrong number of common units in the rows (" ME_FRMT ", " ME_FRMT "): " ME_FRMT " != %s\n",
		i, j, lambda, buf);
}

template<class T>
class C_UncoordinatedGraph_Enumerator : public CPBIBD_Enumerator<T>
{
public:
	CK C_UncoordinatedGraph_Enumerator(const C_InSys<T> *pBIBD, bool matrOwner = false, bool noReplicatedBlocks = false, int treadIdx = -1, uint nCanonChecker = 0) :
		CPBIBD_Enumerator(pBIBD, matrOwner, noReplicatedBlocks, treadIdx, nCanonChecker) {}
protected:
	CK virtual bool TestFeatures(CEnumInfo<T> *pEnumInfo, const CMatrixData<T> *pMatrix, int *pMatrFlags) const;
};

template<class T>
bool C_UncoordinatedGraph_Enumerator<T>::TestFeatures(CEnumInfo<T> *pEnumInfo, const CMatrixData<T> *pMatrix, int *pMatrFlags) const
{
	if (!CPBIBD_Enumerator::TestFeatures(pEnumInfo, pMatrix, pMatrFlags))
		return false;

	return true;
}