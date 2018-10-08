#pragma once
#include "PBIBD_Enumerator.h"

template<class T>
class CIG_Enumerator : public CPBIBD_Enumerator<T>
{
public:
	CK CIG_Enumerator(const C_InSys<T> *pBIBD, const designParam *pParam, unsigned int enumFlags = t_enumDefault, bool firstPath = false, int treadIdx = -1, uint nCanonChecker = 0) :
		CPBIBD_Enumerator<T>(pBIBD, enumFlags, treadIdx, nCanonChecker), m_firstPath(firstPath) {
		const auto lambdaSet = this->getInSys()->GetNumSet(t_lSet);
		const auto nLambd = lambdaSet->GetSize();

		m_pRowIntersections = nLambd > 2 || lambdaSet->GetAt(0) || lambdaSet->GetAt(1) != 1?
							  new T[lenRowIntersection(rowNumb())] : NULL;

		// Allocate memory to store current lambdaA, lambdaB for blocks
		// This set of lambda's is not known, but it cannot have more than k + 1 elements
		const auto len = this->getInSys()->GetK() + 1;
		m_pLambda[0] = new T[2 * len + nLambd];
		m_pLambda[2] = (m_pLambda[1] = m_pLambda[0] + nLambd) + len;

		// Converting the values lamdaB into T format
		const auto &coeffB = pParam->lambdaB();
		for (auto i = nLambd; i--;)
			*(m_pLambda[0] + i) = coeffB[i];
	}

	CK ~CIG_Enumerator()								{ delete[] rowIntersections(); delete[] lambdaBSrc(); }
	CK void CloneMasterInfo(const CEnumerator<T> *p, size_t nRow) override {
		auto pMaster = static_cast<const CIG_Enumerator<T> *>(p);
		copyRowIntersection(pMaster->rowIntersections());

		// No need to clone remaining information from master if it is not yet constructed
		if (nRow < this->getInSys()->GetK())
			return;

		memcpy(lambdaA(), pMaster->lambdaA(), lenLambdas());
	}

	CK inline T *rowIntersections() const				{ return m_pRowIntersections; }
	inline T *lambdaA() const							{ return m_pLambda[1]; }
	inline T *lambdaB() const							{ return m_pLambda[2]; }
protected:
	CK bool TestFeatures(CEnumInfo<T> *pEnumInfo, const CMatrixData<T> *pMatrix, 
		       int *pMatrFlags = NULL, CEnumerator<T> *pEnum = NULL) const override;
	bool makeFileName(char *buffer, size_t lenBuffer, const char *ext) const override;
	CK const char *getObjName() const override			{ return "I-Graph"; }
	CK bool fileExists(const char *path, bool file = true) const override;
	CK bool createNewFile(const char *fName) const override;
	CK bool SeekLogFile() const override				{ return true; }

	CK bool prepareToFindRowSolution() override;
private:
	CK void copyRowIntersection(const T *pntr)			{
		if (rowIntersections())
			memcpy(rowIntersections(), pntr, lenRowIntersection(rowNumb()) * sizeof(*pntr)); 
	}
	inline T *lambdaBSrc() const						{ return m_pLambda[0]; }
	inline bool firstPath() const						{ return m_firstPath; }
	inline size_t lenRowIntersection(T v) const			{ return v * (v - 1) / 2; }
	inline size_t nLambdas() const						{ return lambdaB() - lambdaA(); }
	inline size_t lenLambdas() const                    { return 2 * nLambdas() * sizeof(*lambdaA()); }
	bool CheckConstructedBlocks(T nRow, T k) const;
	bool DefineInterstructForBlocks(size_t nColCurr, T k, const T *pElementNumb, T i, T *lambdaACurrCol) const;

	const bool m_firstPath;
	T *m_pRowIntersections;
	T *m_pLambda[3];
};
