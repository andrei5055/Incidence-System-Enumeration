#pragma once
#include "PBIBD_Enumerator.h"

template<class T>
class CIG_Enumerator : public CPBIBD_Enumerator<T>
{
public:
	CK CIG_Enumerator(const C_InSys<T> *pBIBD, unsigned int enumFlags = t_enumDefault, bool firstPath = false, int treadIdx = -1, uint nCanonChecker = 0) :
		CPBIBD_Enumerator<T>(pBIBD, enumFlags, treadIdx, nCanonChecker), m_firstPath(firstPath) {
		const auto pLambdaSet = this->getInSys()->GetNumSet(t_lSet);
		if (pLambdaSet->GetSize() > 2) {
			const auto v = rowNumb();
			m_pRowIntersections = new T[v * (v - 1) / 2];
		}
		else
			m_pRowIntersections = NULL;
	}

	CK ~CIG_Enumerator()								{ delete[] rowIntersections(); }
protected:
	CK bool TestFeatures(CEnumInfo<T> *pEnumInfo, const CMatrixData<T> *pMatrix, int *pMatrFlags = NULL, CEnumerator<T> *pEnum = NULL) const override;
	bool makeFileName(char *buffer, size_t lenBuffer, const char *ext) const;
	CK const char *getObjName() const					{ return "I-Graph"; }
	CK bool fileExists(const char *path, bool file = true) const override;
	CK bool createNewFile(const char *fName) const override;
	CK bool SeekLogFile() const override				{ return true; }
	CK T *rowIntersections() const						{ return m_pRowIntersections; }
	CK bool prepareToFindRowSolution() override;
private:
	inline bool firstPath() const						{ return m_firstPath; }
	bool CheckConstructedBlocks(T nRow, T k);

	const bool m_firstPath;
	T *m_pRowIntersections;
};
