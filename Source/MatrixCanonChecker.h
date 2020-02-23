#pragma once
#include "ColOrbitManager.h"
#include "matrix.h"

template<class T>
class CMatrixCol : public CColOrbitManager<T>
{
public:
	CC CMatrixCol(const CMatrixData<T> *pMatrix, bool IS_enum, bool matrOwner) :
		CColOrbitManager(pMatrix->maxElement() + 1, pMatrix->rowNumb(), pMatrix->colNumb()) {
		initiateMatrixCol(pMatrix, IS_enum, matrOwner);
	}
	CC CMatrixCol(CMatrixData<T> *pMatrix, T rowNumb, T colNumb, T maxElem, bool IS_enum) :
		CColOrbitManager(maxElem + 1, rowNumb, colNumb) {
		initiateMatrixCol(pMatrix, IS_enum, false);
	}
	CC ~CMatrixCol() {
		if (isMatrOwner())
			delete static_cast<const CMatrix<T> *>(matrix());
	}
	CC inline void initiateMatrixCol(const CMatrixData<T> *pMatrix, bool IS_Enum = true, bool matrOwner = false) {
		m_pMatrix = pMatrix;
		setMatrOwner(matrOwner);
		setIS_Enumerator(IS_Enum);
		setOutFile(NULL);
	}
	CC inline const CMatrixData<T> *matrix() const { return m_pMatrix; }
	inline bool IS_enumerator() const { return m_bIS_Emunerator; }
	inline void closeFile() { if (outFile()) { fclose(outFile()); setOutFile(NULL); } }
	inline CC void setOutFile(FILE *file) { m_pFile = file; }
	inline auto outFile() const { return m_pFile; }
	inline auto outFilePntr() { return &m_pFile; }
protected:
	CC inline void setIS_Enumerator(bool val)	{ m_bIS_Emunerator = val; }
private:
	CC inline void setMatrOwner(bool val)		{ m_bMatrOwner = val; }
	CC inline bool isMatrOwner() const			{ return m_bMatrOwner; }

	const CMatrixData<T> *m_pMatrix;
	bool m_bMatrOwner;
	bool m_bIS_Emunerator;
	FILE *m_pFile;
};

Class2Def(CMatrixCanonChecker) : public Class2(CMatrixCol), public CCanonicityChecker<T>
{
public:
	CC CMatrixCanonChecker(const Class2(CMatrixData) *pMatrix, bool IS_enum, bool matrOwner = false) :
		Class2(CMatrixCol)(pMatrix, IS_enum, matrOwner),
		Class2(CCanonicityChecker)(pMatrix->rowNumb(), pMatrix->colNumb(), pMatrix->maxElement() + 1) {}

	CC CMatrixCanonChecker(Class2(CMatrixData) *pMatrix, S rowNumb, S colNumb, T maxElem, bool IS_enum) :
		Class2(CMatrixCol)(pMatrix, rowNumb, colNumb, maxElem, IS_enum),
		Class2(CCanonicityChecker)(rowNumb, colNumb, maxElem) {}
	CC ~CMatrixCanonChecker() {}
};

Class2Def(CMatrixCanonCheckerGPU) : public Class2(CMatrixCanonChecker)
{
public:
	CC CMatrixCanonCheckerGPU(const Class2(CMatrixData) *pMatrix, bool IS_enum, bool matrOwner = false) :
		Class2(CMatrixCanonChecker)(pMatrix, IS_enum, matrOwner) {}
	CC ~CMatrixCanonCheckerGPU() {}

};
