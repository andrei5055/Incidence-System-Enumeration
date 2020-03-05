#pragma once
#include "ColOrbitManager.h"
#include "matrix.h"

Class2Def(CMatrixCol) : public CColOrbitManager<S>
{
public:
	CC CMatrixCol(const Class2(CMatrixData) *pMatrix, bool IS_enum, bool matrOwner) :
		CColOrbitManager(pMatrix->maxElement() + 1, pMatrix->rowNumb(), pMatrix->colNumb()),
		m_pMatrix(pMatrix), m_bIS_Emunerator(IS_enum), m_bMatrOwner(matrOwner) {
		setOutFile(NULL);
	}
	CC CMatrixCol(Class2(CMatrixData) *pMatrix, S rowNumb, S colNumb, T maxElem, bool IS_enum) :
		CColOrbitManager(maxElem + 1, rowNumb, colNumb),
		m_pMatrix(pMatrix), m_bIS_Emunerator(IS_enum), m_bMatrOwner(false) {
		setOutFile(NULL);
	}
	CC ~CMatrixCol() {
		if (isMatrOwner())???
			delete static_cast<const Class2(CMatrix) *>(matrix());
	}
	CC inline const auto matrix() const			{ return m_pMatrix; }
	inline auto IS_enumerator() const			{ return m_bIS_Emunerator; }
	inline void closeFile()						{ if (outFile()) { fclose(outFile()); setOutFile(NULL); }}
	CC inline void setOutFile(FILE *file)		{ m_pFile = file; }
	inline auto outFile() const					{ return m_pFile; }
	inline auto outFilePntr()					{ return &m_pFile; }
protected:
private:
	CC inline auto isMatrOwner() const			{ return m_bMatrOwner; }

	const Class2(CMatrixData) *m_pMatrix;
	const bool m_bMatrOwner;
	const bool m_bIS_Emunerator;
	FILE *m_pFile;
};

Class2Def(CMatrixCanonChecker) : public Class2(CMatrixCol), public Class2(CCanonicityChecker)
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
