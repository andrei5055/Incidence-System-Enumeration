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

template<class T>
class CMatrixCanonChecker : public CMatrixCol<T>, public CCanonicityChecker<T>
{
public:
	CC CMatrixCanonChecker(const CMatrixData<T> *pMatrix, bool IS_enum, bool matrOwner = false) :
		CMatrixCol<T>(pMatrix, IS_enum, matrOwner),
		CCanonicityChecker<T>(pMatrix->rowNumb(), pMatrix->colNumb(), pMatrix->maxElement() + 1) {}

	CC CMatrixCanonChecker(CMatrixData<T> *pMatrix, T rowNumb, T colNumb, T maxElem, bool IS_enum) :
		CMatrixCol<T>(pMatrix, rowNumb, colNumb, maxElem, IS_enum),
		CCanonicityChecker<T>(rowNumb, colNumb, maxElem) {}
	CC ~CMatrixCanonChecker() {}
};

template<class T>
class CMatrixCanonCheckerGPU : public CMatrixCanonChecker<T>
{
public:
	CC CMatrixCanonCheckerGPU(const CMatrixData<T> *pMatrix, bool IS_enum, bool matrOwner = false) :
		CMatrixCanonChecker<T>(pMatrix, IS_enum, matrOwner) {}
	CC ~CMatrixCanonCheckerGPU() {}

};
