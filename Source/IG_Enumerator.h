#pragma once
#include "PBIBD_Enumerator.h"

template<class T>
class CIncidenceStorage
{
public:
	CIncidenceStorage(T numb, T size) : m_Size(size)	{ m_pIncidence = new T[numb * size]; }
	~CIncidenceStorage()								{ delete[] m_pIncidence; }
	void addIncidence(T elem, T index, int numb)		{ m_pIncidence[m_Size * index + numb] = elem; }
private:
	T *m_pIncidence;
	const T m_Size;
};

template<class T>
class CIG_Enumerator : public CPBIBD_Enumerator<T>
{
public:
	CK CIG_Enumerator(const C_InSys<T> *pBIBD, const designParam *pParam, unsigned int enumFlags = t_enumDefault, bool firstPath = false, int treadIdx = -1, uint nCanonChecker = 0) :
		CPBIBD_Enumerator<T>(pBIBD, enumFlags, treadIdx, nCanonChecker), m_firstPath(firstPath) {
		const auto inSys = this->getInSys();
		const auto lambdaSet = inSys->GetNumSet(t_lSet);
		const auto nLambd = lambdaSet->GetSize();

		m_pRowIntersections = nLambd > 2 || lambdaSet->GetAt(0) || lambdaSet->GetAt(1) != 1?
							  new T[lenRowIntersection(rowNumb())] : NULL;

		// Allocate memory to store current lambdaA, lambdaB for blocks
		// This set of lambda's is not known, but it cannot have more than k + 1 elements

		const auto len = inSys->GetK() + 1;
		m_pLambda[0] = new T[2 * len + nLambd];
		m_pLambda[2] = (m_pLambda[1] = m_pLambda[0] + nLambd) + len;

		// Converting the values lamdaB into T format
		const auto &coeffB = pParam->lambdaB();
		for (auto i = nLambd; i--;)
			*(m_pLambda[0] + i) = coeffB[i];

		m_pElements = new CIncidenceStorage<T>(inSys->rowNumb(), inSys->GetR());
		m_pBlocks = new CIncidenceStorage<T>(inSys->colNumb(), inSys->GetK());
		setElementFlags(NULL);
	}

	CK ~CIG_Enumerator()								{ delete[] rowIntersections(); delete[] lambdaBSrc(); delete [] elementFlags();}
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
	inline void setNumRows(T nRow)						{ m_nNumbRows  = nRow; }
	inline T numRows() const							{ return m_nNumbRows; }
	CK void copyRowIntersection(const T *pntr)			{ 
		if (rowIntersections())
			memcpy(rowIntersections(), pntr, lenRowIntersection(rowNumb()) * sizeof(*pntr)); 
	}
	inline T *lambdaBSrc() const						{ return m_pLambda[0]; }
	inline bool firstPath() const						{ return m_firstPath; }
	inline size_t lenRowIntersection(T v) const			{ return v * (v - 1) / 2; }
	inline size_t nLambdas() const						{ return lambdaB() - lambdaA(); }
	inline size_t lenLambdas() const                    { return 2 * nLambdas() * sizeof(*lambdaA()); }
	bool CheckConstructedBlocks(T nRow, T k, T *pElementNumb);
	bool CheckTransitivityOnConstructedBlocks(T nRow, T k, T r, T *pElementNumb, uchar *pBlockFlags);
	bool CheckOrbits(const CPermutStorage<T> *permRowStorage, T *pRowOrbits = NULL) const;
	bool DefineInterstructForBlocks(size_t nColCurr, T k, const T *pElementNumb, T i, T *lambdaACurrCol) const;
	void FindAllElementsOfBlock(T nRow, size_t nColCurr, int j, T *pElementNumb) const;
	inline void setElementFlags(uchar *pntr)			{ m_pElementFlags = pntr; }
	inline uchar *elementFlags() const					{ return m_pElementFlags; }

	const bool m_firstPath;
	T *m_pRowIntersections;
	T *m_pLambda[3];
	T m_nNumbRows;							// Number of rows where all blocks associated with the first elements were constructed
	uchar *m_pElementFlags;					// Flags to mark the elements, which are already chosen
	CIncidenceStorage<T> *m_pElements;
	CIncidenceStorage<T> *m_pBlocks;
};
