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
	CK CIG_Enumerator(const C_InSys<T> *pBIBD, const designParam *pParam, unsigned int enumFlags = t_enumDefault, bool firstPath = false, int treadIdx = -1, uint nCanonChecker = 0);
	CK ~CIG_Enumerator();
	CK void CloneMasterInfo(const CEnumerator<T> *p, size_t nRow) override;
	CK inline T *rowIntersections() const				{ return m_pRowIntersections; }
	inline T *lambdaA() const							{ return m_pLambda[1]; }
	inline T *lambdaB() const							{ return m_pLambda[2]; }
	inline int *areaWeight() const						{ return m_pAreaWeight; }
protected:
	CK bool TestFeatures(CEnumInfo<T> *pEnumInfo, const CMatrixData<T> *pMatrix, 
		       int *pMatrFlags = NULL, CEnumerator<T> *pEnum = NULL) const override;
	bool makeFileName(char *buffer, size_t lenBuffer, const char *ext) const override;
	CK const char *getObjName() const override			{ return "I-Graph"; }
	CK bool fileExists(const char *path, bool file = true) const override;
	CK bool createNewFile(const char *fName) const override;
	CK bool SeekLogFile() const override				{ return true; }

	CK bool prepareToFindRowSolution() override;
	virtual void reset(T nRow);
private:
	inline void setNumRows(T nRow)						{ m_nNumbRows  = nRow; }
	inline T numRows() const							{ return m_nNumbRows; }
	CK void copyRowIntersection(const T *pntr);
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
	inline void setPermCols(T *pntr)					{ m_pPrmCols = pntr; }
	inline T *permCols() const							{ return m_pPrmCols; }
	inline void setAreaWeight(int *pntr)				{ m_pAreaWeight = pntr; }
	inline bool checkUnitsInThreeAreas() const			{ return areaWeight() != NULL; }
	inline bool strongCheckUnitsInThreeAreas() const	{ return m_bStrongCheck; }
	inline void setStrongCheckUnitsInThreeAreas(bool v) { m_bStrongCheck = v; }
	bool checkThreeAreasUnits(int r, int k, T nrow) const;
	virtual const char* getTopLevelDirName() const		{ return "Semi-symmetric_Graphs"; }

	const bool m_firstPath;
	T *m_pRowIntersections;
	T *m_pLambda[3];
	T m_nNumbRows;							// Number of rows where all blocks associated with the first elements were constructed
	T *m_pPrmCols;							// Memory to keep permitations of columns
	uchar *m_pElementFlags;					// Flags to mark the elements, which are already chosen
	CIncidenceStorage<T> *m_pElements;
	CIncidenceStorage<T> *m_pBlocks;
	int *m_pAreaWeight;
	bool m_bStrongCheck;
};
