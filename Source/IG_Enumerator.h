#pragma once
#include "PBIBD_Enumerator.h"

Class2Def(CIncidenceStorage)
{
public:
	CIncidenceStorage(S numb, S size) : m_Size(size)	{ m_pIncidence = new T[numb * size]; }
	~CIncidenceStorage()								{ delete[] m_pIncidence; }
	void addIncidence(S elem, S index, int numb)		{ m_pIncidence[m_Size * index + numb] = elem; }
private:
	T *m_pIncidence;
	const S m_Size;
};

Class2Def(CIG_Enumerator) : public Class2(CPBIBD_Enumerator)
{
public:
	CK CIG_Enumerator(const InSysPntr pBIBD, const designParam *pParam, unsigned int enumFlags = t_enumDefault, bool firstPath = false, int treadIdx = -1, uint nCanonChecker = 0);
	CK ~CIG_Enumerator();
	CK void CloneMasterInfo(const EnumeratorPntr p, size_t nRow) override;
	CK inline auto rowIntersections() const				{ return m_pRowIntersections; }
	inline auto lambdaA() const							{ return m_pLambda[1]; }
	inline auto lambdaB() const							{ return m_pLambda[2]; }
	inline auto areaWeight() const						{ return m_pAreaWeight; }
protected:
	CK bool TestFeatures(EnumInfoPntr pEnumInfo, const MatrixDataPntr pMatrix,
		       int *pMatrFlags = NULL, EnumeratorPntr pEnum = NULL) const override;
	bool makeFileName(char *buffer, size_t lenBuffer, const char *ext) const override;
	CK const char *getObjName() const override			{ return "I-Graph"; }
	CK bool fileExists(const char *path, bool file = true) const override;
	CK bool createNewFile(const char *fName) const override;
	CK bool SeekLogFile() const override				{ return true; }

	CK bool prepareToFindRowSolution() override;
	virtual void reset(S nRow);
private:
	inline void setNumRows(S nRow)						{ m_nNumbRows  = nRow; }
	inline auto numRows() const							{ return m_nNumbRows; }
	CK void copyRowIntersection(const S *pntr);
	inline S *lambdaBSrc() const						{ return m_pLambda[0]; }
	inline bool firstPath() const						{ return m_firstPath; }
	inline size_t lenRowIntersection(S v) const			{ return v * (v - 1) / 2; }
	inline size_t nLambdas() const						{ return lambdaB() - lambdaA(); }
	inline size_t lenLambdas() const                    { return 2 * nLambdas() * sizeof(*lambdaA()); }
	bool CheckConstructedBlocks(S nRow, S k, S *pElementNumb);
	bool CheckTransitivityOnConstructedBlocks(S nRow, S k, S r, S *pElementNumb, uchar *pBlockFlags);
	bool CheckOrbits(const PermutStoragePntr permRowStorage, S *pRowOrbits = NULL) const;
	bool DefineInterstructForBlocks(size_t nColCurr, S k, const S *pElementNumb, S i, S *lambdaACurrCol) const;
	void FindAllElementsOfBlock(S nRow, size_t nColCurr, int j, S *pElementNumb) const;
	inline void setElementFlags(uchar *pntr)			{ m_pElementFlags = pntr; }
	inline uchar *elementFlags() const					{ return m_pElementFlags; }
	inline void setPermCols(S *pntr)					{ m_pPrmCols = pntr; }
	inline auto permCols() const						{ return m_pPrmCols; }
	inline void setAreaWeight(uint *pntr)				{ m_pAreaWeight = pntr; }
	inline bool checkUnitsInThreeAreas() const			{ return areaWeight() != NULL; }
	inline bool strongCheckUnitsInThreeAreas() const	{ return m_bStrongCheck; }
	inline void setStrongCheckUnitsInThreeAreas(bool v) { m_bStrongCheck = v; }
	bool checkThreeAreasUnits(S r, S k, S nrow) const;
	virtual const char* getTopLevelDirName() const		{ return "Semi-symmetric_Graphs"; }

	const bool m_firstPath;
	S *m_pRowIntersections;
	S *m_pLambda[3];
	S m_nNumbRows;							// Number of rows where all blocks associated with the first elements were constructed
	S *m_pPrmCols;							// Memory to keep permitations of columns
	uchar *m_pElementFlags;					// Flags to mark the elements, which are already chosen
	Class2(CIncidenceStorage) *m_pElements;
	Class2(CIncidenceStorage) *m_pBlocks;
	uint *m_pAreaWeight;
	bool m_bStrongCheck;
};
