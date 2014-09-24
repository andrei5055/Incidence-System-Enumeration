#include "Enumerator.h"

class CRightPartFilter;

class C_InSysEnumerator : public CEnumerator, public CInSysSolver, public CVector
{
public:
	C_InSysEnumerator(const C_InSys *pInSysm, bool matrOwner = false, bool noReplicatedBlocks = false);
	~C_InSysEnumerator();
	virtual VECTOR_ELEMENT_TYPE getX0_3() const                 { return m_x0_3; }
	virtual size_t firstUnforcedRow() const                     { return m_firstUnforcedRow; }
	virtual void setFirstUnforcedRow(size_t rowNum = 0)         { m_firstUnforcedRow = rowNum; }
	virtual size_t *forcibleLambdaPntr() const                  { return m_pForsibleLambda; }
	virtual bool isIS_enumerator() const                        { return true; }
	virtual bool isTDesign_enumerator(size_t t) const			{ return false; }
protected:
	virtual void setX0_3(VECTOR_ELEMENT_TYPE value)             { m_x0_3 = value; }
	virtual bool sortSolutions(CRowSolution *ptr, size_t idx);
	inline CInSysRowEquation *inSysRowEquation() const          { return (CInSysRowEquation *)rowEquation(); }
	virtual bool solutionsForRightSideNeeded(const VECTOR_ELEMENT_TYPE *pRighPart, const VECTOR_ELEMENT_TYPE *pCurrSolution, size_t nRow) const
																{ return true; }
	virtual bool noReplicatedBlocks() const						{ return m_bNoReplBlock; }
	CColOrbit **unforcedOrbits(size_t nRow)	const				{ return getUnforceColOrbPntr() + rank() * nRow; }
private:
	void addForciblyConstructedColOrbit(CColOrbit *pColOrbit, CColOrbit *pPrev, int idx);
	virtual CRowSolution *setFirstRowSolutions();
	virtual size_t MakeSystem();
	virtual CRowSolution *FindSolution(size_t nVar, PERMUT_ELEMENT_TYPE lastRightPartIndex = -1);
	void setVariableLimit(size_t nVar, size_t len, size_t nRowToBuild, size_t k, int colWeight);
	virtual CColOrbit **getUnforceColOrbPntr() const            { return forcibleLambda(currentRowNumb()) ? unforcedColOrbPntr() : NULL; }
	virtual bool checkForcibleLambda(size_t fLambda)            { return true; }
	virtual void resetFirstUnforcedRow();
	virtual void prepareCheckSolutions(size_t nVar)				{}
    virtual bool constructing_t_Design()                        { return false; }
	inline CRightPartFilter *rightPartFilter()                  { return m_pRightPartFilter; }
	inline void setForcibleLambdaPntr(size_t *p)                { m_pForsibleLambda = p; }
	inline void setForcibleLambda(size_t i, size_t v)           { *(forcibleLambdaPntr() + i) = v; }
	inline size_t forcibleLambda(size_t i) const                { return *(forcibleLambdaPntr() + i); }

	VECTOR_ELEMENT_TYPE m_x0_3;
	size_t *m_pForsibleLambda;
	size_t m_firstUnforcedRow;
	CRightPartFilter *m_pRightPartFilter;
	const bool m_bNoReplBlock;
};

