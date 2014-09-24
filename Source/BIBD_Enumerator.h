#include "InsSysEnumerator.h"

class CBIBD_Enumerator : public C_InSysEnumerator
{
public:
	CBIBD_Enumerator(const C_BIBD *pBIBD, bool matrOwner = false, bool noReplicatedBlocks = false) : C_InSysEnumerator(pBIBD, matrOwner, noReplicatedBlocks)   {}
	virtual bool makeJobTitle(char *buffer, int lenBuffer, const char *comment = "") const;
	virtual bool isTDesign_enumerator(size_t t)  const			{ return t <= 2; }

protected:
	virtual bool sortSolutions(CRowSolution *ptr, size_t idx);
	virtual int unforcedElement(const CColOrbit *p, int nRow) const;
	virtual bool solutionsForRightSideNeeded(const VECTOR_ELEMENT_TYPE *pRighPart, const VECTOR_ELEMENT_TYPE *pCurrSolution, size_t nRow) const;
	virtual bool makeFileName(char *buffer, int lenBuffer, const char *ext = NULL) const;
	virtual bool TestFeatures(CEnumInfo *pEnumInfo);
private:
	bool checkChoosenSolution(CRowSolution *pPrevSolution, size_t nRow, size_t usedSolIndex) const;
	virtual bool checkForcibleLambda(size_t fLambda)               { return fLambda == getInSys()->lambda(); }
};

