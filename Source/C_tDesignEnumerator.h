#pragma once
#include "BIBD_Enumerator.h"
#include "Equation.h"

class CIntersectionStorage;

class C_tDesignEnumerator : public CBIBD_Enumerator
{
public:
	C_tDesignEnumerator(const C_tDesign *pBIBD, bool matrOwner = false, bool noReplicatedBlocks = false);
	~C_tDesignEnumerator();
	virtual bool makeJobTitle(char *buffer, int lenBuffer, const char *comment = "") const;
	virtual bool isValidSolution(const VECTOR_ELEMENT_TYPE *pSol) const;
	virtual bool isTDesign_enumerator(size_t t)  const			{ return t <= tDesign()->getT(); }
	PERMUT_ELEMENT_TYPE *getIntersectionParam(const size_t **ppNumb) const;
protected:
	virtual bool makeFileName(char *buffer, int lenBuffer, const char *ext = NULL) const;
	virtual void prepareCheckSolutions(size_t nVar);
    virtual void prepareToTestExtraFeatures();
	virtual void copyInfoFromMaster(const CEnumerator *pMaster);
    virtual bool constructing_t_Design()                        { return true; }
	virtual bool TestFeatures(CEnumInfo *pEnumInfo);
private:
	inline C_tDesign *tDesign() const                           { return static_cast<C_tDesign *>(getInSys()); }
    inline CIntersectionStorage *intersectionStorage() const    { return m_pIntersectionStorage; }
	inline CVariable **varPntr() const							{ return m_pVarPntr; }
	inline CEquation **equPntr() const							{ return m_pEquPntr; }
	inline VECTOR_ELEMENT_TYPE *equIndx() const					{ return m_pEquIndx; }
	inline void resetEquNumb()									{ setEquNumb(0); }
	inline void setEquNumbMax(size_t n)							{ m_nEquNumbMax = n; }
	inline void setEquNumb(size_t val)							{ m_nEquNumb = val; }
	inline size_t equNumb() const								{ return m_nEquNumb; }
	inline size_t equNumbMax() const							{ return m_nEquNumbMax; }
	inline CEquation *equation(size_t idx) const				{ return *(equPntr() + idx); }
	inline CEquation *equation(const CVariable *pVar) const		{ return *(equPntr() + *(equIndx() + pVar->equIdx())); }
	void setEquation(CEquation *pEqu, size_t idx);
	void setEquation(CEquation *pEquFrom, const CEquation *pEquTo);
	void releaseVariables();
	void addEquation(CVariable *pFirstVar, size_t nVar, size_t lambda);
	void allocateMemoryForEquations(size_t nEqu);
	void constructAdditionalEquations(size_t t, size_t nVar);
	int solveAdditionalEquations();

    CIntersectionStorage *m_pIntersectionStorage;
	CVariable **m_pVarPntr;
	CEquation **m_pEquPntr;
	VECTOR_ELEMENT_TYPE *m_pEquIndx;
	size_t m_nEquNumb;
	size_t m_nEquNumbMax;
};

