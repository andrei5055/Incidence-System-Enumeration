#pragma once
#include "BIBD_Enumerator.h"

typedef struct {
	size_t nVar;
	CColOrbit<SIZE_TYPE> *pColOrb;
} OrbToVarMapping;

typedef CContainer<OrbToVarMapping *> COrbToVar;
template<class T> class CIntersectionStorage;

Class2Def(C_tDesignEnumerator) : public  Class2(CBIBD_Enumerator)
#if USE_EXRA_EQUATIONS
	, public CEquSystem, private COrbToVar
#endif
{
public:
	CK C_tDesignEnumerator(const Class2(C_tDesign) *pBIBD, uint enumFlags = t_enumDefault, int treadIdx = -1, uint nCanonChecker = 0);
	CK ~C_tDesignEnumerator();
#if !CONSTR_ON_GPU
	virtual bool makeJobTitle(const designParam *pParam, char *buffer, int lenBuffer, const char *comment = "") const;
#endif
	CK virtual bool isValidSolution(const VECTOR_ELEMENT_TYPE *pSol) const;
	PERMUT_ELEMENT_TYPE *getIntersectionParam(const size_t **ppNumb) const;
protected:
#if !CONSTR_ON_GPU
	virtual bool makeFileName(char *buffer, size_t lenBuffer, const char *ext = NULL) const;
#endif
	virtual CVariableMapping<T> *prepareCheckSolutions(size_t nVar);
	CK virtual void prepareToTestExtraFeatures();
	CK virtual void copyInfoFromMaster(const EnumeratorPntr pMaster);
	CK size_t numLambdas()const override						{ return 1; }
	CK virtual bool TestFeatures(EnumInfoPntr pEnumInfo, const MatrixDataPntr pMatrix, int *pMatrFlags = NULL, const EnumeratorPntr pEnum = NULL) const;
#if USE_EXRA_EQUATIONS
	CK virtual CEquSystem *equSystem()							{ return this; }
	CK virtual bool prepareToFindRowSolution()					{ equSystem()->resetArray(); return true; }
	CK virtual void setColOrbitForCurrentRow(CColOrbit *pColOrb) { m_pColOrbForCurrentRow = pColOrb; COrbToVar::resetArray(); }
	virtual CColOrbit *colOrbitForCurrentRow() const			{ return m_pColOrbForCurrentRow; }
	CK virtual void addColOrbitForVariable(S nVar, CColOrbit *pColOrb);
#endif
private:
	inline auto tDesign() const									{ return static_cast<TDesignPntr>(this->getInSys()); }
	inline auto intersectionStorage() const						{ return m_pIntersectionStorage; }
	virtual const char* getTopLevelDirName() const				{ return "t-designs"; }
	virtual void getEnumerationObjectKey(char* pInfo, int len) const { makeJobTitle(NULL, pInfo, len); }
	virtual void outputTitle(FILE *file) const;

#if USE_EXRA_EQUATIONS	
	CVariableMapping *constructExtraEquations(size_t t, size_t nVar);

	//	void setColOrbIdx(size_t idx)								{ m_nColOrbIdx = idx; }

	CColOrbit *m_pColOrbForCurrentRow;
	//	size_t m_nColOrbIdx;
	//	OrbToVarMapping *m_pOrbVarMapping;
#endif
	CIntersectionStorage<S> *m_pIntersectionStorage;
};

