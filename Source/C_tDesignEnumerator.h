#pragma once
#include "BIBD_Enumerator.h"

typedef struct {
	size_t nVar;
	CColOrbit<SIZE_TYPE> *pColOrb;
} OrbToVarMapping;

typedef CContainer<OrbToVarMapping *> COrbToVar;


Class2Def(C_tDesignEnumerator) : public Class2(CBIBD_Enumerator), public Class2(CIntersection)
#if USE_EXRA_EQUATIONS
	, public CEquSystem, private COrbToVar
#endif
{
public:
	CK C_tDesignEnumerator(const InSysPntr pBIBD, uint enumFlags = t_enumDefault, int treadIdx = -1, uint nCanonChecker = 0);
	CK ~C_tDesignEnumerator();
#if !CONSTR_ON_GPU
	virtual void makeJobTitle(const designParam *pParam, char *buffer, int lenBuffer, const char *comment = "") const;
#endif
	CK virtual bool isValidSolution(const VECTOR_ELEMENT_TYPE *pSol) const;
	T *getIntersectionParam(const T **ppNumb) const {
		return intersectionParam(ppNumb, this->currentRowNumb());
	}
protected:
#if !CONSTR_ON_GPU
	virtual bool makeFileName(char *buffer, size_t lenBuffer, const char *ext = NULL) const;
#endif
	virtual VariableMappingPntr prepareCheckSolutions(size_t nVar) {
		return prepareRowIntersections(this->matrix(), this->currentRowNumb(), tDesign()->GetNumSet(t_lSet)->GetAt(0), tDesign()->getT());
	}
	CK virtual void prepareToTestExtraFeatures() {
		InitIntersection(tDesign()->getT(), this->rowNumb(), tDesign()->GetNumSet(t_lSet));
	}
	CK virtual void copyInfoFromMaster(const EnumeratorPntr pMaster);
	CK T numLambdas() const override							{ return 1; }
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
	virtual const char* getTopLevelDirName() const				{ return "t-designs"; }
	virtual void getEnumerationObjectKey(char* pInfo, int len) const { makeJobTitle(NULL, pInfo, len); }
	virtual const char* getObjNameFormat() const				{ return "%12s:        "; }

#if USE_EXRA_EQUATIONS	
	CVariableMapping *constructExtraEquations(size_t t, size_t nVar);

	//	void setColOrbIdx(size_t idx)								{ m_nColOrbIdx = idx; }

	CColOrbit *m_pColOrbForCurrentRow;
	//	size_t m_nColOrbIdx;
	//	OrbToVarMapping *m_pOrbVarMapping;
#endif
};

