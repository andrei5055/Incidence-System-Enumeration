#pragma once
#include "BIBD_Enumerator.h"

typedef struct {
	size_t nVar;
	CColOrbit<MATRIX_ELEMENT_TYPE> *pColOrb;
} OrbToVarMapping;

typedef CContainer<OrbToVarMapping *> COrbToVar;
template<class T> class CIntersectionStorage;
//template<class T> class C_tDesign;
template<class T> class CEnumerator;

template<class T>
class C_tDesignEnumerator : public CBIBD_Enumerator<T>
#if USE_EXRA_EQUATIONS
	, public CEquSystem, private COrbToVar
#endif
{
public:
	CK C_tDesignEnumerator(const C_tDesign<T> *pBIBD, bool matrOwner = false, bool noReplicatedBlocks = false, int treadIdx = -1, uint nCanonChecker = 0);
	CK ~C_tDesignEnumerator();
#if !CONSTR_ON_GPU
	virtual bool makeJobTitle(char *buffer, int lenBuffer, const char *comment = "") const;
#endif
	CK virtual bool isValidSolution(const VECTOR_ELEMENT_TYPE *pSol) const;
	CK virtual bool isTDesign_enumerator(size_t t) const { return t <= tDesign()->getT(); }
	PERMUT_ELEMENT_TYPE *getIntersectionParam(const size_t **ppNumb) const;
protected:
#if !CONSTR_ON_GPU
	virtual bool makeFileName(char *buffer, size_t lenBuffer, const char *ext = NULL) const;
#endif
	virtual CVariableMapping<T> *prepareCheckSolutions(size_t nVar);
	CK virtual void prepareToTestExtraFeatures();
	CK virtual void copyInfoFromMaster(const CEnumerator<T> *pMaster);
	CK virtual bool constructing_t_Design() { return true; }
	CK virtual bool TestFeatures(CEnumInfo<T> *pEnumInfo, const CMatrixData<T> *pMatrix, int *pMatrFlags = NULL) const;
#if USE_EXRA_EQUATIONS
	CK virtual CEquSystem *equSystem() { return this; }
	CK virtual void prepareToFindRowSolution() { equSystem()->resetArray(); }
	CK virtual void setColOrbitForCurrentRow(CColOrbit *pColOrb) { m_pColOrbForCurrentRow = pColOrb; COrbToVar::resetArray(); }
	virtual CColOrbit *colOrbitForCurrentRow() const { return m_pColOrbForCurrentRow; }
	CK virtual void addColOrbitForVariable(size_t nVar, CColOrbit *pColOrb);
#endif
private:
	inline C_tDesign<T> *tDesign() const { return static_cast<C_tDesign<T> *>(this->getInSys()); }
	inline CIntersectionStorage<T> *intersectionStorage() const { return m_pIntersectionStorage; }

#if USE_EXRA_EQUATIONS	
	CVariableMapping *constructExtraEquations(size_t t, size_t nVar);

	//	void setColOrbIdx(size_t idx)								{ m_nColOrbIdx = idx; }

	CColOrbit *m_pColOrbForCurrentRow;
	//	size_t m_nColOrbIdx;
	//	OrbToVarMapping *m_pOrbVarMapping;
#endif
	CIntersectionStorage<T> *m_pIntersectionStorage;
};

