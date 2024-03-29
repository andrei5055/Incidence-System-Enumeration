#pragma once
#include "PBIBD_Enumerator.h"

Class2Def(CCombBIBD_Enumerator) : public Class2(CombEnumBase)
{
public:
	CCombBIBD_Enumerator(const InSysPntr pBIBD, uint enumFlags = t_enumDefault, int treadIdx = -1, uint nCanonChecker = 0) :
		Class2(CombEnumBase)(pBIBD, enumFlags, treadIdx, nCanonChecker) {
		setR(getInSys()->GetR(lenStabilizer())); 
		m_FirstPartSolutionIdx = new PERMUT_ELEMENT_TYPE[rowNumb()];
		const auto len = rowNumb() * numParts();
		m_bSolutionsWereConstructed = new unsigned char[len];
		memset(m_bSolutionsWereConstructed, 0, len);
	}
	~CCombBIBD_Enumerator();
	const auto *columnPermut() const					{ return m_pColumnPermut; }
protected:
	CK const char* getObjName() const override			{ return getInSys()->objectType() == t_objectType::t_Kirkman_Triple ? "k-System" : "CBIBD"; }
	CK const char* getTopLevelDirName() const override;
	CK RowSolutionPntr setFirstRowSolutions() override;
	CK void CreateForcedRows() override;
	CK T firtstNonfixedRowNumber() const override		{ return 3; }
	CC T lenStabilizer() const override					{ return 1; }
	CK bool createNewFile(const char* fName) const override   { return true; } // To overrite method for regular BIBDs
	CK VectorPntr paramSet(t_numbSetType idx) const override { 
		return (static_cast<Class2(CCombinedBIBD)*>(this->getInSys()))->paramSet(idx); 
	}
	CK T numLambdas() const override					{ return m_nNum_lambdas; }  // Used only in FindSolution which is called separatly for each part
	CK T getLambda(const VectorPntr pLambdaSet, T idx = 0, T numPart = 0) const override { return pLambdaSet->GetAt(numPart); }
#if !CONSTR_ON_GPU
	CK int addLambdaInfo(char *buffer, size_t lenBuffer, const char* pFrmt = NULL, size_t *pLambdaSetSize = NULL) const override;
	CK int getJobTitleInfo(char *buffer, int lenBuffer) const override;
	void getEnumerationObjectKey(char* pKey, int len) const override;
	char* getEnumerationObjectKeyA(char* pKey, int len, const char* pKeyIn = NULL) override;
	const char* getObjNameFormat() const override		{ return "  %14s:      "; }
#endif
	CK bool checkForcibleLambda(T fLambda, T nRows, T numPart) const override {
		if (blockIdx())
			return true;   // We don't need to check it for Kirkman Triple Systems

		const auto lambda = paramSet(t_lSet)->GetAt(numPart);
		return nRows == 2 ? fLambda == lambda : fLambda <= lambda;
	}
	CK CGroupOnParts<T>* makeGroupOnParts(const CCanonicityChecker* owner) override;
	CK void InitGroupOderStorage(const CGroupOnParts<T> *pGroupOnParts)  override {
		if (!m_pGroupOrder && pGroupOnParts)
			m_pGroupOrder = new CGroupOrder<T>;
		if (!m_pGroupOrders && pGroupOnParts)
			m_pGroupOrders = new size_t [pGroupOnParts->numGroups()];
	}
	CK MatrixDataPntr CreateSpareMatrix(const EnumeratorPntr pMaster) override;
	CK void CreateAuxiliaryStructures(EnumeratorPntr pMaster) override;
	CK int compareEnumerationDB_record(const char* record) override;
	CK void resetGroupOrder() override					{ m_pGroupOrder->setGroupOrder(1); }
	CK void incGroupOrder() override					{ m_pGroupOrder->setGroupOrder(m_pGroupOrder->groupOrder() + 1); }
	CK CGroupOrder<T>* extraGroupOrder() const override { return m_pGroupOrder; }
	CK void ConstructedDesignProcessing() const override;
	CK void beforeEnumInfoOutput() const override;
	CK void setDesignParams(designParam* pntr) override { 
		CEnumerator::setDesignParams(pntr); 
		m_bOutputMasters = pntr->outType & t_outputType::t_CombMasters;
	}

private:
	CK void setFirstPartSolutionIndex(PERMUT_ELEMENT_TYPE idx) override { *(m_FirstPartSolutionIdx + currentRowNumb()) = idx; }
	CK PERMUT_ELEMENT_TYPE firstPartSolutionIndex(T nRow) const override { return *(m_FirstPartSolutionIdx + nRow); }
	CK void CreateFirstRow();
	CK void createColumnPermut();
	CK bool outputMaster() const override				{ return m_bOutputMasters; }
	PERMUT_ELEMENT_TYPE* m_FirstPartSolutionIdx;
	size_t *m_pGroupOrders = NULL;				 // orders of group, acting on the parts with the same lambda
	CGroupOrder<T> *m_pGroupOrder = NULL;
	MatrixDataPntr m_pSpareMatrix = NULL;
	CMatrixCanonChecker* m_pCanonChecker = NULL; // used for cannonization of original matrix
	T* m_pColPermut = NULL;						 // Memory, which will keep current permutation of columns
	bool m_bColPermutOwner = false;
	T* m_pColumnPermut = NULL;					 // Permutation of columns (usially calculated by master and used by the threads
	bool m_bOutputMasters = false;
	const char* m_pAdjKeyPrefix;
	char* m_pAdjKey;
	size_t m_lenAdjKey;
	T m_nNum_lambdas = 1;
};

