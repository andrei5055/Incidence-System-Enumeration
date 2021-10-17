#pragma once
#include "BIBD_Enumerator.h"

Class2Def(CCombBIBD_Enumerator) : public Class2(CBIBD_Enumerator)
{
public:
	CCombBIBD_Enumerator(const InSysPntr pBIBD, uint enumFlags = t_enumDefault, int treadIdx = -1, uint nCanonChecker = 0) :
		Class2(CBIBD_Enumerator)(pBIBD, enumFlags, treadIdx, nCanonChecker) {
		setR(getInSys()->GetR(lenStabilizer()));
		m_FirstPartSolutionIdx = new PERMUT_ELEMENT_TYPE[rowNumb()];
		const auto len = rowNumb() * numParts();
		m_bSolutionsWereConstructed = new unsigned char[len];
		memset(m_bSolutionsWereConstructed, 0, len);
	}
	~CCombBIBD_Enumerator()	{
		delete[] m_FirstPartSolutionIdx;
		delete[] m_bSolutionsWereConstructed;
	}
protected:
	CK virtual const char* getObjName() const		{ return "CBIBD"; }
	virtual const char* getTopLevelDirName() const	{ return "Combined_BIBDs"; }
	CK virtual RowSolutionPntr setFirstRowSolutions();
	CK virtual void CreateForcedRows();
	CK virtual T firtstNonfixedRowNumber() const	{ return 3; }
	CC virtual T lenStabilizer() const              { return 1; }
	CK virtual VectorPntr paramSet(t_numbSetType idx) const	{ return (static_cast<Class2(CCombinedBIBD)*>(this->getInSys()))->paramSet(idx); }
	CK virtual size_t numLambdas()					{ return 1; }
	CK virtual T getLambda(const VectorPntr pLambdaSet, T idx = 0, T numPart = 0) const { return pLambdaSet->GetAt(numPart); }
#if !CONSTR_ON_GPU
	CK virtual int addLambdaInfo(char *buffer, size_t lenBuffer, const char* pFrmt = NULL, int* pLambdaSetSize = NULL) const;
	CK virtual int getJobTitleInfo(char *buffer, int lenBuffer) const;
	virtual void getEnumerationObjectKey(char* pKey, int len) const;
	virtual char* getEnumerationObjectKeyA(char* pKey, int len) const;
	virtual const char* getObjNameFormat() const	{ return "  %14s:      "; }
#endif
	CK virtual bool checkForcibleLambda(T fLambda, T nRows, T numPart) const {
		const auto lambda = paramSet(t_lSet)->GetAt(numPart);
		return nRows == 2 ? fLambda == lambda : fLambda <= lambda;
	}
	CK virtual CGroupOnParts<T> * makeGroupOnParts(const EnumeratorPntr owner) {
		auto lanbdaSet = paramSet(t_lSet);
		auto jMax = lanbdaSet->GetSize() - 1;
		CVector<T> lengths;
		T prevLambda = 0;
		uint count;
		uint factorial;
		for (int j = 0; j <= jMax; j++) {
			const auto lambda = lanbdaSet->GetAt(j);
			if (prevLambda == lambda) {
				factorial *= (++count);
				if (j < jMax)
					continue;
				j++;
			}

			if (prevLambda) {
				if (factorial > 1) {
					lengths.AddElement(j-count); 	// index of the first BIBD with the same lambda
					lengths.AddElement(count);		// number of BIBDs with the same lambda
					lengths.AddElement(factorial);  // group order
				}
			}
			prevLambda = lambda;
			factorial = count = 1;
		}
		if (!lengths.GetSize())
			return NULL;

		updateCanonicityChecker(rowNumb(), colNumb());
		return  new CGroupOnParts<T>(owner, lengths, 3);
	}
	CK MatrixDataPntr CreateSpareMatrix(const MatrixDataPntr pMatr);
private:
	CK virtual void setFirstPartSolutionIndex(PERMUT_ELEMENT_TYPE idx)	{ *(m_FirstPartSolutionIdx + currentRowNumb()) = idx; }
	CK virtual PERMUT_ELEMENT_TYPE firstPartSolutionIndex(T nRow) const	{ return *(m_FirstPartSolutionIdx + nRow); }
	void CreateFirstRow(S* pFirstRow=NULL);
	PERMUT_ELEMENT_TYPE* m_FirstPartSolutionIdx;
};

