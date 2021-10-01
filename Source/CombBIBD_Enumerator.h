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
	CK virtual S firtstNonfixedRowNumber() const	{ return 3; }
	CC virtual S lenStabilizer() const              { return 1; }
	CK virtual VectorPntr paramSet(t_numbSetType idx) const	{ return (static_cast<Class2(CCombinedBIBD)*>(this->getInSys()))->paramSet(idx); }
	CK virtual size_t numLambdas()					{ return 1; }
	CK virtual S getLambda(const VectorPntr pLambdaSet, S idx = 0, S numPart = 0) const { return pLambdaSet->GetAt(numPart); }
#if !CONSTR_ON_GPU
	CK virtual int addLambdaInfo(char *buffer, size_t lenBuffer, const char* pFrmt = NULL, int* pLambdaSetSize = NULL) const;
	CK virtual int getJobTitleInfo(char *buffer, int lenBuffer) const;
	virtual void getEnumerationObjectKey(char* pInfo, int len) const;
	virtual const char* getObjNameFormat() const	{ return "  %14s:      "; }
#endif
	CK virtual bool checkForcibleLambda(S fLambda, S nRows, S numPart) const {
		const auto lambda = paramSet(t_lSet)->GetAt(numPart);
		return nRows == 2 ? fLambda == lambda : fLambda <= lambda;
	}
	CK virtual CGroupOnParts<uint> * makeGroupOnParts(const EnumeratorPntr owner) const {
		const auto aaa = this->getInSys()->GetNumSet(t_lSet);
		auto b = paramSet(t_lSet);
		auto jMax = b->GetSize() - 1;
		CVector<uint> lengths;
		S prev = 0;
		uint count;
		uint factorial;
		size_t numNontrivialGroups = 0;
		for (int j = 0; j <= jMax; j++) {
			const auto lambda = b->GetAt(j);
			if (prev == lambda) {
				factorial *= (++count);
				if (j < jMax)
					continue;
			}

			if (prev) {
				if (factorial > 1)
					numNontrivialGroups++;

				lengths.AddElement(count);
				lengths.AddElement(factorial);
			}
			prev = lambda;
			factorial = count = 1;
		}

		return numNontrivialGroups ? new CGroupOnParts<uint>(owner, lengths, numNontrivialGroups) : NULL;
	}
private:
	CK virtual void setFirstPartSolutionIndex(PERMUT_ELEMENT_TYPE idx)	{ *(m_FirstPartSolutionIdx + currentRowNumb()) = idx; }
	CK virtual PERMUT_ELEMENT_TYPE firstPartSolutionIndex(S nRow) const	{ return *(m_FirstPartSolutionIdx + nRow); }
	PERMUT_ELEMENT_TYPE* m_FirstPartSolutionIdx;
};

