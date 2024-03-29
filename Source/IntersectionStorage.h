#pragma once
#include "DataTypes.h"
#include "Vector.h"
#include "matrix.h"

typedef CSimpleArray<MATRIX_ELEMENT_TYPE *> Ct_Storage;

template<class T>
class CPrevRowsIntersection
{
public:
    ~CPrevRowsIntersection()										{ delete[] rowIntersectionPntr(); }

    void init(const T *pNumb, size_t len) {
		m_pNumb = pNumb;
		if (len)
			m_pRowIntersection = new T[len];
	}

	inline T *rowIntersectionPntr() const							{ return m_pRowIntersection; }
    inline const T *numbIntersection() const						{ return m_pNumb; }
private:
    const T *m_pNumb = NULL;
	T *m_pRowIntersection = NULL;
};

template<class T>
class CIntersectionStorage
{
public:
	CIntersectionStorage(T t, T nRow, const CVector<T> *pLambdaSet) {
		m_pRowsIntersection = new CPrevRowsIntersection<T>[nRow];
		const auto len = t - 2;
		T *pNumb = new T[len * nRow];
		// pNumb[0] - number of common blocks for pairs of elements * all previous pairs 
		// pNumb[1] - number of common blocks for triples of elements * all previous triples etc. 
		rowsIntersection(0)->init(pNumb, 0);
		memset(pNumb, 0, len * nRow * sizeof(*pNumb));

		if (pLambdaSet->GetSize() > 1) {
			// For t-designs
			T row = 1;
			pNumb += row * len;
			while (++row < nRow) {
				pNumb += len;
				size_t total = pLambdaSet->GetAt(0) * (pNumb[0] = row - 1);
				for (T i = 1; i < len; i++) {
					pNumb[i] = pNumb[i - 1] * (row - i - 1) / (i + 1);
					total += pLambdaSet->GetAt(i) * pNumb[i];
				}

				rowsIntersection(row)->init(pNumb, total);
			}
		}
		else {
			// For symmetrical design with t_symmetrical_t_cond flag activated
			const auto λ = pLambdaSet->GetAt(0);
			T row = 1;
			pNumb += row * len;
			while (++row < nRow) {
				pNumb += len;
				size_t total = (pNumb[0] = row - 1);
				for (T i = 1; i < len; i++) {
					pNumb[i] = pNumb[i - 1] * (row - i - 1) / (i + 1);
					total += pNumb[i];
				}

				rowsIntersection(row)->init(pNumb, λ * total);
			}
		}
	}

    ~CIntersectionStorage() {
		delete[] rowsIntersection(0)->numbIntersection();
		delete[] rowsIntersection(0);
	}

    CPrevRowsIntersection<T> *rowsIntersection(size_t idx) const      { return m_pRowsIntersection + idx; }
private:
    CPrevRowsIntersection<T> *m_pRowsIntersection;
	Ct_Storage **m_pStorage;
};


Class2Def(CIntersection)
{
public:
	CIntersection()							{}
	CIntersection(T t, T nRows, const CVector<S>* pLambdaSet) {
		InitIntersection(t, nRows, pLambdaSet);
	}
	~CIntersection()						{ delete intersectionStorage(); }
	VariableMappingPntr prepareRowIntersections(const InSysPntr pMatrix, T currRowNumb, T lambda, T t) const;
	T* intersectionParam(const T** ppNumb, size_t row_numb) const {
		const auto* pPrev = intersectionStorage()->rowsIntersection(row_numb);
		*ppNumb = pPrev->numbIntersection();
		return pPrev->rowIntersectionPntr();
	}
protected:
	void InitIntersection(T t, T nRows, const CVector<S>* pLambdaSet) {
		m_pIntersectionStorage = new CIntersectionStorage<T>(t, nRows, pLambdaSet);
	}
	inline auto intersectionStorage() const { return m_pIntersectionStorage; }
private:
	CIntersectionStorage<T>* m_pIntersectionStorage = NULL;
};