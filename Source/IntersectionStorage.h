#pragma once
#include "DataTypes.h"
#include "Vector.h"

typedef CSimpleArray<MATRIX_ELEMENT_TYPE *> Ct_Storage;

class CPrevRowsIntersection
{
public:
    CPrevRowsIntersection()                                         { m_pRowIntersection = NULL; }
    ~CPrevRowsIntersection()										{ delete[] rowIntersectionPntr(); }

    void init(const size_t *pNumb, size_t len) {
		m_pNumb = pNumb;
		if (len)
			m_pRowIntersection = new PERMUT_ELEMENT_TYPE[len];
	}

	inline PERMUT_ELEMENT_TYPE *rowIntersectionPntr() const         { return m_pRowIntersection; }
    inline const size_t *numbIntersection() const					{ return m_pNumb; }
private:
    const size_t *m_pNumb;
	PERMUT_ELEMENT_TYPE *m_pRowIntersection;
};

template<class T>
class CIntersectionStorage
{
public:
	CIntersectionStorage(size_t t, size_t nRow, const CVector<T> *pLambdaSet) {
		m_pRowsIntersection = new CPrevRowsIntersection[nRow];
		const size_t len = t - 2;
		size_t *pNumb = new size_t[len * nRow];
		rowsIntersection(0)->init(pNumb, 0);
		memset(pNumb, 0, len * nRow * sizeof(*pNumb));

		size_t row = 1;
		pNumb += row * len;
		while (++row < nRow) {
			pNumb += len;
			size_t total = pLambdaSet->GetAt(0) * (pNumb[0] = row - 1);
			for (size_t i = 1; i < len; i++) {
				pNumb[i] = pNumb[i - 1] * (row - i - 1) / (i + 1);
				total += pLambdaSet->GetAt(i) * pNumb[i];
			}

			rowsIntersection(row)->init(pNumb, total);
		}
	}

    ~CIntersectionStorage() {
		delete[] rowsIntersection(0)->numbIntersection();
		delete[] rowsIntersection(0);
	}

    CPrevRowsIntersection *rowsIntersection(size_t idx) const      { return m_pRowsIntersection + idx; }
private:
    CPrevRowsIntersection *m_pRowsIntersection;
	Ct_Storage **m_pStorage;
};

template<class T>
class CIntersection
{
public:
	~CIntersection()						{ delete intersectionStorage(); }
protected:
	void InitIntersection(size_t t, size_t nRow, const CVector<T>* pLambdaSet) {
		m_pIntersectionStorage = new CIntersectionStorage<T>(t, nRow, pLambdaSet);
	}
	PERMUT_ELEMENT_TYPE *intersectionParam(const size_t** ppNumb, size_t row_numb) const {
		const auto* pPrev = intersectionStorage()->rowsIntersection(row_numb);
		*ppNumb = pPrev->numbIntersection();
		return pPrev->rowIntersectionPntr();
	}
	inline auto intersectionStorage() const { return m_pIntersectionStorage; }
private:
	CIntersectionStorage<T>* m_pIntersectionStorage = NULL;
};