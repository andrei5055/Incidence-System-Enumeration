#include "IntersectionStorage.h"
#include "Vector.h"

CPrevRowsIntersection::~CPrevRowsIntersection()
{
    delete [] rowIntersectionPntr();
}

void CPrevRowsIntersection::init(const size_t *pNumb, size_t len)
{
    m_pNumb = pNumb;
    if (len)
		m_pRowIntersection = new PERMUT_ELEMENT_TYPE[len];
}

CIntersectionStorage::CIntersectionStorage(size_t t, size_t nRow, const CVector *pLambdaSet)
{
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

CIntersectionStorage::~CIntersectionStorage()
{
    delete [] rowsIntersection(0)->numbIntersection();
    delete [] rowsIntersection(0);
}