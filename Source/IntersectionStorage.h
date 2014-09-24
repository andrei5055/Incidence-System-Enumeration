#include "DataTypes.h"

typedef CSimpleArray<MATRIX_ELEMENT_TYPE *> Ct_Storage;
class CVector;

class CPrevRowsIntersection
{
public:
    CPrevRowsIntersection()                                         { m_pRowIntersection = NULL; }
    ~CPrevRowsIntersection();
    void init(const size_t *pNumb, size_t len);
	inline PERMUT_ELEMENT_TYPE *rowIntersectionPntr() const         { return m_pRowIntersection; }
    inline const size_t *numbIntersection() const					{ return m_pNumb; }
private:
    const size_t *m_pNumb;
	PERMUT_ELEMENT_TYPE *m_pRowIntersection;
};

class CIntersectionStorage
{
public:
	CIntersectionStorage(size_t t, size_t nRow, const CVector *pLambdaSet);
    ~CIntersectionStorage();
    CPrevRowsIntersection *rowsIntersection(size_t idx) const      { return m_pRowsIntersection + idx; }
private:
    CPrevRowsIntersection *m_pRowsIntersection;
	Ct_Storage **m_pStorage;
};