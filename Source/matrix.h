#include "ColOrbits.h"

class CCanonicityChecker;

class CMatrix : public CArrayOfMatrixElements 
{
 public:
 	CMatrix(size_t nRows, size_t nCols = 0) : CArrayOfMatrixElements() { Init(nRows, nCols); }
    CMatrix(const CMatrix *pMatrix, const int *pPermRow, const int *pPermCol);
	virtual ~CMatrix()										{}
	inline size_t colNumb() const							{ return m_nCols; }
	inline size_t rowNumb() const							{ return m_nRows; }
	inline MATRIX_ELEMENT_TYPE *GetRow(size_t nRow) const	{ return (MATRIX_ELEMENT_TYPE *)m_pData + nRow * m_nCols; }
	virtual int maxElement() const							{ return 1; }
	void printOut(FILE *pFile = NULL, size_t nRow = -1, ulonglong matrNumber = -1, const CCanonicityChecker *pCanonCheck = NULL) const;
 protected:
    void Init(size_t nRows, size_t nCols);
    
	size_t m_nRows;
	size_t m_nCols;
};

typedef enum {
        t_rSet,
        t_kSet,
        t_lSet
} t_numbSetType;

class C_InSys : public CMatrix
{
 public:
	C_InSys(int nRows, int nCols);
	C_InSys(const C_InSys *pMaster, size_t nRow);
	~C_InSys();
	inline void AddValueToNumSet(VECTOR_ELEMENT_TYPE value, t_numbSetType type)
                                                    	{ GetNumSet(type)->AddElement(value); }
    inline CVector *GetNumSet(t_numbSetType type) const	{ return *(m_ppNumbSet + type); }
	inline size_t GetK() const							{ return GetNumSet(t_kSet)->GetAt(0); }
	inline size_t GetR() const							{ return colNumb() * GetK() / rowNumb(); }
    inline size_t lambda() const                        { return GetNumSet(t_lSet)->GetAt(0); }
	inline CVector **numbSet() const					{ return m_ppNumbSet; }
private:
	inline void setDataOwner(bool val)					{ m_bDataOwner = val; }
	inline bool isDataOwner() const  					{ return m_bDataOwner; }

	CVector **m_ppNumbSet;
	bool m_bDataOwner;
};

class C_BIBD : public C_InSys
{
 public:
 	C_BIBD(int v, int k, int lambda = 0);
	C_BIBD(const C_BIBD *pMaster, size_t nRow) : C_InSys(pMaster, nRow) {}
protected:
	void Init_BIBD_param(int v, int k, int lambda);
};

class C_tDesign : public C_BIBD
{
public:
	C_tDesign(int t, int v, int k, int lambda);
	C_tDesign(const C_tDesign *pMaster, size_t nRow);
	inline size_t getT() const							{ return m_t; }
	inline size_t lambda() const						{ return GetNumSet(t_lSet)->GetAt(getT() - 2); }
private:
	const size_t m_t;
};