#pragma once
#include "ColOrbits.h"
#include "CanonicityChecker.h"

template<class T>
class CMatrixData {
public:
	CK CMatrixData()						{ setDataOwner(false); }
	CK virtual ~CMatrixData()				{ if (dataOwner()) delete [] m_pData; }
	CC void InitWithData(T nRows, T nCols = 0, T maxElement = 1) {
		Init(nRows, nCols, maxElement, (uchar *)this + sizeof(*this));
	}

	CC inline T colNumb() const				{ return m_nCols; }
	CC inline T rowNumb() const				{ return m_nRows; }
	CC inline T maxElement() const			{ return m_nMaxElement; }
	CC inline T *GetRow(T nRow) const		{ return m_pData + nRow * m_nCols; }
	CC inline void setDataOwner(bool val)	{ m_bDataOwner = val; }
	CK inline T *GetDataPntr() const		{ return m_pData; }
	CK inline size_t lenData() const		{ return m_nLenData; }
	CK void InitTransposed(const CMatrixData<T> *pMatr, int mult = 1) {
		Init(pMatr->colNumb(), pMatr->rowNumb() * mult);
		T *pMatrData = GetDataPntr();
		T *pMatrSrc = pMatr->GetDataPntr();
		for (T i = 0; i < rowNumb(); ++i) {
			if (mult == 1) {
				for (T j = 0; j < pMatr->rowNumb(); ++j)
					*(pMatrData + j) = *(pMatrSrc + j * colNumb());
			}
			else {
				for (T j = 0; j < pMatr->rowNumb(); ++j) {
					const auto val = *(pMatrSrc + j * colNumb());
					for (T k = 0; k < mult; ++k)
						*(pMatrData + j * mult + k) = val;
				}
			}

			pMatrSrc++;
			pMatrData += colNumb();
		}
	}
	CK void InitWithPermutedRows(const CMatrixData<T> *pMatr, T *pPermRows) {
		Init(pMatr->rowNumb(), pMatr->colNumb());
		const auto len = colNumb() * sizeof(T);
		T *pMatrData = GetDataPntr() - len;
		T *pMatrSrc = pMatr->GetDataPntr();
		for (T i = 0; i < rowNumb(); ++i)
			memcpy(pMatrData += len, pMatrSrc + pPermRows[i] * len, len);
	}

protected:

	CC void Init(T nRows, T nCols, T maxElement = 1, T *data = NULL) {
		if (!nCols)
			nCols = nRows;

		m_nRows = nRows;
		m_nCols = nCols;
		m_nMaxElement = maxElement;
		setDataOwner(!data);
		m_nLenData = m_nRows * m_nCols * sizeof(T);
		m_pData = data? data : m_nLenData? new T[m_nRows * m_nCols] : NULL;
	}
private:
	CK inline bool dataOwner()	const		{ return m_bDataOwner; }
	T m_nRows;
	T m_nCols;
	T m_nMaxElement;
	bool m_bDataOwner;
	size_t m_nLenData;
	T *m_pData;
};

template<class T> 
class CMatrix : public CMatrixData<T>
{
 public:
 	CK CMatrix(T nRows = 0, T nCols = 0, T maxElement = 1)	{ this->Init(nRows, nCols, maxElement); }
    CK CMatrix(const CMatrix *pMatrix, const int *pPermRow, const int *pPermCol)
	{
		Init(pMatrix->rowNumb(), pMatrix->colNumb());
		for (auto i = this->rowNumb(); i--;) {
			T *pRow = GetRow(i);
			const T *pRowFrom = pMatrix->GetRow(*(pPermRow + i));
			for (auto j = this->colNumb(); j--;)
				*(pRow + j) = *(pRowFrom + *(pPermCol + j));
		}
	}
	CK virtual ~CMatrix()										{}
	void printOut(FILE *pFile = NULL, T nRow = MATRIX_ELEMENT_MAX, ulonglong matrNumber = UINT64_MAX, const CCanonicityChecker<T> *pCanonCheck = NULL) const
	{
		if (nRow == (T)-1)
			nRow = this->rowNumb();

		auto nCol = this->colNumb();
		char buffer[256], *pBuf = buffer;
		size_t lenBuf = sizeof(buffer);
		if (nCol >= lenBuf - 4)
			pBuf = new char[lenBuf = nCol + 14];

		if (matrNumber > 0) {
			if (pCanonCheck)
				SNPRINTF(pBuf, lenBuf, "\nMatrix # %3llu    |Aut(M)| = %6d:\n", matrNumber, pCanonCheck->groupOrder());
			else
				SNPRINTF(pBuf, lenBuf, "\nMatrix # %3llu\n", matrNumber);
		}
#if PRINT_CURRENT_MATRIX
		else {
			static int cntr;
			char *pBufTmp = pBuf;
			pBufTmp += sprintf_s(pBufTmp, lenBuf, "\n");
			memset(pBufTmp, '=', nCol);
			pBufTmp += nCol;
			sprintf_s(pBufTmp, lenBuf - (pBufTmp - pBuf), " # %d\n", ++cntr);
		}
#endif

		outString(pBuf, pFile);

		// Let's make the symbol table
		char symbols[32], *pSymb = symbols;
		int nMax = this->maxElement() + 1;
		if (nMax > sizeof(symbols))
			pSymb = new char[nMax];

		int i = -1;
		int iMax = nMax <= 9 ? nMax : 9;
		while (++i < iMax)
			*(pSymb + i) = '0' + i;

		i--;
		while (++i < nMax)
			*(pSymb + i) = 'A' + i - 10;

		for (T i = 0; i < nRow; i++) {
			auto pRow = this->GetRow(i);
			for (size_t j = 0; j < nCol; j++)
				*(pBuf + j) = pSymb[*(pRow + j)];

			SNPRINTF(pBuf + nCol, 2, "\n");
			if (pFile)
				fputs(pBuf, pFile);
			else
				std::cout << pBuf;
		}

		if (pBuf != buffer)
			delete[] pBuf;

		if (pSymb != symbols)
			delete[] pSymb;

		if (pCanonCheck && pCanonCheck->groupOrder() > 1)
			pCanonCheck->outputAutomorphismInfo(pFile);
	}

 protected:
};

typedef enum {
        t_rSet,
        t_kSet,
        t_lSet
} t_numbSetType;

template<class T> 
class C_InSys : public CMatrix<T>
{
 public:
	CK C_InSys(int nRows, int nCols, int t) : CMatrix<T>(nRows, nCols), m_t(t) {
		setDataOwner(true);
		// Create 3 sets of vector: R, K and Lambda
		m_ppNumbSet = new CVector<T> *[3];
		for (int i = 0; i < 3; i++)
			m_ppNumbSet[i] = new CVector<T>();
	}

	CK C_InSys(const C_InSys *pMaster, size_t nRow) : CMatrix<T>(pMaster->rowNumb(), pMaster->colNumb()), m_t(pMaster->GetT())  {
		setDataOwner(nRow == 0);
		m_ppNumbSet = pMaster->numbSet();

		// Copy first nRow rows of master's matrix
		const size_t len = this->colNumb() * sizeof(*this->GetRow(0));
		for (T i = 0; i < nRow; i++)
			memcpy(this->GetRow(i), pMaster->GetRow(i), len);
	}

	CK C_InSys(int t = 2) : CMatrix<T>(), m_t(t) { m_ppNumbSet = NULL; }
	CK ~C_InSys() {
		if (!m_ppNumbSet || !isDataOwner())
			return;

		for (int i = 0; i < 3; i++)
			delete m_ppNumbSet[i];

		delete[] m_ppNumbSet;
	}

	CK inline void AddValueToNumSet(VECTOR_ELEMENT_TYPE value, t_numbSetType type)
                                                    		{ GetNumSet(type)->AddElement(value); }
	CK inline CVector<T> *GetNumSet(t_numbSetType t) const	{ return *(m_ppNumbSet + t); }
	CK inline uchar GetT() const							{ return m_t; }
	CK inline T GetK() const								{ return GetNumSet(t_kSet)->GetAt(0); }
	CK inline T GetR() const								{ return this->colNumb() * GetK() / this->rowNumb(); }
	CK inline T lambda() const								{ return GetNumSet(t_lSet)->GetAt(0); }
	CK inline CVector<T> **numbSet() const					{ return m_ppNumbSet; }
	CK inline void setObjectType(t_objectType type)			{ m_objectType = type; }
	CK inline t_objectType objectType() const				{ return m_objectType; }
private:
	CK inline void setDataOwner(bool val)					{ m_bDataOwner = val; }
	CK inline bool isDataOwner() const  					{ return m_bDataOwner; }

	t_objectType m_objectType;
	CVector<T> **m_ppNumbSet;
	bool m_bDataOwner;
	const uchar m_t;
};

template<class T>
class C_BIBD : public C_InSys<T>
{
 public:
	CK C_BIBD(int v, int k, int t, int lambda = 0) : C_InSys<T>(v, lambda * v * (v - 1) / (k * (k - 1)), t) 
													{ InitParam(v, k, lambda); }
	CK C_BIBD(const C_BIBD *pMaster, size_t nRow) : C_InSys<T>(pMaster, nRow) {}
	CK ~C_BIBD() {}
protected:
	CK void InitParam(int v, int k, int lambda) {
										if (!lambda)
											return;

										this->AddValueToNumSet(k, t_kSet);
										this->AddValueToNumSet(lambda * (v - 1) / (k - 1), t_rSet);
										this->AddValueToNumSet(lambda, t_lSet);
									}
};

template<class T>
class C_PBIBD : public C_InSys<T>
{
public:
	CK C_PBIBD(int v, int k, int r, const std::vector<int> &lambdaSet) : C_InSys<T>(v, v * r/k, 2) {
		InitParam(v, k, r, lambdaSet);
	}
	CK C_PBIBD(const C_InSys<T> *pMaster, size_t nRow) : C_InSys<T>(pMaster, nRow) {}
	CK ~C_PBIBD() {}
protected:
	CK void InitParam(int v, int k, int r, const std::vector<int> &lambdaSet) {
		this->AddValueToNumSet(k, t_kSet);
		this->AddValueToNumSet(r, t_rSet);
		for (auto lambda : lambdaSet)
			this->AddValueToNumSet(lambda, t_lSet);
	}
};

template<class T>
class CInconsistentGraph : public C_PBIBD<T>
{
public:
	CK CInconsistentGraph(int v, int k, int r, const std::vector<int> &lambdaSet) :
		C_PBIBD<T>(v, k, r, lambdaSet) {}
	CK CInconsistentGraph(const C_InSys<T> *pMaster, size_t nRow) : C_PBIBD<T>(pMaster, nRow) {}
	CK ~CInconsistentGraph()		{}
};

template<class T>
class C_tDesign : public C_BIBD<T>
{
public:
	CK C_tDesign(int t, int v, int k, int lambda);
	CK C_tDesign(const C_tDesign *pMaster, size_t nRow) : m_t(pMaster->getT()), C_BIBD<T>(pMaster, nRow) {}
	CK ~C_tDesign()											{}
	CK inline T getT() const								{ return m_t; }
	CK inline T lambda() const								{ return this->GetNumSet(t_lSet)->GetAt(getT() - 2); }
private:
	const T m_t;
};