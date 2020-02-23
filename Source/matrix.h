#pragma once
#include "ColOrbits.h"
#include "CanonicityChecker.h"


Class2Def(CMatrixData) {
public:
	CK CMatrixData()						{ setDataOwner(false); }
	CK virtual ~CMatrixData()				{ if (dataOwner()) delete [] m_pData; }
	CC void InitWithData(S nRows, S nCols = 0, T maxElement = 1) {
		Init(nRows, nCols, maxElement, (uchar *)this + sizeof(*this));
	}

	CC inline S colNumb() const				{ return m_nCols; }
	CC inline S rowNumb() const				{ return m_nRows; }
	CC inline T maxElement() const			{ return m_nMaxElement; }
	CC inline T *GetRow(S nRow) const		{ return m_pData + nRow * m_nCols; }
	CC inline void setDataOwner(bool val)	{ m_bDataOwner = val; }
	CK inline T *GetDataPntr() const		{ return m_pData; }
	CK inline size_t lenData() const		{ return m_nLenData; }
	CK void InitTransposed(const CMatrixData *pMatr, uint mult = 1) {
		Init(pMatr->colNumb(), pMatr->rowNumb() * mult);
		auto *pMatrData = GetDataPntr();
		auto *pMatrSrc = pMatr->GetDataPntr();
		for (S i = 0; i < rowNumb(); ++i) {
			if (mult == 1) {
				for (S j = 0; j < pMatr->rowNumb(); ++j)
					*(pMatrData + j) = *(pMatrSrc + j * colNumb());
			}
			else {
				for (S j = 0; j < pMatr->rowNumb(); ++j) {
					const auto val = *(pMatrSrc + j * colNumb());
					for (S k = 0; k < mult; ++k)
						*(pMatrData + j * mult + k) = val;
				}
			}

			pMatrSrc++;
			pMatrData += colNumb();
		}
	}
	CK void InitWithPermutedRows(const CMatrixData *pMatr, S *pPermRows, S rowNumb) {
		Init(rowNumb, pMatr->colNumb());
		const auto len = colNumb() * sizeof(T);
		T *pMatrData = GetDataPntr() - len;
		T *pMatrSrc = pMatr->GetDataPntr();
		for (S i = 0; i < rowNumb; ++i)
			memcpy(pMatrData += len, pMatrSrc + pPermRows[i] * len, len);
	}

	CC void Init(S nRows, S nCols, T maxElement = 1, T *data = NULL) {
		if (!nRows)
			return;

		if (!nCols)
			nCols = nRows;

		m_nRows = nRows;
		m_nCols = nCols;
		m_nMaxElement = maxElement;
		setDataOwner(!data);
		m_nLenData = m_nRows * m_nCols * sizeof(T);
		m_pData = data? data : m_nLenData? new T[m_nRows * m_nCols] : NULL;
	}

	CK void AssignData(T *data)				{ memcpy(GetDataPntr(), data, m_nLenData); }
private:
	CK inline bool dataOwner()	const		{ return m_bDataOwner; }
	S m_nRows;
	S m_nCols;
	T m_nMaxElement;
	bool m_bDataOwner;
	size_t m_nLenData;
	T *m_pData;
};

Class2Def(CMatrix) : public Class2(CMatrixData)
{
 public:
 	CK CMatrix(S nRows = 0, S nCols = 0, T maxElement = 1)	{ this->Init(nRows, nCols, maxElement); }
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
	void printOut(FILE *pFile = NULL, S nRow = ELEMENT_MAX, ulonglong matrNumber = UINT64_MAX, const CanonicityCheckerPntr pCanonCheck = NULL) const
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
			pCanonCheck->outputAutomorphismInfo(pFile, this);
	}

 protected:
};

typedef enum {
	t_lSet,  // Order is important and it used in used in CCombBIBD_Enumerator class
	t_rSet,
	t_kSet
} t_numbSetType;

Class2Def(C_InSys) : public Class2(CMatrix)
{
 public:
	CK C_InSys(int nRows, int nCols, int t) : Class2(CMatrix)(nRows, nCols), m_t(t) {
		setDataOwner(true);	 
		m_ppNumbSet = createParamStorage(t_kSet); // Create 3 sets of vector: Lambda, R, and K 
	}

	CK C_InSys(const C_InSys *pMaster, size_t nRow) : Class2(CMatrix)(pMaster->rowNumb(), pMaster->colNumb()), m_t(pMaster->GetT())  {
		setDataOwner(nRow == 0);
		m_ppNumbSet = pMaster->numbSet();

		// Copy first nRow rows of master's matrix
		const size_t len = this->colNumb() * sizeof(*this->GetRow(0));
		for (T i = 0; i < nRow; i++)
			memcpy(this->GetRow(i), pMaster->GetRow(i), len);
	}

	CK C_InSys(int t = 2) : Class2(CMatrix)(), m_t(t)		{ m_ppNumbSet = NULL; }
	CK ~C_InSys() {
		if (!m_ppNumbSet || !isDataOwner())
			return;

		this->deleteParamStorage(m_ppNumbSet, t_kSet);
	}

	CK inline void AddValueToNumSet(VECTOR_ELEMENT_TYPE value, t_numbSetType type)
                                                    		{ GetNumSet(type)->AddElement(value); }
	CK inline Class1(CVector) *GetNumSet(t_numbSetType t) const	{ return *(m_ppNumbSet + t); }
	CK inline auto GetT() const								{ return m_t; }
	CK inline S GetK() const								{ return GetNumSet(t_kSet)->GetAt(0); }
	CK inline S GetR() const								{ return this->colNumb() * GetK() / this->rowNumb(); }
	CK inline S lambda() const								{ return GetNumSet(t_lSet)->GetAt(0); }
	CK inline VectorPntr *numbSet() const					{ return m_ppNumbSet; }
	CK inline void setObjectType(t_objectType type)			{ m_objectType = type; }
	CK inline t_objectType objectType() const				{ return m_objectType; }
	virtual S rowNumbExt() const							{ return this->rowNumb(); }
protected:
	CK inline bool isDataOwner() const						{ return m_bDataOwner; }
	CK VectorPntr *createParamStorage(int n) const  {
		const auto ppNumbSet = new VectorPntr [n + 1];
		for (int i = 0; i <= n; i++)
			ppNumbSet[i] = new Class1(CVector)();
		return ppNumbSet;
	}
	CK void deleteParamStorage(VectorPntr *ppParam, int n) {
		for (int i = 0; i <= n; i++)
			delete ppParam[i];

		delete[] ppParam;
	}
private:
	CK inline void setDataOwner(bool val)					{ m_bDataOwner = val; }

	t_objectType m_objectType;
	VectorPntr *m_ppNumbSet;
	bool m_bDataOwner;
	const uchar m_t;
};

Class2Def(C_BIBD) : public Class2(C_InSys)
{
 public:
	CK C_BIBD(int v, int k, int t = 2, int lambda = 0) : Class2(C_InSys)(v, lambda * v * (v - 1) / (k * (k - 1)), t)
													{ InitParam(v, k, lambda); }
	CK C_BIBD(const C_BIBD *pMaster, size_t nRow) : Class2(C_InSys)(pMaster, nRow) {}
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

Class2Def(C_PBIBD) : public Class2(C_InSys)
{
public:
	CK C_PBIBD(int v, int k, int r, const std::vector<uint> &lambdaSet) : Class2(C_InSys)(v, v * r/k, 2) {
		InitParam(v, k, r, lambdaSet);
	}
	CK C_PBIBD(const Class2(C_InSys) *pMaster, size_t nRow) : Class2(C_InSys)(pMaster, nRow) {}
	CK ~C_PBIBD() {}
protected:
	CK void InitParam(int v, int k, int r, const std::vector<uint> &lambdaSet) {
		this->AddValueToNumSet(k, t_kSet);
		this->AddValueToNumSet(r, t_rSet);
		for (auto lambda : lambdaSet)
			this->AddValueToNumSet(lambda, t_lSet);
	}
};

Class2Def(CSemiSymmetricGraph) : public Class2(C_PBIBD)
{
public:
	CK CSemiSymmetricGraph(int v, int k, int r, const std::vector<uint> &lambdaSet) :
		Class2(C_PBIBD)(v, k, r, lambdaSet) {}
	CK CSemiSymmetricGraph(const Class2(C_InSys) *pMaster, size_t nRow) : Class2(C_PBIBD)(pMaster, nRow) {}
	CK ~CSemiSymmetricGraph()		{}
};

Class2Def(C_tDesign) : public Class2(C_BIBD)
{
public:
	CK C_tDesign(int t, int v, int k, int lambda);
	CK C_tDesign(const C_tDesign *pMaster, size_t nRow) : m_t(pMaster->getT()), Class2(C_BIBD)(pMaster, nRow) {}
	CK ~C_tDesign()											{}
	CK inline S getT() const								{ return m_t; }
	CK inline S lambda() const								{ return this->GetNumSet(t_lSet)->GetAt(getT() - 2); }
private:
	const S m_t;
};

Class2Def(CCombinedBIBD) : public Class2(C_BIBD)
{
public:
	CK CCombinedBIBD(int v, int k, const std::vector<uint>& lambda);
	CK CCombinedBIBD(const CCombinedBIBD* pMaster, size_t nRow) : Class2(C_BIBD)(pMaster, nRow), m_ppParamSet(pMaster->paramSets()) {}
	CK ~CCombinedBIBD()										{ if (isDataOwner()) this->deleteParamStorage(paramSets(), t_rSet); }
	CK inline VectorPntr *paramSets() const					{ return m_ppParamSet; }
	CK inline VectorPntr paramSet(int idx) const			{ return paramSets()[idx]; }
	virtual S rowNumbExt() const							{ return this->rowNumb() - 1; }
private:
	VectorPntr *m_ppParamSet;
};