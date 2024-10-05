#pragma once
#include "ColOrbits.h"
#include "CanonicityChecker.h"
#if USE_THREADS
#include <mutex>
#endif


Class2Def(CMatrixData) {
public:
	CK CMatrixData()						{ setDataOwner(false); }
	CK virtual ~CMatrixData()				{ 
		releaseData();
		delete partsInfo();
	}
	CC void InitWithData(T nRows, T nCols = 0, S maxElement = 1) {
		Init(nRows, nCols, maxElement, (uchar *)this + sizeof(*this));
	}

	CC inline T colNumb() const				{ return m_nCols; }
	CC inline T rowNumb() const				{ return m_nRows; }
	CC inline S maxElement() const			{ return m_nMaxElement; }
	CC inline S *GetRow(T nRow) const		{ return GetDataPntr() + nRow * colNumb(); }
	CC inline void setDataOwner(bool val)	{ m_bDataOwner = val; }
	CK inline S *GetDataPntr() const		{ return m_pData; }
	CK inline auto lenData() const			{ return m_nLenData; }
	CK inline void setLenData()				{ m_nLenData = sizeof(S) * m_nRows * m_nCols; }
	CK inline void setColNumb(T val)		{ m_nCols = val; }
	CK inline void SetDataPntr(S* pntr)		{ m_pData = pntr; }
	CK void InitTransposed(const CMatrixData *pMatr, uint mult = 1) {
		Init(pMatr->colNumb(), pMatr->rowNumb() * mult);
		auto *pMatrData = GetDataPntr();
		auto *pMatrSrc = pMatr->GetDataPntr();
		for (T i = 0; i < rowNumb(); ++i) {
			if (mult == 1) {
				for (T j = 0; j < pMatr->rowNumb(); ++j)
					*(pMatrData + j) = *(pMatrSrc + j * colNumb());
			}
			else {
				for (T j = 0; j < pMatr->rowNumb(); ++j) {
					const auto val = *(pMatrSrc + j * colNumb());
					for (uint k = 0; k < mult; ++k)
						*(pMatrData + j * mult + k) = val;
				}
			}

			pMatrSrc++;
			pMatrData += colNumb();
		}
	}
	CK void InitWithPermutedRows(const CMatrixData *pMatr, T *pPermRows, T rowNumb) {
		Init(rowNumb, pMatr->colNumb());
		const auto len = colNumb() * sizeof(S);
		S *pMatrData = GetDataPntr() - len;
		S *pMatrSrc = pMatr->GetDataPntr();
		for (T i = 0; i < rowNumb; ++i)
			memcpy(pMatrData += len, pMatrSrc + pPermRows[i] * len, len);
	}

	CC void Init(T nRows, T nCols, S maxElement = 1, S *data = NULL) {
		releaseData();
		if (!nRows)
			return;

		if (!nCols)
			nCols = nRows;

		m_nRows = nRows;
		m_nCols = nCols;
		setMaxElement(maxElement);
		setDataOwner(!data);
		setLenData();
		SetDataPntr(data? data : m_nLenData? new S[m_nRows * m_nCols] : NULL);
	}

	CK inline void AssignData(S *data)			{ memcpy(GetDataPntr(), data, m_nLenData); }
	void printOut(FILE* pFile = NULL, T nRow = ELEMENT_MAX, ulonglong matrNumber = UINT64_MAX, 
				  const CanonicityCheckerPntr pCanonCheck = NULL, ulonglong number = 0, int canonFlag = 1) const;
	CC virtual T numParts() const				{ return 1; }
	CC inline auto *partsInfo() const			{ return m_pPartInfo;  }
	CC inline S* ResetRowPart(T nRow, T idx, uint clean_flags = t_MatrixFlags::t_default_flag) const {
		auto *pRow = GetRow(nRow);
		T len;
		if (!partsInfo() || clean_flags & t_MatrixFlags::t_resetEntireRow)
			len = m_nCols;
		else
			pRow += partsInfo()->GetPartInfo(idx, &len);

		return static_cast<S *>(memset(pRow, 0, len * sizeof(*pRow)));
	}

	CC inline auto stabLengthExt() const			{ return m_nStabExtern; }
	CC inline auto InitPartsInfo(size_t nParts)		{ return m_pPartInfo = new BlockGroupDescr<T>(nParts); }
	CC inline S* GetRow(T nRow, T idx, S* pLen = nullptr) const {
		return GetRow(nRow) + (idx? partsInfo()->GetPartInfo(idx, pLen) : 0);
	}
	CC bool isSimple(bool* pFlag = NULL) const;
	CC void setMaxElement(S val)					{ m_nMaxElement = val; }
private:
	CK inline bool dataOwner()	const				{ return m_bDataOwner; }
	CK inline void releaseData()					{ if (dataOwner()) delete[] GetDataPntr();
													  SetDataPntr(NULL);
	                                                }
	T m_nRows;
	T m_nCols;

	S m_nMaxElement;
	bool m_bDataOwner = true;
	size_t m_nLenData;
	S *m_pData = NULL;
	BlockGroupDescr<T> *m_pPartInfo = NULL;
	mutable S m_nStabExtern;
#if USE_THREADS
	static ulonglong m_matrixCounter;
	static ulonglong GetNextCounter();
public:
	static void ResetCounter()				{ m_matrixCounter = 0; }
	static auto getMatrixCounter()			{ return m_matrixCounter; }
#endif
};


Class2Def(CMatrix) : public Class2(CMatrixData)
{
 public:
 	CK CMatrix(T nRows = 0, T nCols = 0, S maxElement = 1)	{ this->Init(nRows, nCols, maxElement); }
    CK CMatrix(const CMatrix *pMatrix, const int *pPermRow, const int *pPermCol)
	{
		Init(pMatrix->rowNumb(), pMatrix->colNumb());
		for (auto i = this->rowNumb(); i--;) {
			auto *pRow = GetRow(i);
			const auto *pRowFrom = pMatrix->GetRow(*(pPermRow + i));
			for (auto j = this->colNumb(); j--;)
				*(pRow + j) = *(pRowFrom + *(pPermCol + j));
		}
	}
	CK virtual ~CMatrix()										{}
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
	CK C_InSys(T nRows, T nCols, int t = 2, S maxElem = 1) : Class2(CMatrix)(nRows, nCols, maxElem), m_t(t) {
		setDataOwner(true);
		setMaxBlockIntrsection(nRows);
		m_ppNumbSet = createParamStorage(t_kSet); // Create 3 sets of vector: Lambda, R, and K 
	}

	CK C_InSys(const C_InSys *pMaster, T nRows) : Class2(CMatrix)(pMaster->rowNumb(), pMaster->colNumb()), m_t(pMaster->GetT())  {
		DuplicateMasterMatrix(pMaster, nRows);
	}

	CK void DuplicateMasterMatrix(const C_InSys* pMaster, T nRows) {
		setDataOwner(nRows == 0);
		m_ppNumbSet = pMaster->numbSet();
		setMaxBlockIntrsection(nRows);

		// Copy first nRow rows of master's matrix
		const size_t len = this->colNumb() * sizeof(*this->GetRow(0));
		for (auto i = nRows; i--;)
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
	CK inline T GetK() const								{ return GetNumSet(t_kSet)->GetAt(0); }
	CK inline T GetR(T nSuplemRows) const					{ return this->colNumb() * GetK() / (this->rowNumb() - nSuplemRows); }
	CK inline T lambda() const								{ return GetNumSet(t_lSet)->GetAt(0); }
	CK inline VectorPntr *numbSet() const					{ return m_ppNumbSet; }
	CK inline void setObjectType(t_objectType type)			{ m_objectType = type; }
	CK inline auto objectType() const						{ return m_objectType; }
	virtual T rowNumbExt() const							{ return this->rowNumb(); }
	CK inline auto maxBlockIntrsection() const              { return m_maxBlockIntersection; }
	CK inline void setMaxBlockIntrsection(T value)          { m_maxBlockIntersection = value; }
	CK void prepareFirstMatrixRow(int nDay) {
		auto* pNextCol = GetDataPntr();
		const auto numGroups = colNumb() / nDay;
		for (; nDay--; pNextCol += numGroups)
			memset(pNextCol, nDay, numGroups);
	}
	CK void convertToBinaryMatrix(ctchar *pTripleCol, int k, int dayNumb = 0) const {
		const auto pDataSave = GetDataPntr();
		const auto b = colNumb();

		auto* pNextCol = GetDataPntr() + b;
		const auto jMax = rowNumb() - 1;
		memset(pNextCol, 0, jMax * b * sizeof(pNextCol[0]));
		while (dayNumb--) {												// Iterate through all days
			const auto pLastCol = pNextCol + jMax / k;
			for (; pNextCol < pLastCol; pTripleCol += k, pNextCol++)	// Iterate through day groups
				for (int n = 0; n < k; n++)								// Iterate through elements of the group
					*(pNextCol + pTripleCol[n] * b) = 1;
		}
	}
	CK void adjustData(int val) {
		setColNumb(colNumb() - val);
		SetDataPntr(GetDataPntr() + val);
		setLenData();
	}

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
	const T m_t;
	T m_maxBlockIntersection;
};

Class2Def(C_BIBD) : public Class2(C_InSys)
{
 public:
	CK C_BIBD(int v, int k, int t = 2, int lambda = 0) : Class2(C_InSys)(v, lambda * v * (v - 1) / (k * (k - 1)), t)
													{ InitParam(v, k, lambda); }
	CK C_BIBD(const C_BIBD *pMaster, T nRow) : Class2(C_InSys)(pMaster, nRow) {}
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
	CK C_PBIBD(int v, int k, int r, const std::vector<uint> &lambdaSet) : Class2(C_InSys)(v, v * r/k) {
		InitParam(v, k, r, lambdaSet);
	}
	CK C_PBIBD(const Class2(C_InSys) *pMaster, T nRow) : Class2(C_InSys)(pMaster, nRow) {}
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
	CK CSemiSymmetricGraph(const Class2(C_InSys) *pMaster, T nRow) : Class2(C_PBIBD)(pMaster, nRow) {}
	CK ~CSemiSymmetricGraph()		{}
};

Class2Def(C_tDesign) : public Class2(C_BIBD)
{
public:
	CK C_tDesign(int t, int v, int k, int lambda);
	CK C_tDesign(const C_tDesign *pMaster, T nRow) : m_t(pMaster->getT()), Class2(C_BIBD)(pMaster, nRow) {}
	CK ~C_tDesign()											{}
	CK inline auto getT() const								{ return m_t; }
	CK inline auto lambda() const							{ return this->GetNumSet(t_lSet)->GetAt(getT() - 2); }
private:
	const T m_t;
};

Class2Def(CCombinedBIBD) : public Class2(CombDesignBase)
{
public:
	CK CCombinedBIBD(int v, int k, bool kirkmanTriples, const std::vector<uint>& lambda);
	CK CCombinedBIBD(const CCombinedBIBD* pMaster, T nRow) : Class2(CombDesignBase)(pMaster, nRow), m_ppParamSet(pMaster->paramSets()) {
		auto *pPartsInformation = InitPartsInfo(pMaster->numParts());
		pPartsInformation->CopyPartInfo(pMaster->partsInfo());
	}
	CK ~CCombinedBIBD()										{ if (isDataOwner()) this->deleteParamStorage(paramSets(), t_rSet); }
	CK inline VectorPntr *paramSets() const					{ return m_ppParamSet; }
	CK inline VectorPntr paramSet(t_numbSetType idx) const	{ return paramSets()[idx]; }
	virtual T rowNumbExt() const							{ return this->rowNumb() - 1; }
	CC virtual T numParts() const							{ return static_cast<T>(paramSet(t_lSet)->GetSize()); }
protected:
private:
	VectorPntr *m_ppParamSet;
};
