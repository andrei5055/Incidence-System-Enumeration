#include <iostream>
#include <algorithm>
#include "stdafx.h"
#include "matrix.h"


template class CMatrixData<TDATA_TYPES>;
template class C_tDesign<TDATA_TYPES>;
template class CCombinedBIBD<TDATA_TYPES>;

#if USE_THREADS_ENUM
ulonglong CMatrixData<TDATA_TYPES>::m_matrixCounter;

FClass2(CMatrixData, ulonglong)::GetNextCounter() {
	const auto matrixNumber = ++m_matrixCounter;
	return matrixNumber;
}
#endif

FClass2(CMatrixData, bool)::isSimple(bool *pFlag) const {
	const auto nCol = this->colNumb();
	const auto* pData = GetDataPntr();
	for (T j = 1; j < nCol; j++) {
		const auto* pColData = ++pData;
		T i = -1;
		while (++i < this->rowNumb() && *pColData == *(pColData - 1))
			pColData += nCol;

		if (i == this->rowNumb()) {
			if (pFlag)
				*pFlag = j == 1;  // When TRUE, the first 2 columns are identical

			return false;
		}
	}

	return true;
}

FClass2(CMatrixData, void)::printOut(FILE* pFile, T nRow, ulonglong matrNumber, const CanonicityCheckerPntr pCanonCheck, ulonglong number, int canon) const
{
	if (nRow == ELEMENT_MAX)
		nRow = this->rowNumb();
#if USE_THREADS_ENUM
	if (!matrNumber)
		matrNumber = GetNextCounter();
#endif
	const auto nCol = this->colNumb();
	char buffer[256], *pBuf = buffer;
	auto lenBuf = sizeof(buffer);
	if (nCol >= lenBuf - 4)
		pBuf = new char[lenBuf = nCol + 14];

	if (matrNumber > 0) {
		auto* pTmp = pBuf;
		auto len = SNPRINTF(pTmp, lenBuf, "\nMatrix # %3llu", matrNumber);
		if (number)
			len = SNPRINTF(pTmp +=len, lenBuf-=len, ".%llu", number);
#if TEST
		else
			startPrinting = matrNumber >= START_PRINTING_AFTER;
#endif

		if (!canon)
			len = SNPRINTF(pTmp += len, lenBuf -= len, "  *** Non-canonical ***");
		else
		if (canon < 0)
			len = SNPRINTF(pTmp += len, lenBuf -= len, "  *** Cannonisity was not tested ***");

		if (pCanonCheck) {
			len = SNPRINTF(pTmp += len, lenBuf -= len, "    |Aut(M)| = %6zu", pCanonCheck->groupOrder());
			size_t exrtraGroupOrder = pCanonCheck->extraGroupOrder() ? pCanonCheck->extraGroupOrder()->groupOrder() : 0;
			if (exrtraGroupOrder > 1)
				len = SNPRINTF(pTmp += len, lenBuf -= len, "*%zu", exrtraGroupOrder);
		}

		SNPRINTF(pTmp +=len, lenBuf-=len, "\n");
	}
#if PRINT_CURRENT_MATRIX
	else {
		static int cntr;
		char* pBufTmp = pBuf;
		pBufTmp += sprintf_s(pBufTmp, lenBuf, "\n");
		memset(pBufTmp, '=', nCol);
		pBufTmp += nCol;
		sprintf_s(pBufTmp, lenBuf - (pBufTmp - pBuf), " # %d\n", ++cntr);
	}
#endif

	outString(pBuf, pFile);

	// Let's make the symbol table
	char symbols[32], * pSymb = symbols;
	const auto nMax = (partsInfo()? partsInfo()->numParts()  : this->maxElement()) + 1;
	if (nMax > sizeof(symbols))
		pSymb = new char[nMax];

	int i = -1;
	int iMax = nMax <= 9 ? nMax : 9;
	while (++i < iMax)
		*(pSymb + i) = '0' + i;

	i--;
	while (++i < nMax)
		*(pSymb + i) = 'A' + i - 10;

	for (S i = 0; i < nRow; i++) {
		auto pRow = this->GetRow(i);
		for (auto j = nCol; j--;)
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

TDesign()::C_tDesign(int t, int v, int k, int lambda) : Class2(C_BIBD)(v, k, t), m_t(t)
{
	// Define all lambdas: 
	int i = t;
	int *pLambda = new int[t - 2];
	while (--i > 1)
		lambda = (pLambda[i - 2] = lambda) * (v - i) / (k - i);

	// Initiate BIBD's parameter
	this->InitParam(v, k, lambda);

	// Add remaining lambda's to Lambda set
	while (++i < t)
		this->AddValueToNumSet(pLambda[i - 2], t_lSet);

	// Initiate matrix 
	this->Init(v, lambda * v * (v - 1) / (k * (k - 1)));
	delete[] pLambda;
}

CombinedBIBD()::CCombinedBIBD(int v, int k, bool kirkmanTriples, const std::vector<uint>& lambdaInp) : Class2(CombDesignBase)(0, k) {
	int lambda = 0;
	const auto v1 = v - 1;
	const auto k1 = k - 1;
	m_ppParamSet = this->createParamStorage(t_rSet); // Create 2 sets of vector (for Lambda and R of the component of combined BIBD)

	if (!kirkmanTriples) {
		std::vector<uint> lambdaSet(lambdaInp);
		std::sort(lambdaSet.begin(), lambdaSet.end(), std::greater<int>());
		const auto nSubDesigns = lambdaSet.size();
		auto partsInfo = this->InitPartsInfo(nSubDesigns);
		T shift = 0;
		for (size_t i = 0; i < nSubDesigns; ++i) {
			const auto lambdaCurr = lambdaSet[i];
			const auto r = lambdaCurr * v1 / k1;
			const auto b = r * v / k;
			assert(r * k1 == lambdaCurr * v1);
			m_ppParamSet[t_lSet]->AddElement(lambdaCurr);
			m_ppParamSet[t_rSet]->AddElement(r);
			lambda += lambdaCurr;
			partsInfo->SetPartInfo(i, shift, b);
			shift += b;
		}
	}
	else {
		lambda = 1;
		const auto n = (v - 3) / 6;
		const auto nSubDesigns = 3 * n + 1;
		auto partsInfo = this->InitPartsInfo(nSubDesigns);
		const auto b = v / k;
		assert(b == (2 * n + 1));
		for (T i = 0; i < nSubDesigns; ++i) {
			m_ppParamSet[t_lSet]->AddElement(0);
			m_ppParamSet[t_rSet]->AddElement(1);
			partsInfo->SetPartInfo(i, i*b, b);
		}
	}

	Init(v + 1, lambda * v * v1 / (k * k1));
	InitParam(v, k, lambda);
}