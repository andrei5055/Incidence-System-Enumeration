#include <iostream>
#include "stdafx.h"
#include "matrix.h"
#include "CanonicityChecker.h"

using namespace std; 


CMatrix::CMatrix(const CMatrix *pMatrix, const int *pPermRow, const int *pPermCol)
{
    Init(pMatrix->rowNumb(), pMatrix->colNumb());
    for (auto i = rowNumb(); i--;) {
        MATRIX_ELEMENT_TYPE *pRow = GetRow(i);
        const MATRIX_ELEMENT_TYPE *pRowFrom = pMatrix->GetRow(*(pPermRow + i));
        for (auto j = colNumb(); j--;)
            *(pRow+j) = *(pRowFrom + *(pPermCol + j));
    }
}

void CMatrix::Init(size_t nRows, size_t nCols)
{
	if (!nCols)
		nCols = nRows;
        
    m_nRows = nRows;
    m_nCols = nCols;
    SetSize(m_nRows * m_nCols);
}

void CMatrix::printOut(FILE *pFile, size_t nRow, ulonglong matrNumber, const CCanonicityChecker *pCanonCheck) const
{   
	if (nRow == -1)
		nRow = rowNumb();
	
	auto nCol = colNumb();
	char buffer[256], *pBuf = buffer;
    size_t lenBuf = sizeof(buffer);
	if (nCol >= lenBuf - 4)
		pBuf = new char [lenBuf = nCol + 14];

    if (matrNumber > 0) {
		if (pCanonCheck)
			sprintf_s(pBuf, lenBuf, "\nMatrix # %3llu    |Aut(M)| = %6d:\n", matrNumber, pCanonCheck->groupOrder());
		else
			sprintf_s(pBuf, lenBuf, "\nMatrix # %3llu\n", matrNumber);
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
	int nMax = maxElement() + 1;
	if (nMax > sizeof(symbols))
		pSymb = new char [nMax];
	
	int i = -1;
	int iMax = nMax <= 9? nMax : 9;
	while (++i < iMax) 
		*(pSymb + i) = '0' + i;

	i--;
	while (++i < nMax)
		*(pSymb + i) = 'A' + i - 10;

	for (size_t i = 0; i < nRow; i++) {  
		MATRIX_ELEMENT_TYPE *pRow = GetRow(i);
		for (size_t j = 0; j < nCol; j++)
			*(pBuf+j) = pSymb[*(pRow+j)];
		
		strcpy_s(pBuf + nCol, 2, "\n");
		if (pFile)
			fputs(pBuf, pFile);
		else 
			cout << pBuf;
	}	
		
	if (pBuf != buffer)
		delete [] pBuf;
	
	if (pSymb != symbols)
		delete [] pSymb;

	if (pCanonCheck && pCanonCheck->groupOrder() > 1)
		pCanonCheck->outputAutomorphismInfo(pFile);
}

C_InSys::C_InSys(int nRows, int nCols) : CMatrix(nRows, nCols)
{
	setDataOwner(true);
    m_ppNumbSet = new CVector *[3];
    for (int i = 0; i < 3; i++)
        m_ppNumbSet[i] = new CVector();
}

C_InSys::C_InSys(const C_InSys *pMaster, size_t nRow) : CMatrix(pMaster->rowNumb(), pMaster->colNumb())
{
	setDataOwner(nRow == 0);
	m_ppNumbSet = pMaster->numbSet();

	// Copy first nRow rows of master's matrix
	const size_t len = colNumb() * sizeof(*GetRow(0));
	for (size_t i = 0; i < nRow; i++)
		memcpy(GetRow(i), pMaster->GetRow(i), len);
}

C_InSys::~C_InSys()
{
	if (!isDataOwner())
		return;

    for (int i = 0; i < 3; i++)
        delete m_ppNumbSet[i];

	delete [] m_ppNumbSet;
}

C_BIBD::C_BIBD(int v, int k, int lambda) : C_InSys(v, lambda * v * (v-1) / (k * (k - 1)))
{
	Init_BIBD_param(v, k, lambda);
}

void C_BIBD::Init_BIBD_param(int v, int k, int lambda)
{
	if (!lambda)
		return;

	AddValueToNumSet(k, t_kSet);
	AddValueToNumSet(lambda * (v - 1) / (k - 1), t_rSet);
	AddValueToNumSet(lambda, t_lSet);
}

C_tDesign::C_tDesign(int t, int v, int k, int lambda) : m_t(t), C_BIBD(v, k)
{
	// Define all lambdas: 
	int i = t;
	int *pLambda = new int[t - 2];
	while (--i > 1) {
		pLambda[i - 2] = lambda;
		lambda = lambda * (v - i) / (k - i);
	}

	// Initiate BIBD's parameter
	Init_BIBD_param(v, k, lambda);

	// Add remaining lambda's to Lambda set
	while (++i < t)
		AddValueToNumSet(pLambda[i - 2], t_lSet);

	// Initiate matrix 
	Init(v, lambda * v * (v - 1) / (k * (k - 1)));
	delete[] pLambda;
}

C_tDesign::C_tDesign(const C_tDesign *pMaster, size_t nRow) : m_t(pMaster->getT()), C_BIBD(pMaster, nRow)
{

}