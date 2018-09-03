#pragma once
#include "BIBD_Enumerator.h"

template<class T>
class CPBIBD_Enumerator : public CBIBD_Enumerator<T>
{
public:
	CK CPBIBD_Enumerator(const C_InSys<T> *pBIBD, bool matrOwner = false, bool noReplicatedBlocks = false, int treadIdx = -1, uint nCanonChecker = 0) :
		CBIBD_Enumerator<T>(pBIBD, matrOwner, noReplicatedBlocks, treadIdx, nCanonChecker) {}
	CK virtual bool isPBIB_enumerator() const					{ return true; }
protected:
	CK virtual bool checkLambda(size_t lambdaCur) const;
	CK virtual void ReportLamdaProblem(T i, T j, size_t lambda) const;
	CK const char *getObjName() const							{ return "PBIBD"; }
	CK virtual int addLambdaInfo(char *buffer, size_t lenBuffer) const;
};

template<class T>
bool CPBIBD_Enumerator<T>::checkLambda(size_t lambdaCur) const {
	const auto lambdaSet = this->getInSys()->GetNumSet(t_lSet);
	for (size_t i = 0; i < lambdaSet->GetSize(); i++) {
		if (lambdaSet->GetAt(i) == lambdaCur)
			return true;
	}

	return false;
}

template<class T>
int CPBIBD_Enumerator<T>::addLambdaInfo(char *buf, size_t lenBuffer) const {
	const auto lambdaSet = this->getInSys()->GetNumSet(t_lSet);
	auto pBuf = buf;
	for (size_t i = 0; i < lambdaSet->GetSize(); i++) {
		if (pBuf == buf)
			pBuf += SNPRINTF(buf, lenBuffer, "{%d", lambdaSet->GetAt(i));
		else
			pBuf += SNPRINTF(pBuf, lenBuffer - (pBuf - buf), ", %d", lambdaSet->GetAt(i));
	}
	pBuf += SNPRINTF(pBuf, lenBuffer - (pBuf - buf), "}");
	return (int)(pBuf - buf);
}

template<class T>
void CPBIBD_Enumerator<T>::ReportLamdaProblem(T i, T j, size_t lambda) const {
	char buf[128];
	addLambdaInfo(buf, sizeof(buf));
	OUT_STRING(buff, 256, "Wrong number of common units in the rows (" ME_FRMT ", " ME_FRMT "): %zu is not in %s\n",
		i, j, lambda, buf);
}

template<class T>
class CInconsistentGraph_Enumerator : public CPBIBD_Enumerator<T>
{
public:
	CK CInconsistentGraph_Enumerator(const C_InSys<T> *pBIBD, bool matrOwner = false, bool firstPath = false, int treadIdx = -1, uint nCanonChecker = 0) :
		CPBIBD_Enumerator<T>(pBIBD, matrOwner, false, treadIdx, nCanonChecker), m_firstPath(firstPath) {}
protected:
	CK virtual bool TestFeatures(CEnumInfo<T> *pEnumInfo, const CMatrixData<T> *pMatrix, int *pMatrFlags = NULL, CEnumerator<T> *pEnum = NULL) const;
	bool makeFileName(char *buffer, size_t lenBuffer, const char *ext) const;
	CK const char *getObjName() const						{ return "I-Graph"; }
	CK virtual bool fileExists(const char *path, bool file = true) const;
	CK virtual bool createNewFile(const char *fName) const;
	CK bool SeekLogFile() const override					{ return true; }

private:
	inline bool firstPath() const					        { return m_firstPath; }
	const bool m_firstPath;
};

#define PRINT_DEBUG  0
#if PRINT_DEBUG
static bool fff(int nCols, int rowNumb, uchar *pMatr, FILE *file) {
	for (int j = 0; j < nCols; ++j) {
		int sum = 0;
		for (int row = 0; row < rowNumb; ++row)
			sum += *(pMatr + row * nCols + j);

		if (sum != 6) {
			fprintf(file, "Problem with column # %d\n", j);
			return true;
		}
	}
	return false;
}
#endif

template<class T>
bool CInconsistentGraph_Enumerator<T>::fileExists(const char *path, bool file) const {
	const bool retVal = CEnumerator<T>::fileExists(path, file);
	if (!file || !retVal)
		return retVal;

	// The answer is fake. When following statement is false, the file exist,
	// but we don't need the caller know that, becase for Inconsistent graphs
	// all outputs for same order graphs will be in the same file
	return designParams()->logFile != string(path);
}

template<class T>
bool CInconsistentGraph_Enumerator<T>::createNewFile(const char *fName) const {
	if (!fName)
		return firstPath();

	if (designParams()->logFile == fName)
		return false;

	designParams()->logFile = fName;
	return true; 
}

template<class T>
bool CInconsistentGraph_Enumerator<T>::TestFeatures(CEnumInfo<T> *pEnumInfo, const CMatrixData<T> *pMatrix, int *pMatrFlags, CEnumerator<T> *pEnum) const
{
	if (!CPBIBD_Enumerator<T>::TestFeatures(pEnumInfo, pMatrix, pMatrFlags))
		return false;

	if (!this->groupIsTransitive())
		return false;

	*pMatrFlags |= t_trahsitiveGroup;

	const auto nCols = pMatrix->colNumb();
	const auto nRows = pMatrix->rowNumb();
	const auto mult = nCols / nRows;

	// Need to check that transposed matrix is not isomorphic to constructed one
	// CMatrixData<T>  transpMatr;
	CMatrix<T> transpMatr(nCols, nCols, 1);
	transpMatr.InitTransposed(pMatrix, mult);
	const auto rowNumb = transpMatr.rowNumb();
	CMatrixCol<T> matrCol(&transpMatr);

	// For now we will consider only binary incidence systems
	// This part of the program will be a bit more complicated for general case
	assert(matrCol.rankMatr() <= 2);

	matrCol.initiateColOrbits(rowNumb, this->IS_enumerator());
	auto pColIdxMem = new T[nCols];

	const auto colOrbLen = matrCol.colOrbitLen();
	bool flag = false;

	CCanonicityChecker canonChecker(rowNumb, nCols, 2);
#if PRINT_DEBUG
	FOPEN(file, "C:\\Users\\andreii\\Calc\\fff.txt", "w");
	transpMatr.printOut(file, rowNumb, 99, NULL);
	int nMatr = 1;
#endif
	auto pMatr = matrCol.matrix()->GetDataPntr();
	T *pTmp = NULL;		// Memory for reordering the rows and columns of the matrix
	                    // If needed, it will be allocated. 
	T i = 0;
	while (true) {
		for (auto j = nCols; j--;)
			pColIdxMem[j] = j;

		CColOrbit<T> *pColOrbitNext = matrCol.colOrbits()[i];
		while (i < rowNumb) {
			T *pBeg = matrCol.matrix()->GetRow(i);
			CColOrbit<T> *pColOrbit = pColOrbitNext;
			auto pColOrbNext = pColOrbitNext = matrCol.colOrbits()[++i];
			auto pColIdx = pColIdxMem;

			// Loop over the orbits of columns
			while (pColOrbit) {
				// Each column orbits could be splited into sub-orbits
				// For binary incidence systems the algorythm is simple (and it is implemented here)
				// For general case we need to re-order columns of current orbits according to their
				// incidence with current element (corresponding to the current row of the matrix)

				const auto len = pColOrbit->length();
				auto idx = 0;
				auto idxLast = len;
				while (true) {
					// Find first zero from left
					while (idx < idxLast && pBeg[pColIdx[idx]])
						idx++;

					if (idx == idxLast)
						break;

					// Find first non-zero form right
					while (idx < --idxLast && !pBeg[pColIdx[idxLast]]);

					if (idx == idxLast)
						break;

					const auto i1 = pColIdx[idx++];
					const auto i2 = pColIdx[idxLast] - i1;

					auto pTmp = pBeg + i1;
					*pTmp = 1;
					*(pTmp+i2) = 0;

					int j = i;
					while (++j <= rowNumb) {
						auto tmp = *(pTmp += nCols);
						*pTmp = *(pTmp + i2);
						*(pTmp + i2) = tmp;
					}
				}

				pColOrbit = pColOrbit->next();
				if (i < rowNumb) {
					CColOrbit<T> *pNext = pColOrbit ? (CColOrbit<T> *)((char *)pColOrbNext + colOrbLen * len) : NULL;

					// Save the column's orbit information
					if (idx == 0 || idx == len) // Orbit was not splitted
						pColOrbNext->Init(len, pNext);
					else {
						CColOrbit<T> *pNxt = (CColOrbit<T> *)((char *)pColOrbNext + colOrbLen * idx);
						pColOrbNext->Init(idx, pNxt);
						pNxt->Init(len - idx, pNext);
					}

					pColOrbNext = pNext;
					pColIdx += len;
				}
#if PRINT_DEBUG
				if (fff(nCols, rowNumb, pMatr, file))
					break;
#endif
			}
		}

#if PRINT_DEBUG
		transpMatr.printOut(file, rowNumb, -nMatr, NULL);
		if (fff(nCols, rowNumb, pMatr, file))
			break;
#endif
		if (canonChecker.TestCanonicity(rowNumb, &matrCol, t_saveRowPermutations))
			break;  // Matrix is canonized

		// Reorder the rows and columns of the matrix according to  
		// permutations which were found during canonicity testing

		// Define the row number, where non-canonicity was noticed
		i = 0;
		const auto pPermRow = canonChecker.permRow();
		while (pPermRow[i] == i)
			i++;

		if (!pTmp)
			pTmp = new T[nCols * rowNumb];

		// Copy last (rowNum - i) rows of matrix into temporary buffer
		memcpy(pTmp, matrCol.matrix()->GetRow(i), nCols * (rowNumb - i));

#if PRINT_DEBUG
		fprintf(file, "i = %d\n", i);

		for (int k = 0; k < rowNumb; k++)
			fprintf(file, "%2d ", k);

		fprintf(file, "\n");

		auto ppp = pPermRow;
		for (int j = 0; j < 2; ++j) {
			for (int k = 0; k < rowNumb; k++)
				fprintf(file, "%2d ", ppp[k]);

			fprintf(file, "\n");
			ppp = canonChecker.permCol();
		}
#endif
		// Reorder the row and the columns:
		const auto pPermCol = canonChecker.permCol();
		for (T row = i; row < rowNumb; ++row) {
			auto pMatrTo = pMatr + row * nCols;
			const auto pMatrFrom = pTmp + (pPermRow[row] - i) * nCols;
			for (T j = 0; j < nCols; ++j)
				pMatrTo[j] = pMatrFrom[pPermCol[j]];
		}

#if PRINT_DEBUG
		transpMatr.printOut(file, rowNumb, nMatr++, &canonChecker);
		if (fff(nCols, rowNumb, pMatr, file))
			break;
#endif
	}
#if PRINT_DEBUG
	if (file)
		fclose(file);
#endif
	delete[] pColIdxMem;
	delete[] pTmp;

	auto pntr = pMatrix->GetDataPntr();
	if (mult == 1)
		return memcmp(pntr, transpMatr.GetDataPntr(), transpMatr.lenData()) != 0;

	auto pntrTr = transpMatr.GetDataPntr() - nCols;
	for (int i = 0; i < nRows; i++, pntr += nCols) {
		for (int j = 0; j < mult; ++j) {
			if (memcmp(pntr, pntrTr += nCols, nCols * sizeof(T)))
				return true;
		}	
	}

	return false;
}

template<class T>
bool CInconsistentGraph_Enumerator<T>::makeFileName(char *buffer, size_t lenBuffer, const char *ext) const
{
	const auto dirLength = this->getDirectory(buffer, lenBuffer);
	const auto pParam = designParams();
	const auto nVertex = this->rowNumb() * pParam->r / pParam->k;
	SNPRINTF(buffer + dirLength, lenBuffer - dirLength, "%ss_of_order_%d%s", getObjName(), nVertex, ext ? ext : FILE_NAME(""));
	return true;
}