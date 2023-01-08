#include <fstream>      // std::ifstream
#include "BIBD_Enumerator.h"

using namespace std;

template class CBIBD_Enumerator<TDATA_TYPES>;

FClass2(CBIBD_Enumerator, int)::unforcedElement(const CColOrbit<S> *pOrb, int nRow) const
{
	const size_t diffWeight = this->getInSys()->GetK() - pOrb->columnWeight();
	if (diffWeight)
		return diffWeight == this->rowNumb() - nRow ? 1 : -1;

	// All units are there
	return 0;
}

FClass2(CBIBD_Enumerator, bool)::isValidSolution(const VECTOR_ELEMENT_TYPE* pSol) const
{
	// Check if solution is valid (for elimination of invalid solutions)
	// As of today (08/31/2021) we do it only for regular BIBDs OR first part of Combined BIBDs
	if (currentNumPart())
		return true;

	auto currRowNumb = this->currentRowNumb();
	if (currRowNumb <= firtstNonfixedRowNumber())
		return true;

	// For canonical BIBD the number of blocks containing any three elements cannot be
	// bigger than the number of blocks containing first, second and third elements.
	// Let's check it
	const auto lambda = getLambda(paramSet(t_lSet));
	const auto x0_3 = this->getX0_3();
	if (lambda == x0_3)     // Intersection of first three rows is maximal
		return true;		// Nothing to test

	const auto k = this->getInSys()->GetK();
	auto rowNumb = this->rowNumb();
	auto limit = currRowNumb + k + 1;
	if (limit >= rowNumb && currRowNumb + 3 < rowNumb) {
		// Theorem: The number of columns of canonical matrix which are forcible constructed by units
		// cannot be bigger than number of blocks containing first, second and third element.
		// We should start to check this condition on the first row which tested solutions could create
		// first column forcible constructed by units: rowNumb >= this->rowNumb() - k - 1 AND
		// at least three rows need to be constructed: rowNumb < this->rowNumb() - 3
		const auto* pColOrbit = this->colOrbit(currRowNumb);
		const auto* pRowSolution = pSol;
		limit -= rowNumb;
		auto nForcible = forcibleLambda(currRowNumb, 0);
		while (pColOrbit) {
			if (pColOrbit->columnWeight() == limit) {
				// Define the number of new columns that will be enforceable completed by units
				const auto newEnforsed = pColOrbit->length() - *pRowSolution;
				if (newEnforsed && (nForcible += newEnforsed) > x0_3)
					return false;
			}

			pRowSolution++;
			pColOrbit = pColOrbit->next();
		}
	}

	CMatrixCanonChecker::MakeRow(currRowNumb, pSol, false);

	// Define intersection of current row with previous one:
	auto lastRowToCheck = lenStabilizer();
	rowNumb -= lenStabilizer();
	const auto pMatrix = this->matrix();
	const auto *pCurrRow = pMatrix->GetRow(currRowNumb--);
	const auto pPrevRow = pMatrix->GetRow(currRowNumb--);
	const auto partsInfo = pMatrix->partsInfo();
	const auto colNumb = partsInfo? partsInfo->colNumb() : this->colNumb();
	const auto r = partsInfo ? colNumb * k / rowNumb : this->getR();

	S columns[32], *pColumnIdx = columns;
	if (lambda > countof(columns))
		pColumnIdx = new S[lambda];

	// Define "lambda" blocks, which contain both current and previous elements.
	// When doing that, check the necessary and sufficient conditions
	// for the intersection of the first (second), previous and current elements
	S idx = 0;
	S j = 0;
	while (true) {
		if (pCurrRow[j] && pPrevRow[j]) {
			if (idx == x0_3) {					// (x0_3+1)-th common block found
				if (j < r)   					//      among the first r blocks of design
					return false;				// The tested solution can not be used in the canonical matrix

				if (j < 2 * r - lambda) {	//      among the blocks, which contain second, but not first element
					S i = -1;				// Check necessary and sufficient conditions for the
					while (++i < idx) {     // intersection of second, previous and current elements
						if (pColumnIdx[i] >= lambda)
							break;
					}

					if (i == idx || pColumnIdx[i] >= r) // All blocks are amongth first lambda block OR [r+1,...2*r-lambda]
						return false;       // The tested solution can not be used in the canonical matrix
				}
				else {
					// When we are here, the intersection of second, previous and current elements is OK
					if (++lastRowToCheck == currRowNumb) // adjust the limit of the loop below
						return true;					 // there are no untested elements
				}
			}

			pColumnIdx[idx++] = j;
			if (idx == lambda)
				break;
		}

		++j;
	}

	//     for remaining elements:
	do {
		pCurrRow = pMatrix->GetRow(currRowNumb);
		auto j = x0_3;
		for (auto i = lambda; i--;) {
			if (pCurrRow[pColumnIdx[i]]) {
				if (!j--)
					return false;
			}
		}
	} while (--currRowNumb > lastRowToCheck);

	if (pColumnIdx != columns)
		delete[] pColumnIdx;

	return true;
}

FClass2(CBIBD_Enumerator, void)::getEnumerationObjectKey(char *pInfo, int len) const {
	SNPRINTF(pInfo, len, "(%3" _FRMT", %2" _FRMT", %2" _FRMT")",
		this->rowNumb(), this->getInSys()->GetK(), this->getInSys()->lambda());
}

FClass2(CBIBD_Enumerator, void)::CreateForcedRows() {
	CEnumerator::CreateForcedRows();
	if (designParams()->find_all_2_decomp)
		initDesignDB(NULL);
}

FClass2(CBIBD_Enumerator, void)::CreateAuxiliaryStructures(EnumeratorPntr pMaster) {
	if (pMaster && designParams()->find_all_2_decomp) {
		setMaster(pMaster);
		initDesignDB(pMaster);
	}
}

FClass2(CBIBD_Enumerator, void)::initDesignDB(const EnumeratorPntr pMaster, size_t rowAdj) {
	const auto b = matrix()->colNumb();
	const auto v = matrix()->rowNumb() - rowAdj;
	bool flag = designParams()->thread_master_DB;
	flag = pMaster ? flag : !flag || !designParams()->threadNumb;
	setDesignDB(flag ? new CDesignDB((v - 2) * b + LEN_HEADER) : pMaster? pMaster->designDB() : NULL);
}

FClass2(CBIBD_Enumerator, void)::AddMatrixToDB(const CMatrixCanonChecker *pCanonChecker, int rowAdj) const {
	auto* pMatr = pCanonChecker->matrix();
#if USE_MUTEX
	if (sharedDB())
		m_mutexDB.lock();
#endif
	// No need to keep first two rows, they are the same for all master BIBDs
	const auto idx = designDB()->AddRecord(pMatr->GetRow(2), pCanonChecker->groupOrder());
	if (outputMaster()) {
		outBlockTitle();
		pMatr->printOut(this->outFile(), pMatr->rowNumb() - rowAdj, matrix()->getMatrixCounter() + 1, pCanonChecker, idx + 1);
	}
#if USE_MUTEX
	if (sharedDB())
		m_mutexDB.unlock();
#endif
}

#if !CONSTR_ON_GPU
FClass2(CBIBD_Enumerator, bool)::makeFileName(char* buffer, size_t lenBuffer, const char* ext) const
{
	const auto inSys = this->getInSys();
	auto len = this->getDirectory(buffer, lenBuffer);
	len += SNPRINTF(buffer + len, lenBuffer - len, ME_FRMT"_" ME_FRMT"_", inSys->rowNumbExt(), inSys->GetK());
	len += this->addLambdaInfo(buffer + len, lenBuffer - len, ME_FRMT);
	SNPRINTF(buffer + len, lenBuffer - len, "%s", ext ? ext : FILE_NAME(""));
	return true;
}

FClass2(CBIBD_Enumerator, bool)::makeJobTitle(const designParam *pParam, char *buffer, int lenBuffer, const char *comment) const
{
	size_t lambdaSetSize = 0;
	auto len = getJobTitleInfo(buffer, lenBuffer);
	len += addLambdaInfo(buffer + len, lenBuffer - len, "%2" _FRMT, &lambdaSetSize);

	if (pParam->lambdaSizeMax() > lambdaSetSize) {
		auto maxSize = pParam->lambdaSizeMax() - lambdaSetSize;
		auto pBuf = buffer + len;
		pBuf += SNPRINTF(pBuf, lenBuffer - (pBuf - buffer), ")");
		while (maxSize-- > 0)
			pBuf += SNPRINTF(pBuf, lenBuffer - (pBuf - buffer), "   ");
	}
	else
		SNPRINTF(buffer + len, lenBuffer - len, ")%s", comment);

	return true;
}

FClass2(CBIBD_Enumerator, int)::getJobTitleInfo(char *buffer, int lenBuffer) const
{
	const auto v = this->rowNumb();
	const auto b = this->matrix()->colNumb();
	const auto k = this->getInSys()->GetK();
	return SNPRINTF(buffer, lenBuffer, "%s(%3" _FRMT", %3" _FRMT", %2" _FRMT", %2" _FRMT", ", getObjName(), v, b, b * k / v, k);
}

FClass2(CBIBD_Enumerator, int)::addLambdaInform(const Class1(CVector) *lambdaSet, char* buf, size_t lenBuffer, size_t *pLambdaSetSize) const
{
	const auto lambdaNumb = lambdaSet->GetSize();
	if (pLambdaSetSize)
		*pLambdaSetSize = lambdaNumb;

	const auto* pFrmt = "{%2d";
	int len = 0;
	for (size_t i = 0; i < lambdaNumb; i++) {
		len += SNPRINTF(buf + len, lenBuffer - len, pFrmt, lambdaSet->GetAt(i));
		pFrmt = ",%2d";
	}

	return len + SNPRINTF(buf + len, lenBuffer - len, "}");
}

int validLastLine(const char *filename, const char *testLine)
{
	ifstream fin;
	fin.open(filename);
	if (!fin.is_open())
		return -1;

	fin.seekg(-1, ios_base::end);                // go to one spot before the EOF
	bool flag = false;
	while (true) {
		char ch;
		fin.get(ch);                            // Get current byte's data

		if ((int)fin.tellg() <= 1) {            // If the data was at or before the 0th byte
			fin.seekg(0);                       // The first line is the last line
			break;
		}

		if (ch == '\n') {						// If the data was a newline
			if (flag)
				break;
		}
		else
			flag = true;						// first non-end line symbol found

		fin.seekg(-3, ios_base::cur);			// Move to the front of that data, then to the front of the data before it
	}

	char lastLine[256];
	fin.getline(lastLine, countof(lastLine));   // Read the current line
	fin.close();
	return strstr(lastLine, testLine)? 1 : 0;
}

FClass2(CBIBD_Enumerator, bool)::outFileIsValid(const struct stat& info, const char *pFileName) const {
	if (!CEnumerator::outFileIsValid(info))
		return false;

	if (designParams()->objType != t_objectType::t_BIBD || !designParams()->find_all_2_decomp)
		return true;

	if (designParams()->logFile.empty())
		return validLastLine(pFileName, END_OF_FILE);

	return designParams()->logFile != pFileName;
}

FClass2(CBIBD_Enumerator, bool)::outNonCombinedDesigns(designParam* pParam, const CDesignDB* pDesignDB, const char* pOutputInfo) {
	if (!pOutputInfo && !pDesignDB->recNumb()) {
		auto* pDB = pParam->designDB(1);
		if (pDB) {
			delete pDesignDB;
		} else
			pParam->setDesignDB(pDesignDB, 1);

		return true;
	}

	pParam->objType = t_objectType::t_BIBD;
	setDesignParams(pParam);
	const auto flag = setOutputFile();
	pParam->objType = t_objectType::t_CombinedBIBD;
	if (!flag)
		return false;

	const auto v = matrix()->rowNumb();
	const auto b = matrix()->colNumb();
	const size_t r = getR();
	const size_t length[] = { b, r, lambda(), r - lambda() };
	auto* pntr = matrix()->GetRow(0);
	memset(pntr + getR(), 0, (2 * b - r) * sizeof(*pntr));
	for (int j = 1; j <= 3; j++) {
		for (auto i = length[j]; i--;)
			*(pntr + i) = 1;

		pntr += length[j - 1];
	}

	pntr = matrix()->GetRow(2);
	const auto len = pDesignDB->recordLength() - LEN_HEADER;
	const auto nBIBDs = pDesignDB->recNumb();
	fprintf(outFile(), "\n\n%zd BIBD%s which %s NOT combined for ",
		nBIBDs, (nBIBDs > 1 ? "s" : ""), (nBIBDs > 1 ? "are" : "is"));

	if (!pOutputInfo) {
		const auto& lambdaSet = pParam->InterStruct()->lambda();
		fprintf(outFile(), "lambdas = {%d,%d}:\n", lambdaSet[0], lambdaSet[1]);
	}
	else {
		fprintf(outFile(), "any pair of lambdas\n");
	}

	for (size_t i = 0; i < nBIBDs; i++) {
		const auto rec = (const masterInfo*)pDesignDB->getRecord(i);
		memcpy(pntr, (uchar*)rec + LEN_HEADER, len);
		setGroupOrder(rec->groupOrder);
		matrix()->printOut(outFile(), v, i+1, this, rec->designNumber());
	}

	if (!pOutputInfo) {
		auto* pDB = pParam->designDB(1);
		if (pDB) {
			if (pDB->recNumb()) {
				// Calculate intersections with previously obtained DB of BIBDs
				auto* pIntersectionDB = new CDesignDB(pDesignDB->recordLength());
				pIntersectionDB->combineDesignDBs(pDB, pDesignDB, false, true);
				pParam->setDesignDB(pIntersectionDB, 1);
				delete pDB;
			}
			delete pDesignDB;
		}
		else
			pParam->setDesignDB(pDesignDB, 1);

		fclose(outFile());
	}
	else {
		fprintf(outFile(), "\n\n%s\n", pOutputInfo);
		const bool resetMTlevel = pParam->mt_level == 0;
		if (resetMTlevel) {
			// The row number, on which the threads will be launched was not defined.
			// Let's do it here by other parameters
			pParam->mt_level = define_MT_level(pParam);
		}
		compareResults(pParam->enumInfo(), 0, NULL, END_OF_FILE);
		if (resetMTlevel)
			pParam->mt_level = 0;

	}


	return true;
}

#endif
