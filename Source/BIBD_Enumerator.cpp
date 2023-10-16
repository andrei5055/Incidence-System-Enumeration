#include <fstream>      // std::ifstream
#include "BIBD_Enumerator.h"

using namespace std;

template class CBIBD_Enumerator<TDATA_TYPES>;

FClass2(CBIBD_Enumerator)::CBIBD_Enumerator(const InSysPntr pBIBD, uint enumFlags, int treadIdx, uint nCanonChecker) :
	Class2(C_InSysEnumerator)(pBIBD, enumFlags, treadIdx, nCanonChecker) {
	setR(getInSys()->GetR(0));
	const auto k = pBIBD->GetK();
	C_InSysEnumerator::setFirstUnforcedRow(pBIBD->rowNumb() - k);
	setUseLambdaCond(k);

	const auto λ = lambda();
	if (enumFlags & (t_use_3_condition | t_symmetrical_t_cond)) {
		const auto t = (enumFlags & t_symmetrical_t_cond) ? λ + 1 : 3;
		const auto lambdaSet = matrix()->GetNumSet(t_lSet);
		m_pRowIntersection = new CIntersection<T, S>(t, this->rowNumb(), lambdaSet);
		// Let's define some numbers which will only be used for x0_3 = 1
		// number of triplets of elements belonging to 1 or 0 common blocks
		m_Num_3[1] = λ * (k - 2);
		m_Num_3[0] = pBIBD->rowNumb() - m_Num_3[1] - 2;
	}
	
	m_bUseFilterFor_3d_RowSolutions = pBIBD->colNumb() == 3 * getR() - 2 * λ;
}

FClass2(CBIBD_Enumerator, int)::unforcedElement(const CColOrbit<S> *pOrb, int nRow) const
{
	const auto diffWeight = this->getInSys()->GetK() - pOrb->columnWeight();
	if (diffWeight)
		return diffWeight == this->rowNumb() - nRow ? 1 : -1;

	// All units are there
	return 0;
}

FClass2(CBIBD_Enumerator, VariableMappingPntr)::prepareCheckSolutions(size_t nVar) {
	if (!rowIntersection())
		return NULL;

	const auto λ = lambda();
	const auto t = (enumFlags() & t_symmetrical_t_cond) ? λ + 1 : 3;
	return rowIntersection()->prepareRowIntersections(matrix(), this->currentRowNumb(), λ, t);
}

FClass2(CBIBD_Enumerator, bool)::isValidSolution(const VECTOR_ELEMENT_TYPE* pSol) const
{
	// Check if solution is valid (for elimination of invalid solutions)
	// As of today (08/31/2021) we do it only for regular BIBDs OR first part of Combined BIBDs
	if (currentNumPart())
		return true;

	auto currRowNumb = this->currentRowNumb();
	if (currRowNumb < firtstNonfixedRowNumber())
		return true;

	const auto λ = getLambda(paramSet(t_lSet));
	if (currRowNumb == firtstNonfixedRowNumber()) {
		// All solutions for 3-rd row coud be descrided as (x, λ-x, λ-x, r-2*λ + x).
		if (*pSol == λ && useFilterFor_3d_RowSolutions()) {
			// If b = 3 * r - 2 * λ, and the solution with first element equal to λ is used 
			// for row #3, then the element #3 will have (r - 2*λ + 2*x) common blocks with 
			// any elements constructed in accordance with the solution (x, λ-x, λ-x, r - 2*λ + x)
			// Therefore such solution is only could be combined with the solution with x defined 
			// by equation: r - 2*λ + 2*x = λ, which has no integer solutions, when r - λ is odd:
			return (getR() - λ) % 2 == 0;
			// Notes: works for (17, 34, 16, 8, 7)
			// 
			// The only 3rd row solution that can be used with solutions  (λ, 0, 0, r-λ) 
			// is the solution (x=λ-(r-λ)/2, (r-λ)/2, (r-λ)/2, (r-λ)/2) or (x, λ-x, λ-x, λ-x).
			// 
			// It is quite obvious that for the construction of (v-2) remaining rows, we can use only 
			// one "descendant" of the solution (λ, 0, 0, r-λ). Denote by n the number of "descendants"
			// of the solutions (x, λ-x, λ-x, λ-x).
			// 
			// Using obvious equations to construct the completion of all blocks that are equivalent 
			// with respect to the first two elements:
			//    [λ-(r-λ)/2] * n + λ = (k-2)*λ
			//     [(r-λ)/2] * n = (k - 1)*(r - λ)
			//     [(r-λ)/2] * n + (r - λ) = k * (r - λ)
			// 
			// and standard equations for BIBD parameters:
			//   v * r = b * k
			//   r * (k - 1) = λ * (v - 1)
			//
			// We see that the parameters of such a BIBD should only be
			//      v = n + 3
			//      b = 4 * λ + 12 * λ / n
			//      r = 2 * λ + 4 * λ / n
			//      k = 1 + n / 2
			//      λ
			//
			// Since for BIBD k>=3, n>=4 and it should be even. Let's represent n as 
			// n = 2*m + 4, m>=0. Now we have:
			//      v = 2*m + 7                 = 2 * k + 1
			//      b = 4 * λ + 6 * λ / (m + 2) = 2 * λ * (2 * k + 1) / (k - 1)
			//      r = 2 * λ + 2 * λ / (m + 2) = 2 * λ * k / (k - 1)
			//      k = m + 3,
			//      λ
			// 
		}
		return true;
	}
	
	// For canonical BIBD the number of blocks containing any three elements cannot be
	// bigger than the number of blocks containing first, second and third elements.
	// Let's check it
	const auto x0_3 = this->getX0_3();
	if (λ == x0_3)     // Intersection of first three rows is maximal
		return true;   // Nothing to test

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

	const auto lastRow = currRowNumb;
	CMatrixCanonChecker::MakeRow(lastRow, pSol, false);
	OUTPUT_MATRIX(matrix(), outFile(), currentRowNumb() + 1, enumInfo(), -1);

	// Define intersection of current row with previous one:
	auto lastRowToCheck = lenStabilizer();
	rowNumb -= lenStabilizer();
	const auto pMatrix = this->matrix();
	const auto *pCurrRow = pMatrix->GetRow(currRowNumb--);
	const auto pPrevRow = pMatrix->GetRow(currRowNumb--);
	const auto partsInfo = pMatrix->partsInfo();
	const auto colNumb = partsInfo? partsInfo->colNumb() : this->colNumb();
	const auto r = partsInfo ? colNumb * k / rowNumb : this->getR();

	T columns[256], *pColumnIdx = columns;
	if (λ > countof(columns))
		pColumnIdx = new T[λ];

	// Define λ blocks, which contain both current and previous elements.
	// When doing that, check the necessary and sufficient conditions
	// for the intersection of the first (second), previous and current elements
	bool check_remaining_element = false;
	bool retVal = true;
	T idx = 0;
	T j = 0;
	while (true) {
		if (pCurrRow[j] && pPrevRow[j]) {
			if (idx == x0_3) {			// (x0_3+1)-th common block found
				if (j < r) {  			//      among the first r blocks of design
					retVal = false;		// The tested solution can not be used in the canonical matrix
					break;
				}

				if (j < 2 * r - λ) {	//      among the blocks, which contain second, but not first element
					S i = -1;			// Check necessary and sufficient conditions for the
					while (++i < idx) { // intersection of second, previous and current elements
						if (pColumnIdx[i] >= λ)
							break;
					}

					if (i == idx || pColumnIdx[i] >= r) {// All blocks are amongth first λ block OR [r+1,...2*r-λ]
						retVal = false;       // The tested solution can not be used in the canonical matrix
						break;
					}
				}
				else {
					// When we are here, the intersection of second, previous and current elements is OK
					if (++lastRowToCheck == currRowNumb) // adjust the limit of the loop below
						break;					 // there are no untested elements
				}
			}

			pColumnIdx[idx++] = j;
			if (idx == λ) {
				check_remaining_element = true;
				break;
			}
		}

		++j;
	}

	if (retVal && check_remaining_element) {
		//     for remaining elements:
		do {
			pCurrRow = pMatrix->GetRow(currRowNumb);
			auto j = x0_3;
			for (auto i = λ; i--;) {
				if (pCurrRow[pColumnIdx[i]]) {
					if (!j--) {
						retVal = false;
						break;
					}
				}
			}
		} while (retVal && --currRowNumb > lastRowToCheck);
	}

	if (pColumnIdx != columns)
		delete[] pColumnIdx;

	if (!retVal)
		return false;

	if (enumFlags() & t_use_3_condition /*t_symmetrical_t_cond*/) {
		const auto* pCurrRow = this->matrix()->GetRow(lastRow);
		// Get indices of the intersection of the i previous rows (pIntersection) 
		// one of which is the curent last row of the matrix and the other (i-1) correspond to all 
		// CurrentRowNumb()-1, i-1) combinations of previous currentRowNumb()-1 rows 
		// ... and the pointer to the number of the intersections (pNumb) we need to check (*pNumb = 1, when t = 3)
		const T* pNumb;
		auto* pIntersection = rowIntersection()->intersectionParam(&pNumb, lastRow);
		
		T num3[2] = { 0, 0 };
		const T t = 3;// λ + 1;
		for (T i = 2; i < t; i++) {
			const auto λ_prev = λ;
			for (auto k = pNumb[i - 2]; k--; pIntersection += λ_prev) {
				T val = 0;
				for (auto j = λ_prev; j--;) {
					if (*(pCurrRow + *(pIntersection + j)))
						val++;
				}

				if (val > x0_3) {
					// Any matrix constructed using this
					// solution cannot be canonical.
					FOPEN(f, "aaa.txt", "a");
					fprintf(f, "Catched!");
					FCLOSE(f);
					return false;
				}

				if (x0_3 == 1 && ++num3[val] > m_Num_3[val])
					return false; // It really works: see (17, 68, 16, 4, 3)
			}
		}
	}

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
	flag = pMaster ? flag : !flag || !designParams()->create_commonData();
	const auto recordLength = (v - 2) * (designParams()->compressMatrices() ? (b + 7) >> 3 : b);
	setRecordLen(recordLength);
	setDesignDB(flag ? new CDesignDB(recordLength + LEN_HEADER) : pMaster? pMaster->designDB() : NULL);
}

FClass2(CBIBD_Enumerator, void)::AddMatrixToDB(const CMatrixCanonChecker *pCanonChecker, int rowAdj) const {
	auto* pMatr = pCanonChecker->matrix();
#if USE_MUTEX
	if (sharedDB())
		m_mutexDB.lock();
#endif
	// No need to keep first two rows, they are the same for all master BIBDs
	recPtr pRec = pMatr->GetRow(2);
	uchar* memoryAllocated = NULL;
	uchar buffer[512];
	if (designParams()->compressMatrices()) {
		uchar* pntrMatr = recordLen() > sizeof(buffer)
			? memoryAllocated = new uchar[recordLen()]
			: buffer;

		auto* pntrCur = pntrMatr;
		for (int i = 2; i < pMatr->rowNumb(); i++) {
			uchar shift = 0;
			uchar val = 0;
			for (int j = 0; j < pMatr->colNumb(); j++) {
				if (*pRec++)
					val |= 0x80 >> shift;
				if (++shift == 8) {
					*pntrCur++ = val;
					val = shift = 0;
				}
			}

			if (shift)
				*pntrCur++ = val;
		}
		pRec = pntrMatr;
	}

	const auto idx = designDB()->AddRecord(pRec, (DB_INFO_DATA_TYPE)pCanonChecker->groupOrder());
	if (outputMaster()) {
		outBlockTitle();
		pMatr->printOut(this->outFile(), pMatr->rowNumb() - rowAdj, matrix()->getMatrixCounter() + 1, pCanonChecker, idx + 1);
	}

	if (memoryAllocated)
		delete[] memoryAllocated;

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

FClass2(CBIBD_Enumerator, void)::makeJobTitle(const designParam *pParam, char *buffer, int lenBuffer, const char *comment) const
{
	size_t lambdaSetSize = 0;
	auto len = getJobTitleInfo(buffer, lenBuffer);
	const auto flag = getInSys()->objectType() != t_objectType::t_Kirkman_Triple;
	if (flag)
		len += addLambdaInfo(buffer + len, lenBuffer - len, "%2" _FRMT, &lambdaSetSize);

	if (flag && pParam->lambdaSizeMax() > lambdaSetSize) {
		auto maxSize = pParam->lambdaSizeMax() - lambdaSetSize;
		auto pBuf = buffer + len;
		pBuf += SNPRINTF(pBuf, lenBuffer - (pBuf - buffer), ")");
		while (maxSize-- > 0)
			pBuf += SNPRINTF(pBuf, lenBuffer - (pBuf - buffer), "   ");
	}
	else
		SNPRINTF(buffer + len, lenBuffer - len, ")%s", comment);
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

	if (designParams()->objType != t_objectType::t_BIBD || designParams()->find_all_2_decomp != 1)
		return true;

	if (designParams()->logFile.empty())
		return validLastLine(pFileName, END_OF_FILE);

	return designParams()->logFile != pFileName;
}

FClass2(CBIBD_Enumerator, bool)::outNonCombinedDesigns(designParam* pParam, const CDesignDB* pDesignDB, string& outputInfo, bool addInfo) {
	pParam->objType = t_objectType::t_BIBD;
	setDesignParams(pParam);
	auto flag = setOutputFile();
	pParam->objType = t_objectType::t_CombinedBIBD;
	if (!flag)
		return false;

	const auto nBIBDs = pDesignDB->recNumb();
	auto* pDB = pParam->designDB(1);
	const auto& lambdaSet = pParam->InterStruct()->lambda();
	// It's not our first apperance here OR wee will be here at least three times.
	flag = pDB || lambdaSet[0] + pParam->lambdaStep() < lambdaSet[1];
	if (flag)
		fprintf(outFile(), "%s%zd BIBD%s which %s NOT combined for ", pParam->printEmptyLines()? "\n\n" : "",
				nBIBDs, (nBIBDs != 1 ? "s" : ""), (nBIBDs != 1 ? "are" : "is"));

	pParam->setEmptyLines(nBIBDs != 0 || !flag);
	if (addInfo) {
		if (flag)
			fprintf(outFile(), "lambdas = {%d,%d}%s\n", lambdaSet[0], lambdaSet[1], nBIBDs? ":" : "");
	}
	else {
		fprintf(outFile(), "any pair of lambdas%s\n", nBIBDs ? ":" : "");
	}

	size_t nSimpleBIBDs = 0;
	if (nBIBDs) {
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

		bool sameFirst2Columns;
		pntr = matrix()->GetRow(2);
		const auto len = pDesignDB->recordLength() - LEN_HEADER;
		for (size_t i = 0; i < nBIBDs; i++) {
			const auto rec = (const masterInfo<DB_INFO_DATA_TYPE> *)pDesignDB->getRecord(i);
			memcpy(pntr, (uchar*)rec + LEN_HEADER, len);

			if (matrix()->isSimple(&sameFirst2Columns))
				nSimpleBIBDs++;
			else
			if (pParam->printOnlySimpleDesigns()) {
				if (sameFirst2Columns)
					break;
				continue;
			}

			if (flag) {
				setGroupOrder(rec->groupOrder);
				matrix()->printOut(outFile(), v, i + 1, this, rec->designNumber());
			}
		}
	}

	if (addInfo) {
		char str[256];
		sprintf_s(str, "   {%d, %d}: %zd (%zd)", lambdaSet[0], lambdaSet[1], nBIBDs, nSimpleBIBDs);
		outputInfo += str;

		if (pDB) {
			if (!nBIBDs) {
				delete pDesignDB;
				((CDesignDB*)pDB)->resetRecNumb();
			} else {
				if (pDB->recNumb()) {
					// Calculate intersections with previously obtained DB of BIBDs
					auto* pIntersectionDB = new CDesignDB(pDesignDB->recordLength());
					pIntersectionDB->combineDesignDBs(pDB, pDesignDB, false, true);
					pParam->setDesignDB(pIntersectionDB, 1);
					delete pDB;
				}
				delete pDesignDB;
			}
		} else
			pParam->setDesignDB(pDesignDB, 1);

		fclose(outFile());
	}
	else {
		if (pParam->printOnlySimpleDesigns())
			fprintf(outFile(), "\n\n" ONE_LINE_BLOCK "(NOTE: It was required to print only simple non-combined BIBDs) " ONE_LINE_BLOCK "\n");

		fprintf(outFile(), "\n\n%s\n", outputInfo.c_str());

		// Set MT_level used for regular BIBDs enumeration.
		const auto MT_level = pParam->MT_level();
		pParam->set_MT_level(pParam->MT_level(1), 0);
		compareResults(pParam->enumInfo(), 0, NULL, END_OF_FILE);
		// Restore MT_Level
		pParam->set_MT_level(MT_level);
	}

	return true;
}

#define USE_CHECK   0

#if USE_CHECK
void check(const size_t* solution, const size_t* right_part, size_t last_dx = 0) {
	static size_t r_part[2], last;
	if (!solution) {
		last = last_dx;
		for (int i = 0; i < countof(r_part); i++)
			r_part[i] = right_part[i];

		return;
	}

	size_t sum[2] = { solution[0] + solution[1], solution[1] };
	for (int i = 2; i <= last; i++) {
		sum[0] += solution[i];
		sum[1] += i * solution[i];
	}

	for (int i = 0; i < countof(r_part); i++) {
		if (sum[i] + right_part[i] != r_part[i])
			assert(false);
	}
}
#else
#define check(x, y, ...)
#endif
static bool is_valid_solution(const size_t* solution, size_t r2_minus_λ, int λ) {
	// Let's consider the descendants of the solution sol[i] != 0 which has i = {0, 1,... λ} common
	// blocks with the first two elements of the BIBD. Each of them must have j = {λ, λ-1,..., 0}
	// common blocks among both sets of blocks B1 and B2 containing r-λ blocks each.
	// Sometimes, this is impossible for small i, and therefore the solution must be discarded.
	// 
	// For j1, j2, the minimum intersection in B1 (and in B2) will be j1 + j2 - (r-λ).
	// Therefore, the solution should be discarded if 2 * (j1 + j2 - (r-λ)) > λ
	// Last inequality could be rewritten as 2 * (j1 + j2) > 2*r - λ.
	// When j1 = j2 we also should check that sol[i]>=2 for corresponding i = λ - j1. 
	int j = -1;
	while (!solution[++j]);

	int j1 = λ - j;
	// Loop while minimal intersection of two solutions 
	// will be greater than λ
	while (4 * j1 > r2_minus_λ) {
		const auto numDescendand = solution[λ - j1];
		if (numDescendand >= 2) {
			// At least two descendants of that solutions
			// will be in more than λ blocks
			return false;  // invalid solution;
		}

		if (numDescendand) {
			int j2 = j1;
			while (2 * (j1 + --j2) > r2_minus_λ) {
				if (solution[λ - j2])
					return false;
			}
		}

		j1--;
	}
	return true;
}

void solve_system(size_t v, size_t r, size_t k, size_t λ, size_t d) {
	auto len = d;
	auto i = d + 1;
	size_t sol[32], *solution = i < countof(sol)? sol : new size_t [i];
	memset(solution, 0, i * sizeof(solution[0]));
	size_t rightPart[2] = { v - 3, λ * (k - 2) - d };

	check(NULL, rightPart, d);
	const auto r2_minus_λ = 2 * r - λ;
	int numSolution = 0;
	int step = -1;
	size_t val;

#define PRINT_SOL  0
#if PRINT_SOL
	FOPEN(f, "sol.txt", "w");
	char buffer[512];
	char* pBuff = buffer;
	pBuff += sprintf_s(pBuff, countof(buffer), "idx:");
	for (int i = 0; i <= d; i++)
		pBuff += sprintf_s(pBuff, countof(buffer) - (pBuff - buffer), "%4d", i);

	fprintf(f, "%s\n", buffer);

#endif
	while (true) {
		while ((i += step) >= 0 && i <= len) {
			if (step > 0) {
				if (solution[i]) {
					--solution[i];
					rightPart[1] += i;
					rightPart[0]++;
					check(solution, rightPart);
					step = -1;
				}
				else {
					if (i >= len) {       // the current highest coordinate has been set to 0
						break;
						if (--len == 2)   // decrease the number of the last modified coordinate
							break;
						step = -1;
					}
				}

				continue;
			}

			switch (i) {
			case 1:     val = rightPart[1];
				if (rightPart[0] >= val) {
					rightPart[1] = 0;
					rightPart[0] -= solution[1] = val;
					check(solution, rightPart);
				}
				else {
					step = 1;
					if (val = solution[i = 2]) {
						rightPart[1] += val << 1;
						rightPart[0] += val;
						solution[i] = 0;
						check(solution, rightPart);
					}
				}
				break;
			case 0:     solution[0] = rightPart[0];
				rightPart[0] = 0;
				break;
			default:
				val = rightPart[1] / i;
				if (val > rightPart[0])
					val = rightPart[0];

				if (solution[i] = val) {
					rightPart[1] -= val * i;
					rightPart[0] -= val;
					check(solution, rightPart);
				}
			}
/*
			if (!rightPart[0]) {
				if (rightPart[1]) {
					// It is impossible to satisfy the second equation
					// We need to reject the coordinate just assigned
					val = solution[i];
					rightPart[1] += val * i;
					rightPart[0] += val;
					solution[i] = 0;
					check(solution, rightPart);

					while (i < len && !solution[++i]);
					if (i == len)
						break;
				}

				step = 1;
				if (!i--)
					rightPart[i = 0] += solution[0];
			}
			*/
		}


		while (len > 0 && !solution[len])
			len--;

		if (!rightPart[0] && !rightPart[1] && is_valid_solution(solution, r2_minus_λ, (int)λ)) {
			// Solution found
			numSolution++;
#if PRINT_SOL
			pBuff = buffer;
			pBuff += sprintf_s(pBuff, countof(buffer), "%2d: ", numSolution);
			for (int i = 0; i <= d; i++)
				pBuff += sprintf_s(pBuff, countof(buffer) - (pBuff - buffer), "%4zd", solution[i]);

			fprintf(f, "%s\n", buffer);
#endif
		}

		step = 1;
		rightPart[0] += solution[0] + solution[1];
		rightPart[1] += solution[1];
		solution[0] = solution[1] = 0;
		check(solution, rightPart);

		if (len <= 1)   // no reasons to change x{1}, it is the only one in the 2nd equation
			break;

		i = 1;
		while (!solution[++i]);
		i--;
	}

#if PRINT_SOL
	fclose(f);
#endif
	if (solution != sol)
		delete [] solution;
}

FClass2(CBIBD_Enumerator, bool)::useAsRightPart(CRowSolution<TDATA_TYPES>* pRowSol, PERMUT_ELEMENT_TYPE idx) {
	// DEfine by d - index of solutions used for 3-d row
	// System of equations for remaining v-2 rows:
	//         n0 +           n1 +           n2 + ... +           nd = v - 3
	//                        n1 +         2*n2 + ... +         d*nd = λ*(k - 2) - d
	//       λ*n0 +     (λ-1)*n1 +     (λ-2)*n2 + ... +     (λ-d)*nd = (r-λ)*(k-1) - λ + d
	// (r-2*λ)*n0 + (r-2*λ+1)*n1 + (r-2*λ+2)*n2 + ... + (r-2*λ+d)*nd = (b-2*r+λ)*k - (r-2*λ+d)
	// 
	// Actually, we should exclude third and fourth equations, because 
	//    [3] = λ * [1] - [2]  and [4] = (r-2*λ) * [1] + [2] 
	const auto d = pRowSol->solutionIndex();
	if (idx == d)
		return true;

	const auto k = this->getInSys()->GetK();
	const auto λ = this->getInSys()->lambda();
	const auto v = matrix()->rowNumb();
	const auto r = this->getR();
	const auto a = pRowSol->currSolution()[0];
	if (a * (v - 2) == λ * (k - 2))
		return false;

	if (!idx)
		solve_system(v, r, k, λ, d);

	return true; 
}
#endif
