#include <filesystem>
#include <set>
#include <fstream>
#include "TopGun.h"
#include "Table.h"
#include "SRGToolkit.h"

TopGunBase::TopGunBase(const kSysParam& param) : SizeParam(param), 
	m_param(param) {
	m_nRowsOut = param.val[t_nRowsInResultMatrix];
	reserveInputMatrixMemory(nRowsStart(), nMatricesReserved(), param.val[t_orderMatrices]);

	if (m_nRowsOut == 0)
		m_nRowsOut = m_numDays;
}

TopGunBase::~TopGunBase() {
	delete[] m_pInputMatrices;
	delete[] cnt();
	delete m_pMatrixInfo;
	delete[] m_pMatrixPerm;
	delete [] m_pGraphDB;
}

int TopGunBase::readMatrices(int tFolder, int nRows) {
	const std::string mf = *m_param.strVal[t_InputDataFileName];
	if (nRows == 0 && mf.length() && param(t_CBMP_Graph) == 2) {
		if ((m_nMatrices = readInputDataFile(mf, m_pInputMatrices, nMatricesReserved(), m_numDays, m_numPlayers)) < 1) {
			printfRed("*** Can't load input matrices data from file '%s'. Exit\n", mf.c_str());
			return m_nMatrices;
		}
		printf("%d '%s Matrices' (%d,%d,%d) created from data in file '%s'\n", m_nMatrices, matrixType(nRows), m_numPlayers, m_numDays, m_groupSize, mf.c_str());
	}
	else if ((m_nMatrices = loadMatrices(tFolder, nRows)) < 1) {
		printfRed("*** Can't load '%s Matrices'. Exit\n", matrixType(nRows));
		return m_nMatrices;
	}

	const auto totalLen = numMatrices2Process() * inputMatrixSize();
	if (!reallocStorageMemory(&m_pInputMatrices, totalLen, totalLen)) {
		printfRed("*** Memory allocation issue encountered while reading the matrices\n");
		return -1;
	}
	return m_nMatrices;
}
int TopGunBase::readInputDataFile(const std::string& fn, tchar *pData, int nMatricesMax, int nRows, int numColumnsInResult) {
#define LatinSquareData 1
	if (fn.length() == 0)
		return 0;
	std::ifstream mf;
	mf.open(fn, std::ios::in);
	if (!mf)
		return 0;
	std::string digits = "0123456789";
	tchar lineWithData = '*';
	const auto lenMatr = numColumnsInResult * nRows;
	int maxValue = 0, maxValuesOnLine = 0;
	int maxSize = lenMatr * nMatricesMax;
	std::string line;
	size_t start, end;
	bool bDataTypeDefined = false; // input data type must be defined in first input data line, for example: Format=LS
	int currentLineNumber = 0;
	int resultRow = 0;
	int numValuesInResultMatrix = 0, numValuesInResult = 0;
	int numMatricesInResult = 0;
	int inputDataType = 1; // input data type: 1 - Latin Square
	int inputLStype = param(t_inputLatinSquareType); // Latin Square data type can be (1-6) - all permutations of R,C,S
	while (!mf.eof()) {
		if (!getline(mf, line)) {
			mf.close();
			return numMatricesInResult;
		}
		currentLineNumber++;
		trim(line);
		if (line.empty() || line[0] == '#')
			continue;
		if (!bDataTypeDefined) {
			if (line.starts_with("Data=LS"))
				lineWithData = 'X';
			else
				lineWithData = '*';
			inputDataType = LatinSquareData;
			maxValue = numColumnsInResult / 2 - 1;
			maxValuesOnLine = maxValue + 1;
			bDataTypeDefined = true;
			continue;
		}
		if (line[0] != lineWithData && lineWithData != 'X')
			continue;

		start = lineWithData == 'X' ? 0 : 1;
		if (inputDataType == LatinSquareData && numValuesInResultMatrix == 0 && line.find("Latin Square") != std::string::npos) {
			// Temporary test for input LS data saved previously by this program
			// *1(3) : row - Hamiltonian RCS Latin Square(Atomic) :
			const char* lsName[] = { "RCS", "CRS", "SCR", "RSC", "CSR", "SRC" };
			for (int k = 0; k < 6; k++) {
				if (line.find(lsName[k]) != std::string::npos) {
					inputLStype = k + 1;
					break;
				}
			}
			continue;
		}
		auto nCol = line.length();
		for (int k = 0; k < maxValuesOnLine; k++) {
			start = line.find_first_of(digits, start);
			if (start == std::string::npos)
				end = std::string::npos;
			else {
				end = line.find_first_not_of(digits, start);
				if (end == std::string::npos)
					end = line.length();
			}
			if (end <= start) {
				printfRed("*** Error in data (# of values %d is less than expected %d): %s\n", k, maxValuesOnLine, line.c_str());
				printfRed("*** line %d, file: '%s'\n", currentLineNumber, fn.c_str());
				mf.close();
				return 0;
			}
			int iv = stoi(line.substr(start, end - start));
			start = end + 1;
			if (iv < 0 || iv > maxValue) {
				printfRed("*** Error in data (value #%d is above max=%d): %s\n", k + 1, maxValue, line.c_str());
				printfRed("*** line %d, file: '%s'\n", currentLineNumber, fn.c_str());
				mf.close();
				return 0;
			}
			if (inputDataType == LatinSquareData && numValuesInResultMatrix == 0) {
				numValuesInResult += lenMatr;
				if (numValuesInResult > maxSize) {
					printfYellow("*** Warning: input number of matrices exceeds limit(%d): line %d, file: '%s'\n", nMatricesMax, currentLineNumber, fn.c_str());
					mf.close();
					return numMatricesInResult;
				}
				memset(pData, 0, lenMatr);
			}
			if (inputDataType == LatinSquareData) {/**
				if (numValuesInResultMatrix == 0 && iv != 0) {
					printfRed("*** Error in data (first value of LS must be equal 0): %s\n", line.c_str());
					printfRed("*** line %d, file: '%s'\n", currentLineNumber, fn.c_str());
					mf.close();
					return 0;
				} **/
				// inputLStype (1-6):"RCS", "CRS", "SCR", "RSC", "CSR", "SRC" 
				int R = 0, C = 0, S = 0;
				int r = resultRow, c = k, s = iv;
				switch (inputLStype) {
				default:
				case 1: R = r; C = c; S = s; break;
				case 2: R = c; C = r; S = s; break;
				case 3: R = s; C = c; S = r; break;
				case 4: R = r; C = s; S = c; break;
				case 5: R = c; C = s; S = r; break;
				case 6: R = s; C = r; S = c; break;
				}
				if (R * numColumnsInResult + C * 2 + 1 >= lenMatr || S > maxValue) {
					printfRed("*** Error in data (value(s) out of range): %s\n", line.c_str());
					printfRed("*** line %d, file: '%s'\n", currentLineNumber, fn.c_str());
					mf.close();
					return 0;
				}
				*(pData + R * numColumnsInResult + C * 2) = C * 2;
				*(pData + R * numColumnsInResult + C * 2 + 1) = S * 2 + 1;
				numValuesInResultMatrix += 2;
			}
		}
		resultRow++;
		if (numValuesInResultMatrix >= lenMatr) {
			if (numMatricesInResult > 0) {
				if (MEMCMP(pData, pData + lenMatr, lenMatr) == 0) {
					printfYellow("* 2-partite Matrix #%d (created from input file with Latin Squares) is the same as previous\n", numMatricesInResult + 1);
					printfYellow("* line %d, file: '%s'\n", currentLineNumber, fn.c_str());
				}
			}
			pData += lenMatr;
			numMatricesInResult++;
			resultRow = 0;
			numValuesInResultMatrix = 0;
		}
	}
	mf.close();
	return numMatricesInResult;
}
int TopGunBase::loadMatrices(int tFolder, int nRows)
{
	// Matrix file name with folders: StartFolder/ColumnsxRowsxGroupSize[U1FName]/MatrixID.txt
	// StartFolder, Columns, Rows, GroupSize, U1FName - input parameters
	// MatrixID: starts from prefix (P, K, U, PN, KN, UN, PC, KC, UC, PCN, KCN, UCN), then 10 digits and extension ".txt"
	// Examples: 
	//   Logs/16x15x2/PM0123456789.txt 
	//   Logs/16x15x2_4444_88/U0123456789.txt
	//   Logs/27x13x2_999/U0123456789.txt
	//
	const std::string ch(getFileNameAttr(paramPtr()));
	const std::string ext(".txt");
	const auto extLen = ext.length();
	const auto chLen = ch.length();
	const auto nameExampleLength = strlen("0123456789.txt") + chLen;

	int nMatricesFromOneFile, nMatricesAll, nfr;
	nMatricesFromOneFile = nMatricesAll = nfr = 0;
	int nMax = nMatricesMax();
	int nReserved = nMatricesReserved();

	if (!nRows)
		nRows = nRowsStart();

	auto* pMatrixType = matrixType(nRows);

	std::string path_name;
	createFolderAndFileName(path_name, paramPtr(), tFolder, nRows);

	std::set<std::filesystem::path> sorted_by_name;

	for (auto& entry : std::filesystem::directory_iterator(path_name))
		sorted_by_name.insert(entry.path());

	for (auto& filename : sorted_by_name)
	{
		const auto& sfn = filename.string();
		if (sfn.length() <= nameExampleLength)
			continue;

		const std::string name = sfn.substr(sfn.length() - nameExampleLength);
		if (name.substr(0, chLen) != ch
			|| name.substr(name.length() - extLen) != ext)
			continue;

		const auto fnumber = name.substr(chLen, name.length() - chLen - extLen);
		if (fnumber.find_first_not_of("0123456789") != -1)
			continue;

		nMatricesFromOneFile = readInputData(sfn, nRows, nMatricesAll, &m_pInputMatrices, nMax, nReserved, m_pMatrixInfo);
		if (!nMatricesFromOneFile)
		{
			printfRed("Can't load file with '%s Matrices': %s\n", pMatrixType, sfn.c_str());
			break;
		}

		nfr++;
		if (param(t_printMatrices) & 1)
			printf("\n%d %d-rows '%s Matrices' loaded from file %s", nMatricesFromOneFile, nRows, pMatrixType, sfn.c_str());
		nMatricesAll += nMatricesFromOneFile;
		nMax -= nMatricesFromOneFile;
		if (nMax <= 0)
			break;
	}

	if (nMatricesAll)
		printf("\n%d %d-rows '%s Matrices' loaded from %d file(s). Process will start from FirstIndexOfStartMatrices(%d)\n",
			nMatricesAll, nRows, pMatrixType, nfr, param(t_nFirstIndexOfStartMatrices));
	else
		printfRed("*** Can't load '%s-matrices' from folder %s\n", ch.c_str(), path_name.c_str());

	if (m_pMatrixInfo && updateMatrReserved())
		m_pMatrixInfo->updateReservedMatrNumb(nReserved);

	return nMatricesAll;
}

#include <io.h>
void truncate_at_current_pos(FILE* f)
{
	fflush(f);                           // flush stdio buffer

	long pos = ftell(f);                 // current offset
	int fd = _fileno(f);                 // C runtime FD
	HANDLE h = (HANDLE)_get_osfhandle(fd);

	if (h == INVALID_HANDLE_VALUE)
		return;                          // handle error if needed

	SetFilePointer(h, pos, NULL, FILE_BEGIN);
	SetEndOfFile(h);                     // truncate
}

void TopGunBase::outputIntegratedResults(const paramDescr* pParSet, int numParamSet) {
	const char* pResFileName = paramPtr()->strVal[t_ResultsName]->c_str();
	const auto exploreMatrices = param(t_exploreMatrices);
	if (exploreMatrices < -3)
		return;

	std::string IntegratedResults;
	auto nRows = nRowsOut();
	if (!nRows)
		nRows = m_numDays;

	const auto finalReport = pParSet != NULL;
	createFolderAndFileName(IntegratedResults, paramPtr(), t_ResultFolder, nRows, pResFileName);
	FOPEN_F(f, IntegratedResults.c_str(), finalReport || !numParamSet? "w" : "r+");
	
	char buffer[128];
	if (finalReport) {
		const auto nParts = param(t_CBMP_Graph);
		if (nParts <= 1) 
			sprintf_s(buffer, "Factorization of K(%d)", param(t_numPlayers));
		else
			sprintf_s(buffer, "Factorization of K(%dx%d)", nParts, param(t_numPlayers)/nParts);

		setTableTitle(buffer);
		const auto totalMatr = reportResult(f);
		if (exploreMatrices && totalMatr && nRowsOut() == m_numDays) {
			reserveInputMatrixMemory(nRows, (int)totalMatr, 2);
			if (readMatrices(t_ResultFolder, nRows) >= 0) {
				orderAndExploreMatrices(nRowsOut(), 2, exploreMatrices);

				for (int i = 0; i < 2; i++)
					(m_pGraphDB + i)->reportResult(f, false);
			}
		}
		if (!m_reportInfo.empty())
			fprintf(f, m_reportInfo.c_str());
	}
	else {
		if (numParamSet) {
			char *pStr = NULL;
			while (fgets(buffer, sizeof(buffer), f) && !(pStr = strstr(buffer, DATE_TIME_TAG))); 
			
			if (pStr) {
				const auto offset = (long)strlen(pStr) + 3;
				fseek(f, -offset, SEEK_CUR);
			}
			numParamSet = 0;
		}
		else {
			numParamSet = 3;
			pParSet = paramPtr()->pParamDescr;
		}
	}

	const auto allowUndefinedCycles = param(t_allowUndefinedCycles);
	for (int j = 0; j < numParamSet; j++) {
		auto* paramNames = pParSet[j].paramNames;
		const int iMax = pParSet[j].numParams;
		switch (j) {
		case 0:
			fprintf(f, "\n%s\n", SECTION_PARAM_MAIN);
			for (int i = 0; i < iMax; i++)
				fprintf(f, "%30s: %d\n", paramNames[i], param(static_cast<paramID>(i)));
			break;

		case 1:
			fprintf(f, "\n%s\n", SECTION_PARAM_STRINGS);
			for (int i = 0; i < iMax; i++) {
				const auto ptr = paramPtr()->strVal[i];
				if (ptr)
					fprintf(f, "%30s: %s\n", paramNames[i], ptr->c_str());
			}
			break;
		default: {
				auto pntr = paramPtr()->u1fCycles[0];
				if (!pntr || !allowUndefinedCycles && pntr[0] == 1 && pntr[1] == m_numPlayers)
					break;

				fprintf(f, "\n%s\n", SECTION_PARAM_U1F_CONF);
				char buffer[256], *pBuf = buffer;
				const auto lenBuf = countof(buffer);
				SPRINTFS(pBuf, buffer, lenBuf, "%c", '{');
				const auto ngrp = pntr[0];
				pntr++;
				tchar symb;
				for (int i = 0; i < ngrp; i++) {
					SPRINTFS(pBuf, buffer, lenBuf, i ? ",{" : "{");
					int k = 0;
					while (symb = pntr[k])
						SPRINTFS(pBuf, buffer, lenBuf, k++ ? ",%d" : "%d", symb);

					SPRINTFS(pBuf, buffer, lenBuf, "}");
					pntr += MAX_CYCLES_PER_SET;
				}
				SPRINTFS(pBuf, buffer, lenBuf, "}");
				fprintf(f, "%30s: %s\n", paramNames[0], buffer);
			}
		}
	}
	fprintf(f, "\n");

	if (!finalReport) {
		reportResult(f);
		if (!m_reportInfo.empty())
			fprintf(f, m_reportInfo.c_str());
	}

	truncate_at_current_pos(f);
	FCLOSE_F(f);
}

int matrixSize;
ctchar* pStartMatrix;

int compare_matr_fn(const void* pA, const void* pB) {
	return memcmp(pA, pB, matrixSize);
}

int compare_matr_fn_perm(const void* pA, const void* pB) {
	const auto pA_ = pStartMatrix + matrixSize * *(uint*)pA;
	const auto pB_ = pStartMatrix + matrixSize * *(uint*)pB;
	return memcmp(pA_, pB_, matrixSize);
}

int TopGunBase::orderMatrices(int orderMatrixMode) {
	int nDuplicate = 0;
	const auto nMatrices = numMatrices2Process();
	matrixSize = inputMatrixSize();
	if (orderMatrixMode == 2) {
		delete[] m_pMatrixPerm;
		m_pMatrixPerm = new uint[nMatrices];
		for (uint i = 0; i < nMatrices; i++)
			m_pMatrixPerm[i] = i;

		pStartMatrix = inputMatrices();
		std::qsort(m_pMatrixPerm, nMatrices, sizeof(m_pMatrixPerm[0]), compare_matr_fn_perm);
	}
	else {
		auto* pMatrixDst = (tchar*)inputMatrices();
		std::qsort(pMatrixDst, nMatrices, inputMatrixSize(), compare_matr_fn);
		auto* pMatrixSrc = pMatrixDst + inputMatrixSize();
		for (uint i = 1; i < nMatrices; i++, pMatrixSrc += inputMatrixSize()) {
			if (memcmp(pMatrixDst, pMatrixSrc, inputMatrixSize())) {
				pMatrixDst += inputMatrixSize();
				if (nDuplicate)
					memcpy(pMatrixDst, pMatrixSrc, inputMatrixSize());
			}
			else
				nDuplicate++;
		}
	}
	return nDuplicate;
}

void TopGunBase::orderAndExploreMatrices(int nRows, int orderMatrixMode, int exploreMatrices) {
	const auto nDuplicate = orderMatrices(orderMatrixMode);
	printfGreen("%d '%s Matrices' sorted, %d duplicate matrices removed\n", numMatrices2Process(), matrixType(nRows), nDuplicate);
	m_nMatrices -= nDuplicate;
	if (orderMatrixMode != 2)
		return;

	TableAut Result(MATR_ATTR, m_numDays, m_numPlayers, 0, m_groupSize, true, true);
	Result.allocateBuffer(32);
	std::string ResultFile;
	createFolderAndFileName(ResultFile, paramPtr(), t_ResultFolder, nRows, "_OrderedMatrices.txt");

	SRGToolkit* pSRGtoolkit = NULL;
	if (exploreMatrices) {
#if OUT_SRG_TO_SEPARATE_FILE
		std::string ResultSRGFile;
		createFolderAndFileName(ResultSRGFile, paramPtr(), t_ResultFolder, nRows, "_SRG_Type_");
		const auto& srgResFile = ResultSRGFile;
#else
		const auto& srgResFile = ResultFile;
#endif
		pSRGtoolkit = new SRGToolkit(paramPtr(), nRows, srgResFile, exploreMatrices);
		m_pGraphDB = new GraphDB[2]();
		for (int i = 0; i < 2; i++)
			m_pGraphDB[i].setGraphType(i + 1);
	}

	Result.setOutFileName(ResultFile.c_str());
	printfGreen("Saved to a file: \"%s\"\n", ResultFile.c_str());
	for (uint i = 0; i < numMatrices2Process(); i++) {
		const auto idx = m_pMatrixPerm[i];
		const auto groupOrder = (*m_pMatrixInfo->groupOrdersPntr())[idx];
		Result.setGroupOrder(groupOrder);
		Result.setInfo(m_pMatrixInfo->cycleInfo(idx));
		const auto pMatr = inputMatrices() + idx * inputMatrixSize();
		Result.printTable(pMatr, true, false, nRows);
		if (groupOrder > 1)
			Result.printTableInfo(m_pMatrixInfo->groupInfo(idx));

		if (pSRGtoolkit) {
			if (!pSRGtoolkit->exploreMatrix(pMatr, m_pGraphDB, i+1, groupOrder)) {
				delete pSRGtoolkit;
				pSRGtoolkit = NULL;
			}
		}
	}

	if (pSRGtoolkit) {
		pSRGtoolkit->printStat();
		delete pSRGtoolkit;
	}
}

void TopGunBase::allocateMatrixInfoMemory(size_t nMatr, int orderMatrixMode) {
	// orderMatrixMode: 0 - No matrix reordering will be performed.
	//                  1 - Matrix reordering will be performed, but |Aut(M)| will not be needed.
	//                  2 - Matrix reordering will be performed AND |Aut(M)| will be used. 
	delete m_pMatrixInfo;
	m_pMatrixInfo = NULL;
	if (orderMatrixMode == 2) {
		m_pMatrixInfo = new CMatrixInfo((uint)nMatr);
		updateMatrReserved(false);
	}
}
