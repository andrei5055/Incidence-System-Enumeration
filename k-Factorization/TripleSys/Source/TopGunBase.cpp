#include <filesystem>
#include <set>
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
	if ((m_nMatrices = loadMatrices(tFolder, nRows)) < 1)
	{
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
		printf("\n%d %d-rows '%s Matrices' loaded from %d file(s)\n", nMatricesAll, nRows, pMatrixType, nfr);
	else
		printfRed("*** Can't load '%s-matrices' from folder %s\n", ch.c_str(), path_name.c_str());

	if (m_pMatrixInfo && updateMatrReserved())
		m_pMatrixInfo->updateReservedMatrNumb(nReserved);

	return nMatricesAll;
}

void TopGunBase::outputIntegratedResults(const paramDescr* pParSet, int numParamSet, const char* pResFileName) {
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
				const auto offset = (long)strlen(pStr) + 2;
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

	if (!finalReport) {
		reportResult(f);
		if (!m_reportInfo.empty())
			fprintf(f, m_reportInfo.c_str());
	}

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
		pSRGtoolkit = new SRGToolkit(paramPtr(), numPlayers(), nRows, m_groupSize, srgResFile, exploreMatrices);
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
