#include <filesystem>
#include <set>
#include "TopGun.h"

TopGunBase::TopGunBase(const kSysParam& param) : SizeParam(param), 
	m_param(param) {
	m_nRowsOut = param.val[t_nRowsInResultMatrix];
	mStartMatrixSize = m_numPlayers * nRowsStart();
	startMatrix = (tchar*)malloc(nMatricesMax() * mStartMatrixSize);
	if (m_nRowsOut == 0)
		m_nRowsOut = m_numDays;

	if (m_nRowsOut < 2 || m_nRowsOut > m_numDays) {
		printfRed("*** NRowsInResultMatrix(%d) must be either 0 or within the range 2:%d\n",
			nRowsOut(), m_numDays);
		myExit(1);
	}
}

bool TopGunBase::readStartMatrices() {
	if (nRowsStart() < 2 || nRowsStart() > m_numDays)
	{
		printfRed("*** NRowsInStartMatrix(%d) with UseMultiThreading must be in range 2:%d\n",
			nRowsStart(), m_numDays);
		myExit(1);
	}
	if ((nMatrices = getStartMatrices()) < 1)
	{
		printfRed("*** Cant load 'Start Matrices'. Exit\n");
		return false;
	}

	startMatrix = (tchar*)realloc(startMatrix, nMatrices * mStartMatrixSize);
	return true;
}

int TopGunBase::getStartMatrices() const
{
	// Matrix file name with folders: StartFolder/ColumnsxRowsxGroupSize[U1FName]/MatrixID.txt
	// StartFolder, Columns, Rows, GroupSize, U1FName - input parameters
	// MatrixID: starts from prefix (P, K, U, PM, KM, UM), then 10 digits and extension ".txt"
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
	auto* pStartMat = startMatrix;
	int reserved = nMatricesMax();
	std::string path_name;
	createFolderAndFileName(path_name, paramPtr(), t_StartFolder, nRowsStart());

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

		nMatricesFromOneFile = readStartData(sfn, &pStartMat, nMax, reserved);
		if (!nMatricesFromOneFile)
		{
			printfRed("Can't load file with 'Start Matrices': %s\n", sfn.c_str());
			break;
		}

		nfr++;
		printf("\n%d %d-rows 'Start Matrices' loaded from file %s", nMatricesFromOneFile, nRowsStart(), sfn.c_str());
		nMatricesAll += nMatricesFromOneFile;
		pStartMat += nMatricesFromOneFile * mStartMatrixSize;
		nMax -= nMatricesFromOneFile;
		if (nMax <= 0)
			break;
	}

	if (nMatricesAll)
		printf("\n%d %d-rows 'Start Matrices' loaded from %d file(s)\n", nMatricesAll, nRowsStart(), nfr);
	else
		printfRed("*** Can't load '%s-matrices' from folder %s\n", ch.c_str(), path_name.c_str());
	return nMatricesAll;
}

void TopGunBase::outputIntegratedResults(const paramDescr* pParSet, int numParamSet, const char* pResFileName) const {
	std::string IntegratedResults;
	auto nRows = nRowsOut();
	if (!nRows)
		nRows = m_numDays;

	std::string fileName(pResFileName);
	createFolderAndFileName(IntegratedResults, paramPtr(), t_ResultFolder, nRows, &fileName);
	FOPEN_F(f, IntegratedResults.c_str(), "w");
	
	reportResult(f);

	if (!m_reportInfo.empty())
		fprintf(f, m_reportInfo.c_str());

	for (int j = 0; j < numParamSet; j++) {
		auto* paramNames = pParSet[j].paramNames;
		const int iMax = pParSet[j].numParams;
		switch (j) {
		case 0:
			fprintf(f, "\nMain parameters:\n");
			for (int i = 0; i < iMax; i++)
				fprintf(f, "%30s: %d\n", paramNames[i], param(static_cast<paramID>(i)));
			break;

		case 1:
			fprintf(f, "\nParameters of string type:\n");
			for (int i = 0; i < iMax; i++) {
				const auto ptr = paramPtr()->strVal[i];
				if (ptr)
					fprintf(f, "%30s: %s\n", paramNames[i], ptr->c_str());
			}
			break;
		default:
			if (paramPtr()->u1fCycles[0]) {
				fprintf(f, "\nU1F configurations:\n");
				char buffer[128], *pBuf = buffer;
				const auto lenBuf = countof(buffer);
				SPRINTFS(pBuf, buffer, lenBuf, "%c", '{');
				auto pntr = paramPtr()->u1fCycles[0];
				const auto ngrp = pntr[0];
				pntr++;
				tchar symb;
				for (int i = 0; i < ngrp; i++) {
					SPRINTFS(pBuf, buffer, lenBuf, i ? ", {" : "{");
					int k = 0;
					while (symb = pntr[k])
						SPRINTFS(pBuf, buffer, lenBuf, k++ ? ", %d" : "%d", symb);

					SPRINTFS(pBuf, buffer, lenBuf, "}");
					pntr += MAX_UNIFOM_CONF_LENGTH;
				}
				SPRINTFS(pBuf, buffer, lenBuf, "}");
				fprintf(f, "%30s: %s\n", paramNames[0], buffer);
			}
			break;
		}
	}

	FCLOSE_F(f);
}

int matrixSize;

int compare_matr_fn(const void* pA, const void* pB) {
	return memcmp(pA, pB, matrixSize);
}

void TopGunBase::orderMatrices() const {
	matrixSize = mStartMatrixSize;
	std::qsort(pntrStartMatrix(), nMatrices, mStartMatrixSize, compare_matr_fn);
}