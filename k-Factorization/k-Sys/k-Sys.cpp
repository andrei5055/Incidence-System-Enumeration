// k-Sys.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <fstream>
#include <algorithm>
#include <regex>
#include "TopGun.h"

const char* intParamNames[]{
	"nPlayers",
	"GroupSize",
	"UseUniform1Factorization",
	"Use2RowsCanonization",
	"SubmatrixGroupOrderMin",
	"ResultGroupOrderMin",
	"USE_GPU",
	"UseMultiThreading",
	"NThreads",
	"NRowsInStartMatrix",
	"MaxNumberOfStartMatrices",
	"FirstIndexOfStartMatrices",
	"NRowsInResultMatrix",
	"ExpectedResult",
	"UseCheckLinksV",
	"UseRowsPrecalculation",
	"UseSolutionCliquesAfterRow",
	"UseAutForPrecRows",	// Use the Automorphism group of precalculated rows.
	"LastRowSecondPlayer",	// Stop invoking addRow after processing all row solutions for this second player
	"PrintMatrices",
	"SavingMatricesToDisk",
	"MatrixCanonInterval",
	"CheckConstructedMatrices",
	"UseSS",
	"p1f_counter",
	"AutLevelMinDef",	// minimum and maximum numbers of matrix rows for which the 
	"AutLevelMaxDef",	// automorphism groups will be calculated to filter solutions
	"AutLevelMinApp",	// minimum and maximum numbers of matrix rows for which 
	"AutLevelMaxApp",	// solution filtering by automorphism groups will be applied
	"AutDirection",		// Direct (0) or reversed order.  0: G2, G3, ... Gt OR 1: Gt, ..., G3, G2
	"AutGroupsToTest",  // Number of automorphism groups that need to be tested.
	"AutSaveTestedTrs", // 0 - not saving; 1 - saving only while testing on groups ; 
	                    // 2 - saving only while testing paths; 3 - saving while testing on groups OR paths
	"UseImproveMatrix",
	"UseCombinedSolutions",
	"OutAutomorphismGroup",
	"NestedGroups",
	"GridSize",
	"BlockSize",
	"OrderMatrices",
	"AllowMissingCycles", // 1  - Not all cycles need to be represented.
	"Any2RowsConvertToFirst2",
};

const char* strParamNames[]{
	"U1FName",
	"StartFolder",
	"ResultFolder",
	"ImprovedResultFolder",
	"UseBinaryCanonizer",
	"TestName",
	"MatrTest",
};

const char* arrayParamNames[]{
	"U1FCycles"
};

using namespace std;

bool checkInputParam(const kSysParam& param, const char** paramNames);

static size_t getValue(string& tmp, size_t* pPos)
{
	ltrim(tmp);
	if (!tmp.length()) {
		*pPos = string::npos;
		return 0;
	}

	if (tmp.front() != '"' || tmp.back() != '"') {
		size_t posE = tmp.find(' ');
		const size_t posE1 = tmp.find(';');
		if (posE1 != string::npos && posE1 < posE)
			posE = posE1;

		*pPos = posE != string::npos ? posE + 1 : string::npos;
		if (posE != string::npos)
			tmp = tmp.substr(0, posE);
	} else
		tmp = tmp.substr(1, tmp.length() - 2);

	return tmp.length();
}

static int getInteger(const string& str, size_t* pPos) {
	string tmp = str.substr(*pPos + 1);
	if (!getValue(tmp, pPos))
		return 1;

	if (isdigit(tmp[0]))
		return atoi(tmp.c_str());

	return tmp == "YES" ? 1 : 0;
}

static size_t find(const string& str, const char* strValue)
{
	const size_t pos = str.find(strValue);
	return pos != string::npos ? pos + strlen(strValue) : pos;
}

template<typename T>
void setDefaultValue(int* pValue) {
	*pValue = 1;
}

template<typename T>
void setDefaultValue(string** pValue) {
	if (*pValue)
		**pValue = "";
}

template<typename T>
void setDefaultValue(tchar **pValue) {
	printfRed("No default parameter defined for tchar **");
}

template <typename T>
bool setParameter(int *pValue, const string& str, size_t* pos) {
	*pValue = getInteger(str, pos);
	return true;
}

template <typename T>
bool setParameter(string** pValue, const string& str, size_t* pos) {
	string tmp = str.substr(*pos + 1);
	if (getValue(tmp, pos)) {
		if (!*pValue)
			*pValue = new string();

		**pValue = tmp;
	} else
		setDefaultValue<string**>(pValue);

	return true;
}

template <typename T>
bool setParameter(tchar** pValue, const string& str, size_t* pos) {
	auto pStr = str.c_str();
	auto* chptr = strrchr(pStr, '}');
	if (!chptr)
		return false;

	string* groups[10];
	auto last = chptr - pStr;
	chptr = strchr(pStr, '{');
	if (!chptr)
		return false;

	auto first = chptr - pStr;
	string tmp1 = str.substr(first + 1, last - first - 1);
	string tmp = regex_replace(tmp1, regex("\\s"), "");
	std::regex patternA("([0-9]+,)+[0-9]+");
	std::regex patternB("[0-9]+");
	int ngrp = 0;
	while (tmp.length()) {
		pStr = tmp.c_str();
		chptr = strchr(pStr, '{');
		if (!chptr)
			break;

		first = chptr - pStr;
		chptr = strchr(pStr, '}');
		if (!chptr)
			return false;

		last = chptr - pStr;
		assert(ngrp < countof(groups));
		groups[ngrp] = new string(tmp.substr(first + 1, last - first - 1));
		if (!std::regex_match(*groups[ngrp], patternA) && !std::regex_match(*groups[ngrp], patternB))
			return false;

		ngrp++;
		tmp = tmp.substr(last + 1);
	}

	delete[] *pValue;
	const auto len = MAX_CYCLES_PER_SET * ngrp;
	*pValue = new tchar[len + ngrp + 1];
	**pValue = ngrp;
	auto pTo = *pValue + 1;
	memset(pTo, 0, len);
	bool retVal = true;
	int val;
	for (int j = 0; j < ngrp; j++) {
		auto& tmp = *groups[j];
		trim(tmp);
		int k = 0;
		while (tmp.length()) {
			if (isdigit(tmp[0])) {
				pStr = tmp.c_str();
				chptr = strchr(pStr, ',');
				if (chptr) {
					val = atoi(tmp.substr(0, chptr - pStr).c_str());
					tmp = tmp.substr(chptr - pStr + 1);
					ltrim(tmp);
				}
				else {
					val = atoi(tmp.c_str());
					tmp.clear();
				}

				if (k < MAX_CYCLES_PER_SET && val) {
					*(pTo + k++) = val;
				}
				else
					retVal = false;
			}
			else
				retVal = false;
		}

		pTo += MAX_CYCLES_PER_SET;
		delete groups[j];
	}

	return retVal;
}

template <typename T>
static int getParam(const string& str, const char* pKeyWord, T* pValue, size_t* pPos = nullptr) {
	size_t pos = find(str, pKeyWord);
	if (pos == string::npos)
		return 0;		// keyWord was not found

	const auto flg = pos == str.length();
	if (!flg && str[pos] != ' ' && str[pos] != '=')
		return 0;

	string tmp = str.substr(0, pos);
	if (tmp != pKeyWord)
		return 0;

	tmp = str.substr(pos);
	ltrim(tmp);
	if (flg || tmp[0] != '=') {
		// the value is determined by the presence of the keyword
		setDefaultValue<T>(pValue);
		if (!tmp.length())
			return -1;	// there is nothing after keyWords
	}
	else {
		if (!setParameter<T>(pValue, str, &pos)) {
			printfRed("Parameter is not set as expected: '%s'\n", str.c_str());
			exit(1);
		}
	}

	if (pPos)
		*pPos = pos;

	return 1;
}

bool getParameters(ifstream& infile, const paramDescr *par, int nDescr, kSysParam &param, bool &endJob) {
	bool retVal = false;
	string line;
	while (getline(infile, line)) {		// For all the lines of the file
		trim(line);
		size_t pos = line.find("//");
		if (pos != string::npos)
			line = line.substr(0, pos);	// deleting a comment at the end of a line

		if (!line.size() || line[0] == ';')
			continue;					// Skip line if it is a comment OR empty

		if (line[0] == '/' && line[1] == '*') {
			// Commented out part of the input file
			// Let's find the closing 
			do {
				const size_t pos = line.find("*/");
				if (pos != string::npos) {
					line = line.substr(pos + 2);
					break;
				}
			} while (getline(infile, line));

			trim(line);
			if (!line.length())
				continue;

		}

		pos = line.find("=");
		if (pos != string::npos) {
			auto beg = line.substr(0, pos);
			auto end = line.substr(pos + 1);
			rtrim(beg);
			transform(beg.begin(), beg.end(), beg.begin(), ::toupper);
			ltrim(end);
			line = beg + "=" + end;
		}
		else
			transform(line.begin(), line.end(), line.begin(), ::toupper);

		if (line.find("END_JOB") != string::npos) {
			endJob = true;
			break;
		}

		if (line.find("RUN_JOB") != string::npos)
			break;
		
		auto j = nDescr;
		while (j--) {
			auto paramNames = (par+j)->paramNames;
			bool rc = false;
			int i = (par+j)->numParams;
			while (!rc && i--) {
				if (j == 0)
					rc = getParam<int>(line, paramNames[i], param.val + i);
				else
				if (j == 1)
					rc = getParam<string*>(line, paramNames[i], &param.strVal[i]);
				else
					rc = getParam<tchar*>(line, paramNames[i], &param.u1fCycles[i]);
			}

			if (i >= 0) {
				retVal = true;  // at least one parameter has been changed
				break;
			}
		}

		if (j < 0)
			printf("Sorry, the parameter for keyword '%s' was not found\n", line.c_str());
	}

	return retVal;
}
string getUF(const tchar * pU1F) {
	char group[32];
	tchar val;
	int k = 0;
	string str("{");
	auto nu = *pU1F++;
	while (nu--) {
		group[k] = '{';
		for (int i = 0; i < MAX_CYCLES_PER_SET && (val = pU1F[i]); i++) {
			int j = val < 10 ? 1 : val < 100 ? 10 : 100;
			do {
				group[++k] = '0' + val / j;
				val -= (val / j) * j;
			} while (j /= 10);

			group[++k] = ',';
		}
		pU1F += MAX_CYCLES_PER_SET;
		group[k] = '}';
		group[k + 1] = '\0';
		str.append(group);
		group[0] = ',';
		k = 1;
	}
	str.append("}");
	return str;
}
int factorial(int n) {
	return n > 2 ? n * factorial(n - 1) : n;
}
int main(int argc, const char* argv[])
{
	std::cout << "k - Sys 10.61\n";

	paramDescr par[] = { 
		intParamNames, countof(intParamNames), 
		strParamNames, countof(strParamNames),
		arrayParamNames, countof(arrayParamNames)
	};

	paramDescr params[] = {
		intParamNames, countof(intParamNames),
		strParamNames, countof(strParamNames),
		arrayParamNames, countof(arrayParamNames)
	};

	for (int j = 0; j < countof(par); j++) {
		auto* paramNames = par[j].paramNames;
		int i = par[j].numParams;
		auto* upperNames = new char* [i];
		par[j].paramNames = (const char**)upperNames;
		for (; i--;) {
			string line(paramNames[i]);
			transform(line.begin(), line.end(), line.begin(), ::toupper);
			const auto len = line.size() + 1;
			upperNames[i] = new char[len];
			strcpy_s(upperNames[i], len, line.c_str());
		}
	}

	kSysParam param;
	param.pParamDescr = params;
	// Set default integer parameters:
	auto& val = param.val;
	memset(val, 0, sizeof(val));
	val[t_numPlayers] = nPlayers;
	val[t_groupSize] = GroupSize;
	val[t_u1f] = UseUniform1Factorization;
	val[t_use2RowsCanonization] = Use2RowsCanonization;
	val[t_submatrixGroupOrderMin] = SubmatrixGroupOrderMin;
	val[t_resultGroupOrderMin] = ResultGroupOrderMin;
	val[t_useGPU] = USE_GPU;
	val[t_MultiThreading] = UseMultiThreading;
	val[t_numThreads] = NThreads;
	val[t_nRowsInStartMatrix] = NRowsInStartMatrix;
	val[t_nMaxNumberOfStartMatrices] = MaxNumberOfStartMatrices;
	val[t_nFirstIndexOfStartMatrices] = FirstIndexOfStartMatrices;
	val[t_nRowsInResultMatrix] = NRowsInResultMatrix;
	val[t_expectedResult] = ExpectedResult;
	val[t_useCheckLinksV] = UseCheckLinksV;
	val[t_printMatrices] = PrintMatrices;
	val[t_savingMatricesToDisk] = SavingMatricesToDisk;
	val[t_checkConstructedMatrices] = CheckConstructedMatrices;
	val[t_useSS] = UseSS;
	val[t_p1f_counter] = 5000;
	val[t_autLevelMinDef] = val[t_autLevelMin] = NRBase;
	val[t_useImproveMatrix] = UseImproveMatrix;
	val[t_gridSize] = 32;
	val[t_blockSize] = 24;
	val[t_any2RowsConvertToFirst2] = Any2RowsConvertToFirst2;

	// Set default string parameters:
	auto* strVal = param.strVal;
	memset(strVal, 0, t_lastStrParam * sizeof(strVal[0]));
	if (U1FName && strlen(U1FName))
		strVal[t_UFname] = new string(U1FName);

	strVal[t_StartFolder] = new string(StartFolder);
	strVal[t_ResultFolder] = new string(ResultFolder);
	strVal[t_ImprovedResultFolder] = new string(ImprovedResultFolder);

	for (int i = 0; i < countof(param.u1fCycles); i++)
		param.u1fCycles[i] = 0;

	string* testToRun = NULL;
	ifstream* infile = NULL;
	if (argc > 1) {
		infile = new ifstream(argv[1], ifstream::in);
		if (infile && !infile->is_open()) {
			printf("Cannot open file \"%s\"\n", argv[1]);
			delete infile;
			infile = NULL;
		}
		else
			if (argc > 2)
				testToRun = new string(argv[2]);
	}

	if (!infile)
		printf("The precompiled parameters will be used:\n");

	// Job is defined by external file
	vector<string> failedTests;
	auto& testName = strVal[t_testName];
	int numErrors = 0;
	int numTests = 0;
	bool endJob = false;
	do {
		if (infile) {
			if (!getParameters(*infile, par, countof(par), param, endJob))
				break;
		}

		if (!testToRun || (testName && *testToRun == *testName)) {
			const auto numPlayers = val[t_numPlayers];
			const auto groupSize = val[t_groupSize];
			const auto useGPU = val[t_useGPU];

			char buffer[64];
			if (testName)
				sprintf_s(buffer, "\"%s\" ", testName->c_str());
			else
				buffer[0] = 0;

			auto* pAutLevel = val + t_autLevelMinDef;
			for (int i = 2; i < 4; i++, pAutLevel++) {
				if (*pAutLevel < i)
					*pAutLevel = i;

				if (!*(++pAutLevel)) {
					*pAutLevel = val[t_nRowsInResultMatrix];
					if (!*pAutLevel)
						*pAutLevel = (numPlayers - 1) / (groupSize - 1);

					if (i == 2 && val[t_nRowsInResultMatrix] > 2)
						--*pAutLevel;
				}
				else {
					if (*pAutLevel < *(pAutLevel - 1))
						*pAutLevel = *(pAutLevel - 1);
				}
			}

			printfGreen("Test %sis launched with the following parameters:\n", buffer);
			if (val[t_u1f]) {
				if (!val[t_use2RowsCanonization])
					printfYellow(" 1 row canonization %s\n", val[t_groupSize] > 3 || val[t_numPlayers] < 12  ? "" : "(can be slow)");
				else
					printfGreen(" 2 rows canonization\n");

				const auto uf_code = param.u1fCycles[0];
				auto& ufName = strVal[t_UFname];
				if (uf_code) {
					const auto uf = getUF(uf_code);
					if (!ufName) {
						ufName = new string("");
						auto lenCycles = uf_code + 1;
						for (int i = 0; i < *uf_code; i++, lenCycles += MAX_CYCLES_PER_SET) {
							*ufName += '_';
							auto lenCycle = lenCycles;
							while (*lenCycle)
								*ufName += std::to_string(*lenCycle++);
						}
					}
					printfGreen(" %s-matrices with cycles(%s): %s ", getFileNameAttr(&param), ufName->c_str(), uf.c_str());
				}
				else {
					if (!ufName)
						ufName = new string("");

					printfGreen(" %s-matrices with p1f cycles(%d) ", getFileNameAttr(&param), val[t_numPlayers]);
				}
			}
			else {
				printfGreen("%s-matrices ", getFileNameAttr(&param));
			}

			printfGreen("for numPlayers = %d, groupSize = %d using ", numPlayers, groupSize);
			if (!useGPU && val[t_MultiThreading])
				printfGreen("%d CPU threads\n", val[t_numThreads]);
			else
				printfGreen("%cPU\n", useGPU ? 'G' : 'C');

			bool testOK = checkInputParam(param, intParamNames);
			if (testOK) {
				param.groupSizeFactorial = factorial(groupSize);
				TopGunBase* topGun;
				if (!useGPU)
					topGun = new TopGun(param);
				else
					topGun = new TopGunGPU(param);


				topGun->outputIntegratedResults(NULL, 0);
				if (topGun->Run())
					testOK = false;

				topGun->outputIntegratedResults(params, countof(par));
				delete topGun;
			}

			numTests++;
			if (!testOK) {
				numErrors++;
				failedTests.push_back(testName ? *testName : string("NO_NAME TEST"));
			}
		}

		delete testName;
		testName = NULL;
	} while (infile && !endJob);

	if (numErrors) {

		printfRed("The following %d tests failed:\n", numErrors);
		for (int i = 0; i < failedTests.size(); i++)
			printfRed("\"%s\"\n", failedTests[i].c_str());

		wchar_t buffer[64];
		swprintf_s(buffer, L"%d out of %d tests failed", numErrors, numTests);
		speakText(buffer);
	}

	printf("There's nothing left to do.\n");
	speakText(L"There's nothing left to do");

	for (int j = 0; j < countof(par); j++) {
		const auto* upperNames = par[j].paramNames;
		for (auto i = par[j].numParams; i--;)
			delete[] upperNames[i];

		delete[] upperNames;
	}

	for (int j = t_lastStrParam; j--;)
		delete strVal[j];

	for (int j = par[2].numParams; j--;)
		for (int i = 0; i < countof(param.u1fCycles); i++)
			delete [] param.u1fCycles[i];

	delete infile;
	delete testToRun;
	return numErrors;
}

