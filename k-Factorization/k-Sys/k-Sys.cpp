// k-Sys.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <fstream>
#include "TopGun.h"

const char* intParamNames[]{
	"nPlayers",
	"GroupSize",
	"CBMP_Graph",           // Complete Balanced Multi-Partite Graph
	"UseUniform1Factorization",
	"Use2RowsCanonization",
	"UseFastCanonizerForG2", // 0 currently does not work with u1f
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
	"OutAutomorphismGroup",  // 0 - No output, 
	                         // otherwise any combinations of:
	                         //  1,  2 - generators and groups acting on elements
							 //  4,  8 - generators and groups acting on the matrix rows
							 // 16, 32 - generators ogroup acting on the k-sets
	"NestedGroups",
	"GridSize",
	"BlockSize",
	"OrderMatrices",
	"AllowUndefinedCycles", // 1  - allow rows pairs with cycles not defined in input params.
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

bool getParameters(ifstream& infile, const paramDescr* par, int nDescr, kSysParam& param, bool& firstSet, bool& endJob);
bool checkInputParam(const kSysParam& param, const char** paramNames);
string getUF(const tchar* pU1F);

int factorial(int n) {
	return n > 2 ? n * factorial(n - 1) : n;
}

void kSysParam::setup() {
	if (!val[t_CBMP_Graph])
		val[t_CBMP_Graph] = 1;
	groupSizeFactorial = factorial(val[t_groupSize]);
	m_numFactors = (val[t_numPlayers] - (completeGraph() ? 1 : partitionSize())) / (val[t_groupSize] - 1);
}

int main(int argc, const char* argv[])
{
	std::cout << "k - Sys 10.61\n";

	paramDescr params[] = {
		intParamNames, countof(intParamNames),
		strParamNames, countof(strParamNames),
		arrayParamNames, countof(arrayParamNames)
	};

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
	//if (U1FName)// && strlen(U1FName))
	//	strVal[t_UFname] = new string(U1FName);

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
	bool firstSet = true;
	bool endJob = false;
	do {
		if (infile) {
			if (!getParameters(*infile, params, countof(params), param, firstSet, endJob))
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
						*pAutLevel = val[t_CBMP_Graph] > 1 ? numPlayers / groupSize : (numPlayers - 1) / (groupSize - 1);

					if (i == 2 && val[t_nRowsInResultMatrix] > 2)
						--*pAutLevel;
				}
				else {
					if (*pAutLevel < *(pAutLevel - 1))
						*pAutLevel = *(pAutLevel - 1);
				}
			}

			param.setup();
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
						if (param.val[t_allowUndefinedCycles])
							*ufName += "_all";
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
				TopGunBase* topGun;
				if (useGPU)
					topGun = new TopGunGPU(param);
				else
					topGun = new TopGun(param);

				topGun->outputIntegratedResults(NULL, 0);
				if (topGun->Run())
					testOK = false;

				topGun->outputIntegratedResults(params, countof(params));
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

	for (int j = t_lastStrParam; j--;)
		delete strVal[j];

	for (int j = params[2].numParams; j--;)
		for (int i = 0; i < countof(param.u1fCycles); i++)
			delete [] param.u1fCycles[i];

	delete infile;
	delete testToRun;
	delete TopGun::secondRowDB();
	return numErrors;
}

