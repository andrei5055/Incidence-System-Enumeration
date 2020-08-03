// BIBD.cpp : Defines the entry point for the console application.
// 

#include "stdafx.h"
#include "C_tDesignEnumerator.h"
#include "CombBIBD_Enumerator.h"
#include "IG_Enumerator.h"

#include <fstream>
#include <string>      // for getline
#include <algorithm>
#include <iterator>
#include <functional>   // for ptr_fun

#define SDL_MAIN_HANDLED
using namespace std;

int find_T_designParam(int v, int k, int lambda)
{
	int lam = lambda;
	int prevLambda, lambdas[10];
	int i = 1;
	do {
		lambdas[++i - 2] = prevLambda = lambda;
		lambda = prevLambda * (k - i) / (v - i);
	} while (lambda && lambda * (v - i) / (k - i) == prevLambda);

    if (i > 2)
        printf("{%d, %3d, %2d, %2d},\n", i, v, k, lambdas[i - 2]);
    
	return i;
}

// trim from start (in place)
static inline void ltrim(string &s) {
 	s.erase(s.begin(), find_if(s.begin(), s.end(),
		not1(ptr_fun<int, int>(isspace))));
}

// trim from end (in place)
static inline void rtrim(string &s) {
	s.erase(find_if(s.rbegin(), s.rend(),
		not1(ptr_fun<int, int>(isspace))).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(string &s) {
	ltrim(s);
	rtrim(s);
}

typedef enum {
	t_objectTypeStage,
	t_objectParamStage,
	t_operationTypeStage,
	t_operationParamStage,
} t_parsingStage;

typedef enum {
	t_Enumeration,	// default
	t_Automorphism,
	t_Isomorphism,
	t_Canonicity
} t_operationType;

static bool getBIBDParam(const string &paramText, designParam *param, bool BIBD_flag = true)
{
	int i = 0;
	int j = 0;
	int flag = 0;
	char symb;
	do {
		int num = 0;
		while (true) {
			symb = paramText[i++];
			if (symb == ',' || symb == '\0' || j == 2 && flag && symb == '}')
				break;

			if (symb == ' ') {
				if (!num)
					continue;

				break;
			}
			else
			if (j == 2 && symb == '{') {
				// Only lambda can be set with multiple values
				if (flag || num)
					return false;

				flag = 1;
				continue;
			}

			symb -= '0';
			if (symb < 0 || symb > 9)
				return false;

			num = 10 * num + symb;
		}

		auto iStruct = param->InterStruct();
		if (j == 2) {
			iStruct->lambdaPtr()->push_back(num);
			if (symb == '}')
				return true;

			continue;
		}

		if (j++ == 2) {
			if (BIBD_flag) {
				iStruct->lambdaPtr()->push_back(num);
				return true;
			}

			param->r = num;
			continue;
		}

		if (j == 2)
			param->k = num;
		else
			param->v = num;
	} while (symb != '\0');

	return true;
}

static bool getTParam(const string &paramText, designParam *param)
{
	param->t = 2;
	size_t len = paramText.length();
	if (len-- && paramText[len] != '-')
		return true;

	int num = 0;
	int mult = 1;
	while (len--) {
		const char symb = paramText[len] - '0';
		if (symb < 0 || symb > 9)
			break;

		num += mult * symb;
		mult *= 10;
	}

	if (mult == 1)
		return false;

	param->t = num;
	return true;
}

template <typename T, typename S>
bool RunOperation(designParam *pParam, const char *pSummaryFileName, bool FirstPath)
{
	if (pParam->v <= 0)
		return false;

	const string workingDir = pParam->workingDir + pSummaryFileName;
	const char *pSummFile = workingDir.c_str();
	InitCanonInfo(pParam->threadNumb);
	Class2(C_InSys) *pInSys = NULL;
	Class2(C_InSysEnumerator) *pInSysEnum = NULL;
	t_objectType objType = pParam->objType;
	uint enumFlags = pParam->noReplicatedBlocks ? t_noReplicatedBlocks : t_enumDefault;
	if ((pParam->outType & t_GroupGeneratingSet) == t_GroupGeneratingSet)
		enumFlags |= t_outRowOrbits + t_outRowPermute;
	else
	if (pParam->outType & t_GroupOrbits)
		enumFlags |= t_outRowOrbits;

	const auto lambda = pParam->InterStruct()->lambda();
	if (pParam->t <= 2) {
		if (lambda.size() > 1 || objType == t_SemiSymmetricGraph) {
			switch (objType) {
			case t_PBIBD:
				pInSys = new Class2(C_PBIBD)(pParam->v, pParam->k, pParam->r, lambda);
				pInSysEnum = new Class2(CPBIBD_Enumerator)(pInSys, enumFlags);
				break;
			case t_SemiSymmetricGraph:
				enumFlags |= t_outColumnOrbits + t_outStabilizerOrbit + t_colOrbitsConstructed + t_alwaysKeepRowPermute;
				pInSys = new Class2(CSemiSymmetricGraph)(pParam->v, pParam->k, pParam->r, lambda);
				pInSysEnum = new Class2(CIG_Enumerator)(pInSys, pParam, enumFlags, FirstPath);
				break;
			case t_CombinedBIBD:
				pInSys = new Class2(CCombinedBIBD)(pParam->v, pParam->k, lambda);
				pInSysEnum = new Class2(CCombBIBD_Enumerator)(pInSys, enumFlags);
				break;
			default: return false;
			}
		}
		else {
			pInSys = new Class2(C_BIBD)(pParam->v, pParam->k, 2, lambda[0]);
			pInSysEnum = new Class2(CBIBD_Enumerator)(pInSys, enumFlags);
			objType = t_BIBD;
		}
	}
	else {
		pInSys = new Class2(C_tDesign)(pParam->t, pParam->v, pParam->k, lambda[0]);
		pInSysEnum = new Class2(C_tDesignEnumerator)(static_cast<TDesignPntr>(pInSys), enumFlags);
		objType = t_tDesign;
	}

	pInSys->setObjectType(objType);

	char buff[256], buffer[256];
	MAKE_JOB_TITLE(pInSysEnum, pParam, buff, countof(buff));
	cout << buff;
	Class2(CInsSysEnumInfo) enumInfo(buff);
	enumInfo.setDesignInfo(pParam);
	if (FirstPath) {
		FOPEN(outFile, pSummFile, "w");
		enumInfo.outEnumInformation(&outFile, false);
		outString("         BIBDs:                     Canonical:      NRB #:      Constructed:    Run Time (sec):\n", pSummFile);
	}

	const bool resetMTlevel = pParam->mt_level == 0;
	if (resetMTlevel) {
		// The row number, on which the threads will be launched was not defined.
		// Let's do it here by other parameters 
		pParam->mt_level = pInSysEnum->define_MT_level(pParam);
	}

	try {
		if (pInSysEnum->Enumerate(pParam, PRINT_TO_FILE, &enumInfo)) {
			enumInfo.reportResult(buffer, countof(buffer));
			outString(buffer, pSummFile);
			cout << '\r' << buffer;
		}
		else {
			cout << "Some problem durring the enumeration found\n";
			pInSysEnum->closeFile();
		}
	}
	catch (...) {
		pInSysEnum->closeFile();
	}

	if (resetMTlevel)
		pParam->mt_level = 0;

	delete pInSys;
	delete pInSysEnum;
	CloseCanonInfo();
	return true;
}

static size_t find(const string &str, const char *strValue)
{
	const size_t pos = str.find(strValue);
	return pos != string::npos ? pos + strlen(strValue) : pos;
}

static size_t getInteger(const string &str, size_t *pPos) {
	string tmp = str.substr(*pPos + 1);
	ltrim(tmp);
	if (!tmp.length())
		return *pPos = string::npos;

	size_t posE = tmp.find(' ');
	const size_t posE1 = tmp.find(';');
	if (posE1 != string::npos && posE1 < posE)
		posE = posE1;

	*pPos = posE != string::npos ? posE + 1 : string::npos;
	if (posE != string::npos)
		tmp = tmp.substr(0, posE);

	if (isdigit(tmp[0]))
		return atoi(tmp.c_str());

	return tmp == "YES"? 1 : 0;
}

int main(int argc, char * argv[])
{
	_CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_FILE);
	_CrtSetReportFile(_CRT_WARN, _CRTDBG_FILE_STDOUT);
	_CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_FILE);
	_CrtSetReportFile(_CRT_ERROR, _CRTDBG_FILE_STDOUT);
	_CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_FILE);
	_CrtSetReportFile(_CRT_ASSERT, _CRTDBG_FILE_STDOUT);

	const char *pSummaryFile = "EnumSummary.txt";
	cudaDeviceReset();
	/*
	int nDevices;
	cudaError_t err = cudaGetDeviceCount(&nDevices);
	for (int i = 0; i < nDevices; i++) {
	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, i);
	printf("Device Number: %d\n", i);
	printf("  Device name: %s\n", prop.name);
	printf("  ConcurrentManagedAccess: %d\n", prop.concurrentManagedAccess);
	printf("  Memory Clock Rate (KHz): %d\n", prop.memoryClockRate);
	printf("  Memory Bus Width (bits): %d\n", prop.memoryBusWidth);
	printf("  Peak Memory Bandwidth (GB/s): %f\n\n",
	2.0*prop.memoryClockRate*(prop.memoryBusWidth / 8) / 1.0e6);
	}
	*/

	// Limit of the heap size used for memory allocation on GPU	
	cudaDeviceSetLimit(cudaLimitMallocHeapSize, 256 * 1024 * 1024);

	if (argc == 1) {
		printf("Program \"%s\" launched without parameters\n", (const char *)argv[0]);
		abort();
	}

	// Job is defined by external file
	ifstream infile(argv[1], ifstream::in);
	if (!infile.is_open()) {
		printf("Cannot open file: \"%s\"\n", argv[1]);
		abort();
	}

	bool firstRun = true;
	designParam *param = new designParam();
	param->threadNumb = USE_THREADS;
	param->noReplicatedBlocks = false;
	std::string newWorkDir = "./";

	// By default, we enumerating BIBDs 
	t_parsingStage stage = t_objectTypeStage;
	t_objectType objType = t_BIBD;
	t_operationType operType = t_Enumeration;
	uint outType = t_Summary;
	string *pLine = new string();
	string &line = *pLine;
	while (getline(infile, line)) {	// For all the lines of the file
		trim(line);
		if (!line.size() || line[0] == ';' || line[0] == '/' && line[1] == '/')
			continue;				// Skip line if it is a comment OR empty

		transform(line.begin(), line.end(), line.begin(), ::toupper);

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

			transform(line.begin(), line.end(), line.begin(), ::toupper);
		}

		if (line.find("END_JOB") != string::npos)
			break;

		size_t pos = find(line, "WORKING_DIR");
		if (pos != string::npos) {
			newWorkDir = line.substr(pos + 1);
			std::replace(newWorkDir.begin(), newWorkDir.end(), '\\', '/');
			if (newWorkDir.c_str()[newWorkDir.length() - 1] != '/')
				newWorkDir += '/';

			struct stat sb;
			const auto* pWorkDir = newWorkDir.c_str();
			const char* pCause = NULL;
			if (stat(pWorkDir, &sb))
				pCause = "get information about";
			else
				if (!(S_IFDIR & sb.st_mode))
					pCause = "find";
				else
					if (!(S_IREAD & sb.st_mode))
						pCause = "read from";
					else
						if (!(S_IWRITE & sb.st_mode))
							pCause = "write into";

			if (!pCause)
				continue;

			printf("Cannot %s working directory: \'%s\'\n", pCause, pWorkDir);
			break;
		}

		pos = find(line, "THREAD_NUMBER");
		if (pos != string::npos) {
			// Define the number of threads launched to perform task
			const size_t threadNumb = getInteger(line, &pos);
			if (threadNumb == string::npos) {
				printf("Cannot define thread number from: \"%s\"\n", line.c_str());
				printf("Will use the default value: %d\n", USE_THREADS);
				param->threadNumb = USE_THREADS;
			}
			else
				param->threadNumb = threadNumb;

			continue;
		}

		pos = pos = find(line, "THREAD_LEVEL");
		if (pos != string::npos) {
			// Define the row number, where threads will be launched
			const auto mt_level = getInteger(line, &pos);
			if (mt_level == string::npos) {
				printf("Cannot define thread level from: \"%s\"\n", line.c_str());
				printf("Will use the default calculated from object's parameter: v / 2\n");
			}
			else
				param->mt_level = static_cast<int>(mt_level);
		}
		pos = find(line, "NO_REPLICATED_BLOCKS");
		if (pos != string::npos) {
			param->noReplicatedBlocks = getInteger(line, &pos);
			if (pos == string::npos)
				continue;
			line = line.substr(pos);
		}
		// Define a job type
		if (line.find("ENUMERATION") != string::npos)
			operType = t_Enumeration;
		else
		if (line.find("AUTOMORPHISM") != string::npos)
			operType = t_Automorphism;
		else
		if (line.find("ISOMORPHISM") != string::npos)
			operType = t_Isomorphism;
		else
		if (line.find("CANONICITY") != string::npos)
			operType = t_Canonicity;

		if (line.find("PBIBD") != string::npos)
			objType = t_PBIBD;
		else
		if (line.find("COMBINED_BIBD") != string::npos)
			objType = t_CombinedBIBD;
		else
		if (line.find("BIBD") != string::npos)
			objType = t_BIBD;
		else
		if (line.find("T-DESIGN") != string::npos)
			objType = t_tDesign;
		else
		if (line.find("INCIDENCE") != string::npos)
			objType = t_IncidenceSystem;
		else
		if (line.find("SEMI_SYMMETRIC_GRAPH") != string::npos)
			objType = t_SemiSymmetricGraph;

		// Define output type
		pos = find(line, "OUTPUT");
		if (pos != string::npos) {
			outType = t_Summary;
			string output = line.substr(pos);
			if (output.find("NO SUMMARY") != string::npos)
				outType &= (-1 ^ t_Summary);
			else
			if (output.find("SUMMARY") != string::npos)
				outType |= t_Summary;

			if (output.find("ORBIT") != string::npos)
				outType |= t_GroupOrbits;

			if (output.find("GENERATING SET") != string::npos)
				outType |= t_GroupGeneratingSet;

			if (output.find("ALL OBJECTS") != string::npos)
				outType |= t_AllObject;
			else
			if (output.find("ALL TRANSITIVE") != string::npos)
				outType |= t_Transitive;
			else {
				pos = find(output, "GROUP ORDER");
				if (pos != string::npos) {
					string tmp = output.substr(pos);
					if ((pos = tmp.find(">")) != string::npos)
						outType |= t_GroupOrderGT;
					else
					if ((pos = tmp.find("<")) != string::npos)
						outType |= t_GroupOrderLT;

					if (pos != string::npos && tmp[++pos] == '=' ||
						pos == string::npos && (pos = tmp.find("=")) != string::npos) {
						outType |= t_GroupOrderEQ;
					}

					if (!(outType & (t_GroupOrderEQ | t_GroupOrderGT | t_GroupOrderLT))) {
						printf("Cannot define group order from: \"%s\"\n", output.c_str());
						printf("Will use the default output\n");
						outType = t_Summary;
					}
					else {
						if (!(outType & (t_GroupOrderEQ | t_GroupOrderGT)))
							pos++;

						param->grpOrder = (uint)getInteger(tmp, &pos);
					}
				}
			}
		}

		// Define  parameter values
		const size_t beg = line.find("(");
		if (beg == string::npos)
			continue;

		const size_t end = line.find(")");
		if (end == string::npos || end < beg) {
			printf("Wrong expression: '%s'\n", line.substr(end == string::npos? beg : end).c_str());
			break;
		}

		param->outType = outType;
		param->t = 2;
		param->InterStruct()->lambdaPtr()->resize(0);
		size_t from = -1;
		bool BIBD_flag = false;
		switch (objType) {
		case t_tDesign:	if (!getTParam(line.substr(0, beg), param)) {
							from = 0;
							break;
						}

		case t_BIBD:	BIBD_flag = true;
		case t_CombinedBIBD:
		case t_SemiSymmetricGraph:
		case t_PBIBD:	if (!getBIBDParam(line.substr(beg + 1, end - beg - 1), param, BIBD_flag))
							from = beg + 1;

						break;

		case t_IncidenceSystem: 
						break;
		}

		if (from != -1) {
			printf("Can't read parameters: '%s'\n", line.substr(from, end - from + 1).c_str());
			break;
		}

		if (param->workingDir != newWorkDir || firstRun)
			param->workingDir = string(newWorkDir);

		if ((param->objType = objType) == t_SemiSymmetricGraph) {
			int InconsistentGraphs(designParam *pParam, const char *pSummaryFileName, bool firstPath);
			InconsistentGraphs(param, pSummaryFile, firstRun);
		}
		else
			if (!RunOperation<TDATA_TYPES>(param, pSummaryFile, firstRun))
				break;

		firstRun = false;
	}

	delete pLine;
	infile.close();

	delete param;
	_CrtDumpMemoryLeaks();
	return 0;
}
// To 
// cuda - memcheck --leak - check full . / CDTools_GPU.exe
// nvprof CDTools_GPU.exe
// C program to find GCD of two numbers
