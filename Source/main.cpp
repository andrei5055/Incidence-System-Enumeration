// BIBD.cpp : Defines the entry point for the console application.
// 

#include "stdafx.h"
#include "C_tDesignEnumerator.h"

//#include <iostream> 
#include <fstream>
#include <string>      // for getline
//#include <cctype>
#include <algorithm>
//#include <numeric>
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
	t_BIBD,			// default
	t_tDesign,
	t_PBIBD,
	t_InsidenceSystem
} t_objectType;

typedef enum {
	t_Enumeration,	// default
	t_Automorphism,
	t_Isomorphism,
	t_Canonicity
} t_operationType;

static bool getBIBDParam(const string &paramText, designRaram *param)
{
	int i = 0;
	int j = 0;
	while (true) {
		int num = 0;
		while (true) {
			char symb = paramText[i++];
			if (symb == ',' || symb == '\0')
				break;

			if (symb == ' ') {
				if (!num)
					continue;

				break;
			}

			symb -= '0';
			if (symb < 0 || symb > 9)
				return false;

			num = 10 * num + symb;
		}
		if (j++ == 2) {
			param->lambda = num;
			break;
		} else
		if (j == 2)
			param->k = num;
		else
		   param->v = num;
	}

	return true;
}

static bool getTParam(const string &paramText, designRaram *param)
{
	param->t = 2;
	size_t len = paramText.length();
	if (len-- && paramText[len] != '-')
		return true;

	int num = 0;
	int mult = 1;
	while (len--) {
		char symb = paramText[len] - '0';
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

template <class T>
bool RunOperation(designRaram *pParam, const char *pSummaryFileName, bool FirstPath)
{
	if (pParam->v <= 0)
		return false;

	const string workingDir = pParam->workingDir + pSummaryFileName;
	const char *pSummFile = workingDir.c_str();
	InitCanonInfo(pParam->threadNumb);
	C_InSys<T> *pInSys = NULL;
	C_InSysEnumerator<T> *pInSysEnum = NULL;
	if (pParam->t <= 2) {
		pInSys = new C_BIBD<T>(pParam->v, pParam->k, 2, pParam->lambda);
		pInSysEnum = new CBIBD_Enumerator<T>(static_cast<C_BIBD<T> *>(pInSys), false, pParam->noReplicatedBlocks);
	}
	else {
		pInSys = new C_tDesign<T>(pParam->t, pParam->v, pParam->k, pParam->lambda);
		pInSysEnum = new  C_tDesignEnumerator<T>(static_cast<C_tDesign<T> *>(pInSys), false, pParam->noReplicatedBlocks);
	}

	char buff[256], buffer[256];
	MAKE_JOB_TITLE(pInSysEnum, buff, countof(buff));
	cout << buff;
	CInsSysEnumInfo<T> enumInfo(buff);
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
		pInSysEnum->Enumerate(pParam, PRINT_TO_FILE, &enumInfo);
		enumInfo.reportResult(buffer, countof(buffer));
		outString(buffer, pSummFile);
		cout << '\r' << buffer;
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
	size_t pos = str.find(strValue);
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
	ifstream infile;
	infile.open(argv[1]);
	if (!infile.is_open()) {
		printf("Cannot open file: \"%s\"\n", (const char *)argv[1]);
		abort();
	}

	bool firstRun = true;
	designRaram param;
	memset(&param, 0, sizeof(param));
	param.threadNumb = USE_THREADS;
	param.noReplicatedBlocks = false;;
	string newWorkDir = "./";

	// By default, we enumerating BIBDs 
	t_parsingStage stage = t_objectTypeStage;
	t_objectType objType = t_BIBD;
	t_operationType operType = t_Enumeration;
	uint outType = t_Summary;
	string *pLine = new string();
	string &line = *pLine;
	while (getline(infile, line)) {	// For all the lines of the file
		trim(line);
		if (line.size() <= 2 || line[0] == ';' || line[0] == '/' && line[1] == '/')
			continue;				// Skip line if it is a comment OR too short

		transform(line.begin(), line.end(), line.begin(), ::toupper);

		if (line.find("END_JOB") != string::npos)
			break;

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
		
		size_t pos = find(line, "WORKING_DIR");
		if (pos != string::npos) {
			newWorkDir = line.substr(pos + 1);
			if (newWorkDir.c_str()[newWorkDir.length() - 1] != '/')
				newWorkDir += '/';

			continue;
		}

		pos = find(line, "THREAD_NUMBER");
		if (pos != string::npos) {
			// Define the number of threads launched to perform task
			const size_t threadNumb = getInteger(line, &pos);
			if (threadNumb == string::npos) {
				printf("Cannot define thread number from: \"%s\"\n", line.c_str());
				printf("Will use the default value: %d\n", USE_THREADS);
				param.threadNumb = USE_THREADS;
			}
			else
				param.threadNumb = threadNumb;

			continue;
		}

		pos = find(line, "NO_REPLICATED_BLOCKS");
		if (pos != string::npos) {
			param.noReplicatedBlocks = getInteger(line, &pos);
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
		if (line.find("BIBD") != string::npos)
			objType = t_BIBD;
		else
		if (line.find("T-DESIGN") != string::npos)
			objType = t_tDesign;
		else
		if (line.find("INCIDENCE") != string::npos)
			objType = t_InsidenceSystem;

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

						param.grpOrder = (uint)getInteger(tmp, &pos);
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

		param.outType = outType;
		param.t = 2;
		size_t from = -1;
		switch (objType) {
		case t_tDesign:	if (!getTParam(line.substr(0, beg), &param)) {
							from = 0;
							break;
						}

		case t_BIBD:	if (!getBIBDParam(line.substr(beg + 1, end - beg - 1), &param))
							from = beg;

						break;
		case t_PBIBD:

		case t_InsidenceSystem: 
						break;
		}

		if (from != -1) {
			printf("Cannot read parameters: '%s'\n", line.substr(from, end - from + 1).c_str());
			break;
		}

		param.outType = outType;
		if (newWorkDir != param.workingDir || firstRun)
			param.workingDir = newWorkDir;

		if (!RunOperation<MATRIX_ELEMENT_TYPE>(&param, pSummaryFile, firstRun))
			break;

		firstRun = false;
	}

	delete pLine;
	infile.close();

	_CrtDumpMemoryLeaks();
	return 0;
}
// To 
// cuda - memcheck --leak - check full . / CDTools_GPU.exe
// nvprof CDTools_GPU.exe
