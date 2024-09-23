// BIBD.cpp : Defines the entry point for the console application.
// 

#include "stdafx.h"
#include "C_tDesignEnumerator.h"
#include "CombBIBD_Enumerator.h"
#include "IG_Enumerator.h"
#include "TripleSys.h"
#include "designParam.h"

#include <fstream>
#include <string>      // for getline
#include <algorithm>
#include <iterator>

#define SDL_MAIN_HANDLED
using namespace std;

const char *obj_name[] = {
	"BIBD",					// t_BIBD,			- default
	"COMBINED_BIBD",		// t_CombinedBIBD,
	"KIRKMAN_TRIPLE_SYSTEM",// t_Kirkman_Triple
	"TRIPLE_SYSTEM",        // t_TripleSystem - construction by Leo's program
	"T-DESIGN",				// t_tDesign,
	"PBIBD",				// t_PBIBD,
	"INCIDENCE",            // t_IncidenceSystem,
	"SEMI_SYMMETRIC_GRAPH", // t_SemiSymmetricGraph
	"CANON_MATR"			// t_CanonMatr
};

// The order of checking object names used during parameter parsing.
t_objectType idx_obj_type[] = {
	t_objectType::t_PBIBD,
	t_objectType::t_CombinedBIBD,
	t_objectType::t_BIBD,
	t_objectType::t_tDesign,
	t_objectType::t_IncidenceSystem,
	t_objectType::t_SemiSymmetricGraph,
	t_objectType::t_Kirkman_Triple,
	t_objectType::t_TripleSystem,
	t_objectType::t_CanonMatr,
};

int find_T_designParam(int v, int k, int lambda)
{
	int lam = lambda;
	int prevLambda(0), lambdas[10] = {};
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
static inline void ltrim(string& s) {
	s.erase(s.begin(), find_if(s.begin(), s.end(), 
		[](int c) { return !isspace(c); }));
}

// trim from end (in place)
static inline void rtrim(string& s) {
	s.erase(find_if(s.rbegin(), s.rend(),
		[](int c) {return !isspace(c); }).base(), s.end());
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
	char symb(0);
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

		auto lambdaSet = param->InterStruct()->lambdaPtr();
		if (j == 2) {
			lambdaSet->push_back(num);
			if (symb == '}')
				return true;

			continue;
		}

		if (j++ == 2) {
			if (BIBD_flag) {
				lambdaSet->push_back(num);
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

static size_t find(const string &str, const char *strValue)
{
	const size_t pos = str.find(strValue);
	return pos != string::npos ? pos + strlen(strValue) : pos;
}

static int getInteger(const string &str, size_t *pPos) {
	string tmp = str.substr(*pPos + 1);
	ltrim(tmp);
	if (!tmp.length()) {
		*pPos = string::npos;
		return 1;
	}

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

template <typename T>
static int getIntegerParam(const string& str, const char* pKeyWord, T *pValue, size_t* pPos = nullptr) {
	size_t pos = find(str, pKeyWord);
	if (pos == string::npos)
		return 0;		// keyWord was not found

	if (str[pos] != '=') {
		// the value is determined by the presence of the keyword
		*pValue = 1;	
		string tmp = str.substr(pos);
		ltrim(tmp);
		if (!tmp.length())
			return -1;	// there is nothing after keyWords
	} else
		*pValue = static_cast<T>(getInteger(str, &pos));

	if (pPos)
		*pPos = pos;

	return 1;
}

static int getBooleanParam(const string& str, const char* pKeyWord, uint flag, uint *pEnumFlags, size_t* pPos = nullptr) {
	size_t pos = find(str, pKeyWord);
	if (pos == string::npos)
		return 0;		// keyWord was not found

	if (str[pos] != '=') {
		// the value is determined by the presence of the keyword
		*pEnumFlags |= flag;	
		string tmp = str.substr(pos);
		ltrim(tmp);
		if (!tmp.length())
			return -1;	// there is nothing after keyWords
	}
	else {
		if (getInteger(str, &pos))
			*pEnumFlags |= flag;
		else
			*pEnumFlags &= (-1) ^ flag;
	}

	if (pPos)
		*pPos = pos;

	return 1;
}

void getMasterBIBD_param(const string& line, const size_t lineLength, designParam* param) {
	size_t pos = find(line, "THREAD_MASTER_BIBD");
	if (pos != string::npos) {
		const auto val = lineLength > pos ? getInteger(line, &pos) : string::npos;
		param->thread_master_DB = val != string::npos ? static_cast<int>(val) : 1;
	}
}

const char* pSummaryFile = "EnumSummary.txt";

bool parsingPath(string& line, size_t pos, bool isDir = true) {
	line = line.substr(pos + 1);
	std::replace(line.begin(), line.end(), '\\', '/');
	if (isDir && line.back() != '/')
		line += '/';

	struct stat sb;
	const auto* pWorkDir = line.c_str();
	const char* pCause = NULL;
	if (stat(pWorkDir, &sb))
		pCause = isDir? "get information about" : "find";
	else
	if (isDir && !(S_IFDIR & sb.st_mode))
		pCause = "find";
	else
	if (!(S_IREAD & sb.st_mode))
		pCause = "read from";
	else
	if (!(S_IWRITE & sb.st_mode))
		pCause = "write into";

	if (!pCause)
		return true;

	printf("Cannot %s %s: \'%s\'\n", pCause, isDir? "working directory" : "file", pWorkDir);
	return false;
}

int main(int argc, char * argv[])
{
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	_CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_FILE);
	_CrtSetReportFile(_CRT_WARN, _CRTDBG_FILE_STDOUT);
	_CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_FILE);
	_CrtSetReportFile(_CRT_ERROR, _CRTDBG_FILE_STDOUT);
	_CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_FILE);
	_CrtSetReportFile(_CRT_ASSERT, _CRTDBG_FILE_STDOUT);
//	_CrtSetBreakAlloc(194);  // Put here the memory allocation number you want to stop at.

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
	auto lambdaSet = param->InterStruct()->lambdaPtr();

	// By default, we enumerating BIBDs 
	t_parsingStage stage = t_objectTypeStage;
	t_objectType objType = t_objectType::t_BIBD;
	t_operationType operType = t_Enumeration;
	uint outType = t_Summary;
	string *pLine = new string();
	string &line = *pLine;
	string* pWorkDir = new string("./");
	string& newWorkDir = *pWorkDir;
	int use_master_sol = 0;
	int find_master_design = 0;
	int find_all_2_decomp = 0;
	int use_3_element_condition = 0;
	while (getline(infile, line)) {		// For all the lines of the file
		trim(line);
		size_t pos = line.find("//");
		if (pos != string::npos)
			line = line.substr(0, pos);	// deleting a comment at the end of a line

		if (!line.size() || line[0] == ';')
			continue;					// Skip line if it is a comment OR empty

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

		const auto length = line.length();
		pos = find(line, "WORKING_DIR");
		if (pos != string::npos) {
			if (!parsingPath(line, pos))
				break;
			newWorkDir = line;
			continue;
		}

		pos = find(line, "INPUT_FILE");
		if (pos != string::npos) {
			if (!parsingPath(line, pos, false))
				break;

			param->logFile = line;
			continue;
		}

		switch (getIntegerParam(line, "THREAD_NUMBER", &param->threadNumb)) {
		case 0: break;
		case -1:
			printf("Cannot define thread number from: \"%s\"\n", line.c_str());
			printf("Will use the default value: %d\n", USE_THREADS);
			param->threadNumb = USE_THREADS;
		default:
#if !USE_MUTEX
			param->threadNumb = USE_THREADS;
#endif
			continue;
		}

#if WAIT_THREADS
		// the option USE_MASTER_SOLUTIONS coud be used only when master waits for all threads finish on rowMaster() level
		pos = find(line, "USE_MASTER_SOLUTIONS");
		if (pos != string::npos) {
			// When 1, the solutions obtained by master will be used in the threads
			const auto val = length > pos ? getInteger(line, &pos) : string::npos;
			use_master_sol = val != string::npos ? static_cast<int>(val) : 1;
			if (use_master_sol > 1 || use_master_sol < 0) {
				printf("USE_MASTER_SOLUTIONS should be 0 or 1, got %d\n", use_master_sol);
				printf("Will use the default value 1\n");
				use_master_sol = 1;
			}
			param->use_master_sol = use_master_sol;
			continue;
		}
#endif
		if (getIntegerParam(line, "FIND_MASTER_BIBD", &find_master_design)) {
			if (find_master_design) {
				// When 1, the parts of Combined BIBD will be merged and canonical BIBD will be constructed
				// After that we will try to find that "original" BIBD in the ordered list of previously constructed
				// "original" BIBDs. In that way we will find all BIBD which could be split into several parts
				// and all such non-isomorphic splits
				if (find_master_design)
					getMasterBIBD_param(line, length, param);
				continue;
			}
		}

		getMasterBIBD_param(line, length, param);

		size_t mt_level;
		switch (getIntegerParam(line, "THREAD_LEVEL", &mt_level)) {
			case -1:
				printf("Cannot define thread level from: \"%s\"\n", line.c_str());
				printf("Will use the default calculated from object's parameter: v / 2\n");
			case 0: break;
			default: param->set_MT_level(static_cast<int>(mt_level));
		}

		getIntegerParam(line, "SAVE_RESTART_INFO", &param->save_restart_info);
		getIntegerParam(line, "FIND_ALL_2_DECOMPOSITIONS", &find_all_2_decomp);
		getIntegerParam(line, "COMPRESS_MATRICES", &param->m_compress_matrices);
		getIntegerParam(line, "NO_REPLICATED_BLOCKS", &param->noReplicatedBlocks);
		getIntegerParam(line, "USE_THREAD_POOL", &param->m_bUseThreadPool);
		getBooleanParam(line, "USE_SYMMETRICAL_T-CONDITION", t_symmetrical_t_cond, param->enumFlagsPtr());
		getBooleanParam(line, "USE_3-ELEMENT_CONDITION", t_use_3_condition, param->enumFlagsPtr());
		getBooleanParam(line, "UPDATE_RESULTS", t_update_results, param->enumFlagsPtr());

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

		for (int i = 0; i < countof(idx_obj_type); i++) {
			const auto type = idx_obj_type[i];
			if (line.find(obj_name[static_cast<int>(type)]) != string::npos) {
				objType = type;
				pos = line.find("=");
				if (pos != string::npos) {
					param->objSubType = line.substr(pos + 1);
					trim(param->objSubType);
					line.clear();
				}
				break;
			}
		}

		// Define output type
		pos = find(line, "OUTPUT");
		if (pos != string::npos) {
			outType = t_Summary;
			string output = line.substr(pos);
			bool val;
			if (getIntegerParam(output, "ONLY_SIMPLE_DESIGN", &val))
				param->setPrintOnlySimpleDesigns(val);

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

			if (output.find("COMB_MASTERS") != string::npos)
				outType |= t_CombMasters;   // Make output of Master BIBDs for Combined BIBDs
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
		lambdaSet->resize(0);
		size_t from = -1;
		bool BIBD_flag = false;
		switch (objType) {
		case t_objectType::t_tDesign:
						if (!getTParam(line.substr(0, beg), param)) {
							from = 0;
							break;
						}

		case t_objectType::t_BIBD:	BIBD_flag = true;
		case t_objectType::t_CombinedBIBD:
		case t_objectType::t_SemiSymmetricGraph:
		case t_objectType::t_PBIBD:
		case t_objectType::t_Kirkman_Triple:
		case t_objectType::t_TripleSystem:
		case t_objectType::t_CanonMatr:
						if (!getBIBDParam(line.substr(beg + 1, end - beg - 1), param, BIBD_flag))
							from = beg + 1;

						break;

		case t_objectType::t_IncidenceSystem:
						break;
		}

		if (from != -1) {
			printf("Can't read parameters: '%s'\n", line.substr(from, end - from + 1).c_str());
			break;
		}

		param->find_all_2_decomp = 0;
		if (param->workingDir != newWorkDir || firstRun)
			param->workingDir = newWorkDir;

		//		if (param->save_restart_info) {
		//		}
		param->setEnumFlags();

		if ((param->objType = objType) == t_objectType::t_SemiSymmetricGraph) {
			int InconsistentGraphs(designParam * pParam, const char* pSummaryFileName, bool firstPath);
			InconsistentGraphs(param, pSummaryFile, firstRun);
		}
		else
		if (objType == t_objectType::t_CanonMatr) {
			param->LaunchCanonization();
		}
		else {
			if (!param->LaunchEnumeration(pSummaryFile, find_master_design, find_all_2_decomp, use_master_sol, firstRun))
				continue;
		}

		firstRun = false;
	}

	delete pLine;
	delete pWorkDir;
	infile.close();

	delete param;
	_CrtDumpMemoryLeaks();
	return 0;
}

// To 
// cuda - memcheck --leak - check full . / CDTools_GPU.exe
// nvprof CDTools_GPU.exe
