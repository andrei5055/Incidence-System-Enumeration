// BIBD.cpp : Defines the entry point for the console application.
// 

#include "stdafx.h"
#include "C_tDesignEnumerator.h"
#include "CombBIBD_Enumerator.h"
#include "IG_Enumerator.h"
#include "TripleSys/TripleSys.h"

#include <fstream>
#include <string>      // for getline
#include <algorithm>
#include <iterator>
#include <functional>   // for ptr_fun

#define SDL_MAIN_HANDLED
using namespace std;

const char *obj_name[] = {
	"BIBD",					// t_BIBD,			- default
	"COMBINED_BIBD",		// t_CombinedBIBD,
	"KIRKMAN_TRIPLE_SYSTEM",// t_Kirkman_Triples
	"TRIPLE_SYSTEM",        // t_TripleSystem - construction by Leo's program
	"T-DESIGN",				// t_tDesign,
	"PBIBD",				// t_PBIBD,
	"INCIDENCE",            // t_IncidenceSystem,
	"SEMI_SYMMETRIC_GRAPH"  // t_SemiSymmetricGraph
};

// Order of checking object names during parsing of parameters
t_objectType idx_obj_type[] = {
	t_objectType::t_PBIBD,
	t_objectType::t_CombinedBIBD,
	t_objectType::t_BIBD,
	t_objectType::t_tDesign,
	t_objectType::t_IncidenceSystem,
	t_objectType::t_SemiSymmetricGraph,
	t_objectType::t_Kirkman_Triple,
	t_objectType::t_TripleSystem,
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

template <typename T, typename S>
void output_2_decompInfo(designParam *param, const CDesignDB *pDesignDB, string& outputInfo, bool addInfo = false, const char* pSummaryFileName = NULL) {
	const auto& lambda = param->InterStruct()->lambda();
	Class2(C_BIBD) bibd(param->v, param->k, 2, lambda[0] + lambda[1]);
	Class2(CBIBD_Enumerator) bibdEnum(&bibd, t_enumDefault);
	bibdEnum.outNonCombinedDesigns(param, pDesignDB, outputInfo, addInfo);
	if (!addInfo) {
		char buffer[256];
		param->enumInfo()->reportResult(buffer, countof(buffer));
		outString(buffer, pSummaryFileName);
		cout << '\r' << buffer;
	}
}

template <typename T, typename S>
bool RunOperation(designParam *pParam, const char *pSummFile, bool FirstPath, std::string* outInfo) {
	if (pParam->v <= 0) {
		printf("Problem in RunOperation. Number of elements v = %d <= 0", pParam->v);
		return false;
	}

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

	const auto k = pParam->k;
	const auto& lambda = pParam->InterStruct()->lambda();
	if (pParam->t <= 2) {
		if (lambda.size() > 1 || objType == t_objectType::t_SemiSymmetricGraph) {
			bool kirkmanTriples = false;
			static t_objectType obj_types[] = { t_objectType::t_PBIBD, t_objectType::t_CombinedBIBD, t_objectType::t_SemiSymmetricGraph };
			switch (objType) {
			case t_objectType::t_PBIBD:
				pInSys = new Class2(C_PBIBD)(pParam->v, k, pParam->r, lambda);
				pInSysEnum = new Class2(CPBIBD_Enumerator)(pInSys, enumFlags);
				break;
			case t_objectType::t_SemiSymmetricGraph:
				enumFlags |= t_outColumnOrbits + t_outStabilizerOrbit + t_colOrbitsConstructed + t_alwaysKeepRowPermute;
				pInSys = new Class2(CSemiSymmetricGraph)(pParam->v, k, pParam->r, lambda);
				pInSysEnum = new Class2(CIG_Enumerator)(pInSys, pParam, enumFlags, FirstPath);
				break;
			case t_objectType::t_Kirkman_Triple: 
				kirkmanTriples = true;
			case t_objectType::t_CombinedBIBD:
				pInSys = new Class2(CCombinedBIBD)(pParam->v, k, kirkmanTriples, lambda);
				pInSysEnum = new Class2(CCombBIBD_Enumerator)(pInSys, enumFlags + t_useGroupOnParts);
				break;
			default:
				printf("LambdaSet has %zu elements, in that case the object type cannot be \'%s\'\n", lambda.size(), obj_name[+objType]);
				printf("Possible values are:");
				for (int i = 0; i < countof(obj_types); i++)
					printf("\n    %s", obj_name[+obj_types[i]]);

				return false;
			}
		}
		else {
			const auto λ = lambda[0];
			pInSys = new Class2(C_BIBD)(pParam->v, k, 2, λ);
			const auto r = pInSys->GetNumSet(t_rSet)->GetAt(0);
			const int maxNumbCommonElement = (2 * k * λ / r) + (r - k - λ);
			if (maxNumbCommonElement < k) {
				// According the theorem by Connor, maxNumbCommenElement is
				// the maximal number of common elements for two blocks
				enumFlags |= t_noReplicatedBlocks;
				pInSys->setMaxBlockIntrsection(static_cast<T>(maxNumbCommonElement));
			}

			if (λ > 1) {
				if (pParam->enumFlags & t_use_3_condition)
					enumFlags |= t_use_3_condition;

				if ((pParam->enumFlags & t_symmetrical_t_cond) && k == r)
					enumFlags |= t_symmetrical_t_cond;
			}

			if (objType != t_objectType::t_TripleSystem) {
				pInSysEnum = new Class2(CBIBD_Enumerator)(pInSys, enumFlags);
				objType = t_objectType::t_BIBD;
			}
		}
	}
	else {
		pInSys = new Class2(C_tDesign)(pParam->t, pParam->v, k, lambda[0]);
		pInSysEnum = new Class2(C_tDesignEnumerator)(static_cast<TDesignPntr>(pInSys), enumFlags);
		objType = t_objectType::t_tDesign;
	}

	pInSys->setObjectType(objType);

	char buff[256] = {}, buffer[256] = {};
	if (pInSysEnum)
		MAKE_JOB_TITLE(pInSysEnum, pParam, buff, countof(buff));
	else {
		const auto name = obj_name[static_cast<int>(t_objectType::t_TripleSystem)];
		sprintf_s(buff, "%s(%d)", name, pParam->v);
	}
	cout << buff;

	auto* pEnumInfo = new CInsSysEnumInfo<TDATA_TYPES>(buff);
	const bool resetMTlevel = pParam->MT_level() == 0;
	if (pInSysEnum) {
		pEnumInfo->setDesignInfo(pParam);
		if (FirstPath) {
			FOPEN(outFile, pSummFile, "w");
			pEnumInfo->outEnumInformation(&outFile, pInSysEnum->enumFlags(), false);
			outString("         BIBDs:                     Canonical:      NRB #:      Constructed:    Run Time (sec):\n", pSummFile);
		}

		if (resetMTlevel) {
			// The row number, on which the threads will be launched was not defined.
			// Let's do it here by other parameters
			pParam->set_MT_level(pInSysEnum->define_MT_level(pParam));
		}
	}

	try {
		bool retVal = false;
		if (objType != t_objectType::t_TripleSystem) {
			retVal = pInSysEnum->Enumerate(pParam, PRINT_TO_FILE, pEnumInfo);
			if (retVal) {
				pEnumInfo->reportResult(buffer, countof(buffer));
				outString(buffer, pSummFile);
				cout << '\r' << buffer;
			}
		}
		else {
			alldata sys(pParam->v);
			retVal = sys.Run();
		}

		if (!retVal) {
			cout << "\rSome problem was found during the enumeration\n";
			pInSysEnum->closeFile();
		}
	}
	catch (...) {
		pInSysEnum->closeFile();
	}


	delete pInSys;
	if (pParam->find_all_2_decomp) {
		auto* pDesignDB = pInSysEnum->designDB();
		if (pParam->find_master_design) {
			// Compare two databases with BIBDs and make output
			// of the designs which are NOT in second database.
			auto* pComplementDB = new CDesignDB(pDesignDB->recordLength());
			pComplementDB->combineDesignDBs(pParam->designDB(), pDesignDB, true);
			delete pDesignDB;
			delete pEnumInfo;
			output_2_decompInfo<TDATA_TYPES>(pParam, pComplementDB, *outInfo, true);
		}
		else {
			// Saving pEnumInfo for use it after enumeration of combined BIBDs
			pParam->setEnumInfo(pEnumInfo);
			pParam->setDesignDB(pDesignDB);
			// Saving mt_level used for regular BIBDs enumeration
			pParam->set_MT_level(pParam->MT_level(), 1);
			char buffer[265];
			sprintf_s(buffer, "Number of %s's which are NOT combined for the following {lambda_1, lambda_2}:\n  ", buff);
			*outInfo = buffer;
		}

		pInSysEnum->setDesignDB(NULL);
	}
	else {
		delete pEnumInfo;
	}

	delete pInSysEnum;

	if (resetMTlevel)
		pParam->set_MT_level(0);

	CloseCanonInfo();
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

bool designParam::LaunchEnumeration(t_objectType objType, int find_master, int find_all_2_decomp, int use_master_sol, bool& firstRun)
{
	uint iMax = 0;
	uint lambda = 1;
	find_master_design = 0;   // This option for CombBIBDs only
	auto lambdaSet = InterStruct()->lambdaPtr();
	uint baseLambda = 0;
	switch (objType) {
	case t_objectType::t_Kirkman_Triple:
	case t_objectType::t_TripleSystem:
		break;
	default: baseLambda = this->lambda()[0];
	}

	switch (objType) {
	case t_objectType::t_TripleSystem:
	case t_objectType::t_Kirkman_Triple:
		if (this->v % 6 != 3) {
			printf("\nNumber of elements for Kirkman triple system should be equal to 3(mod 6)");
			printf("\nNow it is set to %d", this->v);
			return false;
		}
 
		this->k = 3;
		if (objType == t_objectType::t_TripleSystem) {
			lambdaSet->push_back(1);
			break;
		}
	
		lambdaSet->push_back(0);
		lambdaSet->push_back(r = 1);
	case t_objectType::t_CombinedBIBD:
		find_master_design = find_master;
		this->use_master_sol = 0;
		break;
	case t_objectType::t_BIBD: {
			const auto vMinus1 = this->v - 1;
			const auto kMinus1 = this->k - 1;
			const auto r = vMinus1 * baseLambda / kMinus1;
			if (r * kMinus1 != vMinus1 * baseLambda) {
				printf("\nParameters (v, k, lambda) = (%d, %d, %d) cannot be parameters of BIBD.", this->v, this->k, baseLambda);
				return false;
			}
			this->find_all_2_decomp = find_all_2_decomp ? 2 : 0;
			if (find_all_2_decomp) {
				while (lambda <= (baseLambda >> 1) && (vMinus1 * lambda / kMinus1) * kMinus1 != vMinus1 * lambda)
					lambda++;

				if (lambda > (baseLambda >> 1)) {
					printf("\n BIBD is not a Combined BIBD.");
					return false;
				}

				setLambdaStep(lambda);
				iMax = baseLambda / (2 * lambda);
			}
		}
	}

	const string summaryFile = this->workingDir + pSummaryFile;
	const char* pSummFile = summaryFile.c_str();
	std::string outInfo;
	for (uint i = 0; i <= iMax; i++) {
		if (i) {
			// For all 2-part decomposition search only
			this->objType = t_objectType::t_CombinedBIBD;
			find_all_2_decomp = this->find_master_design = 1;
			this->use_master_sol = 0;
			lambdaSet->resize(0);
			lambdaSet->push_back(i * lambda);
			lambdaSet->push_back(baseLambda - i * lambda);
		}
		if (!RunOperation<TDATA_TYPES>(this, pSummFile, firstRun, &outInfo))
			break;

		firstRun = false;
	}

	if (this->find_all_2_decomp) {
		output_2_decompInfo<TDATA_TYPES>(this, designDB(1), outInfo, false, pSummFile);
		delete enumInfo();
		setEnumInfo(NULL);
	}

	for (int i = 0; i < 2; i++) {
		delete designDB(i);
		setDesignDB(NULL, i);
	}

	this->use_master_sol = use_master_sol;
	find_master_design = find_master;
	logFile = "";
	setLambdaStep(0);
	setEmptyLines();
	return true;
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
			line = line.substr(pos + 1);
			std::replace(line.begin(), line.end(), '\\', '/');
			if (line.back() != '/')
				line += '/';


			struct stat sb;
			const auto* pWorkDir = line.c_str();
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

			if (!pCause) {
				newWorkDir = line;
				continue;
			}

			printf("Cannot %s working directory: \'%s\'\n", pCause, pWorkDir);
			break;
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
		getBooleanParam(line, "USE_SYMMETRICAL_T-CONDITION", t_symmetrical_t_cond, &param->enumFlags);
		getBooleanParam(line, "USE_3-ELEMENT_CONDITION", t_use_3_condition, &param->enumFlags);
		getBooleanParam(line, "UPDATE_RESULTS", t_update_results, &param->enumFlags);

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
				break;
			}
		}

		// Define output type
		pos = find(line, "OUTPUT");
		if (pos != string::npos) {
			outType = t_Summary;
			string output = line.substr(pos);
			bool val;
			if (getIntegerParam(output, "ONLY_SIMPLE_DESSIGN", &val))
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

		if ((param->objType = objType) == t_objectType::t_SemiSymmetricGraph) {
			int InconsistentGraphs(designParam * pParam, const char* pSummaryFileName, bool firstPath);
			InconsistentGraphs(param, pSummaryFile, firstRun);
		}
		else {
			if (!param->LaunchEnumeration(objType, find_master_design, find_all_2_decomp, use_master_sol, firstRun))
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
