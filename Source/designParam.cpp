#include <fstream>
#include "designParam.h"
#include "EnumInfo.h"
#include "DesignDB.h"
#include "Matrix.h"
#include "IG_Enumerator.h"
#include "CombBIBD_Enumerator.h"
#include "C_tDesignEnumerator.h"
#include "k-SysSupport.h"

const char* obj_name[] = {
	"BIBD",					// t_BIBD,			- default
	"COMBINED_BIBD",		// t_CombinedBIBD,
	"KIRKMAN_TRIPLE_SYSTEM",// t_Kirkman_Triple
	"TRIPLE_SYSTEM",        // t_TripleSystem - construction by Leo's program
	"T-DESIGN",				// t_tDesign,
	"PBIBD",				// t_PBIBD,
	"INCIDENCE",            // t_IncidenceSystem,
	"SEMI_SYMMETRIC_GRAPH", // t_SemiSymmetricGraph
	"CANON_MATR",			// t_CanonMatr
	"SEMI_SYMMETRIC_KSYS"	// t_SemiSym_KSystems
};

template <typename T, typename S>
void output_2_decompInfo(designParam* param, const CDesignDB* pDesignDB, std::string& outputInfo, bool addInfo = false, const char* pSummaryFileName = NULL) {
	const auto& lambda = param->InterStruct()->lambda();
	Class2(C_BIBD) bibd(param->v, param->k, 2, lambda[0] + lambda[1]);
	Class2(CBIBD_Enumerator) bibdEnum(&bibd, t_enumDefault);
	bibdEnum.outNonCombinedDesigns(param, pDesignDB, outputInfo, addInfo);
	if (!addInfo) {
		char buffer[256];
		param->enumInfo()->reportResult(buffer, countof(buffer));
		outString(buffer, pSummaryFileName);
		std::cout << '\r' << buffer;
	}
}

template <typename T, typename S>
void PrepareBIBD_Enumeration(const designParam* pParam, Class2(C_InSys)** ppInSys, Class2(C_InSysEnumerator)** ppInSysEnum, t_objectType &objType, uint& enumFlags)
{
	switch (objType) {
	case t_objectType::t_CanonMatr:
	case t_objectType::t_SemiSym_KSystems:
		*ppInSys = new Class2(C_InSys)(pParam->v, pParam->b, 2, pParam->matrixRank);
		*ppInSysEnum = new Class2(C_InSysCanonizator)(*ppInSys, enumFlags);
		return;
	}

	const auto& lambda = pParam->InterStruct()->lambda();
	const auto λ = lambda[0];
	const auto k = pParam->k;
	auto* pInSys = *ppInSys = new Class2(C_BIBD)(pParam->v, k, 2, λ);
	const auto r = pInSys->GetNumSet(t_rSet)->GetAt(0);
	const int maxNumbCommonElement = (2 * k * λ / r) + (r - k - λ);
	if (maxNumbCommonElement < k) {
		// According the theorem by Connor, maxNumbCommonElement is
		// the maximal number of common elements for two blocks
		enumFlags |= t_noReplicatedBlocks;
		pInSys->setMaxBlockIntrsection(static_cast<T>(maxNumbCommonElement));
	}

	if (λ > 1) {
		if (pParam->enumFlags() & t_use_3_condition)
			enumFlags |= t_use_3_condition;

		if ((pParam->enumFlags() & t_symmetrical_t_cond) && k == r)
			enumFlags |= t_symmetrical_t_cond;
	}

	if (objType != t_objectType::t_TripleSystem) {
		*ppInSysEnum = new Class2(CBIBD_Enumerator)(*ppInSys, enumFlags);
		objType = t_objectType::t_BIBD;
	}
}

template <typename T, typename S>
CInsSysEnumInfo<T, S>* createEnumInfo(designParam* pParam, C_InSysEnumerator<T, S>* pInSysEnum = NULL) {
	char buff[256] = {0};
	if (pInSysEnum)
		MAKE_JOB_TITLE(pInSysEnum, pParam, buff, countof(buff));
	else {
		const auto name = obj_name[static_cast<int>(t_objectType::t_TripleSystem)];
		sprintf_s(buff, "%s(%d)", name, pParam->v);
	}
	std::cout << buff;

	auto* pEnumInfo = new CInsSysEnumInfo<TDATA_TYPES>(buff);
	if (pInSysEnum)
		pEnumInfo->setDesignInfo(pParam);

	return pEnumInfo;
}

template <typename T, typename S>
bool RunOperation(designParam* pParam, const char* pSummFile, bool FirstPath, std::string* outInfo) {
	if (pParam->v <= 0) {
		printf("Problem in RunOperation. Number of elements v = %d <= 0", pParam->v);
		return false;
	}

	InitCanonInfo(pParam->threadNumb);
	C_InSys<TDATA_TYPES>* pInSys = NULL;
	C_InSysEnumerator<TDATA_TYPES>* pInSysEnum = NULL;
	t_objectType objType = pParam->objType;
	uint enumFlags = pParam->enumFlags();

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
			PrepareBIBD_Enumeration(pParam, &pInSys, &pInSysEnum, objType, enumFlags);
		}
	}
	else {
		pInSys = new Class2(C_tDesign)(pParam->t, pParam->v, k, lambda[0]);
		pInSysEnum = new Class2(C_tDesignEnumerator)(static_cast<TDesignPntr>(pInSys), enumFlags);
		objType = t_objectType::t_tDesign;
	}

	pInSys->setObjectType(objType);

	auto* pEnumInfo = createEnumInfo(pParam, pInSysEnum);
	const bool resetMTlevel = pParam->MT_level() == 0;
	if (pInSysEnum) {
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

	char buffer[256] = {};
	try {
		bool retVal = false;
		if (objType != t_objectType::t_TripleSystem) {
			retVal = pInSysEnum->Enumerate(pParam, PRINT_TO_FILE, pEnumInfo);
			if (retVal) {
				pEnumInfo->reportResult(buffer, countof(buffer));
				outString(buffer, pSummFile);
				std::cout << '\r' << buffer;
			}
		}
		else {
			//alldata sys(pParam->v);
			//retVal = sys.Run();
		}

		if (!retVal) {
			std::cout << "\rSome problem was found during the enumeration\n";
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
			sprintf_s(buffer, "Number of %s's which are NOT combined for the following {lambda_1, lambda_2}:\n  ", pEnumInfo->strToScreen());
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

bool designParam::LaunchEnumeration(const char *pSummaryFile, int find_master, int find_all_2_decomp, int use_master_sol, bool& firstRun)
{
	uint iMax = 0;
	uint lambda = 1;
	find_master_design = 0;   // This option for CombBIBDs only
	auto lambdaSet = InterStruct()->lambdaPtr();
	uint baseLambda = 0;
	const auto objType = this->objType;
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

	const std::string summaryFile = this->strParam[t_workingDir] + pSummaryFile;
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
	setLogFile();
	setLambdaStep(0);
	setEmptyLines();
	return true;
}

bool writeTable(const std::string& fn, FILE *file, const char *pComment = NULL) {
	// Copying input file into output
	if (!file)
		return false;

	std::ifstream input(fn);
	if (!input)
		return false;

	if (pComment)
		fprintf(file, "%s\n", pComment);

	for (std::string line; getline(input, line);)
		fprintf(file, "%s\n", line.c_str());

	input.close();
	return true;
}

bool designParam::LaunchCanonization() {
	int retVal = 1;
	C_InSys<TDATA_TYPES>* pInSys = NULL;
	C_InSysEnumerator<TDATA_TYPES>* pInSysEnum = NULL;
	if (strParam[t_objSubType] == "K-SYSTEM") {
		const auto nCols = this->v;
		const auto nRows = (nCols - 1) / (k - 1);
		int reservedElement = 1;
		unsigned char* pSm = new unsigned char[nRows * nCols];

		const auto& inputFile = strParam[t_input_file];
		const auto retVal = readTable(inputFile, nRows, nCols, 1, 0, &pSm, reservedElement, 1000);
		if (retVal) {
			this->b = nRows * v / k;
			// An additional matrix row will be used to store the indices 
			// of the columns corresponding to the original matrix rows.
			this->v++;
			uint enumFlags = this->enumFlags();
			matrixRank = nRows - 1;
			PrepareBIBD_Enumeration(this, &pInSys, &pInSysEnum, objType, enumFlags);

			auto* pEnumInfo = createEnumInfo(this, pInSysEnum);
			pEnumInfo->startClock();
			pInSysEnum->assignDesignParam(this);
			std::string comment("Input data: ");
			comment += "\"" + inputFile + "\"\n";
			const auto file = pInSysEnum->outFile();
			writeTable(inputFile, file, comment.c_str());

			pInSys->prepareFirstMatrixRow(nRows);
			pInSys->convertToBinaryMatrix(pSm, k, nRows);
			pInSysEnum->ConstructCanonicalMatrix(k);

			pEnumInfo->setRunTime();
			pEnumInfo->outRunTimeInfo(file, "\n\nCanonization was done in ");
			pEnumInfo->outRunTimeInfo(file);
			fclose(file);
			delete pEnumInfo;
		}

		delete[] pSm;
	}
	else
		return false;  // there is no converter for such subType.

	delete pInSys;
	delete pInSysEnum;
	return retVal > 0;
}

bool designParam::SemiSymByKSystems() {
	const auto nCols = this->v;
	const auto nRows = (nCols - 1) / (k - 1);
	int lenMatr = nRows * nCols;
	tchar* pSm = new tchar[2 * lenMatr];

	auto pMatr = pSm;
	int numMatr = 2;
	int i = 0;
	for (auto file : { t_input_file, t_extraStrParam }) {
		const auto& inputFile = strParam[file];
		if (!readTable(inputFile, nRows, nCols, 1, i++, &pSm, numMatr, numMatr)) {
			delete[] pSm;
			return false;
		}

		pMatr += lenMatr;
	}

	this->v = this->b = 2 * nRows * v / k;
	uint enumFlags = this->enumFlags();
	matrixRank = 2;
	C_InSys<TDATA_TYPES>* pInSys = NULL;
	C_InSysEnumerator<TDATA_TYPES>* pInSysEnum = NULL;
	PrepareBIBD_Enumeration(this, &pInSys, &pInSysEnum, objType, enumFlags);
	auto* pEnumInfo = createEnumInfo(this, pInSysEnum);
	pEnumInfo->startClock();
	pInSysEnum->assignDesignParam(this);
	const auto file = pInSysEnum->outFile();
	pInSys->convertToSemiSymGraph(pSm, nCols, nRows, k);
	pInSysEnum->ConstructCanonicalMatrix(-1);

	pEnumInfo->setRunTime();
	pEnumInfo->outRunTimeInfo(file, "\n\nCanonization was done in ");
	pEnumInfo->outRunTimeInfo(file);
	fclose(file);

    delete[] pSm;

	delete pInSys;
	delete pInSysEnum;
	delete pEnumInfo;
	return true;
}

void designParam::setEnumFlags() {
	m_enumFlags = noReplicatedBlocks ? t_noReplicatedBlocks : t_enumDefault;
	if ((outType & t_GroupGeneratingSet) == t_GroupGeneratingSet)
		m_enumFlags |= t_outRowOrbits + t_outRowPermute;
	else
		if (outType & t_GroupOrbits)
			m_enumFlags |= t_outRowOrbits;
}

LIBRARY_API const char** designParam::objNames() const {
	return obj_name;
}