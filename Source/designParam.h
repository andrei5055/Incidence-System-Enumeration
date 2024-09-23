#pragma once
#include "DataTypes.h"

class CDesignDB;
Class2Def(CInsSysEnumInfo);

class publicParam {
public:
	t_objectType objType = t_objectType::t_BIBD;
	int v = 0;
	int b = 0;
	int k = 0;
	int r = 0;
	int t = 0;
	unsigned char matrixRank = 2;
	uint outType = 0;				// Flags which defines the output information of the task
	uint grpOrder = 0;				// Limits for order of the group of the matrices which will be printed
	size_t threadNumb = 0;			// Number of threads launched to perform task
	int use_master_sol = 0;			// the solutions defined by master will be used (1) or copied (0) by the threads
	int find_master_design = 0;     // find master BIBD, when CombBIBDs are enumerated
	int thread_master_DB = 0;       // Keep Master DB for each thread (no mutex will be used)
	int find_all_2_decomp = 0;      // Find all decompositions of the BIBDs into 2 components
	bool firstMatr = true;			// TRUE, when first matrix of the set was not yet outputted
	bool noReplicatedBlocks = true;	// TRUE, when only block designs with no replicated blocks should be constructed
	bool m_compress_matrices = false; // Use bitwise compression of the matrices stored in the database
	bool m_bUseThreadPool = false;	// allow threads to start threads
	std::string workingDir = "";	// Current working directory name
	std::string logFile = "";		// Used for semi-symmetric graphs and non-combined BIBDs search
	std::string objSubType = "";
	size_t rewindLen = 0;			// Length of the portion of log file, which probably will be rewinded
	int save_restart_info = 0;		// Save restart information that will be used to restart the program.
	std::string restart_info_dir;
	size_t restart_update_unterval = 10 * 60; // default update interval in sec.
	uint m_enumFlags = 0;
};

class designParam : public publicParam {
public:
	designParam() { m_pInterStruct = new CInterStruct(); }
	~designParam() {
		CInterStruct* pntr, * pTmp = InterStruct();
		while (pntr = pTmp) {
			pTmp = pTmp->getNext();
			delete pntr;
		}
	}
	void setEnumFlags();
	inline auto enumFlags() const					{ return m_enumFlags; }
	inline auto enumFlagsPtr()                      { return &m_enumFlags; }
	inline CInterStruct* InterStruct()	const		{ return m_pInterStruct; }
	inline void SetInterStruct(CInterStruct* pntr)	{ m_pInterStruct = pntr; }
	const auto& lambda() const						{ return m_pInterStruct->lambda(); }
	const auto& lambdaA() const						{ return m_pInterStruct->lambdaA(); }
	const auto& lambdaB() const						{ return m_pInterStruct->lambdaB(); }
	inline auto lambdaSizeMax() const				{ return m_lambdaSizeMax; }
	inline void setLambdaSizeMax(size_t val)		{ m_lambdaSizeMax = val; }
	inline void setDesignDB(const CDesignDB* pntr, int idx = 0) { m_pDesignDB[idx] = pntr; }
	inline const CDesignDB* designDB(int idx = 0) const { return m_pDesignDB[idx]; }
	inline auto* enumInfo() const					{ return m_pEnumInfo; }
	inline void setEnumInfo(CInsSysEnumInfo<TDATA_TYPES>* pntr) { m_pEnumInfo = pntr; }
	inline int MT_level(int idx = 0) const			{ return mt_level[idx]; }
	inline void set_MT_level(int val, int idx = 0)	{ mt_level[idx] = val; }
	inline void setLambdaStep(size_t step)			{ m_lambdaStep = step; }
	inline auto lambdaStep() const					{ return m_lambdaStep; }
	inline auto printEmptyLines() const				{ return m_emptyLines; }
	inline void setEmptyLines(bool val = true)		{ m_emptyLines = val; }
	inline auto printOnlySimpleDesigns() const		{ return m_printSimpleDesign; }
	inline void setPrintOnlySimpleDesigns(bool val = true) { m_printSimpleDesign = val; }
	inline auto compressMatrices() const			{ return m_compress_matrices; }
	inline auto useThreadPool() const				{ return m_bUseThreadPool; }
	inline bool create_commonData() const			{ return threadNumb && MT_level() < v; }
	bool LaunchEnumeration(const char *pSummaryFile, int find_master, int find_all_2_decomp, int use_master_sol, bool& firstRun);
	bool LaunchCanonization();

private:
	CInterStruct* m_pInterStruct = NULL;
	int mt_level[2] = { 0, 0 };		// Matrix row number, where the threads will be launched
	size_t m_lambdaSizeMax = 0;		// Maximal number of elements in lambda() (it will be used for formated output)
	size_t m_lambdaStep = 0;        // Step for parameter lambda, used in the loop for non-combined BIBDs search
	const CDesignDB* m_pDesignDB[2] = { NULL, NULL };
	CInsSysEnumInfo<TDATA_TYPES>* m_pEnumInfo = NULL;
	bool m_emptyLines = true;
	bool m_printSimpleDesign = false;
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