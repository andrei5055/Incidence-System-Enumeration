#pragma once
#include <vector>
#include <string>
#include "CudaInfo.h"

#ifdef LIBRARY_EXPORTS
#    define LIBRARY_API __declspec(dllexport)
#else
#    define LIBRARY_API __declspec(dllimport)
#endif

typedef enum {
	t_Summary = 1 << 0,		// default
	t_AllObject = 1 << 1,
	t_GroupOrbits = 1 << 2,
	t_GroupGeneratingSet = t_GroupOrbits + (1 << 3),
	t_Transitive = 1 << 4,
	t_GroupOrderEQ = 1 << 5,
	t_GroupOrderGT = 1 << 6,
	t_GroupOrderLT = 1 << 7,
	t_CombMasters = 1 << 8,
} t_outputType;


enum class t_objectType {
	t_BIBD,			// default
	t_CombinedBIBD,
	t_Kirkman_Triple,
	t_TripleSystem,
	t_tDesign,
	t_PBIBD,
	t_IncidenceSystem,
	t_SemiSymmetricGraph,
	t_CanonMatr
};


typedef enum {
	t_enumDefault			= 0,
	t_IS_enumerator			= 1 << 0,
	t_matrixOwner			= 1 << 1,
	t_noReplicatedBlocks	= 1 << 2,
	t_outColumnOrbits		= 1 << 3,
	t_outStabilizerOrbit	= 1 << 4,	// Keep the orbits of stabilizer of first elements.
	t_colOrbitsConstructed  = 1 << 5,
	t_printTransposedMatrix = 1 << 6,
	t_alwaysKeepRowPermute	= 1 << 7,   // Keep generating set of the Aut(M) acting on the rows of partially constructed matrix M.
	t_outRowPermute			= 1 << 8,   // Keep generating set of the Aut(M) acting on the rows of fully constructed matrix M.
	t_outRowOrbits			= 1 << 9,
	t_useGroupOnParts       = 1 << 10,
	t_symmetrical_t_cond	= 1 << 11,  // Any(λ + 1) elements of symmetric BIBDs(v, v, k, k, λ) simultaneously belong to 0 or 1 block.
	t_use_3_condition		= 1 << 12,  // Use limits on the number of common blocks that any 3 elements belong to. 
	t_update_results        = 1 << 13,  // Updating enumeration results
	t_kSystems              = 1 << 14,  // Enumeration of k-systems
	t_allFlags				= -1
} t_EnumeratorFlags;
#define MATRIX_ELEMENT_TYPE  	unsigned __int8 //uchar
#define SIZE_TYPE				unsigned __int8 //uchar //uint16_t //uchar //uint16_t
#define TDATA_TYPES				SIZE_TYPE, MATRIX_ELEMENT_TYPE 

#define Class1(x)               x<S>
#define Class1Def(x)            template<typename S> class x
#define FClass1(x, ...)			template<typename S> __VA_ARGS__  Class1(x)

#define TFunc2(x, ...)          template<typename T, typename S> __VA_ARGS__ x
#define Class2(x)               x<T,S>
#define Class2Def(x)            TFunc2(x, class)
#define FClass2(x, ...)			TFunc2(Class2(x), __VA_ARGS__)

class CDesignDB;
Class2Def(CInsSysEnumInfo);

class CInterStruct {
public:
	inline CInterStruct(int mult = 1) { setMult(mult); }
	inline auto* Counterparts() const { return m_pCounterparts; }
	inline ~CInterStruct() { delete Counterparts(); }
	inline const auto& lambda() const { return iParam[0]; }
	inline const auto& lambdaA() const { return iParam[1]; }
	inline const auto& lambdaB() const { return iParam[2]; }
	inline auto* lambdaPtr() { return iParam; }
	inline auto* lambdaAPtr() { return iParam + 1; }
	inline auto* lambdaBPtr() { return iParam + 2; }

	inline bool isValid() const { return Counterparts(); }
	inline void InitCounterparts() { m_pCounterparts = new std::vector<CInterStruct*>(); }
	inline void setNext(CInterStruct* pntr) { m_pNext = pntr; }
	inline auto* getNext() const { return m_pNext; }
	inline void setMult(int val) { m_mult = val; }
	inline int mult() const { return m_mult; }

private:
	std::vector<uint> iParam[3];
	std::vector<CInterStruct*>* m_pCounterparts = NULL;
	int m_mult;
	CInterStruct* m_pNext = NULL;
};

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
	LIBRARY_API void setEnumFlags();
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
	LIBRARY_API bool LaunchEnumeration(const char *pSummaryFile, int find_master, int find_all_2_decomp, int use_master_sol, bool& firstRun);
	LIBRARY_API bool LaunchCanonization();
	LIBRARY_API int InconsistentGraphs(const char* pSummaryFileName, bool firstPath);
	LIBRARY_API const char** objNames() const;

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


