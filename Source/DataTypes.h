#pragma once
#include <vector>

#define TEST					  false //true   // Test mode

#ifdef USE_CUDA
#define CONSTR_ON_GPU				0						// 1 - Start using GPU for object construction
#define USE_THREADS					1						// Should be at least 1
#define CANON_ON_GPU                (1 && CONSTR_ON_GPU==0)	// 1 - Start using GPU for canonicity testing
#define NUM_GPU_WORKERS             128
//#include "host_defines.h"
#include "cuda_runtime_api.h"
#define CC __host__ __device__		// CUDA_CALLABLE
#if CONSTR_ON_GPU
#define CK	CC __noinline__
#else
#define CK
#endif
#else
#define USE_THREADS					15 //(TEST? 0 : 15)	// default numer if thread used, If NOT 0, it can be chanaged by THREAD_NUMBER
#define CONSTR_ON_GPU				0
#define CANON_ON_GPU				0
#define NUM_GPU_WORKERS				0
#define CC
#define CK
#endif

#define USE_COL_PERMUT				1
#define WAIT_THREADS				(USE_THREADS && 0)

#if CONSTR_ON_GPU
#define USE_THREADS_ENUM			0 
#define EXIT(x)
#define setOutFile(x)
#define REPORT_PROGRESS(x,...)
#else
#define USE_THREADS_ENUM			USE_THREADS
#if !USE_CUDAFthreads
#define EXIT(x)						
#else
#define EXIT(x)						exit(x)
#endif
#define REPORT_PROGRESS(x, ...)		x->reportProgress(__VA_ARGS__)
#endif

#define MAC							1
#if WIN && MAC == 0
    #define _AFXDLL
    #define _WIN32_WINNT _WIN32_WINNT_MAXVER
    #include <afx.h>
    #include <afxtempl.h>
#else
    #include "ClassArray.h"
    #define CArray  CClassArray
#endif
#include <iostream>

#define OUT_RESULT					1 && (CONSTR_ON_GPU == 0 && (CANON_ON_GPU == 0 || CANON_ON_GPU))

#ifdef WIN


#ifndef USE_CUDA
#ifdef _DEBUG  
	// For Memory leakage detection
	#define CHECK_MEMORY_LEAK		0
#if CHECK_MEMORY_LEAK
	#define _CRTDBG_MAP_ALLOC
	#include <stdlib.h>
	#include <crtdbg.h>
	#define DEBUG_CLIENTBLOCK   new( _CLIENT_BLOCK, __FILE__, __LINE__)
	#define new DEBUG_CLIENTBLOCK
#else
	#define _CrtSetReportMode(x,y)
	#define _CrtSetReportFile(x,y)
	#define _CrtDumpMemoryLeaks()
#endif

#else  
#define DEBUG_CLIENTBLOCK  
#endif // _DEBUG  
#endif

	// Assembly is used for Windows only
    #ifndef X64
        #define USE_ASM				0	// 0 or 1  // 09/07/2014 1 is not working since CCanonicityChecker::checkColOrbit() needs to be updated
    #else
		#if !USE_THREADS				// Assembly code IS NOT thread safe for 64-bit mode
			#define USE_ASM			0	// 0 or 2  // 09/07/2014 1 is not working since CCanonicityChecker::checkColOrbit() needs to be updated
		#else
			#define USE_ASM			0
		#endif
    #endif

#if USE_THREADS
	#define USE_BOOST				0	// Fof OS the usage of boost is not implemented
	#define USE_POOL				(USE_BOOST && 0)
#else
	#define USE_BOOST				0	// For OS the usage of boost is not implemented
#endif

	#define OS "Windows"
#else
    #define OS  "Mac OS"
	#define _CrtSetReportMode(x,y)
	#define _CrtSetReportFile(x,y)
	#define _CrtDumpMemoryLeaks()
  #define USE_ASM				  0
#endif

#define USE_EXRA_EQUATIONS			0   // Construct and use additional equations (for instance, for t-designs)

#define PRINT_RESULTS				(1 && USE_EXRA_EQUATIONS)
#define USE_CANON_GROUP				1	// Reshafle the solutions for (n+1)-th row of matrix M
										// with respect to the Aut G(M(n)) of the first n rows of M
#define USE_STRONG_CANONICITY	(USE_CANON_GROUP && 0)			// Find equivalent solutions for ALL USED canonical solutions of M(n+1)
#define USE_STRONG_CANONICITY_A	(USE_STRONG_CANONICITY && 1)	// Strong canonicity on ALL (not only used) canonical solutions of M(n+1)

#define SLIP_TIME					500  //    in microseconds
#define REPORT_INTERVAL				(5 * 1000 * 100)
#define REPORT_INTERVAL_OBJ_NUMB	10000

#define SOLUTION_STATISTICS			0
#define WRITE_MULTITHREAD_LOG		0

#if TEST
#define PRINT_TO_FILE				0
#define PRINT_SOLUTIONS				1
#define PRINT_CURRENT_MATRIX		1
#define PRINT_PERMUTATION			0 // Output of the current permutation of column's orbit during the canonicity check
#define OUT_PERMUTATION             0 // Output of permutations generated on
									  //   a) matrix rows: when canonicity is checked on the totally constructed matrix:        1
									  //   b) orbits of columns, during the construction of the group for sorting of solutions: 2
									  //   c) (a)  + (b):                                                                       3

#else
#define PRINT_TO_FILE				1	// Write files with the results for each set of parameters
#define PRINT_SOLUTIONS				0
#define PRINT_CURRENT_MATRIX		0
#endif

#define PRINT_CANON_GROUP			0
#define CHECK_CCC					0 // 151315
#define CHECK_CONSTRUCTED			0


#define USE_MUTEX     (USE_THREADS >= 1 && (PRINT_SOLUTIONS || PRINT_PERMUTATION || PRINT_PERMUTATION))
#if USE_MUTEX
#include <mutex>
extern std::mutex out_mutex;
#define MUTEX_LOCK(x)		x.lock()
#define MUTEX_UNLOCK(x)		x.unlock();
#else
#define MUTEX_LOCK(x)
#define MUTEX_UNLOCK(x)
#endif

#if USE_THREADS
	#if USE_BOOST
		#include <boost/thread/thread.hpp>
		#define THREAD_NAME_SPACE boost
	#else
		#define THREAD_NAME_SPACE std
	#endif
#endif

#define FILE_NAME(x)			x".txt"
#define CURRENT_RESULTS			"_dbl"
#define INTERMEDIATE_RESULTS	"_tmp"

#ifdef WIN
	#define OPEN_FILE(x, y, z)	 fopen_s(&x, y, z)
    #define FOPEN(x, y, z)	  	 FILE *x; OPEN_FILE(x, y, z)
#else
	#define sprintf_s(x, y, ...) sprintf(x, __VA_ARGS__)
	#define strcpy_s(x, y, z)    strcpy(x, z)
	#define memcpy_s(x, y, ...)  memcpy(x, __VA_ARGS__)
    #define OPEN_FILE(x, y, z)	 x = fopen(y, z)
	#define FOPEN(x, y, z)	  	 FILE *OPEN_FILE(x, y, z)
#endif

#define FCLOSE(file)			do {if (file) fclose(file);} while(0)

#define SNPRINTF(x, len, ...)	static_cast<size_t>(snprintf(x, len, __VA_ARGS__))
#define SPRINTF(x, ...)			SNPRINTF(x, sizeof(x), __VA_ARGS__)
#define BEG_OUT_BLOCK			"<<<< "		// Marks for beginning and end of the output info, which 
#define END_OUT_BLOCK			">>>> "		// will be skipped during the comparison of output results

#ifndef USE_CUDA
#define THROW() throw;
#else
#define THROW()
#endif


typedef unsigned int		uint;
typedef unsigned short		ushort;
typedef unsigned char		uchar;
typedef unsigned long long	ulonglong;


#define MATRIX_ELEMENT_TYPE  	unsigned __int8 //uchar
#define SIZE_TYPE				unsigned __int8 //uchar //uint16_t //uchar //uint16_t
#define TDATA_TYPES				SIZE_TYPE, MATRIX_ELEMENT_TYPE 
#define ELEMENT_MAX				static_cast<SIZE_TYPE>(-1)

#define Class1(x)               x<S>
#define Class1Def(x)            template<typename S> class x
#define FClass1(x, ...)			template<typename S> __VA_ARGS__  Class1(x)

#define TFunc2(x, ...)          template<typename T, typename S> __VA_ARGS__ x
#define Class2(x)               x<T,S>
#define Class2Def(x)            TFunc2(x, class)
#define FClass2(x, ...)			TFunc2(Class2(x), __VA_ARGS__)

#define MatrixData(...)			FClass2(CMatrixData, __VA_ARGS__)
#define PermutStorage(...)		FClass2(CPermutStorage, __VA_ARGS__)
#define TDesign(...)			FClass2(C_tDesign, __VA_ARGS__)
#define CombinedBIBD(...)		FClass2(CCombinedBIBD, __VA_ARGS__)
#define CanonicityChecker(...)	FClass2(CCanonicityChecker, __VA_ARGS__)


#define EnumInfo				Class2(CEnumInfo)
#define Enumerator				Class2(CEnumerator)
#define GPU_CheckerInfo         Class2(CGPU_CheckerInfo)
#define GPU_CanonChecker		Class2(CGPU_CanonChecker)
#define MatrixCanonCheckerGPU	Class2(CMatrixCanonCheckerGPU)


#define MatrixDataPntr			Class2(CMatrixData) *
#define MatrixPntr				Class2(CMatrix) *
#define MatrixColPntr			Class2(CMatrixCol) *
#define InSysPntr               Class2(C_InSys) *
#define TDesignPntr				Class2(C_tDesign) *
#define RowSolutionPntr			Class2(CRowSolution) *
#define CanonicityCheckerPntr	Class2(CCanonicityChecker) *
#define InSysSolverPntr         Class2(CInSysSolver) *
#define EnumInfoPntr			Class2(CEnumInfo) *
#define EnumeratorPntr			Class2(CEnumerator) *
#define ThreadEnumeratorPntr    Class2(CThreadEnumerator) *
#define PermutStoragePntr		Class2(CPermutStorage) *

#define ColOrbPntr				Class1(CColOrbit) *
#define SimpleArrayPntr			Class1(CSimpleArray) *
#define VectorPntr				Class1(CVector) *

#define _FRMT				"u"
#define ME_FRMT				"%" _FRMT


#define VECTOR_ELEMENT_TYPE  	SIZE_TYPE
typedef CArray<VECTOR_ELEMENT_TYPE, VECTOR_ELEMENT_TYPE> CArrayOfVectorElements;

#define PERMUT_ELEMENT_TYPE  	size_t
typedef CArray<PERMUT_ELEMENT_TYPE, PERMUT_ELEMENT_TYPE> CArraySolutionPerm;
#define PERMUT_ELEMENT_MAX		UINT64_MAX

Class1Def(CSimpleArray) {
public:
	CC inline CSimpleArray(size_t len = 0) : m_nLen(len)	{ m_pElem = len? new S[len] : NULL; }
	CC virtual ~CSimpleArray()								{ delete [] elementPntr(); }
	inline void Init(size_t len, S *pElem)					{ m_pElem = pElem; m_nLen = len; }
	CC inline auto element(size_t idx) const				{ return m_pElem[idx]; }
	inline void setElement(size_t idx, S val)				{ m_pElem[idx] = val; }
	CC inline auto *elementPntr() const						{ return m_pElem; }
	CC inline auto numElement() const						{ return m_nLen; }
	inline auto GetAt(size_t idx)  const					{ return element(idx); }
	inline auto *GetElement(size_t idx)  const				{ return elementPntr() + idx; }
protected:
private:
    S *m_pElem;
    size_t m_nLen;
};

Class1Def(BlockGroupDescr) : public CSimpleArray<S> {
public:
	CC BlockGroupDescr(size_t nParts) : CSimpleArray<S>(nParts << 1) {}
	CC inline auto numParts() const							{ return static_cast<S>(numElement() >> 1); }
	CC inline auto getShift(S idx) const					{ return static_cast<S>(element(idx << 1)); }
	CC inline auto colNumb(S idx = 0) const					{ return static_cast<S>(element((idx << 1) + 1)); }
	CC inline void SetPartInfo(size_t idx, S shift, S len)	{
		setElement(idx << 1, shift);
		setElement((idx << 1) + 1, len);
	}
	CC auto GetPartInfo(S idx, S *pLen) const {
		*pLen = colNumb(idx);
		return getShift(idx);
	}
	CC inline void CopyPartInfo(const BlockGroupDescr<S> *pPartInfo) {
		for (size_t idx = 0; idx < pPartInfo->numElement(); idx++)
			setElement(idx, pPartInfo->GetAt(idx)); 
	}
};

Class1Def(CCounter) : public Class1(CSimpleArray) {
public:
    CC CCounter(size_t len) : Class1(CSimpleArray)(len)  { resetArray(); }
	CC ~CCounter()									{}
	CC inline void resetArray()						{ memset(this->elementPntr(), 0, this->numElement() * sizeof(S)); }
	CC inline void incCounter(int idx)              { ++*(this->elementPntr() + idx); }
private:
};

Class1Def(CContainer) : public Class1(CSimpleArray) {
public:
	CC CContainer(size_t len) : Class1(CSimpleArray)(len){ resetArray(); }
	CC ~CContainer()								{}
    CC inline void resetArray()                     { m_nNumb = 0; }
	CC inline void addElement(S val)                { *(this->elementPntr() + m_nNumb++) = val; }
    CC inline auto numb() const						{ return m_nNumb; }
    inline void incNumb()                           { m_nNumb++; }
	inline auto GetSize() const						{ return numb(); }
private:
	size_t m_nNumb;
};

Class1Def(CMapping) : public Class1(CSimpleArray) {
public:
    CK CMapping(size_t len) : Class1(CSimpleArray)(len<<1)  {}
	CK CMapping(S to, S from, size_t len = 1) : Class1(CSimpleArray)(len << 1)
													{ resetMapping(); addMapping(to, from); }
	CC ~CMapping()									{}
    CK void addMapping(S to, S from);
	void removeMapping(S to);
	void removeLastMapping(size_t n = 1)			{ m_nMapPos -= n << 1; }
	void restoreLastMapping(size_t n = 1)			{ m_nMapPos += n << 1; }
	CK inline void resetMapping()                   { m_nMapPos = 0; }
	inline uint nElement() const					{ return getMapPosition() >> 1; }
	CK inline const S *getMappingPntr() const		{ return this->elementPntr(); }
	CK inline const S *getLastMapping() const		{ return getMappingPntr() + getMapPosition(); }
	const S *findMapping(S to, const S *pTo = NULL, const S *pToLast = NULL) const;
	bool isEmpty() const                            { return !m_nMapPos; }
	inline void adjustElement(int idx, S val = 1)   { *(this->elementPntr() + idx) -= val; }
protected:
	CK inline auto getMapPosition() const           { return m_nMapPos; }
    
    uint m_nMapPos;
};

FClass1(CMapping, void)::addMapping(S to, S from)
{
	assert(m_nMapPos < (this->numElement() << 1));
	auto *pTo = (S *)getLastMapping();
    *pTo = to;
    *(pTo+1) = from;
    m_nMapPos += 2;
}

FClass1(CMapping, const S*)::findMapping(S to, const S *pTo, const S *pToLast) const
{
	// Searching for element to be removed
	if (!pTo)
		pTo = getLastMapping();

	if (!pToLast)
		pToLast = getMappingPntr();

	while ((pTo -= 2) >= pToLast && *pTo != to);

	return pTo >= pToLast ? pTo : NULL;
}

FClass1(CMapping, void)::removeMapping(S to)
{
	// Searching for element to be removed
	auto *pTo = findMapping(to);
	if (!pTo)
		return;

	// Element found
	const int len = getLastMapping() - pTo - 2;
	if (len)
		memcpy(pTo, pTo + 2, len * sizeof(*pTo));

	removeLastMapping();
}

#define countof(x)  (sizeof(x)/sizeof(x[0]))

extern int ccc;
#if CHECK_CCC
	#define CHECK_VAL(ccc)		(ccc > CHECK_CCC)
#else
	#define CHECK_VAL(ccc)		true
#endif

#if CHECK_CONSTRUCTED
	#define START_NUMBER         9
	#define END_NUMBER			91
	#define CHECK_CONSTR(x, y)    enumInfo()->numMatrOfType(t_design_type::t_canonical) >= x && enumInfo()->numMatrOfType(t_design_type::t_canonical) <= y
#else
#define CHECK_CONSTR(x, y)		true
#endif

#define MAKE_OUTPUT()			CHECK_VAL(ccc) && CHECK_CONSTR(START_NUMBER, END_NUMBER)

#if PRINT_SOLUTIONS
	#define OUTPUT_SOLUTION(x, file, nRow, f, nPartFrom, nPartTo)	if (MAKE_OUTPUT()) \
																		{ this->printSolutions(x, file, nRow, f, nPartFrom, nPartTo); }
#else
    #define OUTPUT_SOLUTION(x,...)
#endif

#define OUT_MATRIX(x, y, z, w, v, canon)		{ MUTEX_LOCK(out_mutex); x->printOut(y, z, w, NULL, v->numMatrOfType(t_design_type::t_canonical)+1, canon);  MUTEX_UNLOCK(out_mutex); }
#if PRINT_CURRENT_MATRIX
	#define OUTPUT_MATRIX(x, y, z, v, canon)	if (MAKE_OUTPUT()) \
													OUT_MATRIX(x, y, z, ++ccc, v, canon)
#else
    #define OUTPUT_MATRIX(x,...)
#endif

#if PRINT_PERMUTATION						
	#define OUTPUT_PERMUTATION(x, f, n, p)		MUTEX_LOCK(out_mutex); x->outputPerm(f, n, p);  MUTEX_LOCK(out_mutex); 
#else
    #define OUTPUT_PERMUTATION(x,...)
#endif

#if PRINT_CANON_GROUP
	#define OUTPUT_CANON_GROUP(cond, x, f)		if (cond) x->outputPermutations(f, x->lenPerm()/*numColOrb()*/)
#else
	#define OUTPUT_CANON_GROUP(cond,...)
#endif

typedef enum {
	t_threadUndefined,
	t_threadLaunchFail,
	t_threadRunning,
	t_threadFinished,
	t_threadNotUsed
} t_threadCode;

#if USE_THREADS && WRITE_MULTITHREAD_LOG
void thread_message(int threadIdx, const char *pComment, t_threadCode code, void *pntr = 0);
#else
	#define thread_message(threadIdx, pComment, code, ...)
#endif

size_t outString(const char *str, FILE *file);
size_t outString(const char *str, const char *fileName, const char *mode = "a");

#if PRINT_RESULTS
#define PR_EQU_NUMB				20
#define PRINT_RES_ROW_NUMB		 4
#define PRINT_RES_SOL_NUMB		 1
#define EQUATION_ID_TO_PRINT    -1
#define PRINT_RES_EQU_NUMB		PRINT_RES_ROW_NUMB - 1   // Number of equations to be printed 
														 // for 3-design it should be PRINT_RES_ROW_NUMB - 1,
														 // if you want to print all equestions

#define MY_ID			public:	int myID;
#define MY_ORB_ID		public: int myOrbID;

class CVariableMapping;
class CEquation;
extern CEquation *pMyEquA[PR_EQU_NUMB];
extern int ppp, ggg;
extern size_t printResRowNumb;
extern size_t printResNumVar;

FILE *openOutFile();
void assignPrintResultInfo(int &myID, CEquation *pntr);
void assignOrbID(int &myOrbID);
void printResults(FILE *file, VECTOR_ELEMENT_TYPE *pResult, size_t len, int varIdx);
void printResults(VECTOR_ELEMENT_TYPE *pResult, int lambdaToSplit, int idx, int varIdx);
void printAddEquationVariable(const CEquation *pEqu, VECTOR_ELEMENT_TYPE varIdx, VECTOR_ELEMENT_TYPE valVar, bool addFlg = true);
void printDefinedVariable(VECTOR_ELEMENT_TYPE varID, VECTOR_ELEMENT_TYPE varValue);
void printFinalResultForRightPart(CVariableMapping *pVarValue);
void setPrintResultRowNumber(size_t nRow);
void setPrintResultNumVar(size_t numVar);
#else
// empty macroses
#define MY_ID
#define MY_ORB_ID

#define assignPrintResultInfo(x, y)
#define assignOrbID(myOrbID)
#define printResults(pResult, lambdaToSplit, idx, varIdx)
#define printAddEquationVariable(pEqu, varIdx, valVar, flg)
#define printMapInfo(pVarValue)
#define printDefinedVariable(varID, varValue)
#define printFinalResultForRightPart(pVarValue)
#define setPrintResultRowNumber(nRow)
#define setPrintResultNumVar(nVar)
#endif

enum class t_objectType {
	t_BIBD,			// default
	t_CombinedBIBD,
	t_tDesign,
	t_PBIBD,
	t_IncidenceSystem,
	t_SemiSymmetricGraph
};

template<typename T>
inline int operator + (T val) { return static_cast<int>(val); }

template<typename T>
inline T operator ++ (T& val) { return val = static_cast<T>(+val + 1); }

typedef enum {
	t_enumDefault			= 0,
	t_IS_enumerator			= 1 << 0,
	t_matrixOwner			= 1 << 1,
	t_noReplicatedBlocks	= 1 << 2,
	t_outColumnOrbits		= 1 << 3,
	t_outStabilizerOrbit	= 1 << 4,	// Keep the orbits of stabilizer of first elements
	t_colOrbitsConstructed  = 1 << 5,
	t_printTransposedMatrix = 1 << 6,
	t_alwaysKeepRowPermute	= 1 << 7,   // Keep generating set of the Aut(M) acting on the rows of partially constructed matrix M
	t_outRowPermute			= 1 << 8,   // Keep generating set of the Aut(M) acting on the rows of fully constructed matrix M
	t_outRowOrbits			= 1 << 9,
	t_useGroupOnParts       = 1 << 10,
	t_allFlags				= -1
} t_EnumeratorFlags;

class CInterStruct {
public:
	inline CInterStruct(int mult = 1)		{ setMult(mult); }
	inline auto* Counterparts() const		{ return m_pCounterparts; }
	inline ~CInterStruct()					{ delete Counterparts(); }
	inline const auto &lambda() const		{ return iParam[0]; }
	inline const auto &lambdaA() const		{ return iParam[1]; }
	inline const auto &lambdaB() const		{ return iParam[2]; }
	inline auto *lambdaPtr()				{ return iParam; }
	inline auto *lambdaAPtr()				{ return iParam + 1; }
	inline auto *lambdaBPtr()				{ return iParam + 2; }

	inline bool isValid() const				{ return Counterparts(); }
	inline void InitCounterparts()			{ m_pCounterparts = new std::vector<CInterStruct *>(); }
	inline void setNext(CInterStruct *pntr)	{ m_pNext = pntr; }
	inline auto *getNext() const			{ return m_pNext; }
	inline void setMult(int val)			{ m_mult = val; }
	inline int mult() const					{ return m_mult; }

private:
	std::vector<uint> iParam[3];
	std::vector<CInterStruct *> *m_pCounterparts = NULL;
	int m_mult;
	CInterStruct *m_pNext = NULL;
};


class designParam {
public:
	designParam()									{ m_pInterStruct = new CInterStruct(); }
	~designParam() {
		CInterStruct *pntr, *pTmp = InterStruct();
		while (pntr = pTmp) {
			pTmp = pTmp->getNext();
			delete pntr;
		}
	}
	inline CInterStruct *InterStruct()	const		{ return m_pInterStruct; }
	inline void SetInterStruct(CInterStruct *pntr)	{ m_pInterStruct = pntr; }
	t_objectType objType = t_objectType::t_BIBD;
	int v = 0;
	int k = 0;
	int r = 0;
	int t = 0;
	int mt_level = 0;				// Matrix row number, where the threads will be launched
	uint outType = 0;				// Flags which define the output information of the task
	uint grpOrder = 0;				// Limits for order of the group of the matrices which will be printed
	size_t threadNumb = 0;			// Number of threads launched to perform task
	int use_master_sol = 0;			// the solutions defined by master will be used (1) or copied (0) by the threads
	int find_master_design = 0;     // find master BIBD, when CombBIBDs are enumerated
	int format_master_BIBDs = 0;	// Format of output of information regarding master BIBD decompositions: 0-2
	int thread_master_DB = 0;       // Keep Master DB for each thread (no mutex will be used)
	bool firstMatr = true;			// TRUE, when first matrix of the set was not yet outputted
	bool noReplicatedBlocks = true;	// TRUE, when only block designs with no replicated blocks should be constructed
	std::string workingDir = "";	// Current working directory name
	std::string logFile = "";		//
	size_t rewindLen = 0;			// Length of the portion of log file, which probably will be rewinded
	const auto &lambda() const					{ return m_pInterStruct->lambda(); }
	const auto &lambdaA() const					{ return m_pInterStruct->lambdaA(); }
	const auto &lambdaB() const					{ return m_pInterStruct->lambdaB(); }
	inline auto lambdaSizeMax() const			{ return m_lambdaSizeMax; }
	inline void setLambdaSizeMax(size_t val)	{ m_lambdaSizeMax = val; }
private:
	CInterStruct *m_pInterStruct = NULL;
	size_t m_lambdaSizeMax = 0;		// Maximal number of elements in lambda()
									// (will be used for formated output)
};

#define VAR_1		1
