#pragma once
#include <vector>

#ifdef USE_CUDA
#define CONSTR_ON_GPU				0						// 1 - Start using GPU for object construction
#define USE_THREADS					1						// Should be at least 1
#define CANON_ON_GPU                (1 && !CONSTR_ON_GPU)	// 1 - Start using GPU for canonicity testing
#define NUM_GPU_WORKERS             128
#include "host_defines.h"
#define CC __host__ __device__		// CUDA_CALLABLE
#if CONSTR_ON_GPU
#define CK	CC __noinline__
#else
#define CK
#endif
#else
#define USE_THREADS					3
#define CONSTR_ON_GPU				0
#define CANON_ON_GPU				0
#define NUM_GPU_WORKERS				0
#define CC
#define CK
#endif

#define WAIT_THREADS				(USE_THREADS && 0)

#if CONSTR_ON_GPU
#define USE_THREADS_ENUM			0 
#define EXIT(x)
#define setOutFile(x)
#define REPORT_PROGRESS(x,...)
#else
#define USE_THREADS_ENUM			USE_THREADS
#if !USE_CUDA
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

#define TEST						0   // Test mode
#define PRINT_RESULTS				(1 && USE_EXRA_EQUATIONS)
#define USE_CANON_GROUP				1	// Reshafle the solutions for (n+1)-th row of matrix M
										// with respect to the Aut G(M(n)) of the first n rows of M
#define USE_STRONG_CANONICITY	(USE_CANON_GROUP && 0)			// Find equivalent solutions for ALL USED canonical solutions of M(n+1)
#define USE_STRONG_CANONICITY_A	(USE_STRONG_CANONICITY && 1)	// Strong canonicity on ALL (not only used) canonical solutions of M(n+1)

#define USE_MY_QUICK_SORT			(USE_THREADS > 1 || 1)      // For multithead version we cannot use regular quicksort, since we are using some global variable here

#define SLIP_TIME					500  //    in microseconds
#define REPORT_INTERVAL				(5 * 1000 * 100)
#define REPORT_INTERVAL_OBJ_NUMB	10000

#define PRINT_MATRICES				0
#define SOLUTION_STATISTICS			0
#define WRITE_MULTITHREAD_LOG		0

#if TEST
#define PRINT_TO_FILE				1
#define PRINT_SOLUTIONS				0
#define PRINT_CURRENT_MATRIX		1
#else
#define PRINT_TO_FILE				1	// Write files with the results for each set of parameters
#define PRINT_SOLUTIONS				0
#define PRINT_CURRENT_MATRIX		0
#endif

#define PRINT_PERMUTATION			0   // Output of the current permutation during the canonicity check
#define PRINT_CANON_GROUP			0
#define CHECK_CCC					0 // 151315
#define CHECK_CONSTRUCTED			0


#define USE_MUTEX     (USE_THREADS > 1 && (PRINT_SOLUTIONS || PRINT_PERMUTATION || PRINT_PERMUTATION))
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

#define FCLOSE(file)			if (file) fclose(file)

#define SNPRINTF(x, len, ...)			snprintf(x, len, __VA_ARGS__)
#define SPRINTF(x, ...)			      SNPRINTF(x, sizeof(x), __VA_ARGS__)
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

#define MATRIX_ELEMENT_TYPE  	uchar
#define MATRIX_ELEMENT_IS_BYTE	(MATRIX_ELEMENT_TYPE == uchar)
#if MATRIX_ELEMENT_IS_BYTE
#define _FRMT				"u"
#define ME_FRMT				"%" _FRMT
#define MATRIX_ELEMENT_MAX	UINT8_MAX 
#endif

#define VECTOR_ELEMENT_TYPE  	uchar
#define VECTOR_ACCESS_TYPE		VECTOR_ELEMENT_TYPE
typedef CArray<VECTOR_ELEMENT_TYPE, VECTOR_ACCESS_TYPE> CArrayOfVectorElements;
#define VECTOR_ELEMENT_TYPE_MAX	0xff

#define PERMUT_ELEMENT_TYPE  	size_t
#define PERMUT_ACCESS_TYPE		PERMUT_ELEMENT_TYPE
typedef CArray<PERMUT_ELEMENT_TYPE, PERMUT_ACCESS_TYPE> CArraySolutionPerm;
#define PERMUT_ELEMENT_MAX		UINT64_MAX


template <class T>
class CSimpleArray {
public:
    CC inline CSimpleArray(size_t len) : m_nLen(len){ m_pElem = new T[len]; }
	CC virtual ~CSimpleArray()						{ delete [] elementPntr(); }
	inline void Init(size_t len, T *pElem)			{ m_pElem = pElem; m_nLen = len; }
    CC inline T element(size_t idx) const           { return *(elementPntr() + idx); }
	inline void setElement(size_t idx, T val)		{ *GetElement(idx) = val; }
    CC inline T *elementPntr() const				{ return m_pElem; }
    CC inline size_t numElement() const				{ return m_nLen; }
	inline T GetAt(size_t idx)  const				{ return element(idx); }
	inline T *GetElement(size_t idx)  const			{ return elementPntr() + idx; }
protected:
private:
    T *m_pElem;
    size_t m_nLen;
};

template <class T>
class CCounter : public CSimpleArray<T> {
public:
    CC CCounter(size_t len) : CSimpleArray<T>(len)  { resetArray(); }
	CC ~CCounter()									{}
	CC inline void resetArray()						{ memset(this->elementPntr(), 0, this->numElement() * sizeof(T)); }
	CC inline void incCounter(int idx)              { ++*(this->elementPntr() + idx); }
private:
};

template <class T>
class CContainer : public CSimpleArray<T> {
public:
	CC CContainer(size_t len) : CSimpleArray<T>(len){ resetArray(); }
	CC ~CContainer()								{}
    CC inline void resetArray()                     { m_nNumb = 0; }
	CC inline void addElement(T val)                { *(this->elementPntr() + m_nNumb++) = val; }
    CC inline size_t numb() const                   { return m_nNumb; }
    inline void incNumb()                           { m_nNumb++; }
	inline size_t GetSize() const					{ return numb(); }
private:
	size_t m_nNumb;
};

template <class T>
class CMapping : public CSimpleArray<T> {
public:
    CK CMapping(size_t len) : CSimpleArray<T>(len<<1)  {}
	CK CMapping(T to, T from, size_t len = 1) : CSimpleArray<T>(len << 1)		
													{ resetMapping(); addMapping(to, from); }
	CC ~CMapping()									{}
    CK void addMapping(T to, T from);
	void removeMapping(T to);
	void removeLastMapping(size_t n = 1)			{ m_nMapPos -= n << 1; }
	void restoreLastMapping(size_t n = 1)			{ m_nMapPos += n << 1; }
    CK inline void resetMapping()                   { m_nMapPos = 0; }
	inline uint nElement() const					{ return getMapPosition() >> 1; }
	CK inline const T *getMappingPntr() const       { return this->elementPntr(); }
    CK inline const T *getLastMapping() const		{ return getMappingPntr() + getMapPosition(); }
	const T *findMapping(T to, const T *pTo = NULL, const T *pToLast = NULL) const;
    bool isEmpty() const                            { return !m_nMapPos; }
	inline void adjustElement(int idx, T val = 1)   { *(this->elementPntr() + idx) -= val; }
protected:
	CK inline uint getMapPosition() const           { return m_nMapPos; }
    
    uint m_nMapPos;
};

template <class T>
void CMapping<T>::addMapping(T to, T from)
{
	assert(m_nMapPos < (this->numElement() << 1));
	T *pTo = (T *)getLastMapping();
    *pTo = to;
    *(pTo+1) = from;
    m_nMapPos += 2;
}

template <class T>
const T *CMapping<T>::findMapping(T to, const T *pTo, const T *pToLast) const
{
	// Searching for element to be removed
	if (!pTo)
		pTo = getLastMapping();

	if (!pToLast)
		pToLast = getMappingPntr();

	while ((pTo -= 2) >= pToLast && *pTo != to);

	return pTo >= pToLast ? pTo : NULL;
}

template <class T>
void CMapping<T>::removeMapping(T to)
{
	// Searching for element to be removed
	T *pTo = (T *)findMapping(to);
	if (!pTo)
		return;

	// Element found
	const int len = getLastMapping() - pTo - 2;
	if (len)
		memcpy(pTo, pTo + 2, len * sizeof(T));

	removeLastMapping();
}

#define countof(x)  (sizeof(x)/sizeof(x[0]))

extern int ccc;
extern int nnn;
#if CHECK_CCC
	#define CHECK_VAL(ccc)		(ccc > CHECK_CCC)
#else
	#define CHECK_VAL(ccc)		true
#endif

#if CHECK_CONSTRUCTED
	#define START_NUMBER         9
	#define END_NUMBER			91
	#define CHECK_CONSTR(x, y)    enumInfo()->constrCanonical() >= x && enumInfo()->constrCanonical() <= y
#else
#define CHECK_CONSTR(x, y)		true
#endif

#define MAKE_OUTPUT()			CHECK_VAL(ccc) && CHECK_CONSTR(START_NUMBER, END_NUMBER)

#if PRINT_SOLUTIONS
	#define OUTPUT_SOLUTION(x,file,f)		if (MAKE_OUTPUT()) \
												{ x->printSolutions(file, f); }
#else
    #define OUTPUT_SOLUTION(x,file, f)
#endif

#define OUT_MATRIX(x, y, z, w)				{ MUTEX_LOCK(out_mutex); x->printOut(y, z, w);  MUTEX_UNLOCK(out_mutex); }
#if PRINT_CURRENT_MATRIX
#define OUTPUT_MATRIX(x, y, z)				if (MAKE_OUTPUT()) \
												OUT_MATRIX(x, y, z, ccc)
#else
    #define OUTPUT_MATRIX(x, y, z)
#endif

#if PRINT_PERMUTATION						
	#define OUTPUT_PERMUTATION(x, f, n, p)		MUTEX_LOCK(out_mutex); x->outputPerm(f, n, p);  MUTEX_LOCK(out_mutex); 
#else
    #define OUTPUT_PERMUTATION(x, f, n, p)
#endif

#if PRINT_CANON_GROUP
	#define OUTPUT_CANON_GROUP(cond, x, f)		if (cond) x->outputPermutations(f, x->lenPerm()/*numColOrb()*/)
#else
	#define OUTPUT_CANON_GROUP(cond, x, f)
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

typedef enum {
	t_BIBD,			// default
	t_tDesign,
	t_PBIBD,
	t_InsidenceSystem,
	t_InconsistentGraph
} t_objectType;

typedef enum {
	t_enumDefault			= 0,
	t_IS_enumerator			= 1 << 0,
	t_matrixOwner			= 1 << 1,
	t_noReplicatedBlocks	= 1 << 2,
	t_outColumnOrbits		= 1 << 3,
	t_outStabilizerOrbit	= 1 << 4,	// Keep the orbits of stabilizer of first elements
	t_colOrbitsConstructed  = 1 << 5,
	t_printTransposedMatrix = 1 << 6,
	t_alwaisKeepRowPermute	= 1 << 7,   // Keep forming elements of the Aut(M) acting on the rows of partially constructed matrix M 
	t_allFlags				= -1
} t_EnumeratorFlags;

class CInterStruct {
public:
	inline CInterStruct(int mult = 1)							{ setMult(mult); }
	inline ~CInterStruct()										{ delete Counterparts(); }
	inline const std::vector<int> &lambda() const				{ return iParam[0]; }
	inline const std::vector<int> &lambdaA() const				{ return iParam[1]; }
	inline const std::vector<int> &lambdaB() const				{ return iParam[2]; }
	inline std::vector<int> *lambdaPtr()						{ return iParam; }
	inline std::vector<int> *lambdaAPtr()						{ return iParam + 1; }
	inline std::vector<int> *lambdaBPtr()						{ return iParam + 2; }
	inline std::vector<CInterStruct *> *Counterparts() const	{ return m_pCounterparts; }
	inline bool isValid() const									{ return Counterparts(); }
	inline void InitCounterparts()								{ m_pCounterparts = new std::vector<CInterStruct *>(); }
	inline void setNext(CInterStruct *pntr)						{ m_pNext = pntr; }
	inline CInterStruct *getNext() const						{ return m_pNext; }
	inline void setMult(int val)								{ m_mult = val; }
	inline int mult() const										{ return m_mult; }
private:
	std::vector<int> iParam[3];
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
	t_objectType objType;
	int v;
	int k;
	int r;
	int t;
	int mt_level = 0;		// Matrix row number, where the threads will be launched
	uint outType;			// Flags which define the output information of the task
	uint grpOrder;			// Limits for order of the group of the matrices which will be printed
	size_t threadNumb;		// Number of threads launched to perform task
	bool firstMatr;			// TRUE, when first matrix of the set was not yet outputted
	bool noReplicatedBlocks;// TRUE, when only block designs with no replicated blocks should be constructed
	std::string workingDir; // Current working directory name
	std::string logFile = "";    // 
	size_t rewindLen = 0;   // Length of the portion of log file, which probably will be rewinded
	const std::vector<int> &lambda() const	{ return m_pInterStruct->lambda(); }
	const std::vector<int> &lambdaA() const { return m_pInterStruct->lambdaA(); }
	const std::vector<int> &lambdaB() const { return m_pInterStruct->lambdaB(); }
	inline int lambdaSizeMax() const		{ return m_lambdaSizeMax; }
	inline void setLambdaSizeMax(int val)	{ m_lambdaSizeMax = val; }
private:
	CInterStruct *m_pInterStruct = NULL;
	int m_lambdaSizeMax = 0;// Maximal number of elements in lambda()
	                        // (will be used for formated output)
};

#define VAR_1		1
