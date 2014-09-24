#pragma once

#define MAC                 1
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

#define USE_THREADS					0
#define WAIT_THREADS				(USE_THREADS && 0)

#ifdef WIN
	// For Memory leakage detection
	#define _CRTDBG_MAP_ALLOC
	#include <stdlib.h>

	// Assembly is used for Windows only
    #ifndef X64
        #define USE_ASM				0	// 0 or 1  // 09/07/2014 1 is not working since CCanonicityChecker::checkColOrbit() needs to be updated
    #else
		#if !USE_THREADS				// Assembly code IS NOT thread safe for 64-bit mode
			#define USE_ASM			2	// 0 or 2
		#endif
    #endif

	#define USE_BOOST				0	// Fof OS the usage of boost is not implemented
	#define USE_POOL				(USE_BOOST && 0)
	#define OS "Windows"
#else
    #define OS  "Mac OS"
	#define _CrtSetReportMode(x,y)
	#define _CrtDumpMemoryLeaks()
#endif

#define T_DESIGN_ENUM               1   // Enumerating t-designs for t >= 3
#define TEST						0
#define MT_LEVEL					4
#define USE_CANON_GROUP				1	// Reshafle the solutions for (n+1)-th row of matrix M
										// with respect to the Aut G(M(n)) of the first n rows of M
#define USE_STRONG_CANONICITY	(USE_CANON_GROUP && 1)			// Find equivalent solutions for ALL USED canonical solutions of M(n+1)
#define USE_STRONG_CANONICITY_A	(USE_STRONG_CANONICITY && 1)	// Strong canonicity on ALL (not only used) canonical solutions of M(n+1)

#define NO_REPLICATED_BLOCKS		0
#define USE_MY_QUICK_SORT			(USE_THREADS > 1 || 1)      // For multithead version we cannot use regular quicksort, since we are using some global variable here

#define SLIP_TIME					500  //    in microseconds
#define REPORT_INTERVAL				(5 * 1000 * 100)
#define REPORT_INTERVAL_OBJ_NUMB	10000

#define PRINT_MATRICES				0
#define SOLUTION_STATISTICS			0
#define WRITE_MULTITHREAD_LOG		0

#if TEST
#define PRINT_TO_FILE				0
#define PRINT_SOLUTIONS				0
#define PRINT_CURRENT_MATRIX		0
#else
#define PRINT_TO_FILE				0
#define PRINT_SOLUTIONS				1
#define PRINT_CURRENT_MATRIX		1
#endif
#define PRINT_PERMUTATION			0
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
		#define THREAD_NAME_SPACE boost;
	#else
		#define THREAD_NAME_SPACE std;
	#endif
#endif

#ifdef WIN
	#define FOPEN(x, y, z)	  	 FILE *x; fopen_s(&x, y, z)
#else
	#define sprintf_s(x, y, ...) sprintf(x, __VA_ARGS__)
	#define strcpy_s(x, y, ...)  strcpy(x, __VA_ARGS__)
	#define memcpy_s(x, y, ...)  memcpy(x, __VA_ARGS__)
	#define FOPEN(x, y, z)	  	 FILE *x = fopen(y, z)
#endif

#define FCLOSE(file)			if (file) fclose(file)

#define SPRINTF(x, ...)			sprintf_s(x, sizeof(x), __VA_ARGS__) 

typedef unsigned int		uint;
typedef unsigned short		ushort;
typedef unsigned char		uchar;
typedef unsigned long long	ulonglong;

#define MATRIX_ELEMENT_TYPE  	uchar
#define MATRIX_ACCESS_TYPE		MATRIX_ELEMENT_TYPE
typedef CArray<MATRIX_ELEMENT_TYPE, MATRIX_ACCESS_TYPE> CArrayOfMatrixElements;
#define MATRIX_ELEMENT_IS_BYTE	(MATRIX_ELEMENT_TYPE == uchar)

#define VECTOR_ELEMENT_TYPE  	uchar
#define VECTOR_ACCESS_TYPE		VECTOR_ELEMENT_TYPE
typedef CArray<VECTOR_ELEMENT_TYPE, VECTOR_ACCESS_TYPE> CArrayOfVectorElements;

#define PERMUT_ELEMENT_TYPE  	size_t
#define PERMUT_ACCESS_TYPE		PERMUT_ELEMENT_TYPE
typedef CArray<PERMUT_ELEMENT_TYPE, PERMUT_ACCESS_TYPE> CArraySolutionPerm;


template <class T>
class CSimpleArray {
public:
    inline CSimpleArray(size_t len) : m_nLen(len)   { m_pElem = new T[len]; }
    virtual ~CSimpleArray()                         { delete [] elementPntr(); }
    inline T element(int idx) const                 { return *(elementPntr() + idx); }
    inline void setElement(int idx, T val)          { *(elementPntr() + idx) = val; }
    inline T *elementPntr() const                   { return m_pElem; }
    inline size_t numElement() const                { return m_nLen; }
protected:
private:
    T *m_pElem;
    const size_t m_nLen;
};

template <class T>
class CCounter : public CSimpleArray<T> {
public:
    CCounter(size_t len) : CSimpleArray<T>(len)     { resetArray(); }
	inline void resetArray()                        { memset(CSimpleArray<T>::elementPntr(), 0, CSimpleArray<T>::numElement() * sizeof(T)); }
	inline void incCounter(int idx)                 { ++*(CSimpleArray<T>::elementPntr() + idx); }
private:
};

template <class T>
class CContainer : public CCounter<T> {
public:
    CContainer(size_t len) : CCounter<T>(len)       { resetArray(); }
    inline void resetArray()                        { m_nNumb = 0; }
	inline void addElement(T val)                   { *(CSimpleArray<T>::elementPntr() + m_nNumb++) = val; }
    inline int numb() const                         { return m_nNumb; }
    inline void incNumb()                           { m_nNumb++; }
private:
    int m_nNumb;
};

template <class T>
class CMapping : public CSimpleArray<T> {
public:
    CMapping(size_t len) : CSimpleArray<T>(len<<1)  {}
    void addMapping(T to, T from);
    void resetMapping()                             { m_nMapPos = 0; }
	const T *getMappingPntr() const                 { return CSimpleArray<T>::elementPntr(); }
    const T *getLastMappping() const                { return getMappingPntr() + getMapPosition(); }
protected:
    uint getMapPosition() const                     { return m_nMapPos; }
    bool isEmpty() const                            { return !m_nMapPos; }
    
    uint m_nMapPos;
};

template <class T>
void CMapping<T>::addMapping(T to, T from)
{
    T *pTo = (T *)getMappingPntr() + getMapPosition();
    *pTo = to;
    *(pTo+1) = from;
    m_nMapPos += 2;
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
	#define CHECK_CONSTR(x, y)
#endif

#if PRINT_SOLUTIONS
#define OUTPUT_SOLUTION(x,file,f)			if (CHECK_VAL(ccc) CHECK_CONSTR(START_NUMBER, END_NUMBER)) { MUTEX_LOCK(out_mutex); x->printSolutions(file, f); MUTEX_UNLOCK(out_mutex); }
#else
    #define OUTPUT_SOLUTION(x,file, f)
#endif

#if PRINT_CURRENT_MATRIX
#define OUTPUT_MATRIX(x, y, z)				if (PRINT_CURRENT_MATRIX && CHECK_VAL(ccc)  CHECK_CONSTR(START_NUMBER, END_NUMBER)) { MUTEX_LOCK(out_mutex); x->printOut(y, z, ccc);  MUTEX_UNLOCK(out_mutex); }
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

#if USE_THREADS && WRITE_MULTITHREAD_LOG
	void thread_message(int threadIdx, const char *pComment, t_threadCode code, void *pntr = 0)
#else
	#define thread_message(threadIdx, pComment, code, ...)
#endif

void outString(const char *str, FILE *file);
void outString(const char *str, const char *fileName, const char *mode = "a");

