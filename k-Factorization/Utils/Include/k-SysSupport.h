#pragma once
#include <string>
#include <wtypes.h>

#ifndef USE_CUDA
#ifdef UTILS_EXPORTS
#    define UTIL_LIBRARY __declspec(dllexport)
#else
#    define UTIL_LIBRARY __declspec(dllimport)
#endif

#ifdef K_SYS_LIBRARY_EXPORTS
#    define K_SYS_LIBRARY_API __declspec(dllexport)
#else
#    define K_SYS_LIBRARY_API __declspec(dllimport)
#endif

#define CC
#define CK
#else
#define CC __host__ __device__		// CUDA_CALLABLE
#if CONSTR_ON_GPU
#define CK	CC __noinline__
#else
#define CK
#endif
#define UTIL_LIBRARY
#define K_SYS_LIBRARY_API
#endif

#define AUT			"|Aut(M)| = "
#define MATR_ATTR	"\n\n"##AUT

#define RedText "\x1b[1;31m"
#define GreenText "\x1b[1;32m"
#define YellowText "\x1b[1;33m"
#define ResetTextColor "\x1b[0m"
#define printfRed(fmt, ...) printf(RedText fmt ResetTextColor, __VA_ARGS__)
#define printfGreen(fmt, ...) printf(GreenText fmt ResetTextColor, __VA_ARGS__)
#define printfYellow(fmt, ...) printf(YellowText fmt ResetTextColor, __VA_ARGS__)

void myAssert(int code, const char* file, int line);

typedef const char cchar;
typedef unsigned char tchar;
typedef const tchar  ctchar;
typedef unsigned short ushort;
typedef unsigned int uint;

template<typename T>
CC T* reallocStorageMemory(T** pObjects, size_t lenObj, size_t lenObjPrev = 0, bool bPointers = false) {
	auto* pNewObjMemory = new T[lenObj];
	if (!pNewObjMemory)
		return NULL;

	if (!lenObjPrev)
		lenObjPrev = lenObj >> 1;

	memcpy(pNewObjMemory, *pObjects, lenObjPrev * sizeof(pNewObjMemory[0]));
	if (bPointers) {
		// After pointer reallocation, set to NULL's all unused pointers
		memset(pNewObjMemory + lenObjPrev, 0, (lenObj - lenObjPrev) * sizeof(void *));
	}

	delete[] * pObjects;
	return (*pObjects = pNewObjMemory) + lenObjPrev;
}

// trim from start (in place)
inline void ltrim(std::string& s) {
	s.erase(s.begin(), find_if(s.begin(), s.end(),
		[](int c) { return !isspace(c); }));
}

// trim from end (in place)
inline void rtrim(std::string& s) {
	s.erase(find_if(s.rbegin(), s.rend(),
		[](int c) { return !isspace(c); }).base(), s.end());
}
// trim from both ends (in place)
inline void trim(std::string& s) {
	ltrim(s);
	rtrim(s);
}

class CMatrixInfo {
public:
	CMatrixInfo(uint nMatrReserved) {
		updateReservedMatrNumb(nMatrReserved);
		m_pGroupOrders = new uint[nMatrReserved];
		m_pCycleInfo = new std::string * [nMatrReserved];
		m_pGroupsInfo = new std::string * [nMatrReserved];
		memset(m_pCycleInfo, 0, nMatrReserved * sizeof(m_pCycleInfo[0]));
		memset(m_pGroupsInfo, 0, nMatrReserved * sizeof(m_pGroupsInfo[0]));
	}
	~CMatrixInfo() { 
		delete[] m_pGroupOrders;
		for (auto i = m_nMatrReserved; i--;) {
			if (m_pCycleInfo[i])
				delete m_pCycleInfo[i];
			if (m_pGroupsInfo[i])
				delete m_pGroupsInfo[i];
		}
		delete[] m_pCycleInfo;
		delete[] m_pGroupsInfo;
	}
	inline auto cycleInfo(int idx) const		{ return m_pCycleInfo[idx]->c_str(); }
	inline auto groupInfo(int idx) const		{ return m_pGroupsInfo && m_pGroupsInfo[idx]? m_pGroupsInfo[idx]->c_str() : NULL; }
	inline auto** groupOrdersPntr()				{ return &m_pGroupOrders; }
	inline auto** cycleInfoPntr()				{ return &m_pCycleInfo; }
	inline auto** groupInfoPntr()				{ return &m_pGroupsInfo; }
	inline void updateReservedMatrNumb(int val) { m_nMatrReserved = val; }
private:
	uint* m_pGroupOrders = NULL;
	std::string** m_pCycleInfo = NULL;
	std::string** m_pGroupsInfo = NULL;
	int m_nMatrReserved = 0;
};

UTIL_LIBRARY int readTable(const std::string& fn, int nRows, int nCols, int nmax, int nTotal, tchar** ppSm, int& reservedElement, int nMatricesMax, CMatrixInfo* pMatrixInfos = NULL, char infoSymb = '"');
K_SYS_LIBRARY_API void speakText(LPCWSTR text);

