#pragma once

#ifndef CD_TOOLS
#include "TripleSys.h"

#define TFunc2(x, ...)          template<typename T, typename S> __VA_ARGS__ x
#define Class2(x)               x<T,S>
#define Class2Def(x)            TFunc2(x, class)
#define FClass2(x, ...)			TFunc2(Class2(x), __VA_ARGS__)

#define countof(x)     sizeof(x)/sizeof(x[0])
#else
#include "../DataTypes.h"
#endif

#define USE_2_ROW_CANON 0		// Canonizer will use only 2 rows     

#define CheckerCanon(...)		FClass2(CCheckerCanon, __VA_ARGS__)

typedef enum {
	t_notReady			 = 0,
	t_readyToExplainTxt  = 1,
	t_readyToExplainMatr = 1 << 1,
	t_readyCompletely    = 255,
} t_bResultFlags;

Class2Def(CCheckerCanon) {
public:
	CCheckerCanon(T nRow, T nCol, T groupSize = GroupSize)
		: m_numElem(nCol), m_numElem2(2 * nCol), m_numDaysMax(nRow), m_groupSise(groupSize) {
		m_players = new T[2 * m_numElem];
		m_tmpBuffer = new T[m_numElem + nRow];
		m_pResutMemory = new T[(m_numElem + 1) * nRow];
		initCommentBuffer(256);
	}
	~CCheckerCanon()						{ delete[] m_players;
											  delete [] getTmpBuffer();
											  delete[] resultMemory();
											  resetComments();
	}
	bool CheckCanonicity(const T* result, int nLines, T *bResult=NULL);
	inline auto numDays() const				{ return m_numDays; }
	inline auto comment() const				{ return m_pComment; }
	inline bool improvedResultIsReady(t_bResultFlags flag = t_bResultFlags::t_readyCompletely) const {
											  return (flag & m_bResultFlag) == flag; }
private:
	inline auto groupSize() const			{ return m_groupSise; }
	auto stabiliserLengthExt() const		{ return m_nStabExtern; }
	void setStabiliserLengthExt(T len)		{ m_nStabExtern = len; }
	bool copyTuple(const T* res, T inc = 0, bool doCheck = true) const;
	bool rollBack(T* p_dayRes, T* p_dayIsUsed, int& j, int nDays) const;
	inline void setNumDays(T nDays)			{ m_lenResult = (m_numDays = nDays) * m_numElem; }
	inline void setResultOut(T* pntr)		{ m_pResultOut = pntr; }
	inline auto resultOut() const			{ return m_pResultOut; }
	inline auto getTmpBuffer() const		{ return m_tmpBuffer; }
	inline auto resultMemory() const		{ return m_pResutMemory;  }
	inline auto lenResult()	const			{ return m_lenResult; }
	inline void resetImprovedResultFlag()   { m_bResultFlag = t_bResultFlags::t_notReady; }
	inline void addImproveResultFlags(t_bResultFlags flags) { m_bResultFlag |= flags;  }
	int checkDay_1(const T* result, int iDay, T *pDest);
	bool checkDay(const T* res, T iDay, T numGroup, T* pNumReason);
	void orderigRemainingDays(T daysOK, T groupsOK, T numGroup, T *pDest) const;
	bool permutPlayers4Day(const T* p_players, const T* resDayIn, T numGroup, T* resDayOut) const;
	bool reportTxtError(T* bBuffer, const char* pReason, T* pDays = NULL, T nDays = 2);
	inline void resetComments()				{ delete[] m_pComment; m_pComment = NULL; }
	inline void initCommentBuffer(int len)  { resetComments(); m_pComment = new char[m_nCommentBufferLength = len]; }
	inline auto commentBufferLength() const { return m_nCommentBufferLength; }

	T m_nStabExtern = 0;		// number of first elements of permutation which Canonicity Checker will not move
	T* m_players = NULL;
	const T m_numElem;			// The number of elements that will be the same for all partially constructed objects
								// (it is equal nCol for combinatorial designs or number of players for k-system)
	const T m_numElem2;			// This is twice the number of elements. 
	const T m_numDaysMax;
	const T m_groupSise;
	size_t m_lenResult;
	T m_numDays;
	T* m_pResultOut;
	T* m_tmpBuffer = NULL;		// Buffer uswd for groups and days ordering
	T* m_pResutMemory = NULL;	// Memory allocated to improve results
	unsigned int m_bResultFlag;
	int m_nCommentBufferLength = 0;
	char *m_pComment = NULL;
};


