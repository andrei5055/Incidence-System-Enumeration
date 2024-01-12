#pragma once
#include "TripleSys.h"

#ifndef CD_TOOLS
#define SIZE_TYPE				unsigned char
#define ELEMENT_MAX				static_cast<SIZE_TYPE>(-1)
#define CC
#else
#include "../DataTypes.h"
#endif

#include "GroupOrder.h"

#define TFunc1_(x, ...)			template<typename T> __VA_ARGS__ x
#define Class1_(x)              x<T>
#define Class1Def_(x)           TFunc1_(x, class)
#define FClass1_(x, ...)		TFunc1_(Class1_(x), __VA_ARGS__)

#define IDX_MAX					(ELEMENT_MAX - 1)
#define CheckerCanon(...)		FClass1_(CCheckerCanon, __VA_ARGS__)

typedef enum {
	t_notReady			 = 0,
	t_readyToExplainTxt  = 1,
	t_readyToExplainMatr = 1 << 1,
	t_readyCompletely    = 255,
} t_bResultFlags;

template<typename T>
inline void revert(T* perm, T j, T i) {
	while (++i < --j) perm[i] ^= (perm[j] ^= (perm[i] ^= perm[j]));
}

Class1Def_(CCheckerCanon) : public Class1_(CGroupOrder) {
public:
	CCheckerCanon(T nRow, T nCol, T groupSize = GroupSize)
		: m_numElem(nCol), m_numElem2(2 * nCol), m_numDaysMax(nRow), 
		  m_groupSise(groupSize), m_numGroups(nCol/ groupSize), m_lenRow(m_numElem*sizeof(T)) {
		m_players = new T[5 * nCol];
		m_tmpBuffer = new T[nCol + nRow];
		m_pResultMemory = new T[(nCol + 1) * nRow];
		m_pOrbits = (m_pPermutation = m_players + 2 * nCol) + nCol;
		m_tmpBuffer1 = m_pOrbits + nCol;
		initCommentBuffer(256);
	}
	~CCheckerCanon()						{ delete[] m_players;
											  delete [] tmpBuffer();
											  delete[] resultMemory();
											  resetComments();
											}
	bool CheckCanonicity(const T* result, int nLines, T *bResult=NULL);
	bool CheckPermutations(const T* result, const T* pMatrix, int nRows);
	inline auto numDays() const				{ return m_numDays; }
	inline auto comment() const				{ return m_pComment; }
	inline bool improvedResultIsReady(t_bResultFlags flag = t_bResultFlags::t_readyCompletely) const {
											  return (flag & m_bResultFlag) == flag; }
	inline void setPreordered(bool v = true) { m_bPreordered = v; }
private:
	inline auto numElem() const				{ return m_numElem; }
	inline auto groupSize() const			{ return m_groupSise; }
	auto stabiliserLengthExt() const		{ return m_nStabExtern; }
	void setStabiliserLengthExt(T len)		{ m_nStabExtern = len; }
	bool copyTuple(const T* res, T inc = 0, bool doCheck = true) const;
	bool rollBack(T* p_dayRes, T* p_dayIsUsed, int& j, int nDays) const;
	inline auto setStudiedMatrix(const T* pntr, T nDays) {
		m_pStudiedMatrix = pntr;
		m_lenResult = (m_numDays = nDays) * m_numElem;
	}
	inline auto studiedMatrix() const		{ return m_pStudiedMatrix; }
	inline auto getMatrixRow(T nRow) const  { return studiedMatrix() + nRow * m_numElem; }
	inline void setResultOut(T* pntr)		{ m_pResultOut = pntr; }
	inline auto resultOut() const			{ return m_pResultOut; }
	inline auto tmpBuffer() const			{ return m_tmpBuffer; }
	inline auto resultMemory() const		{ return m_pResultMemory; }
	inline auto lenResult()	const			{ return m_lenResult; }
	inline void resetImprovedResultFlag()   { m_bResultFlag = t_bResultFlags::t_notReady; }
	inline void addImproveResultFlags(t_bResultFlags flags) { m_bResultFlag |= flags;  }
	inline auto numGroups() const			{ return m_numGroups;}
	inline auto lenRow() const				{ return m_lenRow; }
	inline auto destMemory() const			{ return m_pDestMemory; }
	bool checkDay_1(int iDay);
	bool checkDay(T iDay);
	int checkDayCode(int diff, T iDay, const T* secontRow);
	void orderigRemainingDays(T daysOK, T groupsOK, T *pDest) const;
	bool reportTxtError(T* bBuffer, const char* pReason, T* pDays = NULL, T nDays = 2);
	inline void resetComments()				{ delete[] m_pComment; m_pComment = NULL; }
	inline void initCommentBuffer(int len)  { resetComments(); m_pComment = new char[m_nCommentBufferLength = len]; }
	inline auto commentBufferLength() const { return m_nCommentBufferLength; }
	inline auto preordered() const			{ return m_bPreordered; }
	inline void setNumReason(T numReason)	{ m_numReason = numReason; }
	inline auto numReason() const			{ return m_numReason; }
	inline void setReasonParam(T val)       { m_nReasonParam = val; }
	inline auto reasonParam() const			{ return m_nReasonParam; }
	void createDaySequence(T iDay = 1) const;
	bool checkOrderingForDay(T iDay) const;
	bool checkRemainingDays(T iDay, int retVal = -1, const T* pPerm = NULL);
	bool checkPosition1_4(const T* players);
	bool explainRejection(const T* players, T playerPrevID, T playerNewID, T firstDayID = 0, bool doOutput = false, const T* pNewOrder = NULL);
	int orderingMatrix(T nDays, T numGroups, bool expected = true, bool invert = false, const T* permPlayer = NULL);
	void sortTuples(T *players) const;
	inline auto permutation() const			{ return m_pPermutation; }
	inline auto oprbits() const				{ return m_pOrbits; }
	inline void setTrivialPerm(const T* p)  { m_pTrivialPerm = p; }
	inline auto trivialPerm() const			{ return m_pTrivialPerm; }
	inline auto playersPerm(int idx) const  { return m_players + idx * m_numElem; }
	bool checkPermutationOfFirstDayGroups(int numGroups, const T* pCurrentRow, bool useRecording = false);
	bool checkWithGroup(T numElem, int (CCheckerCanon<T>::*func)(const T*, T, const T*), const T* pCurrentRow = NULL);
	int checkPermutationOnGroups(const T* permGroups, T numElem, const T* pCurrentRow);
	int orderingMatrix(const T* permut, T numElem)     {
		return orderingMatrix(0, 0, false, false, permut);
	}
	inline void recordTuples(const T* pTuples, T *pPlayers) const {
		for (T j = 0; j < numElem(); j++)
			pPlayers[pTuples[j]] = j;
	}

	T nextPermutation(T* perm, const T* pOrbits, T nElem, T idx = ELEMENT_MAX, T lenStab = 0);

	T m_nStabExtern = 0;		// number of first elements of permutation which Canonicity Checker will not move
	T* m_players = NULL;
	const T m_numElem;			// The number of elements that will be the same for all partially constructed objects
								// (it is equal nCol for combinatorial designs or number of players for k-system)
	const T m_numElem2;			// This is twice the number of elements. 
	const T m_numDaysMax;
	const T m_groupSise;
	const T m_numGroups;
	T m_numReason;
	T m_nReasonParam;

	const size_t m_lenRow;
	size_t m_lenResult;
	T m_numDays;
	const T* m_pStudiedMatrix = NULL;
	T* m_pResultOut;
	T* m_tmpBuffer = NULL;		// Buffer to be used for groups and days ordering
	T* m_tmpBuffer1 = NULL;     // Buffer to be used in ordering by groups
	T* m_pResultMemory = NULL;	// Memory allocated to improve results
	T* m_pDestMemory = NULL;	// Memory for recorded matrix (equal to bResults OR resultMemory())
	T* m_pPermutation = NULL;
	T* m_pOrbits = NULL;
	const T* m_pTrivialPerm;

	unsigned int m_bResultFlag;
	int m_nCommentBufferLength = 0;
	char *m_pComment = NULL;
	bool m_bPreordered = true;
};


