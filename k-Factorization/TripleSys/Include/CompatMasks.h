#pragma once
#include <vector>
#include "k-SysSupport.h"
#include "Storage.h"

#define SAME_MASK_IDX		0			// We allow the same mask index to be used for three consecutive rows. 
// If 0, we will not apply the acceleration method that analyzes valid solutions for the remaining rows in such situations. 
#define USE_64_BIT_MASK		!USE_CUDA
#define UseSolutionCliques	!USE_CUDA	
										// The graph whose vertices are the remaining solutions must have a maximum 
										// clique whose size is equal to the number of unconstructed rows of the matrix.

typedef long long ll;

#if USE_64_BIT_MASK
typedef ll tmask;
#define SHIFT						6
#else
typedef tchar tmask;
#define SHIFT						3
#endif

#define MASK_BIT(idx)				((tmask)1 << ((idx) & ((1<<SHIFT) - 1)))	
#define IDX(n)						(n + (1<<SHIFT) - 1) >> SHIFT
#define REM(n)						(n % ((tmask)1<<SHIFT))			// remainder from division
#define SET_MASK_BIT(mask, idx)		(mask)[(idx) >> SHIFT] |= MASK_BIT(idx)
#define RESET_MASK_BIT(mask, idx)	(mask)[(idx) >> SHIFT] ^= MASK_BIT(idx)
#define CHECK_MASK_BIT(mask, idx)	((mask)[(idx) >> SHIFT] & MASK_BIT(idx))


#include "CompSolGraph.h"

class alldata;

class CCompatMasks {
public:
	CCompatMasks(const kSysParam* pSysParam, const alldata* pAllData);
	virtual ~CCompatMasks();
	CC inline auto numSolutionTotalB() const				{ return m_numSolutionTotalB; }
	CC inline tmask* getSolutionMask(uint solNumb) const	{ return rowsCompatMasks() + (solNumb + m_solAdj) * lenSolutionMask(); }
	CC inline auto numLongs2Skip(int iRow) const			{ return m_pNumLongs2Skip[iRow]; }
	CC inline const auto rowSolutionMasksIdx() const		{ return m_pRowSolutionMasksIdx; }
	CC inline auto numDaysResult() const					{ return m_numDaysResult; }
	CC inline auto numPreconstructedRows() const			{ return m_numPreconstructedRows; }
	CC inline void reset()									{ m_numObjects = 0; }
	CC inline auto getNumSolution() const					{ return m_numObjects; }
	CC inline const auto rowSolutionMasks() const			{ return m_pRowSolutionMasks; }
	CC inline auto maskTestingCompleted() const				{ return m_pMaskTestingCompleted; }
	CC void releaseSolMaskInfo();
	inline void releaseCompatMaskMemory()					{ delete[] rowsCompatMasks(); }
	void initMaskMemory(uint numSolutions, int lenUsedMask, int numRecAdj = 0, int numSolAdj = 0);
	CC inline auto selectPlayerByMask() const				{ return m_bSelectPlayerByMask; }
	inline auto numMasks() const							{ return m_numMasks; }
	CC inline auto numPlayerSolutionsPtr() const			{ return m_pPlayerSolutionCntr; }
	inline auto lenSolutionMask() const						{ return m_lenSolutionMask; }
	CC void initRowUsage(tmask** ppCompatibleSolutions, bool* pUsePlayersMask) const;
protected:
	inline void setNumSolutions(uint numSol)				{ m_numSolutionTotal = numSol; }
	inline void resetSolutionMask(uint idx) const			{ memset(rowsCompatMasks() + lenSolutionMask() * idx, 0, numSolutionTotalB()); }
	inline auto rowsCompatMasks() const						{ return m_pRowsCompatMasks; }
private:
protected:
	uint m_numSolutionTotal;
	uint m_numSolutionTotalB; // length of one solution mask in bytes
	uint* m_pPlayerSolutionCntr = NULL;
	uint* m_pNumLongs2Skip = NULL; // Pointer to the number of long long's that we don't need to copy for each row.
	// For each row of the matrix, we define two masks, each containing an interval of consecutive bits set to 1
	// These intervals represent the row's first and last sets of solutions that lie outside the separately tested 64-bit intervals.
	tmask* m_pRowSolutionMasks = NULL;
	//  ... and the the set of indices of the long long elements which corresponds to two mask's sets.                           
	uint* m_pRowSolutionMasksIdx = NULL;
	uint m_numObjects;
private:
	const int m_numPreconstructedRows;     // Number of preconstructed matrix rows
	const int m_numDaysResult;
	const bool m_bSelectPlayerByMask;      // Find players by mask of unused players
	uint m_lenSolutionMask;				   // length of one solution mask in tmask units 
	uint m_numMasks;
	int m_solAdj = 0;
	tmask* m_pRowsCompatMasks = NULL;
	bool* m_pMaskTestingCompleted = NULL;
};

class CCompressedMask : public CCompatMasks {
public:
	CCompressedMask(const kSysParam* pSysParam, const alldata* pAllData) : CCompatMasks(pSysParam, pAllData) {}
	~CCompressedMask()						{ releaseSolIndices(); }
	void compressCompatMasks(tmask* pCompSol, uint first, const CCompatMasks* pCompMask);
private:
	inline auto solIndices() const			{ return m_pSolIdx; }
	inline void releaseSolIndices()	const	{ delete[] solIndices(); }

	uint* m_pSolIdx = NULL;					// Indices of the solutions stored in compressed mask
};