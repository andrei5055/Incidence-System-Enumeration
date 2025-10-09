#pragma once
#include "CudaAttributes.h"
#include "k-SysSupport.h"
#include "Table.h"

class CKOrbits : public RowGenerators<tchar> {
public:
	CKOrbits(uint outGroupMask, int numElems, int groupSize, int nRows);
	~CKOrbits();
	void makeGroupOutput(const CRepository<tchar>* pElemInfo, bool outToScreen = false, bool checkNestedGroups = true) override;
	CC void UpdateGroup(ctchar* pSolution) {
#if !USE_CUDA
		encodeSolution(pSolution);
		m_pRowGroup->updateGroup(m_pSolution);
#endif
	}
protected:
	void createTable(ctchar* pSolution) override;
	int createGroup(const CRepository<tchar>* pElemGroup) override {
		return testNestedGroups(pElemGroup, NULL, ((alldata*)pElemGroup)->numDaysResult(), this);
	}
private:
	void encodeSolution(ctchar* pSolution);

	const int m_numElems;
	const int m_numRows;
	const int m_groupSize;
	tchar* m_pTable = NULL;
	tchar* m_pSolution = NULL;
};

