#pragma once
#include "CudaAttributes.h"
#include "k-SysSupport.h"
#include "Table.h"

class CKOrbits : public RowGenerators<tchar> {
public:
	CKOrbits(uint outGroupMask, int numElems, int groupSize, int nRows);
	~CKOrbits();
	int createGroupAndOrbits(const CRepository<tchar>* pElemGroup) override;
	CC void UpdateGroup(ctchar* pSolution) {
#if !USE_CUDA
		encodeSolution(pSolution);
		m_pRowGroup->updateGroup((ctchar*)m_pSolution);
#endif
	}
protected:
	void createTable(ctchar* pSolution) override;
	int createGroup(const CRepository<tchar>* pElemGroup) override {
		return testNestedGroups(pElemGroup, NULL, ((alldata*)pElemGroup)->numDaysResult(), this);
	}
	void createOrbitsSet(const CRepository<tchar>* pElemGroup) override;
	void generatorOutput(bool outToScreen, const char* pErr = NULL) override {
		m_pKOrbGenerators->setName(name());
		m_pKOrbGenerators->makeGroupOutput(NULL, outToScreen, false);
	}
private:
	void encodeSolution(ctchar* pSolution);

	const int m_numElems;
	const int m_numRows;
	const int m_groupSize;
	Generators<ushort>* m_pKOrbGenerators = NULL;
	ushort* m_pTable = NULL;
	ushort* m_pSolution = NULL;
	size_t m_len = 0;
};

