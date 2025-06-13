#include "kOrbits.h"

CKOrbits::CKOrbits(uint outGroupMask, int numElems, int groupSize, int numRows) : 
    m_numElems(numElems), m_numRows(numRows), m_groupSize(groupSize),
    RowGenerators(outGroupMask, numRows * numElems / groupSize) { 
    m_outMask = 16;
    m_sActionOn = "k-sets, |Aut(K)|";

    unsigned int len = numElems;
    int i = groupSize;
    while (--i)
        len *= numElems;

    m_pTable = new tchar[len];
    m_pSolution = new tchar[len = m_numElems / groupSize * numRows];
    auto pntr = getNextObject();
    for (tchar i = 0; i < len; i++)
        pntr[i] = i;
}

CKOrbits::~CKOrbits() {
    delete[] m_pTable;
    delete[] m_pSolution;
}

void CKOrbits::createTable(ctchar* pSolution) {
    tchar id = 0;
    const auto numGroups = m_numElems / m_groupSize;
    for (int i = 0; i < m_numRows; i++) {
        for (int j = 0; j < numGroups; j++) {
            unsigned int idx = *pSolution++;
            for (int k = m_groupSize; --k;) {
                idx *= m_numElems;
                idx += *pSolution++;
            }

            m_pTable[idx] = id++;
        }
    }
}

void CKOrbits::encodeSolution(ctchar* pSolution) {
    tchar id = 0;
    const auto numGroups = m_numElems / m_groupSize;
    for (int i = 0; i < m_numRows; i++) {
        for (int j = 0; j < numGroups; j++) {
            unsigned int idx = *pSolution++;
            for (int k = m_groupSize; --k;) {
                idx *= m_numElems;
                idx += *pSolution++;
            }

            m_pSolution[id++] = m_pTable[idx];
        }
    }
}

void CKOrbits::makeGroupOutput(const CGroupInfo* pElemGroup, bool outToScreen, bool checkNestedGroups) {
    createTable(((alldata*)pElemGroup)->result());
    RowGenerators::makeGroupOutput(pElemGroup, outToScreen, checkNestedGroups);
}

