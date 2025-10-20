#include "kOrbits.h"

CKOrbits::CKOrbits(uint outGroupMask, int numElems, int groupSize, int numRows) : 
    m_numElems(numElems), m_numRows(numRows), m_groupSize(groupSize),
    RowGenerators(outGroupMask, numRows * numElems / groupSize, sizeof(ushort)) {
    m_outMask = 16;
    m_sActionOn = "k-sets, |Aut(K)|";

    size_t len = numElems;
    int i = groupSize;
    while (--i)
        len *= numElems;


    m_pTable = new ushort[m_len = len];
    m_pSolution = new ushort[groupDegree()];
    auto pntr = getNextObject();
    for (i = groupDegree(); i--;)
        pntr[i] = i;

    m_pKOrbGenerators = new Generators<ushort>(0, "\nOrbits and group generators for k-sets", groupDegree());
}

CKOrbits::~CKOrbits() {
    delete[] m_pTable;
    delete[] m_pSolution;
    delete m_pKOrbGenerators;
}

void CKOrbits::createTable(ctchar* pSolution) {
    ushort id = 0;
    const auto numGroups = m_numElems / m_groupSize;
    for (int i = 0; i < m_numRows; i++) {
        for (int j = 0; j < numGroups; j++) {
            unsigned int idx = *pSolution++;
            for (int k = m_groupSize; --k;) {
                idx *= m_numElems;
                idx += *pSolution++;
            }

            ASSERT_IF(idx >= m_len);
            m_pTable[idx] = id++;
        }
    }
}

void CKOrbits::encodeSolution(ctchar* pSolution) {
    ushort id = 0;
    const auto numGroups = m_numElems / m_groupSize;
    for (int i = 0; i < m_numRows; i++) {
        for (int j = 0; j < numGroups; j++) {
            unsigned int idx = *pSolution++;
            for (int k = m_groupSize; --k;) {
                idx *= m_numElems;
                idx += *pSolution++;
            }

            ASSERT_IF(idx >= m_len || id >= groupDegree() || m_pTable[idx] == m_pTable[0]);
            m_pSolution[id++] = m_pTable[idx];
        }
    }
    ASSERT_IF(id != groupDegree());
}

int CKOrbits::createGroupAndOrbits(const CRepository<tchar>* pElemGroup) {
    m_pKOrbGenerators->setOutFileName(outFileName(), false);
    m_pKOrbGenerators->resetOrbits();
    setOrbitsCreated(false);
    return RowGenerators::createGroupAndOrbits(pElemGroup);
}

void CKOrbits::createOrbitsSet(const CRepository<tchar>* pElemGroup) {
    // Adding orbits:
    auto *pOrb = m_pKOrbGenerators->getNextObject();
    const auto grDegree = m_pKOrbGenerators->groupDegree();
    for (auto i = grDegree; i--;)
        pOrb[i] = i;

    // ...  and trivial permutation:
    m_pKOrbGenerators->addObject(pOrb);
    const auto groupOrder = pElemGroup->numObjects();
    for (int i = 1; i < groupOrder; i++) {
        const auto* c = (ushort*)pElemGroup->getObject(i);
        m_pKOrbGenerators->addAutomorphism(grDegree, c, pOrb, true, false, true);
    }

    // In fact, we don't need to do that
    m_pKOrbGenerators->updateGroupOrder(grDegree, pOrb);
}
