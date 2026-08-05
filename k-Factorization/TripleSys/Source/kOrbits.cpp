#include "kOrbits.h"

template <typename Func>
void generatePermutations(int n, Func f)
{
    std::vector<tchar> p(n);
    for (int i = 0; i < n; ++i)
        p[i] = i;

    do {
        f(p);
    } while (std::next_permutation(p.begin(), p.end()));
}

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

    generatePermutations(groupSize, [this](const std::vector<tchar>& p) {
        m_pSymGroup.push_back(p);
    });
}

CKOrbits::~CKOrbits() {
    delete[] m_pTable;
    delete[] m_pSolution;
    delete m_pKOrbGenerators;
}

void CKOrbits::createTable(ctchar* pSolution) {
    // Generate a table of indices for all k-subsets of different elements
    ushort id = 0;
    const auto numGroups = m_numElems / m_groupSize;
    for (int i = 0; i < m_numRows; i++) {
        for (int j = 0; j < numGroups; j++, id++, pSolution += m_groupSize) {
            for (const auto perm : m_pSymGroup) {
                unsigned int idx = pSolution[perm[0]];
                for (int k = 0; ++k < m_groupSize;) {
                    idx *= m_numElems;
                    idx += pSolution[perm[k]];
                }

                ASSERT_IF(idx >= m_len);
                m_pTable[idx] = id;
            }
        }
    }
}

int CKOrbits::createGroupAndOrbits(const CRepository<tchar>* pElemGroup) {
    m_pKOrbGenerators->setOutFileName(outFileName());
    m_pKOrbGenerators->resetOrbits();
    setOrbitsCreated(false);
    return RowGenerators::createGroupAndOrbits(pElemGroup);
}

void CKOrbits::updateKSetGroup(ctchar* pMatrix, ctchar* pElemPerm, ctchar* pRowPerm) {
	const auto nGroups = m_numElems / m_groupSize;
    ushort id = 0;
	for (int i = 0; i < m_numRows; i++) {
		auto const* pKset = pMatrix + pRowPerm[i] * m_numElems;
		for (int j = 0; j < nGroups; j++, pKset += m_groupSize) {
            unsigned int idx = pKset[0];
            for (int k = 0; ++k < m_groupSize;) {
                idx *= m_numElems;
                idx += pKset[k];
            }

            m_pSolution[id++] = m_pTable[idx];
		}
	}
    ASSERT_IF(id != groupDegree());
    m_pRowGroup->updateGroup((ctchar*)m_pSolution);
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
