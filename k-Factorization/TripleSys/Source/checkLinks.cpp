#include "TripleSys.h"

CC CChecklLink::CChecklLink(const SizeParam & sizeParam, const kSysParam* p) : CChecklLinkBase(sizeParam) {
	set_kSysParam(p);
	const auto len = m_numDays * m_numPlayers;
	m_remainder3 = new tchar[m_numPlayers];
	for (int i = 0; i < m_numPlayers; i++)
		m_remainder3[i] = i % 3;
	initArray(&m_pU1Ftable, len);
	initArray(&m_pLinksCopy, m_numPlayers * m_numPlayers);
	initArray(&m_v, m_numPlayers);
	initArray(&m_vo, m_numPlayers);

#if PrintNVminmax
	initArray(&nvmn, len, (char)99);
	initArray(&nvmx, len, (char)(-1));
#endif
}

CC CChecklLink::~CChecklLink() {
	delete[] m_remainder3;
	delete[] m_pU1Ftable;
	delete[] m_pLinksCopy;
	delete[] m_v;
	delete[] m_vo;
#if PrintNVminmax
	delete[] nvmn;
	delete[] nvmx;
#endif
}

#if PrintNVminmax
void CChecklLink::setNV_MinMax(int id, int idx, char nv) {
	idx += id * m_numPlayers;
	if (nvmn[idx] > nv)
		nvmn[idx] = nv;
	if (nvmx[idx] < nv)
		nvmx[idx] = nv;
}
#endif
CC bool CChecklLink::checkLinks(tchar* pLinks, int id, bool printLinksStatTime)
{
	const auto cbmpGraph = !completeGraph();
	int idMin = cbmpGraph ? 3 : (m_numPlayers < 27 ? 2 : 5);
	if (id < idMin)
		return true;

	bool ret = true;
	tchar* lnks = pLinks;

	PrepareIDs();
	int iOffset = 0;// 5;
	int ie = m_numPlayers;
	int i = -1;
	bool bUseTable = m_numPlayers == 24 && cbmpGraph;
	for (int i0 = 0; i0 < ie; i0++) {
		i = (i0 + iOffset) % m_numPlayers;
		const auto ip3 = m_remainder3[i];
		auto* ci = lnks + i * m_numPlayers;
		int nv = 0;
		for (int j = 0; j < m_numPlayers; j++)
		{
			if (ci[j] == unset && i != j) {
				if (cbmpGraph) {
					if (m_remainder3[j] == ip3)
						continue;
				}
				m_v[nv++] = j;
			}
		}
		if (!nv)
			continue;

#ifndef USE_CUDA
		if (nv % 2)
		{
			//continue;
			printf("CheckLinks: error in links table: nv=%d for i=%d\n", nv, i);
			printTableColor("Links", lnks, m_numPlayers, m_numPlayers, m_groupSize);
			exit(1);
		}
#endif
		if (checkLinksV(lnks, m_v, nv, -1, m_vo))
		{
			setNV_MinMax(id, i, nv);
			UpdateFaultsID(i);
		}
		else
		{
			setNV_MinMax(id, i, nv);
			SetLinksFault(i);
			ret = false;
#if 0
			static int a[24];
			a[i]++;
			if (i == 0 && (a[0] % 10000) == 0) {
				for (int k = 0; k < 24; k++)
					printf(" %3d", a[k] * 100 / a[0]);
				printf("\n");
				//if ((a[0] % 1000) == 0)
				//	memset(a, 0, sizeof(a));
			}
			continue;//???
#endif
			break;
		}
	}
	UpdateCntOK(ret)
	return ret;
}
