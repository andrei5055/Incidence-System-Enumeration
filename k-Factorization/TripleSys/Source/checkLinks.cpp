#include "TripleSys.h"

CC CChecklLink::CChecklLink(const SizeParam & sizeParam, const kSysParam* p) : CChecklLinkBase(sizeParam) {
	set_kSysParam(p);
	const auto len = m_numDays * m_numPlayers;
	initArray(&m_pP1Ftable, len);
	initArray(&m_pLinksCopy, m_numPlayers * m_numPlayers);
	initArray(&m_v, m_numPlayers);
	initArray(&m_vo, m_numPlayers);

#if PrintNVminmax
	initArray(&nvmn, len, (char)99);
	initArray(&nvmx, len, (char)(-1));
#endif
}

CC CChecklLink::~CChecklLink() {
	delete[] m_pP1Ftable;
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

CC bool CChecklLink::checkLinks(tchar *pLinks, int id, bool printLinksStatTime)
{
	int idMin = 3;
	if (id < idMin) return true;

	bool ret = true;
	tchar* lnks = pLinks;

	PrepareIDs();
	int iOffset = 5;
	int ie = m_numPlayers;
	for (int i0 = 0; i0 < ie; i0++)
	{
		int i = (i0 + iOffset) % m_numPlayers;
		auto *ci = lnks + i * m_numPlayers;
		int nv = 0;
		for (int j = 0; j < m_numPlayers; j++)
		{
			if (ci[j] == unset && i != j)
				m_v[nv++] = j;
		}
		if (!nv)
			continue;

		//if (nv >= (m_numDays - 1) * 2) // this check makes it faster
		//	continue;

#ifndef USE_CUDA
		if (nv % 2)
		{
			//continue;
			printf("CheckLinks: error in links table: nv=%d for i=%d\n", nv, i);
			printTableColor("Links", lnks, m_numPlayers, m_numPlayers, m_groupSize);
			abort();
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
			break;
		}
	}
	UpdateCntOK(ret)
	return ret;
}
