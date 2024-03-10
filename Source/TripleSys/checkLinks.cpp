#include "TripleSys.h"
#include <iostream>

void alldata::setCheckLinks() {

}

CChecklLink::CChecklLink(int numDays, int numPlayers, int groupSize) :
		SizeParam(numDays, numPlayers, groupSize) {
	const auto len = numDays * numPlayers;
	initArray(&counts, len);
	initArray(&faults, len);
	initArray(&tmfalse, len);
	initArray(&tmok, len);
	initArray(&m_pLinksCopy, numPlayers * numPlayers);
	initArray(&m_v, numPlayers);
	initArray(&m_vo, numPlayers);

#if PrintNVminmax
	initArray(&nvmn, len, (char)99);
	initArray(&nvmx, len, (char)(-1));
#endif
}

CChecklLink::~CChecklLink() {
	delete[] counts;
	delete[] faults;
	delete[] tmfalse;
	delete[] tmok;
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

bool CChecklLink::checkLinks(char *pLinks, int id, bool printLinksStatTime)
{
	bool ret = true;
	const auto len = m_numPlayers * m_numPlayers;

	if (m_numPlayers == 15 && id < 3)
		return true;
	else if (m_numPlayers == 21 && id < 2) //1)
		return true;
	else if (m_numPlayers == 27 && id < 3)
		return true;

	memcpy(m_pLinksCopy, pLinks, len);

	cnt++;

	const auto idx = id * m_numPlayers;
	auto *faults_id = faults + idx;
	auto* counts_id = counts + idx;
#if UseSS == 0
	for (int i0 = 0; i0 < m_numPlayers; i0++)
	{
		int i = (i0 + 5) % m_numPlayers;
#else
	for (int i0 = 0; i0 < m_numPlayers - m_groupSize; i0 += m_groupSize)
	{
		int i = i0;
#endif
		auto *ci = m_pLinksCopy + i * m_numPlayers;
		int nv = 0;

		for (int j = 0; j < m_numPlayers; j++)
		{
			if (ci[j] == unset && i != j)
				m_v[nv++] = j;
		}
		if (nv == 0)
			continue;

		if (nv >= (m_numDays - 1) * 2) // this check makes it faster
			continue;

		if ((nv % 2) != 0)
		{
			//continue;
			printf("CheckLinks: error in links table: nv=%d for i=%d\n", nv, i);
			printTableColor("Links", m_pLinksCopy, m_numPlayers, m_numPlayers);
			abort();
		}

		if (checkLinksV(m_pLinksCopy, m_v, nv, -1, m_vo))
		{
#if UseSS == 0
			goto okplayer;
#endif
			int idd = 0;
			for (int n = 0; n < nv; n += 2)
			{
				const auto a = m_vo[n];
				const auto b = m_vo[n + 1];

				if (a == unset || b == unset)
					abort();
				auto* ca = m_pLinksCopy + a * m_numPlayers;
				if (ci[a] != unset || ci[b] != unset || ca[b] != unset)
					abort();
				idd = -2;// m_numDays - nv / 2 + n;

				auto* cb = m_pLinksCopy + b * m_numPlayers;
				ci[a] = ca[i] = idd;
				ci[b] = cb[i] = idd;
				cb[a] = ca[b] = idd;
			}
			goto okplayer;
		}
	//fltPlayer:
		if (id < 3)
			id = id;
		setNV_MinMax(id, i, nv);
		faults_id[i] |= 1;
		counts_id[i]++;
		cntErr++;
		ret = false;
		break;
	okplayer:
		setNV_MinMax(id, i, nv);
		faults_id[i] |= 2;
	}
	if (ret)
	{
		if (id == 6)
			id = id;
		cntOk++;
		/**/
		if (id > 33)
		{
			printTableColor("CheckLinks Links", m_pLinksCopy, m_numPlayers, m_numPlayers);
			convertLinksToResult(m_pLinksCopy, m_co, m_numPlayers, m_groupSize);
			printTable("CheckLinks Result", m_co, m_numDays, m_numPlayers);
		}
		/**/
	}

	return ret;
}