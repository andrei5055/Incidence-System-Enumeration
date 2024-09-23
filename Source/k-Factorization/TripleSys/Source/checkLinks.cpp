#include "TripleSys.h"
#include <iostream>

#define SetLinksFault(i) \
		faults_id[i] |= 1; counts_id[i]++; m_cntErr++

int CChecklLink::getAllUnlinked(int ic, char* v, int nvMax)
{
	auto* ci = m_pLinksCopy + ic * m_numPlayers;
	char nv = 1;

	for (int j = 0; j < m_numPlayers && nv < nvMax; j++)
	{
		if (ci[j] == unset && ic != j)// && j > ic)
		{
			v[nv] = j;
			nv++;
		}
	}
	if (nv > nvMax - 1)
		abort();
	v[0] = nv - 1;
	v[nv] = unset;
	return nv;
}


CChecklLink::CChecklLink(int numDays, int numPlayers, int groupSize) :
		SizeParam(numDays, numPlayers, groupSize) {
	const auto len = numDays * numPlayers;
	initArray(&m_counts, len);
	initArray(&m_faults, len);
	initArray(&m_tmfalse, len);
	initArray(&m_tmok, len);
	initArray(&m_pLinksCopy, numPlayers * numPlayers);
	initArray(&m_v, numPlayers);
	initArray(&m_vo, numPlayers);

#if PrintNVminmax
	initArray(&nvmn, len, (char)99);
	initArray(&nvmx, len, (char)(-1));
#endif
}

CChecklLink::~CChecklLink() {
	delete[] m_counts;
	delete[] m_faults;
	delete[] m_tmfalse;
	delete[] m_tmok;
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
	char* lnks;
#define LT 0
#if UseSS == 0 && LT == 0
	lnks = pLinks;
#else
	memcpy(m_pLinksCopy, pLinks, len);
	lnks = m_pLinksCopy;
#endif
	m_cnt++;

	const auto idx = id * m_numPlayers;
	auto *faults_id = m_faults + idx;
	auto* counts_id = m_counts + idx;
#if UseSS == 0
	for (int i0 = 0; i0 < m_numPlayers; i0++)
	{
		int i = (i0 + 5) % m_numPlayers;
#else
	for (int i0 = 0; i0 < m_numPlayers - m_groupSize; i0 += m_groupSize)
	{
		int i = i0;
#endif
		auto *ci = lnks + i * m_numPlayers;
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
			printTableColor("Links", lnks, m_numPlayers, m_numPlayers);
			abort();
		}

		if (checkLinksV(lnks, m_v, nv, -1, m_vo))
		{
#if UseSS == 0 && LT == 0
			goto okplayer;
#endif
#if LT
			if (i0 > 1)
				goto okplayer;
#endif
			int idd = 0;
			for (int n = 0; n < nv; n += 2)
			{
				const auto a = m_vo[n];
				const auto b = m_vo[n + 1];

				if (a == unset || b == unset)
					abort();
				auto* ca = lnks + a * m_numPlayers;
				if (ci[a] != unset || ci[b] != unset || ca[b] != unset)
					abort();
				idd = -2;// m_numDays - nv / 2 + n;

				auto* cb = lnks + b * m_numPlayers;
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
		m_cntErr++;
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
		m_cntOk++;
		/**
		if (id > 3)
		{
			printTableColor("CheckLinks Links", lnks, m_numPlayers, m_numPlayers);
			convertLinksToResult(lnks, m_co, m_numPlayers, m_groupSize);
			printTable("CheckLinks Result", m_co, m_numDays, m_numPlayers);
		}
		**/
	}

	return ret;
}
#if 0
#define Stat_checkLinksT Stat
#else
#define Stat_checkLinksT(a, b, c) 
#endif
bool CChecklLink::checkLinksT(char* pLinks, int id, bool printLinksStatTime)
{
	Stat_checkLinksT("all checkLinksT", 0, true);
	if (id < 2)//4,5,6
	{
		Stat_checkLinksT("id < 1", 1, true);
		return true;
	}

	//return true;
	const auto len = m_numPlayers * m_numPlayers;
	m_cnt++;
	const auto idx = id * m_numPlayers;
	auto* faults_id = m_faults + idx;
	auto* counts_id = m_counts + idx;
	char v[256];
	int nv, nv2 = 0, nv3 = 0, nvMax = sizeof(v);
	for (int i0 = 0; i0 < m_numPlayers; i0 += 3)
	//int i0 = 0;
	{
		int i = i0;
		memcpy(m_pLinksCopy, pLinks, len);
		if ((nv = getAllUnlinked(i, v, nvMax)) < 3
			|| (nv2 = getAllUnlinked(i + 1, v + nv, nvMax - nv)) < 3
			|| (nv3 = getAllUnlinked(i + 2, v + nv + nv2, nvMax - nv - nv2)) < 3
			)
			return true;
		if ((nv & 1) == 0 || (nv2 & 1) == 0 || (nv3 & 1) == 0)
		{
			printTable("Links", m_pLinksCopy, m_numPlayers, m_numPlayers);
			printf("i=%d nv=%d %d %d\n", i, nv, nv2, nv3);
			return false;
			abort();
		}
		int nvAll = nv + nv2 + nv3;
		v[nvAll + 1] = 0;
		if (!checkLinksTR(v, nvAll, 0, -1))
		{
			SetLinksFault(0);
			Stat_checkLinksT("false", 2, true);
			Stat_checkLinksT("i>=3 false", 3, i>=3);
			return false;
		}
	}
	faults_id[0] |= 2;
	m_cntOk++;
	return true;
}
