#include "TripleSys.h"
#if ReportCheckLinksData
int CReportCheckLinksData::getAllUnlinked(int ic, tchar* v, int nvMax) const
{
	auto* ci = m_pLinksCopy + ic * m_numPlayers;
	char nv = 1;
	for (int j = 0; j < m_numPlayers && nv < nvMax - 1; j++)
	{
		if (ci[j] == unset && ic != j)
		//if (ci[j] == unset)
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

bool CReportCheckLinksData::checkLinksTR(const tchar* v, int nvAll, int nv, int ind) const
{
	if (nvAll <= 0)
		return true;
	if (nv <= 0)
	{
		if (ind != -1)
			v++;
		nv = *v++;
		nvAll--;
		ind = -1;
	}
	tchar t[256];
	if (ind == -1)
		memcpy(t, v, nvAll);
	else if (ind == 0)
		memcpy(t, v + 1, nvAll);
	else
	{
		memcpy(t, v, ind);
		if (nvAll > ind)
			memcpy(t + ind, v + ind + 1, nvAll - ind);
	}
	auto* lnkt0 = m_pLinksCopy + t[0] * m_numPlayers;
	for (int i = 1; i < nv; i++)
	{
		auto* lnk = lnkt0 + t[i];
		if (*lnk == unset)
		{
			*lnk = 1;
			if (checkLinksTR(t + 1, nvAll - 2, nv - 2, i - 1))
				return true;
			*lnk = unset;
		}
	}
	return false;
}
bool CReportCheckLinksData::checkLinksT(const tchar* pLinks, int id)
{
	id--;
	Stat_checkLinksT("all checkLinksT", 0, true);
	if (id < 2)//4,5,6
	{
		Stat_checkLinksT("id<1", 1, true);
		return true;
	}
	//return true;
	const auto len = m_numPlayers * m_numPlayers;
	m_cnt++;
	const auto idx = id * m_numPlayers;
	auto* faults_id = m_faults + idx;
	auto* counts_id = m_counts + idx;
#define MAX_V_GR 4
	tchar v[MAX_PLAYER_NUMBER * MAX_V_GR];
	int nv = 0, nv1 = 0, nv2 = 0, nv3 = 0, nvMax = sizeof(v);
	for (int i0 = 0; i0 < m_numPlayers; i0+=1)
		//int i0 = 0;
	{
		int i = i0;
		if (nv == 0)
		    memcpy(m_pLinksCopy, pLinks, len);
		if ((nv1 = getAllUnlinked(i, v + nv, nvMax - nv)) < 3)
			break;
		if ((nv1 & 1) == 0)
		{/**
			printTable("Links", m_pLinksCopy, m_numPlayers, m_numPlayers);
			printf("i=%d nv=%d %d\n", i, nv, nv1);
			//return true;
			abort();
			**/
			break;
		}
		nv += nv1;
		if (nv > m_numPlayers * MAX_V_GR)
		{
			int nvAll = nv;
			v[nvAll + 1] = 0;
			if (!checkLinksTR(v, nvAll, 0, -1))
			{
				SetLinksFault(i);

				Stat_checkLinksT("false(all)", 2, true);

				return false;
			}
			nv = 0;
			faults_id[i] |= 2;
		}
	}
	m_cntOk++;
	return true;
}
#endif