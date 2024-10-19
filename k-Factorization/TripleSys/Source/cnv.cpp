#include "TripleSys.h"

const char *getFileNameAttr(const kSysParam* param, const char** uf) {
	static const char *fhdr;
	if (uf)
		*uf = "";

	fhdr = param->val[t_p1f] ? "P" : "K";
	if (param->val[t_u1f]) {
		fhdr = "U";
		if (uf)
			*uf = param->strVal[t_UFname]->c_str();
	}
	if (param->val[t_nestedGroups])
		return *fhdr == 'K' ? "KM" : *fhdr == 'P' ? "PM" : "UM";
	return fhdr;
}

CC bool checkSet(const tchar* tr, int nt)
{
	int n = 1<<tr[0];
	for (int i = 1; i < nt; i++)
	{
		const int m = 1<<tr[i];
		if (n & m)
			return false;
		n |= m;
	}
	return true;
}

CC void alldata::cnvInit()
{
	m_cnvMode = 0;
	memcpy(m_groups, result(), m_groupSize);

	const auto groupSizeMinus1 = m_groupSize - 1;
	auto* groups_ig = m_groups;
	for (int i = 1; i < m_groupSizeFactorial; i++)
	{
		groups_ig += m_groupSize;
		memcpy(groups_ig, groups_ig - m_groupSize, m_groupSize);
		do {
			auto j = m_groupSize;
			while (groups_ig[--j] >= groupSizeMinus1)
				groups_ig[j] = 0;

			groups_ig[j]++;
		} while (!checkSet(groups_ig, m_groupSize));
	}
}
CC int alldata::cnvCheckKm1(ctchar* tr, int nrows, tchar* pOrbits)
{
	int ret = 1;
	tchar ttr1[MAX_PLAYER_NUMBER], *ttr = (tchar *)tr;
	const auto * res = result();
	tchar n = 0, day = 0;
#if !USE_CUDA
	if (m_cnvMode) {
		cnvPrintAuto(tr, nrows);
		return 1; // print only
	}
#endif
	if (!m_TrInd && !groupOrder())
	{
		//printf("cnvCheckStart\n");
		updateGroup(res);
		day = 1;
	}
	for (; day < m_NumDaysToTransform; day++)
	{
		if (day) {
			ttr = ttr1;
			const auto* resn = result(n = m_DayIdx[day]);
			for (int i = m_numPlayers; i--;)
			{
				ttr[resn[i]] = tr[i];
			}
		}

		const auto icmp = (this->*m_pProcessMatrix)(res, ttr, nrows, n);
#if !USE_CUDA
		//TestkmProcessMatrix(nrows, n, tr, ttr, icmp);
#endif
		Stat_cnvCheckKm1("cmp(2)", 4, icmp == 2);
		Stat_cnvCheckKm1("cmp(all)", 5, true);

		if (icmp == 0) {
#if NEED_TO_DEBUG
			if (incGroupOrder())
				orbits()->UpdateOrbits(ttr);
#else
			//printTable("ttr", ttr, 1, 21, 3);
			updateGroup(ttr);
#endif
#if USE_EQUAL
#if 0
			if (pOrbits) {
				if (!ret) {
					;
				}
				else
					memcpy(pOrbits, ttr, m_numPlayers);
			}
#endif
			ret = 0;
			if (day) 
				m_DayIdx[day--] = m_DayIdx[--m_NumDaysToTransform];
#endif
			continue;
		}

		if (icmp < 0) {
#if PRINT_TRANSFORMED
			printTransformed(nrows, m_numPlayers, tr, ttr, res, m_Km, n, nLoops, m_finalKMindex);
#endif
			ret = -1;
			break;
		}
	}
	Stat_cnvCheckKm1("can(all)", 0, true);
	Stat_cnvCheckKm1("(-1)", 1, ret == -1);
	return ret;
}
CC bool alldata::cnvCheckKm(ctchar* tr, ctchar* tg, int nrows)
{
	auto* trmk = m_trmk;
	for (int i = 0; i < m_nGroups; i++, trmk += m_groupSize)
	{
		const auto itr = tr[i];
		const auto* pGroups = m_groups + tg[i] * m_groupSize;
		for (int j = 0; j < m_groupSize; j++)
		{
			trmk[j] = itr + pGroups[j];
		} 
	}
	const bool ret = cnvCheckKm1(m_trmk, nrows, NULL) >= 0;

	m_TrInd++;
	return ret;
}
CC bool alldata::cnvCheckTgNew(ctchar* tr, int nrows)
{
	tchar tg[MAX_GROUP_NUMBER];
	bool ret = true;
	memset(tg, 0, m_nGroups);
	while(1)
	{
		if (!cnvCheckKm(tr, tg, nrows))
		{
			ret = false;
			break;
		}
		int i = m_nGroups;
		while (i-- && ++tg[i] >= m_groupSizeFactorial)
		{
			tg[i] = 0;
		}
		if (i < 0)
			break;
	}
	return ret;
}
CC bool alldata::cnvCheckNew(int iMode, int nrows, bool useAutomorphisms)
{
	m_TrInd = 0;
	resetGroupOrder();
	auto* pTestedTRs = testedTRs();
	if (pTestedTRs)
		pTestedTRs->resetGroupOrder();

	if (useAutomorphisms && utilizeGroups(nrows)) {
		int step, lastVal, j, nGroupsTested = 0;
		int i = setupIteratorByGroups(&lastVal, &step);
		int playerIndexCycle = nrows * m_numPlayers - m_groupSize - 1;
		while (i != lastVal) {
			auto* m_pRowGroup = rowGroup(i += step);
			if (!m_pRowGroup || (j = m_pRowGroup->groupOrder()) < 1)
				continue;
			
			const auto nRowsToTest = nrows - i;
			if (nRowsToTest < 1)
				continue;

			const auto* pMatrToTest = result(i);
			auto* cmpTr = m_pRowGroup->getObject();
			while (j--) {
				if (pTestedTRs && pTestedTRs->isProcessed(cmpTr))
				{
					cmpTr += m_numPlayers;
					continue;
				}
				m_playerIndex = playerIndexCycle - i * m_numPlayers;
				const auto cmp = kmProcessMatrix(pMatrToTest, cmpTr, nRowsToTest);
				if (cmp < 0) {
					m_playerIndex += i * m_numPlayers;
					return false;
				}
				if (!cmp)
					updateGroup(cmpTr);
				cmpTr += m_numPlayers;
			}

			if (++nGroupsTested == param(t_autGroupNumb))
				break;
		}
	}

	m_playerIndex = nrows * m_numPlayers - m_groupSize - 1;

	if (param(t_nestedGroups) > 1 && nrows > 2)
	{
		updateGroup(result(0));
		if (groupOrder() >= param(t_nestedGroups)) {
			saveGroup(*this, nrows);
			if (nrows < m_numDaysResult)
				return true;
		}
		else
			return false;
	}

	bool ret = true;
	if (m_p1f && iMode >= 0) { // leo && !param(t_u1f)) {
		ret = (this->*m_pCheckP1F)(nrows);
	}
	else {
		tchar a[MAX_GROUP_NUMBER + 1], p[MAX_GROUP_NUMBER + 1];
		initCheckByGroup(nrows, 1);
		m_cnvMode = iMode;

		// change "ng = m_nGroups;" below (for "fast version" of cnfCheckNew) to: 
		//int ng = m_nGroups > 4 ? 8 : m_nGroups; // check permutations for the first ng groups only
		int ng = m_nGroups; // "full check"

		// Head Permutations Using a Linear Array Without Recursion by Phillip Paul Fuchs
		for (auto i = m_nGroups; i--;)   // initialize arrays; a[N] can be any type
		{
			a[i] = (p[i] = i) * m_groupSize;   // a[i] value is not revealed and can be arbitrary
		}
		if (!cnvCheckTgNew(a, nrows))
		{
			ret = false;
		}
		else
		{
			p[m_nGroups] = ng;//m_nGroups; // p[N] > 0 controls iteration and the index boundary for i
			auto i = 1;   // setup first swap points to be 1 and 0 respectively (i & j)
			while (i < ng)//m_nGroups)
			{
				p[i]--;             // decrease index "weight" for i by one
				const auto j = (i % 2) ? p[i] : 0;   // IF i is odd then j = p[i] otherwise j = 0
				SWAP(a[i], a[j]);
				if (!cnvCheckTgNew(a, nrows))
				{
					ret = false;
					break;////
				}
				i = 0;              // reset index i to 1 (assumed)
				while (!p[++i])     // while (p[i] == 0)
				{
					p[i] = i;       // reset p[i] zero value
				}
			}
		}
	}

	if (ret)
		saveGroup(*this, nrows);

	return ret;
}
