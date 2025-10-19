#include "TripleSys.h"

const char *getFileNameAttr(const kSysParam* param, const char** uf) {
	static const char *fhdr;
	if (uf)
		*uf = "";

	fhdr = "K";
	const auto& val = param->val;
	if (val[t_u1f]) {
		const auto u1fCycles = param->u1fCycles[0];
		if (!val[t_allowUndefinedCycles] &&
			(!u1fCycles || (u1fCycles[0] == 1 && u1fCycles[1] == val[t_numPlayers])))
			fhdr = "P";
		else
			fhdr = "U";
	}
	if (uf)
		*uf = param->strVal[t_UFname] ? param->strVal[t_UFname]->c_str() : "";

	if (!param->completeGraph()) {
		if (val[t_nestedGroups] > 0)
			return *fhdr == 'K' ? "KCN" : *fhdr == 'P' ? "PCN" : "UCN";
		else
			return *fhdr == 'K' ? "KC" : *fhdr == 'P' ? "PC" : "UC";
	}
	else if (val[t_nestedGroups] > 0)
		return *fhdr == 'K' ? "KN" : *fhdr == 'P' ? "PN" : "UN";
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
	if (m_TrInd == 100000000)
	{
		printf("%d ", m_TrInd);
		m_TrInd = 1;
	}
	//return 1;
	//bool first = false;
	if (!m_TrInd && !orderOfGroup())
	{
		//first = true;
		//printf("cnvCheckStart\n");
		updateGroup(res);
		day = 1;
	}
	const auto cbmpGraph = !completeGraph();
	for (; day < m_NumDaysToTransform; day++)
	{
		if (day) {
			ttr = ttr1;
			const auto* resn = result(n = m_DayIdx[day]);
			if (cbmpGraph) {
				for (int i = 0; i < m_numPlayers; i += m_groupSize)
				{
					for (int j = 0; j < m_groupSize; j++) {
						auto r1 = resn[i + j];
						ttr[r1] = tr[i + m_groupSizeRemainder[r1]];
					}
				}
			}
			else {
				for (int i = m_numPlayers; i--;)
				{
					ttr[resn[i]] = tr[i];
				}
			}
		}

		const auto icmp = (this->*m_pProcessMatrix)(res, ttr, nrows, n, NULL);
		/**
		if (first)
		{
			printf("d=%d icmp=%d ", day, icmp);
			printTable("result", res, nrows, m_numPlayers, m_groupSize);
			printTable("tr\n", ttr, 1, m_numPlayers, m_groupSize);
			printTable("processed", m_Km, nrows, m_numPlayers, m_groupSize);
		}**/
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
			ret = -1;
#if !USE_CUDA //PRINT_TRANSFORMED
			if (param(t_printMatrices) & 16) {
				printTransformed(nrows, m_numPlayers, m_groupSize, tr, ttr, res, m_Km, n, nLoops, m_finalKMindex);
				//continue;
			}
#endif
			if (m_ignoreCanonizationMinus1)
				continue; // Calculate |Aut| and minimum player index to comeback for all such tr's 
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
	const auto step = completeGraph() ? 1 : 0;
	for (int i = 0; i < m_nGroups; i++, trmk += m_groupSize, tg += step)
	{
		const auto itr = tr[i];
		const auto* pGroups = m_groups + *tg * m_groupSize;
		for (int j = 0; j < m_groupSize; j++)
		{
			trmk[j] = itr + pGroups[j];
		} 
	}
	const bool ret = cnvCheckKm1(m_trmk, nrows, NULL) >= 0;

	m_TrInd++;
	return ret;
}
CC bool alldata::cnvCheckTgNew(ctchar* tr, int nrows, int ngroups)
{
	tchar tg[MAX_GROUP_NUMBER];
	bool ret = true;
	memset(tg, 0, m_nGroups);
	const auto ng = completeGraph() ? m_nGroups : 1;

	while(1)
	{
		if (!cnvCheckKm(tr, tg, nrows))
		{
			ret = false;
			if (!m_ignoreCanonizationMinus1)
				break;
		}
		int i = ng;
		while (i-- && ++tg[i] >= m_groupSizeFactorial)
		{
			tg[i] = 0;
		}
		if (i < m_nGroups - ngroups)
			break;
	}
	return ret;
}
CC bool alldata::cnvCheckNew(int iMode, int nrows, bool useAutomorphisms)
{
	const int playerIndexCycle = nrows * m_numPlayers - m_groupSize - 1;
	m_TrInd = 0;
	resetGroupOrder();
	m_playerIndex = playerIndexCycle;

	if (m_precalcMode != eCalculateRows) {
		if (useAutomorphisms && param(t_autGroupNumb) && utilizeGroups(nrows)) {
			int step, lastVal, j, nGroupsTested = 0;
			bool bRet = true;
			int i = setupIteratorByGroups(&lastVal, &step);
			while (i != lastVal) {
				auto* m_pRowGroup = rowGroup(i += step);
				if (!m_pRowGroup || (j = m_pRowGroup->numObjects()) < 1)
					continue;
				const auto i1 = i;
				const auto nRowsToTest = nrows - i;
				if (nRowsToTest < 1)
					continue;

				const auto* pMatrToTest = result(i);

				while (j--) {
					auto* cmpTr = m_pRowGroup->CStorage<tchar>::getObject(j);
					m_playerIndex = playerIndexCycle - i * m_numPlayers;
					const auto cmp = kmProcessMatrix(pMatrToTest, cmpTr, nRowsToTest);
						m_playerIndex += i * m_numPlayers;
					if (cmp < 0) {
						bRet = false;
						if (!m_ignoreCanonizationMinus1)
						return false;
					}
					else if (!cmp)
						updateGroup(cmpTr);
				}
				i = i1;
				if (++nGroupsTested == param(t_autGroupNumb))
					break;
			}
			if (!bRet)
				return false;
		}

		//???m_playerIndex = playerIndexCycle;

		if (m_precalcMode != eCalculateRows && param(t_nestedGroups) > 1 && nrows > 2)
		{
			updateGroup(result(0));
			if (orderOfGroup() >= param(t_nestedGroups)) {
				saveGroup(*this, nrows);
				if (nrows < numDaysResult())
					return true;
			}
			else
				return false;
		}
	}
	return canonizator(iMode, nrows);
}
CC bool alldata::canonizator(int iMode, int nrows)
{
	bool ret = true;
	initCheckByGroup(nrows, 1);
	m_cnvMode = iMode;
	/**/
	if (m_pCheckFunc && iMode >= 0) {
		ret = (this->*m_pCheckFunc)(nrows, nrows);
	}
	else /**/{
		tchar a[MAX_GROUP_NUMBER + 1], p[MAX_GROUP_NUMBER + 1];
		// change "ng = m_nGroups;" below (for "fast version" of cnfCheckNew) to: 
		//int ng = m_nGroups > 4 ? 3 : m_nGroups; // check permutations for the first ng groups only
		int ng = m_nGroups; // "full check"

		// Head Permutations Using a Linear Array Without Recursion by Phillip Paul Fuchs
		for (auto i = m_nGroups; i--;)   // initialize arrays; a[N] can be any type
		{
			a[i] = (p[i] = i) * m_groupSize;   // a[i] value is not revealed and can be arbitrary
		}
		if (!cnvCheckTgNew(a, nrows, ng))
			ret = false;
		if (ret || m_ignoreCanonizationMinus1)
		{
			p[m_nGroups] = ng;//m_nGroups; // p[N] > 0 controls iteration and the index boundary for i
			auto i = 1;   // setup first swap points to be 1 and 0 respectively (i & j)
			while (i < ng)//m_nGroups)
			{
				p[i]--;             // decrease index "weight" for i by one
				const auto j = (i % 2) ? p[i] : 0;   // IF i is odd then j = p[i] otherwise j = 0
				SWAP(a[i], a[j]);
				if (!cnvCheckTgNew(a, nrows, ng)) {
					ret = false;
					if (!m_ignoreCanonizationMinus1)
						break;
				}
				i = 0;              // reset index i to 1 (assumed)
				while (!p[++i])     // while (p[i] == 0)
				{
					p[i] = i;       // reset p[i] zero value
				}
			}
		}
	}

	if (ret) {
		//m_playerIndex = 0;
		if (m_precalcMode != eCalculateRows || nrows < param(t_useRowsPrecalculation))
			saveGroup(*this, nrows);
	}
	return ret;
}
