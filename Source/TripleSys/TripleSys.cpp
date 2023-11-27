// TripleSys.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <windows.h>
#include "TripleSys.h"

alldata::alldata(int numPlayers, int groupSize, bool useCheckLinks) :
		SizeParam((numPlayers - 1) / (groupSize - 1), numPlayers, groupSize),
		m_np2(numPlayers * numPlayers),
		m_np3(m_np2 * numPlayers),
		m_nGroups(numPlayers / groupSize) {
	m_nLenResults = m_numDays * numPlayers;
	maxResult = new char[m_nLenResults];
	m_pResults = new char[m_nLenResults];
	selPlayers = new char[m_numPlayers];
	tmpPlayers = new char[m_numPlayers];
	indexPlayer = new char[m_numPlayers];
	m_pLinks = new char[m_np2];
	index6 = new char[m_numDays];
	if (groupSize == 3 && useCheckLinks)
		m_pCheckLink = new CChecklLink(m_numDays, m_numPlayers);

	Init();
}

alldata::~alldata() {
	delete[] maxResult;
	delete[] m_pResults;
	delete[] selPlayers;
	delete[] tmpPlayers;
	delete[] indexPlayer;
	delete[] m_pLinks;
	delete[] index6;
	delete[] bg;
	delete m_pCheckLink;
}

void alldata::Init() {
	int lastRep = GetTickCount();
	memset(links(), unset, m_numPlayers * m_numPlayers);
	memset(result(), 0, m_nLenResults);
	maxDays = -1;
	nLoops = 0;
	noResult = false;
	iDay = 0;
	iPrevResult = unset;
}

bool alldata::Run() {
	int iTime = GetTickCount();
	while (nLoops < LoopsMax)
	{
		while (iDay < m_numDays)
		{
			if (iPrevResult >= 0)
			{
				iDay = iPrevResult;
				initPrevDayProcess();
				iPrevResult = unset;
				if (iDay < 0)
				{
					noResult = true;
					goto finish;
				}
			}
			else
				initCurrentDay();

			while (iPlayer < m_numPlayers)
			{
				int iPlayerNumber = 0;
				if (iPlayer < 0)
				{
					int k = -1;
					if (iDay <= 0)
					{
						noResult = true;
						goto finish;
					}
					else
					{
						memset(result(iDay), 0, m_numPlayers);
						iDay--;
						iPrevResult = iDay;
						goto runday;
					}
					abort();
				}

				if (iPlayer >= m_numPlayers)
					abort();

				if (UseLastSix && iPlayer == m_numPlayers - 6)
				{
					iPlayerNumber = processLastSix();
					if (iPlayer == m_numPlayers)
						goto nextDay;
				}
				else
					iPlayerNumber = getNextPlayer();
				
				if (iPlayerNumber >= m_numPlayers)
				{
					//printf("\n%d  go back", iDay);
					getPrevPlayer();
					continue;
				}
				indexPlayer[iPlayer] = iPlayerNumber;
				tmpPlayers[iPlayer] = iPlayerNumber;
				if ((iPlayer + 1) < m_numPlayers)
					indexPlayer[(iPlayer + 1)] = 0;
#if UseBitmask
				bmask = bmask & (~(1 << iPlayerNumber));
#endif
				selPlayers[iPlayerNumber] = iPlayer;
				if (PrintGroupStat)
				    setStat();
				iPlayer++;
			}
nextDay:
			if (PrintGroupStat)
				printStat();
			if (noResult)
				break;

			memcpy(result(iDay), tmpPlayers, m_numPlayers);

			if (maxDays < iDay)
			{
				/**/
				printf("day %d  Time = %d\n", iDay, GetTickCount() - iTime);
				printTable("Result", result(), m_numDays, m_numPlayers, 0, GroupSize, true);
				//printTable("Links", links[0], m_numPlayers, m_numPlayers);
				/**/
				maxDays = iDay;
				memcpy(maxResult, result(0), m_nLenResults);
			}
			if (iDay > 0)
			{
				const int iCheck = check1(result(), iDay + 1);
				if (iCheck == 0)
				{
					// get new matrix
					iPrevResult = iDay - 1;
					continue;
				}
			}
			iDay++;
		runday:;
		}
	finish:

		nLoops++;

		if (noResult)
		{
			printf("no more results\n");
			break;
		}

		//report result
		printf("Result %d, Time = %d\n", nLoops, GetTickCount() - iTime);
		//printTable("Links", links(), m_numPlayers, m_numPlayers);
		printTable("Result table", result(), m_numDays, m_numPlayers, 0, GroupSize, true);
		if (nLoops >= LoopsMax)
			break;
		iDay--;
		iPrevResult = iDay;
	}
	printf("Total time = %d ms\n", GetTickCount() - iTime);
	if (nLoops == 1)
	{
		if (memcmp(maxResult, result(), m_nLenResults) != 0)
			printTable("'Maximum days' Result", maxResult, m_numDays, m_numPlayers);
		printTable("Links", links(), m_numPlayers, m_numPlayers);
		convertLinksToResult(links());
		if (memcmp(m_co, result(), m_nLenResults) != 0)
			printTable("Result from link (different than result)", m_co, m_numDays, m_numPlayers);
	}

	return true;
}

void alldata::initCurrentDay()
{
	iPlayer = 0;
	memset(indexPlayer, 0, m_numPlayers);
	memset(selPlayers, unset, m_numPlayers);
	memset(tmpPlayers, unset, m_numPlayers);
#if UseBitmask
	bmask = (1 << m_numPlayers) - 1;
#endif
	//checkbmask(selPlayers, bmask);
}


bool alldata::setLinksAndDevCounts(char* p, int ip, char iset)
{
	const int i = ip % GroupSize;
	if (i == 0)
		return true;
	//	if (iset != 1 && p[ip] < 0)
	//		return true;
	char bset = iset == 1 ? iDay : unset;
	char i1 = p[ip];
	//	if (i1 < 0 || i1 >= m_numPlayers)
	//		i1 = i1;
	auto* linkPtr = links(i1);
	if (bset != unset)
	{
		for (int j = 1; j <= i; j++)
		{
			char i2 = p[ip - j];
			//		if (i2 < 0 || i2 >= m_numPlayers)
			//			i2 = i2;
			if (linkPtr[i2] != unset)
				return false;

		}
	}
	for (int j = 1; j <= i; j++)
	{
		const char i2 = p[ip - j];
		linkPtr[i2] = *(links(i2) + i1) = bset;
	}
	return true;
}

void alldata::initPrevDayProcess()
{
	if (m_pCheckLink)
	{
		if (iDay < 0)
			return;
	}
	else
	{
		if (iDay < 1)  // keep first line (first day result)
		{
			iDay = -1;
			return;
		}
	}

	auto* const pRes = result(iDay);
	memcpy(indexPlayer, pRes, m_numPlayers);
	memcpy(tmpPlayers, pRes, m_numPlayers);

	if (UseLastSix)
	{
		iPlayer = m_numPlayers - GroupSize * 2;
		indexPlayer[iPlayer] = index6[iDay]; // only one of six indices used, last 5 values of array indexPlayers not used
	}
	else
		iPlayer = m_numPlayers - GroupSize - 1;

	if (iPlayer < 0)
	{
		iDay = -1;
		return;
	}

	int ind = indexPlayer[iPlayer];
	if (ind < 0)
		abort();
	else
	{
		memset(selPlayers, 1, m_numPlayers);
#if UseBitmask
		bmask = 0;
#endif

		for (int j = iPlayer; j < m_numPlayers; j++)
		{
			if (*(pRes+j) < 0)
				abort();
			if (!setLinksAndDevCounts(pRes, j, unset))
				abort();
			int k = tmpPlayers[j];
#if UseBitmask
			bmask = bmask | (1 << k);
#endif
			tmpPlayers[j] = selPlayers[k] = unset;
			indexPlayer[j] = 0;
		}
	}
	//checkbmask(s);
	if (UseLastSix && iPlayer == m_numPlayers - 6)
		index6[iDay] = ind + 1;
	indexPlayer[iPlayer] = ind + 1;
}

bool alldata::initStartValues(const char* ivcb, bool printStartValues)
{
	char *iv = m_co; // We can use existing array m_co
	int v;
	int ind = 0;
	int id = iDay = 0;
	memset(iv, unset, m_nLenResults);

	for (ind = 0; ; ind++)
	{
		while (*ivcb == ' ' || *ivcb == '\n')
		{
			if (*ivcb++ == '\n')
			{
				if (id >= m_numDays - 1)
					goto doneInit;
				ind = 0;
				id++;
			}
		}

		if (sscanf_s(ivcb, "%d", &v) != 1)
			break;
		if (ind >= m_numPlayers)
		{
			ind = 0;
			if (++id >= m_numDays)
				break;
		}
		*(iv + id * m_numPlayers + ind) = (char)v;
		while ((*ivcb >= '0' && *ivcb <= '9') || *ivcb == '-')
			ivcb++;
	}
doneInit:
	if (ind <= 0 && id <= 0)
		return false;

	char* iv_id = iv;
	auto* res = result(0);
	for (int i = 0; i < id; i++, iv_id += m_numPlayers, res += m_numPlayers)
	{
		iDay = i;
		for (int j = 0; j < m_numPlayers; j++)
		{
			const auto ivId = iv_id[j];
			if (ivId == unset)
			{
				printf("Init: value for day %d position %d not devined\n", i, j);
				printTable("Initial result", result(0), m_numDays, m_numPlayers);
				exit(0);
			}

			*(res+j) = ivId;
			if (!setLinksAndDevCounts(res, j, 1))
			{
				printf("Init: value of %d (for day %d position %d) already devined in links table\n",
					ivId, i, j);
				printTable("Initial result", result(0), m_numDays, m_numPlayers);
				exit(0);
			}
		}
		index6[i] = getLastSixIndex(res);
	}

	iDay = id;
	if (printStartValues) {
		printTable("Result", result(), numDays(), numPlayers());
		printTable("Links", links(0), numPlayers(), numPlayers());
	}
	return true;
}

void alldata::getPrevPlayer()
{
	if (iPlayer >= m_numPlayers)
		abort();

	indexPlayer[iPlayer] = 0;
	while (--iPlayer >= 0)
	{
		int iPlayerNumber;
		if (!setLinksAndDevCounts(tmpPlayers, iPlayer, unset))
			abort();
		if (tmpPlayers[iPlayer] < 0 || tmpPlayers[iPlayer] >= m_numPlayers)
			abort();
		iPlayerNumber = tmpPlayers[iPlayer];
		//checkbmask(selPlayers, bmask);
		selPlayers[iPlayerNumber] = unset;
#if UseBitmask
		bmask = bmask | (1 << iPlayerNumber);
#endif
		//checkbmask(selPlayers, bmask);
		tmpPlayers[iPlayer] = unset;

		if (UseLastSix && iPlayer >= m_numPlayers - 6)
			abort();

		if (iPlayer >= m_numPlayers - GroupSize || iPlayerNumber >= m_numPlayers)
		{
			indexPlayer[iPlayer] = 0;
			continue;
		}

		iPlayerNumber++;
		indexPlayer[iPlayer] = iPlayerNumber;
		break;
	}
}

int alldata::getNextPlayer()
{
	int iPlayerNumber = indexPlayer[iPlayer];
	int m0 = iPlayer % GroupSize;
	int m1 = m0 == 0 ? GroupSize : 1;
	int ifixedPlayer = -1;

	if (!m_pCheckLink)
	{
		if (iDay == 0)
		{
			if (iPlayerNumber > iPlayer)
				return m_numPlayers;
			iPlayerNumber = iPlayer;
		}
		else if (GroupSize == 3)
		{
			if (iDay == 1 && iPlayer > 0 && iPlayer < 3)
			{
				ifixedPlayer = iPlayer * 3;
				if (ifixedPlayer >= iPlayerNumber)
				{
					tmpPlayers[iPlayer] = ifixedPlayer;
					if (setLinksAndDevCounts(tmpPlayers, iPlayer, 1))
						return ifixedPlayer;
					tmpPlayers[iPlayer] = unset;
				}
				return m_numPlayers;
			}
			if (iPlayer < 7 && m0 == 0 && iDay > 0)
			{
				ifixedPlayer = iPlayer / 3;
				if (ifixedPlayer >= iPlayerNumber)
					return ifixedPlayer;
				return m_numPlayers;
			}
			else if (iDay == 3 && iPlayer == 2 && iPlayerNumber <= result(iDay - 1)[iPlayer])
				iPlayerNumber = result(iDay - 1)[iPlayer] + 1;

			/** attemt to sort days by the second columnt **/
			else if (iPlayer == 1 && iDay > 0)
			{
				if (iPlayerNumber <= result(iDay - 1)[iPlayer])
					iPlayerNumber = result(iDay - 1)[iPlayer] + 1;
			}
			/**/
		}
	}

	if (iPlayer >= m1)
	{
		// the following 2 lines makes signinficant change in speed
		if (iPlayerNumber <= tmpPlayers[iPlayer - m1])
			iPlayerNumber = tmpPlayers[iPlayer - m1] + 1;


		/** p7 > p4 in the Second day
		4 <= z1 <= 11, 5 <= z2 <= 14,
		and if z1 <= 8, then z2 <= 11 **/
		if (GroupSize == 3 && iDay == 1)
		{
			if (iPlayer == 4)
			{
				if (iPlayerNumber < 4)
					iPlayerNumber = 4;
				if (iPlayerNumber > 11)
					return m_numPlayers;
			}
			else if (iPlayer == 7)
			{
				if (iPlayerNumber <= tmpPlayers[4])
					iPlayerNumber = tmpPlayers[4] + 1;
				if (iPlayerNumber > 11 && tmpPlayers[4] <= 8)
					return m_numPlayers;
				if (iPlayerNumber > 14)
					return m_numPlayers;
			}
		}
	}

	if (m0 == 0)
	{
		//new
		if (iPlayerNumber > iPlayer)
			return m_numPlayers;
#if UseBitmask
		int msk = bmask;
		int minBit = msk ^ (msk - 1); // last bit and '1' in all bits below
		int iBit = 1 << iPlayerNumber;

		for (; iPlayerNumber < m_numPlayers; iPlayerNumber++, iBit += iBit)
		{
			if (selPlayers[iPlayerNumber] != unset)
				continue;
			/**/
			if (minBit < iBit)
			{
				iPlayerNumber = m_numPlayers;
				break;
			}/**/
			//tmpPlayers[iPlayer] = iPlayerNumber;
			break;
		}
#else
		int firstNotSel = 0;
		for (int i = 0; i < m_numPlayers; i++)
		{
			if (selPlayers[i] == unset)
			{
				firstNotSel = i;
				break;
			}
		}
		if (iPlayerNumber > firstNotSel)
			return m_numPlayers;
		else
			return firstNotSel;
#endif
	}
	else
	{
		char i0 = iPlayer > 0 ? tmpPlayers[iPlayer - 1] : 0;
		char i1 = iPlayer > 1 ? tmpPlayers[iPlayer - 2] : 0;
		char* l0 = links(i0);
		char* l1 = links(i1);
		char day = iDay;
		for (; iPlayerNumber < m_numPlayers; iPlayerNumber++)
		{
			if (selPlayers[iPlayerNumber] != unset)
				continue;
			if (m_numPlayers) { // GroupSize == 3
				if (l0[iPlayerNumber] != unset)
					continue;
				char* li = links(iPlayerNumber);
				if (m0 == 2)
				{
					if (l1[iPlayerNumber] != unset)
						continue;
					l1[iPlayerNumber] = li[i1] = day;
				}
				l0[iPlayerNumber] = li[i0] = day;

				/**/
				if (m0 == 2 && m_pCheckLink)
				{
					if (!m_pCheckLink->checkLinks(links(), iDay))
					{
						li[i0] = li[i1] = unset;
						l0[iPlayerNumber] = l1[iPlayerNumber] = unset;
						continue;
					}
				}
				/**/
				break;
			} else {
				tmpPlayers[iPlayer] = iPlayerNumber;
				if (setLinksAndDevCounts(tmpPlayers, iPlayer, 1))
				{
					break;
				}
				tmpPlayers[iPlayer] = unset;
			}
		}
	}

	return iPlayerNumber;
}

#define sl(a, b) *(links(a)+b) = *(links(b)+a) = iDay;
#define ul(a, b) *(links(a)+b) = *(links(b)+a) = unset
#define l(a, b) (*(links(a)+b) == unset)
#define process6(cs, v0, v1, v2, v3, v4, v5) { \
			if (l(v0, v1) && l(v1, v2) && l(v0, v2) && l(v3, v4) && l(v4, v5) && l(v3, v5)) \
			{ cact = cs; r[0] = v0; r[1] = v1; r[2] = v2; r[3] = v3; r[4] = v4; r[5] = v5; break; }}
#define index6(cs, v0, v1, v2, v3, v4, v5)  \
			if ((vt[0] == v0) && (vt[1] == v1) && (vt[2] == v2) && (vt[3] == v3) && (vt[4] == v4) && (vt[5] == v5)) \
			    return cs

void alldata::getLastSix(char* v)
{
	int j = 6;
	for (int i = m_numPlayers; --i >= 0;)
	{
		if (selPlayers[i] == unset)
		{
			if (--j < 0)
				abort();
			v[j] = i;
		}
	}
	if (j != 0)
		abort();
}

int alldata::getLastSixIndex(const char *resDay)
{
	char v[6], vt[6];
	memcpy(v, resDay + m_numPlayers - 6, 6);
	memcpy(vt, v, 6);
	for (int j = 1; j < 6; j++)
	{
		for (int i = 0; i < 6 - j; i++)
		{
			if (v[i + 1] < v[i])
			{
				char tmp = v[i];
				v[i] = v[i + 1];
				v[i + 1] = tmp;
			}
		}
	}
	index6(0, v[0], v[1], v[2], v[3], v[4], v[5]);
	index6(1, v[0], v[1], v[3], v[2], v[4], v[5]);
	index6(2, v[0], v[1], v[4], v[2], v[3], v[5]);
	index6(3, v[0], v[1], v[5], v[2], v[3], v[4]);
	index6(4, v[0], v[2], v[3], v[1], v[4], v[5]);
	index6(5, v[0], v[2], v[4], v[1], v[3], v[5]);
	index6(6, v[0], v[2], v[5], v[1], v[3], v[4]);
	index6(7, v[0], v[3], v[4], v[1], v[2], v[5]);
	index6(8, v[0], v[3], v[5], v[1], v[2], v[4]);
	index6(9, v[0], v[4], v[5], v[1], v[2], v[3]);
	abort();
}

int alldata::processLastSix()
{
	char v[6];
	char cact = 0;
	const auto ip = (unsigned char)iPlayer;
	char* r = tmpPlayers + ip;
	char c = *(indexPlayer + ip);
	if (c > 9)
		return m_numPlayers;
	else if (c < 0)
		abort();
	getLastSix(v);
	cact = -1;
	while (c <= 10)
	{
		switch (c)
		{
		case 0: process6(0, v[0], v[1], v[2], v[3], v[4], v[5])
		case 1: process6(1, v[0], v[1], v[3], v[2], v[4], v[5])
		case 2: process6(2, v[0], v[1], v[4], v[2], v[3], v[5])
		case 3: process6(3, v[0], v[1], v[5], v[2], v[3], v[4])
		case 4: process6(4, v[0], v[2], v[3], v[1], v[4], v[5])
		case 5: process6(5, v[0], v[2], v[4], v[1], v[3], v[5])
		case 6: process6(6, v[0], v[2], v[5], v[1], v[3], v[4])
		case 7: process6(7, v[0], v[3], v[4], v[1], v[2], v[5])
		case 8: process6(8, v[0], v[3], v[5], v[1], v[2], v[4])
		case 9: process6(9, v[0], v[4], v[5], v[1], v[2], v[3])
		default:
			return m_numPlayers;
		};
		if (cact < 0)
			return m_numPlayers;
		sl(r[0], r[1]);
		sl(r[0], r[2]);
		sl(r[1], r[2]);
		sl(r[3], r[4]);
		sl(r[3], r[5]);
		sl(r[4], r[5]);
		if (!m_pCheckLink || m_pCheckLink->checkLinks(links(), iDay))
			break;
		ul(r[0], r[1]);
		ul(r[0], r[2]);
		ul(r[1], r[2]);
		ul(r[3], r[4]);
		ul(r[3], r[5]);
		ul(r[4], r[5]);
		c = cact + 1;
	}
	index6[iDay] = cact;
	indexPlayer[ip] = cact;
	return iPlayer = m_numPlayers;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
