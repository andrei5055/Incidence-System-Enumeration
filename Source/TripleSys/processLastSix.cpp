#include <iostream>
#include "TripleSys.h"

#define sl(a, b) *(links(a)+b) = *(links(b)+a) = iDay;
#define ul(a, b) *(links(a)+b) = *(links(b)+a) = unset
#define l(a, b) (*(links(a)+b) == unset)
#define process6(cs, v0, v1, v2, v3, v4, v5) { \
			if (l(v0, v1) && l(v1, v2) && l(v0, v2) && l(v3, v4) && l(v4, v5) && l(v3, v5)) \
			{ cact = cs; r[0] = v0; r[1] = v1; r[2] = v2; r[3] = v3; r[4] = v4; r[5] = v5; break; }}
#define index6(cs, v0, v1, v2, v3, v4, v5)  \
			if ((vt[0] == v0) && (vt[1] == v1) && (vt[2] == v2) && (vt[3] == v3) && (vt[4] == v4) && (vt[5] == v5)) \
			    return cs

void alldata::getUnselected(char* v, int nv)
{
	int j = nv;
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

int alldata::getLastSixIndex(const char* resDay)
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
	getUnselected(v, 6);
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
		if (!m_bCheckLinkV || m_pCheckLink->checkLinks(links(), iDay))
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
