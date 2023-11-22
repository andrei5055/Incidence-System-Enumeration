
#include <iostream>
#include "TripleSys.h"
void printResult(alldata* s)
{
	for (int j = 0; j < nDays; j++)
	{
		for (int i = 0; i < nPlayers; i = i + GroupSize)
		{
			for (int k = 0; k < GroupSize; k++)
				printf(" %2d", (int)s->result[j][i + k]);
			printf(" ");
		}
		printf("\n");
	}
	printf("\n");
}
void init(alldata* s, const char** initdata)
{
	char p[nPlayers];
	for (s->iDay = 0; s->iDay < nDays; s->iDay++)
	{
		const char* data = initdata[s->iDay];
		int n = 0;
		if (data[0] == '\0')
			break;
		int check = 0;
		for (int i = 0; n < nPlayers; i++)
		{
			if (data[i] == '\0')
				break;

			if (data[i] >= '0' && data[i] <= '9')
				p[n] = data[i] - '1';
			else if (data[i] >= 'a' && data[i] <= 'f')
				p[n] = data[i] - 'a' + 9;
			else if (data[i] >= 'A' && data[i] <= 'F')
				p[n] = data[i] - 'A' + 9;
			else
				continue;
			check = check | (1 << p[n]);
			if (!setLinksAndDevCounts(s, p, n, 1))
			{
				printf("Error in 'setLinksAndDevCounts' iDay=%d pos=%d player=%d d=%s\n",
					s->iDay, n, p[n], data);
			}
			n++;
		}
		if (check != (1 << nPlayers) - 1)
		{
			printf("Error in 'setlinks' iDay=%d x=%8.8x n=%d d=%s\n", s->iDay, check, n, data);
		}
		memcpy(s->result[s->iDay], p, sizeof(p));
		//memcpy(&s->bklinks[0], &s->links, sizeof(s->links));
		//memcpy(&s->bklinks[1], &s->links, sizeof(s->links));
	}
}
char adjust(int id, int ip)
{
	if (ip == 0)
		return 0;
	else if (ip < GroupSize)
		return id * (GroupSize - 1) + ip;
	else if (ip == GroupSize)
		return ((id > 0) ? 1 : GroupSize);

	return unset;
}

void adjust0(alldata* s)
{
	for (int m = 0; m < nDays; m++)
	{
		for (int i = 0; i < nPlayers; i++)
		{
			if ((s->preSet[m][i] = adjust(m, i)) != unset)
			{
				if (!setLinksAndDevCounts(s, s->preSet[m], i, 1))
					abort();
			}
		}
	}
}
bool setLinksAndDevCounts(alldata* s, char* p, int ip, char iset)
{
	char tmp[(nPlayers - 1) / 2 + 1];
	char tmp2[GroupSize];
	int i = ip % GroupSize;
	if (i == 0)
		return true;
	//	if (iset != 1 && p[ip] < 0)
	//		return true;
	char bset = iset == 1 ? 1 : 0;
	char i1 = p[ip];
	//	if (i1 < 0 || i1 >= nPlayers)
	//		i1 = i1;
	memset(tmp, 0, sizeof(tmp));
	for (int j = 1; j <= i; j++)
	{
		char i2 = p[ip - j];
		//		if (i2 < 0 || i2 >= nPlayers)
		//			i2 = i2;
		if (s->links[i1][i2] == bset)
			return false;
		char idev = i1 < i2 ? i2 - i1 : i1 - i2;
		if (idev > (nPlayers / 2))
			idev = nPlayers - idev;
		//		if (idev < 0 || idev > nPlayers/2+1)
		//			idev = 0;
		tmp[idev] += iset;
		char count = s->devCounts[idev] + tmp[idev];
		if (iset == 1 && count > nPlayers)
			return false;
		else if (iset != 1 && count < 0)
			return false;
		tmp2[j] = idev;

	}
	for (int j = 1; j <= i; j++)
	{
		char i2 = p[ip - j];
		s->links[i1][i2] = bset;
		s->links[i2][i1] = bset;
		s->devCounts[tmp2[j]] += iset;
	}
	return true;
}
void checkbmask(alldata *s)
{
	for (int i = 0; i < nPlayers; i++)
	{
		if (s->selPlayers[i] != unset && (s->bmask & (1 << i)) != 0)
			i = i;
		else if (s->selPlayers[i] == unset && (s->bmask & (1 << i)) == 0)
			i = i;
	}
}
void initPrevDayProcess(alldata *s)
{
	s->iPlayer = nPlayers - GroupSize;
	s->maxPlayers = s->iPlayer;
	memcpy(s->indexPlayer, s->result[s->iDay], sizeof(s->indexPlayer));
	memcpy(s->tmpPlayers, s->result[s->iDay], sizeof(s->tmpPlayers));
	memset(s->selPlayers, 1, sizeof(s->selPlayers));
	s->bmask = 0;
	for (int j = s->iPlayer; j < nPlayers; j++)
	{
		if (!setLinksAndDevCounts(s, s->result[s->iDay], j, unset))
			abort();
		int k = s->tmpPlayers[j];
		s->bmask = s->bmask | (1 << k);
		s->selPlayers[k] = unset;
		s->tmpPlayers[j] = unset;
	}
	//checkbmask(s);
	s->indexPlayer[s->iPlayer]++;
	// temporary


	//memset(s->indexPlayer, 0, sizeof(s->indexPlayer));
	//s->indexPlayer[4] = s->result[s->iDay][4];

}
void initCurrentDayProcess(alldata* s)
{
	s->iPlayer = 0;
	s->maxPlayers = 0;
	memset(s->indexPlayer, 0, sizeof(s->indexPlayer));
	memset(s->selPlayers, unset, sizeof(s->selPlayers));
	memset(s->tmpPlayers, unset, sizeof(s->tmpPlayers));
	s->bmask = (1 << nPlayers) - 1;
	for (int j = 0; j < nPlayers; j++)
	{
		int k = s->preSet[s->iDay][j];
		if (k != unset)
		{
			s->selPlayers[k] = j;
			s->bmask = s->bmask ^ (1 << k);
		}
	}
	//checkbmask(s->selPlayers, s->bmask);

}