#include <iostream>
#include <assert.h>
#include "TripleSys.h"

void kmSortAllRowsWith2PlayersInGroup(char* mo, char* mi, int nr, int nc)
{
	char* r1 = mi;
	for (int i = 0; i < nr; i++, r1+=nc)
	{
		int ip1 = r1[1] - 1;
		if (ip1 < 0 || ip1 >= nr)
			abort();
		char* r2 = mo + ip1 * nc;
		int ip2 = r2[1];
		memcpy(r2, r1, nc);
	}
}
void kmSortOneRowWith2PlayersInGroup(char* r, int nr, int nc)
{
	for (int i = 1; i < nr; i++)
	{
		bool ex = true;
		char* r1 = r;
		char* r2 = r1 + nc;
		for (int j = 0; j < nr - i; j++, r1 += nc, r2 += nc)
		{
			if (*r1 > *r2)
			{
				ex = false;
				register short int t = *((short int*)r1);
				*((short int*)r1) = *((short int*)r2);
				*((short int*)r2) = t;
			}
		}
		if (ex)
			break;
	}
}
void kmSortRowsByOneValue(char* r, int nr, int nc, int ip)
{
	char t[33];
	if (nc > sizeof(t))
		abort();
	for (int i = 1; i < nr; i++)
	{
		bool ex = true;
		char* r1 = r;
		char* r2 = r1 + nc;
		for (int j = 0; j < nr - i; j++, r1 += nc, r2 += nc)
		{
			if (*(r1 + ip) > *(r2 + ip))
			{
				ex = false;
				memcpy(t, r1, nc);
				memcpy(r1, r2, nc);
				memcpy(r2, t, nc);
			}
		}
		if (ex)
			break;
	}
}

void kmFullSort2(char* mo, char* mi, int nr, int nc)
{
	int nrnc = nr * nc;
	for (int i = 0; i < nrnc; i += 2)
	{
		char t = mi[i];
		if (t > mi[i + 1])
		{
			mi[i] = mi[i + 1]; mi[i + 1] = t;
		}
	}
	for (int i = 0; i < nr; i++)
	{
		char* km = mi + i * nc;
		kmSortOneRowWith2PlayersInGroup(km, nc / 2, 2);
	}
	kmSortAllRowsWith2PlayersInGroup(mo, mi, nr, nc);
}
void kmFullSort(char* mi, int nr, int nc, int gs)
{
	int nrnc = nr * nc;
	for (int i = 0; i < nrnc; i += gs)
	{
		char* mii = mi + i;
		for (int j = 1; j < gs; j++)
		{
			bool ex = true;
			for (int k = 0; k < gs - j; k++)
			{
				if (mii[k] > mii[k + 1])
				{
					ex = false;
					char t = mii[k]; mii[k] = mii[k + 1]; mii[k + 1] = t;
				}
			}
			if (ex)
				break;
		}
	}
	for (int i = 0; i < nr; i++)
	{
		char* km = mi + i * nc;
		kmSortRowsByOneValue(km, nc / gs, gs, 0);
	}
	kmSortRowsByOneValue(mi, nr, nc, 1);
}
void kmTranslate(char* mo, char* mi, char* tr, int nr, int nc)
{
	if (tr[0] < 0)
	printTable("2", tr, 1, nc, 0, 2);
	for (int i = 0; i < nr * nc; i++)
	{
		mo[i] = tr[mi[i]];
	}
}
