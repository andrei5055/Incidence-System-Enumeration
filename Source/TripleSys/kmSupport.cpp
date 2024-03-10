#include <iostream>
#include <assert.h>
#include "TripleSys.h"

void kmSortRowsByOneValue(char* r, int nr, int nc, int ip)
{
	char t[33];
	if (nc > sizeof(t))
		abort();
	for (int i = 1; i < nr; i++)
	{
		char* r1 = r;
		char* r2 = r1 + nc;
		for (int j = 0; j < nr - i; j++, r1+=nc, r2+=nc)
		{
			if (*(r1+ip) > *(r2+ip))
			{
				memcpy(t, r1, nc);
				memcpy(r1, r2, nc);
				memcpy(r2, t, nc);
			}
		}
	}
}

void kmFullSort(char* mo, int nr, int nc, int gs)
{
	int nrnc = nr * nc;
	if (gs == 2)
	{
		for (int i = 0; i < nrnc; i += gs)
		{
			char t = mo[i];
			if (t > mo[i + 1])
			{
				mo[i] = mo[i + 1]; mo[i + 1] = t;
			}
		}
	}
	else
	{
		for (int i = 0; i < nrnc; i += gs)
		{
			char* moi = mo + i;
			for (int j = 1; j < gs; j++)
			{
				for (int k = 0; k < gs - j; k++)
				{
					if (moi[k] > moi[k + 1])
					{
						char t = moi[k]; moi[k] = moi[k + 1]; moi[k + 1] = t;
					}
				}
			}
		}
	}

	for (int i = 0; i < nr; i++)
	{
		char* km = mo + i * nc;
		kmSortRowsByOneValue(km, nc / gs, gs, 0);
	}
	kmSortRowsByOneValue(mo, nr, nc, 1);
}
void kmTranslate(char* mo, char* mi, char* tr, int nr, int nc)
{
	for (int i = 0; i < nr * nc; i++)
	{
		mo[i] = tr[mi[i]];
	}
}
