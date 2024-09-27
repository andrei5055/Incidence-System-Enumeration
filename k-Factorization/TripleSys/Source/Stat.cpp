#include "TripleSys.h"

static sLongLong __a[16] = { 0 };
static const char* __t[16] = { 0 };
void StatAdd(const char* t, int ind, bool bAdd)
{
	if (ind < 0 || ind >= countof(__t))
		abort();
	__t[ind] = t;
	if (bAdd)__a[ind]++;
}
void StatEnd(bool bReset, const char* cHdr, int iHdr, bool bPrint)
{
	if (bPrint)
	{
		bool noStat = true;
		for (int i = 0; i < countof(__t); i++)
		{
			if (__t[i])
				noStat = false;
		}
		printf(GreenText);
		printf("%s=%02d", cHdr, iHdr);
		if (!noStat)
		{
			for (int i = 0; i < countof(__t); i++)
			{
				if (__t[i])
				{
					noStat = false;
					printf("|%8zd:%s", __a[i], __t[i]);
				}
			}
			if (bReset)
			{
				memset(__a, 0, sizeof(__a));
				memset(__t, 0, sizeof(__t));
			}
		}
		printf("|\n");
		printf(ResetTextColor);
	}
}
