#include <iostream>
#include <assert.h>
#include "TripleSys.h"

#define absstat 0

#define SetTaTb(_ind) va=tr[mii[_ind]]; vb=tr[mii[_ind+1]]; \
		if (va > vb){ta[vb] = vb;tb[vb] = va;} \
			else{ta[va] = va;tb[va] = vb;}
#define SetTwoGroups2(_ind) \
		if (ta[_ind] != unset) \
		{mon[m] = ta[_ind]; mon[m + 1] = tb[_ind]; m += 2;} \
		if (ta[_ind + 1] != unset) \
		{mon[m] = ta[_ind + 1]; mon[m + 1] = tb[_ind + 1]; m += 2;}
#define SetFirstTwoGroups2() \
		mon[m] = ta[0]; mon[m + 1] = tb[0]; \
		mon[m + 2] = ta[1]; mon[m + 3] = tb[1]; m += 4;

void kmSortRowsByOneValue(char* mo, char* mi, char np, char nr, char nc, char ip)
{
	char* pmi[28];
	if (np > (sizeof(pmi)/sizeof(pmi[0])))
		abort();
	char n = np;
	memset(pmi, 0, n * sizeof(pmi[0]));
	char* mic = mi;
	char nrows = nr;
	while (nrows-- > 0)
	{
		char iv = *(mic + ip);
		//if (iv >= n || iv < 0 || pmi[iv] != NULL)
		//	abort();
		pmi[iv] = mic;
		mic += nc;
	}
	switch (nc)
	{
		case 2: {
			short int* mos = (short int*)mo;
			for (char i = 0; i < n; i++)
			{
				short int* mis = (short int*)(pmi[i]);
				if (mis != NULL)
				{
					*mos = *mis;
					mos++;
				}
			}
			break;
		}
		case 3: {
			char* moc = mo;
			for (char i = 0; i < n; i++)
			{
				mic = pmi[i];
				if (mic != NULL)
				{
					*moc = *mic;
					*(moc + 1) = *(mic + 1);
					*(moc + 2) = *(mic + 2);
					moc += 3;
				}
			}
			break;
		}
		case 4: {
			int* moi = (int*)mo;
			for (char i = 0; i < n; i++)
			{
				int* mii = (int*)(pmi[i]);
				if (mii != NULL)
				{
					*moi = *mii;
					moi++;
				}
			}
			break;
		}
		default: {
			char* moc = mo;
			for (char i = 0; i < n; i++)
			{
				mic = pmi[i];
				if (mic != NULL)
				{
					memcpy(moc, mic, nc);
					moc += nc;
				}
			}
		}
	}
}
void kmSortGroups2(char* mi, int nr, int nc)
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
}
void kmSortGroups(char* mi, int nr, int nc, int gs)
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
}
void kmSortGroups3(char* mi, int nr, int nc)
{
	char ng = nr * nc / 3;
	char* mii = mi;
	for (int i = 0; i < ng; i++, mii += 3)
	{
		char a = *mii, b = *(mii + 1), c = *(mii + 2);
		if (a > b) {
			if (b > c) { *mii = c; *(mii + 2) = a; }                  // c b a 
			else if (a < c) { *mii = b; *(mii + 1) = a; }             // b a c
			else { *mii = b; *(mii + 1) = c; *(mii + 2) = a; }        // b c a
		} else { 
			if (a > c) { *mii = c; *(mii + 1) = a; *(mii + 2) = b; }  // c a b
			else if (b > c) { *(mii + 1) = c; *(mii + 2) = b; }       // a c b
		}
	}
}
void kmFullSort(char* tmp, char* mi, int nr, int nc, int gs)
{
	if (gs == 3)
		kmSortGroups3(mi, nr, nc);
	else if (gs == 2)
	    kmSortGroups2(mi, nr, nc);
	else
		kmSortGroups(mi, nr, nc, gs);
	char* mii = mi;
	char* moi = tmp;
	for (int i = 0; i < nr; i++, moi += nc, mii += nc)
	{
		kmSortRowsByOneValue(moi, mii, nc, nc / gs, gs, 0);
	}
	// result of the loop above is in tmp, send it back to mi
	kmSortRowsByOneValue(mi, tmp, nc, nr, nc, 1);
}
int kmTranslateAndCheckOneRow2(char* mo, char* mi, int mind, char* ta, char* tb, char* tr, int nr, int nc, int ir)
{
	char* mii = mi + mind * nc;
	memset(ta, unset, nc);
	char va, vb;
	switch (nc)if (nc == 12)
	{
	case 12:
		SetTaTb(0);
		SetTaTb(2);
		SetTaTb(4);
		SetTaTb(6);
		SetTaTb(8);
		SetTaTb(10);
		break;
	case 10:
		SetTaTb(0);
		SetTaTb(2);
		SetTaTb(4);
		SetTaTb(6);
		SetTaTb(8);
		break;
	default:
		for (int j = 0; j < nc; j += 2)
		{
			SetTaTb(j);
		}
	}
	char iRow = tb[0];
	int iret = 3;
#if 1
	if (iRow == 2)
	{
		char diff = tb[1] - mi[nc + 3];
		if (diff != 0 && ir < 2)
		{
			iret = diff > 0 ? 1 : -1;
			return iret;
		}
		iret = 2;
	}
#endif

	char* mon = mo + (iRow - 1) * nc;
	unsigned char m = 0;
	switch (nc)
	{
	case 10:
		SetTwoGroups2(0)
		SetTwoGroups2(2)
		SetTwoGroups2(4)
		SetTwoGroups2(6)
		SetTwoGroups2(8)
		break;
	case 12:
		SetTwoGroups2(0)
		SetTwoGroups2(2)
		SetTwoGroups2(4)
		SetTwoGroups2(6)
		SetTwoGroups2(8)
		SetTwoGroups2(10)
		break;
	default:
		for (int j = 0; j < nc; j+=2)
		{
			SetTwoGroups2(j)
		}
	}
	return iret;
}
int kmTranslateAndSort2(char* mo, char* mi, char* tr, int nr, int nc, int ir, char* rind, int ind)
{
	int nall = nr * nc;
	int iRet;
	char ta[16*2], tb[16*2];
	bool bProc2 = false, bProc3 = false;
	if (nc > 16)
		abort();
	tb[0] = 1; // to make compiler happy
#if absstat
	static int a, b, c, d, e, z[13];
	a++;
	if ((a % 1000000) == 0)
	{
		printf("%d fast=%d second=%d no2=%d left=%d", a, b, c, d, e);
		for (int i = 0; i < nc - 1; i++)
			printf(" %d:%d", i, z[i]);
		printf("\n");
	}
#endif
	if (nr < nc - 1)
	{
		for (int i = nc; i < nc * (nc - 1); i += nc)
			mo[i] = unset;
	}
	char r2ind = rind[ind];
	iRet = 3;
	if (r2ind < nr)
	    iRet = kmTranslateAndCheckOneRow2(mo, mi, r2ind, ta, tb, tr, nr, nc, ir);
	//if (nr < nc - 1)
	//	return iRet < 0 ? -1 : 1;
	if (iRet < 3)
	{
		//printf("\n 2\n");
#if absstat
		b++;
#endif
		if (iRet < 2)
		    return iRet;
		iRet = memcmp(mo + nc, mi + nc, nc);
		if (iRet != 0 && ir < 2)
			return iRet;
		bProc2 = true;
	}
	for (int i = 0; i < nr; i++)
	{
		iRet = kmTranslateAndCheckOneRow2(mo, mi, i, ta, tb, tr, nr, nc, ir);
		char iRow = tb[0];
		if (bProc2 && iRow == 2)
			continue;
		/**
		if (i == 0)
			printf("\n%d ", nr);
		printf(" %2d:%2d", i, iRow);
		if (i == nr - 1)
			printf("\n");
		if (iRow < 1)
			abort();
		**/
		if (iRow == 2)
		{
#if absstat
			c++;
			if (rind[ind] == i && iRet < 2)
				i = i;
#endif
			rind[ind] = i;
			if (iRet < 2)
			{
#if absstat
				//c++;
#endif
				return iRet;
			}
		}
#if 1
		if (iRow <= nr)
		{
			switch (iRow) {
				case 2: bProc2 = true; iRet = memcmp(mo + nc, mi + nc, bProc3 ? nc * 2 : nc); break;
				case 3: bProc3 = true;
					if (bProc2)
					{
						iRet = memcmp(mo + nc, mi + nc, nc * 2);
						break;
					}
					else
						continue;
				default: continue;
			}
			if (iRet == 1 || (iRet != 0 && ir < 2))
			{
#if absstat
				z[i]++;
#endif
				return iRet;
			}
		}	
#endif
	}
	if (!bProc2 && ir < 2)
	{
#if absstat
		d++;
#endif
		return 1;
	}
	if (nr < nc - 1)
	{
		int k = nc;
		int ncncm1 = nc * (nc - 1);
		if (ir < 2)
		{
			for (int i = nc; i < ncncm1; i += nc)
			{
				if (mo[i] != unset)
				{
					iRet = memcmp(mo + k, mo + i, nc);
					if (iRet != 0)
						return iRet;
					k += nc;
				}
			}
			return 0;
		}
		for (int i = nc; i < ncncm1; i += nc)
		{
			if (mo[i] != unset)
			{
				if (i != k)
					memcpy(mo + k, mo + i, nc);
				k += nc;
			}
		}
	}
#if 0
	for (int i = 1; i < nr; i++)
	{
		int ic = memcmp(mo + nc*i, mi + nc*i, nc);
		if (ic != 0)
		{
			//if (mo[nc*i + 1] == i+1)
			   z[i]++;
			return ic;
			
		}
	}
	return 0;
#endif
#if absstat
	e++;
#endif
	return memcmp(mo + nc, mi + nc, nall - nc);
	//z[nr - 1]++;
	return 1;
}
void kmTranslate(char* mo, char* mi, char* tr, int nr, int nc)
{
	int nall = nr * nc;
	for (int i = 0; i < nall; i++)
	{
		mo[i] = tr[mi[i]];
	}
}
