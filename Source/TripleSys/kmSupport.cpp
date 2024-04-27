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

void kmSortGroupsByFirstValue(char* mo, char* mi, char nr, char nc, char np)
{
	char* pmi[28];
	if (np > (sizeof(pmi) / sizeof(pmi[0])))
		abort();
	memset(pmi, 0, np * sizeof(pmi[0]));
	char* mic = mi;

	for (char ir = 0; ir < nr; ir++)
	{
		char iv = *mic;
		//if (iv >= nr || iv < 0 || pmi[iv] != NULL)
		//	abort();
		pmi[iv] = mic;
		mic += nc;
	}
	switch (nc)
	{
		case 2: {
			short int* mos = (short int*)mo;
			for (char i = 0; i < np; i++)
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
			for (char i = 0; i < np; i++)
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
			for (char i = 0; i < np; i++)
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
			for (char i = 0, j = 0; i < np; i++)
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
void kmSortRowsBy2ndValue(char* mo, char* mi, char nr, char nc, char* tm)
{
	char* pmi[28];
	if (nc > (sizeof(pmi) / sizeof(pmi[0])))
		abort();
	memset(pmi, 0, nc * sizeof(pmi[0]));
	char* mic = mi;

	for (char ir = 0; ir < nr; ir++)
	{
		char iv = *(mic + 1);
		//if (iv >= nr || iv < 0 || pmi[iv] != NULL)
		//	abort();
		tm[iv] = ir;
		pmi[iv] = mic;
		mic += nc;
	}
	char* moc = mo;
	for (char i = 0, j = 0; i < nc; i++)
	{
		mic = pmi[i];
		if (mic != NULL)
		{
			tm[j++] = tm[i];
			memcpy(moc, mic, nc);
			moc += nc;
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
int kmProcessMatrix(char* mo, char* mi, char* tmp, int nr, int nc, int gs, char* tr, int* pDayMax)
{
	if (tr != NULL)
		kmTranslate(mo, mi, tr, nr, nc);
	else 
		memcpy(mo, mi, nr * nc);
	if (gs == 3)
		kmSortGroups3(mo, nr, nc);
	else if (gs == 2)
	    kmSortGroups2(mo, nr, nc);
	else
		kmSortGroups(mo, nr, nc, gs);
	char* cii = mo;
	char* coi = tmp;
	for (int i = 0; i < nr; i++, coi += nc, cii += nc)
	{
		kmSortGroupsByFirstValue(coi, cii, nc / gs, gs, nc);
	}
	// result of the loop above is in tmp, sort and send it to mo
	char tm[32];
	if (nc > sizeof(tm))
		abort();
	kmSortRowsBy2ndValue(mo, tmp, nr, nc, tm);
	if (pDayMax == NULL)
		return memcmp(mo + nc, mi + nc, nc * nr - nc);
	int dayMax = tm[0];
	int icmp = 0;
	coi = mo;
	cii = mi;
	for (int i = 1; i < nr; i++)
	{   // start from 2nd row
		coi += nc;
		cii += nc;
		icmp = memcmp(coi, cii, nc);
		switch (icmp) 
	    {
			case -1: *pDayMax = dayMax > tm[i] ? dayMax : tm[i]; return -1;
			case 0: dayMax = dayMax > tm[i] ? dayMax : tm[i]; break;
			case 1: return 1;
		}
	}
	return 0;
}
int kmProcessOneRow2(char* mo, char* mi, int mind, char* ta, char* tb, char* tr, int nr, int nc)
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
	int iret;
	if (iRow == 2)
	{
		char diff = tb[1] - mi[nc + 3];
		if (diff != 0)
		{
			iret = diff > 0 ? 1 : -1;
			return iret;
		}
	}
	//if (iRow > nr)
	//	return 1;
	int ncc = (iRow - 1) * nc;
	char* mon = mo + ncc; // needed for macros below
	unsigned char m = 0;
	switch (nc)
	{
	case 10:
		SetTwoGroups2(0)
		//SetFirstTwoGroups2()
		SetTwoGroups2(2)
		SetTwoGroups2(4)
		SetTwoGroups2(6)
		SetTwoGroups2(8)
		break;
	case 12:
		SetTwoGroups2(0)
		//SetFirstTwoGroups2()
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
	if (iRow == 2)
		return memcmp(mo + ncc, mi + ncc, nc);
	return 0;
}
int kmProcessMatrix2(char* mo, char* mi, char* tr, int nr, int nc, char* rind, int ind, int* pDayMax)
{
	int iRet;
	char ta[16], tb[16], tm[16];
	bool bProc2 = false;
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
	memset(tm, unset, nc);
	tm[1] = ind;
	char rowMax = ind;
	char r2ind = rind[ind];
	if (r2ind >= nr)
		r2ind = 1;
	iRet = kmProcessOneRow2(mo, mi, r2ind, ta, tb, tr, nr, nc);
    char iRow = tb[0];
	if (iRow == 2)
	{
#if absstat
		b++;
#endif
		if (iRet != 0)
		{
			if (iRet < 0 && pDayMax != NULL)
				*pDayMax = std::max(r2ind, rowMax);
			return iRet;
		}
		bProc2 = true;
	}
	tm[iRow] = r2ind;

	for (int i = 0; i < nr; i++)
	{
		if (i == r2ind || i == ind)
			continue;
		iRet = kmProcessOneRow2(mo, mi, i, ta, tb, tr, nr, nc);
		iRow = tb[0];
		if (iRow == 2)
		{
#if absstat
			c++;
			if (rind[ind] == i)
				i = i;
#endif
			rind[ind] = i;
			if (iRet != 0)
			{
#if absstat
				//c++;
#endif
				if (iRet < 0 && pDayMax != NULL)
					*pDayMax = std::max((char)i, rowMax);
				return iRet;
			}
			bProc2 = true;
		}
		else if (iRow <= nr) // cant return if iRow > nr
		{
			if(bProc2 == true && iRow == 3)
			{
				iRet = memcmp(mo + nc * 2, mi + nc * 2, nc);
				if (iRet == -1)
				{
					if (pDayMax != NULL)
						*pDayMax = std::max(std::max((char)i, rowMax), tm[2]);
					return -1;
				}
				else if (iRet == 1)
					return 1;
			}
		}
		tm[iRow] = i;
	}
	if (!bProc2)
	{
#if absstat
		d++;
#endif
		return 1;
	}
	rowMax = std::max(rowMax, tm[2]);
	for (int i = 3; i <= nr; i++)
	{
		if (tm[i] == unset)
			return 1;
		iRet = memcmp(mo + nc * (i - 1), mi + nc * (i - 1), nc);
		switch (iRet)
		{
		case -1: if (pDayMax != NULL)*pDayMax = std::max(tm[i], rowMax); return -1;
		case 0: rowMax = std::max(tm[i], rowMax); break;
		case 1: return 1;
		}
	}
	return 0;
}
void kmTranslate(char* mo, char* mi, char* tr, int nr, int nc)
{
	int nall = nr * nc;
	for (int i = 0; i < nall; i++)
	{
		mo[i] = tr[mi[i]];
	}
}
