#include <iostream>
#include <assert.h>
#include "TripleSys.h"

#define SetUpToTwoGroups2(i, m, j, ta, tb) \
		if (ta[i] != unset) \
		{m[j] = ta[i]; m[j + 1] = tb[i]; j += 2;} \
		if (ta[i + 1] != unset) \
		{m[j] = ta[i + 1]; m[j + 1] = tb[i + 1]; j += 2;}
#define SetUpToThreeGroups3(i, m, j, ta, tb, tc) \
		if (ta[i] != unset) \
		{m[j] = ta[i]; m[j + 1] = tb[i]; m[j + 2] = tc[i]; j += 3;} \
		if (ta[i + 1] != unset) \
		{m[j] = ta[i + 1]; m[j + 1] = tb[i + 1]; m[j + 2] = tc[i + 1]; j += 3;} \
		if (ta[i + 2] != unset) \
		{m[j] = ta[i + 2]; m[j + 1] = tb[i + 2]; m[j + 2] = tc[i + 2]; j += 3;}

#define sort3bytes(a, b, c) \
	if (a > b)std::swap(a, b); if (a > c)std::swap(a, c); if (b > c)std::swap(b, c)

#define sort3bytesInMemory(p) \
{ char a = *p, b = *(p + 1), c = *(p + 2); sort3bytes(a, b, c); *p = a; *(p + 1) = b; *(p + 2) = c; } 

#define SetTaTb(i, m, ta, tb, a, b) a=tr[m[i]]; b=tr[m[i+1]]; \
		if (a > b)std::swap(a, b); ta[a] = a; tb[a] = b

#define SetTaTbTc(i, m, ta, tb, tc, a, b, c) a=tr[m[i]]; b=tr[m[i+1]]; c=tr[m[i+2]]; \
		sort3bytes(a, b, c); \
		ta[a] = a; tb[a] = b; tc[a] = c

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
		sort3bytesInMemory(mii);
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
	//if (pDayMax == NULL)
	//	return memcmp(mo + nc, mi + nc, nc * nr - nc);
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
			case -1: if (pDayMax)*pDayMax = (dayMax > tm[i]) ? dayMax : tm[i]; return -1;
			case 0: dayMax = dayMax > tm[i] ? dayMax : tm[i]; break;
			case 1: if (coi[1] != cii[1]) return 2; return 1;
		}
	}
	return 0;
}
int kmProcessOneRow2(char* mo, char* mi, int mind, char* ta, char* tb, char* tr, int nr, int nc)
{
	char* mii = mi + mind * nc;
	memset(ta, unset, nc);
	char va, vb;
	switch (nc)
	{
	case 12:
		SetTaTb(0, mii, ta, tb, va, vb);
		SetTaTb(2, mii, ta, tb, va, vb);
		SetTaTb(4, mii, ta, tb, va, vb);
		SetTaTb(6, mii, ta, tb, va, vb);
		SetTaTb(8, mii, ta, tb, va, vb);
		SetTaTb(10, mii, ta, tb, va, vb);
		break;
	case 10:
		SetTaTb(0, mii, ta, tb, va, vb);
		SetTaTb(2, mii, ta, tb, va, vb);
		SetTaTb(4, mii, ta, tb, va, vb);
		SetTaTb(6, mii, ta, tb, va, vb);
		SetTaTb(8, mii, ta, tb, va, vb);
		break;
	default:
		for (int j = 0; j < nc; j += 2)
		{
			SetTaTb(j, mii, ta, tb, va, vb);
		}
	}
	char row2ndValue = tb[0];
	int iret;
	if (row2ndValue == 2)
	{
		char diff = tb[1] - mi[nc + 3];
		if (diff != 0)
		{
			iret = diff > 0 ? 1 : -1;
			return iret;
		}
	}
	//if (row2ndValue > nr)
	//	return 1;
	int ncc = (row2ndValue - 1) * nc;
	char* mon = mo + ncc; // needed for macros below
	unsigned char m = 0;
	switch (nc)
	{
	case 10:
		SetUpToTwoGroups2(0, mon, m, ta, tb)
		SetUpToTwoGroups2(2, mon, m, ta, tb)
		SetUpToTwoGroups2(4, mon, m, ta, tb)
		SetUpToTwoGroups2(6, mon, m, ta, tb)
		SetUpToTwoGroups2(8, mon, m, ta, tb)
		break;
	case 12:
		SetUpToTwoGroups2(0, mon, m, ta, tb)
		SetUpToTwoGroups2(2, mon, m, ta, tb)
		SetUpToTwoGroups2(4, mon, m, ta, tb)
		SetUpToTwoGroups2(6, mon, m, ta, tb)
		SetUpToTwoGroups2(8, mon, m, ta, tb)
		SetUpToTwoGroups2(10, mon, m, ta, tb)
		break;
	default:
		for (int j = 0; j < nc; j += 2)
		{
			SetUpToTwoGroups2(j, mon, m, ta, tb)
		}
	}
	if (row2ndValue == 2)
		return memcmp(mo + ncc, mi + ncc, nc);
	return 0;
}
int kmProcessOneRow3(char* mo, char* mi, int mind, char* ta, char* tb, char* tc, char* tr, int nr, int nc)
{
	char* mii = mi + mind * nc;
	memset(ta, unset, nc);
	char va, vb, vc;
	for (int j = 0; j < nc; j += 3)
	{
		SetTaTbTc(j, mii, ta, tb, tc, va, vb, vc);
	}
	char row2ndValue = tb[0];
	int iret;
	if (row2ndValue == 3)
	{
		char diff = (tc[0] == mi[nc + 2]) ? tb[1] - mi[nc + 4] : tc[0] - mi[nc + 2];
		if (diff != 0)
		{
			iret = diff > 0 ? 1 : -1;
			return iret;
		}
	}
	//if (row2ndValue > nr)
	//	return 1;
	int ncc = (row2ndValue - 1) * nc;
	char* mon = mo + ncc; // needed for macros below
	unsigned char m = 0;
	for (int j = 0; j < nc - 2; j += 3)
	{
		SetUpToThreeGroups3(j, mon, m, ta, tb, tc)
	}
	if (row2ndValue == 3)
		return memcmp(mo + ncc, mi + nc, nc);
	return 0;
}
int kmProcessMatrix2(char* mo, char* mi, char* tr, int nr, int nc, char* rind, int ind, int* pDayMax)
{
	int iRet;
	char ta[16], tb[16], tm[16];
	bool bProc2 = false;
	if (nc > sizeof(ta))
		abort();
	tb[0] = 1; // to make compiler happy
	memset(tm, unset, nc);
	tm[1] = ind;
	char rowMax = ind;
	char r2ind = rind[ind];
	if (r2ind >= nr)
		r2ind = 1;
	iRet = kmProcessOneRow2(mo, mi, r2ind, ta, tb, tr, nr, nc);
	char row2ndValue = tb[0];
#if 0
	Stat("1tr", 10, true);
	Stat("2nd row first", 11, row2ndValue == 2);
	Stat("fast ret", 12, row2ndValue == 2 && iRet != 0);
	/*
	without bitmask and tr calculation only for last day:
		iDay=9: kmProcessMatrix2:3825958 2nd row first:3138093 fast returns:3104670
		Thread 1: 396 non-isomorphic matrices (10,9,2) created
		Thread execution time = 8831 ms

	with "bitmask" and tr calculations for last 2 days:
		iDay=9: kmProcessMatrix2:3967449 2nd row first:2682637 fast returns:2651996
		Thread 1: 396 non-isomorphic matrices (10,9,2) created
		Thread execution time = 13148 ms
	*/
#endif
	if (row2ndValue == 2)
	{
		if (iRet != 0)
		{
			if (iRet < 0 && pDayMax != NULL)
				*pDayMax = std::max(r2ind, rowMax);
			return iRet;
		}
		bProc2 = true;
	}
	if (row2ndValue < 0 || row2ndValue >= nc)
		abort();
	tm[row2ndValue] = r2ind;

	for (int i = 0; i < nr; i++)
	{
		if (i == r2ind || i == ind)
			continue;
		iRet = kmProcessOneRow2(mo, mi, i, ta, tb, tr, nr, nc);
		row2ndValue = tb[0];
		if (row2ndValue == 2)
		{
			rind[ind] = i;
			if (iRet != 0)
			{
				if (iRet < 0 && pDayMax != NULL)
					*pDayMax = std::max((char)i, rowMax);
				return iRet;
			}
			bProc2 = true;
		}
		else if (row2ndValue <= nr) // cant return if row2ndValue > nr
		{
			if (bProc2 == true && row2ndValue == 3)
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
		tm[row2ndValue] = i;
	}
	if (!bProc2)
	{
		return 2;
	}
	rowMax = std::max(rowMax, tm[2]);
	for (int i = 3; i <= nr; i++)
	{
		if (tm[i] == unset)
			return 2;
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
int kmProcessMatrix3(char* mo, char* mi, char* tr, int nr, int nc, char* rind, int ind, int* pDayMax)
{
	int iRet;
	char ta[33], tb[33], tc[33], tm[33];
	bool bProc2 = false;
	if (nc > sizeof(ta))
		abort();
	tb[0] = 1; // to make compiler happy
	memset(tm, unset, nc);
	tm[1] = ind;
	char rowMax = ind;
	char r2ind = rind[ind];
	if (r2ind >= nr)
		r2ind = 1;
	iRet = kmProcessOneRow3(mo, mi, r2ind, ta, tb, tc, tr, nr, nc);
	char row2ndValue = tb[0];
#if 0
	Stat("kmProcessMatrix3", 0, true);
	Stat("2nd row first", 1, row2ndValue == 3);
	Stat("fast returns", 2, row2ndValue == 3 && iRet != 0);
#endif
	if (row2ndValue == 3) // 2nd row
	{
		if (iRet != 0)
		{
			if (iRet < 0 && pDayMax != NULL)
				*pDayMax = std::max(r2ind, rowMax);
			return iRet;
		}
		bProc2 = true;
	}
	tm[row2ndValue] = r2ind;

	for (int i = 0; i < nr; i++)
	{
		if (i == r2ind || i == ind)
			continue;
		iRet = kmProcessOneRow3(mo, mi, i, ta, tb, tc, tr, nr, nc);
		row2ndValue = tb[0];
		if (row2ndValue < 0 || row2ndValue >= nc)
			abort();
		if (row2ndValue == 3) // 2nd row
		{
			rind[ind] = i; // update index for next time
			if (iRet != 0)
			{
				if (iRet < 0 && pDayMax != NULL)
					*pDayMax = std::max((char)i, rowMax);
				return iRet;
			}
			bProc2 = true;
		}
		else if (row2ndValue <= nr) // check for row 3
		{
			if (bProc2 == true && row2ndValue == 4) // 3rd row and 2nd row was processed cmp was equal 0
			{
				iRet = memcmp(mo + nc * 3, mi + nc * 2, nc);
				if (iRet == -1)
				{
					if (pDayMax != NULL)
						*pDayMax = std::max(std::max((char)i, rowMax), tm[3]);
					return -1;
				}
				else if (iRet == 1)
					return 1;
			}
		}
		tm[row2ndValue] = i;
	}
	if (!bProc2)
	{
		return 2;
	}
	rowMax = std::max(rowMax, tm[3]);
	char* moi = mo + nc * 3;
	char* mii = mi + nc * 2;
	// following loop checs all rows starting from row 3 
	for (int i = 4, j = 2; i < nc - 1 && j < nr; i++, moi+= nc)
	{
		if (tm[i] == unset)
			continue;
		if (*(moi + 1) != *(mii + 1))
			return 2;
		iRet = memcmp(moi, mii, nc);
		switch (iRet)
		{
		case -1: if (pDayMax != NULL)*pDayMax = std::max(tm[i], rowMax); return -1;
		case 0: rowMax = std::max(tm[i], rowMax); break;
		case 1: return 1;
		}
		j++;
		mii += nc;
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
