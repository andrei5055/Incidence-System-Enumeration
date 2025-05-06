#pragma once
#include "TripleSys.h"

#define sort3bytes(a, b, c) \
	if (a > b) SWAP(a, b); if (a > c) SWAP(a, c); if (b > c) SWAP(b, c)

#define U1F() \
	if ((k = t2[k1 = t1[k]]) == unset) \
		return i; \
	t2[k] = t2[k1] = unset

#define U1FN() \
	if ((k = t2[k1 = t1[k]]) != unset) \
		return i

#define U1FT(i) ro[ro[ri[i]] = ri[i + 1]] = ri[i]

//#define U1FT3(i) ro[ri[i]] = ri[i+1]; ro[ri[i+1]] = ri[i+2]; ro[ri[i+2]] = ri[i]
#define U1FT3(i) ro[ri[i]] = i; ro[ri[i+1]] = i+1; ro[ri[i+2]] = i+2

#define U1FCheck6(a, b, c, d, e) if (a != unset && \
	((a == b && a == c) || (a == d && a == e))) \
		goto nextPlayer6

#define P3Cycle1(m, a, idx) \
	for (int a = 0; a < 3; a++) \
	{ \
		const tchar tmp = t2d3[v[idx] = res1[m + a]]; \
		if (tmp >= m_nGroups || us[tmp] == unset) \
			continue; \
		us[tmp] = unset
#define P3Cycles3() \
	P3Cycle1(0, a, 0); \
	P3Cycle1(3, b, 1); \
	P3Cycle1(6, c, 2)
#define P3Cycles5() \
	P3Cycles3(); \
	P3Cycle1(9, d, 3); \
	P3Cycle1(12, e, 4)
#define P3Cycles7() \
	P3Cycles5(); \
	P3Cycle1(15, f, 5); \
	P3Cycle1(18, g, 6)
#define P3Cycles9() \
	P3Cycles7(); \
	P3Cycle1(21, h, 7); \
	P3Cycle1(24, o, 8)
#define P3Cycle1End() \
	us[tmp] = 0
#define P3Cycles3EndAndReturn(ret) \
	P3Cycle1End(); } \
	P3Cycle1End(); } \
	P3Cycle1End(); } \
	return ret
#define P3Cycles5EndAndReturn(ret) \
	P3Cycle1End(); } \
	P3Cycle1End(); } \
	P3Cycles3EndAndReturn(ret)
#define P3Cycles7EndAndReturn(ret) \
	P3Cycle1End(); } \
	P3Cycle1End(); } \
	P3Cycles5EndAndReturn(ret)
#define P3Cycles9EndAndReturn(ret) \
	P3Cycle1End(); } \
	P3Cycle1End(); } \
	P3Cycles7EndAndReturn(ret)
#define P3CyclesCheck(trc)	\
			ncycles = p3Cycles(trc, ncr, t1, t2, v, res1, res2); \
			if (cyclesNotOk(ncr, ncycles, trc->length)) \
				return ncycles ? -1: 0;


