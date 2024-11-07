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
#define P3CyclesCheck()	\
			ncycles = p3Cycles(ncr, t1, t2, v, res1, res2); \
			if (cyclesNotOk(ncr, ncycles, m_TrCycles.length)) \
				return ncycles;

// expected second row
static tchar _expected2ndRow3uf_21[] = {
//1: | Aut(M) | = 1728, Cycles:18(6:6:9) 36(9:12)
   0,3,6, 1,4,7, 2,5,8, 9,12,15,10,13,18,11,16,19,14,17,20,
//2:| Aut(M) | = 64, Cycles:16(9:12) 16(21)
   0,3,6, 1,4,7, 2,5,9, 8,10,12,11,15,18,13,16,19,14,17,20,
//3:| Aut(M) | = 16, Cycles:2(6:6:9) 8(6:15) 4(9:12) 16(21)
   0,3,6, 1,4,7, 2,5,9, 8,12,15,10,13,18,11,16,19,14,17,20,
//4:| Aut(M) | = 32, Cycles:24(21)
   0,3,6, 1,4,7, 2,9,12, 5,10,13, 8,15,18,11,16,19,14,17,20,
//5:| Aut(M) | = 48, Cycles:24(21)
   0,3,6, 1,4,7, 2,9,12, 5,10,15, 8,11,18,13,16,19,14,17,20,
//6:| Aut(M) | = 4, Cycles:6(6:15) 4(9:12) 16(21)
   0,3,6, 1,4,7, 2,9,12, 5,10,15, 8,13,18,11,16,19,14,17,20,
//7:| Aut(M) | = 12, Cycles:3(6:6:9) 12(6:15) 12(21)
   0,3,6, 1,4,9, 2,7,10, 5,12,15, 8,13,18,11,16,19,14,17,20,
//8:| Aut(M) | = 28, Cycles:7(6:6:9) 7(6:15) 7(9:12) 10(21)
   0,3,6, 1,4,9, 2,7,12, 5,10,15, 8,13,18,11,16,19,14,17,20,
//9:| Aut(M) | = 4, Cycles:1(6:6:9) 7(6:15) 5(9:12) 14(21)
   0,3,6, 1,4,9, 2,7,12, 5,10,15, 8,16,18,11,13,19,14,17,20,
//10:| Aut(M) | = 4, Cycles:2(6:6:9) 4(6:15) 6(9:12) 14(21)
   0,3,6, 1,4,9, 2,7,12, 5,13,15, 8,16,18,10,14,19,11,17,20,
//11:| Aut(M) | = 8, Cycles:1(6:6:9) 4(6:15) 4(9:12) 16(21)
   0,3,6, 1,4,9, 2,7,12, 5,15,18, 8,16,19,10,13,17,11,14,20,
//12:| Aut(M) | = 336, Cycles:24(21)
   0,3,6, 1,9,12, 2,15,18, 4,10,16, 5,13,19, 7,11,20, 8,14,17,
   };

