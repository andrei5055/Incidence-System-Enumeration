#include <iostream>
#include "TripleSys.h"
const char* ms =
//Matrix #   5 | Aut(M) | = 1
"11111110000000000000000000000000000"
"10000001111110000000000000000000000"
"10000000000001111110000000000000000"
"01000001000001000001111000000000000"
"01000000100000100000000111100000000"
"00100001000000100000000000011110000"
"00100000010000010001000100000001100"
"00010000100000001000100000010001010"
"00010000010000000100010010001000001"
"00001000001001000000000010000100110"
"00001000000100010000001001010000001"
"00000100001000000010010001000011000"
"00000100000010001001000000100100001"
"00000010000100000010100000101000100"
"00000010000010000100001100000010010"
;
bool s2k(const char* s, char *lnk, int nr, int nc)
{
	char* m = s == NULL ? (char*)ms : (char*)s;
	memset(lnk, unset, nr * nr);
	for (int i = 0; i < nc; i++)
	{
		int n = 0;
		char t[3];
		for (int j = 0; j < nr; j++)
		{
			char v = m[j * nc + i];
			if (v == '1')
			{
				if (n == 3)
				{
					printf("Error: more than 3 values in group (column) %d\n", i);
					return false;
				}
				t[n] = j;
				n++;
			}
		}

		if (n != 3)
		{
			printf("Error: less than 3 (%d) values in group (column) %d\n", n, i);
			return false;
		}
		if (lnk[t[0] * nr + t[1]] != unset ||
			lnk[t[0] * nr + t[2]] != unset ||
			lnk[t[1] * nr + t[2]] != unset)
		{
			printf("Error: one of the pair in group (column) %d (%d,%d,%d) already defined\n",
				i, t[0], t[1], t[2]);
			return false;
		}
		int d = i % ((nr - 1) / 2);
		lnk[t[0] * nr + t[1]] = d;
		lnk[t[0] * nr + t[2]] = d;
		lnk[t[1] * nr + t[2]] = d;
		lnk[t[1] * nr + t[0]] = d;
		lnk[t[2] * nr + t[0]] = d;
		lnk[t[2] * nr + t[1]] = d;
	}
	char* ln = lnk;
	for (int i = 0; i < nr; i++, ln += nr)
	{
		int n = 0;
		char t[100];
		memset(t, 0, nr);
		for (int j = 0; j < nr; j++)
		{
			char v = ln[j];
			if (v != unset)
			{
				t[v]++;
				if (t[v] > 2)
				{
					printf("Error: player %d present more than 1 time in day %d\n", i, v);
					return false;
				}
				if (t[v] == 2)
					n++;
			}
		}
		if (n != (nr - 1) / 2)
		{
			printf("Error: player %d missed one or more day(s)\n", i);
			return false;
		}
	}
	return true;
}