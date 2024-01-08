#include <iostream>
#include "TripleSys.h"

char lastError[256];

bool checkMatrix1(char* lnks, int ln, int nc, const char* matrix, int i1, int i2, bool printError, int* errLine, int* errGroup, int* dubLine)
{
	const char* c = matrix + ln * nc;
	char a = c[i1], b = c[i2];
	if (a < 0 || a >= nc ||
		b < 0 || b >= nc ||
		a == b)
	{
		printf("CheckMatrix: matrix values (%d %d) in group %d on line %d are incorrect, job aborted\n", a, b, i1 / 3, ln);
		abort();
	}
	char dLine = lnks[a * nc + b];
	if (dLine != unset)
	{
		sprintf_s(lastError, "CheckMatrix: matrix pair (%d %d) in group %d on line %d already defined in line %d\n", a, b, i1 / 3, ln, dLine);
		if (printError)
			printf(lastError);

		*errLine = ln;
		*errGroup = i1 / 3;
		*dubLine = dLine;
		return false;
	}
	lnks[a * nc + b] = lnks[b * nc + a] = ln;
	return true;
}
bool _CheckMatrix(const char* matrix, int nl, int nc, char *lnks, bool printError, int* errLine, int* errGroup, int* dubLine)
{
	if (nl < 0 || nl > 10 || (nc != 15 && nc != 21))
	{
		printf("CheckMatrix: incorrect parameters (nlines=%d ncolumns=%d), job aborted\n", nl, nc);
		abort();
	}

	memset(lnks, unset, nc * nc);

	for (int j = 0; j < nl; j++)
	{
		for (int i = 0; i < nc; i = i + 3)
		{
			if (!checkMatrix1(lnks, j, nc, matrix, i,     i + 1, printError, errLine, errGroup, dubLine) ||
				!checkMatrix1(lnks, j, nc, matrix, i,     i + 2, printError, errLine, errGroup, dubLine) ||
				!checkMatrix1(lnks, j, nc, matrix, i + 1, i + 2, printError, errLine, errGroup, dubLine))
			{
				return false;
			}
		}
	}
	return true;
}
bool alldata::CheckMatrix(const char* matrix, int nl, int nc, bool printError, int* errLine, int* errGroup, int* dubLine)
{
	return _CheckMatrix(matrix, nl, nc, links(), printError, errLine, errGroup, dubLine);
}