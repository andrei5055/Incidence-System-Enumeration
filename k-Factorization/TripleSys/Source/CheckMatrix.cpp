#include "TripleSys.h"

char lastError[256];

bool checkMatrix1(tchar* lnks, int ln, int nc, int gs, const tchar* matrix, int i1, int i2, bool printError, int* errLine, int* errGroup, int* dubLine)
{
	const auto* c = matrix + ln * nc;
	auto a = c[i1], b = c[i2];
	if (a == unset || a >= nc ||
		b == unset || b >= nc ||
		a == b)
	{
		printf("CheckMatrix: matrix values (%d %d) in group %d on line %d are incorrect, job aborted\n", a, b, i1 / gs, ln);
		abort();
	}
	tchar dLine = lnks[a * nc + b];
	if (dLine != unset)
	{
		sprintf_s(lastError, "CheckMatrix: matrix pair (%d %d) in group %d on line %d already defined in line %d\n", a, b, i1 / gs, ln, dLine);
		if (printError)
			printf(lastError);

		*errLine = ln;
		*errGroup = i1 / gs;
		*dubLine = dLine;
		return false;
	}
	lnks[a * nc + b] = lnks[b * nc + a] = ln;
	return true;
}
bool _CheckMatrix(const tchar* matrix, int nrows, int nc, int gs, tchar *lnks, bool printError, int* errLine, int* errGroup, int* dubLine, int nr)
{
	if (gs <= 1 || (nr < nrows))
	{
		printf("CheckMatrix: incorrect parameters (nrows=%d, ncolumns=%d, group size=%d), job aborted\n", nrows, nc, gs);
		abort();
	}

	memset(lnks, unset, nc * nc);

	for (int j = 0; j < nrows; j++)
	{
		for (int i = 0; i < nc; i = i + gs)
		{
			for (int k = 0; k < gs - 1; k++)
			{
				for (int m = k + 1; m < gs; m++)
				{
					if (!checkMatrix1(lnks, j, nc, gs, matrix, i + k, i + m, printError, errLine, errGroup, dubLine))
					{
						return false;
					}
				}
			}
		}
	}
	return true;
}
bool alldata::CheckMatrix(const tchar* matrix, int nl, int nc, int gs, bool printError, int* errLine, int* errGroup, int* dubLine)
{
	return _CheckMatrix(matrix, nl, nc, gs, links(), printError, errLine, errGroup, dubLine, sysParam()->numFactors());
}