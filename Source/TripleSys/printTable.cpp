#include <iostream>
#include "TripleSys.h"

void printTableColor(char const* name, const char *c, int nl, int nc, int ns, int np, bool makeString)
{
	printf("%s:\n", name);
	for (int j = 0; j < nl; j++)
	{
		if (makeString) printf("\"");
		for (int i = 0; i < nc; i++)
		{
			char v = c[j * nc + i];
			if (np > 0 && (i % np) == 0 && i > 0)
				printf(" ");
			if (v >= 0 && v < 67)
			{
				printf("\x1b[38;5;%dm%2d", 28 + v * 3, v);
				printf("\x1b[0m");
			}
			else
				printf("%2d", v);
		}
		if (j + 1 >= nl || ns <= 0 || ((j + 1) % ns) == 0)
			makeString ? printf(" \"\n") : printf("\n");
		else
			printf(" ");
	}
}
void printTable(char const* name, const char *c, int nl, int nc, int ns, int np, bool makeString)
{
	if (strlen(name) != 0)
		printf("%s:\n", name);
	for (int j = 0; j < nl; j++)
	{
		if (makeString) printf("\"");
		for (int i = 0; i < nc; i++)
		{
			if (np > 0 && (i % np) == 0 && i > 0)
				printf(" ");
			if (c[j * nc + i] == -1)
				printfGreen(" %3d", c[j * nc + i]);
			else
				printf(" %3d", c[j * nc + i]);
		}
		if (j + 1 >= nl || ns <= 0 || ((j + 1) % ns) == 0)
			makeString ? printf(" \"\n") : printf("\n");
		else
			printf(" ");
	}
}
void printTable(char const* name, const int *c, int nl, int nc, int ns, int np, bool makeString, double scale)
{
	printf("%s:\n", name);
	for (int j = 0; j < nl; j++)
	{
		if (makeString) printf("\"");
		for (int i = 0; i < nc; i++)
		{
			const auto v = c[j * nc + i];
			if (np > 0 && (i % np) == 0 && i > 0)
				printf(" ");
			scale > 0.0 ? printf(" %4.1f", v * scale) : printf(" %4d", v);
		}
		if (j + 1 >= nl || ns <= 0 || ((j + 1) % ns) == 0)
			makeString ? printf(" \"\n") : printf("\n");
		else
			printf(" ");
	}
}
void printTable(char const* name, const double *c, int nl, int nc, int ns, int np, bool makeString, double scale)
{
	printf("%s:\n", name);
	for (int j = 0; j < nl; j++)
	{
		if (makeString) printf("\"");
		for (int i = 0; i < nc; i++)
		{
			const auto v = c[j * nc + i] * scale;
			if (np > 0 && (i % np) == 0 && i > 0)
				printf(" ");
			if (v != 0)
				printfGreen(" %4.1f", v);
			else
			    printf(" %4.1f", v);
		}
		if (j + 1 >= nl || ns <= 0 || ((j + 1) % ns) == 0)
			makeString ? printf(" \"\n") : printf("\n");
		else
			printf(" ");
	}
}
void printThreadsStat(double *cntTotal, double *cnt, double dNumMatrices, int nrows, int mStep, int nThreads, clock_t iTime, bool bPrintSetup)
{
	double sum1 = 0.0, sum2 = 0;
	for (int i = 0; i < nThreads; i++)
	{
		int j = i * 2;
		double d = cnt[j];
		d = ((d < 0) ? -1.0 - d : d);
		sum1 += cntTotal[j] + d;
		sum2 += cntTotal[j + 1] + cnt[j + 1];
	}
	if (bPrintSetup)
		printf("\n\nSetup: %.0f %d-rows matrices, step=%d, %d threads\n", dNumMatrices, nrows, mStep, nThreads);
	else
		printf("\n");
	printf("%.0f non-isomorphic matrices (%d,%d,%d) were selected from %.0f generated, time=%d\nThread:Non-isomorphic/Generated Matrices:",
		sum1, nPlayers, (nPlayers - 1) / (GroupSize - 1), GroupSize, sum2, clock() - iTime);
	for (int i = 0; i < nThreads; i++)
	{
		if ((i % 10) == 0)
			printf("\n");
		int j = i * 2;
		double d = cnt[j];
		d = ((d < 0) ? -1.0 - d : d);
		printf(" %d:%.0f/%.0f", i + 1, cntTotal[j] + d, cntTotal[j + 1] + cnt[j + 1]);
	}
	printf("\n");
}
