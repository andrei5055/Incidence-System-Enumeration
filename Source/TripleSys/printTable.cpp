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
void printThreadsStat(sLongLong* cntTotal, sLongLong* cnt, int nMatrices, int nProccesed, 
	int nrowsStart, int nrowsOut, int nThreads, clock_t iTime, bool bPrintSetup)
{
	sLongLong sum1 = 0, sum2 = 0;
	for (int i = 0; i < nThreads; i++)
	{
		int j = i * 2;
		sLongLong d = cnt[j];
		d = ((d < 0) ? -1 - d : d);
		sum1 += cntTotal[j] + d;
		sum2 += cntTotal[j + 1] + cnt[j + 1];
	}
	printf("T=%dsec: %d (from %d) matrices (%dx%d) processed by %d threads. %zd matrices (%dx%d) generated\n", 
		(clock() - iTime) / 1000, nProccesed, nMatrices, nPlayers, nrowsStart, nThreads, sum1, nPlayers, nrowsOut);

	if (bPrintSetup)
	{
		printf("Thread:Matrices generated");
		for (int i = 0; i < nThreads; i++)
		{
			if ((i % 12) == 0)
				printf("\n");
			int j = i * 2;
			sLongLong d = cnt[j];
			d = ((d < 0) ? -1 - d : d);
			printf(" %d:%zd", i + 1, cntTotal[j] + d);
		}
		printf("\n");
	}
}
void printTransformed(int nrows, int ncols, const char* tr, const char* ttr, const char* pImatr, const char* pTmatr, int numRow, sLongLong nLoops, int finalKMindex)
{
	if (nLoops)
		printf("Calculated Matrix %zd can be improved. See below (nKm=%d row=%d)\n", nLoops, finalKMindex, numRow);
	
	printTable("Tr source", tr, 1, ncols);
	printTable("Tr actual", ttr, 1, ncols);
	printTable("Original", pImatr, nrows, ncols);
	printTable("Translated", pTmatr, nrows, ncols);
}
