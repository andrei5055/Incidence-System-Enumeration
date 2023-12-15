#include <iostream>
#include "TripleSys.h"

void printTableColor(char const* name, const char *c, int nl, int nc, int ns, int np, bool makeString)
{
	printf("%s:\n", name);
	for (int j = 0; j < nl; j++)
	{
		if (makeString) printf(" \" ");
		for (int i = 0; i < nc; i++)
		{
			char v = c[j * nc + i];
			if (np > 0 && (i % np) == 0)
				printf(" ");
			if (v >= 0 && v < 67)
			{
				printf("\x1b[38;5;%dm %3d", 28 + v * 3, v);
				printf("\x1b[0m");
			}
			else
				printf(" %3d", v);
		}
		if (j + 1 >= nl || ns <= 0 || ((j + 1) % ns) == 0)
			makeString ? printf(" \\n\"\n") : printf("\n");
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
		if (makeString) printf(" \" ");
		for (int i = 0; i < nc; i++, c++)
		{
			if (np > 0 && (i % np) == 0)
				printf(" ");
			if (*c == -1)
				printfGreen(" %3d", *c);
			else
				printf(" %3d", *c);
		}
		if (j + 1 >= nl || ns <= 0 || ((j + 1) % ns) == 0)
			makeString ? printf(" \\n\"\n") : printf("\n");
		else
			printf(" ");
	}
}

void printTable(char const* name, const int *c, int nl, int nc, int ns, int np, bool makeString, double scale)
{
	printf("%s:\n", name);
	for (int j = 0; j < nl; j++)
	{
		if (makeString) printf(" \" ");
		for (int i = 0; i < nc; i++, c++)
		{
			if (np > 0 && (i % np) == 0)
				printf(" ");
			scale > 0.0 ? printf(" %4.1f", *c * scale) : printf(" %3d", *c);
		}
		if (j + 1 >= nl || ns <= 0 || ((j + 1) % ns) == 0)
			makeString ? printf(" \\n\"\n") : printf("\n");
		else
			printf(" ");
	}
}
void printTable(char const* name, const double *c, int nl, int nc, int ns, int np, bool makeString, double scale)
{
	printf("%s:\n", name);
	for (int j = 0; j < nl; j++)
	{
		if (makeString) printf(" \" ");
		for (int i = 0; i < nc; i++, c++)
		{
			const auto v = *c;
			if (np > 0 && (i % np) == 0)
				printf(" ");
			if (v != 0)
				printfGreen(" %4.1f", v);
			else
			    printf(" %4.1f", v);
		}
		if (j + 1 >= nl || ns <= 0 || ((j + 1) % ns) == 0)
			makeString ? printf(" \\n\"\n") : printf("\n");
		else
			printf(" ");
	}
}
