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
			if (v == 1)
				printfRed(" %3d", v);
			else if (v == 2)
				printfGreen(" %3d", v);
			else if(v == 3)
				printfYellow(" %3d", v);
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
	printf("%s:\n", name);
	for (int j = 0; j < nl; j++)
	{
		if (makeString) printf(" \" ");
		for (int i = 0; i < nc; i++)
		{
			if (np > 0 && (i % np) == 0)
				printf(" ");
			if (c[j * nc + i] == -1)
				printfGreen(" %3d", c[j * nc + i]);
			else
				printf(" %3d", c[j * nc + i]);
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
		for (int i = 0; i < nc; i++)
		{
			const auto v = c[j * nc + i];
			if (np > 0 && (i % np) == 0)
				printf(" ");
			scale > 0.0 ? printf(" %4.1f", v * scale) : printf(" %3d", v);
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
		for (int i = 0; i < nc; i++)
		{
			const auto v = c[j * nc + i] * scale;
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