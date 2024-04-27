#include <iostream>
#include <fstream>
#include <filesystem>
#include "TripleSys.h"

void createStartFolderAndFileName(char* fn, size_t fns, const char* folder, const char* fileNameFmt, int np, int nr, int gs)
{
	namespace fs = std::filesystem;
	if (folder != NULL && strlen(folder) > 0)
	{
		if (fns > 0 && fn != NULL)
		{
			int cnt = sprintf_s(fn, fns, "%s", folder);
			if (fs::create_directories(fn))
			{
				printf("Error: Cant create folder %s\n", fn);
				exit(1);
			}
			sprintf_s(fn + cnt, fns - cnt, "SM_%dx%dx%d.txt", np, nr, gs);
		}
	}
}
void createFolderAndFileName(char* fn, size_t fns, const char* folder, const char* fileNameFmt, int np, int nr, int gs, int iset)
{
	namespace fs = std::filesystem;
	if (folder != NULL && strlen(folder) > 0)
	{
		if (fns > 0 && fn != NULL)
		{
			int cnt = sprintf_s(fn, fns, "%s%dx%dx%d/", folder, np, nr, gs);
			if (fs::create_directories(fn))
			{
				printf("Error: Cant create folder %s\n", fn);
				exit(1);
			}
			sprintf_s(fn + cnt, fns - cnt, fileNameFmt, iset);
		}
	}
}
int readStartData(char* fn, char* sm, int nmax, int np, int nr, int gs)
{
	int nmm = 0, npm = 0, nrm = 0, gsm = 0;
	if (strlen(fn) == 0)
		return 0;
	std::ifstream mf;
	mf.open(fn, std::ios::in);
	if (!mf)
		return 0;
#if 0
	mf >> nmm >> npm >> nrm >> gsm;
	if (nmm < 1 || nmm > nmax)
	{
		printf("*** Value of number of files (%d) in the first line of %s\nis out of range (1:%d)\n", nmm, fn, nmax);
		exit(1);
	}
	if (npm != np || nrm != nr || gsm != gs)
	{
		printf("*** First line of file %s must be:\n%d %d %d %d\n", fn, npm, np, nr, gs);
		exit(1);
	}
	mf.ignore(256, '\n');
	int i = 0;
	int iCheck = 0;
	int nl = 1;
	while (i < nmm)
	{
		for (int j = 0; j < nr; j++) {
			for (int k = 0; k < np; k++) {
				int iv = 0;
				mf >> iv;
				if (mf.eof() || (k == 0 && iv != 0) || (gs == 2 && k == 1 && iv != j + 1))
				{
					printf("Error in data line %d\n", nl);
					exit(0);
				}
				*(sm++) = iv;
			}
			nl++;
		}
		i++;
		mf >> iCheck;
		nl++;
		if (mf.eof() || iCheck != i && iCheck != -i)
		{
			printf("*** Error in file %s line=%d:\nMatrix number %d, expected %d. Exit\n", fn, nl, iCheck, i);
			exit(1);
		}
	}
	if (iCheck != -nmm)
	{
		printf("*** Error in file %s :\nLast Matrix number (at the end of the file)=%d, expected %d. Exit\n", 
			fn, iCheck, -nmm);
		exit(1);
	}
#else
	int i = 0;
	int nl = 0;
	while (i < nmax)
	{
		char fb;
		mf >> fb;
		if (i == 1097)
			i = i;
		if (mf.eof())
			break;
retry:	if (fb != '\"')
		{
			mf.ignore(256, '\n');
			nl++;
		}
		else
		{
	        char* sms = sm;
			for (int j = 0; j < nr; j++) {
				if (j > 0)
					mf >> fb;
				if (mf.eof() || fb != '\"' || fb == '\0' || fb == '\n')
				{
					printf("Error in data\n");
					sm = sms;
					goto retry;
					exit(0);
				}
				for (int k = 0; k < np; k++) {
					int iv = 0;
					mf >> iv;
					if (mf.eof() || (k == 0 && iv != 0))
					{
						printf("Error in data line %d, file: %s\n", nl, fn);
						exit(0);
					}
					*(sm++) = iv;
				}
				nl++;
				mf.ignore(256, '\n');
			}
			i++;
		}
	}
#endif
	mf.close();
	return i;
}

void saveStartData(char* fn, char* sm, int nm, int np, int nr, int gs)
{
	std::ofstream mf;
	mf.open(fn, std::ios::out);
	if (!mf)
	{
		printf("*** Cant open file %s to save start %d-rows matrices\n", fn, nr);
		exit(1);
	}
	mf << nm << ' ' << np << ' ' << nr << ' ' << gs;
	mf << " : Number of matrices, Number of players, Number of rows, Group size" << std::endl;
	for (int i = 0; i < nm; i++) {
		for (int j = 0; j < nr; j++) {
			for (int k = 0; k < np; k++) {
				mf << (int)(*(sm++)) << ' ';
			}
			mf << std::endl;
		} 
		mf << ((i + 1 == nm) ? -nm : i + 1) << std::endl;
	}
	mf.close();
}
