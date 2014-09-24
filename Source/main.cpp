// BIBD.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "C_tDesignEnumerator.h"
#include "EnumInfo.h"


typedef struct {
#if T_DESIGN_ENUM
    int t;
#endif
	int v;
	int k;
	int lambda;
} tDesign_param;


int find_T_designParam(int v, int k, int lambda)
{
	static int cntr; cntr++; int lam = lambda;
	int prevLambda, lambdas[10];
	int i = 1;
	do {
		lambdas[++i - 2] = prevLambda = lambda;
		lambda = prevLambda * (k - i) / (v - i);
	} while (lambda && lambda * (v - i) / (k - i) == prevLambda);

    if (i > 2)
        printf("{%d, %3d, %2d, %2d},\n", i, v, k, lambdas[i - 2]);
    
	return i;
}

int _tmain(int argc, _TCHAR* argv[])
{
	_CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_FILE);
	_CrtSetReportFile(_CRT_WARN, _CRTDBG_FILE_STDOUT);
	_CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_FILE);
	_CrtSetReportFile(_CRT_ERROR, _CRTDBG_FILE_STDOUT);
	_CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_FILE);
	_CrtSetReportFile(_CRT_ASSERT, _CRTDBG_FILE_STDOUT);

    tDesign_param bibd_param[] = {
#if T_DESIGN_ENUM
		{3, 8, 4, 3},
		{0, 0, 0, 0},
		{4, 11, 5, 3},
		{0, 0, 0, 0},
        #include "t_design_param.txt"
#else
		{ 6, 3, 2 },
		{ 0, 0, 0 },
        #include "BIBD_param.txt"
            { 0, 0, 0 },
		{9, 4, 6},
    { 0, 0, 0 },
//		{ 19, 9, 4 },
		{ 25, 5, 1 },
		{ 25, 9, 3 },
		{ 29, 7, 3 },
//        {16, 6, 4},
//        {15, 3, 2},

        
  //      {13, 3, 2},
  //      {13, 4, 1},
  //      {13, 4, 2},
  //      {13, 6, 5},
//        {8, 4, 12},
//        { 6, 3, 2 },
//		{9, 3, 2},
//		{9, 4, 3},
		{ 0, 0, 0 },
//#include "TestSet.txt"
//		{ 18, 6, 5 },
//		{ 7, 3, 9 },
//		{8, 4, 12},
//		{ 9, 3, 2 },  //            36 (13)  [4]
//		{ 6, 3, 2 },  //
//		{ 9, 3, 3 },  //+       22,521 (332) [2], [3], [4], [5]
//		{ 9, 3, 4 },
//#include "V_9.txt"
//#include "V_10.txt"
#include "BIBD_param.txt"
#endif
	};
    
	const char *pSummaryFile = "EnumSummary.txt";
	outString("         BIBDs:                     Canonical:      NRB #:      Constructed:    Run Time (sec):\n", pSummaryFile, "w");

	char buff[256], buffer[256];
	int i;
	for (i = 0; i < countof(bibd_param); i++) {
        tDesign_param *pParam = bibd_param + i;
		if (pParam->v <= 0)
			break;


#if T_DESIGN_ENUM
//        find_T_designParam(pParam->v, pParam->k, pParam->lambda);
//        continue;        
		C_tDesign bibd(pParam->t, pParam->v, pParam->k, pParam->lambda);
        #define DesignEnumerator    C_tDesignEnumerator
#else
		C_BIBD bibd(pParam->v, pParam->k, pParam->lambda);
        #define DesignEnumerator    CBIBD_Enumerator
#endif
		DesignEnumerator BIBD_Enum(&bibd, false, NO_REPLICATED_BLOCKS);
		BIBD_Enum.makeJobTitle(buff, countof(buff));
		std::cout << buff;
		CInsSysEnumInfo enumInfo(buff);
		try {
			BIBD_Enum.Enumerate(PRINT_TO_FILE, &enumInfo);
			enumInfo.reportResult(buffer, countof(buffer));
			outString(buffer, pSummaryFile);
			std::cout << "\xd" << buffer;
		}
		catch (...)  {
			BIBD_Enum.closeFile();	
			break;
		}
	}

	_CrtDumpMemoryLeaks();
	return i;
}

