#pragma once
#include "k-SysSupport.h"

#define USE_GPU	  0
#define nPlayers  21
#define GroupSize 3

// The maximum numbers of players and groups per day for which the program is compiled.
#define MAX_PLAYER_NUMBER			33
#define MAX_GROUP_NUMBER			12
#define MAX_3PF_SETS			    218 // for 15 we need 13, for 21 - 40(54?), for 27 we need it to be 217 (probably)

#define MAX_CYCLE_SETS		30
#define MAX_CYCLES_PER_SET	 8
#define USE_GROUP_4_2_ROWS   0  // The use of the Aut(M) of the 2-row matrix

#define unset ((tchar)(-1))

typedef enum {
	// indices of parameters with the integer type values
	t_numPlayers,			 
	t_groupSize,
	t_u1f,
	t_use2RowsCanonization,
	t_submatrixGroupOrderMin,
	t_resultGroupOrderMin,
	t_useGPU,
	t_MultiThreading,
	t_numThreads,
	t_nRowsInStartMatrix,
	t_nMaxNumberOfStartMatrices,
	t_nFirstIndexOfStartMatrices,
	t_nRowsInResultMatrix,
	t_expectedResult,
	t_useCheckLinksV,
	t_useRowsPrecalculation,
	t_useSolutionCliquesAfterRow,
	t_useAutForPrecRows,
	t_lastRowSecondPlayer,
	t_printMatrices,
	t_savingMatricesToDisk,
	t_matrixCanonInterval,
	t_checkConstructedMatrices,
	t_useSS,
	t_p1f_counter,
	t_autLevelMinDef,
	t_autLevelMaxDef,
	t_autLevelMin,
	t_autLevelMax,
	t_autDirection,
	t_autGroupNumb,
	t_autSaveTestedTrs,
	t_useImproveMatrix,
	t_useCombinedSolutions,
	t_outAutomorphismGroup,
	t_nestedGroups,
	t_gridSize,
	t_blockSize,
	t_orderMatrices,
	t_allowMissingCycles,
	t_any2RowsConvertToFirst2,
	t_lastParam,
	// indices of parameters with the string type values
	t_UFname = 0,
	t_StartFolder,
	t_ResultFolder,
	t_ImprovedResultFolder,
	t_binaryCanonizer,
	t_testName,
	t_matrTest, // parameter, which will define the test to be launched for constructed matrix or matrices 
	t_lastStrParam
} paramID;

typedef struct {
	const char** paramNames;
	int numParams;
} paramDescr;

typedef struct {
	int val[t_lastParam];
	int groupSizeFactorial;
	std::string* strVal[t_lastStrParam];
	tchar* u1fCycles[1];
	paramDescr* pParamDescr;
} kSysParam;


