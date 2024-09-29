#pragma once
#include "k-SysSupport.h"

#define USE_GPU	  0
#define nPlayers  16
#define GroupSize 3

// The maximum numbers of players and groups per day for which the program is compiled.
#define MAX_PLAYER_NUMBER			27
#define MAX_GROUP_NUMBER			9
#define MAX_UNIFOM_CONF_LENGTH		4
#define MAX_3PF_SETS			    218 // for 15 we need 13, for 21 - 40(54?), for 27 we need it to be 217 (probably)
#define MAX_3PF_SECOND_ROWS			19 // need 9-18+1? for 21 (not p1f first two rows)

#define MAX_CYCLE_SETS 30
#define MAX_CYCLES_PER_SET 6

#define unset ((tchar)(-1))

typedef enum {
	// indices of parameters with the integer type values
	t_numPlayers,			 
	t_groupSize,
	t_u1f,
	t_p1f,
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
	t_useCheckLinksH,
	t_useCheckLinksV,
	t_printMatrices,
	t_savingMatricesToDisk,
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
	t_outAutomorphismGroup,
	t_nestedGroups,
	t_binaryCanonizer,
	t_gridSize,
	t_blockSize,
	t_lastParam,
	// indices of parameters with the string type values
	t_UFname = 0,
	t_StartFolder,
	t_ResultFolder,
	t_ImprovedResultFolder,
	t_testName,
	t_matrTest, // parameter, which will define the test to be launched for constructed matrix or matrices 
	t_lastStrParam
} paramID;

typedef struct {
	int val[t_lastParam];
	int groupSizeFactorial;
	std::string* strVal[t_lastStrParam];
	tchar* u1f[1];
} kSysParam;


