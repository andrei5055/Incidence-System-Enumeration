#pragma once
#include <vector>
#include "k-SysSupport.h"

#define USE_GPU	  0
#define nPlayers  21
#define GroupSize 3

// The maximum numbers of players and groups per day for which the program is compiled.
#define MAX_PLAYER_NUMBER			92
#define MAX_GROUP_NUMBER			(MAX_PLAYER_NUMBER / 2)

#define MAX_CYCLE_SETS		32
#define MAX_CYCLES_PER_SET	(MAX_PLAYER_NUMBER / 4)
#define USE_GROUP_4_2_ROWS   0  // The use of the Aut(M) of the 2-row matrix

#define unset ((tchar)(-1))

typedef enum {
	t_nonregular = 0,
	t_regular = 1 << 0,
	t_srg = t_regular + (1 << 1),
	t_4_vert = t_srg + (1 << 2),
	t_complete = -1
} t_graphType;

typedef enum {
	// indices of parameters with the integer type values
	t_numPlayers,
	t_numPlayersMax,
	t_groupSize,
	t_CBMP_Graph,			// Complete Balanced Multi-Partite Graph
	t_generateMatrixExample,
	t_u1f,
	t_use2RowsCanonization,
	t_useFastCanonizerForG2,
	t_submatrixGroupOrderMin,
	t_resultGroupOrderMin,
	t_useGPU,
	t_MultiThreading,
	t_numThreads,
	t_nRowsInStartMatrix,
	t_nRowsInResultMatrix,
	t_nMaxNumberOfStartMatrices,
	t_nFirstIndexOfStartMatrices,
	t_expectedResult,
	t_useCheckLinksV,
	t_useRowsPrecalculation,
	t_useSolutionCliquesAfterRow,
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
	t_out_CSV_file,
	t_outAutomorphismGroup,
	t_nestedGroups,
	t_gridSize,
	t_blockSize,
	t_orderMatrices,
	t_allowUndefinedCycles,
	t_any2RowsConvertToFirst2,
	t_exploreMatrices,
	t_semiSymmetricGraphs,
	t_rejectCycleLength,
	t_test,
	t_keepPrevResult,
	t_saveLatinSquareType,
	t_inputLatinSquareType,
	t_v4,
	t_v4Row,
	t_lastParam,
	// indices of parameters with the string type values
	t_UFname = 0,
	t_StartFolder,
	t_ResultFolder,
	t_ImprovedResultFolder,
	t_binaryCanonizer,
	t_testName,
	t_matrTest, // parameter, which will define the test to be launched for constructed matrix or matrices
	t_CSV_FileName,
	t_ResultsName,
	t_InputDataFileName,
	t_lastStrParam,
} paramID;

typedef enum {
	t_printNothing = 0,
	t_printExploringSRG = 1 << 6, // report the process of SRG's exploring 
} printFlags;

typedef struct {
	const char** paramNames;
	int numParams;
	std::vector<int>* m_pTmpParamStorage;
} paramDescr;

class kSysParam {
public:
	kSysParam()				{ Init(NULL); }
	kSysParam(const kSysParam* param) {
		*this = *param;
		Init(param->u1fCycles[0]);
	}
	~kSysParam()			{ 
		if (m_b_dataOwner) {
			delete[] u1fCycles[0];
			u1fCycles[0] = NULL;
		}
	}
	CC inline auto partiteNumb() const		{ return val[t_CBMP_Graph]; }
	CC inline auto completeGraph() const	{ return partiteNumb() == 1; }
	inline auto partitionSize() const		{ return val[t_numPlayers] / partiteNumb(); }
	CC inline auto numFactors() const		{ return m_numFactors; }
	void setup();

	int val[t_lastParam];
	int groupSizeFactorial;
	int m_numFactors;
	std::string* strVal[t_lastStrParam];
	tchar* u1fCycles[1];
	paramDescr* pParamDescr;
private:
	void Init(ctchar* pU1fCycles) {
		pParamDescr = NULL;
		memset(strVal, 0, sizeof(strVal));
		if (m_b_dataOwner = (pU1fCycles != NULL)) {
			const auto len = 1 + *pU1fCycles * MAX_CYCLES_PER_SET;
			u1fCycles[0] = new tchar[len];
			memcpy(u1fCycles[0], pU1fCycles, len);
		}
		else
			u1fCycles[0] = NULL;
	}
	bool m_b_dataOwner;
};

#define SECTION_PARAM_MAIN		"Main parameters:"
#define SECTION_PARAM_STRINGS	"Parameters of string type:"
#define SECTION_PARAM_U1F_CONF	"U1F configurations:"
