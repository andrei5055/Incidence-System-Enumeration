#pragma once
#include <thread>
#include "TripleSys.h"

class RowDB : public CStorageSet<tchar>, public kSysParam {
public:
	RowDB(const kSysParam& param) : CStorageSet<tchar>(10, param.val[t_numPlayers]), kSysParam(param) {}
	bool isValid(const kSysParam* param) {
		ctchar* p0 = param->u1fCycles[0];
		ctchar* p1 = u1fCycles[0];
		if ((p0 && !p1) || (!p0 && p1))
			return false;

		if (p0 && p1 && (*p0 != *p1 || memcmp(p0+1, p1+1, (*p0) * MAX_CYCLES_PER_SET))) {
			return false;
		}

		const paramID IDs[] = { t_numPlayers, t_groupSize, t_use2RowsCanonization, t_u1f, t_CBMP_Graph,
			t_allowUndefinedCycles, t_any2RowsConvertToFirst2, t_binaryCanonizer };
		for (int i = 0; i < countof(IDs); i++) {
			const auto id = IDs[i];
			if (param->val[id] != val[id])
				return false;
		}
		return true;
	}
};

class TopGunBase : public SizeParam, public MatrixDB {
public:
	TopGunBase(const kSysParam& param);
	virtual ~TopGunBase()					{ 
		delete[] m_pInputMatrices;
		delete[] cnt();
		delete m_pMatrixInfo;
		delete[] m_pMatrixPerm;
	}
	int virtual Run() = 0;
	void K_SYS_LIBRARY_API outputIntegratedResults(const paramDescr *pParSet = NULL, int numParamSet = 0, const char* pResults = "_Results.txt");
	inline auto numPlayers() const			{ return m_numPlayers; }
	inline auto groupSize() const			{ return m_groupSize; }
	inline auto groupSizeFactorial() const	{ return m_groupSizeFactorial; }
	inline auto numMatrices2Process() const	{ return m_nMatrices; }
	inline auto nRowsStart() const			{ return param(t_nRowsInStartMatrix); }
	inline auto nRowsOut() const			{ return m_nRowsOut; }
	inline sLongLong* cnt() const			{ return m_cnt; }
	inline const auto* inputMatrices() const { return m_pInputMatrices; }
	inline auto* paramPtr() const			{ return &m_param; }
protected:
	inline int param(paramID id) const		{ return m_param.val[id]; }
	bool readStartMatrices(const char* pMatrixType = "Start");
	void InitCnt(size_t nThrds) { m_cnt = new sLongLong[2 * nThrds]; memset(m_cnt, 0, 2 * nThrds * sizeof(long long)); }
	inline auto nMatricesMax() const		{ return param(t_nMaxNumberOfStartMatrices); }
	inline auto nMatricesReserved() const	{ return nMatricesMax() > 100000 ? 100000 : nMatricesMax(); } // Memory will be reserved for the specified number of matrices before the reading process begins.
	inline auto inputMatrixSize() const		{ return m_nInputMatrixSize; }
	int orderMatrices(int orderMatrixMode);

	int m_nRowsOut;
	int m_nInputMatrixSize;
	tchar* m_pInputMatrices = NULL;
	uint* m_pMatrixPerm = NULL;
	CMatrixInfo* m_pMatrixInfo = NULL;	// Information about loaded matrices: |Aut(M)|, cycle's, group's
	int m_nMatrices;
	sLongLong *m_cnt = NULL;
	const kSysParam m_param;
	std::string m_reportInfo;
private:
	inline void setInputMatrixSize(int size){ m_nInputMatrixSize = size; }
	int loadMatrices(const char* pMatrixType);
	int readStartData(const std::string& fn, int nTotal, tchar** ppSm, int nm, int &reservedElem, CMatrixInfo *pMatrixInfo = NULL) const {
		return readTable(fn, nRowsStart(), numPlayers(), nm, nTotal, ppSm, reservedElem, nMatricesMax(), pMatrixInfo);
	}
};

class TopGun : public TopGunBase {
public:
	K_SYS_LIBRARY_API TopGun(const kSysParam& param);

	K_SYS_LIBRARY_API ~TopGun();
	int K_SYS_LIBRARY_API Run();
	inline static K_SYS_LIBRARY_API auto secondRowDB()		{ return m_pSecondRowsDB; }
private:
	sLongLong printThreadsStat(int nMatrices, int nProcessed, const clock_t& iTime, bool bPrintSetup);
	void myTemporaryCheck();
	void startThread(int iTask, int iTaskId, eThreadStartMode iMode = eCalcResult, CRowStorage* pRowStorage = NULL);
	void threadStopped(int iTask);
	void waitAllThreadFinished();
	sLongLong *m_cntTotal = NULL;
	sLongLong dNumMatrices[2];
	bool *threadActive = NULL;
	int mLinksSize;
	clock_t cTime = 0, rTime = 0, mTime = 0, iTime = 0;
	ctchar* mstart = NULL, *mfirst = NULL;
	int numThreads;
	int m_iMatrix;
	int m_iPrintCount;
	std::vector<std::thread> threads;//(NThreads);
	static RowDB* m_pSecondRowsDB;
};

class TopGunGPU : public TopGunBase {
public:
	K_SYS_LIBRARY_API TopGunGPU(const kSysParam& param) : TopGunBase(param)		{}
	K_SYS_LIBRARY_API ~TopGunGPU();
	K_SYS_LIBRARY_API int Run();
	inline int gridSize() const		{ return param(t_gridSize); }
	inline int blockSize() const	{ return param(t_blockSize); }  // 32 x 28 also works
};
