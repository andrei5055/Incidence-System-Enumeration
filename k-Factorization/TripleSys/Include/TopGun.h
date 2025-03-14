#pragma once
#include <thread>
#include "TripleSys.h"

class TopGunBase : public SizeParam, public MatrixDB {
public:
	TopGunBase(const kSysParam& param);
	virtual ~TopGunBase()					{ 
		delete[] startMatrix;
		delete[] cnt();
		delete[] m_pMatrixAutOrder;
		delete[] m_pMatrixPerm;
	}
	int virtual Run() = 0;
	void K_SYS_LIBRARY_API outputIntegratedResults(const paramDescr *pParSet = NULL, int numParamSet = 0, const char* pResults = "_Results.txt") const;
	inline auto numPlayers() const			{ return m_numPlayers; }
	inline auto groupSize() const			{ return m_groupSize; }
	inline auto groupSizeFactorial() const	{ return m_groupSizeFactorial; }
	inline auto numMatrices2Process() const	{ return nMatrices; }
	inline auto nRowsStart() const			{ return param(t_nRowsInStartMatrix); }
	inline auto nRowsOut() const			{ return m_nRowsOut; }
	inline sLongLong* cnt() const			{ return m_cnt; }
	inline auto* pntrStartMatrix() const	{ return startMatrix; }
	inline auto* paramPtr() const			{ return &m_param; }
protected:
	inline int param(paramID id) const		{ return m_param.val[id]; }
	bool readStartMatrices();
	void InitCnt(size_t nThrds) { m_cnt = new sLongLong[2 * nThrds]; memset(m_cnt, 0, 2 * nThrds * sizeof(long long)); }
	inline auto nMatricesMax() const		{ return param(t_nMaxNumberOfStartMatrices); }
	inline auto nMatricesReserved() const	{ return nMatricesMax() > 100000 ? 100000 : nMatricesMax(); } // Memory will be reserved for the specified number of matrices before the reading process begins.
	void orderMatrices(int orderMatrixMode);

	int m_nRowsOut;
	int mStartMatrixSize;
	tchar* startMatrix;
	uint* m_pMatrixAutOrder = NULL;
	uint* m_pMatrixPerm = NULL;
	int nMatrices;
	sLongLong *m_cnt = NULL;
	const kSysParam m_param;
	std::string m_reportInfo;
private:
	int getStartMatrices();
	int readStartData(const std::string& fn, int nTotal, tchar** ppSm, int nm, int &reservedElem, uint ** ppGroupOrders = NULL) const {
		return readTable(fn, nRowsStart(), numPlayers(), nm, nTotal, ppSm, reservedElem, nMatricesMax(), ppGroupOrders);
	}
};

class TopGun : public TopGunBase {
public:
	K_SYS_LIBRARY_API TopGun(const kSysParam& param);

	K_SYS_LIBRARY_API ~TopGun();
	int K_SYS_LIBRARY_API Run();
private:
	sLongLong printThreadsStat(int nMatrices, int nProcessed, const clock_t& iTime, bool bPrintSetup);
	void myTemporaryCheck();
	void startThread(int iTask, eThreadStartMode iMode = eCalcResult, bool bOnlyStart = false, CRowStorage* pRowStorage = NULL);
	void threadStopped(int iTask);
	void waitAllThreadFinished();
	sLongLong *m_cntTotal = NULL;
	sLongLong dNumMatrices[2];
	bool *threadActive = NULL;
	int mLinksSize;
	clock_t cTime = 0, rTime = 0, mTime = 0, iTime = 0;
	tchar* mstart = NULL, *mfirst = NULL;
	int numThreads;
	int m_iMatrix;
	int m_iPrintCount;
	int m_iTaskSeq;
	std::vector<std::thread> threads;//(NThreads);
	CStorageSet<tchar>* m_pSecondRowsDB = NULL;
};

class TopGunGPU : public TopGunBase {
public:
	K_SYS_LIBRARY_API TopGunGPU(const kSysParam& param) : TopGunBase(param)		{}
	K_SYS_LIBRARY_API ~TopGunGPU();
	K_SYS_LIBRARY_API int Run();
	inline int gridSize() const		{ return param(t_gridSize); }
	inline int blockSize() const	{ return param(t_blockSize); }  // 32 x 28 also works
};
