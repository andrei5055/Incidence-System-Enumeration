#pragma once
#include <thread>
#include <vector>
#include "TripleSys.h"
#define LIBRARY_EXPORTS 1
#ifdef LIBRARY_EXPORTS
#    define LIBRARY_API __declspec(dllexport)
#else
#    define LIBRARY_API __declspec(dllimport)
#endif
class TopGun : private SizeParam {
public:
	LIBRARY_API TopGun(int numPlayers, int groupSize);

	LIBRARY_API ~TopGun();
	bool LIBRARY_API Run();
private:
	int getStartMatrices();
	void myTemporaryCheck();
	void startThread(int iTask);
	void threadStopped(int iTask);
	void waitAllThreadFinished();
	sLongLong m_cnt[NThreads * 2];
	sLongLong m_cntTotal[NThreads * 2];
	sLongLong dNumMatrices[2];
	char startMatricesFileName[256];
	char startMatricesFullFileName[256];
	bool threadActive[NThreads];
	int nMatricesMax;
	int nMatrices;
	int nRowsStart = NRowsInStartMatrix;
	int nRowsOut = NRowsInResultMatrix;
	int mStartMatrixSize = nPlayers * nRowsStart;
	int mLinksSize = nPlayers * nPlayers;
	clock_t cTime, rTime, mTime, iTime;
	char* startMatrix;
	char* mStartLinks;
	char* mstart;
	int numThreads;
	int iMatrix;
	int iPrintCount;
	int iTaskSeq;
	std::vector<std::thread> threads;//(NThreads);
};
void RunThread(int threadNumber, int iMode,
	char* mstart0, char* mstart, int nRowsStart, int nRowsOut, sLongLong* pcnt, bool bPrint);