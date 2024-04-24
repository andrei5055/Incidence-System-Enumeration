#pragma once
#include "TripleSys.h"
#include <thread>
#include <vector>
class TopGun : private SizeParam {
public:
	TopGun(int numPlayers, int groupSize);

	~TopGun();
	bool Run();
private:
	int getStartMatrices();
	void myTemporaryCheck();
	void startThread(int iTask);
	void threadStopped(int iTask);
	void waitAllThreadFinished();
	sLongLong cnt[NThreads * 2];
	sLongLong cntTotal[NThreads * 2];
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