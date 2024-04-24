#include "TopGun.h"
#include "data.h"

TopGun::TopGun(int numPlayers, int groupSize) :
	SizeParam((numPlayers - 1) / (groupSize - 1), numPlayers, groupSize) {

	dNumMatrices[0] = nMatricesMax = MaxNumberOfStartMatrices;
	startMatrix = (char*)malloc(nMatricesMax * mStartMatrixSize);
	mStartLinks = new char[nMatricesMax * mLinksSize];
	if (startMatrix == NULL || mStartLinks == NULL)
	{
		printf("Not enough memory for initial %d-rows %d matrices. Exit\n", nRowsStart, nMatricesMax);
		exit(1);
	}
	memset(cntTotal, 0, sizeof(cntTotal));
	memset(cnt, 0, sizeof(cnt));
	memset(threadActive, false, sizeof(threadActive));
	iMatrix = 0;
	iPrintCount = 0;
	iTaskSeq = 0;
}

TopGun::~TopGun() {
	free(startMatrix);
	delete[] mStartLinks;
}

bool TopGun::Run()
{
	iTime = clock();

	if (nRowsOut == 0)
		nRowsOut = m_numDays;

	if (!UseMultiThreading)
	{
		if (nRowsOut < 2 || nRowsOut > m_numDays)
		{
			printf("Number of rows(%d) for result matrices must be in range 2:%d\n",
				nRowsOut, m_numDays);
			exit(1);
		}
		alldata sys(m_numPlayers);
		sys.initStartValues(ivc);// can be used to start from previous result
		sys.Run(0, eCalcResult, (char*)0, (char*)0, nRowsStart, nRowsOut, NULL, true);
	}
	else
	{
		if (nRowsStart < 2 || nRowsStart > m_numDays)
		{
			printf("Number of rows(%d) in 'Start' matrices (input for threads) must be in range 2:%d\n", 
				nRowsStart, m_numDays - 1);
			exit(1);
		}
		if (getStartMatrices() < 1)
		{
			printf("Cant get/calculate 'Start' matrices for threads. Exit\n");
			exit(1);
		}

		char* tsm = (char*)realloc(startMatrix, nMatrices * mStartMatrixSize);
		if (tsm != NULL)
			startMatrix = tsm;

		myTemporaryCheck();
		delete[] mStartLinks;
		mStartLinks = NULL;


		//exit(0);

		numThreads = (nMatrices > NThreads) ? NThreads : nMatrices;
		threads.resize(numThreads);

		int nThreadsRunning = 1;
		int nMatricesProc = 0;
		mstart = startMatrix;
		mTime = clock() - iTime;

		printfYellow("\nMultithread Matrices Calculation started (time=%dsec)\n", mTime / 1000);

		cTime = clock();
		while (nThreadsRunning > 0)
		{
			int iTask = 0;
			nThreadsRunning = 0;
			for (auto& t : threads)
			{
				if (!threadActive[iTask])
				{
					if (iMatrix < nMatrices)
					{
						//printTable("Input matrix", mstart, nRowsStart, m_numPlayers, 0, m_groupSize, true);
						startThread(iTask);
						nThreadsRunning++;
					}
				}
				else
				{
					nThreadsRunning++;
					if (cnt[iTask * 2] >= 0)
					{
						// thread finished
						nMatricesProc++;
						t.join();
						threadStopped(iTask);
					}
				}
				iTask++;
			}
			if (nThreadsRunning == 0)
				break;

			std::this_thread::sleep_for(std::chrono::milliseconds(5));

			if (clock() - cTime > 20000)
			{
				printThreadsStat(cntTotal, cnt, nMatrices, nMatricesProc, nRowsStart, nRowsOut, numThreads, iTime, ((iPrintCount++) % 10) == 0);
				cTime = clock();
			}
		}
		waitAllThreadFinished();
		rTime = clock() - iTime;
		printThreadsStat(cntTotal, cnt, nMatrices, iMatrix, nRowsStart, nRowsOut, numThreads, iTime, true);
		printf("Total time=%dms (include prep time=%dms)\n", rTime, mTime);
	}
	std::cout << "\7" << std::endl; // play sound
	return 0;
}
