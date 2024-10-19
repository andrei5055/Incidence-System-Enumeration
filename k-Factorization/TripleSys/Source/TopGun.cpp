#include "TopGun.h"
#include "data.h"

TopGun::TopGun(const kSysParam& param) : TopGunBase(param) {
	numThreads = param.val[t_numThreads];

	dNumMatrices[0] = nMatricesMax();
	mLinksSize = numPlayers() * numPlayers();
	if (!startMatrix)
	{
		printfRed("*** Not enough memory for initial %d-rows %d matrices. Exit\n", nRowsStart(), nMatricesMax());
		myExit(1);
	}
	if ((m_p1f || param.val[t_u1f]) && m_groupSize > 3)
	{
		printfRed("*** u1f and/or p1f cannot be used with 'GroupSize'=%d. Exit\n", m_groupSize);
		myExit(1);
	}

	if (this->param(t_MultiThreading)) {
		m_cntTotal = new sLongLong[2 * numThreads];
		memset(m_cntTotal, 0, 2 * numThreads * sizeof(m_cntTotal[0]));
		threadActive = new bool[numThreads];
		memset(threadActive, false, numThreads * sizeof(threadActive[0]));
	}

	m_iTaskSeq = m_iPrintCount = m_iMatrix = 0;
}

TopGun::~TopGun() {
	delete[] m_cntTotal;
	delete[] threadActive;
}

int TopGun::Run()
{
	iTime = clock();
	sLongLong resultMatr = 0;
	if (!param(t_MultiThreading))
	{
		alldata sys(*this, paramPtr());
		sys.initStartValues(ivc);// can be used to start from previous result
		resultMatr = sys.Run(1, eCalcResult, NULL, NULL, nRowsStart(), nRowsOut(), NULL, &m_reportInfo);
		transferMatrixDB(sys.matrixDB());
	}
	else
	{
		if (!readStartMatrices())
			myExit(1);

		//myTemporaryCheck();

		const auto firstIndexOfStartMatrices = param(t_nFirstIndexOfStartMatrices);
		if (nMatrices <= firstIndexOfStartMatrices)
		{
			printfRed("*** Value of FirstIndexOfStartMatrices(%d) must be from 0 to number of 'Start matrices'(%d). Exit\n", 
				firstIndexOfStartMatrices, nMatrices);
			myExit(1);
		}

		if (numThreads > nMatrices - firstIndexOfStartMatrices)
			numThreads = nMatrices - firstIndexOfStartMatrices;

		threads.resize(numThreads);

		InitCnt(numThreads);
		int nThreadsRunning = 1;
		int nMatricesProc = m_iMatrix = firstIndexOfStartMatrices;
		mstart = startMatrix + m_iMatrix * mStartMatrixSize;
		mTime = clock() - iTime;

		printfYellow("\nMultithread Matrices Calculation started (time=%dsec)\n", mTime / 1000);

		cTime = clock();
		while (nThreadsRunning > 0 || m_iMatrix < nMatrices)
		{
			int iTask = 0;
			nThreadsRunning = 0;
			for (auto& t : threads)
			{
				if (!threadActive[iTask])
				{
					if (m_iMatrix < nMatrices)
					{
						//printTable("Input matrix", mstart, nRowsStart, m_numPlayers, 0, m_groupSize, true);
						startThread(iTask);
						nThreadsRunning++;
					}
				}
				else
				{
					nThreadsRunning++;
					if (m_cnt[iTask * 2] >= 0)
					{
						//printf("thread %d ended, %zd matrices processed\n", iTask, m_cnt[iTask * 2]);
						// thread finished
						nMatricesProc++;
						t.join();
						threadStopped(iTask);
					}
				}
				iTask++;
			}

			std::this_thread::sleep_for(std::chrono::milliseconds(5));

			if (clock() - cTime > 20000)
			{
				printThreadsStat(nMatrices, nMatricesProc, iTime, ((m_iPrintCount++) % 10) == 0);
				cTime = clock();
			}
		}
		waitAllThreadFinished();
		rTime = clock() - iTime;
		resultMatr = printThreadsStat(nMatrices, m_iMatrix, iTime, true);
		const auto str = std::format("Total time={}ms (including prep time={}ms)\n", rTime, mTime);
		printf(str.c_str());
		m_reportInfo += str;
	}

	const auto expectedResult = param(t_expectedResult);
	const auto code = expectedResult >= 0 && expectedResult != resultMatr ? 1 : 0;
	if (code)
		printfRed("*** Discrepancy Between Expected and Actual Number of Constructed Matrices: (%d != %lld)\n", expectedResult, resultMatr);

	reportEOJ(code);
	return code;
}
