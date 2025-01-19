#include "TopGun.h"
#include "Table.h"
#include "data.h"

TopGun::TopGun(const kSysParam& param) : TopGunBase(param) {
	numThreads = param.val[t_numThreads];

	dNumMatrices[0] = nMatricesMax();
	mLinksSize = numPlayers() * numPlayers();
	if (!startMatrix)
	{
		printfRed("*** Not enough memory for initial %d-rows %d matrices. Exit\n", nRowsStart(), nMatricesReserved());
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
	bool bUsePm = param(t_MultiThreading) == 2 && param(t_useRowsPrecalculation);
	if (!param(t_MultiThreading))
	{
		alldata sys(*this, paramPtr());
		sys.initStartValues(ivc);// can be used to start from previous result
		resultMatr = sys.Run(1, eCalcResult, NULL, NULL, nRowsStart(), NULL, &m_reportInfo);
		transferMatrixDB(sys.matrixDB());
	}
	else
	{
		if (!readStartMatrices())
			myExit(1);

		//myTemporaryCheck();
		const auto orderMatrixMode = param(t_orderMatrices);
		if (orderMatrixMode) {
			orderMatrices(orderMatrixMode);
			printfGreen("%d 'Start Matrices' sorted\n", nMatrices);
			if (orderMatrixMode == 2) {
				TableAut Result("|Aut(M)|", m_numDays, m_numPlayers, 0, m_groupSize, true, true);
				Result.allocateBuffer(32);
				std::string ResultFile;
				createFolderAndFileName(ResultFile, paramPtr(), t_ResultFolder, nRowsStart(), "_OrderedMatrices.txt");
				for (int i = 0; i < nMatrices; i++) {
					const auto idx = m_pMatrixPerm[i];
					const auto pMatr = pntrStartMatrix() + idx * mStartMatrixSize;
					Result.setGroupOrder(m_pMatrixAutOrder[idx]);
					Result.printTable(pMatr, true, ResultFile.c_str(), false, nRowsStart());
				}
				printfGreen("They are saved to a file: \"%s\"\n", ResultFile.c_str());
				reportEOJ(0);
				return 0;
			}
		}
		const auto firstIndexOfStartMatrices = param(t_nFirstIndexOfStartMatrices);
		if (nMatrices <= firstIndexOfStartMatrices)
		{
			printfRed("*** Value of FirstIndexOfStartMatrices(%d) must be from 0 to number of 'Start Matrices'(%d). Exit\n",
				firstIndexOfStartMatrices, nMatrices);
			myExit(1);
		}
		if (!param(t_MultiThreading) == 1) {
			if (numThreads > nMatrices - firstIndexOfStartMatrices)
				numThreads = nMatrices - firstIndexOfStartMatrices;
		}
		threads.resize(numThreads);

		InitCnt(numThreads);
		int nThreadsRunning = 1;
		int nMatricesProc = m_iMatrix = firstIndexOfStartMatrices;
		mstart = startMatrix + m_iMatrix * mStartMatrixSize;
		mTime = clock() - iTime;

		printfYellow("\nMultithread Matrices Calculation started (time=%dsec)\n", mTime / 1000);

		cTime = clock();

		if (bUsePm)
		{
			int icode = 0;
			while (m_iMatrix < nMatrices) {
				alldata sys(*this, paramPtr());
				if (!sys.Run(1, eCalculateRows, mstart, mstart, nRowsStart(), NULL, &m_reportInfo)) {
					printfYellow("*** Number of pre-calculated solutions is 0 for matrix %d\n", m_iMatrix + 1);
				}
				else {
					for (int iTask = 0; iTask < numThreads; iTask++)
						startThread(iTask, eCalculateMatrices, true, sys.RowStorage());

					while (1) {

						std::this_thread::sleep_for(std::chrono::milliseconds(5));

						int iTask = 0;
						nThreadsRunning = 0;
						for (auto& t : threads)
						{
							if (threadActive[iTask])
							{
								nThreadsRunning++;
								if (m_cnt[iTask * 2] >= 0)
								{
									printf("thread %d ended, %zd matrices processed\n", iTask, m_cnt[iTask * 2]);
									// thread finished
									t.join();
									threadStopped(iTask);
								}
							}
							iTask++;
						}
						if (clock() - cTime > 20000)
						{
							printThreadsStat(nMatrices, nMatricesProc, iTime, ((m_iPrintCount++) % 10) == 0);
							cTime = clock();
						}
						if (!nThreadsRunning)
							break;
					}
					waitAllThreadFinished();

					for (int iTask = 0; iTask < numThreads; iTask++)
						threadStopped(iTask);

					nMatricesProc++;

					if (clock() - cTime > 20000)
					{
						printThreadsStat(nMatrices, nMatricesProc, iTime, ((m_iPrintCount++) % 10) == 0);
						cTime = clock();
					}
				}
				mstart += mStartMatrixSize;
				m_iMatrix += 1;
			}
		}
		else
		{
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
