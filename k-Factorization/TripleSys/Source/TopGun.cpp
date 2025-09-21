#include "TopGun.h"
#include "data.h"
#include <set>
#include <filesystem>

RowDB* TopGun::m_pSecondRowsDB = NULL;

TopGun::TopGun(const kSysParam& param) : TopGunBase(param) {
	numThreads = param.val[t_numThreads];

	dNumMatrices[0] = nMatricesMax();
	mLinksSize = numPlayers() * numPlayers();
	if (!inputMatrices())
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

	m_iPrintCount = m_iMatrix = 0;
}

TopGun::~TopGun() {
	delete[] m_cntTotal;
	delete[] threadActive;
	reportEOJ(m_errCode);
}
void TopGun::deleteOldFiles()
{
	namespace fs = std::filesystem;
	std::string resultFile, resultFolder;
	std::string fName = "0.txt";
	const auto numRows = param(t_nRowsInResultMatrix) ? param(t_nRowsInResultMatrix) : m_numDays;
	createFolderAndFileName(resultFolder, paramPtr(), t_ResultFolder, numRows);
	createFolderAndFileName(resultFile, paramPtr(), t_ResultFolder, numRows, fName);
	fs::path filePath = resultFile;
	std::string fileMask = filePath.filename().string();
	const auto i0pos = fileMask.find('0');
	if (i0pos < 1)
		return;

	for (auto& entry : fs::directory_iterator(fs::path(resultFolder))) {
		filePath = entry.path();
		std::string fn = filePath.filename().string();
		if (fn.size() > 4 && !fn.compare(0, i0pos + 1, fileMask, 0, i0pos + 1) && !fn.compare(fn.size() - 4, 4, ".txt")) {
			std::error_code ec;
			bool removed = fs::remove(filePath, ec);

			if (removed) {
				// File was successfully removed
			}
			else {
				// File removal failed, check ec for details
				std::cerr << "Error removing file: " << ec.message() << std::endl;
			}
		}
	}

	//printf("%s\n", resultFile.c_str());
	return;
}
int TopGun::Run()
{
	m_errCode = 0;
	iTime = clock();
	sLongLong resultMatr = 0;
	const auto orderMatrixMode = param(t_orderMatrices);

	if (orderMatrixMode < 2 && (m_groupSize <= 3 || param(t_CBMP_Graph) > 1) && m_use2RowsCanonization) {
		if (m_pSecondRowsDB && !m_pSecondRowsDB->isValid(paramPtr())) {
			delete m_pSecondRowsDB;
			m_pSecondRowsDB = NULL;
		}

		if (!m_pSecondRowsDB) {
			alldata sys(*this, paramPtr(), 1);
			m_pSecondRowsDB = new RowDB(*paramPtr());
			// Use:
			// a) m_pSecondRowsDB->addObject(pSecondRow)
			// b) m_pSecondRowsDB->getObject(int idx)  for accessing record # idx 
			// c) m_pSecondRowsDB->numObjects()        when you need to get the number of records (second rows) in DB
			//
			resultMatr = sys.Run(1, eCalcSecondRow, m_pSecondRowsDB, NULL, NULL, 0, NULL, &m_reportInfo);
			if (resultMatr == 0) {
				printfRed("*** Cannot create second row(s) with these parameters. Exit\n");
				myExit(1);
			}
		}
	}
	if (orderMatrixMode || param(t_MultiThreading)) {
		if (readMatrices() < 0)
			myExit(1);
		const auto exploreMatrices = param(t_exploreMatrices);
		if (!param(t_nFirstIndexOfStartMatrices) && !exploreMatrices)
			deleteOldFiles();

		//myTemporaryCheck();
		if (orderMatrixMode) {
			orderAndExploreMatrices(nRowsStart(), orderMatrixMode, exploreMatrices > 1);
			if (orderMatrixMode == 2)
				return 0;
		}
		const auto firstIndexOfStartMatrices = (uint)param(t_nFirstIndexOfStartMatrices);
		const auto nMatrices = numMatrices2Process();
		if (nMatrices && nMatrices <= firstIndexOfStartMatrices)
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
		mfirst = inputMatrices();
		mstart = mfirst + m_iMatrix * inputMatrixSize();
		mTime = clock() - iTime;

		printfYellow("\nMultithread Matrices Calculation started (time=%dsec)\n", mTime / 1000);

		cTime = clock();

		const bool bUseMultiThread2 = param(t_MultiThreading) == 2 && param(t_useRowsPrecalculation);
		if (bUseMultiThread2)
		{
			int icode = 0;
			while (m_iMatrix < nMatrices) {
				alldata sys(*this, paramPtr());
				if (!sys.Run(1, eCalculateRows, m_pSecondRowsDB, mstart, mfirst, nRowsStart(), NULL, &m_reportInfo)) {
					if (param(t_printMatrices))
						printfYellow("*** Number of pre-calculated solutions is 0 for matrix %d\n", m_iMatrix + 1);
				}
				else {
					for (uint iTask = 0; iTask < numThreads; iTask++)
						startThread(iTask, m_iMatrix * numThreads + iTask + 1, eCalculateMatrices, sys.RowStorage());

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
									//printf("thread %d ended, %zd matrices processed\n", iTask, m_cnt[iTask * 2]);
									// thread finished
									t.join();
									threadStopped(iTask);
								}
							}
							iTask++;
						}
						if (!nThreadsRunning)
							break;
					}
					waitAllThreadFinished();

					for (uint iTask = 0; iTask < numThreads; iTask++)
						threadStopped(iTask);

				}
				nMatricesProc++;
				if (clock() - cTime > 20000)
				{
					printThreadsStat(nMatrices, nMatricesProc, iTime, ((m_iPrintCount++) % 10) == 0);
					cTime = clock();
				}
				mstart += inputMatrixSize();
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
							startThread(iTask, m_iMatrix + 1);
							mstart += inputMatrixSize();
							m_iMatrix += 1;
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
	else {
		deleteOldFiles();
		alldata sys(*this, paramPtr());
		sys.initStartValues(MatrixFromDatah);// can be used for testing to start from matrix selected in data.h
		resultMatr = sys.Run(1, eCalcResult, m_pSecondRowsDB, NULL, NULL, nRowsStart(), NULL, &m_reportInfo);
		transferMatrixDB(sys.matrixDB());
	}

	const auto expectedResult = param(t_expectedResult);
	m_errCode = expectedResult >= 0 && expectedResult != resultMatr ? 1 : 0;
	if (m_errCode)
		printfRed("*** Discrepancy Between Expected and Actual Number of Constructed Matrices: (%d != %lld)\n", expectedResult, resultMatr);
	else
		finalizeSemiM();
	return m_errCode;
}

#include <mutex>

std::mutex mtxLinks; // The mutex to protect the shared resource
CStorageIdx<tchar>** mpLinks = NULL;
CStorageIdx<tchar>** mShLinks = NULL;
CStorageIdx<tchar>** mShLinks2 = NULL;
int SemiPhase = 0;

void TopGun::finalizeSemiM()
{
	const auto nSemiM = mShLinks ? mShLinks[0]->numObjects() : -1;
	if (nSemiM < 0)
		return;
	SemiPhase++;
}
