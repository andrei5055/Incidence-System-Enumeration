#include <thread>
#include "TripleSys.h"
#include "main.h"
#include <vector>
void RunThread(int threadNumber, int iMode, int improveResult, 
	char* mstart0, char* mstart, char* mstop, int nRows, int mStep, double* pcnt, bool bPrint)
{
	alldata sys(nPlayers);
	sys.Run(threadNumber, iMode, improveResult, mstart0, mstart, mstop, nRows, mStep, pcnt, bPrint);
}
int main()
{
	double cnt[NThreads*2];
	double cntTotal[NThreads*2];
	double dNumMatrices[2], dNumMatricesMax;
	bool threadActive[NThreads];
	int mStep = StepForStartMatrices;
	int nRows = NRowsInStartMatrix;
	clock_t cTime, rTime, mTime, iTime = clock();
//	_CrtSetBreakAlloc(174);  // Put here the memory allocation number you want to stop at.
	dNumMatrices[0] = dNumMatricesMax = MaxNumberOfStartMatrices;
	std::cout << "TripleSys 12.2\n";

	if (nRows == 0)
	{
		alldata sys(nPlayers);
		sys.initStartValues(ivc);// can be used to start from previous result
		sys.Run(0, eCalcResult, ImproveResults, (char*)0, (char*)0, (char*)0, 0, 0, NULL, true);
	}
	else
	{
		RunThread(0, eCalcNumberOfMatrices, ImproveResults, (char*)0, (char*)0, (char*)0, nRows, 0, dNumMatrices, false);
		if (dNumMatrices[0] > dNumMatricesMax)
		{
			printf("Number of initial %d-rows matrices (%.0f) > limit (%.0f). Exit\n", nRows, dNumMatrices[0], dNumMatricesMax);
			exit(0);
		}
		printf("Number of initial %d-rows matrices=%.0f, limit=%.0f\n", nRows, dNumMatrices[0], dNumMatricesMax);
		//mStep = (int)((dNumMatrices[0] + NThreads - 1) / NThreads);
		//mStep = mStep / 2;
		//if (mStep < 1)
		mStep = StepForStartMatrices;
		double allMatricesSize = (dNumMatrices[0] / mStep + 1) * nPlayers * nRows;
		if (allMatricesSize > MaxMemoryForStartMatrices)
		{
			printf("Memory needed for initial %d-rows matrices (%.0f) > limit (%.0f). Exit\n", 
				nRows, allMatricesSize, MaxMemoryForStartMatrices);
			exit(0);
		}
		printf("Memory needed for initial %d-rows matrices=%.0f, limit=%.0f\n",
			nRows, allMatricesSize, MaxMemoryForStartMatrices);
		int iAllMatricesSize = (int)((dNumMatrices[0] + mStep - 1) / mStep) * nPlayers * nRows;
		char* startMatrix = new char[iAllMatricesSize];
		char* endMatrix = new char[iAllMatricesSize];

		RunThread(0, eCalcStartStop, ImproveResults, (char*)0, startMatrix, endMatrix, nRows, mStep, dNumMatrices, false);
		mTime = clock() - iTime;
		printf("%.0f %d-rows matrices needed\n", dNumMatrices[0], nRows);
		/**/
		memset(cntTotal, 0, sizeof(cntTotal));
		memset(cnt, 0, sizeof(cnt));
		memset(threadActive, false, sizeof(threadActive));
		int numThreads = NThreads;
		if ((double)numThreads > dNumMatrices[0])
			numThreads = (int)dNumMatrices[0];
		std::vector<std::thread> threads(numThreads);
		int nb = 0;
		cTime = clock();
		int iPrintCount = 0;
		for (double dj = 0; dj < dNumMatrices[0]; dj += mStep, nb++) {
			while (1)
			{
				int iTask = 0;
				for (auto& t : threads) {
					if (!threadActive[iTask])
					{
						cnt[iTask * 2] = -1;
						size_t ind = nb * nPlayers * nRows;
						threads[iTask] = std::thread{ RunThread, iTask + 1, eCalcResult, ImproveResults,
							startMatrix + ind, startMatrix + ind, endMatrix + ind, nRows, 0, &cnt[iTask * 2], false };
						threadActive[iTask] = true;
						break;
					}
					else if (cnt[iTask * 2] >= 0)
					{
						t.join();
						cntTotal[iTask * 2] += cnt[iTask * 2];
						cntTotal[iTask * 2 + 1] += cnt[iTask * 2 + 1];
						cnt[iTask * 2 + 1] = cnt[iTask * 2 + 1] = 0;
						threadActive[iTask] = false;
						printf("thread %d ended\n", iTask + 1);
					}
					iTask++;
				}
				if (iTask < numThreads)
					break;
				std::this_thread::sleep_for(std::chrono::milliseconds(5));

				if (clock() - cTime > 20000)
				{
					printThreadsStat(cntTotal, cnt, dNumMatrices[0], nRows, mStep, numThreads, iTime, ((iPrintCount++) % 10) == 0);
					cTime = clock();
				}
			}
		}
		int i = 0;
		for (auto& t : threads) {
			if (threadActive[i])
			{
				t.join();
				if (cnt[i * 2] < 0)
					printf("\n*** Internal error in thread %d (%.0f) ***\n", i + 1, cnt[i * 2]);
			}
			i++;
		}
		rTime = clock() - iTime;
		printThreadsStat(cntTotal, cnt, dNumMatrices[0], nRows, mStep, numThreads, iTime, true);
		printf("Total time=%d (include prep time=%d)\n", rTime, mTime);
	}
	cout << "\7" << endl; // play sound
	return 0;
}
