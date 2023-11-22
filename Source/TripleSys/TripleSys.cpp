#include <iostream>
#include "TripleSys.h"

void alldata::Init() {
	memset(result, 0, sizeof(result));
	memset(links, 0, sizeof(links));
	memset(bklinks, 0, sizeof(bklinks));
	memset(devCounts, 0, sizeof(devCounts));
	memset(selPlayers, 0, sizeof(selPlayers));
	memset(tmpPlayers, 0, sizeof(tmpPlayers));
	memset(indexPlayer, 0, sizeof(indexPlayer));
	memset(&preSet, unset, sizeof(preSet));
	maxDays = -1;
	maxDaysPlayers = -1;
	nLoops = 0;
	maxPlayers = 0;
	noResult = false;
	adjust0(this);
	iDay = 0;
}

bool alldata::Run(int numPlayers) {
	if (numPlayers != nPlayers) {
		std::cout << "\nProgram needs to be recompiled for nPlayers = " << numPlayers << "\n";
		return true;
	}

	std::cout << "\n";

	Init();

	bool iPrevResult = false;
	while (nLoops < 8)
	{
		if (nLoops > 0)
		{
			iDay = nDays - 1;
			iPrevResult = true;
		}
		else
		{
			iDay = 0;
			iPrevResult = false;
		}
		//int iStartDay = iDay;
		while (iDay < nDays)
		{
		runday:
			if (iDay == 0)
				iDay = iDay;

			if (iPrevResult)
			{
				initPrevDayProcess(this);
				iPrevResult = false;
			}
			else
				initCurrentDayProcess(this);
			while (iPlayer < nPlayers)
			{
				int i = 0;

				if (iPlayer < 0)
				{
					int k = -1;
					if (iDay <= 0)
					{
						noResult = true;
						goto finish;
					}
					iDay--;
					iPrevResult = true;
					goto runday;
				}

				if (iPlayer < nPlayers)
				{
					if ((i = preSet[iDay][iPlayer]) == unset)
					{
						i = indexPlayer[iPlayer];
						/*

						if (iDay == 0 && iPlayer == 4 && i <= iPrevResult)
							i = iPrevResult + 1;
							*/

							//if (iDay == nDays - 1 && iPlayer == (nPlayers - GroupSize) && i <= iPrevResult)
							//	i = iPrevResult + 1;

						int m1 = (iPlayer % GroupSize) == 0 ? GroupSize : 1;
						//int m1 = (iPlayer % GroupSize) == 0 ? 0 : 1;
						if (iPlayer >= m1 && m1 > 0)
						{
							if (i <= tmpPlayers[iPlayer - m1])
								i = tmpPlayers[iPlayer - m1] + 1;
						}
						for (; i < nPlayers; i++)
						{
							//if (iPlayer == 6)
							//	i = i;
							if (selPlayers[i] != unset)
								continue;

							//checkbmask(selPlayers, bmask);
							/*
							if (iDay > 1 && iPlayer > 8 && i <= tmpPlayers[iPlayer-GroupSize])
							{
								i = nPlayers;
								break;
							}
							*/
							if (m1 == GroupSize && (bmask ^ (bmask - 1)) <= (1 << i))
							{
								i = nPlayers;
								break;
							}

							tmpPlayers[iPlayer] = i;
							if (setLinksAndDevCounts(this, tmpPlayers, iPlayer, 1))
							{
								//if (links[1][5] != unset)
								//	i = i;
								break;
							}
							tmpPlayers[iPlayer] = unset;
							//memset(indexPlayer + iPlayer + 1, 0, sizeof(indexPlayer) - iPlayer - 1);
						}
					}
				}
				if (i >= nPlayers)
				{

					if (iPlayer >= nPlayers)
						abort();
					indexPlayer[iPlayer] = 0;
					//for (int k = iPlayer; k < nPlayers; k++)
					//	indexPlayer[k] = 0;

					for (iPlayer--; iPlayer >= 0; iPlayer--)
					{
						if (preSet[iDay][iPlayer] != unset)
							continue;
						if (!setLinksAndDevCounts(this, tmpPlayers, iPlayer, unset))
							abort();
						//
						// memcpy(&links, &bklinks[iPlayer], sizeof(links));
						if (tmpPlayers[iPlayer] < 0 || tmpPlayers[iPlayer] >= nPlayers)
							abort();
						i = tmpPlayers[iPlayer];
						//checkbmask(selPlayers, bmask);
						selPlayers[i] = unset;
						bmask = bmask | (1 << i);
						//checkbmask(selPlayers, bmask);
						tmpPlayers[iPlayer] = unset;
						result[iDay][iPlayer] = unset;
						if (iPlayer >= nPlayers - GroupSize)
							continue;
						i++;
						indexPlayer[iPlayer] = i;
						break;
					}
					continue;
				}
				//if (iPlayer == 4 && iDay == 2)
				//	i = i;

				{
					//c1++;
					if (iDay == 88)
					{
						static int c = 0;
						c++;
						//if (c <= 100)
						{
							static FILE* fd = NULL;
							if (fd == NULL)
								fopen_s(&fd, "testnew.txt", "w");
							if (fd)
							{
								fprintf(fd, "%2d %2d %2d %2d  ", c, iDay, iPlayer, i);
								if ((c % 5) == 0)
									fprintf(fd, "\n");
							}
						}
					}
				}
				indexPlayer[iPlayer] = i;
				tmpPlayers[iPlayer] = i;
				result[iDay][iPlayer] = i;
				if (iPlayer + 1 < nPlayers)
					indexPlayer[iPlayer + 1] = 0;
				//if ((iPlayer%3)>0&&setLinksAndDevCounts(&sys, tmpPlayers, iPlayer, 1))
				//	i=i;

				//checkbmask(selPlayers, bmask);
				bmask = bmask & (~(1 << i));
				selPlayers[i] = iPlayer;
				//checkbmask(selPlayers, bmask);
				iPlayer++;
				//memcpy(&bklinks[iPlayer], &links, sizeof(links));
			}
			//				printf("%d\n", iDay);
						//memcpy(&bklinks[0], &links, sizeof(links));
						//memcpy(&bklinks[1], &links, sizeof(links));
			memcpy(result[iDay], tmpPlayers, sizeof(tmpPlayers));
			if (noResult)
				break;
			//			printResult(&sys);

			memcpy(result[iDay], tmpPlayers, sizeof(tmpPlayers));
			if (maxDays < iDay)
			{
				maxDays = iDay;
				memcpy(&maxResult, &result, sizeof(maxResult));
			}

			int iCheck = check1(result, iDay + 1);
			if (iCheck == 0)
			{
				// get new matrix
				iPrevResult = 1;
				continue;
			}
			iDay++;
		}
		//check result
	finish:
		nLoops++;
		printf("result %d\n", nLoops);
		printResult(this);
		if (noResult)
			break;
		//iPrevResult = result[6][4];
		//if (iDay >= nDays - 1 && iPlayer >= nPlayers - 1)
		//	break;
	}

	return iPrevResult;
}
