#include <iostream>
#include <windows.h>
#include "TripleSys.h"

void alldata::setStat()
{
	if (bg == NULL)
	{
		const auto len = m_np3 * m_nGroups * m_numDays;
		bg = new unsigned char[len];
		if (bg == NULL)
		{
			printf("no memory (%d)\n", len);
			exit(0);
		}
		memset(bg, 0, len);
	}
	if ((iPlayer % GroupSize) == (GroupSize - 1))
	{
		int n = (iDay * m_nGroups + (iPlayer / GroupSize)) * m_np3 +
			tmpPlayers[iPlayer - 2] * m_np2 +
			tmpPlayers[iPlayer - 1] * m_numPlayers +
			tmpPlayers[iPlayer];
		if (bg[n] < 255)
			bg[n]++;
	}
}

void alldata::printStat() const
{
	static int lastRep = GetTickCount();
	if (bg == NULL)
		return;

	if (GetTickCount() - lastRep > 60000)
	{
		lastRep = GetTickCount();
		printf("\n");
		for (int k = 0;k < m_numDays; k++)
		{
			for (int n = 0; n < m_nGroups; n++)
			{
				int cnt = 0;
				int j = (k * m_nGroups + n) * m_np3;
				for (int m = 0; m < m_np3; m++)
				{
					if (bg[j + m] != 0)
						cnt++;
				}
				printf(" %6d", cnt);
			}
			printf("\n");
		}
		printf("\n");
		for (int k = 0;k < m_numDays; k++)
		{
			for (int n = 0; n < m_nGroups; n++)
			{
				char cnt[2][3] = { { -1, -1, -1 }, {22, 22, 22} };
				////memcpy(cnt[0], &result[k][n * 3], 3);
				//memcpy(cnt[1], &result[k][n * 3], 3);

				int j = (k * m_nGroups + n) * m_np3;
				for (int m = 0; m < m_np3; m++)
				{
					if (bg[j + m] != 0)
					{
						char v[3];
						v[0] = m / m_np2;
						v[1] = (m / m_numPlayers) % m_numPlayers;
						v[2] = m % m_numPlayers;
						if (k == 7 && n == 0 && v[2] == 0)
							m = m;
						for (int l = 0; l < 3;l++)
						{
							if (cnt[0][l] < v[l])
								cnt[0][l] = v[l];
							if (cnt[1][l] > v[l])
								cnt[1][l] = v[l];
						}
					}
				}
				if (cnt[0][0] == -1)
					memset(cnt[0], 0, sizeof(cnt));

				for (int l = 0; l < 3;l++)
				{
					printf("%2d", cnt[1][l]);
				}
				if (memcmp(cnt[0], cnt[1], 3) == 0)
					printf("       ");
				else
				{
					printf("-");
					for (int l = 0; l < 3; l++)
					{
						printf("%2d", cnt[0][l]);
					}
				}
				printf(" ");
			}
			printf("\n");
		}
		memset(bg, 0, m_np3 * m_nGroups * m_numDays);
		printTable("Max Result", maxResult, m_numDays, m_numPlayers);
		printTable("Result", result(), m_numDays, m_numPlayers);
	}
}