#include <iostream>
#include "TripleSys.h"
void swap(char* lnk, int np, int ind1, int ind2, int ind3)
{
	for (int i = 0; i < np; i++)
	{
		char a = lnk[ind2 * np + i];
		lnk[ind2 * np + i] = lnk[ind3 * np + i];
		lnk[ind3 * np + i] = a;
	}
	for (int i = 0; i < np; i++)
	{
		char a = lnk[ind2 + np * i];
		lnk[ind2 + np * i] = lnk[ind3 + np * i];
		lnk[ind3 + np * i] = a;
	}
}
void alldata::sortLinks()
{
	return;
	int npj = m_numPlayers;
	memcpy(m_lo, links(), m_numPlayers * m_numPlayers);
	printTableColor("Original Links", m_lo, m_numPlayers, m_numPlayers);
	for (int i = 0; i < m_numPlayers; i++)
	{
		int lastCol = unset;
		for (int j = 0; j < npj; j++)
		{
			if (m_lo[i * m_numPlayers + j] >= unset)
			{
				if (lastCol == unset)
				    lastCol = j;
			}
			else if (i != j && lastCol != unset)
			{
				swap(m_lo, m_numPlayers, i, lastCol, j);
				j = lastCol;
				lastCol = unset;
			}
		} 
		for (int j = 0; j < npj; j++)
		{
			if (m_lo[i * m_numPlayers + j] != unset)
			{
				npj = j;
				break;
			}
		}
	}
    printTableColor("Sorted Links", m_lo, m_numPlayers, m_numPlayers);
	convertLinksToResult(m_lo);
	printTable("Result from Sorted Links", m_co, m_numDays, m_numPlayers);
}
