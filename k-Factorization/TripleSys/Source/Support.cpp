#include "TripleSys.h"
#include <thread>
#define RPC_NO_WINDOWS_H
#include <sapi.h>
void speakText(LPCWSTR text)
{
	ISpVoice* pVoice = NULL;

	if (FAILED(::CoInitialize(NULL)))
	{
		printfRed("*** Error in COM init\n\7");
		exit(1);
	}
	//HRESULT hr = CoCreateInstance(CLSID_SpVoice, NULL, CLSCTX_ALL, IID_ISpVoice, (void**)&pVoice);
	HRESULT hr = CoCreateInstance(CLSID_SpVoice, NULL, CLSCTX_ALL, IID_PPV_ARGS(&pVoice));
	if (SUCCEEDED(hr))
	{
		pVoice->SetVolume(30);
		hr = pVoice->Speak(text, 0, NULL);
		pVoice->Release();
		pVoice = NULL;
	}
	::CoUninitialize();
}
void reportEOJ(int code)
{
	if (code)
		printfRed("Exit with error(%d). See error explanation above.\n\n", code);
	else
		printfGreen("End of job(%d)\n\n", code);
	if (code)
	{
		speakText(L"Oops. k task error");
		//speakText(L"<text>Oops.</text><PRON SYM = 'k eh ey - s ih s'/> task error"); //Oops. k sys task error
	}
	else
	{
		speakText(L"k task exit ok");
		//speakText(L"<PRON SYM=\"k eh ey - s ih s\"></PRON>"); //k sys task ended
	}
	printf("\7");
}
void myAssert(int code, const char* file, int line) {
	printfRed("Problem in line: %d of file: %s. Code = %d\n", line, file, code);
	assert(!code);
}
void myExit(int code)
{
	reportEOJ(code);
	exit(code);
}
CC bool alldata::initCurrentDay()
{
	iPlayer = 0;
	updateIndexPlayerMinMax();
	m_firstNotSel = 0;
	memset(selPlayers, unset, m_numPlayers);
	memset(tmpPlayers, unset, m_numPlayers);
	return true;
}
CC bool SizeParam::setLinksForOnePlayer(tchar id, tchar* lnk, tchar* p, int ip, tchar v) const
{
	const int i = ip % m_groupSize;
	p += ip;
	if (i)
	{
		auto* linkPtr = lnk + v * m_numPlayers;
		for (int j = 1; j <= i; j++)
		{
#if 0
			auto i2 = *(p - j);
			if (v == 13 && i2 == 6)
				j = j;
#endif
			if (linkPtr[*(p - j)] != unset)
				return false;
		}
		for (int j = 1; j <= i; j++)
		{
			const auto i2 = *(p - j);
			linkPtr[i2] = *(lnk + i2 * m_numPlayers + v) = id;
		}
	}
	*p = v;
	return true;
}

CC bool alldata::unsetLinksForOnePlayer(ctchar* p, int ip) const
{
	const int i = ip % m_groupSize;
	if (i == 0)
		return true;
	const auto i1 = *(p += ip);
	auto* linkPtr = links(i1);
	for (int j = 1; j <= i; j++)
	{
		const auto i2 = *(p - j);
		linkPtr[i2] = *(links(i2) + i1) = unset;
	}
	return true;
}
CC bool alldata::initPrevDay()
{
	bPrevResult = false;
	if (iDay >= 0 && iDay < m_numDays)
	{
		memset(result(iDay), 0, m_numPlayers);
	}

	iDay--;
	updateIndexPlayerMinMax();
	m_firstNotSel = 0;

	if (m_bCheckLinkV)
	{
		if (iDay < 1) //0)
			return false;
	}
	else
	{
		if (iDay < 1)  // keep first line (first day result)
		{
			iDay = -1;
			return false;
		}
	}
	auto* const pRes = result(iDay);
	memcpy(indexPlayer, pRes, m_numPlayers);
	memcpy(tmpPlayers, pRes, m_numPlayers);

	const auto ind = indexPlayer[iPlayer = iPlayerIni];
	ASSERT(ind == unset);

	for (int j = 0; j < m_numPlayers; j++)
	{
		ASSERT(*(pRes + j) == unset);
		const auto k = tmpPlayers[j];
		if (j < iPlayer)
		{
			selPlayers[k] = j;
		}
		else
		{
			const auto unsetLinkOK = unsetLinksForOnePlayer(pRes, j);
			ASSERT(!unsetLinkOK);
			tmpPlayers[j] = selPlayers[k] = unset;
			indexPlayer[j] = m_indexPlayerMin[j];
		}
	}

	indexPlayer[iPlayer] = ind + 1;
	return true;
}

CC void alldata::getPrevPlayer()
{
	ASSERT(iPlayer > m_numPlayers);
	if (iPlayer < m_numPlayers)
	    indexPlayer[iPlayer] = m_indexPlayerMin[iPlayer];
	while (--iPlayer >= 0)
	{
		const auto iPlayerNumber = tmpPlayers[iPlayer];
		const auto unsetLinkOk = unsetLinksForOnePlayer(tmpPlayers, iPlayer);
		ASSERT(!unsetLinkOk);
		tmpPlayers[iPlayer] = selPlayers[iPlayerNumber] = unset;
		if (m_firstNotSel > iPlayerNumber)
			m_firstNotSel = iPlayerNumber;

		if (iPlayer >= m_numPlayers - m_groupSize || iPlayerNumber >= m_numPlayers)
		{
			indexPlayer[iPlayer] = m_indexPlayerMin[iPlayer];
			continue;
		}

		indexPlayer[iPlayer] = iPlayerNumber + 1;
		return;
	}
}
CC void alldata::updateIndexPlayerMinMax()
{
	if (iDay > 0)
	{
		switch (m_groupSize){
			case 2:
			if (m_numPlayers >= 4)
			{
				if (m_useRowsPrecalculation == eCalculateRows && iDay == param(t_useRowsPrecalculation))
					m_indexPlayerMin[1] = m_indexPlayerMax[1] = m_secondPlayerInRow4;
				else
					m_indexPlayerMin[1] = m_indexPlayerMax[1] = param(t_CMP_Graph) ? iDay * 2 + 1 : iDay + 1;
				m_indexPlayerMax[2] = 1;
			}
			break;
			case 3:
			if (m_numPlayers >= 9)
			{
				if (m_useRowsPrecalculation == eCalculateRows && iDay == param(t_useRowsPrecalculation))
					m_indexPlayerMin[1] = m_indexPlayerMax[1] = m_secondPlayerInRow4;
				else {
					int ip1;
					if (param(t_CMP_Graph)) {
						ip1 = result(iDay - 1)[1] + 1; // must be the same for (param(t_CMP_Graph) == 0) ?
						for (; ip1 < m_numPlayers; ip1++)
							if ((ip1 % 3) && (links()[ip1] == unset))
								break;
						m_indexPlayerMin[1] = ip1;
						m_indexPlayerMax[1] = m_numPlayers - m_groupSize - (m_numDays - iDay);
					}
					else {
						ip1 = 3;
						for (; ip1 < m_numPlayers; ip1++)
							if (links()[ip1] == unset)
								break;
						m_indexPlayerMin[1] = m_indexPlayerMax[1] = ip1;
					}
				}
				m_indexPlayerMin[2] = m_indexPlayerMin[1] + 1;
				m_indexPlayerMax[3] = 1;
				m_indexPlayerMax[6] = 2;
			}
			break;
			default: {
				int ip1 = m_groupSize + iDay - 1;
				if (param(t_CMP_Graph)) {
					ip1 = result(iDay - 1)[1] + 1; // must be the same for (param(t_partiteGraph) == 0) ?
					for (; ip1 < m_numPlayers; ip1++)
						if ((ip1 % m_groupSize) && (links()[ip1] == unset))
							break;
					m_indexPlayerMin[1] = ip1;
					m_indexPlayerMax[1] = m_numPlayers - m_groupSize - (m_numDays - iDay);
				}
				else {
					if (m_numPlayers != m_groupSize * m_groupSize)
					{
						for (; ip1 < m_numPlayers; ip1++)
							if (links()[ip1] == unset)
								break;
					}
					m_indexPlayerMin[1] = m_indexPlayerMax[1] = ip1;
				}
				m_indexPlayerMin[2] = m_indexPlayerMin[1] + 1;
				for (int i = 1; i < m_groupSize; i++) {
					m_indexPlayerMin[m_groupSize * i] = m_indexPlayerMax[m_groupSize * i] = i;
				}

				break;
			}
		}
	}
	memcpy(indexPlayer, m_indexPlayerMin, m_numPlayers);
}
CC void alldata::setArraysForLastRow(int nrows)
{
	memcpy(indexPlayer, result(nrows - 1), m_numPlayers);
	memcpy(tmpPlayers, result(nrows - 1), m_numPlayers);
	for (int i = 0; i < m_numPlayers; i++)
		selPlayers[tmpPlayers[i]] = i;
	iPlayer = m_numPlayers;
}
