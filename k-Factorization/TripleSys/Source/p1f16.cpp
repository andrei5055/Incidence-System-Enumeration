#include <iostream>
#include "TripleSys.h"
inline void unsetCell(int* cr, int* cb, tchar* pv, ctchar ir, ctchar ic, ctchar np)
{
	ctchar v = *pv;
	if (v == unset)
		return;
	if (v < ic)
		return;
	int iv = ~(1 << v);
	ASSERT(ic >= np);
	ASSERT(ir >= np);
	cb[ic] &= iv;
	cr[ir] &= iv;
	*pv = unset;
	iv = ~(1 << ic);
	cb[v] &= iv;
	cr[ir] &= iv;
	*(pv + v - ic) = unset;
}
inline void setCell(int* cr, int* cb, tchar* pv, ctchar v, ctchar ir, ctchar ic, ctchar np)
{
	ASSERT(v >= np);
	if (v < ic)
		return;
	int iv = 1 << v;
	cb[ic] |= iv;
	cr[ir] |= iv;
	*pv = v;
	pv = pv + v - ic;
	iv = 1 << ic;
	cb[v] |= iv;
	cr[ir] |= iv;
	*pv = ic;
}
void unsetLsCells(int* cr, int* cb, tchar** lsr, ctchar irs, ctchar ire, ctchar ics, ctchar np)
{
	tchar js = ics < 1 ? 1 : ics; // cb[0] not supported (no modification to first column of LS)
	tchar je = np;
	for (int i = irs; i <= ire; i++) {
		ASSERT(i >= np || np == unset || js >= np);
		tchar* p = lsr[i];
		for (int j = js; j < np; j++)
			unsetCell(cr, cb, p + j, i, j, np);
		js = 1;
	}
}
void setLsCells(int* cr, int* cb, tchar** lsr, ctchar irs, ctchar ire, ctchar np)
{
	tchar irStart = irs;
	if (irs == 0) {
		for (int j = 1; j < np; j++)
			cb[j] = 1 << j;
		irStart = 1;
	}
	for (tchar i = irStart; i <= ire; i++) {
		ASSERT(i >= np || np == unset);
		tchar* p = lsr[i];
		for (int j = 0; j < np; j++) // cb[0] not supported (no modification to first column of LS)
			setCell(cr, cb, p + j, *(p + j), i, j, np);
	}
}
void initLsTables(tchar* ls, int* cr, int* cb, tchar** lsr, ctchar nRowsExists, ctchar np)
{
	memset(cb, 0, np * sizeof(cb[0]));
	for (int i = 0; i < np; i++) {
		lsr[i] = ls + i * np;
		cr[i] = 1 << i;
		if (i >= nRowsExists) {
			memset(lsr[i], unset, np);
			setCell(cr, cb, lsr[i], i, i, 0, np);
		}
	}
	memset(cb, 0, np * sizeof(cb[0]));
	setLsCells(cr, cb, lsr, 0, nRowsExists - 1, np);
}
ctchar checkAllRowsP1f(tchar** lsr, ctchar ir, ctchar np) {
	ctchar* t0 = lsr[ir];
	//t0 = lsr[3];
	tchar kMin = np;
	for (int m = 1; m < ir; m++) {
		ctchar* t1 = lsr[m];
		tchar k = 0;
		tchar kMax = 0;
		for (int n = 2; n < np; n += 2) {
			tchar kn = t0[k];
			if (kMax < k)
				kMax = k;
			if (!(k = t1[kn])) {
				if (kMin > kMax)
					kMin = kMax;
				break;
			}
		}
	}
	return kMin;
}
void resultFromNeighbors(tchar* res, tchar** lsr, ctchar nr, ctchar nRes, ctchar np) {
	//printTable("p1f full matrix", ls, nRes, np, 0);
	/**/
	tchar* resi = res + nr * np;
	for (int i = nr; i < nRes; i++, resi += np) {
		tchar k = 0;
		for (int j = 0; j < np && k < np / 2; j++) {
			tchar v = *(lsr[i + 1] + j);
			ASSERT(v >= np);
			if (v > k && j < v) {
				resi[k * 2] = j;
				resi[k * 2 + 1] = v;
				k++;
			}
		}
	}
}
bool alldata::p1f16()
{

	m_cTime = m_iTime = clock();
	for (int i = 0; i < m_numDaysResult; i++)
		m_rowTime[i] = m_cTime - m_iTime;
	tchar df[MAX_PLAYER_NUMBER * MAX_PLAYER_NUMBER];
	tchar* res[MAX_PLAYER_NUMBER];
	tchar ls[MAX_PLAYER_NUMBER * MAX_PLAYER_NUMBER];
	tchar* lsr[MAX_PLAYER_NUMBER];
	int cb[MAX_PLAYER_NUMBER], cr[MAX_PLAYER_NUMBER];
	ctchar np = m_numPlayers;
	int nm = 0;
	int nrMax = 0;
	tchar ic = 1;
	tchar nr = iDay + 1; // number of rows in result
	tchar ir = nr + 1;   // index of row in ls to setup
	memset(df, 0, sizeof(df));
	for (int i = 0; i < np; i++)
		res[i] = result(i);
	memcpy(ls, result(), np); // fill first row of ls with 0, 1, 2, ...
	for (int i = 0; i < nr; i++) 
		u1fSetTableRow(ls + np * (i + 1), result(i), true); // fill nr rows of ls with data on the basis of result
	initLsTables(ls, cr, cb, lsr, nr + 1, np);
	printTable("lsr", ls, np, np, 0);
/*0  1   2  3   4  5   6  7   8  9
  1  0   3  2   5  4   7  6   9  8
  2  4   0  6   1  8   3  9   5  7
  3  5   7  0   9  1   8  2   6  4 */ //0 3  3 2  2 7  7 6  6 8  8 9  9 4  4 5  5 1  1 0
	
	bool bSetMode = true;
	tchar iStart = unset;
	while (1) {
		tchar* pv = lsr[ir] + ic;
		tchar ip;
		if (iStart < ic) {// && *(lsr[ir] + iStart) == ic)
			ip = iStart;
		}
		else {
			ip = iStart + 1;
			int msk = cr[ir] | cb[ic];
			for (; ip < np; ip++) {
				int bp = 1 << ip;
				if (bp & msk)
					continue;
				if (ir == ic && ip != 0) {
					ip = np;
					break;
				}
				if ((ic == 1 && ip == *lsr[ir - 1] && *(lsr[ir]) == *(lsr[ir - 1] + 1)))
					continue;
				if (ip < ic && *pv != ic)
					continue;
				setCell(cr, cb, pv, ip, ir, ic, np);
				bSetMode = true;
				break;
			}
		}
		//printTable("", lsr[ir], 1, ic + 1, 0);
		if (ip >= np)
			bSetMode = false;
		if (bSetMode) {
			ic++;
			pv++;
			iStart = *pv;
			if (ic < np) {
				continue;
			}
			else {
				ctchar irOrg = ir;
				ic = 1;
				ctchar kMin = checkAllRowsP1f(lsr, ir, np);
				if (kMin < np) {
					tchar* t0 = lsr[ir];
					ic = kMin + 1;
					for (tchar k = ic; k < np; k++) {
						unsetCell(cr, cb, t0 + k, ir, k, np);
					}
					//printTable("", lsr[ir], 1, ic, 0);
					goto notP1F;
				}
				ir++;
				if (nrMax < ir || ir >= m_numDaysResult) {
					nrMax = ir;
					//printfGreen("\nnm=%d nr=%d ", nm, nrMax);
				}

				//if (ir <= m_numDaysResult - 7) { // 3 min 6 matrices
				if (ir <= m_numDaysResult) { // 15 sec 6 matrices, 9:40 - 12, 59min - 15
					pv = lsr[ir] + ic;
					iStart = *pv;
					continue;
				}
				ctchar nRes = ir - 1;
				resultFromNeighbors(result(), lsr, nr, nRes, np);
				memcpy(neighbors(), lsr[1], np * nRes);
				/**
				for (int i = 0; i < nRes; i++) {
					for (int j = 0; j < np - 1; j += 2)
						df[i * (np / 2) + j / 2] = result(i)[j + 1] - result(i)[j];
				}
				printTable("df", df, ir - 1, np / 2, 1);
				**/
				m_lastRowWithTestedTrs = 0;
				//printTable("nb", neighbors(), nRes, np, m_groupSize);
				//printTable("rm", result(), nRes, np, m_groupSize);
				bool bRet = cnvCheckNew(0, nRes, false);
				m_lastRowWithTestedTrs = 0;
				if (bRet) {
					if (irOrg == m_numDaysResult) {
						nm++;
						m_cTime = clock();
						for (int i = nr; i < nRes; i++)
							m_rowTime[i] = m_cTime - m_iTime;
						printResultWithHistory("Result", nRes);
						//printTable("nb", neighbors(), nRes, np, vm_groupSize);
						//printTable("rm", result(), nRes, np, m_groupSize);
						printfRed("nm=%d ir=%d OK\n", nm, ir);
					}
					else {
						pv = lsr[ir] + ic;
						iStart = *pv;
						continue;
					}
				}
				else {
					ctchar iv = m_playerIndex % np;
					tchar irn = (m_playerIndex - iv) / np + 1;
					ctchar* t0 = lsr[irn];
					tchar icn = result(irn - 1)[iv];
					if (icn < t0[icn])
						icn = t0[icn];
					if (irn >= m_numDaysResult) {
						unsetLsCells(cr, cb, lsr, m_numDaysResult, m_numDaysResult, 1, np);
						irn = m_numDaysResult;
						icn = 1;
					}
					else {
						ASSERT(irn > irOrg || (irn == irOrg && icn >= np));
						unsetLsCells(cr, cb, lsr, irn, irOrg, icn, np);
					}
					ic = icn;
					ir = irn;
				}
			}
		}
	notP1F:
		bSetMode = false;
		if (ic == 1) {
			ir--;
			if (ir <= nr) {
				printfRed("\nnmatrices=%d\n", nm);
				break;
			}
			ic = np - 1;
		}
		else
			ic--;
		pv = lsr[ir] + ic;
		iStart = *pv;
		unsetCell(cr, cb, pv, ir, ic, np);
	}
	return false;
}