#pragma once
#include "Table.h"

class TableLS : public TableAut
{
public:
	TableLS(char const* name, int nl, int nc, int ns, int np, bool makeString = false, bool outCntr = false, int outLS = 0, bool bCBMP = false) :
		TableAut(name, nl, nc, ns, np, makeString, outCntr),
		m_outLS(outLS), m_bCBMP(bCBMP) {
		m_isLSCreated = m_isSymmetricLS = m_isAtomicLS = false;
		m_latinSquares = NULL;
		if (outLS) {
			m_nColumns = (!bCBMP && np == 2) ? nc - 1 : ((!bCBMP && np == 3) ? nc : nc / np);
			m_latinSquares = new tchar[m_nColumns * m_nColumns * 6];
			m_tmpLS = new tchar[m_nColumns * m_nColumns];
		}
	}
	~TableLS() {
		if (m_latinSquares) {
			delete[] m_latinSquares;
			delete[] m_tmpLS;
		}
	}

	virtual void outLS(ctchar* c, int nl, FILE *f, bool outToScreen) {
		if (m_outLS < 1 || m_outLS > 6 || m_np > 3)
			return;

		/* 2 or 3 partite data V from each group converted to triples r, c, s.
			2-partite: r, c, s : r - even value in group / 2, c - (odd value in group - 1) / 2, s - row of group.
			3-partite: r, c, s : r = V/3 (if V%3==0), c = (V-1)/3 (if V%3==1), s=(V-2)/3 (if V%3==2)
			Then we generate up to 6 conjugates(or parastrophes) of the original Latin square.
			The six possible permutations of the coordinates(row, column, symbol) are:
			(r, c, s) : The original Latin square.
			(c, r, s) : The transpose of the Latin square(swapping rows and columns).
			(s, c, r) : A different Latin square derived from the original.
			(r, s, c) : Another derived Latin square.
			(c, s, r) : Another derived Latin square.
			(s, r, c) : Another derived Latin square.
		*/
		const char* lsName[] = { "RCS", "CRS", "SCR", "RSC", "CSR", "SRC" };
		const char startLineLS[] = {'*', ' ', '\0'};
		bool rowHamiltonian[6];
		memset(rowHamiltonian, false, 6);
		int nRows = nl;
		if (!m_bCBMP && m_np == 3)
			nRows = m_nColumns;
		else if (nRows > m_nColumns)
			nRows = m_nColumns;
		int lsSize = m_nColumns * m_nColumns;
		tchar* pls[6];
		memset(m_latinSquares, 0, lsSize * 6);
		for (int k = 0; k < 6; k++)
			pls[k] = m_latinSquares + k * lsSize;
		for (int i = 0; i < nl; i++) {
			auto* ci = c + i * m_nc;
			for (int j = 0; j < m_nc; j += m_np) {
				tchar v[3];
				if (m_bCBMP) {
					for (int m = 0; m < m_np; m++)
						v[ci[j + m] % m_np] = ci[j + m] / m_np;
				}
				else {
					for (int m = 0; m < m_np; m++)
						v[m] = ci[j + m];
				}
				if (m_np == 2)
					v[2] = i;
				tchar r = v[0];
				tchar c = v[1];
				tchar s = v[2];
				if (!m_bCBMP && m_np == 2) {
					if (r >= nRows)
						continue;
					if (c >= m_nColumns)
						c = r;
				}
				if (r >= m_nColumns || c >= m_nColumns || s >= m_nColumns)
					continue;
				int nPerm = (!m_bCBMP && m_np == 3) ? 6 : 1;
				for (int iPerm = 0; iPerm < nPerm; iPerm++) {
#define setLS(ls, a, b, c) ls[(a) * m_nColumns + (b)] = c
					setLS(pls[iPerm], r, c, s); // : The original Latin square.
					setLS(pls[(iPerm + 1) % 6], c, r, s); // : The transpose of the Latin square(swapping rows and columns).
					setLS(pls[(iPerm + 2) % 6], s, c, r); // : A different Latin square derived from the original.
					setLS(pls[(iPerm + 3) % 6], r, s, c); // : Another derived Latin square.
					setLS(pls[(iPerm + 4) % 6], c, s, r); // : Another derived Latin square.
					setLS(pls[(iPerm + 5) % 6], s, r, c); // : Another derived Latin square.
				}
				if (!m_bCBMP && m_np == 2) {
					SWAP(r, c); // create simmetrica to diagonal parts of matrices
					setLS(pls[0], r, c, s); // : The original Latin square.
					setLS(pls[1], c, r, s); // : The transpose of the Latin square(swapping rows and columns).
					setLS(pls[2], s, c, r); // : A different Latin square derived from the original.
					setLS(pls[3], r, s, c); // : Another derived Latin square.
					setLS(pls[4], c, s, r); // : Another derived Latin square.
					setLS(pls[5], s, r, c); // : Another derived Latin square.
				}
			}
		}
		if (!m_bCBMP && m_np == 3) {
			for (int k = 0; k < 6; k++) {
				for (int j = 0; j < nRows; j++)
					*(pls[k] + j * m_nColumns + j) = j;
			}
		}
		m_isLSCreated = true;
		m_isSymmetricLS = true;
		m_isAtomicLS = true;
		for (int k = 0; k < 6; k++) {
			sortLS(pls[k], m_tmpLS, nRows, m_nColumns);
			if (!(rowHamiltonian[k] = isRowHamiltonian(pls[k], m_tmpLS, nRows, m_nColumns)))
				m_isAtomicLS = false;
			if (k > 0 && MEMCMP(pls[0], pls[k], lsSize))
				m_isSymmetricLS = false;
		}
		char buffer[1024], *pBuf = buffer;
		for (int ki = m_outLS - 1; ki < m_outLS; ki++) {
			int k = (ki <= 0 || ki >= 6) ? 0 : ki;
		//for (int k = 0; k < 6; k++) {
			//if (!(m_outLS & (1 << k)))
			//	continue;
			pBuf = buffer;
			if (rowHamiltonian[k])
				SPRINTFD(pBuf, buffer, "%c %5zd(%d): row-Hamiltonian ", '*', m_cntr, k + 1);
			else
				SPRINTFD(pBuf, buffer, "%c %5zd(%d): Not row-Hamiltonian ", '*', m_cntr, k + 1);
			SPRINTFD(pBuf, buffer, "%s Latin Square%s", lsName[k], m_isAtomicLS ? "(Atomic, " : "(Not Atomic, ");
			SPRINTFD(pBuf, buffer, "%s:\n", m_isSymmetricLS ? "Totally symmetric)" : "Not Totally symmetric)");
			_printf(f, outToScreen, buffer);
			outMatrix(pls[k], nRows, m_nColumns, 0, m_ns, f, false, outToScreen, startLineLS, -1, NULL, false);
		}
	}

	virtual bool isLSCreated() const	{ return m_isLSCreated; }
	inline bool isSymmetricLS() const	{ return m_isSymmetricLS; }
	inline bool isAtomicLS() const		{ return m_isAtomicLS; }
private:
	const int m_outLS;
	const bool m_bCBMP;
	bool m_isLSCreated;
	bool m_isSymmetricLS;
	bool m_isAtomicLS;
	int m_nColumns;
	tchar* m_tmpLS;
	tchar* m_latinSquares;
};
