#include <fstream>
#include "k-SysSupport.h"

int readTable(const std::string& fn, int nRows, int nCols, tchar** pSm, int nmax, int reservedElement, char infoSymb) {
	if (fn.length() == 0)
		return 0;
	std::ifstream mf;
	mf.open(fn, std::ios::in);
	if (!mf)
		return 0;

	const auto checkFirstElem = nRows > 0;
	int i(0), nl(0);
	auto sm = *pSm;
	while (i < nmax)
	{
		char fb;
		mf >> fb;
		if (mf.eof())
			break;

		if (fb == infoSymb) {
			auto* sms = sm;
			int iv = 0;
			int j = 0;
			while (true) {
				for (int k = 0; k < nCols; k++) {
					mf >> iv;
					if (mf.eof() || checkFirstElem && (k == 0 && iv))
					{
						printfRed("*** Error in data line %d, file: %s\n", nl, fn.c_str());
						return 0;
					}
					*(sm++) = iv;
				}
				nl++;
				mf.ignore(256, '\n');
				mf >> fb;
				++j;
				if (nRows > 0) {
					if (j == nRows)
						break;

					if (mf.eof() || fb != '\"')
					{
						printfRed("*** Error in data\n");
						sm = sms;
						goto retry;
					}
				}
				else {
					if (mf.eof() || fb != '\"') {
						mf.close();
						return j;
					}

					if (j == reservedElement)
						sm = reallocStorageMemory(pSm, nCols * (reservedElement <<= 1));
				}

			}
			i++;
			continue;
		}

	retry:
		mf.ignore(256, '\n');
		nl++;
	}
	mf.close();
	return i;
}
