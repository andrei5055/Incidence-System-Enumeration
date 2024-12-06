#include <fstream>
#include "k-SysSupport.h"

int readTable(const std::string& fn, int nRows, int nCols, tchar** pSm, int nmax, int& reservedElement, uint** ppGroupOrders, char infoSymb) {
	if (fn.length() == 0)
		return 0;
	std::ifstream mf;
	mf.open(fn, std::ios::in);
	if (!mf)
		return 0;

	const auto checkFirstElem = nRows > 0;
	int i(0), nl(0);
	auto sm = *pSm;
#if 0
	std::string line;
	size_t start, end;
	while (i < nmax && !mf.eof()) {
		int j = 0;
		while (getline(mf, line)) {
			nl++;
			trim(line);
			if (line[0] == infoSymb) {
				if (!j++) {
					if (i++ == reservedElement)
						sm = reallocStorageMemory(pSm, nCols * nRows * (reservedElement <<= 1));
				}

				int iv;
				start = end = 0;
				std::string dl("\" ");
				for (int k = 0; k < nCols; k++) {
					start = line.find_first_not_of(dl, end);
					end = line.find(" ", start);
					iv = stoi(line.substr(start, end - start));
					if (checkFirstElem && (k == 0 && iv)) {
						printfRed("*** Error in data line %d, file: %s\n", nl, fn.c_str());
						return 0;
					}
					*(sm++) = iv;
				}

				if (j == nRows)
					break;
			}
			else {
				if (!j) {
					start = line.find("|Aut(M)| = ", 0);
					if (start != -1) {
						start += 11;
						end = line.find(",", start);
						auto order = stoi(line.substr(start, end - start));
						start += 0;
					}
				}
				else {
					printfRed("*** Error in data\n");
					return 0;
				}
			}
		}
	}
#else
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
#endif
	mf.close();
	return i;
}
