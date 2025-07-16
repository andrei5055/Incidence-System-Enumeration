#include <fstream>
#include "k-SysSupport.h"

int readTable(const std::string& fn, int nRows, int nCols, int nmax, int nUsed, tchar** ppSm, 
	int& reservedElement, int nMatricesMax, CMatrixInfo * pMatrixInfo, char infoSymb) {
	if (fn.length() == 0)
		return 0;

	std::ifstream mf;
	mf.open(fn, std::ios::in);
	if (!mf)
		return 0;

	const auto checkFirstElem = nRows > 0;
	const auto lenMatr = nCols * nRows;
	tchar* sm = *ppSm + lenMatr * nUsed;
	auto** ppGroupOrders = pMatrixInfo ? pMatrixInfo->groupOrdersPntr() : NULL;
	auto pGroupOrders = ppGroupOrders ? *ppGroupOrders + nUsed : NULL;
	std::string*** ppInfo[2] = { NULL };
	std::string** pInfo[2] = { NULL };
	if (ppGroupOrders) {
		for (auto k = sizeof(pInfo) / sizeof(pInfo[0]); k--;) {
			ppInfo[k] = k? pMatrixInfo->groupInfoPntr() : pMatrixInfo->cycleInfoPntr();
			if (ppInfo[k])
				pInfo[k] = *ppInfo[k] + nUsed;
		}
	}

	auto** ppCycleInfo = ppGroupOrders ? pMatrixInfo->cycleInfoPntr() : NULL;
	std::string line;
	size_t start, end;
	int iPrev = 0;
	int i(0), nl(0);
	while (i < nmax && !mf.eof()) {
		int j = 0;
		while (getline(mf, line)) {
			nl++;
			if (!j && i + nUsed >= reservedElement) {
				const auto prevReserved = reservedElement;
				reservedElement <<= 1;
				if (reservedElement > nMatricesMax)
					reservedElement = nMatricesMax;
				if (i + nUsed == reservedElement) {
					mf.close();
					return i;
				}
				if (!(sm = reallocStorageMemory(ppSm, lenMatr * reservedElement, lenMatr * prevReserved))) {
					printfRed("*** Failed to allocate memory for %d matrices while reading the file: %s\n", reservedElement, fn.c_str());
					mf.close();
					return 0;
				}

				if (pGroupOrders) {
					if (!reallocStorageMemory(ppGroupOrders, reservedElement, prevReserved)) {
						printfRed("*** Failed to allocate memory for %d group orders while reading the file: %s\n", reservedElement, fn.c_str());
						mf.close();
						return 0;
					}

					pGroupOrders = *ppGroupOrders + nUsed;
				}

				for (int k = sizeof(pInfo)/sizeof(pInfo[0]); k--;) {
					if (pInfo[k]) {
						if (!reallocStorageMemory(ppInfo[k], reservedElement, prevReserved, true)) {
							printfRed("*** Failed to allocate memory for %d elements of %s info while reading the file: %s\n", reservedElement, (k? "groups" : "cycles"), fn.c_str());
							mf.close();
							return 0;
						}

						pInfo[k] = *ppInfo[k] + nUsed;
					}
				}
			}
			trim(line);
			if (line.empty())
				continue;

			if (line[0] == infoSymb) {
				int iv;
				start = end = 0;
				std::string dl("\" ");
				for (int k = 0; k < nCols; k++) {
					start = line.find_first_not_of(dl, end);
					end = line.find(" ", start);
					iv = stoi(line.substr(start, end - start));
					if (checkFirstElem && (k == 0 && iv)) {
						printfRed("*** Error in data line %d, file: %s\n", nl, fn.c_str());
						mf.close();
						return 0;
					}
					*(sm++) = iv;
				}

				if (++j == nRows) {
					iPrev = i++;
					j = 0;
				}
			}
			else {
				if (!j) {
					if (pGroupOrders) {
						// Collecting the group order information
						start = line.find(AUT, 0);
						if (start != -1) {
							start += strlen(AUT);
							end = line.find(",", start);
							pGroupOrders[i] = stoi(line.substr(start, end - start));
							if (pInfo[0]) {
								// Collecting the cycle information
								line = line.substr(end + 1);
								ltrim(line);
								if (!pInfo[0][i])
									pInfo[0][i] = new std::string(line);
								else
									*pInfo[0][i] = line;
							}
						}
						else {
							// Group's Info found
							if (pInfo[1]) {
								// We are ready to store it...
								if (line.find("Orbits") == 0)
									line = "\n" + line;

								if (!pInfo[1][iPrev])
									pInfo[1][iPrev] = new std::string(line);
								else
									*pInfo[1][iPrev] += "\n" + line;
							}
						}
					}
				}
				else {
					printfRed("*** Error in data\n");
					mf.close();
					return 0;
				}
			}
		}
	}

	if (pGroupOrders && pInfo[1][iPrev])
		*pInfo[1][iPrev] += "\n";

	mf.close();
	return i;
}
