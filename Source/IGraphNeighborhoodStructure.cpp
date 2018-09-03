#include "DataTypes.h"
#include "matrix.h"
#include <iostream> 
#include <stdio.h>
#include <numeric>

using namespace std;

static unsigned GCD(unsigned u, unsigned v) {
	while (v != 0) {
		unsigned r = u % v;
		u = v;
		v = r;
	}
	return u;
}

int primeNumb[] = { 2, 3, 5, 7, 11, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53 };

static int minDivider(int k) {
	const int iMax = countof(primeNumb);
	int i = -1;
	while (++i < iMax) {
		if (!(k % primeNumb[i]))
			return primeNumb[i];
	}

	return -1;
}

static bool isPrime(int k) {
	const int iMax = countof(primeNumb);
	int i = -1;
	while (++i < iMax && primeNumb[i] < k);

	return k == primeNumb[i];
}

template <class T>
bool RunOperation(designRaram *pParam, const char *pSummaryFileName, bool FirstPath);

int InconsistentGraphs(designRaram *pParam, const char *pSummaryFileName, bool firstPath)
{
    int nVertex = pParam->v;
	if (nVertex < 10) {
		cout << "There are no inconsistent graphs for " << nVertex << " vertices";
		return 0;
	}

	const auto kMax = nVertex - 2;
	const int len = kMax + 1;
	int *mult = new int[len * 3];
	int *step = mult + len;
	int *val = step + len;
	for (int k = 3; k < kMax; ++k) {
		// Theorem 1: k % (a(k) + 1) == 0
		//    Consequence 1: max(a(k)) == k / min(divider(k)) - 1;
		//    Consequence 2: if is prime number, a(k) == 0;
		// Theorem 2: a(k-1) == 0,  
		//    If NOT, a(k-1) = k, and the bipartite will be vertex-transitive

		int aMin = 2 * k - nVertex;
		if (aMin < 0)
			aMin = 0;

		const bool kIsPrime = isPrime(k);
		int aMax = k - (kIsPrime ? 2 : 1);
		const int iMax = aMax - aMin;
		for (int i = 0; i <= iMax; ++i) {
			mult[i] = aMin;
			step[i] = aMin ? k / GCD(k, aMin) : 1;
			val[i] = 0;
			if (++aMin == k - 1)
				aMin++; // We skip a(k - 1) because of the Theorem 2           
		}

		int sum0 = nVertex - 1;
		int sum1 = k * (k - 1);
		int i = iMax;
		int num = kIsPrime ? sum0 - sum0%step[i] : minDivider(k) - 1;
		cout << "K = " << k << endl;

		bool printMult = true;
		while (true) {
			if (sum0 >= num && sum1 >= num * mult[i]) {
				sum0 -= num;
				sum1 -= num * mult[i];
				val[i] = num;
				if (!sum0 && !sum1) {
					// Set of parameters is constructed
					// WE also could try to use following
					// Theorem 3: If aMin == 0 && val[0] > 1, then nVertex/2 >= 3 * k - 2 * b,
					//            where b - max index of a(i) != 0.
					// Prof: follows from nVertex/2 - 2 * k >= k - 2 * b
					char buffer[256], *pBuf;
					if (printMult) {
						printMult = false;
						pBuf = buffer;
						for (int j = 0; j <= iMax; ++j)
							pBuf += sprintf_s(pBuf, countof(buffer) - (pBuf - buffer), "%2d ", mult[j]);

						cout << buffer << endl;
					}

					pBuf = buffer;
					for (int j = 0; j <= iMax; ++j)
						pBuf += sprintf_s(pBuf, countof(buffer) - (pBuf - buffer), "%2d ", val[j]);

					cout << buffer << endl;

					pParam->lambda.resize(0);
					pParam->r = k;
					int jMax = iMax;
					int n = val[jMax];
					if (mult[iMax] == k && n) {
						++n;
						--jMax;
					}
					else
						n = 1;

					pParam->k = k / n;
					pParam->v = nVertex / n;

					for (int j = 0; j <= jMax; ++j) {
						if (val[j])
							pParam->lambda.push_back(mult[j]);
					}

					if (!RunOperation<MATRIX_ELEMENT_TYPE>(pParam, pSummaryFileName, firstPath))
						return 0;

					firstPath = false;
				}

				if (i && sum0) {
					num = sum0 - sum0 % step[--i];
					continue;
				}
			}
			else {
				if (num > step[i]) {
					num -= step[i];
					continue;
				}

				i++;
			}

			do {
				const int numb = val[i] ? (i ? step[i] : val[i]) : 0;
				if (numb) {
					sum0 += numb;
					sum1 += numb * mult[i];
					val[i] -= numb;
					if (i)
						break;
				}
			} while (++i <= iMax);

			if (i > iMax)
				break;

			num = sum0 - sum0 % step[--i];
		}
	}

	return 1;
}