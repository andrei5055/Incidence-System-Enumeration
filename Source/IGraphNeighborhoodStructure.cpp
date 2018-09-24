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
bool RunOperation(designParam *pParam, const char *pSummaryFileName, bool FirstPath);

bool IntersectionArrayIsValid(int nVertex, int k, const int *pVal, const int *pMult, int iMax, int *pUsedFlags) {
	const bool doubleFlag = pMult[iMax] == k && pVal[iMax];
	const int a = doubleFlag? pVal[iMax] + 1 : 1;

	// Theorem 4: If we do have identical elements (forming the groups), which belong to the same k blocks
	// then:
	//    a) v and k are divisible by size of the group a = (pVal[iMax] + 1)
	//    b) k / a >= the minimal number of different groups 
	if (doubleFlag) {
		if (nVertex % a || k % a)
			return false;

		int nGroups = 1;
		for (int i = 0; i < iMax; ++i) {
			if (pVal[i] % a)
				return false;

			const auto b = pVal[i] * pMult[i];
			if (!b)
				continue;

			if (b % a)
				return false;

			nGroups++;
		}

		if (a * nGroups > k)
			return false;
	}

	// Teorem 5:  if pMult[i] != 0, 
	//   a) nVertex * pVal[i] % 2 == 0
	//   b) k * [(pMult[i] * pVal[i]) / k] % 2 == 0
	// Prof: (a) Calculate the number of full bipartite graph K(2, i) for 1 <= i <= k
	//    with two vertices in one part and i vertex in another one
	//		 (b) do the same for K(2, i) for 1 <= i <= k with the fixrd vertex
	const auto iMin = pMult[0] ? 0 : 1;
	if (nVertex % 2) {
		for (int i = iMin; i <= iMax; ++i)
			if (pVal[i] % 2)
				return false;
	}

	if (k % 2) {
		for (int i = iMin; i <= iMax; ++i)
			if ((pVal[i] * pMult[i]) % 2)
				return false;
	}
//	return true;
	// We also could try to use following
	// Theorem 3: If iMin == 0 && val[0] > 1, then nVertex >= 3 * k - 2 * z,
	//            where i - max index of a(i) != 0.
	// Prof: follows from nVertex/2 - 2 * k >= k - 2 * i
	//            e1        e2 
	//             k    |    k    | v - 2k  |       
	// e1:   11111...111                    |
	// e2:   00000...000111...1111          |
	// e3:   1..110..0001..1100000???...????|
	//       e1&e3      e2&e3         m     |
	//       <------------- v ------------->|
	//

	if (!pMult[0] && pVal[0]) {
		// There are two elements (e1, e2) which do not belong to common block
		// Rewrite possible intersections into compressed array
		int idxMax = iMax - (doubleFlag ? 1 : 0);
		int *pTmp = pUsedFlags + idxMax + 1;

		// Maximal number of elements which have no common blocks with e2 
		int nUsed0 = pVal[0] - a;
		int i = 0;
		if (nUsed0 > 0) {
			// There are at least 2 elements which do not belong to any
			// of the blocks of which contains the current element e2
			// (it means that e2&e3 could be empty)
			pTmp[i++] = 0;
		}

		// Compress remaining lambda's
		for (int j = 1; j <= idxMax; ++j) {
			if (pVal[j])
				pTmp[i++] = j;
		}

		// Number of elements in compressed array
		idxMax = i;
		// Number of blocks which contain none of (e1, e2) 
		const int m = nVertex - 2 * k;

		// Flags, which will mark the indices of the intersections 
		// which coud be used as the intersections of the element e2
		memset(pUsedFlags, 0, i * sizeof(pUsedFlags[0]));

		// Loop over the number of possible common blocks of e1 and e3
		while (i--) {			
			const auto e1_e3 = pMult[pTmp[i]];
			const auto remK = k - e1_e3;  // remaining number of units which are not in e1   
			for (int j = 0; j <= i; ++j) {
				const auto def = remK - pMult[pTmp[j]];
				if (def < 0)	// Current and all following intersections 
					break;		// cannot be used

				// Check if the deficite of units in row e3 could be replenished in last m block
				if (def > m)
					continue;  // it couldn't

				// The deficit of units row e3 could be replenished 
				// in blocks which do not contain e1 or e2
				if (def && def == m) {
					// e2 and e3 have no common blocks
					// It is NOT possible to use this replenishment if we don't have
					// enough elements which do not have common blocks with e2 
					// To define that we should compare nUsed0 with a
					if (nUsed0 < a)
						continue;   // no luck

					if (!pUsedFlags[i] && (j == i || remK < pMult[pTmp[j]])) {
						// The i-th intersection was not yet used AND we do have last 
						// chance to do it here. In that case we should not only compare 
						// nUsed0 with the	# of elements, which have the same intersection
						// e1_e3 (similar to e3, none of them can have common blocks with e2),
						// but also modify nUsed0: 
						nUsed0 -= pVal[pTmp[i]];
						if (nUsed0 < 0)
							return false;  // i-th intersection could not be used
					}
				}

				// Mark both intersection as "used"
				pUsedFlags[i] = pUsedFlags[j] = 1;
			}

			if (!pUsedFlags[i])
				return false;   // We cannot find any any valid pair for current intersection e1&e3
		} 
	}

	return true;
}

int InconsistentGraphs(designParam *pParam, const char *pSummaryFileName, bool firstPath)
{
    int nVertex = pParam->v;
	if (nVertex < 10) {
		cout << "There are no inconsistent graphs for " << nVertex << " vertices";
		return 0;
	}

	const auto kMax = nVertex - 2;
	const int len = kMax + 1;
	int *mult = new int[len * 5];
	int *step = mult + len;
	int *val = step + len;
	int *buffer = val + len;
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
		int num = kIsPrime ? sum0 - sum0%step[i] : k / minDivider(k) - 1;
		cout << "K = " << k << endl;

		bool printMult = true;
		int paramIdx = 0;
		while (true) {
			if (sum0 >= num && sum1 >= num * mult[i]) {
				sum0 -= num;
				sum1 -= num * mult[i];
				val[i] = num;
				if (!sum0 && !sum1 && IntersectionArrayIsValid(nVertex, k, val, mult, iMax, buffer)) {
					// Set of parameters is constructed
					char buffer[256], *pBuf;
					if (printMult) {
						printMult = false;
						pBuf = buffer;
						pBuf += sprintf_s(pBuf, countof(buffer) - (pBuf - buffer), "     ");
						for (int j = 0; j <= iMax; ++j)
							pBuf += sprintf_s(pBuf, countof(buffer) - (pBuf - buffer), "%2d ", mult[j]);

						cout << buffer << endl;
					}

					pBuf = buffer;
					pBuf += sprintf_s(pBuf, countof(buffer) - (pBuf - buffer), "%3d: ", ++paramIdx);
					for (int j = 0; j <= iMax; ++j)
						pBuf += sprintf_s(pBuf, countof(buffer) - (pBuf - buffer), "%2d ", val[j]);

					cout << buffer << endl;

					const auto iStruct = pParam->InterStruct();
					iStruct->lambdaPtr()->resize(0);
					iStruct->lambdaAPtr()->resize(0);
					iStruct->lambdaBPtr()->resize(0);
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

					// Define possible values for 
					//   (a) mutual intersection of any two matrix rows;
					//   (b) number of elements, which belong exactly to i common blocks
					//   (c) number of elements, which belong to some block AND exactly to i common blocks 
					for (int j = 0; j <= jMax; ++j) {
						if (!val[j])
							continue;

						iStruct->lambdaPtr()->push_back(mult[j]);
						iStruct->lambdaAPtr()->push_back(val[j] / n);
						iStruct->lambdaBPtr()->push_back(val[j] * mult[j] / k / n);
					}

					if (true/*k >= 4  && val[0] == 8 && val[2] == 4*/) {
						if (!RunOperation<MATRIX_ELEMENT_TYPE>(pParam, pSummaryFileName, firstPath))
							return 0;

						firstPath = false;
					}
				}

				if (i && sum0) {
					num = sum0 - sum0 % step[--i];
					continue;
				}
			}
			else {
				if (num >= step[i]) {
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
/*

Main lists: valid parameters accordint to Teorem 4 and 5
  ! eliminated by Theorem 3
v = 10
K = 3
    0  1
1:  3  6
K = 4
    0  1  2  4
1:  4  0  4  1
2:  0  8  0  1
3:  3  0  6  0
4:  1  4  4  0
K = 5
K = 6
    2  3  4  6
1:  0  8  0  1
2:  3  0  6  0
3:  0  6  3  0
K = 7

v = 12
K = 3
    0  1
1:  5  6
K = 4
    0  1  2  4
1:  6  0  4  1
2:  2  8  0  1
3:  5  0  6  0
4:  3  4  4  0
5:  1  8  2  0
K = 5
    0  1  2  3
1:  1  0 10  0
K = 6
     0  1  2  3  4  6
 1:  3  0  0  6  0  2
 2:  0  0  9  0  0  2
 3:  4  0  0  0  6  1 !
 4:  2  0  0  8  0  1 
 5:  0  0  6  4  0  1
 6:  3  0  0  2  6  0 !
 7:  2  0  3  0  6  0 !
 8:  2  0  0  6  3  0 !
 9:  1  0  3  4  3  0 
10:  0  0  6  2  3  0 
11:  1  0  0 10  0  0
12:  0  0  3  8  0  0
K = 7
K = 8
    4  5  6  8
1:  8  0  0  3
2:  6  0  4  1
3:  2  8  0  1
K = 9
    6  7  9
1:  9  0  2

v = 14
K = 3
    0  1
1:  7  6
K = 4
    0  1  2  4
1:  8  0  4  1
2:  4  8  0  1
3:  7  0  6  0
4:  5  4  4  0
5:  3  8  2  0
6:  1 12  0  0
K = 5
    0  1  2  3
1:  3  0 10  0
K = 6
     0  1  2  3  4  6
 1:  6  0  0  0  6  1
 2:  4  0  0  8  0  1 !
 3:  0  6  0  6  0  1
 4:  2  0  6  4  0  1 
 5:  0  0 12  0  0  1
 6:  5  0  0  2  6  0
 7:  4  0  3  0  6  0
 8:  1  6  0  0  6  0
 9:  4  0  0  6  3  0
10:  3  0  3  4  3  0
11:  0  6  0  4  3  0
12:  2  0  6  2  3  0
13:  1  0  9  0  3  0
14:  3  0  0 10  0  0 !
15:  2  0  3  8  0  0 !
16:  1  0  6  6  0  0
17:  0  0  9  4  0  0
K = 7
K = 8
    2  3  4  5  6  8
1:  4  0  4  0  4  1
2:  0  8  0  0  4  1
3:  4  0  0  8  0  1
4:  0  0 12  0  0  1
K = 9
K = 10
K = 11

*/