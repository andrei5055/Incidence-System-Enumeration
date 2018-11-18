#include "DataTypes.h"
#include <numeric>
#include <cmath>

#define COUT(buf)	cout << buf << endl;

using namespace std;
static int lll = 0;
static unsigned GCD(unsigned u, unsigned v) {
	while (v != 0) {
		unsigned r = u % v;
		u = v;
		v = r;
	}
	return u;
}

int primeNumb[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53 };

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

		assert(pTmp - pUsedFlags + i <= lll);
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
				const auto e2_e3 = pMult[pTmp[j]];
				const auto def = remK - e2_e3;
				if (def < 0)	// Current and all following intersections 
					break;		// cannot be used

				// Check if the deficite of units in row e3 could be replenished in last m block
				if (def > m)
					continue;  // it couldn't

				if (!e2_e3) {
					// The deficit of units row e3 could be replenished 
					// in blocks which do not contain e1 or e2
					if (def && def == m) {
						// e2 and e3 have no common blocks
						// It is NOT possible to use this replenishment if we don't have
						// enough elements which do not have common blocks with e2 
						// To define that we should compare nUsed0 with a
						if (nUsed0 < a)
							continue;   // no luck

						if (!pUsedFlags[i] && (j == i || remK < pMult[pTmp[j+1]])) {
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

static size_t calcSum(const std::vector<int> &lambdaA, const std::vector<int> &lambda) {
	size_t sum = 0;
	for (auto i = lambdaA.size(); i--;)
		sum += lambdaA[i] * lambda[i] * lambda[i];

	return sum;
}

#define NEW_CondB	1
#if NEW_CondB
bool checkCondB(const CInterStruct *iStruct1, const CInterStruct *iStruct2, int *indices, int vMinus2K, int k)
{
	// Skip index which corresponds to the presence of duplicated elements
	const auto &lambdaA = iStruct1->lambdaA();
	const auto &lambdaB = iStruct1->lambdaB();
	auto size = lambdaA.size();
	int adj = lambdaB[size - 1];
	if (lambdaA[size - 1] == adj++)
		size--;
	else
		adj = 1;

	// Calculate minimal value of the counterpart's lambda (intersection with the first block)
	int idx = 0; 
	int lambdaMin = adj;
	for (auto i = size; i--;) {
		const int n = (lambdaB[i] << 1) - lambdaA[i];
		if (n <= 0)
			continue;

		indices[idx++] = static_cast<int>(i);
		if (adj > 1)
			lambdaMin += adj - n % adj;
		else
			lambdaMin += n;
	}

	const auto &lambda = iStruct2->lambda();
	const auto jMax = lambda.size();
	size_t j = lambda[0]? 0 : 1;
	if (lambda[j] < lambdaMin) {
//		cout << "CCC:  lambda[" << j << "] = " << lambda[j] << "  lambdaMin = " << lambdaMin << endl;
		return false;
	}
/*
	if (idx && lambda[j] == lambdaMin) {
		cout << "BBB: lambdaMin = " << lambdaMin << " adj = " << adj << "  ";
		if (iStruct1 == iStruct2)
			cout << "+++";
		else
			cout << "---";

		for (int i = idx; i--;) {
			const auto j = indices[i];
			cout << "  (" << iStruct1->lambda()[j] << ", " << iStruct1->lambdaA()[j] << ", " << iStruct1->lambdaB()[j] << ")";
		}

		cout << endl;
	}
*/
	if (adj == 1)
		return true;

	while (j < jMax) {
		const auto lambdaTmp = lambda[j++];
		if (lambdaTmp % adj) 
			return false;
	}

	return true;
}
#else
bool checkCondB(const CInterStruct *iStruct1, const CInterStruct *iStruct2, int *indices, int vMinus2K, int k)
{
	// Skip index which corresponds to the presence of duplicated elements
	const auto &lambdaA = iStruct1->lambdaA();
	const auto &lambdaB = iStruct1->lambdaB();
	auto size = lambdaA.size();
	int adj = lambdaB[size - 1];
	if (lambdaA[size - 1] == adj++)
		size--;
	else
		adj = 1;

	const auto &lambda = iStruct2->lambda();
	for (auto j : lambda) {
		if (!j)
			continue;
		
		if (j % adj)
			return false;

		const auto jm = j - adj;
		for (auto i = size; i--;) {
			const int n = lambdaA[i] - (lambdaB[i] << 1) + jm;
			if (n < 0) {
				cout << "CCC:  lambda[x] = " << j << "  i = " << i << "  jm = " << jm << "(" << lambdaA[i] << ", " << lambdaB[i] << ")" << endl;
				return false;
			}

			if (n > 0)
				continue;

			// When for some index i we are here, stronger condition 
			// (without addition of jm) could be checked for remaining indices
			for (auto l = size; l--;) {
				if (l == i)
					continue;

				if (lambdaA[l] - (lambdaB[l] << 1) < 0) {
					cout << "Hura!" << endl;
					return false;
				}
			}
		}
	}

	if (!lambdaB[0] && vMinus2K < 0) {
		// In that case a(0) > 0
		const auto j = lambda[0] ? lambda[0] : lambda[1];
		if (vMinus2K + j < lambdaA[0])
			return false;  // No place for a(0) elements 
	}

	return true;
}
#endif

void CheckIntersections(CInterStruct *iStructA, int *pIndices, int k, int v)
{
	// Function which check the conditions from Theorem 5.4 from Andrei Ivanov thesis
	// (a) S(i^2 * a(i)) is a constant which is the same for transposed matrix
	// (b) a(i) - 2*b(i) + j >= 1    (for any i and any j which A(j) != 0

	// Save pointer for possible deletion of the intersection parameters
	const int vMinus2K = v - 2 * k;
	CInterStruct *iStructBase = iStructA;
	CInterStruct *iStructPrev = iStructBase;
	while (iStructA = iStructPrev->getNext()) {
		const auto sum1 = calcSum(iStructA->lambdaA(), iStructA->lambda());
		bool flag = checkCondB(iStructA, iStructA, pIndices, vMinus2K, k);
		CInterStruct *iStructB = iStructA;
		do {
			if (flag || iStructB != iStructA && 
				        sum1 == calcSum(iStructB->lambdaA(), iStructB->lambda()) &&
						checkCondB(iStructA, iStructB, pIndices, vMinus2K, k) &&
						checkCondB(iStructB, iStructA, pIndices, vMinus2K, k)) {
				if (!iStructA->Counterparts())
					iStructA->InitCounterparts();

				iStructA->Counterparts()->push_back(iStructB);
			}

			flag = false;
		} while (iStructB = iStructB->getNext());

		if (!iStructA->isValid()) {
			// This structure is not a valid structure
			// Check if it is in some other structure's list of counterparts
			CInterStruct *iStructTmp = iStructBase;
			while ((iStructTmp = iStructTmp->getNext()) != iStructA) {
				const auto pntr = iStructTmp->Counterparts();
				if (!pntr)
					continue;

				int idx = static_cast<int>(pntr->size());
				while (idx-- && (*pntr)[idx] != iStructA);
				if (idx >= 0)
					break;
			}

			if (iStructTmp == iStructA) {
				// It is not. We can delete it
				iStructPrev->setNext(iStructA->getNext());
				delete iStructA;
				continue;
			}
		}

		iStructPrev = iStructA;
	}
}

void printSolution(bool &printMult, const CInterStruct *iStruct, const int *mult, int iMax, int paramIdx, int *val)
{
	char buffer[256], *pBuf;
	if (printMult) {
		printMult = false;
		pBuf = buffer;
		pBuf += sprintf_s(pBuf, countof(buffer) - (pBuf - buffer), "     ");
		for (int j = 0; j <= iMax; ++j)
			pBuf += sprintf_s(pBuf, countof(buffer) - (pBuf - buffer), "%2d ", mult[j]);

		pBuf += sprintf_s(pBuf, countof(buffer) - (pBuf - buffer), " Sum(a(i)*i^2):");
		COUT(buffer);
	}

	// Restore full solution of the system of equations
	const auto pMult = iStruct->lambda().data();
	auto const pVal = iStruct->lambdaA().data();
	memset(val, 0, sizeof(val[0]) * (iMax + 1));
	int k = 0;
	const auto kMax = iStruct->lambda().size();
	for (int i = 0; i <= iMax; ++i) {
		if (mult[i] != pMult[k])
			continue;

		val[i] = pVal[k];
		if (++k == kMax)
			break;
	}

	pBuf = buffer;
	pBuf += sprintf_s(pBuf, countof(buffer) - (pBuf - buffer), "%3d: ", paramIdx);
	int sumX = 0;
	for (int i = 0; i <= iMax; ++i) {
		sumX += mult[i] * mult[i] * val[i];
		pBuf += sprintf_s(pBuf, countof(buffer) - (pBuf - buffer), "%2d ", val[i]);
	}

	pBuf += sprintf_s(pBuf, countof(buffer) - (pBuf - buffer), "  %5d", sumX);
	COUT(buffer);
}

bool CheckFolkmanConditions(int v)
{
	/*
	THEOREM 5. (Jon Folkman (1967)) Let v be a positive integer. There are no admissible graphs
		with v points, if v satisfies one of the following conditions :
	(3.1) v is odd;
	(3, 2) v = 2p or 2p^2, where p is prime;
	(3.3) v < 30 and 4 does not divide v;
	(3.4) v < 20.
	*/
	if (v % 2 || v < 20 || v < 30 && v % 4)
		return false;

	v /= 2;
	if (isPrime(v))
		return false;

	const int res = static_cast<int>(sqrt(static_cast<float>(v)) + 0.000000001);
	return (res * res != v) || !isPrime(res);
}

int InconsistentGraphs(designParam *pParam, const char *pSummaryFileName, bool firstPath)
{
	if (!CheckFolkmanConditions(pParam->v))
		return 0;

	const int nVertex = pParam->v /= 2;
//	return 0;
	const auto kMax = nVertex - 2;
	const int len = kMax + 1;
	int *mult = new int[len * 6];
	int *step = mult + len;
	int *val = step + len;
	int *bufferTmp = val + len;
	int *pIndices = bufferTmp + 2 * len;
	lll = len;
	char buffer[256];
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
		sprintf_s(buffer, sizeof(buffer), "K = %d", k);
		COUT(buffer);

		CInterStruct *iStruct = NULL;
		bool printMult = true;
		int paramIdx = 0;
		while (true) {
			if (sum0 >= num && sum1 >= num * mult[i]) {
				sum0 -= num;
				sum1 -= num * mult[i];
				val[i] = num;

				if (!sum0 && !sum1 && IntersectionArrayIsValid(nVertex, k, val, mult, iMax, bufferTmp)) {
					int n = val[iMax];
					if (mult[iMax] == k && n)
						++n;
					else
						n = 1;

					auto iStructTmp = new CInterStruct(n);
					if (iStruct)
						iStruct->setNext(iStructTmp);
					else
						pParam->InterStruct()->setNext(iStructTmp);

					iStruct = iStructTmp;

					// Define possible values for 
					//   (a) mutual intersection of any two matrix rows;
					//   (b) number of elements, which belong exactly to i common blocks
					//   (c) number of elements, which belong to some block AND exactly to i common blocks 
					for (int j = 0; j <= iMax; ++j) {
						if (!val[j])
							continue;

						iStruct->lambdaPtr()->push_back(mult[j]);
						iStruct->lambdaAPtr()->push_back(val[j]);
						iStruct->lambdaBPtr()->push_back(val[j] * mult[j] / k);
					}

//					printSolution(printMult, iStruct, mult, iMax, ++paramIdx, val);
				}

				if (i && sum0) {
					num = sum0 - sum0 % step[--i];
					continue;
				}
			}
			else {
				if (num >= step[i]) {
					const int diff = num * mult[i] - sum1;
					if (diff > step[i] * mult[i])
						num -= diff / (step[i] * mult[i]) * step[i];
					else
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


		// To do the formated output, define maximal size of lambda sets
		int maxSize = 0;
		iStruct = pParam->InterStruct();
		CheckIntersections(iStruct, pIndices, k, nVertex);

		CInterStruct *iStructCurr = iStruct->getNext();
		if (!iStructCurr)
			continue;	// There are no valid solutions for current k

		printMult = true;
		paramIdx = 0;
		while (iStructCurr) {
			printSolution(printMult, iStructCurr, mult, iMax, ++paramIdx, val);
			if (maxSize < iStructCurr->lambda().size())
				maxSize = static_cast<int>(iStructCurr->lambda().size());

			iStructCurr = iStructCurr->getNext();
		}

		pParam->setLambdaSizeMax(maxSize);

		// Launching the enumeration for all constructed parameters
		pParam->r = k;
		iStructCurr = iStruct->getNext();
		while (iStructCurr) {
			if (iStructCurr->isValid()) {
				auto n = iStructCurr->mult();
				pParam->k = k / n;
				pParam->v = nVertex / n;

				// Copy current set of parameters into place where from they will be used
				const auto jMax = iStructCurr->lambda().size() - (n > 1 ? 1 : 0);
				iStruct->setMult(n);
				iStruct->lambdaPtr()->resize(jMax);
				iStruct->lambdaAPtr()->resize(jMax);
				iStruct->lambdaBPtr()->resize(jMax);
				for (int j = 0; j < jMax; ++j) {
					(*iStruct->lambdaPtr())[j] = iStructCurr->lambda()[j];
					(*iStruct->lambdaAPtr())[j] = iStructCurr->lambdaA()[j] / n;
					(*iStruct->lambdaBPtr())[j] = iStructCurr->lambdaB()[j] / n;
				}

				if (true/* k >= 4 *//*&& val[0] == 8 && val[2] == 4*/) {
					if (!RunOperation<MATRIX_ELEMENT_TYPE>(pParam, pSummaryFileName, firstPath))
						return 0;

					firstPath = false;
				}
			}

			// Remove just used intersection 
			CInterStruct *iStructNext = iStructCurr->getNext();
			delete iStructCurr;
			iStruct->setNext(iStructCurr = iStructNext);
		}
	}

	return 1;
}


/*

Main lists: 
  K = 3
  0  1  Sum(a(i)*i^2):
  1:  3  6       6
  K = 4
  0  1  2  4  Sum(a(i)*i^2):
  1:  4  0  4  1      32
  2:  0  8  0  1      24
  3:  3  0  6  0      24
  4:  1  4  4  0      20
  K = 5
  K = 6
  2  3  4  6  Sum(a(i)*i^2):
  1:  0  8  0  1     108
  2:  3  0  6  0     108
  3:  0  6  3  0     102
  K = 7
  =====================
  K = 3
  0  1  Sum(a(i)*i^2):
  1:  5  6       6
  K = 4
  0  1  2  4  Sum(a(i)*i^2):
  1:  6  0  4  1      32
  2:  2  8  0  1      24
  3:  5  0  6  0      24
  4:  3  4  4  0      20
  5:  1  8  2  0      16
  K = 5
  0  1  2  3  Sum(a(i)*i^2):
  1:  1  0 10  0      40
  K = 6
  0  1  2  3  4  6  Sum(a(i)*i^2):
  1:  3  0  0  6  0  2     126
  2:  0  0  9  0  0  2     108
  3:  2  0  0  8  0  1     108
  4:  0  0  6  4  0  1      96
  5:  1  0  3  4  3  0      96
  6:  0  0  6  2  3  0      90
  7:  1  0  0 10  0  0      90
  8:  0  0  3  8  0  0      84
  K = 7
  K = 8
  4  5  6  8  Sum(a(i)*i^2):
  1:  8  0  0  3     320
  2:  6  0  4  1     304
  3:  2  8  0  1     296
  K = 9
  6  7  9  Sum(a(i)*i^2):
  1:  9  0  2     486
  ==============================
  K = 3
  0  1  Sum(a(i)*i^2):
  1:  7  6       6
  K = 4
  0  1  2  4  Sum(a(i)*i^2):
  1:  8  0  4  1      32
  2:  4  8  0  1      24
  3:  7  0  6  0      24
  4:  5  4  4  0      20
  5:  3  8  2  0      16
  6:  1 12  0  0      12
  K = 5
  0  1  2  3  Sum(a(i)*i^2):
  1:  3  0 10  0      40
  K = 6
  0  1  2  3  4  6  Sum(a(i)*i^2):
  1:  0  6  0  6  0  1      96
  2:  2  0  6  4  0  1      96
  3:  0  0 12  0  0  1      84
  4:  4  0  3  0  6  0     108
  5:  1  6  0  0  6  0     102
  6:  4  0  0  6  3  0     102
  7:  3  0  3  4  3  0      96
  8:  0  6  0  4  3  0      90
  9:  2  0  6  2  3  0      90
  10:  1  0  9  0  3  0      84
  11:  1  0  6  6  0  0      78
  12:  0  0  9  4  0  0      72
  K = 7
  K = 8
  2  3  4  5  6  8  Sum(a(i)*i^2):
  1:  4  0  4  0  4  1     288
  2:  0  8  0  0  4  1     280
  3:  4  0  0  8  0  1     280
  4:  0  0 12  0  0  1     256
  K = 9
  K = 10
  K = 11

*/