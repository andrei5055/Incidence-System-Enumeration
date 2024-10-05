#include "CDTools.h"
#include "matrix.h"
#include "InsSysEnumerator.h"

void *createCanonizer(int v, int k) {
	const auto numDays = (v - 1) / (k - 1);
	const auto b = numDays * v / k;

	auto pInSys = new C_InSys<tchar, tchar> (v+1, b, 2, numDays - 1);
	uint enumFlags = t_enumDefault; // t_outRowOrbits + t_outRowPermute
	auto pInSysCanon = new C_InSysCanonizator<tchar, tchar>(pInSys, enumFlags);
	pInSys->prepareFirstMatrixRow(numDays);
	return pInSysCanon;
}

void releaseCanonizer(void* pCanon) {
	if (!pCanon)
		return;

	auto* pCanonizer = static_cast<C_InSysCanonizator<tchar, tchar>*>(pCanon);
	delete pCanonizer->matrix();
	delete pCanonizer;
}

ctchar *runCanonizer(void* pCanon, ctchar* pMatrix, int k, int dayNumb) {
	auto* pCanonizer = static_cast<C_InSysCanonizator<tchar, tchar>*>(pCanon);
	auto pInSys = ((C_InSys <tchar, tchar>*)(pCanonizer->matrix()));
	
	const auto v = pInSys->rowNumb() - 1;
	const auto b = pInSys->colNumb();
	const auto dayNumbMax = (v - 1) / (k - 1);
	if (dayNumb <= 0)
		dayNumb = dayNumbMax;

	auto dayAdj = dayNumbMax - dayNumb;
	if (dayAdj) {
		pInSys->adjustData(dayAdj *= v / k);
		pCanonizer->adjustData(dayNumb, pInSys->colNumb());
	}

	pInSys->convertToBinaryMatrix(pMatrix, k, dayNumb);
	const auto retVal = pCanonizer->CanonizeMatrix(k, NULL, dayNumb);
	static int ccc = 0;
#if 0
	FILE* f;
	if (!fopen_s(&f, "C:\\Users\\16507\\source\\repos\\Incidence-System-Enumeration\\aaa.txt", "a")) {
		pInSys->printOut(f, 16, ++ccc);
		fclose(f);
	}
#endif
	if (dayAdj) {
		pInSys->adjustData(-dayAdj);
		pCanonizer->adjustData(dayNumbMax, pInSys->colNumb());
	}

	return retVal;
}