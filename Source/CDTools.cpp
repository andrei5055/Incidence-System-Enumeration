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

	auto numClasses = dayNumb;
	auto dayAdj = dayNumbMax - dayNumb;
	if (dayAdj) {
		numClasses = 0;
		pInSys->adjustData(dayAdj *= v / k);
		pCanonizer->adjustData(dayNumb, pInSys->colNumb());
	}

	pInSys->convertToBinaryMatrix(pMatrix, k, dayNumb);
	const auto retVal = pCanonizer->CanonizeMatrix(k, NULL, numClasses);
	if (dayAdj) {
		// Reverting previously adjusted parameters
		pInSys->adjustData(-dayAdj);
		pCanonizer->adjustData(dayNumbMax, pInSys->colNumb());
	}

	return retVal;
}