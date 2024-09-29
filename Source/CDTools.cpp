#include "CDTools.h"
//#include "DataTypes.h"
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

ctchar *runCanonizer(void* pCanon, ctchar* pMatrix, int k) {
	auto* pCanonizer = static_cast<C_InSysCanonizator<tchar, tchar>*>(pCanon);
	auto pInSys = pCanonizer->matrix();
	const auto numDays = (pInSys->rowNumb() - 1) / (k - 1);
	pInSys->convertToBinaryMatrix(pMatrix, numDays, k);
	return pCanonizer->CanonizeMatrix(-1);
}