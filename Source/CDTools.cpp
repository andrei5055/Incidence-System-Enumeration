#include "CDTools.h"
//#include "DataTypes.h"
#include "matrix.h"
#include "InsSysEnumerator.h"

/*template<typename T, typename S>
C_InSysCanonizator<T, S>* 
*/
void *createCanonizator(int v, int lenGroup) {
	auto pInSys = new C_InSys<tchar, tchar> (v+1, v, 2, 7);
	return new C_InSysCanonizator<tchar, tchar>(pInSys, 1);
}