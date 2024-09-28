#pragma once
#include "k-SysSupport.h"
/*#include "DataTypes.h"

template<typename T, typename S>
class C_InSysCanonizator;

template<typename T, typename S>
LIBRARY_API C_InSysCanonizator<T, S>* createCanonizator(int v, int lenGroup);
*/
LIBRARY_API void* createCanonizator(int v, int lenGroup);