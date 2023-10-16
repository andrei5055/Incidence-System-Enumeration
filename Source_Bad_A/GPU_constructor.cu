#include "stdafx.h"
#include "InsSysEnumerator.h"
#include "EnumInfo.h"

#if CONSTR_ON_GPU

template<typename T, typename S>
__global__ void DummyFuncGPU(S v, S k, S lambda = 0) {
}

#include "InSysSolver.cpp"
#include "BIBD_Enumerator.cpp"
//#include "RowSolution.cpp"
#include "VariableMapping.cpp"
#include "InSysRowEquation.cpp"
#include "InSysEnumerator.cpp"


#endif