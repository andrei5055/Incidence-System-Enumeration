#include "stdafx.h"
#include "Enumerator.h"
#include "EnumInfo.h"
#include "InsSysEnumerator.h"

#if CONSTR_ON_GPU

template<class T>
__global__ void DummyFuncGPU(int v, int k, int lambda = 0) {
}

#include "InSysSolver.cpp"
#include "BIBD_Enumerator.cpp"
//#include "RowSolution.cpp"
#include "VariableMapping.cpp"
#include "InSysRowEquation.cpp"
#include "InSysEnumerator.cpp"


#endif