
#include "stdafx.h"
#include "EnumInfo.h"
//#include "BIBD_Enumerator.h"

#if CONSTR_ON_GPU

#include "Enumerator.cpp"

template<class T>
__global__ void EnumerateGPU(int v, int k, int lambda = 0) {
#if T_DESIGN_ENUM
	//        find_T_designParam(pParam->v, pParam->k, pParam->lambda);
	//        continue;        
	C_tDesign bibd(pParam->t, pParam->v, pParam->k, pParam->lambda);
#define DesignEnumerator    C_tDesignEnumerator<MATRIX_ELEMENT_TYPE>
#else
	C_BIBD<MATRIX_ELEMENT_TYPE> bibd(v, k, lambda);
#define DesignEnumerator    CBIBD_Enumerator<MATRIX_ELEMENT_TYPE>
#endif
	DesignEnumerator BIBD_Enum(&bibd, false, NO_REPLICATED_BLOCKS);
	char buff[256] = "AAAAAAAAAAAAAAAAA", buffer[256] = "TEST TEST";
	MAKE_JOB_TITLE(BIBD_Enum, buff, countof(buff));
//	sprintf(buff, "AAAAAAA");
	printf(buff);
	CInsSysEnumInfo<MATRIX_ELEMENT_TYPE> enumInfo(buff);
	const int mt_level = BIBD_Enum.define_MT_level(v);
	BIBD_Enum.Enumerate(mt_level, PRINT_TO_FILE, &enumInfo);
///	enumInfo.reportResult(buffer, countof(buffer));
//	outString(buffer, pSummaryFile);
	printf("\xd%s", buffer);
}

bool Enumerate(designRaram *pParam)
{
	EnumerateGPU<MATRIX_ELEMENT_TYPE> <<<1, 1>>> (pParam->v, pParam->k, pParam->lambda);
	CudaCheckError();
	cudaDeviceSynchronize();
	return true;
}

int cudaMemCmp(const unsigned char* left, const unsigned char* right, size_t length)
{
	int result = 1;
	while (result && (length-- > 0))
		result &= (left[length] == right[length]);

	return result;
}

#endif
