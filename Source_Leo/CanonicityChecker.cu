//
//  CanonicityChecker.cu
//  BIBD_Mac
//
//  Created by Andrei Ivanov on 2/10/18.
//  Copyright (c) 2018 Andrei Ivanov. All rights reserved.
//

#include "CanonicityChecker.cpp"

#if !CONSTR_ON_GPU
GroupOrderInfo *CGPU_CheckerInfo<TDATA_TYPES>::m_pOrderInfo;
size_t CGPU_CheckerInfo<TDATA_TYPES>::m_nCPUthreads;
COrderInfo **CGPU_CheckerInfo<TDATA_TYPES>::m_ppOrderInfo;
CTimerInfo CGPU_CheckerInfo<TDATA_TYPES>::m_timer;

__global__ void AssignCheckerGlobal(CMatrixCanonCheckerGPU<TDATA_TYPES> **ppCheckers, uint checkerIdx, CMatrixData<TDATA_TYPES> *pMatrix,
	SIZE_TYPE nRows, SIZE_TYPE nCols, MATRIX_ELEMENT_TYPE maxElem, bool IS_enum, const size_t *pColOrbInfo)
{
	auto pChecker = ppCheckers[checkerIdx];
	if (!pChecker) {
#if USE_OLD_CODE
		pChecker = ppCheckers[checkerIdx] = new CMatrixCanonChecker<TDATA_TYPES>(pMatrix, nRows, nCols, maxElem, IS_enum);
#else
		pMatrix->InitWithData(nRows, nCols, maxElem);
		pChecker = ppCheckers[checkerIdx] = new CMatrixCanonCheckerGPU <TDATA_TYPES>(pMatrix);
		pChecker->initiateColOrbits(nRows, 0, pMatrix->partsInfo(), IS_enum, NULL);
		if (!checkerIdx)
			pChecker->setEnumInfo(new CInsSysEnumInfo<TDATA_TYPES>());
#endif
	}

	// Restore ColOrbit information on GPU
	pChecker->restoreColOrbitInfo(nRows, pColOrbInfo);
}

template<>
bool AssignChecker<TDATA_TYPES>(CMatrixCanonCheckerGPU <TDATA_TYPES> **ppCheckers, CMatrixData<TDATA_TYPES> **pAllMatrData, uint checkerIdx,
	const CEnumerator<TDATA_TYPES> * pCanonChecker, cudaStream_t stream
#if TRACE_CUDA_FLAG
					, int myID
#endif
					)
{
	const auto *pMatrix = pCanonChecker->matrix();
	auto nRows = pMatrix->rowNumb();

	// Converting ColOrbit information to GPU
	const size_t nColOrb = pCanonChecker->copyColOrbitInfo(nRows);

	const auto pColOrbInfoBeg = pCanonChecker->CanonCheckerGPU()->ColOrbitData(t_CPU);
	auto pColOrbitDataCPU = pCanonChecker->CanonCheckerGPU()->ColOrbitData(t_GPU);
	CudaSafeCall(cudaMemcpyAsync(pColOrbitDataCPU, pColOrbInfoBeg,
		nColOrb * sizeof(pColOrbInfoBeg[0]), cudaMemcpyHostToDevice, stream));

	// Copying the matrix
	const cudaError_t err = cudaMemcpyAsync((char *)pAllMatrData[checkerIdx] + sizeof(CMatrixData<TDATA_TYPES>),
								pMatrix->GetDataPntr(), pMatrix->lenData(), cudaMemcpyHostToDevice, stream);
	TRACE_CUDA("  assignChecker 3: pAllMatrData[checkerIdx = %d] = %p  err = %d%s\n", checkerIdx, pAllMatrData[checkerIdx], err, err != cudaSuccess ? " - ERROR" : "")
	CudaSafeCallRet(err, false);

	AssignCheckerGlobal <<<1, 1, 0, stream>>> (ppCheckers, checkerIdx, pAllMatrData[checkerIdx],
									nRows, pMatrix->colNumb(), pMatrix->maxElement(), 
									pCanonChecker->IS_enumerator(), pColOrbitDataCPU);
	TRACE_CUDA("  assignChecker 4: err = %d%s  pMatrix = %p\n", cudaGetLastError(), cudaGetLastError() != cudaSuccess ? " (ERROR)" : "", pMatrix)
	CudaCheckError();
	return true;
}

template<typename T, typename S>
__global__ void MakeCopyGroupInfoGlobal(MatrixCanonCheckerGPU **ppCheckers, GroupOrderInfo *orderInfo, int CPU_threadIdx, COrderInfo *pOrderInfo) {
	auto *pEnumInfo = ppCheckers[0]->enumInfo();
	const int iMax = pEnumInfo->GetSize();
	orderInfo += CPU_threadIdx;
	if (!pOrderInfo) {
		if (orderInfo->nOrdersMax < iMax) {
			// Already allocated memory is not enough
			orderInfo->nOrders = -iMax;
			return;
		}

		orderInfo->nOrders = iMax;
		pOrderInfo = orderInfo->pOrderInfo;
	}
	else {
		orderInfo->pOrderInfo = pOrderInfo;
		orderInfo->nOrders = orderInfo->nOrdersMax = iMax;
	}

	for (int i = 0; i < iMax; i++)
		memcpy(pOrderInfo + i, pEnumInfo->GetAt(i), sizeof(COrderInfo));
}

template<>
void MakeCopyGroupInfo<TDATA_TYPES>(CMatrixCanonCheckerGPU<TDATA_TYPES>**ppCheckers, GroupOrderInfo *orderInfo, int CPU_threadIdx, cudaStream_t stream, COrderInfo *pOrderInfo) {
	MakeCopyGroupInfoGlobal <<<1, 1, 0, stream>>> (ppCheckers, orderInfo, CPU_threadIdx, pOrderInfo);
	CudaCheckError();
}

template<typename T, typename S>
__global__ void ResetEnumInfoGlobal(MatrixCanonCheckerGPU **ppCheckers) {
	ppCheckers[0]->enumInfo()->resetEnumInfo();
}

template<>
void ResetEnumInfoGlobal<TDATA_TYPES>(CMatrixCanonCheckerGPU<TDATA_TYPES>**ppCheckers, cudaStream_t stream) {
	ResetEnumInfoGlobal <<<1, 1, 0, stream>>> (ppCheckers);
	CudaCheckError();
}

//template<typename T, typename S>
__global__ void ReleaseCheckersGlobal(CMatrixCanonCheckerGPU<TDATA_TYPES> **ppCheckers) {
	delete ppCheckers[blockIdx.x];
}

void ReleaseCheckers(CMatrixCanonCheckerGPU<TDATA_TYPES> **ppCheckers, uint numCheckers, cudaStream_t stream)
{
	ReleaseCheckersGlobal <<<numCheckers, 1, 0, stream>>> (ppCheckers);
	CudaCheckError();
}

template<class T, class S>
__global__ void TestCanonicity(MatrixCanonCheckerGPU **ppCheckers, uchar *pCanonFlag, uint *pGroupInfo, uint nObj, int *pMatrFlag, bool noReplicatedBlocks)
{
	auto objIdx = threadIdx.x;
	auto pChecker = ppCheckers[objIdx];
	const auto *pMatrix = pChecker->matrix();
	pCanonFlag[objIdx] = pChecker->TestCanonicity(pMatrix->rowNumb(), pChecker, t_saveRowPermutations);
	if (pCanonFlag[objIdx]) {
		if (pChecker->groupIsTransitive())
			pCanonFlag[objIdx] |= 0x02;

		pGroupInfo[objIdx] = pChecker->groupOrder();
	}
	
	__syncthreads();
	if (objIdx == 0) {
		auto *pEnumInfo = pChecker->enumInfo();
		for (size_t i = 0; i < nObj; i++) {
			if (!pCanonFlag[i])
				continue;

			auto pChecker = ppCheckers[i];
			if (noReplicatedBlocks && pEnumInfo->constructedAllNoReplBlockMatrix())
				pEnumInfo->setNoReplBlockFlag(false);
			else
				pEnumInfo->updateConstrCounters(pMatrFlag[i], pChecker->groupOrder(), pChecker->groupIsTransitive());
		}
	}
}

void TestCanonicity(uint nMatr, CMatrixCanonCheckerGPU<TDATA_TYPES> **ppCheckers, uchar *pCanonFlag, uint *pGroupInfo, int *pMatrFlag, cudaStream_t stream)
{
	TestCanonicity <<<1, nMatr, 0, stream>>> (ppCheckers, pCanonFlag, pGroupInfo, nMatr, pMatrFlag, false);
	CudaCheckError();
}
#endif
