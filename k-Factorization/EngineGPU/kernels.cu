#ifdef USE_CUDA
#include "cuda_runtime_api.h"
__device__ void memcpy_gpu(unsigned char* to, const unsigned char* from);

#include "TopGun.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include "TripleSys.cpp"
#include "endOfRowPrecalculation.cpp"
#include "precalculatedSolutions.cpp"
#include "cnvPrecalcRowsCompCheck.cpp"
#include "alldata.cpp"
#include "p1fCheck.cpp"
#include "p1fSupport.cpp"
#include "cnv.cpp"
#include "cnvSupport.cpp"
#include "cnvSupport2.cpp"
#include "cnvSupport3.cpp"
#include "cnvSupport4.cpp"
#include "kmSupport.cpp"
#include "checkLinks.cpp"
#include "checkLinksV.cpp"
#include "Support.cpp"
#include "initStartValues.cpp"
#include "processOneDay.cpp"
#include "getNextPlayer.cpp"
#include "checkPlayer1.cpp"
#include "checkPlayer3.cpp"
#include "checkCurrentResult.cpp"
#include "getUnselected.cpp"
#include "RowStorage.cpp"
#include "RowUsage.cpp"

#include "CheckMatrix.cpp"
#include "TrCycles.cpp"
#include "cycles.cpp"
#include "TrRepo.cpp"

#define GLOBAL_SYNCHRONIZE  1
#if GLOBAL_SYNCHRONIZE
#define CUDA_SYNCHRONIZE(cudaStatus)    cudaStatus = cudaDeviceSynchronize();   \
                                        gpuErrchk(cudaStatus, Error)
#define KERNEL_SYNCHRONIZE()
#else
#define CUDA_SYNCHRONIZE(cudaStatus)
#define KERNEL_SYNCHRONIZE()            __syncthreads();
#endif

#define PROFILING   1
#if PROFILING
typedef struct Timer {
    clock_t timeStamp;
    std::string comment;
    void SetTime(clock_t t, const char* pComment = "") {
        timeStamp = t;
        comment = pComment;
    }
} Timer;

#define TIME_STAMP_T(t, comment)    time[timeIdx++].SetTime(clock(), comment)
#else
#define TIME_STAMP_T(t, comment)
#endif
#define TIME_STAMP(comment)         TIME_STAMP_T(clock(), comment)

__global__ void initKernel(alldata** proc, const SizeParam p, const kSysParam* pSysParam, const int j)
{
    const int treadIndex = threadIdx.x + blockIdx.x * blockDim.x;
    proc[treadIndex] = new alldata(p, pSysParam, UseCheckLinksT, ImproveResults, CreateImprovedMatrix);
    KERNEL_SYNCHRONIZE();
}

__global__ void releaseKernel(alldata** proc)
{
    const int treadIndex = threadIdx.x + blockIdx.x * blockDim.x;
    delete proc[treadIndex];
    KERNEL_SYNCHRONIZE();
}

__global__ void enumerateKernel(alldata** pProc, sLongLong* pMatrNumb, int nMatrixToProcess, 
                                tchar* pMstart, int nRowsStart)
{
    const int treadIndex = threadIdx.x + blockIdx.x * blockDim.x;
    const auto numThreads = blockDim.x * gridDim.x;
    auto proc = pProc[treadIndex];
    const auto lenMatr = nRowsStart * proc->numPlayers();

#if 0
    int matrIdx = (nMatrixToProcess / numThreads) * numThreads + treadIndex - 1; // treadIndex;
    printf(" numThreads = %d  nMatrixToProcess = %d  matrIdx = %d\n", numThreads, nMatrixToProcess, matrIdx);
    if (treadIndex > 0)
        return;
    if (treadIndex < 10)
        printf(" matrIdx = %d\n", matrIdx);
#else
    int matrIdx = treadIndex;
    pMatrNumb[treadIndex] = 0;
#endif
    while (matrIdx < nMatrixToProcess) {
        auto* mstart = pMstart + matrIdx * lenMatr;
        auto* mfirst = pMstart;
        pMatrNumb[treadIndex] = proc->Run(treadIndex, eCalcResult, NULL, mstart, mfirst, nRowsStart);
        matrIdx += numThreads;
    }
    KERNEL_SYNCHRONIZE()
}

__device__ void memcpy_gpu(unsigned char* to, const unsigned char* from) {
    const int treadIndex = threadIdx.x + blockIdx.x * blockDim.x;
    to[treadIndex] = from[treadIndex];
}

#define gpuErrchk(ans, label) { if (!gpuAssert((ans), __FILE__, __LINE__)) goto label; }

inline bool gpuAssert(cudaError_t code, const char* file, int line)
{
    if (code == cudaSuccess)
        return true;
    
    fprintf(stderr, "GPU assert: %s %s %d\n", cudaGetErrorString(code), file, line);
    return false;
}

__declspec(dllexport) sLongLong runWithCuda(const TopGunGPU* pTG)
{
    const auto sTime = clock();
#if PROFILING
    Timer time[20];
    int timeIdx = 0;
    TIME_STAMP("");
#endif

    auto ssTime = sTime;
    auto seTime = sTime;
    const int numSets = 2816; // 2304;//  768;
    const int nn = 1;//3;

    const auto numPlayers = pTG->numPlayers();
    const auto groupSize = pTG->groupSize();
    const auto groupSizeFactorial = pTG->groupSizeFactorial();
    const auto gridSize = pTG->gridSize();
    const auto blockSize = pTG->blockSize();
    const auto matrToProceed = pTG->numMatrices2Process();
    tchar* inMatrices = NULL;
    const size_t lenMatr = pTG->numPlayers() * pTG->nRowsStart() * sizeof(tchar);
    const size_t inputLenTotal = matrToProceed * lenMatr;
    tchar* pntrTo;
    auto* pntrFrom = pTG->inputMatrices();
    kSysParam *paramPtr;
    alldata** dev_proc = NULL;
    sLongLong* matr_res = NULL;

    unsigned int nProc = gridSize * blockSize;
    unsigned int n_proc = nn * nProc;
    const auto len = n_proc * sizeof(alldata*);
    const auto lenParam = sizeof(*pTG->paramPtr());
    const auto lenResult = n_proc * sizeof(matr_res[0]);

    printf("Running program on %d GPU threads:  <<< %d, %d >>>\n", n_proc, gridSize, blockSize);

    cudaDeviceReset();

    sLongLong total = 0;
    sLongLong* res_cpu = new sLongLong[n_proc];

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaError_t cudaStatus = cudaSetDevice(0);
    gpuErrchk(cudaStatus, Error);

#define UNCREASE_GPU_STACK_SIZE 0
#if UNCREASE_GPU_STACK_SIZE
    {
        const size_t newStackSize = 4096;// 2048;
        size_t stackSize;
        cudaStatus = cudaDeviceGetLimit(&stackSize, cudaLimitStackSize);
        gpuErrchk(cudaStatus, Error)

        cudaStatus = cudaDeviceSetLimit(cudaLimitStackSize, newStackSize);
        gpuErrchk(cudaStatus, Error);
        printf(" Stack size was increased from %zd to %zd bytes.\n", stackSize, newStackSize);
    }
#endif
    // Allocate GPU buffers for n_proc processes
    cudaStatus = cudaMalloc(&dev_proc, len + lenParam);
    gpuErrchk(cudaStatus, Error);

    // Copying parameters to GPU memory
    
    paramPtr = (kSysParam *)((char *)(dev_proc) + len);
    cudaMemcpy(paramPtr, pTG->paramPtr(), lenParam, cudaMemcpyHostToDevice);
    gpuErrchk(cudaStatus, Error);

    TIME_STAMP("Preparing for initKernel");
    // Launch a kernel on the GPU.
    for (int j = 0; j < nn; j++) {
        initKernel << <gridSize, blockSize >> > (dev_proc + j * nProc, *pTG, paramPtr, j);
        TIME_STAMP("initKernel launch");
        // Check for any errors launching the kernel
        cudaStatus = cudaGetLastError();
        gpuErrchk(cudaStatus, Error);
        TIME_STAMP("initKernel check error");
    }

    CUDA_SYNCHRONIZE(cudaStatus)
    TIME_STAMP("initKernel synchronize");
        // Allocate memory for enumeration results
    cudaStatus = cudaMalloc(&matr_res, lenResult);
    gpuErrchk(cudaStatus, Error);

    cudaStatus = cudaMalloc(&inMatrices, inputLenTotal * numSets);
    gpuErrchk(cudaStatus, Error);

#if 1
    pntrTo = inMatrices;
    for (int i = 0; i < numSets; i++) {
        cudaMemcpy(pntrTo, pntrFrom, inputLenTotal, cudaMemcpyHostToDevice);
        gpuErrchk(cudaStatus, Error);
        pntrTo += inputLenTotal;
    }
#else
    
    for (int i = 0; i < matrToProceed; i++) {
        //printTable("Matr:", pntrFrom, 3, 12, 0);
        cudaMemcpy(pntrTo, pntrFrom, lenMatr, cudaMemcpyHostToDevice);
        gpuErrchk(cudaStatus, Error);
        pntrTo += lenMatr;
        pntrFrom += lenMatr;
    }
#endif
    ssTime = clock();
    TIME_STAMP_T(ssTime, "Preparing for enumerateKernel");
    for (int j = 0; j < nn; j++) {
        enumerateKernel << < gridSize, blockSize >> > (dev_proc + j * nProc, matr_res + j * nProc, matrToProceed * numSets,
            inMatrices, pTG->nRowsStart());

        // Check for any errors launching the kernel
        TIME_STAMP("enumerateKernel launch");
        cudaStatus = cudaGetLastError();
        gpuErrchk(cudaStatus, Error);
        TIME_STAMP("enumerateKernel check error");
    }

    CUDA_SYNCHRONIZE(cudaStatus);
    seTime = clock();
    TIME_STAMP_T(seTime, "enumerateKernel synchronize");

    for (int j = 0; j < nn; j++) {
        releaseKernel << < gridSize, blockSize >> > (dev_proc + j * gridSize * blockSize);
        TIME_STAMP("releaseKernel launch");
        // Check for any errors launching the kernel
        cudaStatus = cudaGetLastError();
        gpuErrchk(cudaStatus, Error);
        TIME_STAMP("releaseKernel check error");
    }
    CUDA_SYNCHRONIZE(cudaStatus);
    TIME_STAMP("releaseKernel synchronize");

    // Copy output vector from GPU buffer to host memory.
    cudaStatus = cudaMemcpy(res_cpu, matr_res, lenResult, cudaMemcpyDeviceToHost);
    gpuErrchk(cudaStatus, Error);
    TIME_STAMP("Starting Number of matrices calculation");
    for (auto i = n_proc; i--;)
        total += res_cpu[i];

    printf("Total %zd\n", total);
    TIME_STAMP("Number of matrices calculation");
Error:
    cudaFree(dev_proc);
    cudaFree(matr_res);
    cudaFree(inMatrices);
    delete[] res_cpu;

    TIME_STAMP("Release CUDA memory");
    printf("Total %zd matrices constructed, build time = %.2f - %.2f (sec)\n", total, (clock() - sTime)/1000., (seTime-ssTime)/1000.);
#if PROFILING
    for (int i = 1; i < timeIdx; i++) {
        printf("Time #%2d = %6.2f (s) spent in %s\n", i, (time[i].timeStamp - time[i - 1].timeStamp) / 1000., time[i].comment.c_str());
    }
#endif
    return cudaStatus == cudaSuccess ? total : -1;
}

__declspec(dllexport) int closeGPU() {
    // cudaDeviceReset must be called before exiting in order for profiling and
    // tracing tools such as Nsight and Visual Profiler to show complete traces.
    if (cudaDeviceReset() == cudaSuccess)
        return 0;

    fprintf(stderr, "cudaDeviceReset failed!");
    return 1;
}

CC bool alldata::improveMatrix(int improveResult, tchar* bResults, const int lenResult, tchar** pbRes1) {
    return true;
}
#endif