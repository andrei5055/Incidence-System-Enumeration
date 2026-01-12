#ifdef USE_CUDA

#include "TopGun.h"

int TopGunGPU::Run()
{
    __declspec(dllimport) sLongLong runWithCuda(const TopGunGPU * pntr);

    if (readMatrices() < 0)
        myExit(1);

    const auto total = runWithCuda(this);
    if (total < 0) {
        return -1;
    }

    return 0;
}

TopGunGPU::~TopGunGPU() {
    __declspec(dllimport) int closeGPU();
    closeGPU();
}
#endif