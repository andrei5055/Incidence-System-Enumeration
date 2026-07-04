#include "k16A2.h"
#include <algorithm>
#include <cstring>
#include <cstdio>
#include <iostream>
#include <utility>
#include <mutex>
#include <intrin.h>
#include <set>

// File-scope N-constants for free functions in this TU (mirror K16A2's in-class values).
static constexpr int NP = 16;
static constexpr int NM = NP - 1;
static constexpr int NHALF = NP / 2;
static constexpr int NEDGES = NP * (NP - 1) / 2;

K16A2::K16A2(const FactorParams& factParam, int fixed3RowsIndex, int kThreads, const unsigned char* first3Rows,
    ResultCallback callback, void* cbClassPtr, bool bPrint) : KBase<Mask16_C>(factParam, bPrint) {
    for (int t = 0; t < 256; t++) {
        thread_root_idx[t] = 0x7fffffff;
        roots_done[t] = 0;
    }

    init(fixed3RowsIndex, kThreads, first3Rows, callback, cbClassPtr);
}

void K16A2::init(int fixed3RowsIndex, int kThreads, const unsigned char* first3Rows, ResultCallback callback, void* cbClassPtr) {
    if (kThreads > 256) kThreads = 256;
    this->thread_buffers.clear();
    for (int i = 0; i < kThreads; i++) this->thread_buffers.push_back(std::make_unique<ThreadLocalBuffers>());
    if (m_bPrint) {
        std::cout << "Init: K16A2 solver (Path Progress version), Compiled: " << __DATE__ << " " << __TIME__ << ", MS Compiler: " << _MSC_FULL_VER << ", kThreads: " << kThreads << std::endl;
    }
    this->fixed3RowsIndex = fixed3RowsIndex;
    this->kThreads = kThreads;
    this->resultCallback = callback;
    this->cbClass = cbClassPtr;
    start_time = std::chrono::steady_clock::now();
    last_print_time = start_time;
    bTimeSet = false;
    call_counter = 0;

    for (int s = 0; s < m_nFixedRows; s++) {
        memcpy(fixedRows[s].src, first3Rows + s * NP, NP);
        for (int i = 0; i < NP; i += 2) {
            const auto a = fixedRows[s].src[i]; const auto b = fixedRows[s].src[i + 1];
            fixedRows[s].adj[a] = b; fixedRows[s].adj[b] = a;
        }
    }
    for (int t = 0; t < 256; t++) {
        thread_root_idx[t] = 0x7fffffff;
        roots_done[t] = 0;
    }
    KBase::init();
}

bool K16A2::addRow(int rowNum, const unsigned char* source) {
    printf("not supported\n");
    exit(1);
    return false;
}



K16A2::CycleUnion K16A2::find_cycles(const uint8_t* adj1, const uint8_t* adj2) {
    CycleUnion cu; memset(&cu, 0, sizeof(cu));
    bool visited[NP] = { false };
    for (int i = 0; i < NP; ++i) {
        if (visited[i]) continue;
        int curr = i;
        int cycle_pos = 0;
        do {
            visited[curr] = true;
            cu.cycles[cu.count][cycle_pos++] = curr;
            curr = adj1[curr];
            visited[curr] = true;
            cu.cycles[cu.count][cycle_pos++] = curr;
            curr = adj2[curr];
        } while (curr != i);
        cu.lens[cu.count] = cycle_pos;
        cu.count++;
    }
    return cu;
}

void K16A2::solve(int mode) {
    runExhaustiveSearch();
}
