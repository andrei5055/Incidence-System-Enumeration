// ============================================================================
// k20a2.cpp — minimal K20 host: construct, then dispatch solve() to the
// representative-method classifier (runRepresentativeMethod, in k20a2rep.cpp).
// No cyclic engine; the rep method is unseeded so the fixed starter rows are
// ignored. Mirrors the REP_ORDER[S] dispatch used by K16A2/K18A2.
// ============================================================================
#include "k20a2.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <atomic>
#include <vector>

K20A2::K20A2(const FactorParams& factParam, int fixed3RowsIndex, int kThreads,
             const unsigned char* first3Rows, ResultCallback callback, void* cbClassPtr, bool bPrint)
    : KBase<Mask20_C>(factParam, bPrint) {
    this->fixed3RowsIndex = fixed3RowsIndex;
    this->kThreads = (kThreads > 256) ? 256 : kThreads;
    this->resultCallback = callback;
    this->cbClass = cbClassPtr;
    (void)first3Rows;   // rep method is unseeded; fixed rows are not used
    if (m_bPrint) {
        std::cout << "Init: K20A2 solver (representative method only), Compiled: "
                  << __DATE__ << " " << __TIME__ << ", MS Compiler: " << _MSC_FULL_VER
                  << ", kThreads: " << this->kThreads << std::endl;
    }
    KBase::init();
}

bool K20A2::addRow(int /*iRow*/, const unsigned char* /*source*/) {
    printf("K20A2::addRow not supported (rep method is unseeded)\n");
    return false;
}

void K20A2::solve(int /*mode*/) {
    runExhaustiveSearch();
}

void K20A2::runExhaustiveSearch() {
#if K20_USE_REP_METHOD
    // Run the classifier ONCE per process (the job may invoke the solver per start
    // matrix, but this search is unseeded so it must run a single time).
    static std::atomic<bool> repDone{ false };
    bool expected = false;
    if (repDone.compare_exchange_strong(expected, true)) {
        // Orders to classify. Priority: env REP_ORDERS (comma list) > REP_ORDER
        // (single) > K20_REP_ORDER default. All emit into the SAME result pipeline,
        // which canonicalizes the union -> one combined, deduplicated classification.
        // Per-order HARVEST TARGET via "order:target" (e.g. "3:7,19"): that order STOPS
        // after `target` distinct classes (skips the proof-of-exhaustion tail) -- complete
        // only if target == true count, else a lower bound. No ':' = run to completion.
        std::vector<std::pair<int,int>> orders;   // (order, target); target 0 = full run
        if (const char* envs = std::getenv("REP_ORDERS")) {       // full job: comma list
            const char* p = envs;
            while (*p) {
                int v = atoi(p), tgt = 0;
                const char* q = p; while (*q && *q != ',' && *q != ':') q++;   // scan this token
                if (*q == ':') tgt = atoi(q + 1);                              // optional :target
                if (v > 1) orders.push_back({ v, tgt });
                while (*p && *p != ',') p++; while (*p == ',') p++;
            }
        } else if (const char* env1 = std::getenv("REP_ORDER")) { // single-order run (also "order:target")
            int v = atoi(env1), tgt = 0; const char* q = env1; while (*q && *q != ':') q++; if (*q == ':') tgt = atoi(q + 1);
            orders.push_back({ v, tgt });
        } else {
            orders.push_back({ K20_REP_ORDER, 0 });
        }
        for (auto& ot : orders)
            runRepresentativeMethod(ot.first, ot.second);   // private member; uses kThreads + m_bPrint
    }
#endif
}
