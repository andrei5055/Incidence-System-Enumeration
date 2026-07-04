#ifndef K14A2_H
#define K14A2_H
// ============================================================================
// k14a2.h — minimal K14 host for the representative-method classifier.
//
// K14 is a VALIDATION ORACLE for the K18/K20 method (its full P1F count is known
// from this project's own solver: 23 P1Fs, |Aut| {1:2,2:3,3:5,4:1,6:5,12:5,84:1,
// 156:1}), and a structurally better one than K16 because N=14 == 2 (mod 4) like
// K18 -- so it has the fixed-point-free even types (4^3.2^1, 2^7) K16 lacks.
//
// Like K20 there is NO cyclic engine here: K14 is supported ONLY through the
// unseeded representative method, fully self-contained in k14a2rep.cpp (own N=14
// constants). This class is just the KSolver shell the dispatcher constructs:
// solve() -> runExhaustiveSearch() -> runRepresentativeMethod(order). N=14
// (NM=13, NEDGES=91). Order = K14_REP_ORDER, overridable via env REP_ORDER /
// REP_ORDERS (with optional :target harvest). K14_REP_SYMCAP bounds centralizer
// enumeration (matches K16/K18/K20).
// ============================================================================
#include <cstdint>
#include <new>
#include <malloc.h>
#include "kBase.h"

#define K14_USE_REP_METHOD 1          // always on: the rep method is the only K14 path
#define K14_REP_ORDER      4          // default automorphism order to classify
#define K14_REP_SYMCAP     2000000LL  // enumerate C(alpha) only when |C(alpha)| <= this

struct Mask14_C { uint64_t m[2] = { 0, 0 }; };  // 128-bit (>=91 edges); unused by rep

class K14A2 : public KBase<Mask14_C> {
public:
    // N-derived constants (single knob: NP), used by the dispatcher's FactorParams.
    static constexpr int NP     = 14;          // players / vertices
    static constexpr int NM     = NP - 1;      // 13 factors in a 1-factorization
    static constexpr int NFIXED = 3;           // fixed starter rows (unused by rep)
    static constexpr int M_MAX  = 131072;      // pool cap (unused by rep; sizes nothing it touches)

    // Guard custom aligned operator new/delete from MemoryLeak.h's Debug `#define new DBG_NEW`.
#pragma push_macro("new")
#undef new
    void* operator new(size_t size) {
        void* p = _aligned_malloc(size, 64);
        if (!p) throw std::bad_alloc();
        return p;
    }
    void operator delete(void* p) {
        _aligned_free(p);
    }
#if defined(_DEBUG) && defined(_MSC_VER)
    void* operator new(size_t size, int blockType, const char* filename, int linenumber) {
        void* p = _aligned_malloc_dbg(size, 64, filename, linenumber);
        if (!p) throw std::bad_alloc();
        return p;
    }
    void operator delete(void* p, int blockType, const char* filename, int linenumber) {
        _aligned_free(p);
    }
#endif
#pragma pop_macro("new")

    K14A2(const FactorParams& factParam, int fixed3RowsIndex, int kThreads,
          const unsigned char* first3Rows, ResultCallback callback,
          void* cbClassPtr = NULL, bool bPrint = true);

    void solve(int mode = 0) override;
    bool addRow(int iRow, const unsigned char* source) override;

private:
    void runExhaustiveSearch();
    void runRepresentativeMethod(int order, int target = 0);   // defined in k14a2rep.cpp; target>0 = harvest mode

    int kThreads = 1;
    ResultCallback resultCallback = NULL;
    void* cbClass = NULL;
};

#endif // K14A2_H
