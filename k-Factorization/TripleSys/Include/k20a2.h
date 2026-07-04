#ifndef K20A2_H
#define K20A2_H
// ============================================================================
// k20a2.h — minimal K20 host for the representative-method classifier.
//
// Unlike K16A2/K18A2 there is NO cyclic search engine here: K20 is supported
// ONLY through the unseeded orderly-generation representative method, which is
// fully self-contained in k20a2rep.cpp (its own N=20 constants, zero deps on
// this class beyond the few members below). This class is just the KSolver
// shell the dispatcher constructs: solve() -> runExhaustiveSearch() ->
// runRepresentativeMethod(order). N=20 (NM=19, NEDGES=190).
//
// The automorphism order is K20_REP_ORDER, overridable at runtime via env
// REP_ORDER (single) or REP_ORDERS (comma list). K20_REP_SYMCAP bounds
// centralizer enumeration (matches K16/K18).
// ============================================================================
#include <cstdint>
#include <new>
#include <malloc.h>
#include "kBase.h"

#define K20_USE_REP_METHOD 1          // always on: the rep method is the only K20 path
#define K20_REP_ORDER      4          // default automorphism order to classify
#define K20_REP_SYMCAP     2000000LL  // enumerate C(alpha) only when |C(alpha)| <= this

struct Mask20_C { uint64_t m[4] = { 0, 0, 0, 0 }; };  // 256-bit (>=190 edges); unused by rep

class K20A2 : public KBase<Mask20_C> {
public:
    // N-derived constants (single knob: NP), used by the dispatcher's FactorParams.
    static constexpr int NP     = 20;          // players / vertices
    static constexpr int NM     = NP - 1;      // 19 factors in a 1-factorization
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

    K20A2(const FactorParams& factParam, int fixed3RowsIndex, int kThreads,
          const unsigned char* first3Rows, ResultCallback callback,
          void* cbClassPtr = NULL, bool bPrint = true);

    void solve(int mode = 0) override;
    bool addRow(int iRow, const unsigned char* source) override;

private:
    void runExhaustiveSearch();
    void runRepresentativeMethod(int order, int target = 0);   // defined in k20a2rep.cpp; target>0 = harvest mode

    int kThreads = 1;
    ResultCallback resultCallback = NULL;
    void* cbClass = NULL;
};

#endif // K20A2_H
