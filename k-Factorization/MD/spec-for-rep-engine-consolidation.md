# Spec — consolidate the K14/K16/K18/K20 representative engines into one source

*Drafted 2026-06-28. Target attempt: next week (~week of 2026-07-06). Status: PROPOSAL,
not started. No code changed yet.*

## 1. Motivation

The four classifier sources

* `TripleSys/Source/k14a2rep.cpp`
* `TripleSys/Source/k16A2rep.cpp`   *(note the capital `A` in the filename — naming drift)*
* `TripleSys/Source/k18a2rep.cpp`
* `TripleSys/Source/k20a2rep.cpp`

are meant to be the same N-generic engine differing only in `N`. Today
`k14a2rep ≡ k16A2rep ≡ k18a2rep` (verified normalized-identical), but the design invites
**accidental divergence**, which has bitten us repeatedly:

* a parallelism use-after-scope bug had to be hand-mirrored across files (commit `c74f00f`);
* the all-covers `decomposeMissingEdges` fix had to be mirrored K16+K18;
* this session a `REP_PIND` experiment landed in `k18` only (since reverted).

And `k20a2rep.cpp` has silently fallen behind: it is the **pre-CGT** version (no
`cgtEngine.h`; no `coverGen`/`splitNode`/`schreierDedup`/work-queue parallelism). So the
"four copies" model already failed to keep them equal.

**Goal:** one source of truth for the engine body, so divergence is *structurally
impossible* and `k20` is brought up to parity for free.

## 2. Design decision — TEMPLATE, not a runtime class

A single class with `N` as a **runtime** value (sizing arrays at start) is the wrong
mechanism:

* `N` is currently `constexpr` (`constexpr int N = 18;`). The hot path
  (`is_perfect`, `applyAlpha`, `genM`, `cover`, `coverGen`, `canonKey`) runs billions of
  times per run on fixed-size stack arrays (`std::array<uint8_t,N>`, `int color[N][N]`,
  `path_end[NM][N]`). Runtime `N` forces heap allocation (≈100× slower — measured on the
  scratch solver) or `N_MAX=20` fixed arrays with de-optimized loops. Unacceptable on the
  long runs.
* It also isn't only numbers that differ: each class is `KxxA2 : public KBase<MaskNN_C>`,
  and **the mask type `MaskNN_C` is N-dependent**. One runtime class cannot cleanly be
  `KBase<one_mask>`.

**Chosen shape:** a `template<int N, class Mask>` shared engine body (header), instantiated
once per size by a thin per-N wrapper. This keeps `constexpr` performance AND single-source.
Precedent already in the tree: `cgtEngine.h` is `template<int N>`.

## 3. Proposed structure

```
TripleSys/Include/cgtEngine.h            (unchanged; already template<int N>)
TripleSys/Include/repEngine.h    (NEW)   template<int N, class Mask, /*traits*/> — THE engine:
                                          all of the current anonymous-namespace contents
                                          (Match/Perm, edge tables, buildCentralizer,
                                          buildGenerators, RepShared, RepWorker, cover,
                                          coverGen, splitNode, schreierDedup, enumTypes)
                                          plus runRepresentativeMethodImpl(order,target,...).

TripleSys/Source/k14a2rep.cpp    becomes a ~5-line stub:
TripleSys/Source/k16a2rep.cpp      #include "repEngine.h"
TripleSys/Source/k18a2rep.cpp      void K18A2::runRepresentativeMethod(int o,int t){
TripleSys/Source/k20a2rep.cpp        rep::run<18, Mask18_C>(o, t, kThreads, m_bPrint,
                                          K18_REP_SYMCAP, "[K18-REP]", resultCallback, ...);
                                      }
```

The wrapper keeps each class's required entry point and `KBase<MaskNN_C>` typing (integration
rules, §5); the engine body lives once in `repEngine.h`.

Rename `k16A2rep.cpp` → `k16a2rep.cpp` while here, to end the casing drift (update the
`.vcxproj`).

## 4. What differs per N — and how each becomes a parameter

| difference | today | after |
|---|---|---|
| vertex count | `constexpr int N` per file | template param `N` |
| base mask type | `KBase<MaskNN_C>` | template param `Mask` (or trait) |
| centralizer cap | `KNN_REP_SYMCAP` macro | arg / trait `symCap` (env `REP_SYMCAP` still overrides) |
| print tag | `"[KNN-REP]"` literal | arg `const char* tag` |
| default order | `KNN_REP_ORDER` | unchanged (read in the support file's dispatch) |
| result pipeline | `resultCallback`, `cbClass` | passed through from the wrapper |
| worker count / print gate | `kThreads`, `m_bPrint` | passed through from the wrapper |

Everything else (the 700-line body) is identical and moves verbatim.

## 5. Integration rules to honor (unchanged)

* The ONLY external entry point stays `KxxA2::runRepresentativeMethod`, called from
  `KxxA2::runExhaustiveSearch` (in `kXXa2support.cpp` / `k14a2.cpp` / `k20a2.cpp`). The
  `REP_ORDER[S]` env dispatch (with the `:target` counter — see those files) stays in the
  support layer and is NOT moved.
* Nothing new gets external linkage — the engine header's helpers stay in an unnamed
  namespace or are `static`/`inline`/templated.
* Worker count from `KxxA2::kThreads`; every print gated by `m_bPrint` (forwarded as
  `g_bprint`). Note: the current file-scope globals (`g_nodes`, `g_emits`, `g_stop`,
  `g_harvest`, progress state) must become **per-instantiation** (e.g. members of a
  per-`<N>` struct or `inline` template statics) so the four sizes don't share state.
  ⚠️ This is the subtle part of the refactor.

## 6. Validation plan — regression oracles (MUST reproduce byte-for-byte)

Run before and after; counts and full `|Aut|` histograms must match exactly.

* **K14** (`REP_ORDERS=2,3,7,13`): total `|Aut|>1` = 21 =
  `{2:3,3:5,4:1,6:5,12:5,84:1,156:1}`; per-order 2→16, 3→17, 7→1, 13→1.
* **K16** (`REP_ORDERS=3,5,7`): `|Aut|>2` = 30 = `{3:19,5:5,7:4,14:1,15:1}` (Gill–Wanless).
* **K18**: order-4 → 179, order-8 → 30, order-17 → 2.
* **K20**: order-19 → 7 = `{19:3,57:1,171:2,342:1}`.

Extra check unlocked by consolidation: after k20 inherits CGT, re-run any k20 order that
previously used the slow non-CGT path and confirm identical counts.

## 7. Risks & cautions

* **Concurrent code.** The worker pool + shared work-queue (`splitNode`, `QCAP`) is exactly
  where the prior refactor introduced a severe undercount (`7af0ce0`: use-after-scope on
  `queue`/`qmtx`/`active`). Moving globals to per-instantiation state is the same class of
  hazard. Review the lifetime of every shared variable vs the threads that read it.
* **Shared mutable globals → per-N state.** Easy to miss one; a missed one silently
  cross-contaminates two sizes' counts. The K14/K16 oracles catch it only if those sizes are
  actually exercised in the same process — so also test a multi-size run in one process.
* **Timing.** Do it as a deliberate branch, ideally AFTER the paper's numbers (§5 of
  `MD/paper_P1F_K18_K20.md`) are locked, so a refactor can't perturb reported counts.
* **Scope creep.** Do NOT fold in new features (e.g. the K10-pinning idea) during the
  refactor. Refactor first, prove parity, then build new things on the single source.

## 8. Step-by-step plan

1. Branch from current `perf/parallel-small-L`.
2. Capture golden outputs: run the §6 oracles on the current binaries, save logs.
3. Extract `k18a2rep.cpp`'s body into `repEngine.h` as `template<int N, class Mask>`;
   convert file-scope globals to per-instantiation state.
4. Reduce `k18a2rep.cpp` to the thin wrapper; build; reproduce K18 oracles exactly.
5. Convert `k16` and `k14` to wrappers; reproduce their oracles. Rename `k16A2rep.cpp`.
6. Convert `k20` to the wrapper (it now gains CGT). Reproduce K20 order-19 = 7; spot-check a
   previously-non-CGT order.
7. Run a single process exercising ≥2 sizes to catch shared-state bleaks.
8. Diff all four wrappers — they should be near-identical 5-liners.
9. Update `MD/paper_P1F_K18_K20.md` §6 reproducibility note (one engine, instantiated per N).

## 9. Acceptance criteria

* One engine header; four wrappers, each only differing by `N`, `Mask`, `symCap`, tag.
* All §6 oracle counts + histograms reproduced exactly.
* No measurable hot-path slowdown vs the `constexpr` baseline (template keeps it).
* `k20` now uses the CGT engine (parity with the others).
* Naming drift fixed (`k16a2rep.cpp`).

## References
* Current engines: `TripleSys/Source/k{14,16,18,20}*rep.cpp`; shared `cgtEngine.h`.
* Oracles & status: `MD/order2_landscape_K18_K20.md`, `MD/paper_P1F_K18_K20.md`.
* Prior refactor hazard: parallelism bug fixed in commit `c74f00f` (see memory
  `triplesys-decompose-speedup`).
