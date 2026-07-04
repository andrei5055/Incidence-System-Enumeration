# Spec: `K18A2Old` — starter-based solver for K18 L=4 (and other aut orders)

Status: DRAFT for review. No code written yet. Author: porting analysis 2026-06-24.

## 1. Motivation

The cyclic `K18A2` (pool-materialize + `decomposeMissingEdges`) is **stuck at L=4**: the
complement is dense (`num_colors = NM−L = 13`, `AvgPool ≈ 224k`), the exact-cover phase
dominates (~67%), and most time is spent proving non-existence for rejected automorphisms.
It is effectively unbounded at L=4 (and OOMs at L=2).

The **starter method** used by `K16A2Old` (3-row start matrices → candidate pool →
arc-consistency prune → clique/exact-cover search) does **all of K16 incl. L=2 in ~10 min**
because it never materializes a matching pool and the arc-consistency filter collapses the
candidate set hard.

### Measured time estimate for K18 L=4 (calibrated, not guessed)

| | starters | candidates/starter | work units | time | threads |
|---|---|---|---|---|---|
| K16A2Old (measured) | 1,647 | 60,000 | 9.88×10⁷ | 10 min | 10 |
| K18 L=4 (target) | ~20,000 | ~80,000 | 1.60×10⁹ | **~2.7 h baseline** | 10 |

- Throughput: 9.88×10⁷ / 600 s ≈ **165k candidate-rows/s** (end-to-end: pool build + arc-consistency + clique).
- Scale ratio 16.2× → 2.7 h at equal per-unit cost.
- K18 is deeper per unit (14 rows to place vs 12; 14 arc-consistency slots vs 12), apply **~1.2–2×** → **realistic ~3–5 h** wall-clock on the same 10-thread box.
- Tail risk to ~half a day only if order-4 clique trees branch much harder than K16's prime cases.
- 10 threads already saturate ~6 cores → no further parallel headroom on this machine.

**Conclusion: tractable in hours, vs. "never" for the cyclic method at L=4.**

## 1b. BENCHMARK FINDINGS (2026-06-25) — drive the architecture below

Built + ran an L=4 head-to-head (`g_benchL4` prototype, since reverted): baseline
`decomposeMissingEdges` vs. the same + an arc-consistency (unit-propagation) pre-pass.
Result over ~170 L=4 α's: **speedup 0.99–1.01×, AC early-rejects 0, mismatches 0.**

- **AC is sound but useless at L=4:** no complement edge is covered by a unique orbit, so
  unit propagation never fires; rejection-infeasibility is combinatorially deeper. → a cheap
  surgical patch to `decompose` will NOT yield the needed order-of-magnitude.
- Live phase split this run: `en/gr/cv ≈ 12% / 48% / 40%` (AvgPool ≈ 226k), **both group
  and cover dominate**, and the pool is rebuilt on *every* α.
- **Decisive insight:** the cyclic-decompose pool is **α-SPECIFIC** (complement of the
  α-built G) → cannot be amortized. K16A2Old's pool is **α-INDEPENDENT** (rows compatible
  with the 3 fixed rows) → built **once per starter**; α only regroups it into orbits.
  That once-per-starter (vs once-per-α) amortization is the real lever — confirming the
  starter engine, not decompose tweaks.

## 1c. HARD CONSTRAINT: order-4 α source (4 ∤ 18)

K16A2Old generates candidate automorphisms from the **dihedral of a Hamiltonian cycle**
(`get_transformations`). An N-cycle's rotations have order dividing N. **18 = 2·3², and
4 ∤ 18**, so an 18-cycle has NO order-4 rotation → K16A2Old's α-source produces **zero**
order-4 candidates for K18. Order-4 K18 automorphisms (types 4⁴·2¹, 4⁴·1², 4³·2³, …) must
come from elsewhere. The cyclic `K18A2` **already** enumerates exactly these via
`CycleBacktrackState` / `buildPermutation` + `generate_remaining_cycles` (proven; it's what
the 30-core run uses). → **reuse that generator; do NOT port `get_transformations`.**

## 2. REVISED APPROACH (2026-06-25): amortized completion inside K18A2

Because (a) the lever is the α-independent pool amortized across α's, and (b) order-4 α's
must come from K18A2's existing generator, the clean design is **NOT a separate clone** but
an **alternative completion path inside `K18A2`**:

1. **Keep** K18A2's proven order-4 α-generator (`CycleBacktrackState`) — unchanged.
2. **Build once per starter** (the 3 fixed rows) an α-INDEPENDENT candidate-factor pool =
   all perfect-matching rows pairwise-Hamiltonian with the 3 fixed rows (K16A2Old's
   `addRow`/`global_pool` logic, ported to NP=18, scalar), then **arc-consistency prune** it.
3. **Per α from the generator:** group the pruned pool into α-orbits, then **orbit clique**
   exact-cover (K16A2Old's `process_permutation`/`find_clique`, minus prime targeting) to
   emit every α-invariant P1F completion. This **replaces** the per-α `decomposeMissingEdges`
   (which rebuilds the 226k complement pool every call).
4. Gate behind a flag (e.g. `UseKSolve` bit or a member) so the existing decompose path and
   the L≥8=21 gate stay available for cross-validation.

This reuses the proven generator + canonizer, sidesteps 4∤18, and captures the amortization
win the benchmark identified. (A separate `K18A2Old` class is rejected: it would need its own
order-4 α-source, which 4∤18 makes awkward, and the cyclic method — not K16A2Old, which is
broken — is the real K18 oracle.)

### scalar/port notes still apply (from §3 below)
`Mask18_C` (256-bit) already covers 153 edges; reuse `find_cycles`/`is_perfect_scalar`
(NP=18-correct) and SKIP K16A2Old's `__m128i` 16-vertex kernels. The vestigial scaffolding
in `k18a2.h` (`PackedAdj`/`State`/`SearchContext`/`internal_solve`/`addRow`) can be the
skeleton to fill in (esp. `addRow` stub + undefined `internal_solve`).

## 2-OLD (superseded). Original approach: separate class, exact mirror

~~Create `K18A2Old` as a near-exact NP=18 mirror of `K16A2Old`~~ — superseded by §2 above.

- The vestigial `K16A2Old`-derived scaffolding currently sitting in `k18a2.h`
  (`PackedAdj`/`Factor`/`State`/`SearchContext`/`internal_solve`/`addRow`) is a leftover
  from when `K18A2` was first cloned from the old engine. It is **dead** (`addRow` is a
  stub, `internal_solve` is undefined). It can be filled in (see §2); remove it in a separate
  cleanup so `K18A2` is purely the cyclic method.

(Alternative — implement the old engine inside `K18A2` by filling in the stubs — is
rejected: it muddies the working class and breaks oracle symmetry.)

## 3. Porting deltas: `K16A2Old` → `K18A2Old`

### 3a. Constants
| K16 | value | K18 | value |
|---|---|---|---|
| `K16_N` | 16 | `K18_N` | 18 |
| `K16_MATCH` | 15 | `K18_MATCH` | 17 |
| `K16_FIXED` | 3 | `K18_FIXED` | 3 |
| `K16_SEARCH` | 12 (=MATCH−FIXED) | `K18_SEARCH` | 14 |
| `K16_M_MAX` | 60,000 | `K18_M_MAX` | **131,072** (≥80k + headroom) |
| edges | 120 (=16·15/2) | edges | **153** (=18·17/2) |
| `K16_WORDS` | derived | `K18_WORDS` | recompute from M_MAX |

The `s*2+7` slot labeling (K16A2Old.cpp:356) and any hard-coded `16`/`120`/depth bounds
must be re-derived from the K18 constants.

### 3b. Edge mask — already solved
`K16A2Old` uses `Mask256_C { uint64_t m[2] }` (128 bits) with `eid<64→m[0]`, `eid<128→m[1]`.
153 edges overflow 128 bits. **Reuse the existing `Mask18_C { uint64_t m[4] }` (256-bit)**
and add the third word everywhere a mask is built or AND-ed:
- build: `else if (eid < 192) m[2] |= (1ULL << (eid-128));`
- disjointness (`is_perfect_packed`, `fixedEdgesMask` test): add `|| (a.m[2] & b.m[2])`.
Pervasive but purely mechanical. `KBase<Mask18_C>` already makes base-class mask ops generic.

### 3c. SIMD width — the real hazard, resolved by scalar reuse
`K16A2Old`'s `apply_perm_16`, `check_cycle_simd`, `get_fast_sorted`, `pack_factor_adj` use
`__m128i` (16 bytes) to hold the 16-vertex adjacency/permutation. **18 vertices do not fit
`__m128i`**, and `_mm256_shuffle_epi8` cannot do a cross-128-bit-lane byte permute without
AVX-512 VBMI.

**Resolution: drop these SIMD kernels for K18 and reuse the already-correct scalar
equivalents** that exist for NP=18:
- `find_cycles(adj1, adj2)` — already implemented in `k18a2.cpp:65` for NP=18.
- `is_perfect_scalar(adj1, adj2)` — already used by the cyclic decompose.
- `apply_perm` (scalar) — already in `K18A2`.

The SIMD cycle/perm checks are a **minority** of the work; the hot path is the edge-mask
AND (stays fast via `Mask18_C`). The ~3–5 h estimate came from K16A2Old *with* SIMD, so
scalar K18 kernels land us toward the higher end of the bracket — acceptable. AVX2/AVX-512
re-vectorization is an **optional later optimization**, not a blocker.

### 3d. `get_fast_sorted` / `gfs_table`
The 16-byte pshufb lookup tables are K16-specific. Replace with a scalar sorted-factor
build for K18 (or skip if the sorted-factor canonicity path isn't needed — see §4).

## 4. RESOLVED (2026-06-24): full enumeration at run time, classify afterward

**Decision: Option (a).** Starters are **generic** (not order-4-seeded). The clique engine
enumerates **all** P1Fs reachable from the 20k starters (every automorphism order) at run
time; L=4 (and L=3, L=5, …) are extracted afterward by the existing `|Aut|` classification
(KSolve / CheckCanon pipeline). **No order-4 constraint inside the engine; no prime→composite
generalization of `get_transformations` needed** — the automorphism transforms are used only
for **isomorph rejection / canonicity** during search (order-agnostic), which is ported
faithfully from `K16A2Old` as-is. This removes the only real risk; the port is now
clone + reparameterize.

Open follow-up (non-blocking): confirm the 20k-starter set actually *covers* all of K18
(i.e. the starters are the canonical 3-row reps), and that the result/dedup volume is
tractable — verify empirically on the validation run (§6).

### (historical) original blocking questions — now resolved above

`K16A2Old` targets `target_primes = {7,5,3,2}` (k16a2Old.cpp:418) and uses
`get_transformations` / `get_transformations_general` for automorphism-based canonicity
reduction — the framework is **prime-order-oriented**. L=4 is order **4 (composite)**.

For K18, *unrestricted* enumeration is infeasible (the whole reason the cyclic method
exists), so `K18A2Old` **must restrict to P1Fs admitting an order-4 automorphism**. Three
ways this could be intended:

1. **Order-4 starters + α-invariant clique** (most likely): the ~20,000 starters are
   already order-4-seeded, and the clique search assembles α-invariant completions — i.e.
   `K18A2`'s order-4 framing, but driven by the K16A2Old clique engine instead of
   pool-decompose. Need: confirm the starter scheme and how α constrains candidate rows.
2. Add `4` to a `target_orders` list and generalize `get_transformations` to composite
   order 4 (verify it doesn't assume primality).
3. Enumerate broadly from starters and filter by `|Aut|` post-hoc (rejected for K18 —
   too many P1Fs).

**Questions for you:**
- Where do the ~20,000 starters come from — are they already order-4-constrained?
- Does each candidate row need to be α-orbit-compatible at generation time, or is α only
  enforced during the clique assembly?

## 5. Integration points

- **Dispatch** (`TripleSys.cpp:361`): mirror the K16 branch (lines 349–359). Add, under the
  existing `param(t_useKSolve) & (32+64)` + `m_numPlayers == 18` block:
  ```cpp
  if (param(t_useKSolve) & 64)
      m_pKSolver = new K18A2Old(factParam, ...);
  else
      m_pKSolver = new K18A2(factParam, ...);   // existing cyclic path
  ```
  i.e. `UseKSolve=64` → starter method, `UseKSolve=32` → cyclic (unchanged).
- **FactorParams**: `{ K18_N, K18_MATCH, K18_FIXED, K18_M_MAX, m_numDaysResult }`.
- **vcxproj / .filters**: add `k18a2Old.cpp` + `k18a2Old.h` (mirror the K16A2Old entries).
- **Config**: new `NewTests/param18a2Old.txt` with `UseKSolve=64`, `nPlayers=18`,
  `U1FCycles={{18}}`, pointing at the K18 3-row starter folder; new `run18a2Old.bat`.
- **ODR caution**: K16A2Old has no anonymous-namespace globals, but watch for any
  file-scope helpers shared with `k16a2Old.cpp` — wrap in an anonymous namespace if names
  collide (the `FlatMatchingMap` ODR crash precedent).

## 6. Validation plan

1. **Cross-check against the cyclic method (primary):** run `K18A2Old` and cyclic `K18A2`
   on the **same L where the cyclic method completes** (L≥8 → known **21** canonical, or a
   single mid-L like L=8) and confirm identical canonical sets. This validates the starter
   engine independently of K16A2Old's known issues.
2. **K16 oracle (secondary):** `K16A2Old` is the intended oracle, **but memory records it as
   currently broken** for p=3/5/7 (finds almost nothing; only p=2 works → 59). Either fix
   K16A2Old first, or rely on the cyclic-method cross-check (#1) plus the new `K16A2`
   results ({3:18, 5:5, 7:4, 14:1, 15:1}).
3. Only after #1 passes, trust `K18A2Old`'s L=4 output as the production result.

## 7. Effort & risk

- **3b (edge mask m[2])**: low effort, mechanical, low risk.
- **3c (scalar kernels)**: low effort (reuse existing), low risk; mild perf cost.
- **3a/3d/5 (constants, tables, wiring)**: low–medium effort, low risk.
- **§4 (order-4 targeting)**: **medium–high, the only real risk** — correctness of the
  automorphism restriction. Resolve the open questions first.
- **§6.1 validation**: must pass before trusting results.

Rough sequencing: (1) confirm §4 → (2) clone + reparameterize (§3a–3d) → (3) wire dispatch
+ build (§5) → (4) validate L≥8=21 (§6.1) → (5) production L=4 run (~3–5 h).
