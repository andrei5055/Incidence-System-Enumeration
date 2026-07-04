# K18 2^9 engine — status & log (as of 2026-06-27)

*Session goal after the parity theorems: build a counter for the one remaining
fully-fixed-point-free case, K18 type 2^9 (n=9, "L=2" / order-2). This logs what was
built, what was tried, what was learned, and the open options. Branch
`perf/parallel-small-L`.*

---

## 1. Settled results (committed, solid — the real deliverables)

- **Order-2 parity theorem** ([order2_parity_theorem.md](order2_parity_theorem.md)): a
  fpf-involution (2^n) invariant P1F of K_{2n} exists only for n odd (2n≥6). ⇒
  **K16 2^8 = ∅, K20 2^10 = ∅** by proof. K4 lone exception.
- **Order-4k parity theorem** ([order4_parity_theorem.md](order4_parity_theorem.md)):
  a P1F admitting an automorphism of **order divisible by 4** exists only for n odd ⇒
  K16, K20 admit none of order 4,8,12,…. Crux lemma: an involution automorphism of a
  P1F fixes ≤2 vertices. Both fold into the paper §5.3.
- Net: of the fully-fpf over-cap cases, **K16 and K20 are done by argument**; **K18 2^9
  (n=9 odd) is the only live one** the theorem cannot decide.

## 2. The engine that works for K ≤ 14 (committed: `scratchpad/engine.cpp`)

Counts **C(α)-orbits** of α-invariant P1Fs (α = 2^n) by:
- forced-order exact cover over α-orbits of factors (diagonals appear as α-fixed edges);
- **during-search isomorph rejection**: dedup partials by their C(α)-canonical form
  (provably loses no complete-solution orbit);
- **IR canonical form** (`canon_ir.cpp` logic, folded into engine.cpp): individualization-
  refinement over the n blocks under C(α)=C2≀S_n, **no group enumeration** (scales past
  the brute group-min that's impossible at n=9). Refinement driven by the block-pair
  invariant pat[i][j] = multiset over factors of edge-count(0/1/2) between blocks i,j.
  Sheets resolved by an invariant-forced selection (replaces a 2^n brute);
- **root symmetry break**: force the single canonical first fixed factor F0 (all first
  fixed factors are C(α)-equivalent — WLOG);
- **128-bit-hash partial dedup** (replaced 300-byte strings ⇒ ~40× less memory,
  removed the OOM wall);
- **PIN filter** via env `PINS="x:y,…"` ("if neighbor of 0 is x then neighbor of 1 is y")
  + `DUMP=k` to print orbit-rep row-maps. ORDER_CAP tunable (argv[2]).

**Validated (matches brute group-min):** K6=1, K8=0, K10=2, K12=0. **K14 = 7
C(α)-orbits** (reaches exactly 7 leaves; tight dedup). Pin demo on K10: full row-map pin
→ 14 nodes (was 115), finds exactly the pinned orbit.

> NOTE: 7 is the **C(α)-orbit** count, not iso-classes (iso ≤ 7). The orbit→iso-class
> conversion (full-S_2n dedup of the few reps) is **still unbuilt**. The paper's
> "K14 order-2 = 16" counts ALL involution types (2^7 AND 2^6 1^2) up to iso — a
> different quantity, not directly comparable to our 7.

## 3. What was tried for K18 2^9 and why each stalled

| attempt | result | why it failed for K18 |
|---|---|---|
| generic dedup engine (string seen-set) | correct K≤14 | **OOM**: ~580 B/node, seen-set grows ~1:1 with nodes |
| brute within-cell canon | — | (n−1)!=40320 orderings at symmetric partials → 15min stuck on one node |
| IR canon | K14 205s→147s | fixes the canon explosion, but K18 node-count still huge |
| F0 root symmetry break | K10 230→115 nodes | only fixes level 1; levels 2–3 still symmetric |
| **rigid-structure engine** (`scratchpad/rigid.cpp`) | bounded mem (6MB) | Phase 1 = 9 fixed factors via precomputed **compat bitmasks**; Phase 2 = swapped pairs. **Slower than generic** (K14 >1300s) — diagonals-first ordering hurts dedup; abandoned for speed, kept as artifact |
| 128-bit hash dedup (in generic engine) | K14=7, bounded mem | solves OOM but not the tree size; K18 still ~600 nodes/s, ~3% pruning |
| pin-sharding for timing estimate | inconclusive | pins were placed on **late** factors (0~4,0~6,0~8); my lex-edge search order means those don't prune the tree-top, so ×15ⁿ extrapolation invalid. 4–5 pins still didn't finish a shard in 60–90s |

**Key K18 signals (generic engine, unpinned):** ~600–630 nodes/s, **only ~3% pruning**
(vs 24% at K14), 0 complete solutions in 131k nodes. The tree ≈ the full (huge) search.

## 4. The |Aut|=272 result the user provided (`scratchpad/check_result.cpp`)

A real K18 P1F, |Aut|=272 (the **AGL(1,17)/order-17** family). Tool reproduced |Aut|=272
(validates the Aut-finder) and proved **all 17 involutions are type 2^8 1^2 (2 fixed
points), ZERO fixed-point-free (2^9)**. So it is **not 2^9-invariant** — cannot validate
the L=2 engine, and the pin `4:8` from it doesn't correspond to anything my 2^9 engine
searches. Consistent with the AGL structure (x↦−x+b fixes b/2 and ∞) and our ≤2-fixed lemma.
Weakly suggestive that **K18 2^9 may be empty** (the most symmetric K18 P1F has no fpf
involution) — but one example proves nothing.

## 5. Honest assessment & open options

- **K18 2^9 is not tractable in the from-scratch engine** within a useful budget; all
  evidence points to a very large search. The pin technique that gives 10× in the
  **main program** (row-by-row build → pinning early rows prunes the subtree) does not
  transfer to this engine's lex-edge search structure.
- **Best vehicle for K18 2^9 = the main program + pin interface** (proven 10×, sharding
  works by design). Use pins to shard + extrapolate there, or to verify against any
  future known 2^9 result.

**Open next steps (pick up here):**
1. **Math (highest value):** try to prove **K18 2^9 = ∅** (or find an example) by argument,
   like K16/K20. The |Aut|=272 evidence + the rigid structure (9 fixed factors w/ 1
   diagonal each + 4 swapped pairs, all proven) are the handles.
2. Corrected pin-sharding in engine.cpp: pin the **earliest** factors (0~2, 0~3, …) to get
   a valid timing estimate from this engine (quick experiment).
3. Build the **orbit→iso-class conversion** (full-S_2n dedup of orbit reps) so K14 etc. give
   publishable iso-class counts.
4. Optional: a **2^8 1^2 engine variant** to verify against the |Aut|=272 result at K18
   size and open that mixed-involution case.

## 6. Artifacts (all in `scratchpad/`, committed)

- `engine.cpp` — the K≤14 counter (IR canon + during-search hash-dedup + F0 break + PINS/DUMP).
- `canon_ir.cpp` — standalone IR canonical-form validator (matches brute on K10).
- `rigid.cpp` — rigid-structure engine (compat bitmasks; bounded mem but slow; reference).
- `parity_check.cpp`, `order4_check.cpp` — theorem validators.
- `check_result.cpp` — Aut-group / involution-type analyzer for a given P1F.

Build (no nauty; pure MSVC): `vcvars64.bat && cl /O2 /EHsc /std:c++17 <file>.cpp`.
Run: `engine.exe N` (N=2n); env `PINS`, `DUMP`, argv[2]=ORDER_CAP.
