> # 🟢 RESOLVED 2026-06-29 — K18 2⁹ is NON-EMPTY
>
> The question this whole document chases is settled, and the answer is **non-empty** —
> the reformulation's *first* reading ("no obstruction ⇒ likely non-empty") was correct.
> A worked example was sitting in the harvested classification all along:
> `Logs_18-17-A2_Result/Complete_graphs/18/18x17x2_18/P0000000047.txt`, a P1F with
> **|Aut(M)| = 16** carrying two fixed-point-free involutions of type 2⁹ (e.g.
> α = (0 1)(2 4)(3 5)(6 8)(7 9)(10 12)(11 13)(14 16)(15 17)). Independently verified
> (`scratchpad/verify47.ps1`): 17 perfect matchings, 153 edges, all 136 factor-pairs are
> single 18-cycles (genuine *perfect* 1-factorization), and α permutes the 17 factors as
> **9 fixed + 4 swapped pairs** — exactly the forced structure derived below.
>
> **Paper census (complete for |Aut| ≥ 3, file `P_Results_0000001-0000531.txt`):** of the
> 531 K18 P1F classes {|Aut| 3:351, 4:144, 8:22, 16:12, 17:1, 272:1}, **exactly one**
> (the |Aut|=16 class P47) admits a 2⁹ involution. The 352 odd-order classes (|Aut| 3, 17)
> have no involution at all (Lagrange); |Aut|=272 = AGL(1,17) has 17 involutions, all of
> type 2⁸1². (The |Aut| = 2 stratum is not enumerated — intractable — so we count existence,
> not the number of 2⁹-admitting classes.) ⇒ the prove-∅ route, the ~10¹⁶-node hunt, and the
> qsearch/engine emptiness machinery below are now **moot for this cell.** The depth-15 stall
> was solver weakness, not an obstruction.
>
> *(Everything below is the 2026-06-28 record, retained for the structure theory it derived,
> which the example confirms.)*

---

# K18 2⁹ — quotient reformulation & strategic reassessment (2026-06-28)

*Goal of the session: prove K18 2⁹ = ∅ by argument. Outcome: the attempt produced a
clean reformulation of the problem over the K9 quotient that shows **emptiness is almost
certainly FALSE** — the structure is satisfiable at every level a parity/counting argument
can see, and the n=3,5,7 precedent (all non-empty) extends naturally. The high-value move
flips from "prove empty" to "construct an example" — and the reformulation makes that
search small.*

---

## The forced rigid structure (recap, from [order2_parity_theorem.md](order2_parity_theorem.md))

α = fpf involution of type 2⁹ on K18. Blocks B_i = {a_i,b_i}, i=1..9; diagonals
d_i = {a_i,b_i}. For **n=9 odd** the structure of any α-invariant P1F is FORCED:

- **9 fixed factors**, each containing **exactly one diagonal** (d_F=1, forced by
  (⋆): d_F+d_{F'}∈{0,2} and d_F≡n≡1 mod 2 ⇒ every d_F=1, and Σd_F=n ⇒ exactly n fixed).
- **4 swapped pairs** {G, αG} (s = (n−1)/2 = 4), each containing **no diagonal**.

## The K9 quotient (the new part)

Project K18 onto K9 by collapsing each block B_i → vertex i. An off-diagonal edge of K18
between blocks i,k comes in two **α-orbits**:
- **par** = {a_ia_k, b_ib_k}  (Z₂ voltage 0)
- **cross** = {a_ib_k, b_ia_k} (Z₂ voltage 1)

So K18's off-diagonal edges = a **Z₂-voltage double cover of K9** (36 edges × 2 orbits =
72 orbits) plus the 9 diagonals (the fibre edges). Every orbit lies in exactly one factor.

**Fixed factor i** = 1 diagonal (exposes block i) + 4 doubled edges = a perfect matching
of K9∖{i}, each matched edge carrying one orbit (par/cross). The 9 fixed factors:
- have a **forced bijection** factor i ↔ exposed vertex i (since each contains a distinct
  single diagonal);
- their 9 near-perfect-matchings are 8-regular on K9 with 4×9 = 36 edges = |E(K9)|, so
  **they partition E(K9)** — i.e. they form a **1-factorization of K10** on vertices
  {1,…,9,∞}, where matching i pairs ∞ with i (the exposed vertex). Each K9-edge is used by
  **exactly one** fixed factor, consuming **one** of its two orbits.

**Swapped pair** {G,αG}: G is a diagonal-free perfect matching of K18, so every block has
both vertices matched outward ⇒ G projects to a **2-factor C of K9** (2-regular). αG = the
voltage-flip of G, projects to the **same** C. The pair consumes, per edge of C, the
**other** orbit (one orbit per edge, 9 orbits per pair). Perfectness forces G∪αG to be a
single Hamiltonian 18-cycle ⇒ by the Z₂-lift rule **C is a single 9-cycle with an odd
number of cross edges**. The 4 swapped 2-factors partition the *leftover* orbit of every
K9-edge (36 edges) into **4 Hamiltonian cycles of K9** — i.e. a **Walecki-type
Hamiltonian decomposition** of K9 (exists: K_{2m+1} = m Ham cycles, m=4).

## Reformulation (exact)

A 2⁹-invariant P1F of K18 ⟺ a triple (F, Q, σ) with:
- **F** = a 1-factorization of K10 on {1,…,9,∞}  → the 9 fixed factors;
- **Q** = a decomposition of K9 into 4 Hamiltonian cycles  → the 4 swapped pairs;
- **σ : E(K9) → Z₂** = a sign per edge (par/cross) used by its fixed factor; the swapped
  pair on that edge takes the complementary orbit. Each swapped 9-cycle must have **odd**
  Σσ (so it lifts to one 18-cycle, not two 9-cycles);

**subject to PERFECTNESS**: all C(17,2)=136 pairwise factor-unions Hamiltonian. The
fixed–fixed unions and within-pair unions are already partly handled by the structure; the
binding constraints are the **fixed–swapped** and **swapped–swapped (cross-pair)** unions,
which depend on σ.

## Why emptiness is almost certainly FALSE

1. **No parity/counting obstruction remains.** Orbit counts, diagonal counts, (⋆), the
   ≤2-fixed lemma, and the per-pair Hamiltonicity are ALL satisfiable (F exists, Q exists
   via Walecki, σ has odd-weight 9-cycles available). Everything a parity argument can see
   is green. The ≤2-fixed lemma and order-2 parity theorem are already spent — they killed
   n even and the k≤7 mixed types; they say nothing against n=9.
2. **Precedent.** n=3 (K6, 4 P1Fs), n=5 (K10, 576), n=7 (K14, exist) are ALL non-empty,
   each with the identical reformulation shape (1-factorization of K_{n+1} + Walecki K_n +
   signs). n=9 is the next odd term; nothing distinguishes it structurally.
3. **The only contrary evidence is weak**: the |Aut|=272 (AGL(1,17)) K18 P1F happens to
   have no 2⁹ involution (its involutions are all 2⁸1²). That is one P1F out of many and
   says nothing about whether a *different* P1F carries a 2⁹ automorphism. Compute "not
   plateauing" only means the naive tree is huge — expected, since perfectness is global.

**Conclusion: do not invest in an emptiness proof — the truth is very likely "non-empty."**
The win is to **construct** one, which also definitively settles the cell.

## Recommended construction route (small search)

Search the triple (F, Q, σ) instead of the K18 tree — orders of magnitude smaller:
1. Pick Q = a fixed Walecki Hamiltonian decomposition of K9 (4 cycles). (Canonical choice.)
2. Pick F = a 1-factorization of K10 (start with the GK10 / standard one; 396 iso classes
   total, but try structured ones first).
3. Solve for σ ∈ Z₂^36 making all 136 unions Hamiltonian — this is a constraint problem
   over 36 bits with strong local structure (each Hamiltonicity = a Z₂-lift connectivity
   condition). Even brute 2³⁶ is large, but the per-cycle odd-weight + per-union conditions
   prune hard; better: set it up as a SAT/linear-algebra-over-GF(2)-plus-connectivity
   search, or seed σ from a known K18 P1F's coset structure.
4. Alternatively / as an oracle: take the user's known K18 P1Fs and test directly whether
   ANY admits a 2⁹ automorphism (relabel search) — if the AGL(1,17) one is the only P1F in
   hand, generate more (e.g. all P1Fs with |Aut| a multiple of 3,5,… already harvested by
   the rep engine) and screen their involution types for a 2⁹.

## Validation TODO (cheap, do first)

Re-derive the reformulation count on K6/K10/K14 and check it reproduces parity_check's
counts (K6=4, K10=576). If the (F,Q,σ) model reproduces those exactly, it is trustworthy
for K18 and we run the construction search. parity_check.cpp etc. are gone from scratchpad
(session-temp) — would need a small rebuild of just the (F,Q,σ) enumerator (much smaller
than the old engine).

## Enumerator built + validated (2026-06-28) — `tools/qsearch_quotient.cpp`

A clean, fast α-invariant P1F enumerator over the block structure (specialized candidate
generators: α-symmetric matchings for fixed factors, α-disjoint matchings for swapped
pairs; lowest-edge exact cover; pairwise-Hamiltonicity pruning; n-odd d=1 pruning).
Build: `tools/build_qsearch.bat`. Run: `qsearch_quotient <n> <COUNT|FIND> [nodecap]`,
env `PIN1=1` (root symmetry break), `SEED=k` (randomized restart), `PROG=N` (progress).

**Validated exactly** against the (lost) `parity_check.cpp`:
- K6 (n=3) COUNT = **4**, K8 (n=4) = **0**, K10 (n=5) = **576**  ✓ (labeled α-inv P1Fs)
- K14 (n=7) FIND = example found at node 16 (7 fixed + 3 swapped pairs) — n=7 NON-EMPTY ✓
- **PIN1 cross-check:** K10 pinned COUNT = 48, and 48 × 12 = 576 (12 = orbit of canonical
  first fixed factors) ⇒ the root pin is exactly WLOG. So a PIN1 exhaustive K18 COUNT = 0
  would *prove* 2⁹ = ∅; a PIN1 FIND is complete for existence.

## K18 2⁹ empirical result (2026-06-28): NO example yet; consistent depth-15/17 stall

- Deterministic + ~180 randomized restarts + 4 long PIN1 randomized FIND workers: **no 2⁹
  P1F found.** Node rate ~3k/s (per-node candidate-gen cost dominates).
- The search **repeatedly reaches depth 15 of 17 factors and never closes the last two.**
  This near-miss is stable across thousands of branches and all random seeds.
- **Reassessment.** This *tempers* the "almost certainly non-empty" reading above. Three
  independent engines (user's rep engine, the old engine.cpp, and now qsearch) all fail to
  produce a 2⁹ example, and the depth-15 wall looks like a genuine endgame obstruction —
  evidence now leans **"empty or very rare,"** consistent with the |Aut|=272 / AGL(1,17)
  datum (no 2⁹ involution). NOT yet conclusive: the exhaustive PIN1 COUNT has not
  terminated, and "didn't find" ≠ "doesn't exist."
- **Running (background, may take many hours):** (1) exhaustive PIN1 COUNT K18 — if it
  terminates at 0 → emptiness PROVED; if solCount>0 → an example exists. (2) 4 randomized
  PIN1 FIND workers — first to complete exhibits an example.

## NEW THEOREM (2026-06-28, proved + computationally verified): fixed factors = P1F of K_{n+1}

> **Theorem.** In any 2^n-invariant P1F of K_{2n} (n odd), the n fixed factors together with
> a new point ∞ (∞ ↔ the block each fixed factor exposes via its diagonal) form a **perfect
> 1-factorization of K_{n+1}.**

*Proof.* F_i ∪ F_k (two fixed factors) is α-invariant, so α acts on this 18-cycle as a
fixed-point-free element of D_{2n} fixing the two diagonals d_i,d_k — a reflection through two
antipodal edges. Any alternating cycle in the quotient M_i ∪ M_k (no diagonal) would lift to a
proper 2-regular subgraph of the single cycle — impossible. So M_i ∪ M_k is a Hamiltonian path
i→k; with ∞ this is a Hamiltonian cycle of K_{n+1}, i.e. the K_{n+1} 1-factorization is perfect. ∎

**Verified** (`VERIFYFIX`): all 4 K6 solutions → P1F of K4, all 576 K10 solutions → P1F of K6,
**0 violations**.

**Corollary — NO σ-parity obstruction.** Tracing the lift of F_i ∪ F_k: forward pass hits
(v_t, 1+S_t), return pass hits (v_t, S_t) — opposite sheets for ANY signs. So F_i ∪ F_k is a
single 18-cycle *independent of σ*. Hence the fixed factors constrain σ not at all, and the only
σ-parity constraint (swapped: Σ over each 9-cycle even) is freely satisfiable. **The clean parity
lever that killed n-even does not exist for n=9.** Any emptiness must come from the
non-symmetric cross-constraints (fixed∪swapped, swapped∪swapped) — no short proof in sight, and
n=3,5,7 satisfy them. **Updated assessment: the weight of structural evidence now leans
NON-EMPTY-but-rare**, and the decisive move is the accelerated search the theorem unlocks
(restrict fixed factors to a P1F of K10), not an emptiness proof.

## Accelerated search built (2026-06-28) — pin fixed = P1F of K10, search σ + swapped

`qsearch_quotient.cpp` env `ACCEL=1`: generates a random perfect 1-factorization of K_{n+1}
(the pinned fixed structure, per the theorem), restricts fixed-factor candidates to that
matching (signs σ still searched via Hamiltonicity), and covers the rest with swapped pairs.
- **Validated:** K10 and K14 examples found **instantly** via the accelerated route.
- **K18:** still **no example.** Two failure modes both observed: (a) random σ per pinned P1F
  almost never completes (σ is a 2³⁶ space, good σ is rare — must be *searched*, not sampled);
  (b) the proper σ-search per pinned K10 1-factorization stalls at **depth 14/17** at ~4k
  nodes/s and does not close. The matching-pin removes the ~105× fixed-matching branching but
  the σ+swapped tree remains large and the per-node cost dominates.

**Consolidated empirical verdict (5 independent methods).** User's rep engine, old engine.cpp,
full qsearch (depth 15/17), random (P1F,σ) sampling (≈40k, 0 hits), and pinned σ-search (depth
14/17) — **all fail to produce a K18 2⁹ P1F, each stalling 1–3 factors short.** This is real
(if soft) evidence for **empty or extremely rare**, consistent with |Aut|=272/AGL(1,17). It is
NOT a proof — none of these runs *exhausted*; they are time-bounded.

## Highest-value next steps (in priority)
0. **Settle it.** Two viable routes given the above: (i) **perf pass** on qsearch (bitset edge
   sets + incremental union-Hamiltonicity instead of re-walking) for a ~10–100× node-rate gain,
   then an exhaustive pinned run over all K10 P1F iso-classes = either an example or a real
   emptiness proof; (ii) **overnight randomized ACCEL sweep** (heavy-tailed: many K10 P1Fs ×
   seeds, shallow per-P1F cap) — best shot at exhibiting a rare example cheaply.
1. **Mine the depth-15 stall into a proof.** The 15→17 endgame obstruction is the candidate
   for an n=9-specific argument: characterize what the last 2 factors must be (a swapped
   pair closing the leftover orbits) and show the residual graph can't be a Hamiltonian
   swapped pair. This is the most likely route to settle 2⁹ = ∅ rigorously.
2. **Constrained reformulation search**: fix a canonical Walecki Q + a structured K10
   1-factorization F, search σ ∈ Z₂³⁶ with the per-cycle even-weight + union-Hamiltonicity
   constraints (SAT/GF(2)+connectivity). Much smaller than the K18 tree; could find a rare
   example the tree-search misses, or exhaust a slice.
3. **Screen existing K18 P1Fs**: if the user's rep engine harvested P1Fs with various |Aut|,
   test each for a 2⁹ involution directly (relabel search) — a positive instantly settles it.

## Pointers
- Tool: `tools/qsearch_quotient.cpp` (+ `tools/build_qsearch.bat`).
- Structure source: [order2_parity_theorem.md](order2_parity_theorem.md) (n-odd table).
- Prior verdict (compute exhausted): [TBD.md](TBD.md).
- Full engine log: [k18_engine_status.md](k18_engine_status.md).
