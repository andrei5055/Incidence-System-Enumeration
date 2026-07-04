# Order‑2 (involution) automorphisms of P1Fs of K18 and K20 — status & plan

*Session record, 2026‑06‑28 (branch `perf/parallel-small-L`). Companion to
`MD/order2_parity_theorem.md`, `MD/order4_parity_theorem.md`, `MD/k18_engine_status.md`.*

## 0. Context

For the representative method, classifying all P1Fs with `|Aut|>1` reduces (Cauchy) to a
sweep over automorphisms of each prime order. Orders ≥ 3 are done for K18 and K20. **Order 2
is the gap.** An involution automorphism of a P1F fixes ≤ 2 vertices (proved lemma, the crux
of the order‑2 parity theorem), and the fixed‑point count is even, so on `N=2n` vertices only
two cycle types can occur:

* `2^{n-1} 1^2` — exactly **2 fixed points** ("2‑fixed" case)
* `2^{n}`       — **fixed‑point‑free** (fpf)

(K18: `2^8 1^2`, `2^9`. K20: `2^9 1^2`, `2^10`.)

## 1. The one‑fixed‑factor lemma (proved + validated)

**Lemma.** A P1F of `K_{2n}` with an involution automorphism `α` fixing exactly two
vertices `u,v` has **exactly one** `α`‑fixed factor, and it is forced to be the **diagonal
factor** `D = {u,v} ∪ {the n−1 diagonals}` (the `n` α‑fixed edges, which form one perfect
matching).

*Proof.* If two distinct factors `F,G` are both fixed, then `H = F∪G` is a Hamiltonian cycle
(perfectness) with `α(H)=H`. `H` spans all vertices, so `α|_H` is a cycle‑involution fixing
exactly `u,v` ⇒ the **vertex‑reflection** through `u,v` (they are antipodal on `H`). On an even
cycle the unique 2‑edge‑colouring gives the two edges at each fixed vertex opposite colours,
and a vertex‑reflection swaps them ⇒ it **swaps** `F` and `G` — contradicting both fixed. So
≤ 1 fixed factor. Parity: `α` permutes the `2n−1` (odd) factors, non‑fixed in pairs ⇒ the
number of fixed factors is odd ⇒ ≥ 1 ⇒ exactly 1. Each α‑fixed edge lies in the unique fixed
factor, and there are `n` of them = a full matching ⇒ that factor is `D`. ∎

**Validation** (`scratchpad/fix2.cpp`, three independent modes — enumerate‑all‑and‑filter,
enumerate‑α‑invariant‑directly, pin‑D): all agree and #fixed‑factors is always 1:

| K6 | K8 | K10 | K12 |
|----|----|-----|-----|
| 2  | 16 | 96  | 1536 |

Contrast: fpf involutions have *no* fixed vertices, the reflection argument doesn't bite, and
they keep **many** fixed factors (for K18 `2^9`: 9 fixed factors forming a P1F of K10). That is
why the fpf case is the hard one.

## 2. Existence is settled for the 2‑fixed case — Galois construction

The 2‑fixed case is **not** an emptiness question (no mod‑4 obstruction: K16 `2^7 1^2` found by
search, K18 `2^8 1^2` = AGL(1,17); `16≡0`, `18≡2 mod 4`, both non‑empty).

**Galois P1F `GK_{p+1}`** (`p` prime): vertices `GF(p) ∪ {∞}`; factor
`F_k = {∞,k} ∪ { {k+i, k−i} : i=1..(p−1)/2 }`. The involution `x ↦ −x` fixes `0` and `∞`,
swaps `x ↔ −x` ⇒ type `2^{(p−1)/2} 1^2`, and `Aut ⊇ AGL(1,p)` (order `p(p−1)`).

* `p=17` → **K18 `2^8 1^2`**, AGL(1,17), order 272  (matches the project's `|Aut|=272` result)
* `p=19` → **K20 `2^9 1^2`**, AGL(1,19), order 342  ← settles the K20 gap's existence
* `p=13` → K14 `2^6 1^2` (sanity)

**Verified** (`scratchpad/gk.cpp`): each is a genuine P1F (all factor‑pairs Hamiltonian) and
`x↦−x` is an automorphism of the stated type. **So K20 `2^9 1^2 ≠ ∅` — by construction, no
search.**

## 3. What the compute attempts established (negative results)

* **Scratch finder** (`scratchpad/fix4.cpp`, randomized restart, no symmetry‑breaking,
  ~16M nodes/s): K14 `<1M` nodes, K16 `~14M`, **K18 `2^8 1^2` not found in 4B nodes** (it *is*
  non‑empty — solver limitation). K20 out of reach.
* **Idea (a): pin D in the rep engine** (prototyped in `k18a2rep.cpp` behind env `REP_PIND`):
  builds, structurally correct (`rawTasks 1`), but **does not make it tractable** — frontier 211k+
  tasks and growing, no plateau, `classes=0`. ROOT CAUSE: pinning `D` does **not** shrink the
  centralizer — `D` is fixed by *all* of `C(α)` (`|C|≈2e7` for K18 `2^8 1^2`, `3.7e8` for K20
  `2^9 1^2`), so per‑node `setwiseStab` dedup still faces the whole group. **Dead end.** The pin
  code was **reverted** afterward — `k14a2rep`/`k16A2rep`/`k18a2rep` are normalized‑identical
  again and `REP_PIND` no longer exists in the tree. `k20a2rep.cpp` lacks the CGT engine entirely
  (pre‑existing pre‑CGT version); not worth porting given this.

**Conclusion:** enumeration of the 2‑fixed classes is intractable; existence is settled by §2.

## 4. Complete order‑2 landscape

| case            | type        | status        | by |
|-----------------|-------------|---------------|----|
| K18 `2^8 1^2`   | 2 fixed     | **non‑empty** | `GK_18` = AGL(1,17) |
| K18 `2^9`       | fpf         | **OPEN**      | — |
| K20 `2^9 1^2`   | 2 fixed     | **non‑empty** | `GK_20` = AGL(1,19) |
| K20 `2^{10}`    | fpf         | **empty**     | order‑2 parity theorem |

**Lone undecided order‑2 case anywhere: K18 `2^9` (fpf).** It is boxed in: no construction
reaches it (`GK` only ever yields the 2‑fixed type), the parity theorem does not apply
(`18≡2 mod 4`; K14 fpf being non‑empty shows there is no parity obstruction to exploit), and
brute force is ~`10^16` nodes.

Effect on the classification claims:
* **K20 "all aut>1": OPEN** — `2^9 1^2` is non‑empty so cannot be dropped, and its enumeration
  is intractable. Have: all aut≥3 + `2^{10}=∅` + `2^9 1^2` non‑empty (≥1 class, count unknown).
* **K18 "all aut>2": have all order ≥ 3.** Order‑2‑only classes (`2^8 1^2` non‑empty, `2^9`
  open) are the residue.

## 5. Plan — the one remaining shot at K18 `2^9`: the K10‑quotient antipodality argument

For a `2^9`‑invariant P1F of K18 the 9 fixed factors form a **P1F of K10** (quotient on the 9
pairs), plus 4 swapped pairs. Because `n=9` is **odd**, antipodal rotation *swaps* the two
1‑factors of a fixed‑pair 18‑cycle, so it is excluded; the only colour‑preserving fpf
cycle‑involution left is an **edge‑reflection**, which fixes 2 edges = 2 diagonals that must be
**antipodal** on that 18‑cycle. Each of the `C(9,2)=36` fixed‑factor pairs imposes one such
antipodality constraint, and every constraint is a statement purely about the (finite,
enumerable) set of **P1Fs of K10**.

**Goal:** translate the 36 antipodality constraints into a property of P1Fs of K10 and test it
against the complete list of K10 P1Fs. If none satisfies all 36, that is an **emptiness proof
for K18 `2^9` with no `10^16` search** (it runs over a handful of K10 objects). If some do, feed
the surviving fixed parts (+ the 4 swapped pairs) to a small targeted solver. Either way it
decides the lone open case off the back of K10's small classification, not the K18 tree.

Backstop if the argument is inconclusive: SAT/exact‑cover on the K10‑reduced instance.

### 5.1 K10 is UNIQUE (established 2026‑06‑28, `scratchpad/k10p1f.cpp`)

**K10 has exactly ONE perfect 1‑factorization up to isomorphism** (864 labeled with the first
factor fixed; 1 class. K8 also = 1, as known). So **the fixed part of every K18 `2^9` P1F is, up
to relabeling the 9 pairs, the SAME object** — the unique K10 P1F (= `GK_{10}`).

What this buys: the K18 `2^9` decision reduces to **"can the unique K10 P1F be lifted+completed
to a K18 `2^9` P1F?"** Edge/orbit budget: K18 has 9 diagonal‑orbits + 72 swapped‑orbits = 81
orbits. The 9 fixed factors are the lifted K10 factors — each = 1 diagonal + 4 FULL swapped
orbits (α‑invariance forces whole orbits), using 36 swapped‑orbits; the other 36 swapped‑orbits
form the 4 swapped factor‑pairs. The remaining freedom is: (i) for each of the 36 K10‑edges
`{x,y}`, WHICH of the two pair‑(x,y) orbits goes into the fixed factor ("signs"), and (ii) the 4
swapped pairs built from the complementary 36 orbits — all subject to perfectness. The within‑
pair flips (`2^9`) and `Aut(GK_10)` act on the 36 signs, cutting the search.

**Decisive experiment (next):** pin the unique K10 P1F (one canonical labeled lift) and
EXHAUSTIVELY search its completions (signs + 4 swapped pairs) with the flip/Aut symmetry
reduction. Exhausts empty ⇒ **K18 `2^9 = ∅` proven**; finds one ⇒ **non‑empty**. This is a
single pinned fixed part (NOT the `10^16` PIN1 search, NOT the random‑sample hunt) — its size is
the open question, to be measured.

### 5.2 Exhaustive fpf counter built + validated; naive version hits the wall (2026‑06‑28)

`scratchpad/fixfpf.cpp` — clean standalone exhaustive enumerator of α‑invariant P1Fs for the
fpf involution (exact cover by α‑orbit blocks + incremental Hamiltonian prune). **Validated
exact: K6=4, K8=0, K10=576** (matches the project's known values), and confirms **K12 `2^6`=0**
(consistent with the parity theorem, `12≡0 mod 4`). BUT with NO symmetry reduction the LABELED
count explodes via the centralizer (`2^n·n!`): **K14 `2^7` does not even FIND a first solution
in 60s** (it is non‑empty — the 16‑result), and K16/K18 are far beyond. So the naive completer
is correct but cannot reach K18. The K10‑uniqueness lever pays off ONLY with pin‑whole‑fixed‑
part + flip/Aut symmetry‑breaking — the substantial remaining build, payoff uncertain (qsearch
already estimated exhaustive K18 `2^9` ~`10^16` and pivoted to a no‑proof random hunt). Artifacts:
`scratchpad/fixfpf.cpp` (validated), `scratchpad/k10p1f.cpp` (K10 unique).

## Artifacts
* `scratchpad/fix2.cpp` — lemma validator (3 modes).
* `scratchpad/fix4.cpp` — randomized‑restart finder (allocation‑free).
* `scratchpad/gk.cpp`   — Galois `GK_{p+1}` constructor + P1F/automorphism verifier.
* (pin‑D experiment in `k18a2rep.cpp` behind env `REP_PIND` was reverted; siblings k14/k16/k18 normalized‑identical.)
