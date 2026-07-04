# Accelerating the order-2 (fixed-point-free involution) case — two approaches

*Design note. Targets the one open bottleneck of the representative method: order $2$,
cycle type $2^n$ on $K_{2n}$ (the fixed-point-free involution). Both ideas were
prompted by a "split / combine" and an "FFT-like change of basis" intuition; this note
makes each concrete, states what is solid vs. open, and gives a validation plan against
the **known** answer $K_{14}$ order-$2 = 16$.*

---

## 0. The problem and why it is dense

Fix the fixed-point-free involution $\alpha=2^n$ on $V=\{0,\dots,2n-1\}$; it pairs the
vertices into $n$ blocks $B_i=\{x_i,x_i'\}$ with $\alpha:x_i\leftrightarrow x_i'$. The
centralizer is the hyperoctahedral group
$$C(\alpha)=C_2\wr S_n,\qquad |C(\alpha)|=2^n\,n!,$$
the $S_n$ permuting blocks and the $C_2^n$ flipping inside blocks. For $K_{18}$,
$|C(\alpha)|=2^9\cdot 9!\approx1.9\times10^8$ — far past the enumeration cap, so the
orderly-generation dedup runs on a weak generator-only path and the canonical search
tree never plateaus. Crucially, **exact per-node stabilizers and full-group transporter
dedup were already tried and the tree stayed huge** — so the obstruction is the size of
the *canonical tree*, not weak isomorph rejection. That rules out "better dedup" and
points to changing the *object* (Approach A) or the *counting principle* (Approach B).

Throughout, recall the rep method already decomposes an $\alpha$-invariant
factorisation into **$\alpha$-orbits of factors**: an orbit is either a single
$\alpha$-*fixed* factor ($\alpha F=F$) or an $\alpha$-*swapped* pair $\{F,\alpha F\}$.

---

## 1. Edge orbits under $\alpha$ (shared by both approaches)

The $\binom{2n}{2}=n(2n-1)$ edges fall into three $\alpha$-orbit types:

| type | representative | size | count |
|---|---|---|---|
| **diagonal** | $\{x_i,x_i'\}$ | 1 (α-fixed) | $n$ |
| **parallel** | $\{x_i,x_j\},\{x_i',x_j'\}$ | 2 | $\binom n2$ |
| **crossed** | $\{x_i,x_j'\},\{x_i',x_j\}$ | 2 | $\binom n2$ |

Check: $n + 2\binom n2 + 2\binom n2 = n + 2n(n-1) = n(2n-1)$. ✓ The two off-diagonal
orbits for a block-pair $\{i,j\}$ are the **parallel** and **crossed** lifts of the
single quotient edge $ij\in K_n$.

---

## 2. Approach A — quotient by the involution (shrink the object)

**Idea.** Work on the $n$-vertex quotient $K_n = K_{2n}/\alpha$ (super-vertices =
blocks), so the search ranges over small labelled objects on $n$ points instead of
matchings of $2n$ points.

### 2.1 What an $\alpha$-fixed factor becomes

An $\alpha$-fixed perfect matching $M$ ($\alpha M=M$) is a union of edge-orbits.
Writing $d$ = #diagonal orbits and $e$ = #off-diagonal orbits used,
$2d+4e=2n\Rightarrow d+2e=n$, so:

> **An $\alpha$-fixed factor $\;\longleftrightarrow\;$ a partial matching of $K_n$ with
> $d\equiv n\pmod 2$ isolated vertices, together with a binary label
> (parallel / crossed) on each matched edge.**

The $d$ isolated quotient-vertices are the blocks closed by a diagonal edge; the $e$
matched quotient-edges are off-diagonal orbits, each carrying its parallel/crossed bit.
For $n$ odd ($K_{14}$, $n=7$) every fixed factor has $d\ge 1$ diagonal edges. This
object has $O(2^{n/2}\cdot(n{-}1)!!)$ possibilities — vastly fewer than the matchings of
$K_{2n}$.

### 2.2 What an $\alpha$-swapped pair becomes

For $\{F,\alpha F\}$, $F$ is an arbitrary perfect matching and $\alpha F$ its image; the
pair jointly covers a set of edge-orbits (each off-diagonal orbit split one edge to $F$,
one to $\alpha F$; no diagonal edge can appear — it would force $\alpha F=F$). So a
swapped pair $\longleftrightarrow$ **an arbitrary perfect matching of $K_n$ lifted by a
choice, per quotient-edge, of which of its two endpoints' parities $F$ takes** — i.e. an
*oriented/■-decorated* perfect matching of $K_n$. (This is the orbit the rep method
already enumerates; the quotient just renames it.)

### 2.3 The factorisation / perfectness condition

This is the real work and the genuinely open part. Two global constraints must be
pushed through the quotient:

1. **Cover.** Every diagonal orbit once, every parallel/crossed orbit once, across all
   factor-orbits. The $n$ diagonal edges are distributed among the $\alpha$-fixed
   factors only; their placement is highly constrained (each fixed factor takes
   $d\equiv n\bmod 2$ of them).
2. **Perfectness (the coupling).** $F_i\cup F_j$ must be a single Hamiltonian
   $2n$-cycle. This does **not** localise to the quotient — a Hamiltonian cycle of
   $K_{2n}$ projects to a closed walk on $K_n$ that uses each super-vertex twice, with
   the parallel/crossed labels deciding whether the lift is one $2n$-cycle or two
   $n$-cycles. Working out the label-parity rule that makes the lift a single cycle is
   the crux; it is a $\mathbb F_2$ (sign) condition on the labels along cycles of the
   quotient walk, closely related to the orientability/"holonomy" of the $C_2$-bundle
   $V\to V/\alpha$.

**Status.** Solid: the object-level reduction (§2.1–2.2) — the search domain genuinely
shrinks from $K_{2n}$ to labelled $K_n$. This is the same structure underlying *starters*
and the $GK_{2p}$ construction, and **why** $2n\equiv2\pmod4$ behaves specially ($n$ odd
forces $d\ge1$ diagonal edge per fixed factor).

### 2.4 Perfectness rule — DERIVED + VALIDATED (2026-06-27)

$V\to V/\alpha$ is a $C_2$-bundle. Encode each vertex as (block $i$, sheet $s\in\{0,1\}$),
$\alpha$ flipping $s$. Every edge has a **sheet-flip** $\sigma$: diagonal and crossed
flip ($\sigma{=}1$), parallel does not ($\sigma{=}0$). For two $\alpha$-fixed factors
$F,F'$ encoded per block as `DIAG` or `(partner j, flip b)`, walk the union from $(0,0)$,
alternating $F,F'$ steps — `DIAG`$\to(i,1{-}s)$, `(j,b)`$\to(j,\,s\oplus b)$ — for $n$
double-steps; **$F\cup F'$ is a single $2n$-cycle iff the walk returns to $(0,0)$ only at
the end.** This is the monodromy of the double cover, made into an $n$-vertex computation.

**Validated** (`scratchpad/quotient_proto.cpp`): on *every* pair of $\alpha$-fixed factors
for $K_8,K_{10},K_{12},K_{14}$ the quotient rule equals true $K_{2n}$ perfectness — **0
mismatches** (K14: 1{,}696{,}506 pairs, incl. all the diagonal-edge cases). So perfectness
is genuinely *local to $K_n$*. This is the foundational piece and the part with the
$n$-odd diagonal subtlety; it works.

### 2.5 What is already proven vs. the one remaining frontier

> **UPDATE (2026-06-27): the frontier canonical form is BUILT and validated.**
> `scratchpad/canon_ir.cpp` implements a fast $C(\alpha)=C_2\wr S_n$ canonical form by
> **individualization-refinement over the $n$ blocks** (no group enumeration — so it
> scales to $n=9$ where $|C(\alpha)|\approx1.9\times10^8$ is unbuildable). Refinement is
> driven by a genuine $C(\alpha)$-invariant: for each block-pair $(i,j)$, the multiset
> over factors of the edge-count ($0/1/2$) between blocks $i,j$. Sheet ($C_2$) flips are
> resolved by min-serialization. **Validated:** on $K_{10}$ it reproduces the brute
> group-min partition exactly ($2$ orbits, identical classes); invariance
> $\mathrm{canon}(gP)=\mathrm{canon}(P)$ holds on $400/400$ random trials at both $n=5$
> and $n=9$; avg $\sim1.15$ ms/call at $n=9$ (synthetic objects, no group needed).
> **Remaining to reach $K_{18}$:** (1) wire it into orderly generation with canonical
> augmentation so $|X|$ is never enumerated (the actual speed lever), validate the ladder
> $K_{10}\!\to\!2$, $K_{14}\,2^7\!\to\!16$ iso, $K_{16}\,2^8\!\to\!0$ (now a *proven*
> free check), then run $K_{18}\,2^9$; (2) if the per-node cost bites, optimize canon_ir
> (per-block sheet selection instead of the $2^n$ brute, automorphism pruning to cut
> within-cell orderings). The hard conceptual piece is done; what's left is the
> generation harness + the tractability question for the $2^9$ tree.

> **UPDATE 2 (2026-06-27): orderly-generation harness built (`scratchpad/engine.cpp`),
> profiled, and the next blocker pinned down.** The harness = forced-order exact cover
> over $\alpha$-orbits + **during-search isomorph rejection** (dedup partials by their
> $C(\alpha)$-canonical form — provably loses no complete-solution orbit). Validated:
> $K_{10}\!\to\!2$ and $K_{12}\!\to\!0$ (correct), with a **48× tree reduction**
> ($11018\!\to\!230$ nodes at $K_{10}$) from the partial dedup. Two canon optimizations
> landed: invariant-forced sheet selection (replaces the $2^n$ sheet brute — most blocks'
> orientation is forced by a flip-independent vertex invariant, only genuinely-symmetric
> blocks are brute'd) and an ordering-cap tunable.
> **But $K_{14}$ does not finish.** Profiling shows the bottleneck is the canonical
> form's **brute within-cell ordering enumeration** for *symmetric partials*: at shallow
> nodes (few factors placed) the unplaced blocks are genuinely symmetric, so refinement
> leaves big cells $\Rightarrow$ up to $(n{-}1)!$ orderings per canon call. Capping the
> orderings (skip dedup when too symmetric) just moves the cost: it removes pruning
> exactly at the shallow nodes where it matters most, and the tree re-explodes. **The
> real fix is automorphism-pruned IR** (nauty-quality canonical labelling: target-cell
> individualization + refinement + automorphism detection to collapse the $(n{-}1)!$
> orderings), which is the proper next build. Also still pending: the
> $C(\alpha)$-orbit $\to$ isomorphism-class conversion (full-$S_{2n}$ dedup of the orbit
> reps).

> **UPDATE 3 (2026-06-27): root symmetry break added; engine now completes $K_{14}$.**
> All first fixed factors are $C(\alpha)$-equivalent ($S_{n-1}$ + sheet flips act
> transitively), so the engine forces the single canonical $F_0$ (diagonal block $0$,
> blocks $1..n{-}1$ paired $(1,2)(3,4)\dots$ parallel) — WLOG, and it skips the $(n{-}1)!$
> canon calls otherwise spent collapsing them. **Engine validated against the brute
> group-min on $K_6{=}1,K_8{=}0,K_{10}{=}2,K_{12}{=}0$** (all match) and now runs through
> $K_{14}$: **$7$ $C(\alpha)$-orbits** ($521{,}598$ nodes, $205$ s; reaches exactly $7$
> leaves — tight dedup). This $7$ is the $C(\alpha)$-orbit count (the iso-class count is
> $\le7$ and needs the still-pending orbit$\to$class conversion; it is *not* the paper's
> "$K_{14}$ order-$2=16$", which counts all involution types $2^7$ and $2^6 1^2$ up to
> iso). The $F_0$ break helps but $K_{14}$ is still $205$ s — the symmetric-partial canon
> cost (levels $2$–$3$) remains the wall, so $K_{18}$ needs the faster canon (link nauty
> or hand-roll automorphism-pruned IR) regardless.

Pulling the two prototypes together, **every *correctness* piece of Approach A is now
validated** — for swapped factors and general-pair perfectness the quotient walk is just
the ordinary vertex walk (`is_perfect`), so nothing new was needed there:

| piece | validated by | result |
|---|---|---|
| perfectness from quotient (α-fixed pairs, diagonals) | `quotient_proto` | 0 mismatches, K8–K14 |
| full enumeration of α-invariant P1Fs (both factor types) | `burnside_proto` `countInv(⟨α⟩)` | K10 = 576, K8 = 0 (cross-checked vs direct) |
| $C(\alpha)$-orbit dedup of α-invariant P1Fs | `burnside_proto` direct-orbit count | K10 = 2 ( = Burnside) |

So the quotient representation, its perfectness, and the $C(\alpha)$-orbit reduction are all
*correct*. **The single remaining obstacle is purely performance:** the dedup as prototyped
is `enumerate-all` $\times$ brute-group ($O(|X|\cdot|C(\alpha)|)$), which is exactly why the
over-cap $2^9$ case is hard. Making it fast = an **efficient canonical form of the labelled
$K_n$ quotient object under $C_2\wr S_n$** (so the per-level orderly-generation dedup is one
canonical-form call, not an orbit walk). That is standard graph-canonicalisation /
partition-refinement (nauty / McKay) territory — a real implementation (or link an existing
canonicaliser), and it **converges with the "canonical augmentation" idea (§3 of the
options)**: both the quotient route and the McKay route need the same engine — a fast
canonical form of small objects under $C_2\wr S_n$. That engine is the frontier; everything
upstream of it is proven. Build it, validate the orbit/class counts against $K_{10}$ and
$K_{14}\,2^7$, then run $K_{16}\,2^8$ (smallest over-cap case the rep method cannot finish).

---

## 3. Approach B — Burnside / cycle-index over $C(\alpha)$ (change the counting)

**Idea.** Replace one huge isomorph-rejecting search by a group-average of cheap
fixed-structure counts — the rigorous form of "compute in a symmetry-adapted basis,"
the non-abelian analogue of using the Fourier basis for convolution.

### 3.1 The average

Let $X$ = set of $\alpha$-invariant P1Fs. By Cauchy–Frobenius,
$$\#\{C(\alpha)\text{-orbits on }X\}=\frac1{|C(\alpha)|}\sum_{g\in C(\alpha)}\mathrm{Fix}(g),
\qquad \mathrm{Fix}(g)=\#\{F\in X: gF=F\}=\#\{\langle\alpha,g\rangle\text{-invariant P1Fs}\}.$$

Two accelerations make this practical:

- **Sum over conjugacy classes, not elements.** $\mathrm{Fix}$ is a class function, so
  $\sum_g\mathrm{Fix}(g)=\sum_{[g]}|[g]|\,\mathrm{Fix}(g)$. $C_2\wr S_7$ has only a few
  hundred conjugacy classes (indexed by pairs of partitions / "signed cycle types"), so
  there are a few hundred terms, not $645{,}120$.
- **Most terms vanish or are tiny.** $\mathrm{Fix}(g)$ counts P1Fs invariant under the
  *larger* group $\langle\alpha,g\rangle$; for most $g$ that group is too rich to admit
  any P1F, so $\mathrm{Fix}(g)=0$. Each nonzero term is itself a small rep-method run
  (one representative of the bigger cycle type) — the expensive dense tree is replaced
  by many cheap constrained counts. This is exactly the cycle-index / Pólya structure.

### 3.2 The honest subtlety (orbits ≠ isomorphism classes)

The average gives $C(\alpha)$-*orbits* of $\alpha$-invariant P1Fs. That is **not**
directly the number of isomorphism classes admitting a $2^n$ automorphism: two
$\alpha$-invariant P1Fs can be isomorphic via a $g\notin C(\alpha)$ that sends $\alpha$
to a *different* $2^n$-involution of the target's automorphism group. A class $C$
contributes, to the orbit count, the number of $C(\alpha)$-classes of $2^n$-involutions
in $\mathrm{Aut}(C)$; it is $1$ exactly when $\alpha$ is the unique such involution up to
$\mathrm{Aut}(C)$-conjugacy, and $>1$ otherwise (a finite correction). The clean,
complete bookkeeping that yields the full $|\mathrm{Aut}|$-**distribution** (not just a
total) is the **table of marks** of $C(\alpha)$: for each conjugacy class of subgroups
$H\le C(\alpha)$ count P1Fs with stabiliser exactly $H$, then read off classes by
stabiliser type. So Approach B cleanly gives a *related* number fast; matching the
rep method's iso-class histogram needs the table-of-marks refinement.

### 3.3 Validation plan ($K_{14}$, where everything is checkable)

$K_{14}$, $\alpha=2^7$: $|C(\alpha)|=2^7\cdot7!=645{,}120$ — small enough to enumerate
$C(\alpha)$ *and* to run the rep method to completion (it already gives the proven
order-$2=16$). Concretely:

1. Enumerate the conjugacy classes of $C_2\wr S_7$ (signed cycle types).
2. For each class rep $g$, build $\langle\alpha,g\rangle$ and count its invariant P1Fs
   with the existing engine (set the representative cycle type of $\langle\alpha,g\rangle$).
3. Burnside-average $\Rightarrow$ $C(\alpha)$-orbit count; check it reproduces $16$
   after the §3.2 involution-multiplicity correction (most $K_{14}$ order-$2$ classes
   have $\alpha$ as their unique $2^7$-involution, so the correction is small and
   explicit).
4. If it matches and is fast, repeat the *count* (total only) for $K_{16}$ $2^8$ and
   $K_{18}$ $2^9$ — the cases the direct method cannot finish — as an independent
   result, then add the table-of-marks pass for the histogram.

**Status.** Solid: the average and the conjugacy-class/zero-term accelerations. Needs
care: the orbits-vs-classes correction (§3.2) and the table of marks for the
distribution. Best **near-term experiment** because it is self-contained and checkable
against a known number.

### 3.4 Prototype results (2026-06-27) — machinery VALIDATED

A standalone prototype (`scratchpad/burnside_proto.cpp`) computes, for a fixed-point-free
$\alpha=2^{N/2}$, both sides independently: `O_direct` (enumerate all $\alpha$-invariant
P1Fs, group by $C(\alpha)$) and `O_burnside` $=\frac1{|C(\alpha)|}\sum_{g}\mathrm{Fix}(g)$
with each $\mathrm{Fix}(g)=\#\langle\alpha,g\rangle$-invariant P1Fs computed by a *separate*
exact-cover (not read off the orbit set — so agreement is a real test).

| $K_N$ | total P1Fs (H=id) | $\alpha$-inv. (labeled) | `O_direct` | `O_burnside` | machinery |
|---|---|---|---|---|---|
| $K_8$ $(2^4)$ | 960 ($=40320/42$, confirms $\lvert\mathrm{Aut}\rvert{=}42$) | 0 | 0 | 0 | ✓ (no fpf involution) |
| $K_{10}$ $(2^5)$ | 90720 | 576 | **2** | $7680/3840=$ **2** | ✓ **nonzero, exact** |

**$K_{10}$ is the decisive check:** the independent per-element `Fix(g)` average reproduces
the direct orbit count (2) on a genuine nonzero instance — the Burnside *algorithm* is
correct. It also **confirms §3.2 empirically**: $K_{10}$ has a *single* P1F up to
isomorphism, yet `O_direct = 2`, so that P1F carries two $C(\alpha)$-inequivalent
$2^5$-involutions — the orbits-vs-classes correction is real and necessary to turn the
Burnside number into an isomorphism-class count.

The conjugacy-class form (group $C(\alpha)$ by signed cycle type — one $\mathrm{Fix}$ per
class $\times$ size) was then implemented and **re-validated on $K_{10}$**: 36 classes,
only 8 with $\mathrm{Fix}>0$, $O_\text{burnside(cc)}=2=O_\text{direct}$. So both the
algorithm and its scalable conjugacy-class form are correct.

### 3.5 Prototype verdict (2026-06-27) — Burnside is correct but **NOT a speed-up**

Pushing to $K_{14}$ exposed the fatal flaw. The Burnside sum's **identity term**
$\mathrm{Fix}(\mathrm{id})=|X|=\#\{\alpha\text{-invariant P1Fs}\}$ (the full *un-reduced*
set) is exactly the quantity orderly generation **avoids** by quotienting the search by
$C(\alpha)$. Measured: $K_{14}$ $2^7$'s $\mathrm{Fix}(\mathrm{id})$ did not finish in
$>85$ s of un-reduced enumeration, versus the rep method's **7.3 s** for the *whole*
$C(\alpha)$-reduced classification of the same type. Every Burnside/table-of-marks
formulation needs $\mathrm{Fix}$ at the trivial subgroup $=|X|$, so there is no escape:
the dominant term *is* the original bottleneck, and it is strictly **worse** than the rep
method (which never enumerates $|X|$). For $K_{18}$ $2^9$, where even the $C(\alpha)$-
*reduced* search does not terminate, the un-reduced $|X|$ is hopeless.

**Conclusion.** Approach B is a sound, validated *cross-check* (and the $K_{10}$ result —
1 isomorphism class but 2 $C(\alpha)$-orbits — is a clean illustration of the
orbits-vs-classes subtlety), but it does **not** accelerate the dense case. The real
lever is **Approach A (quotient by the involution)** — shrinking the *object* from
$K_{2n}$ to labelled $K_n$ — because that reduces $|X|$ itself, the one thing Burnside
cannot. Effort should go there; Burnside is retired as an acceleration route.

---

## 3.6 STRUCTURAL FINDING (2026-06-27) — fpf-involution P1Fs exist only for $2n\equiv2\pmod4$

> **NOW PROVED (2026-06-27).** The conjecture below is a **theorem** — see
> [order2_parity_theorem.md](order2_parity_theorem.md). For $n$ even ($2n\ge8$) there is
> **no** $2^n$-invariant P1F; $K_4$ is the lone exception. Consequently
> **$K_{16}\,2^8=\varnothing$ and $K_{20}\,2^{10}=\varnothing$ by argument, not compute**,
> leaving $K_{18}\,2^9$ as the only live fully-fpf case. The crux is a one-line lemma: the
> antipodal rotation of a $2n$-cycle swaps its two $1$-factors iff $n$ is odd, so for $n$
> even no swapped factor pair can exist and the all-fixed diagonal count cannot close.

Building the canonicaliser (`scratchpad/canon_proto.cpp`) surfaced a parity law that
**reshapes the whole order-2 question.** Counting *fixed-point-free*-involution-invariant
P1Fs (type $2^n$) directly:

| $K_N$ | $n$ | $2n \bmod 4$ | fpf-involution P1Fs |
|---|---|---|---|
| $K_8$  | 4 | 0 | **0** |
| $K_{10}$ | 5 | 2 | 576 (2 $C(\alpha)$-orbits) |
| $K_{12}$ | 6 | 0 | **0** |
| $K_{14}$ | 7 | 2 | yes (rep method, type $2^7$) |

**They occur only when $2n\equiv 2\pmod4$ ($n$ odd)** — the *same* parity as the paper's
order-4 observation (§5.3 of the manuscript), strongly suggesting a single underlying
obstruction. The consequence is large:

> **Theorem (PROVED 2026-06-27, [order2_parity_theorem.md](order2_parity_theorem.md)).**
> $K_{16}$ has no $2^8$-invariant P1F and $K_{20}$ has no $2^{10}$-invariant P1F. The two
> *intractable fully-fixed-point-free* over-cap types we could not search are **empty**,
> leaving $K_{18}\,2^9$ ($2n=18\equiv2$) as the only genuinely non-empty fully-fpf case.

This is the cleanest possible kind of progress: instead of a heroic search to (probably)
prove $0$, a **parity theorem** settles $K_{16}\,2^8$ and $K_{20}\,2^{10}$ outright. The
proof (the order-2 analogue of the order-4 parity, both at $2n\equiv2\pmod4$) runs through
the diagonal-edge / monodromy structure of §2.4 exactly as anticipated: every diagonal
sits in an $\alpha$-fixed factor (D1), the per-factor count $d\equiv n\pmod2$ (D2), and
the union-cycle automorphism — the antipodal rotation swaps the two $1$-factors of a
$2n$-cycle **iff $n$ is odd** — kills every swapped factor pair when $n$ is even, after
which the diagonal count cannot close. It does **not** fully close order-2 (the mixed
types like $K_{18}\,2^8 1^2$, with fixed points, are separate), but it removes the worst
fully-fpf cases by argument rather than compute.

## 4. Recommendation (updated after the §3.5 prototype verdict)

- **Approach B: validated, but retired as a speed-up.** The prototype proved the
  algorithm correct (K10) and is a useful independent cross-check, but its identity term
  re-introduces the exact bottleneck it was meant to avoid (§3.5). It cannot beat the rep
  method on the dense case. Keep it only as a small-case verifier.
- **Approach A is the path — foundation now VALIDATED (§2.4).** The quotient perfectness
  rule (incl. the $n$-odd diagonal edges) reproduces true perfectness with 0 mismatches
  through $K_{14}$. Perfectness is local to $K_n$. Next: encode swapped factors, run
  orderly generation in the quotient (where $C(\alpha)=C_2\wr S_n$ is a small natural
  symmetry), and test whether that cheap exact dedup dissolves the over-cap $2^9$ wall —
  validating the full count against $K_{10}$ and $K_{14}\,2^7$ first.
- Until Approach A lands, order-$2$ stays the stated open problem in the manuscript, and
  the deliverable remains $|\mathrm{Aut}|\ge 3$ via the prime/order-4 sweeps.
