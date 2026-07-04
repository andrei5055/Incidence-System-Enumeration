# The order-2 parity theorem — PROVED (2026-06-27)

*Settles the two intractable fully-fixed-point-free over-cap cases by argument instead
of compute: $K_{16}\,2^8=\varnothing$ and $K_{20}\,2^{10}=\varnothing$. Companion to
[the acceleration sketch](order2_acceleration_sketch.md) §3.6, which conjectured this
from data ($K_8{=}0,K_{10}{=}576,K_{12}{=}0,K_{14}{=}\text{yes}$).*

---

## Statement

> **Theorem.** Let $\alpha$ be a fixed-point-free involution (cycle type $2^n$) on the
> vertices of $K_{2n}$. If $2n\ge 6$ and $K_{2n}$ has an $\alpha$-invariant perfect
> $1$-factorisation (P1F), then $n$ is **odd** (equivalently $2n\equiv2\pmod4$).
>
> Contrapositive (the useful direction): **for $n$ even and $2n\ge8$, $K_{2n}$ has no
> $\alpha$-invariant P1F.** The case $2n=4$ is a genuine exception ($K_4$, where the
> unique P1F is $\alpha$-invariant with cycle type $2^2$).

**Consequences.** $K_{16}$ ($n=8$) has **no** $2^8$-invariant P1F and $K_{20}$ ($n=10$)
has **no** $2^{10}$-invariant P1F — the two over-cap types the representative method
could not search. The only genuinely non-empty fully-fpf case is $K_{18}\,2^9$
($n=9$ odd), which the theorem does **not** resolve and which remains the one real
target for the canonicaliser.

---

## Setup and notation

$\alpha$ pairs the $2n$ vertices into $n$ **blocks** $B_i=\{a_i,b_i\}$ with
$\alpha:a_i\leftrightarrow b_i$. Call the $n$ edges $\{a_i,b_i\}$ the **diagonals**;
these are exactly the $\alpha$-**fixed edges** (an edge $\{u,v\}$ has $\alpha\{u,v\}=\{u,v\}$
iff $\alpha u=v$, since $\alpha$ is fixed-point-free).

A P1F of $K_{2n}$ has $2n-1$ factors (perfect matchings), pairwise unions Hamiltonian.
$\alpha$ permutes the factor set; being an involution, it splits them into **fixed**
factors ($\alpha F=F$) and **swapped** pairs $\{G,\alpha G\}$ ($\alpha G\ne G$). Let
$k$ = #fixed factors and $s$ = #swapped pairs, so $k+2s=2n-1$; in particular **$k$ is
odd** (so $k\ge1$).

Two elementary facts used below:

- **(D1) Diagonals lie only in fixed factors.** If a diagonal $d\in G$ then
  $d=\alpha d\in\alpha G$; if $G$ were swapped, $G$ and $\alpha G$ would share the edge
  $d$, contradicting edge-disjointness. So every diagonal sits in a fixed factor, and
  $\sum_{\text{fixed }F} d_F=n$ where $d_F$ = #diagonals in $F$. **Swapped factors
  contain no diagonal.**
- **(D2) Per-factor diagonal parity.** A fixed factor is a disjoint union of
  $\alpha$-edge-orbits: a diagonal (covers $1$ block) or an off-diagonal pair (covers
  $2$ blocks). Covering all $n$ blocks gives $d_F+2e_F=n$, so $d_F\equiv n\pmod2$.

---

## The key lemma

> **Lemma (parity of the antipodal rotation).** Let $H$ be a $2n$-cycle with its two
> perfect matchings $M,M'$ (the alternate-edge $1$-factors). The antipodal rotation
> $\rho$ (rotation by $n$) maps $M\mapsto M'$ **iff $n$ is odd**; for $n$ even it fixes
> each of $M,M'$ setwise.

*Proof.* Index the vertices $0,1,\dots,2n-1$ around $H$; then
$M=\{(2t,2t{+}1)\}$, $M'=\{(2t{+}1,2t{+}2)\}$, and $\rho:i\mapsto i+n\pmod{2n}$.
$\rho$ sends the edge $(2t,2t{+}1)$ to $(2t{+}n,\,2t{+}1{+}n)$. If $n$ is even the first
coordinate $2t+n$ is even, so the image is in $M$; if $n$ is odd it is odd, so the image
is in $M'$. $\square$

This is the whole parity phenomenon. Everything else is bookkeeping.

---

## Proof of the theorem

Assume $n$ is **even** and an $\alpha$-invariant P1F exists; we derive $2n\le4$.

### Step 1 — no swapped pairs.

Suppose $\{G,\alpha G\}$ is a swapped pair. Their union $H=G\cup\alpha G$ is a Hamiltonian
$2n$-cycle (perfectness, as $G\ne\alpha G$). $\alpha$ maps $H$ to itself, so
$\alpha|_H\in\mathrm{Aut}(H)=D_{2n}$, and it is a fixed-point-free involution (it is the
restriction of the globally fixed-point-free $\alpha$, and $H$ spans every vertex).

The fixed-point-free involutions of $D_{2n}$ are the antipodal rotation $\rho$ (no fixed
edge) and the $n$ edge-midpoint reflections (each fixing **exactly two** edges). By
**(D1)** the matchings $G,\alpha G$ contain no diagonal, so $H$ has no $\alpha$-fixed
edge; hence $\alpha|_H$ fixes no edge and must be $\rho$.

But $\alpha$ *swaps* the two $1$-factors of $H$ (namely $G$ and $\alpha G$), so $\rho$
swaps $M,M'$ — which by the Lemma forces **$n$ odd**, contradicting $n$ even. Therefore
$s=0$: **every factor is $\alpha$-fixed**, and $k=2n-1$.

### Step 2 — the diagonal count cannot close.

All $2n-1$ factors are fixed. For any two fixed factors $F,F'$, their union $H'$ is a
Hamiltonian $2n$-cycle and $\alpha|_{H'}$ is again a fixed-point-free element of
$D_{2n}$ — so it is $\rho$ (no fixed edge) or an edge-midpoint reflection (two fixed
edges). Its fixed edges are precisely the diagonals lying in $F\cup F'$. Hence

$$d_F+d_{F'}\in\{0,2\}\qquad\text{for every pair of (fixed) factors.}\tag{$\star$}$$

By **(D1)** $\sum_F d_F=n\ge4$, so some factor has $d_F\ge2$. Among $2n-1\ge7$ factors
carrying only $n$ diagonals in total, a factor with $d_F=0$ certainly exists (factors
with a diagonal number at most $n/2<2n-1$). Pairing the $d_F\ge2$ factor with a
$0$-factor in $(\star)$ forces $d_F=2$; pairing two factors that each had $\ge2$ would
give a sum $\ge4$, so **at most one** factor has $d_F=2$ and all others have $d_F=0$.
Then $\sum_F d_F\le2<n$ — contradiction.

Hence no $\alpha$-invariant P1F exists for $n$ even with $2n\ge8$. $\blacksquare$

*(The argument also locates the lone exception: $n=2$ makes $\sum d_F=2$ achievable by a
single factor with $d_F=2$ — the all-diagonals matching $\alpha$ itself — which is the
$K_4$ P1F, cycle type $2^2$.)*

---

## What this gives, and the matching structure for $n$ odd

The two steps are independent levers, both flowing from the Lemma:

| $n$ | Lemma says | swapped pairs $s$ | structure |
|---|---|---|---|
| even | $\rho$ fixes both $1$-factors | $s=0$ (Step 1) forces all-fixed; then $(\star)$ + count fails | **no P1F** ($2n\ge8$) |
| odd  | $\rho$ swaps the two $1$-factors | $s\ge0$ allowed | P1Fs exist (e.g. $K_6,K_{10},K_{14}$) |

The complementary structure for $n$ odd is forced and matches the data: with swapped
pairs allowed, $(\star)$ together with $\sum d_F=n$ (and $d_F\equiv n\equiv1\pmod2$, so
each $d_F$ is odd $\ge1$) gives $d_F=1$ for every fixed factor and exactly $k=n$ fixed
factors. Verified by direct enumeration (`scratchpad/parity_check.cpp`):

```
K4  (n=2): 1 P1F,  k=3=2n-1,  shape {0,0,2}     (the n=2 exception)
K6  (n=3): 4 P1Fs, k=3=n,     shape {1,1,1}
K8  (n=4): 0
K10 (n=5): 576,    k=5=n,     shape {1,1,1,1,1}
K12 (n=6): 0
```

`k==1` never occurs (Step 1 already rules out the "diagonal matching as the sole fixed
factor" scenario, since that needs $n-1$ swapped pairs). All shapes match the theorem.

---

## Remarks

- **Same parity as the order-4 observation** (paper §5.3). The order-4 analogue —
  $4^k$-invariant P1Fs only for $2n\equiv2\pmod4$ — is the natural next target; the
  same two ingredients should apply (the involution $\alpha^2$ of type $2^{2k}$ is a
  fixed-point-free involution, so the Lemma constrains its swapped factor pairs, and the
  order-4 elements further refine the diagonal/orbit count). Left for a follow-up.
- **The canonicaliser is now needed only for $K_{18}\,2^9$.** With $K_{16}\,2^8$ and
  $K_{20}\,2^{10}$ proved empty, the Approach-A engine has a single live target rather
  than three.
- **Scope.** The theorem concerns the *fully* fixed-point-free type $2^n$. Mixed types
  with fixed points (e.g. $K_{18}\,2^8 1^2$) are separate and unaffected.
