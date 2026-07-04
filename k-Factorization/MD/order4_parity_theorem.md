# Parity for automorphisms of order divisible by 4 — PROVED (2026-06-27)

*Settles paper open problem #1 (§5.3) and generalises it: a P1F of $K_{2n}$ admits an
automorphism whose order is **divisible by $4$** (orders $4,8,12,16,\dots$) only if
$2n\equiv2\pmod4$. The proof reduces to the
[order-2 parity theorem](order2_parity_theorem.md) via the involution $\beta^{\,\mathrm{ord}/2}$,
plus one new general lemma about involution automorphisms and a parity count. Companion
to that note; same underlying mechanism (a 2-torsion obstruction that fires when $n$ is
even). The order-$4$ statement is the smallest instance and is presented first.*

---

## Statement

> **Theorem (order-4 parity).** For $2n\ge8$, if a perfect $1$-factorisation of
> $K_{2n}$ admits an automorphism of order $4$, then $n$ is **odd** (equivalently
> $2n\equiv2\pmod4$). In particular $K_{16}$ and $K_{20}$ admit **no** order-$4$
> automorphism of any P1F.

This was the empirical "Observation" of paper §5.3
($\#\{K_{14}\}=5,\#\{K_{16}\}=0,\#\{K_{18}\}=179,\#\{K_{20}\}=0$); it is now a theorem.

---

## A general lemma (the new ingredient)

> **Lemma (involutions fix $\le2$).** Let $\gamma\ne\mathrm{id}$ be an involution
> automorphism of a P1F of $K_{2n}$, $n\ge2$. Then $\gamma$ fixes **at most $2$**
> vertices.

*Proof.* $\gamma$ permutes the $2n-1$ factors in orbits of size $1$ or $2$. As $2n-1$
is odd, the number of $\gamma$-fixed factors is odd, hence $\ge1$. A **$\gamma$-invariant
pair** of factors exists: if there are $\ge2$ fixed factors, take two of them; otherwise
the single fixed factor leaves $2n-2\ge2$ factors in $\gamma$-swapped pairs, take one.
For such a pair $\{F,F'\}$ the union $H=F\cup F'$ is a Hamiltonian $2n$-cycle
(perfectness, $F\ne F'$) and $\gamma(H)=H$, so $\gamma|_H\in\mathrm{Aut}(H)=D_{2n}$.
It is a non-identity involution (if $\gamma|_H=\mathrm{id}$ it fixes every vertex of $H$,
i.e. all $2n$ vertices, so $\gamma=\mathrm{id}$). Every involution in $D_{2n}$ — the
antipodal rotation ($0$ fixed vertices) or a reflection ($0$ or $2$ fixed vertices) —
fixes $0$ or $2$ vertices of $H$. Since $V(H)$ is all $2n$ vertices, $\gamma$ fixes
$\le2$ vertices globally. $\square$

The bound is **tight**: type $2^{n-1}1^2$ involutions ($f=2$) do occur, e.g.
$K_8\,2^31^2$ ($16$ invariant P1Fs) and $K_{10}\,2^41^2$ ($96$) — verified below.

---

## Proof of the theorem

Let $\beta$ be an order-$4$ automorphism of a P1F of $K_{2n}$, with cycle type
$4^a2^b1^c$ ($a\ge1$, $4a+2b+c=2n$). Put $\alpha=\beta^2$: an involution of type
$2^{2a}1^{f}$ with $f=2b+c$ fixed points. Note $c=2n-4a-2b$ is even, so **$f$ is even**,
and $2n=4a+f$.

Assume $n$ is **even**; we derive a contradiction.

1. **$f\le2$** by the Lemma applied to $\gamma=\alpha$ (which is a non-identity
   involution automorphism: $\beta$ has order $4$, so $\alpha=\beta^2\ne\mathrm{id}$).
   With $f$ even, $f\in\{0,2\}$.

2. **$f=2$ is impossible for $n$ even.** Then $2n=4a+2$, i.e. $n=2a+1$ is odd —
   contradicting $n$ even. (Pure arithmetic; no P1F needed.)

3. **$f=0$ is impossible for $n$ even.** Then $\beta$ has type $4^a$ and
   $\alpha=\beta^2$ is the **fixed-point-free** involution of type $2^{2a}=2^{n}$ with
   $n=2a$. The P1F is $\beta$-invariant, hence $\alpha$-invariant. But $n$ is even and
   $2n\ge8$, so the [order-2 parity theorem](order2_parity_theorem.md) says **no**
   $\alpha$-invariant P1F exists — contradiction.

Both cases fail, so no order-$4$ automorphism exists when $n$ is even ($2n\ge8$);
equivalently, an order-$4$ automorphism forces $n$ odd. $\blacksquare$

The fully fixed-point-free order-$4$ types ($K_{16}\,4^4$, $K_{20}\,4^5$) are exactly the
$f=0$ branch and die directly through $\beta^2$; the types with fixed points die by the
Lemma + arithmetic ($f=2$) or again by $\beta^2$.

---

## Generalization: any order divisible by 4

The same three steps work verbatim for every order $\equiv0\pmod4$, with $\beta^2$
replaced by the central involution $\beta^{\,\mathrm{ord}/2}$.

> **Theorem (order-$4k$ parity).** For $2n\ge8$, if a P1F of $K_{2n}$ admits an
> automorphism whose order is divisible by $4$, then $n$ is **odd**. In particular
> $K_{16}$ and $K_{20}$ admit no P1F automorphism of order $4,8,12,16,\dots$.

*Proof.* Let $\beta$ have order $4s$ and put $\gamma=\beta^{2s}$ — an involution
($\gamma\ne\mathrm{id}$ since $\mathrm{ord}(\gamma)=4s/\gcd(4s,2s)=2$). $\gamma$ fixes a
point of an $L$-cycle of $\beta$ iff $L\mid 2s$; a **moved** point lies in a cycle with
$L\mid4s,\ L\nmid2s$, i.e. $v_2(L)=v_2(4s)\ge2$, so $4\mid L$. Hence the number of moved
points is $\equiv0\pmod4$, and the fixed count is $f=2n-(\text{moved})\equiv2n\pmod4$.
The Lemma gives $f\le2$, and $f$ is even (both $2n$ and the moved count are), so
$f\in\{0,2\}$.

- $f=2$: $2n\equiv2\pmod4$, i.e. $n$ odd.
- $f=0$: $\gamma$ is the fixed-point-free involution $2^n$ with $2n\equiv0\pmod4$
  ($n$ even); the P1F is $\gamma$-invariant, so the order-2 theorem forbids it.

Either way $n$ is odd. $\square$

This subsumes the order-$4$ theorem ($s=1$, $\gamma=\beta^2$) and the requested orders
$8$ ($s=2$, $\gamma=\beta^4$) and $12$ ($s=3$, $\gamma=\beta^6$). Note the dichotomy also
pins the cycle structure: $\gamma$'s moved cycles are $4$-divisible, so e.g. an order-$8$
automorphism has **no $4$-cycles** ($\beta^4$ would fix them, $\ge4$ fixed points) and an
order-$12$ one has **no $3$- or $6$-cycles** — confirmed in the validation below.

**Scope of the method.** Orders $\equiv2\pmod4$ ($6,10,\dots$) do **not** reduce this
way: there $\beta^{\,\mathrm{ord}/2}$'s moved-point count is only $\equiv2\pmod4$, the
$f=0$ branch can land on $n$ odd, and the residual case is the *mixed* involution type
$2^{n-1}1^2$ — exactly the still-open order-2-with-fixed-points problem. Odd orders have
no involution power and are unconstrained by parity (e.g. order $3$ occurs at $K_{16}$).
So "order divisible by $4\Rightarrow n$ odd" is the full reach of the $\beta^2$/order-2
lever.

---

## Validation (`scratchpad/order4_check.cpp`)

Enumerating $\langle g\rangle$-invariant P1Fs by exact cover over $\langle g\rangle$-edge-
orbits, for $g$ of the stated cycle type:

```
order-2 sanity:   K8  2^4    -> 0      K10 2^5    -> 576
LEMMA (f>=4):     K8  2^2 1^4-> 0      K8  2 1^6  -> 0
                  K10 2^2 1^6-> 0      K12 2^4 1^4-> 0
bound tight (f=2):K8  2^3 1^2-> 16     K10 2^4 1^2-> 96
ORDER-4, n EVEN:  K8  4^2    -> 0      K8  4 2^2  -> 0    K8 4 2 1^2 -> 0   K8 4 1^4 -> 0
                  K12 4^3    -> 0      K12 4^2 2^2-> 0    K12 4^2 1^4-> 0
ORDER-4, n ODD:   K14 4^3 1^2-> 1152   K10 4^2 1^2-> 16   K10 4^2 2 -> 16   K6 4 1^2 -> 2
ORDER-8, n EVEN:  K8  8 -> 0   K16 8^2 -> 0   K16 8 4^2 -> 0   K16 8 4 2^2 -> 0
ORDER-8, n ODD:   K18 8^2 1^2 -> 1792   K18 8^2 2 -> 1280
ORDER-12,n EVEN:  K12 12 -> 0   K16 12 4 -> 0   K20 12 4^2 -> 0   K16 12 2^2 -> 0
                  K18 12 3^2 -> 0   K18 12 6 -> 0    (3-/6-cycle types forbidden)
ORDER-12,n ODD:   K14 12 1^2 -> 24
```

Every $n$-even type of order $4,8,12$ gives $0$ — including the lemma-forbidden
structures (order-8 with $4$-cycles, order-12 with $3$- or $6$-cycles, or too many fixed
points). The orders genuinely occur for $n$ odd (so the theorem is not vacuous):
$K_{14}\,4^31^2=1152$ (paper: $5$ iso-classes), $K_{18}\,8^21^2=1792$, $K_{14}\,12\,1^2=24$.

---

## Remarks

- **Both paper parity observations are now theorems.** Order-2
  ([here](order2_parity_theorem.md)) and order-4 (this note); §5.3's "structural proof …
  left open" is closed.
- **The Lemma is of independent interest** — it is a statement about *all* involution
  automorphisms of P1Fs, not just $\beta^2$: such an involution has cycle type $2^n$
  (fixed-point-free) or $2^{n-1}1^2$ (exactly two fixed points). Nothing with $\ge4$
  fixed points can be an automorphism of a perfect $1$-factorisation.
- **Latin-square reading.** Via the atomic-Latin-square correspondence (P1F of $K_{2n}$
  $\leftrightarrow$ atomic Latin square of order $2n-1$), this is the same parity family
  as cyclic-group complete mappings (transversals of $\mathbb Z_m$ exist iff $m$ odd):
  in every case a $2$-torsion / summation obstruction blocks the even case.
- **Higher orders — done.** The generalization above proves *all* orders divisible by
  $4$ at once (validated for $8$ and $12$). Orders $\equiv2\pmod4$ remain open via the
  mixed involution type $2^{n-1}1^2$; odd orders are unconstrained.
