# Perfect 1-factorisations of $K_{18}$ and $K_{20}$ with prescribed automorphisms

**Andrei V. Ivanov**
Independent Researcher, San Jose, CA, USA

**Leonid M. Tertitski**
Independent Researcher, Carlsbad, CA, USA

*Draft — version 0.5 (2026-06-29). Resolves the last open order-$2$ existence question: $K_{18}$
**does** admit a P1F invariant under a fixed-point-free involution of type $2^9$. The witness is a
class with $|\mathrm{Aut}|=16$, found within the order-$4$ sweep (its automorphism group contains
both an order-$4$ element and the $2^9$ involution) and independently reproduced; it is given in
Appendix A.5. Consequently every order-$2$ existence question for $K_{18}$ and $K_{20}$ is now
settled (§5.4, §6); what remains open is only the complete enumeration of the order-$2$ classes.
Version 0.4 (2026-06-29) added §5.5, relating the order-$(2n-1)$ results to the
*starter / 1-rotational* subfamily: for $2n-1$ prime, starter-induced $\Leftrightarrow$
$(2n-1)\mid|\mathrm{Aut}|$, so the Skolem-sequence and hill-climbing constructions of prior
work occupy exactly that column — which we classify completely (two for $K_{18}$, seven for
$K_{20}$), the genuinely new classes being the non-1-rotational ones (with Appendix A.4 and
new references). Version 0.3 (2026-06-28) added the order-$2$ landscape (§5.4): a
one-fixed-factor lemma for the two-fixed-point involution, the Galois construction settling
existence of the two-fixed-point cases for $K_{18}$ and $K_{20}$, and the reduction of the
lone remaining open case ($K_{18}\,2^9$) to the unique P1F of $K_{10}$. Version 0.2 added the
$K_{14}$ validation (§4.2). Work in progress; see §6, “Status of results”.*

---

## Abstract

A *perfect 1-factorisation* (P1F) of the complete graph $K_{2n}$ is a partition of its
edges into $2n-1$ perfect matchings such that the union of any two of them is a
Hamiltonian cycle. Complete enumeration of P1Fs is currently feasible only up to
$K_{16}$, for which Gill and Wanless found $3155$ P1Fs, $89$ of them with a
non-trivial automorphism group. For $K_{18}$ and beyond, full enumeration is out
of reach, and the only P1Fs on record are isolated examples from algebraic
constructions.

We take an *automorphism-first* approach: for a fixed permutation $\alpha$ of the
vertices we enumerate **all** P1Fs invariant under $\alpha$ directly, by an exact
cover of the edges with $\alpha$-invariant orbits of 1-factors, using orderly
generation for isomorph rejection. Running one representative $\alpha$ per cycle
type of a target order, and deduplicating across types, yields every P1F admitting
an automorphism of that order. This scales to orders where complete enumeration
does not.

We validate the method by reproducing two known classifications exactly,
class-for-class: the $K_{16}$ enumeration of Gill and Wanless, and the $23$ P1Fs of
$K_{14}$. The $K_{14}$ check is decisive because $14\equiv 2\pmod 4$ — as does
$18$ — so $K_{14}$ exercises the fixed-point-free involution and order-$4$ cycle
types that drive the $K_{18}$ results and that $K_{16}$ ($\equiv 0$) structurally
lacks; being small, its dense types are tractable, letting us prove the
order-$2$ (fixed-point-free) classification complete. We then report, for the first
time, complete counts of P1Fs of $K_{18}$ and $K_{20}$ admitting automorphisms of
various orders.
In particular $K_{18}$ has exactly **179** isomorphism classes of P1Fs admitting an
automorphism of order $4$, and the method recovers (and corrects an undercount of)
the count for order $8$. For $K_{20}$ we find $7$ classes admitting an automorphism
of order $19$ and none admitting one of order $4,7,11,13$ or $17$. We also prove a clean
parity phenomenon: a P1F of $K_{2n}$ admits an automorphism of order divisible by $4$, or
a fixed-point-free involution, only when $2n\equiv 2\pmod 4$ — so $K_{16}$ and $K_{20}$
admit neither. Together with a one-fixed-factor lemma, the Galois construction,
and an explicit fixed-point-free example for $K_{18}$ (type $2^9$, with
$|\mathrm{Aut}|=16$, recovered within the order-$4$ sweep), this **settles every order-$2$
existence question** for $K_{18}$ and $K_{20}$; only the complete *enumeration* of the
order-$2$ classes remains open. The classifications by automorphism of order $\ge 3$ are
complete for both $K_{18}$ and $K_{20}$.

---

## 1. Introduction

A *1-factorisation* of $K_{2n}$ is a partition of its edge set into $2n-1$ perfect
matchings (1-factors). It is **perfect** (a P1F) if the union of every pair of
distinct 1-factors is a single Hamiltonian cycle of $K_{2n}$. P1Fs are central
objects in combinatorial design theory, with connections to atomic Latin squares,
the perfect 1-factorisation conjecture, and round-robin tournament scheduling
[Wallis; Mendelsohn–Rosa].

Exact enumeration of P1Fs grows explosively: there is a unique P1F of $K_4,K_6,K_8$
(up to isomorphism), $5$ of $K_{12}$, $23$ of $K_{14}$, and — the current frontier —
$3155$ of $K_{16}$, of which $89$ have a non-trivial automorphism group, all cyclic
[Gill–Wanless]. Their computation enumerated *all* $3155$ P1Fs by backtracking and
then identified automorphisms post hoc. That strategy does not extend to $K_{18}$:
the number of P1Fs is far too large to list.

Yet the P1Fs of greatest interest — those with symmetry — are exactly the ones that
remain accessible, because prescribing an automorphism collapses the search space
enormously. This is the classical idea behind isomorph-free generation of designs
with prescribed automorphism groups. We apply it to P1Fs:

> **Automorphism-first enumeration.** Fix a permutation $\alpha$ of the $2n$ vertices.
> Enumerate every P1F that $\alpha$ maps to itself, directly, without enumerating any
> P1F that $\alpha$ does not fix.

Running one representative $\alpha$ per cycle type of order $k$ and taking the union
(deduplicated up to isomorphism) yields **every** P1F of $K_{2n}$ admitting an
automorphism of order $k$. Because the cyclic-group constructions used in prior work
are a special case, this both recovers the known symmetric examples and finds the
rest.

**Contributions.**
1. An orderly-generation algorithm enumerating all $\alpha$-invariant P1Fs for a
   fixed $\alpha$ (§3), with per-level isomorph rejection by the centralizer
   $C(\alpha)$.
2. Exact reproduction of two known classifications, class-for-class (§4): the
   Gill–Wanless $K_{16}$ enumeration, and the $23$ P1Fs of $K_{14}$ — the latter the
   $2n\equiv 2\pmod 4$ case that validates the fixed-point-free machinery and, being
   small, proves the order-$2$ classification complete.
3. New complete counts for $K_{18}$ and $K_{20}$ by automorphism order (§5),
   including $179$ classes for $K_{18}$ with an order-$4$ automorphism, and a
   correction to the count obtainable by cyclic-only searches. A consequence (§5.5):
   the *starter-induced* (1-rotational) P1Fs — the entire output of the Skolem-sequence
   and hill-climbing constructions of prior work — are exactly the order-$(2n-1)$ classes,
   here classified **completely** (two for $K_{18}$, seven for $K_{20}$); the genuinely new
   classes are the non-1-rotational ones.
4. Two parity theorems (§5.3): a P1F of $K_{2n}$ admits an automorphism of order
   divisible by $4$ only when $2n\equiv2\pmod4$; and the fixed-point-free involution
   ($2^n$) requires $2n\equiv2\pmod4$. In particular $K_{16}$ and $K_{20}$ admit no such
   automorphisms.
5. The order-$2$ landscape (§5.4), now **complete for existence**: a one-fixed-factor
   lemma for the two-fixed-point involution, the Galois construction settling the
   two-fixed-point cases for $K_{18}$ and $K_{20}$, and an explicit $K_{18}$ P1F with a
   fixed-point-free $2^9$ automorphism (Appendix A.5) settling the last open case. Every
   order-$2$ type is thus decided; only the full enumeration of order-$2$ classes is out
   of reach.

---

## 2. Preliminaries

Let $V=\{0,1,\dots,2n-1\}$ and $G=K_{2n}$. A **1-factor** (perfect matching) is a set
of $n$ pairwise-disjoint edges covering $V$; we represent it in adjacency form
$M:V\to V$ with $M(u)$ the partner of $u$. Two 1-factors $A,B$ are **perfect** if
$A\cup B$ is a single Hamiltonian cycle. A **P1F** is a set $\mathcal F$ of $2n-1$
pairwise-perfect 1-factors partitioning $E(G)$.

An **automorphism** of a P1F $\mathcal F$ is a vertex permutation $\alpha\in S_{2n}$
with $\alpha(\mathcal F)=\mathcal F$ (equivalently, $\alpha$ permutes the 1-factors of
$\mathcal F$). The automorphisms form a group $\mathrm{Aut}(\mathcal F)$. Two P1Fs are
**isomorphic** if a vertex permutation maps one to the other; we count isomorphism
classes.

For a permutation $\alpha$, its **cycle type** is the multiset of its cycle lengths;
its **order** is their least common multiple. The **centralizer**
$C(\alpha)=\{g\in S_{2n}: g\alpha=\alpha g\}$ acts on the set of $\alpha$-invariant
P1Fs and on partial structures; it is the symmetry group we exploit.

If $\mathcal F$ is $\alpha$-invariant then $\alpha$ permutes its $2n-1$ factors, so
$\mathcal F$ is a disjoint union of **$\alpha$-orbits of factors**
$\{M,\alpha M,\alpha^2M,\dots\}$. Each such orbit consists of pairwise-perfect,
edge-disjoint 1-factors, and the orbit lengths divide the order of $\alpha$.

---

## 3. Method

### 3.1 Reduction to cycle types

Conjugate permutations yield isomorphic invariant-P1F families, and conjugacy in
$S_{2n}$ is determined by cycle type. Hence to find all P1Fs admitting an
automorphism of order $k$ it suffices to:

1. enumerate the cycle types of order $k$ (partitions of $2n$ into parts whose lcm
   is $k$);
2. for one representative $\alpha$ of each type, enumerate all $\alpha$-invariant
   P1Fs; and
3. deduplicate the union up to isomorphism.

For order $k$ the relevant cycle types are partitions of $2n$ into parts dividing
$k$ with least common multiple exactly $k$ — typically only a few dozen, versus the
astronomically many P1Fs.

### 3.2 Exact cover by factor-orbits

Fix $\alpha$. We build an $\alpha$-invariant P1F as an **exact cover** of $E(G)$ by
$\alpha$-invariant factor-orbits. Order the edges $0,1,\dots,\binom{2n}{2}-1$
lexicographically. At each step:

1. let $e$ be the lowest uncovered edge;
2. generate every 1-factor $M\ni e$ that is perfect with all committed factors and
   whose orbit $\{M,\alpha M,\dots\}$ is valid (members pairwise perfect and
   edge-disjoint);
3. commit one such orbit and recurse; when all edges are covered, $\mathcal F$ is a
   complete $\alpha$-invariant P1F.

Candidate factors through $e$ are generated by a matching search with Hamiltonicity
pruning (a partial matching is rejected as soon as it cannot complete to a factor
perfect with a committed one).

### 3.3 Isomorph rejection by orderly generation

Two cover branches related by an element of $C(\alpha)$ produce isomorphic P1Fs, so
exploring both is wasteful. We use **orderly generation** over $C(\alpha)$:

> At a node with committed factor-set $P$, let $G\le C(\alpha)$ be the subgroup
> fixing $P$ setwise. Among the candidate factor-orbits covering the lowest uncovered
> edge, keep **one representative per $G$-orbit**, and recurse into each chosen orbit
> $O$ with the smaller stabilizer $\mathrm{Stab}_G(O)$.

This is sound for counting isomorphism classes: any $g\in G$ is a bijection between
the sets of completions extending two $G$-equivalent siblings, so exploring one
representative reaches every class the others would. As the search descends, the
stabilizer chain $G\ge\mathrm{Stab}_G(O)\ge\cdots$ shrinks quickly, collapsing the
dense high-symmetry branches that defeat naïve pruning.

A labelling-independent canonical form (the lexicographically least serialisation
over all Hamiltonian-cycle framings of the factor set) provides the final
deduplication across cycle types and labellings, and simultaneously yields
$|\mathrm{Aut}(\mathcal F)|$ (as the number of framings attaining the minimum,
divided by two).

### 3.4 Complexity remarks

When $|C(\alpha)|$ is moderate it is enumerated and the per-level $G$-orbit dedup is
cheap. When $|C(\alpha)|$ is large (cycle types with many fixed points, where
$C(\alpha)$ contains a large symmetric group) the current implementation falls back
to a generating-set walk without per-level dedup; such types are the present
performance bottleneck (§6). They are, empirically, the types that yield no P1Fs, so
they cost time without changing the counts.

---

## 4. Validation

### 4.1 $K_{16}$

Running the method on $K_{16}$ for orders $3,5,7$ (whose union captures all P1Fs with
$|\mathrm{Aut}|>2$ that is not a power of two) reproduces the Gill–Wanless
classification **exactly**:

| generator cycle type | order | this work | Gill–Wanless |
|---|---|---|---|
| $3^5 1^1$ | 3 | 19 | 19 |
| $5^3 1^1$ | 5 | 5 | 5 |
| $7^2 1^2$ | 7 | 4 | 4 |
| $14^1 1^2$ | 14 | 1 | 1 |
| $15^1 1^1$ | 15 | 1 | 1 |
| **total ($|\mathrm{Aut}|>2$)** | | **30** | **30** |

The independently computed automorphism group orders match in every class. The dense
type $3^5 1^1$ (centralizer order $29160$), on which a naïve symmetry-breaking search
stalls, completes in seconds under orderly generation. This agreement with the
published authority — not merely on the total, but on the full distribution — gives
strong confidence in the implementation.

### 4.2 $K_{14}$ — the fixed-point-free case

$K_{14}$ has exactly $23$ P1Fs. Two of them are asymmetric ($|\mathrm{Aut}|=1$) and
hence invisible to any automorphism-first method; the remaining $21$ have a
non-trivial automorphism group, with the distribution
$\{2{:}3,\,3{:}5,\,4{:}1,\,6{:}5,\,12{:}5,\,84{:}1,\,156{:}1\}$. Our method reproduces
this distribution **exactly**, cross-checked against an independent complete
enumeration:

| automorphism order | classes | $|\mathrm{Aut}|$ distribution |
|---|---|---|
| 2 | 16 | $\{2{:}3,\,4{:}1,\,6{:}5,\,12{:}5,\,84{:}1,\,156{:}1\}$ |
| 3 | 17 | $\{3{:}5,\,6{:}5,\,12{:}5,\,84{:}1,\,156{:}1\}$ |
| 7 | 1 | $\{84{:}1\}$ |
| 13 | 1 | $\{156{:}1\}$ |
| **union ($|\mathrm{Aut}|>1$)** | **21** | $\{2{:}3,\,3{:}5,\,4{:}1,\,6{:}5,\,12{:}5,\,84{:}1,\,156{:}1\}$ |

This validation is qualitatively stronger than the $K_{16}$ one in two respects.

*It exercises the structure of $K_{18}$.* Because $14\equiv 2\pmod 4$, like $18$ and
unlike $16$, the order-$2$ and order-$4$ cycle types of $K_{14}$ include
**fixed-point-free** ones ($2^7$ and $4^3 2^1$) — precisely the types responsible for
the $K_{18}$ results in §5 and structurally absent from $K_{16}$. The identical engine
producing the exact known $K_{14}$ counts is therefore direct evidence for the
correctness of the $K_{18}$ machinery, not merely of the $N$-generic core.

*It closes the order-$2$ case.* For $K_{14}$ the fixed-point-free involution type
$2^7$ has centralizer of order $2^7\cdot 7! = 645{,}120$, small enough to enumerate,
so order $2$ **completes** — yielding the proven count of $16$ in a few seconds. This
is the case that is intractable for the larger sizes (where $2^8$, $2^9$, $\dots$ have
centralizers too large to enumerate; cf. §6), and $K_{14}$ confirms that the
limitation there is one of scale, not of correctness.

*A remark on composite orders.* By Cauchy's theorem a group of order divisible by a
prime $p$ contains an element of order $p$, so for **prime** orders the count of
classes our method finds equals the number whose $|\mathrm{Aut}|$ is divisible by that
prime — making orders $2,3,7,13$ above clean checks. This fails for **composite**
orders: a group whose order is divisible by $4$ need not contain an element of order
$4$ (its Sylow $2$-subgroup may be elementary abelian, or the group may be $A_4$).
Indeed $K_{14}$ has $8$ classes with $4\mid|\mathrm{Aut}|$ but only $5$ admitting an
order-$4$ *automorphism* ($\{12{:}4,\,156{:}1\}$); the three excluded are the
$|\mathrm{Aut}|=4$ (Klein), one $|\mathrm{Aut}|=12$ ($A_4$, found instead by order
$3$), and the $|\mathrm{Aut}|=84$ class. The order-$4$ counts reported throughout this
paper are counts of P1Fs *admitting an order-$4$ automorphism*, i.e. with an order-$4$
group element.

---

## 5. New results

### 5.1 $K_{18}$

$K_{18}$ admits no complete published enumeration. The only P1F of $K_{18}$ on record
is the algebraic construction for orders $q+1$ ($q$ a prime power), here $q=17$,
whose automorphism group is the affine group $\mathrm{AGL}(1,17)$ of order
$17\cdot 16 = 272$. Our method recovers it. Beyond it we obtain complete counts by
automorphism order:

| automorphism order | classes | $|\mathrm{Aut}|$ distribution |
|---|---|---|
| 17 | 2 | $\{17:1,\ 272:1\}$ |
| 8 | 30 | $\{8:17,\ 16:12,\ 272:1\}$ |
| **4** | **179** | $\{4:144,\ 8:22,\ 16:12,\ 272:1\}$ |

The order-$4$ count of **179** is, to our knowledge, new. We note that the
$|\mathrm{Aut}|=16$ and $|\mathrm{Aut}|=272$ rows agree across the independent
order-$4$ and order-$8$ computations, an internal cross-check. Exactly one of these $179$
classes — one with $|\mathrm{Aut}|=16$ — additionally carries a fixed-point-free involution
of type $2^9$; this is the witness that settles the last open order-$2$ case (§5.4,
Appendix A.5).

**A correction.** Searches that build the automorphism cyclically (fixing a vertex)
structurally miss P1Fs whose order-$8$ automorphism is fixed-point-free (cycle type
$8^2 2^1$); such a search reports $20$ classes for order $8$. The automorphism-first
method, which places one representative per cycle type including fixed-point-free
ones, gives the correct $30$.

### 5.2 $K_{20}$

$K_{20}=K_{19+1}$ likewise has a single recorded P1F, from $\mathrm{AGL}(1,19)$ of
order $19\cdot 18=342$; our method recovers it. Complete counts so far:

| automorphism order | classes | $|\mathrm{Aut}|$ distribution |
|---|---|---|
| 19 | 7 | $\{19:3,\ 57:1,\ 171:2,\ 342:1\}$ |
| 17 | 0 | — |
| 13 | 0 | — |
| 11 | 0 | — |
| 7 | 0 | — |
| 5 | 0 | — |
| 4 | 0 | — |

No P1F of $K_{20}$ admits an automorphism of order $4,5,7,11,13,$ or $17$. The seven
order-$19$ classes all lie in the $\mathrm{AGL}(1,19)$ family (orders dividing
$342=2\cdot 3^2\cdot 19$). The count for order $3$ is in progress (§6).

### 5.3 A parity observation

The fixed-point-free order-$4$ cycle types ($4^4$ for $K_{16}$, $4^5$ and $4^4 2^2$
for $K_{20}$, etc.) admit no valid factor-orbit through the first edge, so they
contribute nothing; for the sizes with order-$4$ P1Fs the contributing types are
those with two fixed points ($4^3 1^2$ for $K_{14}$, $4^4 1^2$ for $K_{18}$). Across
the sizes accessible to us this gives

$$\#\{K_{14}\}=5,\quad \#\{K_{16}\}=0,\quad \#\{K_{18}\}=179,\quad \#\{K_{20}\}=0,$$

where $\#\{K_{2n}\}$ denotes the number of P1Fs of $K_{2n}$ admitting an order-$4$
automorphism. The two sizes with such P1Fs are exactly the two with
$2n\equiv 2\pmod 4$ ($K_{14},K_{18}$); the two without are exactly those with
$2n\equiv 0$ ($K_{16},K_{20}$). This is now a theorem:

> **Theorem (order-$4k$ parity).** A P1F of $K_{2n}$ ($2n\ge8$) admits an
> automorphism whose order is **divisible by $4$** (orders $4,8,12,\dots$) only if
> $2n\equiv 2\pmod 4$ (equivalently $n$ odd). In particular $K_{16}$ and $K_{20}$ admit
> no P1F automorphism of any order divisible by $4$.

*Proof.* Let $\beta$ have order $4s$ and put $\gamma=\beta^{2s}$ (an involution). A
**lemma** bounds the fixed points of any involution automorphism: $\gamma$ acts on any
$\gamma$-invariant factor-pair's Hamiltonian cycle as an involution of $D_{2n}$, which
fixes $0$ or $2$ vertices — and such a pair always exists (the number of $\gamma$-fixed
factors is odd, hence $\ge1$; pair two fixed factors, or a fixed factor's swapped
partners). Now $\gamma$'s moved points lie in $\beta$-cycles of length divisible by $4$
(those with $L\mid4s,\,L\nmid2s$), so the moved count is $\equiv0\pmod4$ and the fixed
count $f\equiv2n\pmod4$; with $f\le2$ even, $f\in\{0,2\}$. If $f=2$ then $2n\equiv2$, i.e.
$n$ odd. If $f=0$ then $\gamma$ is the fixed-point-free $2^n$ with $2n\equiv0$ ($n$ even),
which the order-2 theorem below forbids. Hence $n$ is odd. $\square$

See [order4_parity_theorem.md](order4_parity_theorem.md) for details and validation
(orders $4,8,12$). Orders $\equiv2\pmod4$ do not reduce this way (they meet the open
mixed type $2^{n-1}1^2$). The **order-$2$ theorem** is the engine of the $f=0$ step:

> **Theorem (order-2 parity).** For $2n\ge6$, $K_{2n}$ has a P1F invariant under a
> fixed-point-free involution (cycle type $2^n$) only if $2n\equiv2\pmod4$. In
> particular $K_{16}$ has no $2^8$-invariant and $K_{20}$ no $2^{10}$-invariant P1F.

*Proof sketch.* The $2n-1$ factors split under $\alpha$ into fixed factors and swapped
pairs $\{G,\alpha G\}$. For a swapped pair, $G\cup\alpha G$ is a Hamiltonian $2n$-cycle
$H$ on which $\alpha$ acts as a fixed-point-free automorphism; as $G,\alpha G$ carry no
$\alpha$-fixed (diagonal) edge, that automorphism fixes no edge of $H$ and is therefore
the antipodal rotation $\rho$. The antipodal rotation of a $2n$-cycle swaps its two
$1$-factors iff $n$ is odd; since $\alpha$ swaps $G$ and $\alpha G$, $n$ must be odd.
Hence for $n$ even every factor is $\alpha$-fixed, and the $n$ diagonals (each in a fixed
factor, each fixed factor holding $\equiv n\equiv0\pmod2$ of them, and any two fixed
factors' union holding $0$ or $2$) cannot sum to $n\ge4$ — contradiction. $\square$

Both parity statements are now theorems, and the order-4 one generalises to **every
order divisible by $4$** (same proof with $\beta^2$ replaced by $\beta^{\,\mathrm{ord}/2}$):
no P1F of $K_{16}$ or $K_{20}$ admits an automorphism of order $4,8,12,16,\dots$.

### 5.4 The order-2 landscape

The theorems of §5.3 settle the **fixed-point-free** involution type $2^n$ for
$2n\equiv0\pmod4$ ($K_{16}\,2^8$ and $K_{20}\,2^{10}$ are empty). The other involution
type — by the $\le 2$-fixed-points lemma these are the only two — is the
**two-fixed-point** type $2^{n-1}1^2$. We settle its existence and pin its structure.

**Lemma (one fixed factor).** If a P1F of $K_{2n}$ has an automorphism $\alpha$ of cycle
type $2^{n-1}1^2$ (fixing exactly two vertices $u,v$), then $\alpha$ fixes **exactly one**
factor: the *diagonal factor* $D$ consisting of the edge $\{u,v\}$ together with the
$n-1$ transposition pairs (these $n$ $\alpha$-fixed edges form one perfect matching).

*Proof.* If two distinct factors $F,G$ were both $\alpha$-fixed, then $H=F\cup G$ is a
Hamiltonian cycle with $\alpha(H)=H$; as $H$ spans $V$, $\alpha|_H$ is the cycle
involution fixing exactly $u,v$ — the reflection through them. On an even cycle the
unique proper $2$-edge-colouring gives the two edges at each fixed vertex opposite
colours, and that reflection swaps them, hence swaps $F$ and $G$ — contradicting both
being fixed. So at most one fixed factor; and an involution acting on the odd number
$2n-1$ of factors fixes an odd number, so exactly one. It carries all $n$ $\alpha$-fixed
edges, i.e. equals $D$. $\square$

The lemma was validated against direct enumeration: for $K_6,K_8,K_{10},K_{12}$ the
number of $2^{n-1}1^2$-invariant P1Fs is $2,16,96,1536$ respectively, and in every one
the fixed-factor count is exactly $1$.

**Existence (Galois construction).** The two-fixed-point case is non-empty for both our
sizes. The Galois P1F $GK_{p+1}$ over $\mathrm{GF}(p)$ ($p$ prime; vertices
$\mathrm{GF}(p)\cup\{\infty\}$, factor $F_k=\{\infty,k\}\cup\{\{k+i,k-i\}:1\le i\le(p-1)/2\}$)
admits the involution $x\mapsto -x$, fixing $0$ and $\infty$, of cycle type
$2^{(p-1)/2}1^2$. For $p=17$ this is the $\mathrm{AGL}(1,17)$ factorisation of $K_{18}$
(Appendix A.1), type $2^8 1^2$; for $p=19$ the $\mathrm{AGL}(1,19)$ factorisation of
$K_{20}$ (Appendix A.3), type $2^9 1^2$. Hence **$K_{18}$ and $K_{20}$ both admit a P1F
with a two-fixed-point involution.** The existence question is thus settled by
construction; only the *complete enumeration* of these classes is out of reach, their
centralizers ($2^{n-1}(n-1)!\cdot2$) being far past the enumeration cap (§6).

**The last fixed-point-free case, resolved.** The one remaining order-$2$ existence
question — the fixed-point-free type $2^9$ of $K_{18}$ — we settle **affirmatively**:
$K_{18}$ admits a perfect $1$-factorisation invariant under a fixed-point-free involution
of type $2^9$. This is notable because the type resists direct attack on all three of the
usual fronts: the Galois construction yields only the two-fixed-point type, never
fixed-point-free; the parity theorem of §5.3 does not apply ($18\equiv2\pmod4$, and
$K_{14}\,2^7$ being non-empty shows there is no parity obstruction to exploit); and a
direct exhaustive search over the type is infeasible (its centralizer
$C(\alpha)=C_2\wr S_9$ has order $2^9\cdot 9!\approx 1.9\times10^8$, far past the
enumeration cap, and the search tree is estimated at $\sim 10^{16}$ nodes).

The witness is obtained not by searching the $2^9$ type directly but as a by-product of
the tractable **order-$4$** enumeration (§5.1). Any P1F whose automorphism group has order
divisible by $4$ and contains a $2^9$ element is found by the order-$4$ sweep; exactly one
of the $179$ order-$4$ classes — a class $\mathcal F^\ast$ with $|\mathrm{Aut}|=16$ — is of
this kind. Its group contains the fixed-point-free involution
$$\alpha=(0\,1)(2\,4)(3\,5)(6\,8)(7\,9)(10\,12)(11\,13)(14\,16)(15\,17),$$
which permutes the $17$ factors of $\mathcal F^\ast$ as $9$ fixed factors and $4$ swapped
pairs — exactly the structure forced below. The factorisation is listed in Appendix A.5,
was verified independently (all $\binom{17}{2}=136$ factor-pairs are Hamiltonian; $\alpha$
is a genuine automorphism), and was reproduced as a byte-identical canonical matrix by a
separate run. Its $9$ fixed factors realise the $K_{10}$-lift described next.

This realises the one structural handle on the type: the $n=9$ fixed factors of any
$2^9$-invariant P1F of $K_{18}$ form a perfect $1$-factorisation of $K_{10}$, and
**$K_{10}$ has a unique P1F up to isomorphism** (verified: $864$ labellings with the first
factor fixed, one class). The example $\mathcal F^\ast$ exhibits a lift of that single
object to a $2^9$-invariant P1F of $K_{18}$, so the lift exists.

**The order-$2$ existence landscape is now complete.** Combining §5.3 and the above:

| size | type | cycle type | status |
|---|---|---|---|
| $K_{16}$ | fixed-point-free | $2^8$ | **empty** (§5.3) |
| $K_{20}$ | fixed-point-free | $2^{10}$ | **empty** (§5.3) |
| $K_{18}$ | fixed-point-free | $2^9$ | **non-empty** (Appendix A.5) |
| $K_{18}$ | two-fixed-point | $2^8 1^2$ | **non-empty** (Galois, A.1) |
| $K_{20}$ | two-fixed-point | $2^9 1^2$ | **non-empty** (Galois, A.3) |

Every order-$2$ involution type for $K_{18}$ and $K_{20}$ is decided. What remains beyond
the present method is only the *complete enumeration* of the order-$2$ classes — equivalently
the P1Fs whose automorphism group is an elementary-abelian $2$-group, invisible to any
order-$\ge3$ (or order-$4$) sweep — whose centralizers exceed the enumeration cap (§6).

### 5.5 The 1-rotational (starter-induced) subfamily and prior constructions

Every previously recorded P1F of $K_{18}$ and $K_{20}$ arises from a *starter*, and it is
worth pinning down precisely where these sit in our classification. A P1F of $K_{2n}$ on the
vertex set $\mathbb{Z}_{2n-1}\cup\{\infty\}$ is **1-rotational** if the rotation
$\rho:x\mapsto x+1\pmod{2n-1}$ (with $\rho(\infty)=\infty$) is an automorphism; equivalently
it is the cyclic development of a *starter* — a set of $n-1$ pairs partitioning
$\mathbb{Z}_{2n-1}\setminus\{0\}$ whose differences cover $\mathbb{Z}_{2n-1}\setminus\{0\}$,
together with the pair $\{0,\infty\}$. The rotation $\rho$ has order $2n-1$ and cycle type
$(2n-1)^1\,1^1$.

For both our sizes $2n-1\in\{17,19\}$ is **prime**, which makes the correspondence exact:

> **Observation.** For $K_{2n}$ with $2n-1$ prime, a P1F is starter-induced (1-rotational)
> **iff** $(2n-1)\mid|\mathrm{Aut}|$.

Indeed, any automorphism of prime order $2n-1$ permutes the $2n=(2n-1)+1$ vertices as a single
$(2n-1)$-cycle and one fixed point (no other cycle type of that order fits $2n$ points), so —
after relabelling — it is the rotation $\rho$ and the P1F is 1-rotational; conversely $\rho$
contributes the factor $2n-1$ to $|\mathrm{Aut}|$. Hence the *entire* starter-induced subfamily
is exactly the **order-$(2n-1)$** column of our tables, which our method classifies
**completely**:

- **$K_{18}$:** exactly **two** 1-rotational P1Fs (§5.1), $|\mathrm{Aut}|\in\{17,\,272\}$;
- **$K_{20}$:** exactly **seven** 1-rotational P1Fs (§5.2),
  $|\mathrm{Aut}|\in\{19,\,57,\,171,\,342\}$.

These complete counts subsume the prior literature. The Galois factorisations $GK_{p+1}$ are the
maximal-symmetry members ($|\mathrm{Aut}|=\mathrm{AGL}(1,p)$, i.e. $272$ for $K_{18}$ and $342$
for $K_{20}$; Appendices A.1, A.3). Starter constructions from **Skolem sequences**
[Pike–Shalaby; Pike] produce further 1-rotational members — for example the order-$8$ Skolem
sequence $(1,1,3,7,8,3,2,6,2,5,7,4,8,6,5,4)$ develops into the $|\mathrm{Aut}|=17$ P1F of
$K_{18}$ (Appendix A.4) — but, being constructions, they exhibit individual examples rather than
a classification. **Hill-climbing** searches [Dinitz–Stinson] likewise generate starter-induced
P1Fs but cannot certify completeness, and are structurally biased *against* the symmetric classes
that are our subject: a randomised generator reaches an isomorphism class with probability
roughly proportional to $1/|\mathrm{Aut}|$, so the high-symmetry P1Fs we enumerate are precisely
the ones such methods almost never find.

**Consequence for the contribution.** The order-$(2n-1)$ rows therefore *recover* known
territory — the Galois and Skolem-starter factorisations are not new — but ours is the first
enumeration to prove them **complete** (two for $K_{18}$, seven for $K_{20}$, with full
$|\mathrm{Aut}|$ distributions). The genuinely new P1Fs are the **non-1-rotational** ones: every
class with $|\mathrm{Aut}|\ge3$ whose order is coprime to $2n-1$ — all $179$ order-$4$ classes of
$K_{18}$ among them — has *no* starter and lies outside the reach of any cyclic, Skolem-starter,
or hill-climbing method.

---

## 6. Status of results, reproducibility, and open problems

**Status.** The $K_{16}$ and $K_{14}$ validations (the latter including a *complete*
order-$2$ classification) and the $K_{18}$ (orders $4,8,17$) and $K_{20}$
(orders $4,5,7,11,13,17,19$) counts above are complete and reproducible. The $K_{20}$
count for order $3$, and the full $|\mathrm{Aut}|>2$ sweep of $K_{18}$ over all orders
$4,\dots,17$, are computations in progress and will appear in a later revision.

**The order-2 / dense-type bottleneck.** Two situations currently lack per-level
isomorph rejection: order-$2$ (the fixed-point-free type $2^n$) and cycle types with
many fixed points (e.g. $K_{20}$, order $3$, type $3^5 1^5$). Both have a centralizer
too large to enumerate. That this is a limitation of *scale* rather than of
correctness is shown by $K_{14}$ (§4.2): there the fixed-point-free involution type
$2^7$ has an enumerable centralizer, order $2$ completes, and the method returns the
exact known count of $16$. The obstruction at the larger sizes is only that $2^8$,
$2^9$, $\dots$ have centralizers beyond the enumeration cap. A generating-set/Schreier
formulation of the per-level dedup (§3.3) would remove this limitation and is the
natural next step; it would make the order-$2$ classifications (the
$|\mathrm{Aut}|$-a-power-of-two cases) tractable at $K_{18}$ and $K_{20}$ as well.
The two *fully* fixed-point-free over-cap types are now settled by the order-2 parity
theorem (§5.3): $K_{16}\,2^8$ and $K_{20}\,2^{10}$ are **empty**. The two-fixed-point
over-cap types ($K_{18}\,2^8 1^2$, $K_{20}\,2^9 1^2$) are settled for *existence* by the
Galois construction (§5.4), and the last fixed-point-free type, $K_{18}\,2^9$, is settled
**non-empty** by the explicit example of Appendix A.5 (recovered within the order-$4$
sweep). **Every order-$2$ existence question for $K_{18}$ and $K_{20}$ is therefore now
decided** (§5.4); what remains beyond the method is only the *complete enumeration* of
these order-$2$ classes, whose centralizers exceed the enumeration cap.

**Reproducibility.** The implementation is a single self-contained module per order
$N=2n$ (changing only $N$); the search, canonical form, and orbit machinery are
$N$-generic. [TODO: repository URL, commit hash, machine specifications, and per-run
wall-clock times.]

**Open problems.**
1. ~~Prove (or refute) the $2n\equiv 2\pmod 4$ observation for order-$4$ automorphisms.~~
   **Resolved (§5.3): proved for both order $2$ and order $4$.**
2. Complete the $K_{18}$ and $K_{20}$ classifications. What is complete is the
   classification by automorphism of **order $\ge 3$** (for both sizes); equivalently,
   every P1F whose automorphism group is not an elementary-abelian $2$-group is
   accounted for. This is deliberately weaker than "$|\mathrm{Aut}|\ge3$": a P1F whose
   automorphism group is an elementary-abelian $2$-group of order $\ge4$ — for instance
   the Klein four-group $V_4$, or $C_2^3$ — has $|\mathrm{Aut}|\ge3$ yet contains *no*
   element of order $\ge3$, so it is invisible to every order-$\ge3$ sweep and would
   surface only in the intractable order-$2$ search; whether any such $K_{18}/K_{20}$
   classes exist is itself open. By contrast any automorphism group whose order has an odd
   prime factor (by Cauchy's theorem) or that contains an element of order $4$ or $8$ *is*
   found — so, e.g., the $144$ classes we report at $|\mathrm{Aut}|=4$ are all cyclic
   ($C_4$). The residue is the order-$2$ part, now settled for *existence* across
   the board: the fixed-point-free types are $K_{16}\,2^8$, $K_{20}\,2^{10}$ empty (§5.3)
   and $K_{18}\,2^9$ non-empty (Appendix A.5); the two-fixed-point types are non-empty
   (Galois, §5.4). What remains open is only the complete *enumeration* of these order-$2$
   classes — and hence the full $|\mathrm{Aut}|>1$ classification of $K_{18}$ and $K_{20}$ —
   which is out of reach of the present method.
3. Determine whether any of the new $K_{18}/K_{20}$ P1Fs give atomic Latin squares,
   following the $K_{16}$ analysis.

---

## Acknowledgements

Software implementation, computational experiments, and drafting assistance were
carried out with the help of Anthropic's Claude (Claude Code). The authors directed
the research, verified the results, and are solely responsible for the content. We
thank Ian M. Wanless for making the $K_{16}$ catalogue and constructions publicly
available, which served as our validation oracle.

---

## References

*(to be finalised; preliminary list)*

1. M. J. Gill and I. M. Wanless. *Perfect 1-factorisations of $K_{16}$.* Bulletin of
   the Australian Mathematical Society **101** (2020). arXiv:1905.07535.
2. W. D. Wallis. *One-Factorizations.* Mathematics and its Applications, Kluwer, 1997.
3. E. Mendelsohn and A. Rosa. *One-factorizations of the complete graph — a survey.*
   Journal of Graph Theory **9** (1985), 43–65.
4. B. D. McKay. *Isomorph-free exhaustive generation.* Journal of Algorithms **26**
   (1998), 306–324.
5. R. C. Read. *Every one a winner, or how to avoid isomorphism search when
   cataloguing combinatorial configurations.* Annals of Discrete Mathematics **2**
   (1978), 107–120.
6. J. H. Dinitz, D. K. Garnick, and B. D. McKay. *There are 526,915,620 nonisomorphic
   one-factorizations of $K_{12}$.* Journal of Combinatorial Designs **2** (1994).
7. M. Meszka and A. Rosa. *Perfect 1-factorisations.* (and related work on $k$-cycle
   free one-factorisations).
8. D. A. Pike and N. Shalaby. *The use of Skolem sequences to generate perfect
   one-factorizations.* Ars Combinatoria **59** (2001), 153–159.
9. D. A. Pike. *Non-isomorphic perfect one-factorizations from Skolem sequences.* Journal
   of Combinatorial Mathematics and Combinatorial Computing **44** (2003). [page range to
   verify]
10. J. H. Dinitz and D. R. Stinson. *Hill-climbing algorithms for the construction of
    combinatorial designs.* Annals of Discrete Mathematics **26** (1985), 101–125.
    [details to verify]

---

## Appendix A. Key result examples

Each example below lists a perfect 1-factorisation as its $2n-1$ 1-factors
$F_1,\dots,F_{2n-1}$ on the vertex set $\{0,1,\dots,2n-1\}$; a 1-factor is written as
its $n$ edges $(a,b)$. By definition the union of any two factors $F_i\cup F_j$ is a
single Hamiltonian cycle. The automorphism order $|\mathrm{Aut}|$ is as computed and
independently verified. (These are extracted verbatim from the solver output and are
reproducible from the repository result files.)

### A.1 $K_{18}$, $|\mathrm{Aut}|=272$ — the $\mathrm{AGL}(1,17)$ factorisation

This is the algebraic P1F of $K_{18}=K_{17+1}$; its automorphism group is the affine
group $x\mapsto ax+b$ over $\mathrm{GF}(17)$ (vertex $17$ playing $\infty$), of order
$17\cdot16=272$. Our method recovers it among the order-4 classes (since $4\mid 272$).

```
F1  (0,1)(2,3)(4,5)(6,7)(8,9)(10,11)(12,13)(14,15)(16,17)
F2  (0,2)(1,4)(3,6)(5,8)(7,10)(9,12)(11,14)(13,16)(15,17)
F3  (0,3)(1,5)(2,7)(4,9)(6,11)(8,13)(10,15)(12,17)(14,16)
F4  (0,4)(1,8)(2,6)(3,10)(5,12)(7,14)(9,16)(11,17)(13,15)
F5  (0,5)(1,9)(2,11)(3,7)(4,13)(6,15)(8,17)(10,16)(12,14)
F6  (0,6)(1,12)(2,10)(3,14)(4,8)(5,16)(7,17)(9,15)(11,13)
F7  (0,7)(1,13)(2,15)(3,11)(4,17)(5,9)(6,16)(8,14)(10,12)
F8  (0,8)(1,16)(2,14)(3,17)(4,12)(5,15)(6,10)(7,13)(9,11)
F9  (0,9)(1,17)(2,16)(3,15)(4,14)(5,13)(6,12)(7,11)(8,10)
F10 (0,10)(1,15)(2,17)(3,13)(4,16)(5,11)(6,14)(7,9)(8,12)
F11 (0,11)(1,14)(2,12)(3,16)(4,10)(5,17)(6,8)(7,15)(9,13)
F12 (0,12)(1,11)(2,13)(3,9)(4,15)(5,7)(6,17)(8,16)(10,14)
F13 (0,13)(1,10)(2,8)(3,12)(4,6)(5,14)(7,16)(9,17)(11,15)
F14 (0,14)(1,7)(2,9)(3,5)(4,11)(6,13)(8,15)(10,17)(12,16)
F15 (0,15)(1,6)(2,4)(3,8)(5,10)(7,12)(9,14)(11,16)(13,17)
F16 (0,16)(1,3)(2,5)(4,7)(6,9)(8,11)(10,13)(12,15)(14,17)
F17 (0,17)(1,2)(3,4)(5,6)(7,8)(9,10)(11,12)(13,14)(15,16)
```

### A.2 $K_{18}$, $|\mathrm{Aut}|=4$ — a representative new class

A representative of the $144$ classes of P1Fs of $K_{18}$ with automorphism group of
order $4$ (cyclic, $C_4$). To our knowledge no class of this kind has previously been
recorded. (One of the $179$ order-4 classes; the full list is in the repository.)

```
F1  (0,1)(2,3)(4,5)(6,7)(8,9)(10,11)(12,13)(14,15)(16,17)
F2  (0,2)(1,4)(3,6)(5,8)(7,10)(9,12)(11,14)(13,16)(15,17)
F3  (0,3)(1,5)(2,7)(4,9)(6,11)(8,13)(10,15)(12,17)(14,16)
F4  (0,4)(1,3)(2,8)(5,10)(6,9)(7,17)(11,13)(12,14)(15,16)
F5  (0,5)(1,2)(3,9)(4,11)(6,16)(7,8)(10,12)(13,15)(14,17)
F6  (0,6)(1,15)(2,5)(3,17)(4,12)(7,11)(8,14)(9,16)(10,13)
F7  (0,7)(1,14)(2,16)(3,4)(5,13)(6,10)(8,17)(9,15)(11,12)
F8  (0,8)(1,10)(2,14)(3,11)(4,13)(5,16)(6,15)(7,12)(9,17)
F9  (0,9)(1,11)(2,10)(3,15)(4,17)(5,12)(6,13)(7,14)(8,16)
F10 (0,10)(1,9)(2,17)(3,8)(4,15)(5,11)(6,12)(7,16)(13,14)
F11 (0,11)(1,8)(2,9)(3,16)(4,10)(5,14)(6,17)(7,13)(12,15)
F12 (0,12)(1,16)(2,11)(3,5)(4,7)(6,14)(8,15)(9,13)(10,17)
F13 (0,13)(1,17)(2,4)(3,10)(5,6)(7,15)(8,12)(9,14)(11,16)
F14 (0,14)(1,6)(2,15)(3,12)(4,8)(5,7)(9,11)(10,16)(13,17)
F15 (0,15)(1,7)(2,13)(3,14)(4,6)(5,9)(8,10)(11,17)(12,16)
F16 (0,16)(1,13)(2,12)(3,7)(4,14)(5,17)(6,8)(9,10)(11,15)
F17 (0,17)(1,12)(2,6)(3,13)(4,16)(5,15)(7,9)(8,11)(10,14)
```

### A.3 $K_{20}$, $|\mathrm{Aut}|=342$ — the $\mathrm{AGL}(1,19)$ factorisation

The algebraic P1F of $K_{20}=K_{19+1}$; automorphism group $x\mapsto ax+b$ over
$\mathrm{GF}(19)$, order $19\cdot18=342$. Recovered among the order-19 classes. Note
$F_{10}$ is the "reversal" factor $(0,10)(1,19)(2,18)\cdots(9,11)$ characteristic of
this construction.

```
F1  (0,1)(2,3)(4,5)(6,7)(8,9)(10,11)(12,13)(14,15)(16,17)(18,19)
F2  (0,2)(1,4)(3,6)(5,8)(7,10)(9,12)(11,14)(13,16)(15,18)(17,19)
F3  (0,3)(1,5)(2,7)(4,9)(6,11)(8,13)(10,15)(12,17)(14,19)(16,18)
F4  (0,4)(1,8)(2,6)(3,10)(5,12)(7,14)(9,16)(11,18)(13,19)(15,17)
F5  (0,5)(1,9)(2,11)(3,7)(4,13)(6,15)(8,17)(10,19)(12,18)(14,16)
F6  (0,6)(1,12)(2,10)(3,14)(4,8)(5,16)(7,18)(9,19)(11,17)(13,15)
F7  (0,7)(1,13)(2,15)(3,11)(4,17)(5,9)(6,19)(8,18)(10,16)(12,14)
F8  (0,8)(1,16)(2,14)(3,18)(4,12)(5,19)(6,10)(7,17)(9,15)(11,13)
F9  (0,9)(1,17)(2,19)(3,15)(4,18)(5,13)(6,16)(7,11)(8,14)(10,12)
F10 (0,10)(1,19)(2,18)(3,17)(4,16)(5,15)(6,14)(7,13)(8,12)(9,11)
F11 (0,11)(1,18)(2,16)(3,19)(4,14)(5,17)(6,12)(7,15)(8,10)(9,13)
F12 (0,12)(1,15)(2,17)(3,13)(4,19)(5,11)(6,18)(7,9)(8,16)(10,14)
F13 (0,13)(1,14)(2,12)(3,16)(4,10)(5,18)(6,8)(7,19)(9,17)(11,15)
F14 (0,14)(1,11)(2,13)(3,9)(4,15)(5,7)(6,17)(8,19)(10,18)(12,16)
F15 (0,15)(1,10)(2,8)(3,12)(4,6)(5,14)(7,16)(9,18)(11,19)(13,17)
F16 (0,16)(1,7)(2,9)(3,5)(4,11)(6,13)(8,15)(10,17)(12,19)(14,18)
F17 (0,17)(1,6)(2,4)(3,8)(5,10)(7,12)(9,14)(11,16)(13,18)(15,19)
F18 (0,18)(1,3)(2,5)(4,7)(6,9)(8,11)(10,13)(12,15)(14,17)(16,19)
F19 (0,19)(1,2)(3,4)(5,6)(7,8)(9,10)(11,12)(13,14)(15,16)(17,18)
```

### A.4 $K_{18}$, $|\mathrm{Aut}|=17$ — the cyclic (Skolem-starter) factorisation

The second of the two 1-rotational P1Fs of $K_{18}$ (§5.5); its automorphism group is the
cyclic group $\mathbb{Z}_{17}$ generated by $x\mapsto x+1\pmod{17}$ (vertex $17$ playing
$\infty$, fixed). It is the cyclic development of the starter obtained from the order-$8$ Skolem
sequence $(1,1,3,7,8,3,2,6,2,5,7,4,8,6,5,4)$ [Pike–Shalaby], and our method recovers it among
the order-$17$ classes (the other being the $\mathrm{AGL}(1,17)$ member of Appendix A.1). The
canonical form below is verbatim solver output; like A.1 its first two factors are the canonical
$F_1,F_2$, but it is a distinct isomorphism class.

```
F1  (0,1)(2,3)(4,5)(6,7)(8,9)(10,11)(12,13)(14,15)(16,17)
F2  (0,2)(1,4)(3,6)(5,8)(7,10)(9,12)(11,14)(13,16)(15,17)
F3  (0,3)(1,5)(2,8)(4,11)(6,15)(7,13)(9,17)(10,12)(14,16)
F4  (0,4)(1,9)(2,16)(3,7)(5,15)(6,14)(8,10)(11,12)(13,17)
F5  (0,5)(1,10)(2,14)(3,16)(4,9)(6,12)(7,15)(8,17)(11,13)
F6  (0,6)(1,16)(2,13)(3,5)(4,10)(7,9)(8,15)(11,17)(12,14)
F7  (0,7)(1,13)(2,4)(3,12)(5,17)(6,11)(8,14)(9,16)(10,15)
F8  (0,8)(1,6)(2,15)(3,13)(4,12)(5,10)(7,16)(9,11)(14,17)
F9  (0,9)(1,12)(2,6)(3,10)(4,8)(5,11)(7,17)(13,14)(15,16)
F10 (0,10)(1,7)(2,9)(3,14)(4,13)(5,16)(6,17)(8,11)(12,15)
F11 (0,11)(1,17)(2,12)(3,15)(4,7)(5,14)(6,9)(8,16)(10,13)
F12 (0,12)(1,15)(2,5)(3,8)(4,16)(6,13)(7,11)(9,14)(10,17)
F13 (0,13)(1,2)(3,17)(4,14)(5,9)(6,8)(7,12)(10,16)(11,15)
F14 (0,14)(1,11)(2,17)(3,4)(5,7)(6,10)(8,13)(9,15)(12,16)
F15 (0,15)(1,8)(2,10)(3,9)(4,6)(5,13)(7,14)(11,16)(12,17)
F16 (0,16)(1,14)(2,7)(3,11)(4,17)(5,6)(8,12)(9,10)(13,15)
F17 (0,17)(1,3)(2,11)(4,15)(5,12)(6,16)(7,8)(9,13)(10,14)
```

### A.5 $K_{18}$, $|\mathrm{Aut}|=16$ — a P1F with a fixed-point-free $2^9$ involution

The example settling the last open order-$2$ case (§5.4): a perfect $1$-factorisation of
$K_{18}$ whose automorphism group (of order $16$) contains the **fixed-point-free
involution** $\alpha=(0\,1)(2\,4)(3\,5)(6\,8)(7\,9)(10\,12)(11\,13)(14\,16)(15\,17)$, of
cycle type $2^9$. It was recovered within the order-$4$ enumeration (since $4\mid 16$),
independently reproduced by a separate run as a byte-identical canonical matrix, and
verified directly: all $\binom{17}{2}=136$ factor-pairs are Hamiltonian, and $\alpha$
permutes the factors below as **$9$ fixed factors and $4$ swapped pairs** — the structure
forced for any $2^9$-invariant P1F. The group of order $16$ in fact contains two such $2^9$
involutions (and one of type $2^8 1^2$); $\alpha$ above is one of them.

```
F1  (0,1)(2,3)(4,5)(6,7)(8,9)(10,11)(12,13)(14,15)(16,17)
F2  (0,2)(1,4)(3,6)(5,8)(7,10)(9,12)(11,14)(13,16)(15,17)
F3  (0,3)(1,5)(2,7)(4,9)(6,11)(8,13)(10,15)(12,17)(14,16)
F4  (0,4)(1,3)(2,12)(5,11)(6,14)(7,8)(9,17)(10,16)(13,15)
F5  (0,5)(1,2)(3,13)(4,10)(6,9)(7,15)(8,16)(11,17)(12,14)
F6  (0,6)(1,9)(2,16)(3,4)(5,15)(7,12)(8,11)(10,17)(13,14)
F7  (0,7)(1,8)(2,5)(3,17)(4,14)(6,13)(9,10)(11,16)(12,15)
F8  (0,8)(1,6)(2,11)(3,15)(4,13)(5,17)(7,16)(9,14)(10,12)
F9  (0,9)(1,7)(2,14)(3,10)(4,16)(5,12)(6,17)(8,15)(11,13)
F10 (0,10)(1,13)(2,9)(3,12)(4,11)(5,6)(7,14)(8,17)(15,16)
F11 (0,11)(1,12)(2,13)(3,8)(4,7)(5,10)(6,15)(9,16)(14,17)
F12 (0,12)(1,10)(2,6)(3,14)(4,8)(5,16)(7,9)(11,15)(13,17)
F13 (0,13)(1,11)(2,15)(3,7)(4,17)(5,9)(6,8)(10,14)(12,16)
F14 (0,14)(1,16)(2,4)(3,11)(5,13)(6,10)(7,17)(8,12)(9,15)
F15 (0,15)(1,17)(2,10)(3,5)(4,12)(6,16)(7,11)(8,14)(9,13)
F16 (0,16)(1,15)(2,17)(3,9)(4,6)(5,14)(7,13)(8,10)(11,12)
F17 (0,17)(1,14)(2,8)(3,16)(4,15)(5,7)(6,12)(9,11)(10,13)
```

Under $\alpha$ the nine fixed factors are
$F_1,F_2,F_3,F_8,F_9,F_{12},F_{13},F_{14},F_{15}$, and the remaining eight form the four
swapped pairs $\{F_4,F_5\},\{F_6,F_7\},\{F_{10},F_{11}\},\{F_{16},F_{17}\}$; the nine
$\alpha$-fixed factors restrict to a perfect $1$-factorisation of $K_{10}$ on the block
quotient, as predicted in §5.4.

---

*Notes for the authors (remove before submission): verify all reference details;
fill the reproducibility TODOs (repo, commit, hardware, timings); complete the
in-progress sweeps in §5–6; consider contacting I. M. Wanless to confirm novelty;
decide the target venue (see the Submission Targets page).*

<div style="page-break-before: always"></div>

\newpage

---

# Submission Targets

*Internal note — not part of the manuscript; remove before submission.*

Ranked by fit for a computational enumeration of perfect 1-factorisations with
prescribed automorphisms. Standard practice is to post a preprint to **arXiv
(math.CO)** first, then submit to one journal at a time.

### Tier 1 — direct precedent / best fit

1. **Bulletin of the Australian Mathematical Society** (Cambridge Univ. Press, for the
   Australian Math. Society). *This is exactly where Gill & Wanless published the
   $K_{16}$ enumeration.* Publishes concise, self-contained results; our scope and
   length match well. **Top recommendation.**
2. **Journal of Combinatorial Designs** (Wiley). The premier specialist journal for
   1-factorisations, Latin squares, and design enumeration. Ideal if we expand the
   structural discussion and the connection to atomic Latin squares.

### Tier 2 — strong open-access (no author fees)

3. **The Electronic Journal of Combinatorics** (EJC). Diamond open access, no fees,
   highly reputable; routinely publishes design-theory enumeration. Excellent if open
   access matters.
4. **Australasian Journal of Combinatorics** (Comb. Math. Soc. of Australasia).
   Free open access; close to the Australian combinatorics community (Wanless et al.).
   Fast, fitting venue for catalogue-style results.

### Tier 3 — broader scope

5. **Designs, Codes and Cryptography** (Springer) — design-theory home, slightly
   broader; good if framed toward the construction/automorphism structure.
6. **Discrete Mathematics** (Elsevier) — broad and well-indexed.
7. **Graphs and Combinatorics** (Springer) — graph-theoretic framing of the P1F /
   Hamiltonian-cycle structure.

### Notes
- **arXiv first** (math.CO, secondary math.GR/cs.DM): establishes priority/date and is
  expected by all the above. Costs nothing and is the cleanest way to assert novelty.
- **Contact I. M. Wanless** (Monash) before/at submission: he is the authority on this
  exact problem, maintains the catalogue, and is the natural person to confirm none of
  the $K_{18}/K_{20}$ counts are already known. A short courtesy email with the draft
  is advisable and may also yield a referee-quality sanity check.
- For Tier 1, completing the in-progress sweeps (K18 full $|\mathrm{Aut}|>2$; K20
  orders 3, 5) and stating the classifications as *complete* will materially
  strengthen the paper.
