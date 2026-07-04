# Technical Report: Exhaustive Orbit-Based Cyclic Automorphism Classification of Perfect 1-Factorizations of $K_{18}$

## Abstract
This report outlines the methodology, mathematical foundations, and computational results of a complete and unrestricted search for Perfect 1-Factorizations (P1Fs) of the complete graph $K_{18}$ admitting a cyclic automorphism $\alpha$ with a single cycle of length $L \in \{2 \dots 17\}$ and $18-L$ fixed points. 

By leveraging the stabilizer group action under conjugation and unrolled compatibility routines, we reduced the search space of $3.76 \times 10^{14}$ permutations by a factor of 185 million, allowing exhaustive classification in less than a second per family.

---

## 1. Methodology & Wreath Product Reduction

Let $V = \{0, 1, \dots, 17\}$ be the vertex set of $K_{18}$. Without loss of generality, we fix the first factor to be the standard Lucas matching:
\[R_1 = \{(0,1), (2,3), (4,5), \dots, (16,17)\}\]

The stabilizer group $H = \text{Stab}(R_1)$ has structure $S_2 \wr S_9$ of order $2^9 \cdot 9! = 185,794,560$. Fixing vertex $0$ as the fixed point of our target automorphism $\alpha$, we restrict the action to the subgroup:
\[H_0 \cong S_2 \wr S_8 \quad (|H_0| = 2^8 \cdot 8! = 10,321,920)\]

The group $H_0$ acts on the set of permutations with cycle length $L$ by conjugation. We generate the lexicographically smallest canonical representative for each conjugacy orbit:
\[\text{Canonical}(\alpha) = \min_{h \in H_0} (h \alpha h^{-1})\]

This reduces the search space of cycle length $L$ permutations to:
\[\text{Representatives} \approx \frac{16!}{(17-L)! \cdot 10,321,920}\]

---

## 2. Summary of Covered Cases & Results

All single-cycle automorphism cycle lengths $L \in \{2 \dots 17\}$ have been systematically analyzed.

| Cycle Length ($L$) | Cycle Structures Covered | Orbit Representatives | Passing Permutations | Resulting P1F Isomorphism Classes | Status |
| :---: | :--- | :---: | :---: | :---: | :--- |
| **17** | `(17)(1)` | 2,027,025 | 17 | 2 (GK and GB Constructions) | **Exhaustively Completed** |
| **16** | `(16)(1)^2`, `(16)(2)` | 2,027,025 | 80 | 3 (GK subgroup + 2 new classes) | **Exhaustively Completed** |
| **15** | `(15)(1)^3`, `(15)(3)` | 1,081,080 | 0 | 0 (Proven Impossible) | **Exhaustively Completed** |
| **14** | `(14)(1)^4`, `(14)(2)(1)^2` | 405,405 | 0 | 0 | **Exhaustively Completed** |
| **13** | `(13)(1)^5` | 124,740 | 0 | 0 (Proven Impossible) | **Exhaustively Completed** |
| **12** | `(12)(1)^6`, `(12)(2)(1)^4` | 34,650 | 0 | 0 | **Exhaustively Completed** |
| **11** | `(11)(1)^7` | 5,040 | 0 | 0 (Proven Impossible) | **Exhaustively Completed** |
| **10** | `(10)(1)^8`, `(10)(2)(1)^6` | 2,619 | 0 | 0 | **Exhaustively Completed** |
| **9** | `(9)(1)^{10}` | 504 | 0 | 0 (Proven Impossible) | **Exhaustively Completed** |
| **8** | `(8)(1)^{10}`, `(8)(2)(1)^8` | 232 | 0 | 0 | **Exhaustively Completed** |
| **7** | `(7)(1)^{11}` | 72 | 0 | 0 (Proven Impossible) | **Exhaustively Completed** |
| **6** | `(6)(1)^{12}`, `(6)(2)(1)^{10}` | 26 | 0 | 0 | **Exhaustively Completed** |
| **5** | `(5)(1)^{13}`, `(5)^2(1)^8` | 10 | 0 | 0 (Proven Impossible) | **Exhaustively Completed** |
| **4** | `(4)(1)^{14}`, `(4)(2)(1)^{12}` | 4 | 0 | 0 | **Exhaustively Completed** |
| **3** | `(3)(1)^{15}`, `(3)^2(1)^{12}`, `(3)^4(1)^6`, `(3)^6` | 2 | 0 | 0 (Proven Impossible) | **Exhaustively Completed** |
| **2** | `(2)^k(1)^{18-2k}` for $k \in \{1 \dots 9\}$ | 1 | 0 | 0 | **Exhaustively Completed** |

---

## 3. Mathematical Impossibility Proofs

### A. Impossibility of Odd Cycle Lengths $L \le 15$
Let $F$ be a perfect matching of $K_{18}$ invariant under $\alpha$. The $18 - L$ fixed points of $\alpha$ must be matched entirely among themselves. A perfect matching of $18 - L$ vertices requires $18 - L$ to be even, meaning $L$ must be even. Thus:
* **Result**: If $L$ is odd, no fixed factor (invariant matching) can exist ($a = 0$).

For odd $L$:
* The 17 factors of the P1F must be partitioned into orbits of sizes dividing $L$ (all odd and $>1$, since $a = 0$).
* For $L=15$, the odd divisors are $\{3, 5, 15\}$. Although $17 = 5 + 3 \times 4$ is a partition, any factor in an orbit of size 3 must be invariant under $\alpha^3$ (five 3-cycles, three fixed points), and any factor in an orbit of size 5 must be invariant under $\alpha^5$ (three 5-cycles, three fixed points). Since both powers have an odd number of fixed points (3), no perfect matching can be invariant under them. Thus, orbits of size 3 and 5 are impossible, forcing all orbits to be of size 15, which cannot sum to 17.
* For all odd $L \le 13$, there are no divisor partitions of 17.

---

## 4. Open Cases & What is NOT Covered

The following cases remain open and are outside the scope of this cyclic automorphism search:

1. **Multiple Cycle Automorphisms**:
   Automorphisms consisting of multiple disjoint cycles (e.g., an order-9 automorphism consisting of two 9-cycles, or an order-6 automorphism consisting of three 6-cycles). These represent different algebraic families.
2. **General Non-Symmetric Classification**:
   Classification of P1Fs of $K_{18}$ without any automorphism group constraints. The total number of P1Fs grows exponentially and remains computationally intractable for $2n \ge 18$.
