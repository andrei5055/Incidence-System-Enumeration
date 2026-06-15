# Mathematical Case Studies: Cycle Lengths $\{17, 16, 15, 14, 13\}$ for $K_{18}$ P1F

This document provides the detailed mathematical reasoning, cycle structures, orbit partitioning, and symmetry analysis for the cycle lengths $\{17, 16, 15, 14, 13\}$.

---

## 1. Cycle Length 17 (Order-17 Cyclic Automorphisms)

### A. Cycle Structure
* **Structure**: $(17)(1)$ — One 17-cycle and one 1-cycle (1 fixed point). This is the only possible cycle structure of order 17 on 18 vertices.

### B. Mathematical Families
Exhaustive search and canonization confirmed exactly **two non-isomorphic cyclic P1Fs** for this order:
1. **GK (Georges-Kotzig) Construction (Automorphism group size: 272)**:
   * The automorphism group is isomorphic to the 1D affine group $GA(1, 17) \cong \mathbb{Z}_{17} \rtimes \mathbb{Z}_{16}$ over the finite field $\mathbb{F}_{17}$.
   * It consists of all linear transformations of the form $x \mapsto ax + b \pmod{17}$ with $\infty$ fixed.
2. **GB (Bryant-Maenhaut-Wanless) Construction (Automorphism group size: 17)**:
   * Discovered in the 2000s, this is another infinite family of cyclic P1Fs of complete graphs $K_{p+1}$ (for prime $p$).
   * Its automorphism group is just the cyclic shift group $\mathbb{Z}_{17}$ of order 17.

---

## 2. Cycle Length 16 (Order-16 Automorphisms)

### A. Parity and Orbit Partitioning
Any orbit of factors under an order-16 automorphism must have a size dividing 16 (i.e. $\{1, 2, 4, 8, 16\}$). Since the total number of factors in $K_{18}$ is 17:
* The only possible partition of factor orbit sizes is:
  $$17 = 16 + 1$$
* This guarantees that any valid P1F must consist of exactly **one cyclic orbit of 16 factors** ($G[0] \dots G[15]$) and **one fixed factor** ($F_f$) invariant under the automorphism.

### B. Matching of the Fixed Points
Let the two fixed points of the automorphism be $0$ and $f$. 
* In the fixed factor $F_f$, **$0$ and $f$ must be matched to each other** (i.e. $\{0, f\} \in F_f$).
* *Proof*: If $0$ were matched to some vertex $u$ in the 16-cycle, then invariance of $F_f$ would require $0$ to be matched to $\alpha(u)$. Since $F_f$ is a matching, this requires $\alpha(u) = u$, which is a contradiction because all vertices in the 16-cycle have orbits of size 16.

### C. Completeness of Compatibility Checks
The solver verifies the P1F using:
1. **Internal Compatibility**: Checks $G[0] \cup G[j]$ is a single 18-cycle for $j \in \{1 \dots 8\}$. Since $G[a] \cup G[b] \cong G[0] \cup G[b-a \pmod{16}]$, checking the first 8 differences covers all pairs in the cycle.
2. **Fixed Factor Compatibility**: Checks $G[0] \cup F_f$ is perfect. Since $F_f$ is invariant under $\alpha$:
   $$G[k] \cup F_f = \alpha^k(G[0]) \cup \alpha^k(F_f) = \alpha^k(G[0] \cup F_f) \cong G[0] \cup F_f$$
   This ensures that $F_f$ forms a perfect union cycle with all 16 factors.

### D. Canonical Representative Count (5 Matrices)
The 80 passing permutations canonize to exactly 5 matrices:
* **1 matrix of size 272**: This is the GK construction (which naturally contains a 16-cycle subgroup).
* **4 matrices of size 16**: These represent exactly **two unique isomorphism classes**, with each class split into a "forward" and "reverse" (inverse) representation of the cycle.

---

## 3. Cycle Length 15 (Order-15 Cyclic Automorphisms)

### A. Cycle Structures
* **Type A**: $(15)(1)(1)(1)$ (one 15-cycle, three 1-cycles).
* **Type B**: $(15)(3)$ (one 15-cycle, one 3-cycle).

### B. Impossibility Proof
1. **No Fixed Factors**: For any odd order automorphism $\alpha$, any edge in a fixed factor must connect two fixed points of $\alpha$ (since odd order groups contain no 2-cycles). Since $18 - L \le 3 < 18$, the fixed points cannot cover all 18 vertices, so **no fixed factors can exist ($a = 0$)**.
2. **Orbit Partition**: The 17 factors must be partitioned into orbits of sizes dividing 15 (divisors $\{3, 5, 15\}$). The only partition of 17 is:
   $$17 = 5 + 3 \times 4$$
   This requires the existence of orbits of size 3 and size 5.
3. **Parity Contradiction**:
   * Any factor in an orbit of size 3 is invariant under $\alpha^3$ (order 5, odd).
   * Any factor in an orbit of size 5 is invariant under $\alpha^5$ (order 3, odd).
   * Since $\alpha^3$ and $\alpha^5$ have odd order, they cannot have fixed factors. This makes orbits of size 3 and 5 mathematically impossible. Thus, no such P1F can exist.

---

## 4. Cycle Length 14 (Order-14 Automorphisms)

### A. Cycle Structures
* **Type A**: $(14)(1)(1)(1)(1)$ (4 fixed points) — *Exhaustively tested (0 results).*
* **Type B**: $(14)(2)(1)(1)$ (2 fixed points) — *Unexplored.*
* **Type C**: $(14)(2)(2)$ (0 fixed points) — *Unexplored.*

### B. Feasibility Analysis
Unlike $L=15$, cycle length 14 is **mathematically possible** because we can write 17 as a partition of divisors of 14 ($\{1, 2, 7, 14\}$):
* For example, $17 = 14 + 2 + 1$ (one orbit of size 14, one of size 2, and one fixed factor). In the fixed factor, the 4 fixed points of Type A can be matched among themselves in 2 edges, which is valid since 4 is even.
* The search on the 405,405 orbit representatives of Type A yielded **exactly 0 results**, proving computationally that no P1Fs exist for Type A. Types B and C remain open.

---

## 5. Cycle Length 13 (Order-13 Cyclic Automorphisms)

### A. Impossibility Proof
1. **No Fixed Factors**: The number of fixed points is $18 - 13 = 5$. Since 5 is odd, the fixed points cannot be matched among themselves, so **no fixed factors can exist ($a = 0$)**.
2. **Orbit Partition**: The 17 factors must be partitioned into orbits whose sizes divide 13. The only divisor of 13 (excluding 1) is 13.
3. **Partition Contradiction**: It is impossible to write 17 as a sum of multiples of 13:
   $$17 \neq 13 \times k$$
   Thus, no cyclic P1F of cycle length 13 can ever exist.
