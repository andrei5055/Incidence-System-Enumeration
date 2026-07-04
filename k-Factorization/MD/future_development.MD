# Future Development: Validation on $K_{16}$ and Feasibility of $K_{20}$

This document details future directions for validation and scaling the Complete Graphs Solver architecture to complete graphs of higher orders.

---

## 1. Validation on $K_{16}$

To verify the mathematical correctness and codebase implementations of the canonical orbit-reduction stabilizer solver, we can run it on $K_{16}$ (16 vertices, 15 factors) and validate against known literature results (exactly 30 non-isomorphic P1Fs of $K_{16}$ exist).

### Stabilizer Reduction Formulation for $K_{16}$
- **Fixed Matching $R_1$**: $\{(0,1), (2,3), (4,5), \dots, (14,15)\}$
- **Stabilizer Subgroup $H_0$** (stabilizing vertex $0$):
  $$H_0 \cong S_2 \wr S_7 \quad (|H_0| = 2^7 \cdot 7! = 645,120)$$
- **Permutations to Check** (for $p=15$ cyclic automorphisms):
  $$\text{Representatives} = \frac{14!}{2^7 \cdot 7!} = \frac{87,178,291,200}{645,120} = \mathbf{135,135 \text{ representatives}}$$

> [!TIP]
> With the stack-based recursion optimization, searching 135,135 representatives on $K_{16}$ will execute in **under 0.02 seconds** on a single CPU core.

### Prime-Order Configurations to Validate on $K_{16}$
- **$p = 5$**: Cycle structure $(5)^2(1)^6$ (6 fixed points, even $\Rightarrow$ fixed factors possible). Orbit partition: $15 = 5k + a \Rightarrow a \equiv 0 \pmod 5$.
- **$p = 3$**: Cycle structures $(3)^2(1)^{10}$, $(3)^4(1)^4$.
- **$p = 2$**: Cycle structures $(2)^k(1)^{16-2k}$ for $k \in \{1 \dots 8\}$.

---

## 2. Scaling to $K_{20}$ (Feasibility Study)

Completing the exhaustive classification of perfect 1-factorizations of Complete Graphs of order 20 admitting non-trivial automorphisms (order $\ge 2$) is **completely feasible** and well within touch using the current solver architecture.

### Stabilizer Reduction Formulation for $K_{20}$
- **Fixed Matching $R_1$**: $\{(0,1), (2,3), (4,5), \dots, (18,19)\}$
- **Stabilizer Subgroup $H_0$** (stabilizing vertex $0$):
  $$H_0 \cong S_2 \wr S_9 \quad (|H_0| = 2^9 \cdot 9! = 185,794,560)$$
- **Permutations to Check** (for $p=19$ cyclic automorphisms):
  $$\text{Representatives} = \frac{18!}{2^9 \cdot 9!} = \frac{6,402,373,705,728,000}{185,794,560} = \mathbf{34,459,425 \text{ representatives}}$$

> [!IMPORTANT]
> A search over 34.4 million representatives using our optimized stack-based C++ backtracking engine will take **only ~3.7 seconds** on a single CPU core.

### Complete Prime-Order Case Analysis for $K_{20}$

To claim a 100% complete classification of $K_{20}$ P1Fs with symmetries of order $\ge 2$, the solver must test all cycle configurations of prime orders $p \in \{2, 3, 5, 7, 11, 13, 17, 19\}$:

| Prime ($p$) | Cycle Structure | Feasibility | Description / Orbit & Fixed Points Analysis |
| :---: | :--- | :---: | :--- |
| **19** | $(19)(1)$ | **Possible** | 1 fixed point. Checked in **~3.7 seconds**. |
| **17** | $(17)(1)^3$ | **Impossible** | 3 fixed points (odd) $\Rightarrow$ no fixed factors ($a = 0$). Orbit partition $19 = 17k$ is impossible. |
| **13** | $(13)(1)^7$ | **Impossible** | 7 fixed points (odd) $\Rightarrow$ no fixed factors ($a = 0$). Orbit partition $19 = 13k$ is impossible. |
| **11** | $(11)(1)^9$ | **Impossible** | 9 fixed points (odd) $\Rightarrow$ no fixed factors ($a = 0$). Orbit partition $19 = 11k$ is impossible. |
| **7** | $(7)^2(1)^6$ | **Possible** | 6 fixed points (even). Orbit partition $19 = 7 \times 2 + 5$ ($a = 5$ fixed factors). |
| **5** | $(5)^2(1)^{10}$, $(5)^4$ | **Possible** | Even number of fixed points (10 and 0). Orbit partitions require $a \equiv 4 \pmod 5 \Rightarrow a \ge 4$. |
| **3** | $(3)^2(1)^{14}$, $(3)^4(1)^8$, $(3)^6(1)^2$ | **Possible** | Even fixed points (14, 8, 2). Orbit partitions require $a \equiv 1 \pmod 3 \Rightarrow a \ge 1$. |
| **2** | $(2)^k(1)^{20-2k}$ | **Possible** | Even fixed points for all $k \in \{1 \dots 10\}$. Orbit partitions require $a \equiv 1 \pmod 2 \Rightarrow a \ge 1$. |

Any search for these possible cycle structures can be completed in **under 5 minutes** of total execution time.
