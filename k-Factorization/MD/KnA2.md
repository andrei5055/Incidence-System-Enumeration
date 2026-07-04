# Specifications & Concepts: Complete Graph Solvers (K16 & K18 Implementation)

This document specifies the mathematical algorithms, combinatorics, search-space optimizations, and actual implementation details for the complete graph solvers ($K_{16}$ exact cover solver and $K_{18}$ symmetry solver).

---

## 1. Project Goal & Target Rationale

### The $K_{16}$ Verification Solver (`K16A2`)
The [K16A2](file:///C:/TripleSys/Tr71/TripleSys/Source/k16A2.cpp) class (defined in [k16A2.h](file:///C:/TripleSys/Tr71/TripleSys/Include/k16A2.h)) is configured and compiled as a highly optimized **$K_{16}$ verification solver**. All parameters, cycle lengths, and SSE/AVX registers in `K16A2` are tuned to $N=16$ (15 factors, 12 search slots).

### The $K_{18}$ Solver (`K18A2`)
The [K18A2](file:///C:/TripleSys/Tr71/TripleSys/Source/k18a2.cpp) class (defined in [k18a2.h](file:///C:/TripleSys/Tr71/TripleSys/Include/k18a2.h)) implements an order-17 cyclic automorphism symmetry solver. Rather than running a full exact cover backtracking tree, it checks if starter matrices can be completed to a Perfect 1-Factorization (P1F) invariant under a player automorphism $\alpha$ of prime order $p = 17$.

---

## 2. Mathematical Symmetries

A perfect 1-factorization consists of exactly $2n-1$ rows (1-factors). Under a player automorphism $\alpha$ of prime order $p$:
* **$K_{16}$ Symmetries ($N_{factors}=15$)**: Partitioned into rings under prime orders $p \in \{2, 3, 5, 7\}$.
* **$K_{18}$ Symmetries ($N_{factors}=17$)**: Verified under cyclic automorphism of prime order $p = 17$, where all 17 factors belong to a single orbit (ring) of size 17 under $\alpha$ (mapping $R_1 \to R_2 \to R_3 \to \dots \to R_{17} \to R_1$).

---

## 3. The Pair Transition Optimization

Specifying a single row transition (e.g., $R_1 \to R_2$) requires checking permutations, which is computationally expensive.

By specifying a **pair transition** (such as $(R_1, R_2) \to (R_2, R_3)$, $(R_1, R_3) \to (R_3, R_2)$, or $(R_2, R_3) \to (R_3, R_1)$):
1. Because the factorization is perfect, the unions of the pairs (e.g. $R_1 \cup R_2$ and $R_2 \cup R_3$) form Hamiltonian cycles of length $2n$.
2. An automorphism $\alpha$ mapping one union cycle to the other must map one cycle to the other while preserving edge adjacency.
3. This restricts the player permutations $\alpha$ to **exactly $2n = 36$** candidates per transition. This accelerates the permutation search step by over **90 million times**.

---

## 4. K18A2 Solver Algorithm Workflow

Instead of backtracking, the $K_{18}$ solver verifies the order-17 symmetry cyclically:

### Phase 1: Transition Permutation Generation
Given a starter $(R_1, R_2, R_3)$, find the candidate permutations $\alpha$ for each of the three active pair transition configurations:
* $(R_1, R_2) \to (R_2, R_3)$
* $(R_1, R_3) \to (R_3, R_2)$
* $(R_2, R_3) \to (R_3, R_1)$

### Phase 2: Linear Chain Verification
For each candidate permutation $\alpha$ of order 17, we generate the factor chain $R_k = \alpha(R_{k-1})$ for $k = 4 \dots 17$:
1. **Uniqueness of Vertex 0**: Verify that the neighbor of vertex 0 in the new factor ($R_k[0]$) has not been used by any previous factor in the chain. If it is used, prune the transition immediately.
2. **Pool Lookup**: Perform an $O(1)$ lookup using a custom `std::array<uint8_t, 18>` key in `f_map_unordered` to verify if the factor exists in the candidate pool. If not, prune immediately.
3. **Pairwise Compatibility**: Verify that the generated factor is pairwise compatible and perfect cycle-wise with all previously generated factors in the chain using the inlined `is_perfect_packed()`. If not, prune immediately.

### Phase 3: Validation and Callback
If all 17 factors are successfully generated, the solution is sorted into canonical order and reported via `resultCallback`.
