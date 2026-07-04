# Automorphism Search Subsets & Classification Options for $K_{18}$ P1F

This document outlines the recommended subsets of cycle lengths to target for a publishable classification result of parallel 1-factorizations (P1F) of the complete graph $K_{18}$.

---

## 1. The Core Bottleneck: Reconstruction of Missing Matchings

Although smaller cycle lengths ($L$) have very few starting permutations (orbit representatives) to check, they are computationally much slower.

The computation time is dominated by the function `decomposeMissingEdges`, which must find $17 - L$ pairwise disjoint perfect matchings in the complement graph:

* **$L = 17$**: $17 - 17 = 0$ missing matchings. The search is extremely fast per representative (simple compatibility check).
* **$L = 16$**: $17 - 16 = 1$ missing matching. Very fast.
* **$L = 14$**: $17 - 14 = 3$ missing matchings. Highly feasible.
* **$L = 12$**: $17 - 12 = 5$ missing matchings. Feasible and runs quickly.
* **$L = 10$**: $17 - 10 = 7$ missing matchings. Feasible with the Exact Cover/MRV optimization in Tt4.
* **$L \le 8$**: $\ge 9$ missing matchings. The complement graph has millions of possible perfect matchings, and finding a clique of $\ge 9$ mutually disjoint matchings becomes a massive exact cover problem, causing timeouts.

---

## 2. Recommended Classification Subsets for Publication

To obtain a mathematically clean and publishable result, you can restrict the search space to a specific cyclic automorphism order. Below are the three recommended options:

### Option A: Cyclic Automorphisms of Order $L \ge 12$ (Highly Recommended)
* **Goal**: Classify all parallel 1-factorizations of $K_{18}$ admitting an automorphism of order $\ge 12$.
* **Cycle Lengths to Run**: $L \in \{12, 14, 16, 17\}$.
* **Feasibility**: Very High. The maximum number of missing matchings is 5. The exact cover solver will complete quickly.
* **Result Presentation**: *"Classification of parallel 1-factorizations of $K_{18}$ admitting a cyclic automorphism of order at least 12."*

### Option B: Cyclic Automorphisms of Order $L \ge 10$
* **Goal**: Classify all parallel 1-factorizations of $K_{18}$ admitting an automorphism of order $\ge 10$.
* **Cycle Lengths to Run**: $L \in \{10, 12, 14, 16, 17\}$.
* **Feasibility**: High. $L=10$ requires finding 7 missing matchings. With the MRV exact cover optimization in Tt4, this is expected to complete without timing out.
* **Result Presentation**: *"Classification of parallel 1-factorizations of $K_{18}$ admitting a cyclic automorphism of order at least 10."*

### Option C: Large Cyclic Automorphisms ($L \ge 6$) with a Timeout
* **Goal**: Classify parallel 1-factorizations admitting cyclic automorphisms of order $\ge 6$, noting that some small cycle lengths were truncated due to search space complexity.
* **Cycle Lengths to Run**: $L \in \{6, 8, 10, 12, 14, 16, 17\}$.
* **Feasibility**: Moderate. Requires setting a timeout (e.g., 30 minutes) for the slower cases ($L=8$ and $L=6$).
* **Result Presentation**: A complete classification for $L \ge 10$, combined with partial/truncated search statistics for $L = 6, 8$.
