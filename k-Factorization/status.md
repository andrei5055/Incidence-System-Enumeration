# Project Status & Architecture: Complete Graphs Solver (Tr71)

This document describes the current architecture, mathematical models, algorithm implementations, and optimizations of the complete graphs solver system.

---

## 1. Component Architecture & Configuration

The project implements two main active solver classes designed to search for Perfect 1-Factorizations (P1F) of complete graphs $K_{2n}$:

### A. K16A2 Solver
* **Definition**: Defined in [k16A2.h](file:///C:/TripleSys/Tr71/TripleSys/Include/k16A2.h) and implemented in [k16A2.cpp](file:///C:/TripleSys/Tr71/TripleSys/Source/k16A2.cpp).
* **Target**: Complete graph $K_{16}$ ($16$ vertices, $15$ factors).
* **Search Method**: Exact cover backtracking search (Knuth's Algorithm X) with static MRV ordering, look-ahead pruning, and SSE/AVX registers.
* **Symmetries**: Checks for prime-order player automorphisms $p \in \{2, 3, 5, 7\}$.

### B. K18A2 Solver
* **Definition**: Defined in [k18a2.h](file:///C:/TripleSys/Tr71/TripleSys/Include/k18a2.h) and implemented in [k18a2.cpp](file:///C:/TripleSys/Tr71/TripleSys/Source/k18a2.cpp).
* **Target**: Complete graph $K_{18}$ ($18$ vertices, $17$ factors).
* **Search Method**: Order-17 cyclic automorphism symmetry validation. For a starter (3 fixed rows), the solver precalculates valid transition chains in `init()`. It then filters candidates incrementally as they are added via `addRow()`, and reports solutions directly in `solve()`.
* **Symmetries**: Prime-order player automorphisms $p = 17$.

---

## 2. Core Search Optimizations for K18A2

The $K_{18}$ solver utilizes multiple low-level and algorithmic optimizations to accelerate validation:

### A. Unrolled Cycle Union Validation (`is_perfect_scalar`)
* The cycle validation routine `is_perfect_scalar` is unrolled into exactly 9 sequential lookup pairs and moved into [k18a2.h](file:///C:/TripleSys/Tr71/TripleSys/Include/k18a2.h) as a `FORCE_INLINE` function, completely eliminating loop overhead.

### B. Unrolled Permutation Application (`apply_perm_18`)
* Permutation application `apply_perm_18` is unrolled into 18 direct assignments and moved into [k18a2.h](file:///C:/TripleSys/Tr71/TripleSys/Include/k18a2.h) as a `FORCE_INLINE` function.

### C. Allocation-Free Lookups (`f_map_unordered`)
* The look-up table `f_map_unordered` uses `std::array<uint8_t, 18>` as its key with a custom hash function, avoiding dynamic heap allocation during candidate factor validation.

### D. Precalculated Chain Generation (in `init`)
* Candidate automorphism mappings of order 17 are pre-generated on initialization for three transition configurations:
  * $(R_1, R_2) \to (R_2, R_3)$
  * $(R_1, R_3) \to (R_3, R_2)$
  * $(R_2, R_3) \to (R_3, R_1)$
* The entire factor chain $R_4 \dots R_{17}$ is generated and checked for vertex 0 neighbor uniqueness and pairwise compatibility (`is_perfect_packed`) on initialization, keeping only valid chains.

### E. Incremental Candidate Pruning (in `addRow`)
* Incoming candidates for row `rowNum` are compared against expected factors for all active precalculated chains.
* If a candidate does not match any active chain, it is skipped (not inserted into pool), reducing database sizes.
* If a new row index starts but no chains remain active (the previous row had zero matching candidates), the solver returns `false` early to abort the current starter.

### F. Compact Rejection Statistics
* The solver tracks candidates rejected at each step in `rejected_candidates[18]`.
* Compact single-line rejection statistics are output on early return in `addRow` or upon completion in `solve()` (e.g., `Rejections: 4:320 5:0 ...`), avoiding multi-line table printing.

---

## 3. Build & Execution Status

* **Build Method**: MSBuild via [build.bat](file:///C:/TripleSys/Tr71/build.bat) compiles the solution under `Release x64` configuration.
* **Target Executable**: `k-Sys.exe` is generated in `x64\Release` and executed.
* **Run Scripts**: Solver parameters are defined in files like [param18a2.txt](file:///C:/TripleSys/Tr71/NewTests/param18a2.txt) and executed using batch scripts like [run18a2.bat](file:///C:/TripleSys/Tr71/NewTests/run18a2.bat) from the [NewTests](file:///C:/TripleSys/Tr71/NewTests/) folder.

---

## 4. Classification Coverage & Limitations (Order-16 Automorphisms)

* **Exhaustive Coverage (Type A & Type B)**: The solver is 100% complete for all cyclic automorphisms of order 16 on $K_{18}$. This includes:
  - **Type A** (cycle structure $(16)(1)^2$, with two fixed points): Searched using Search Type 1.
  - **Type B** (cycle structure $(16)(2)$, with the two non-cycle vertices transposed): Searched using Search Type 2.
* **Resulting Classes**: By searching both structures, all valid cyclic automorphisms of order 16 are exhaustively covered. The search identified the GK construction (group size 272) and two new isomorphism classes of group size 16.

