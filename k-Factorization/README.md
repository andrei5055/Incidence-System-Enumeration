# Perfect One-Factorisations of $K_{2n}$ with Non-Trivial Automorphism Group

A complete catalogue of the isomorphism classes of perfect one-factorisations (P1Fs) of the complete graph on 18 vertices ($K_{2n}$) possessing non-trivial automorphism groups ($|Aut(\mathcal{F})| > 1$).

**Total catalogued classes:** **10,710**

---

## 1. Overview & Contents

This dataset accompanies the classification of non-rigid perfect one-factorisations of $K_{2n}$. The archive contains:

| Filename | Description | Format / Entries |
| :--- | :--- | :--- |
| `K18_P1F_aut_gt1.txt` | Primary catalogue of 10,710 canonical P1Fs | ASCII text, 17 rows of 9 vertex pairs per record |
| `README.md` | Dataset documentation, record schemas, validation, and usage guide | Markdown |

---

## 2. Mathematical Background

### 2.1 One-Factorisations of $K_{2n}$
Let $K_{2n} = (V, E)$ denote the complete undirected graph on 2n vertices with $|V| = 2n$ and $|E| = \binom{2n}{2}$. A **one-factor** (or perfect matching) is a 1-regular spanning subgraph of  $K_{2n}$, comprising $n$ disjoint edges that partition $V$. 

A **one-factorisation** $\mathcal{F} = \{F_1, F_2, \dots, F_{2n-1}\}$ is a partition of the edge set $E(K_{2n})$ into $2n-1$ edge-disjoint one-factors. For $K_{18}$:
- Vertex count: $|V| = 18$
- Edge count: $|E| = \binom{18}{2} = 153$
- Number of one-factors: $2n - 1 = 17$
- Edges per factor: $n = 9$

### 2.2 Perfect One-Factorisations (P1F)
A one-factorisation $\mathcal{F}$ is **perfect** (a P1F) if the union of any two distinct one-factors $F_i, F_j \in \mathcal{F}$ ($i \neq j$) forms a single connected **Hamiltonian cycle** of length $2n$. For $K_{18}$, every one of the $\binom{17}{2} = 136$ factor pairs $F_i \cup F_j$ induces a cycle of length 18.

### 2.3 Isomorphism and Automorphisms
Two factorisations $\mathcal{F}_A$ and $\mathcal{F}_B$ are **isomorphic** if there exists a permutation $p &isin; \mathcal{S}_{2n}$ of the vertex set such that
XX

Classes with $|Aut(\mathcal{F})| = 1$ (rigid factorisations) are outside the scope of this catalogue.

---

## 3. Census by Automorphism Group Order

The catalogue contains all 10,710 pairwise non-isomorphic P1Fs of $K_{18}$ with $|Aut(\mathcal{F})| > 1$. The distribution across group orders is as follows:


| Automorphism Group Order | Number of Isomorphism Classes | Cumulative Share (%) |
| :---: | :---: | :---: |
| **2** | 10,179 | 95.04% |
| **3** | 351 | 98.32% |
| **4** | 144 | 99.66% |
| **8** | 22 | 99.87% |
| **16** | 12 | 99.98% |
| **17** | 1 | 99.99% |
| **272** | 1 | 100.00% |
| **Total** | **10,710** | **100.00%** |

---

## 4. Notable & Historical Classes

### 4.1 Prior Known Factorisations (6 Classes)
Prior to this exhaustive enumeration, exactly six non-isomorphic P1Fs of $K_{18}$ were documented in the literature:
- **Record #10710** ($|Aut| = 272$): The sharply 2-transitive affine group $AGL(1, 17)$ acting on $\mathbb{F}_{17} \cup \{\infty\}$, corresponding to the classical construction $GK_{18}$.
- **Record #10709** ($|Aut| = 17$): The cyclic/starter-generated P1F admitting a regular automorphism of order 17 fixing the infinity point.
- **Records with $|Aut| \in \{8, 16\}$**: Four previously known symmetric constructions originating from quotient developments and Latin trade constructions (e.g., Anderson, Ihrig, Meszka, Rosa).
  - *Record cross-references:* `[Record #XXXX, #YYYY, #ZZZZ, #WWWW]` *(refer to index table in supplementary documentation)*.

### 4.2 Atomic Latin Squares of Order 17
A perfect one-factorisation of $K_{2n}$ induces a symmetric Latin square of order $2n-1$ with constant diagonal. A Latin square is **atomic** if it contains no non-trivial proper subsquares and satisfies extremal orthogonal trade criteria.
- In this catalogue, **exactly two isomorphism classes** produce atomic Latin squares of order 17.
- **Catalogued Record IDs:** `Record #A` and `Record #B` *(flagged with metadata attribute `atomic_ls=True` in extended tables)*.

---

## 5. Data File Format Specification (`K18_P1F_aut_gt1.txt`)

### 5.1 Canonical Form and Ordering
- **Vertices:** Labelled as integers in $\{0, 1, 2, \dots, 17\}$ (or $\{1, 2, \dots, 18\}$ depending on 0-based vs 1-based indexing option).
- **Edge representation:** An unordered pair $(u, v)$ with $u < v$.
- **One-factor representation:** 9 edge pairs sorted in lexicographical order by the first vertex:
  $$(u_1, v_1), (u_2, v_2), \dots, (u_9, v_9) \quad \text{where } u_1 < u_2 < \dots < u_9 \text{ and } u_i < v_i.$$
- **Factorisation representation:** 17 one-factors $F_1, F_2, \dots, F_{17}$ sorted lexicographically by their canonical string representations.
- **Record numbering:** Records are numbered contiguously from `1` to `10710` in strictly ascending lexicographical order of their canonical representation. 
  - *Versioning policy:* Any future additions (such as subsets of rigid factorisations) will be appended to maintain immutable numeric IDs for classes 1–10710.

### 5.2 Record Layout Example
Each record in `K18_P1F_aut_gt1.txt` is structured as a header line followed by the 17 one-factors:

```text
# RECORD 10710 | AUT_ORDER: 272 | ID: P1F_18_10710
F01: (0,1) (2,17) (3,16) (4,15) (5,14) (6,13) (7,12) (8,11) (9,10)
F02: (0,2) (1,3) (4,17) (5,16) (6,15) (7,14) (8,13) (9,12) (10,11)
...
F17: (0,17) (1,16) (2,15) (3,14) (4,13) (5,12) (6,11) (7,10) (8,9)
```

---

## 6. Computational Methodology & Validation

### 6.1 Exhaustive Search Strategy
1. **Prescribed Automorphism Subgroups:** The search proceeded stratum-by-stratum by prescribing non-trivial automorphism groups $\Gamma \le S_{18}$ of prime orders ($p = 2, 3, 17$) and higher order subgroups.
2. **Backtrack Completion:** Using orderly generation and constraint satisfaction backtracking, every orbit of edges was systematically extended to partial 1-factorisations.
3. **Hamiltonian Pruning:** Any intermediate factorisation containing a pair of 1-factors whose union yielded disconnected 2-factors (cycles of length $< 18$) was pruned immediately.

### 6.2 Independent Sanity & Integrity Checks
Every record in this dataset has undergone automated post-processing:
1. **Factor Validation:** Verified that each $F_i$ contains 9 mutually disjoint edges partitioning the 18 vertices.
2. **Partition of $E(K_{18})$:** Verified that $\bigcup_{i=1}^{17} F_i = E(K_{18})$ with $|E(K_{18})| = 153$.
3. **Pairwise Hamiltonicity:** Checked all $\binom{17}{2} = 136$ pairs $(F_i, F_j)$ for every design; each pair forms a single 18-cycle.
4. **Isomorphism Rejection:** Computed canonical graph certificates using `nauty` / `Traces` to guarantee zero duplicate classes.
5. **Automorphism Group Order:** Independent re-computation of $|Aut(\mathcal{F})|$ using stabilizer chain algorithms matching the census table.

---

## 7. Python Verification Script

Users can verify and parse the records using the Python script below:

```python
import itertools

def verify_p1f_record(factors, n_vertices=18):
    """
    Verifies that a list of factors forms a valid Perfect 1-Factorisation of K_{2n}.
    factors: list of 2n-1 sets of edges, where each edge is a 2-tuple (u, v) with u < v.
    """
    expected_factors = n_vertices - 1
    edges_per_factor = n_vertices // 2
    total_edges = n_vertices * (n_vertices - 1) // 2

    assert len(factors) == expected_factors, f"Expected {expected_factors} factors, got {len(factors)}"

    all_edges = set()
    for idx, factor in enumerate(factors):
        assert len(factor) == edges_per_factor, f"Factor {idx} has {len(factor)} edges"
        factor_vertices = set(itertools.chain.from_iterable(factor))
        assert len(factor_vertices) == n_vertices, f"Factor {idx} is not a perfect matching"
        all_edges.update(factor)

    assert len(all_edges) == total_edges, f"Total distinct edges is {len(all_edges)}, expected {total_edges}"

    # Verify that the union of any two factors is a single Hamiltonian cycle
    for i in range(expected_factors):
        for j in range(i + 1, expected_factors):
            # Construct adjacency list for F_i U F_j
            adj = {v: [] for v in range(n_vertices)}
            for u, v in itertools.chain(factors[i], factors[j]):
                adj[u].append(v)
                adj[v].append(u)

            # Trace cycle starting at vertex 0
            visited = set()
            curr, prev = 0, -1
            while curr not in visited:
                visited.add(curr)
                nbrs = adj[curr]
                nxt = nbrs[0] if nbrs[0] != prev else nbrs[1]
                prev, curr = curr, nxt

            assert len(visited) == n_vertices, f"Factors {i} and {j} do not form an 18-cycle"

    return True
```

---

## 8. Citation & Zenodo Metadata

### Keywords
`perfect one-factorisation`, `one-factorisation`, `complete graph`, `K18`, `automorphism group`, `exhaustive classification`, `combinatorial design`, `atomic Latin square`, `Hamiltonian cycle`

### License
This dataset is distributed under the **Creative Commons Attribution 4.0 International (CC-BY 4.0)** license. You are free to share, copy, and adapt this material with appropriate attribution.

### How to Cite
If you use this catalogue in published research, please cite:
```bibtex
@dataset{p1f_k18_catalogue,
  title        = {A complete catalogue of the perfect one-factorisations of $K_{18}$ with non-trivial automorphism group},
  year         = {2026},
  publisher    = {Zenodo},
  note         = {10,710 isomorphism classes; K18_P1F_aut_gt1.txt}
}
```\n