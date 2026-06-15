# Implementation Plan: Finalizing K18 Cyclic Automorphism Solver

This document outlines the exact steps needed to complete the exhaustive classification of cyclic automorphisms for $K_{18}$ using our current backtracking technique, optimizing it to run in under 20 seconds, and adding case-by-case progress printouts.

---

## 1. Goal & Context
The user has decided **not** to run the identity case ($L=1$, which requires unconstrained search). 
To cover all other mathematically possible cyclic automorphism configurations (including disjoint cycles like $(16)(2)$, $(14)(2)^2$, and the multi-cycle case $(3)^6$), we need to:
1. Eliminate C++ `std::vector` heap allocations in the recursive step to prevent billions of `malloc`/`free` overheads.
2. Enable Search Type 2 and disjoint cycle structures for all even cycle lengths $L < 17$ and for the $L=3$ multi-cycle case.
3. Print progress updates as each case starts.

---

## 2. Step-by-Step Instructions

### Step 1: Optimize Recursion to Stack-Allocated Arrays
In [k18a2support.cpp](file:///c:/TripleSys/Tr71/TripleSys/Source/k18a2support.cpp), rewrite the recursive helper function:
```cpp
void K18A2::CycleBacktrackState::generate_remaining_cycles(int start_idx, const std::vector<uint8_t>& rem, bool* rem_used, uint8_t* alpha_p, int L)
```
to use stack-allocated buffers rather than `std::vector` allocations.

#### Proposed Stack-Based Signature & Logic:
Instead of dynamically passing/creating `std::vector` at every recursive level, pass fixed-size arrays or track index ranges on the stack. For instance:
* `uint8_t rem_nodes[18]` (the out-of-cycle vertices)
* `int rem_size` (the number of out-of-cycle vertices)
* `bool rem_used[18]`
* Stack-allocated variables `int unused_indices[18]`, `uint8_t cycle_nodes[18]`, and `int perm[18]` inside the helper instead of heap vectors.

```cpp
// Example of stack-allocated recursion helper:
void K18A2::CycleBacktrackState::generate_remaining_cycles_stack(
    int start_idx, 
    const uint8_t* rem, 
    int rem_size, 
    bool* rem_used, 
    uint8_t* alpha_p, 
    int L) 
{
    int first_unused = -1;
    for (int i = start_idx; i < rem_size; i++) {
        if (!rem_used[i]) {
            first_unused = i;
            break;
        }
    }
    
    if (first_unused == -1) {
        if (checkPermutationPassed(alpha_p, rem_used, L)) {
            saveAlpha(alpha_p);
        }
        return;
    }
    
    uint8_t v0_rem = rem[first_unused];
    rem_used[first_unused] = true;
    
    for (int d = 1; d <= rem_size; d++) {
        if (L % d != 0) continue;
        
        if (d == 1) {
            alpha_p[v0_rem] = v0_rem;
            generate_remaining_cycles_stack(first_unused + 1, rem, rem_size, rem_used, alpha_p, L);
        } else {
            int unused_indices[18];
            int unused_count = 0;
            for (int j = first_unused + 1; j < rem_size; j++) {
                if (!rem_used[j]) {
                    unused_indices[unused_count++] = j;
                }
            }
            
            if (unused_count >= d - 1) {
                uint8_t cycle_nodes[18];
                cycle_nodes[0] = v0_rem;
                
                int perm[18];
                bool perm_used[18] = { false };
                
                // Perform backtrack arrangement without std::function/heap
                auto arrange_stack = [&](auto& self_fn, int depth) -> void {
                    if (depth == d - 1) {
                        for (int k = 0; k < d - 1; k++) {
                            cycle_nodes[k + 1] = rem[unused_indices[perm[k]]];
                        }
                        for (int k = 0; k < d - 1; k++) {
                            alpha_p[cycle_nodes[k]] = cycle_nodes[k + 1];
                            rem_used[unused_indices[perm[k]]] = true;
                        }
                        alpha_p[cycle_nodes[d - 1]] = cycle_nodes[0];
                        
                        generate_remaining_cycles_stack(first_unused + 1, rem, rem_size, rem_used, alpha_p, L);
                        
                        for (int k = 0; k < d - 1; k++) {
                            rem_used[unused_indices[perm[k]]] = false;
                        }
                        return;
                    }
                    for (int j = 0; j < unused_count; j++) {
                        if (!perm_used[j]) {
                            perm_used[j] = true;
                            perm[depth] = j;
                            self_fn(self_fn, depth + 1);
                            perm_used[j] = false;
                        }
                    }
                };
                arrange_stack(arrange_stack, 0);
            }
        }
    }
    rem_used[first_unused] = false;
}
```

### Step 2: Remove Restrictions on Disjoint Cycle Searches
1. **Remove `L <= 10 && d > 1` restriction**:
   Delete line 117 in [k18a2support.cpp](file:///c:/TripleSys/Tr71/TripleSys/Source/k18a2support.cpp):
   ```cpp
   if (L <= 10 && d > 1) continue; // Prevent combinatorial explosion...
   ```
2. **Enable Search Type 2 for all even L**:
   In `searchCycleLength` (line 335 and 352), change the check:
   ```cpp
   if (L >= 12 && L <= 16 && L % 2 == 0)
   ```
   to:
   ```cpp
   if (L >= 2 && L <= 16 && L % 2 == 0)
   ```

### Step 3: Add L = 3 and L = 2 to the Search List
In `runExhaustiveSearch()` (line 225), update `cycle_lengths` array to include `3` (for the $(3)^6$ structure):
```cpp
int cycle_lengths[] = { 17, 16, 15, 14, 13, 12, 10, 8, 6, 4, 3, 2 };
```

### Step 4: Add Case-by-Case Progress Printouts
In `searchCycleLength`, print status updates to standard output before running the backtrackers:
```cpp
void K18A2::searchCycleLength(int L, std::set<std::vector<uint8_t>>& unique_results, CycleLengthStats& stats) {
    printf("-> Entering Case: L = %d, Search Type 1 (Fixed out-of-cycle points)\n", L);
    fflush(stdout);
    
    // ... run state1 backtrack ...

    if (L >= 2 && L <= 16 && L % 2 == 0) {
        printf("-> Entering Case: L = %d, Search Type 2 (Transposed out-of-cycle points)\n", L);
        fflush(stdout);
        
        // ... run state2 backtrack ...
    }
    
    // ... process automorphisms and record results ...
}
```

### Step 5: Compile & Verify
1. Run `build.bat` from root workspace directory to compile:
   ```cmd
   build.bat
   ```
2. Run the test script:
   ```cmd
   NewTests\run18a2.bat
   ```
3. Verify the final ASCII table output contains all 17 cases, is 100% complete, and displays the correct counts for `Unique Reps` and `Isom Classes` (confirming the GB, GK, and any new classes).
