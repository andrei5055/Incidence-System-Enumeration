# Implementation Plan: Reverting Starter Mapping Constraints, Adding Missing Cases, and Fixing Table Printing in K16A2

This plan addresses:
1. Reverting the incorrect starter mapping constraints in `K16A2` that caused the 3 valid results for $L=3$ to be rejected.
2. Adding the remaining cycle lengths (13, 11, and 7) to make the cyclic solver search exhaustive for $K_{16}$.
3. Fixing a table printing bug where $L=8$ (which is `NHALF` for $K_{16}$) was incorrectly classified and suppressed as "Theoretically Impossible".

Any code changes will include a comment explaining why the change is specific to $K_{16}$ ($NP=16$, $NM=15$) compared to $K_{18}$ ($NP=18$, $NM=17$).

## Proposed Changes

### [K16A2 Solver]

#### [MODIFY] [k16A2support.cpp](file:///c:/TripleSys/Tt4/TripleSys/Source/k16A2support.cpp)
1. **Revert `validateDefinedOrbits`**:
   - Remove the K16-specific starter mapping constraints (Constraint 1 and Constraint 2) from `K16A2::CycleBacktrackState::validateDefinedOrbits`.
   - Restore the clean orbit compatibility check that matches `K18A2`.
   - Add a comment explaining that this revert makes `K16A2` match the core `K18A2` algorithm (which does not enforce starter-specific transitions during cycle backtrack).
2. **Correct remaining cycle limits in `searchCycleLength`**:
   - Replace `L <= 16` with `L <= 14` on lines 900, 939, 977, and 978.
   - Add a comment explaining that `14` is specific to $K_{16}$ since the maximum even cycle length is $NM-1 = 14$ (whereas it is $16$ for $K_{18}$).
3. **Add Cycle Lengths (13, 11, 7) and fix table printing**:
   - Update `cycle_lengths` array in `runExhaustiveSearch()` to include `13, 11, 7`.
   - Remove hardcoded classification checks for `L == 11 || L == NHALF || L == 7` and `L == 15 || L == 13` from the table printing loop, so actual stats are printed for all run cases.
   - Add comments explaining that $L=13, 11, 7$ are included for $K_{16}$ because $K_{16}$ has prime automorphism orders of 7, 5, 3, 2.

---

## Verification Plan

### Automated Tests
- Run `build.bat nopause` to compile the release build.
- Run `run16a2.bat nopause` to execute the solver.
- Verify that:
  - The 3 valid results for $L=3$ are successfully found and recorded in the summary table.
  - All cycle lengths from 15 down to 2 are run and printed in the table.
  - The stats for $L=8$ are printed correctly (not marked theoretically impossible).
  - No compilation warnings or errors are introduced.
