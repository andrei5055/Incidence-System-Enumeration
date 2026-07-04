# TBD — K18 order-2 (2^9) next steps (saved 2026-06-27)

## Verdict: K18 2^9 is intractable for BOTH compute engines

- **k18a2rep (main rep engine — the strong one: exact per-node Stab_C(P) dedup,
  parallel, bounded memory, harvest):** order-2 is intractable. Per `run18a2full.bat`
  (verified 2026-06-27): the two dense over-cap involution types **2^8·1²** and **2^9**
  do **not** plateau (2^8·1² alone expands a 223k-task frontier with no convergence).
  Running `REP_ORDER=2` → the 7 small types finish in seconds at +0, then it **hangs** on
  2^8·1² / 2^9 (kill manually). So a complete L=2 count will NOT finish in 24h.
- **From-scratch engine (`scratchpad/engine.cpp`):** also can't do K18 2^9 (~600 nodes/s,
  ~3% pruning, tree ≈ full search). See [k18_engine_status.md](k18_engine_status.md).

If even the exact-stabilizer rep engine doesn't plateau on 2^9, the search tree is
astronomical. **Compute route is exhausted for 2^9.**

## The math has already reduced order-2 to exactly two cases

Our **≤2-fixed-points lemma** (an involution automorphism of a P1F fixes ≤2 vertices —
[order4_parity_theorem.md](order4_parity_theorem.md)) proves every order-2 type
2^k·1^(18−2k) with k≤7 (≥4 fixed points) is **EMPTY**. This certifies the engine's
"+0 in seconds" for the 7 small types — and extends it. So the ONLY undecided order-2
types for K18 are:
- **2^8·1²** (f=2 fixed points) — the mixed case, open.
- **2^9** (f=0, fully fpf) — the open case our parity theorem leaves for n=9 odd.

Everything else (order-2) is settled by theorem.

## HIGHEST-VALUE NEXT STEP: prove K18 2^9 = ∅ (or find an example) by ARGUMENT

This is the only route left and the handles are strong:
1. The ≤2-fixed lemma + the order-2 parity theorem already killed every other order-2 type.
2. The user's most-symmetric K18 result (|Aut|=272 = AGL(1,17)) has **no** 2^9 involution
   (all involutions are 2^8·1²) — see `scratchpad/check_result.cpp`. Weakly suggests
   2^9 may be empty.
3. The rigid structure is fully pinned: a 2^9-invariant K18 P1F = exactly **9 α-fixed
   factors (one diagonal each) + 4 swapped pairs** (proved). Look for an obstruction in
   this structure (a counting / monodromy / parity argument specific to n=9, analogous to
   how the antipodal-rotation lemma killed n even).
4. If 2^9 turns out non-empty, the goal becomes constructing/exhibiting one example.

Also worth settling by the same kind of argument: **2^8·1²** (the last mixed case).

## OPTIONAL fast experiment: existence probe

Run k18a2 once with `REP_ORDER=2:1` (harvest — stop after the first class). If a 2^9 (or
2^8·1²) P1F exists, it may surface one quickly → then we'd HAVE an example to study/pin.
If it churns, more evidence toward empty. (Caveat: prior "hangs" were full runs, not
harvest=1 — confirm whether harvest changes it.)

## Pointers
- Full session log: [k18_engine_status.md](k18_engine_status.md)
- Theorems: [order2_parity_theorem.md](order2_parity_theorem.md),
  [order4_parity_theorem.md](order4_parity_theorem.md)
- Tools: `scratchpad/{engine,canon_ir,rigid,parity_check,order4_check,check_result}.cpp`
