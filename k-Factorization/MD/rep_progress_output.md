# Reading the representative-method print-outs (K16 / K18 / K20)

This document explains every line the automorphism-**representative** classifier prints
(`[K16-REP]` / `[K18-REP]` / `[K20-REP]` and the indented `[rep]` progress lines), with
real examples, and shows how to tell at a glance whether a run is **ramping toward
completion**, has hit an **intractable plateau**, or is **finding-but-not-terminating**.

All output is produced by `runRepresentativeMethod` in `k16A2rep.cpp` / `k18a2rep.cpp` /
`k20a2rep.cpp`. It is gated by `m_bPrint` (the job's print flag); the periodic `[rep]`
line is emitted **at most once every 30 s**. Worker count = `KThreads`.

---

## 1. The line types

A run prints, in order:

```
[K18-REP] order=4 : 20 cycle types, kThreads=32, symCap=2000000          (A) run header
[K18-REP] HARVEST MODE: stop after 10 distinct classes (...)             (B) only if a target is set
[K18-REP]  type 4 4 4 4      rawTasks 194 uniq 194 reps 194 stab0 3 (dedup 0.1s) -> cover...   (C) per-type, BEFORE the cover
    [rep] o4 type 5/20 "4 4 4 4 " reps 12/194  elapsed=30s nodes=500736 (16684/s) classes=8    (D) every 30s DURING a cover
[K18-REP]  type 4 4 4 4      |C|=6144    rawTasks 194     reps 194     new +131 (cum 131) 15.9s (E) per-type, AFTER the cover
...                                                                       (C/D/E repeat per cycle type)
[K18-REP] |Aut| histogram: 4:144 8:22 16:12 272:1                        (F) final, all types done
[K18-REP] order=4 TOTAL DISTINCT CLASSES = 179  (357.0s)                 (G) final total
```

### (A) Run header
```
[K18-REP] order=<order> : <N> cycle types, kThreads=<workers>, symCap=<cap>
```
- **order** — the automorphism order being classified this pass.
- **N cycle types** — how many cycle types of that order exist (each is searched in turn).
- **symCap** — centralizer-enumeration cap. A type with `|C(alpha)| <= symCap` uses the fast
  *enumerable* path; a larger one uses the *over-cap* path (see §2). Override with `REP_SYMCAP`.

### (B) Harvest-mode banner (only with an `order:target`)
```
[K18-REP] HARVEST MODE: stop after <target> distinct classes (COMPLETE only if <target> == true count, else a LOWER BOUND)
```
Printed when you pass `REP_ORDER=2:60` / `REP_ORDERS=2:60,4,5,7`. The run stops as soon as
`target` distinct classes are found and prints `HARVESTED ...` instead of `TOTAL ...` (see §3).

### (C) Per-type **pre-cover** line
```
[K18-REP]  type <typestr>  rawTasks <R> uniq <U> reps <reps> stab0 <S> (dedup <t>s) -> cover...
```
- **typestr** — the cycle type, e.g. `4 4 4 4` = 4⁴ (+ implicit fixed points to fill N), `2 2 2 2 2 2 2 2 2` = 2⁹.
- **rawTasks R** — raw first-orbit matchings generated at the root.
- **uniq U** — distinct alpha-orbits among them.
- **reps** — representative tasks after centralizer de-duplication. **This is the number of
  independent root subtrees / parallel work units.** Few reps + big subtrees = the hard case.
- **stab0 S** — size of the first rep's stabilizer (enumerable: number of group elements;
  over-cap: number of generators). A large `stab0` on an enumerable type means heavy per-node
  symmetry work (e.g. K16 2⁷·1² has stab0 ≈ 1.29 M).
- **dedup t** — seconds spent collecting + de-duplicating before the cover starts.

### (D) The 30-second `[rep]` progress line  *(K16 / K18 — the rich form)*
```
    [rep] o<order> type <i>/<N> "<typestr>" reps <done>/<created>  elapsed=<s>s nodes=<n> (<rate>/s) classes=<c>
```
This is the line you watch. Fields:
- **o\<order\>** — current order.
- **type i/N** — currently on cycle type *i* of *N*.
- **"typestr"** — the current type.
- **reps done/created** — **the key completion gauge** (meaning depends on the path, see §2):
  - *enumerable type:* `created` = number of root reps; `done` = root subtrees finished.
  - *over-cap type:* `created` = work-queue tasks (partial covers) spawned **so far**;
    `done` = tasks expanded/completed so far. `created` keeps growing while the tree is explored.
- **elapsed** — seconds since the whole run started (across all types, not just this one).
- **nodes** — cumulative exact-cover nodes over all workers.
- **(rate/s)** — nodes per second since the previous `[rep]` line (throughput / how busy the cores are).
- **classes** — distinct classes found so far this pass. ⚠️ This is an **approximate** cross-worker
  count (it can slightly over-count duplicates found by two workers at once); the **exact** number
  is the final histogram/total. For trend-watching the approximation is fine.

> **K20 prints a simpler line** (older engine): `    [rep] elapsed=<s>s nodes=<n> classes=<c>` —
> no `type i/N` and **no `reps done/created`**. For K20 use `nodes`, `rate`, and `classes`
> trends plus the per-type (C)/(E) lines; the reps gauge is K16/K18 only.

### (E) Per-type **post-cover** summary
```
[K18-REP]  type <typestr>  |C|=<c>  rawTasks <R>  reps <reps>  new +<k> (cum <C>)  <t>s
```
- **|C|** — `|C(alpha)|` if it was enumerated, or **`0` = over-cap** (the over-cap flag; that type
  ran on the generator / work-queue path). **`|C|=0` is the single most important field for
  spotting a potentially intractable type.**
- **new +k** — distinct classes this type contributed.
- **cum C** — cumulative distinct across types so far this order.
- **t** — seconds this type's cover took. A type that prints this line **finished** (proved complete).

### (F)/(G) Final lines
```
[K18-REP] |Aut| histogram: 4:144 8:22 16:12 272:1
[K18-REP] order=4 TOTAL DISTINCT CLASSES = 179  (357.0s)
```
- **histogram** — `|Aut|:count` over the distinct classes (e.g. 144 classes with |Aut|=4, etc.).
- **TOTAL DISTINCT CLASSES** — the proven, de-duplicated count for this order.
- In harvest mode this last line is replaced by
  `order=<o> HARVESTED <n> classes (target <t> reached -- NOT proven complete) (<s>s)`.

---

## 2. Enumerable vs over-cap — why the `reps` gauge means two different things

| | enumerable type (`\|C\| > 0`) | over-cap type (`\|C\|=0`) |
|---|---|---|
| symmetry | full `C(alpha)` materialized | only generators + BSGS |
| fan-out | one worker per root rep | shared **work queue** of partial covers |
| `reps created` | fixed = number of root reps | **grows** as the tree is split |
| `reps done` | root subtrees finished | queue tasks expanded |
| "done" reaching "created" | type finishes | type finishes |

Knowing which path a type is on (from `|C|=...` in line E, or from whether `created` grows)
tells you how to read `reps done/created`.

---

## 3. Reading the trend — ramp, plateau, find-but-no-terminate

Watch four things across consecutive `[rep]` lines: **`classes`**, **`reps done`**, **`reps created`**,
and the **`(rate/s)`**. The combinations map to four situations.

### 3a. Healthy completion (the common, good case)
Most types finish in well under a second and you only ever see the (C)/(E) pair, e.g. K18 order-17:
```
[K18-REP]  type 17             rawTasks 17 uniq 17 reps 17 stab0 17 (dedup 1.1s) -> cover...
[K18-REP]  type 17             |C|=17       rawTasks 17      reps 17      new +2 (cum 2) 0.0s
[K18-REP] |Aut| histogram: 17:1 272:1
[K18-REP] order=17 TOTAL DISTINCT CLASSES = 2  (1.1s)
```
**Signal:** the type prints its (E) summary quickly → it terminated. No `[rep]` lines at all means
every type finished inside the 30 s window. Order-4 (the slowest tractable K18 order) finished in
357 s with `TOTAL DISTINCT CLASSES = 179`.

### 3b. Ramp toward completion (a long but finite type)
`reps done` climbs steadily toward `reps created`, `created` stops growing much, and the type
eventually prints its (E) line. This is what a big-but-tractable over-cap type looks like: the gap
`created − done` **shrinks** over time. If `done/created` is closing, let it run.

### 3c. **Plateau = intractable** (the pattern to recognize and abort)
`classes` is flat, `reps done` crawls, and `reps created` **keeps outrunning it with no sign of
catching up**. Real example — K18 order-2, type 2⁹, on a 32-core PC:
```
[rep] o2 type 9/9 "2 2 2 2 2 2 2 2 2 " reps   14/205691  elapsed=30s    nodes=500736    (16684/s) classes=0
[rep] o2 type 9/9 "2 2 2 2 2 2 2 2 2 " reps   53/205691  elapsed=15109s nodes=384692224 (25393/s) classes=0
```
In **4.2 hours**: `done` went 14 → 53, `created` stuck at 205 691 (the queue hit its cap and is
draining locally), `classes` never left 0, rate steady ~25 k/s. **Read:** the cores are 100% busy
(good rate) but the tree is astronomically larger than the work done — `done` will never reach
`created`. This type will not terminate; harvest it or skip it. `classes=0` here does **not** mean
"the answer is 0" — it means the search hasn't reached a single *complete* cover yet.

### 3d. Finding but not terminating (productive, still can't finish)
`classes` **grows** but `reps done` stays stuck (often `0/<small>`), because the few root subtrees
are each enormous. Real example — K16 order-2, type 2⁷·1² (enumerable, stab0 ≈ 1.29 M):
```
[rep] o2 type 7/8 "2 2 2 2 2 2 2 " reps 0/4  elapsed=90s  nodes=15381504  (175920/s) classes=30
[rep] o2 type 7/8 "2 2 2 2 2 2 2 " reps 0/4  elapsed=210s nodes=36212736  (160244/s) classes=44
[rep] o2 type 7/8 "2 2 2 2 2 2 2 " reps 0/4  elapsed=750s nodes=126758912 (135313/s) classes=58
```
**Read:** `classes` rising (8→30→44→58) ⇒ the search **is** finding real classes. But `reps 0/4`
after 12 min ⇒ not one of the 4 root subtrees has finished, so it cannot *prove* it found them all.
This is the case where a **harvest target** is the right tool: `REP_ORDER=2:60` (or a safe `2:55`)
collects the classes and stops, instead of grinding forever on the exhaustion proof.

---

## 4. Quick decision guide

| What you see across `[rep]` lines | Meaning | Action |
|---|---|---|
| Type prints its `(E)` summary line | terminated, proven complete | nothing — read `new +k` / `cum` |
| `done` rising **toward** `created`, gap shrinking | ramping to completion | let it run |
| `classes` flat (e.g. 0), `created` ≫ `done` and gap **widening** | intractable plateau (e.g. 2⁹) | abort / harvest / skip the type |
| `classes` **rising** but `done` stuck at `0/小` | finding, can't terminate | use a harvest target `order:target` |
| `(rate/s)` near 0 or collapsing | cores idle / starved | check `KThreads`, check it's the fixed build |
| `nodes` rising, `classes` rising, then `HARVESTED` printed | target reached | result is a **lower bound** (or exact if target = known count) |

**Rules of thumb**
- `|C|=0` in the `(E)` line = over-cap type — the candidates for the slow plateau cases.
- The single best "is it finishing?" signal is the **`reps done/created` trend**, not elapsed time.
- `classes` in `[rep]` is approximate; trust the final **histogram / TOTAL** for the exact count.
- A type that has printed its `(E)` summary is *done*; only the type named in the live `[rep]`
  line is still running.

---

## 5. Worked example — a full order-4 K18 pass (abridged; inner counts illustrative)

The shape of a complete, tractable pass (verified outer numbers: `+131`, `+48`, `cum 179`, the
`2 4 4 4 4` holdout taking ~224 s, total 357 s; the `reps`/`dedup`/intermediate values are
illustrative of the format):

```
[K18-REP] order=4 : 20 cycle types, kThreads=32, symCap=2000000
[K18-REP]  type 4 4 4 4      rawTasks ... reps 194 stab0 ... (dedup ...s) -> cover...
[K18-REP]  type 4 4 4 4      |C|=12288   ... reps 194  new +131 (cum 131) 15.9s     <- big enumerable type, finished
... (18 more types, most new +0, sub-second) ...
[K18-REP]  type 2 4 4 4 4    ... -> cover...                                         <- 4^4*2^1, the slow holdout
    [rep] o4 type 20/20 "2 4 4 4 4 " reps ... classes=...                           <- watched until it finished
[K18-REP]  type 2 4 4 4 4    |C|=12288   ... new +48 (cum 179) 224.0s               <- terminated after 224s
[K18-REP] |Aut| histogram: 4:144 8:22 16:12 272:1
[K18-REP] order=4 TOTAL DISTINCT CLASSES = 179  (357.0s)
```
Here type `4 4 4 4` (= 4⁴·1²) and `2 4 4 4 4` (= 4⁴·2¹) are the two heavy enumerable types; the
rest are trivial. The final histogram is the deliverable. Contrast with order-2, whose `2⁸·1²` and
`2⁹` types never reach an `(E)` line (§3c) — that is the signature of the intractable order.

---

## 6. Notes

- **Cadence / cost:** the `[rep]` line is throttled to one print per 30 s and flushed every ~1024
  nodes per worker, so it is essentially free and never floods the log.
- **K20 difference:** K20 still uses the older inline progress line
  (`[rep] elapsed=.. nodes=.. classes=..`, no `type`/`reps` fields). Use the per-type `(C)/(E)`
  lines and the `nodes`/`classes` trend there.
- **Harvest mode** (`order:target`) is documented in the batch headers and `MD` notes; in the
  output it shows the `HARVEST MODE` banner (B) and the `HARVESTED ...` final line instead of
  `TOTAL DISTINCT CLASSES`. A harvested count is a lower bound unless `target` equals the true count.
