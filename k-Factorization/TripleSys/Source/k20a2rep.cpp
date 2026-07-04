// =============================================================================
// k20a2rep.cpp — automorphism-REPRESENTATIVE method for classifying K20 Perfect
// 1-Factorizations (P1Fs) by a given automorphism order. N=20 sibling of
// k16A2rep.cpp / k18a2rep.cpp (the engine is N-generic; this file differs only
// in N and names). K20 has no cyclic engine at all -- this rep method is its
// only solver path.
//
// Method: unseeded and complete. For ONE representative permutation alpha per
// cycle TYPE of the target order, it finds ALL alpha-invariant P1Fs by an exact
// cover of the 190 edges of K20 with alpha-invariant factor-orbits, then
// deduplicates across types with a labeling-independent canonical key. Fixed-
// point-free cycle types (which cyclic generators structurally miss) are
// included, so the classification is complete.
//
// The engine is validated on the smaller cases by the same code: K16 |Aut|>2 =
// {3:19,5:5,7:4,14:1,15:1}=30 (Wanless catalogue), K18 order-4 = 179, order-17
// -> 2. Known K20 anchors: order-19 -> 7 classes {19:3,57:1,171:2,342:1}
// (342 = AGL(1,19)); the order-4 parity theorem (MD/order4_parity_theorem.md)
// says K20 admits NO automorphism of order divisible by 4, so order-4 must
// return 0. Isomorph rejection is by ORDERLY GENERATION: at every cover level
// the candidate factor-orbits covering the lowest uncovered edge are
// deduplicated by the current centralizer stabilizer G, and each chosen orbit O
// recurses with its sub-stabilizer Stab_G(O). Sound for class counting because
// any g in G = Stab(partial cover) is a class-preserving bijection between
// sibling subtrees; combined with the shrinking stabilizer chain this collapses
// the dense high-centralizer types that a naive lex-leader prune cannot
// terminate.
//
// Two symmetry-reduction paths, by whether C(alpha) fits the enumeration cap:
//   * ENUMERABLE (|C(alpha)| <= cap): C(alpha) is materialized; cover() carries
//     the stabilizer as indices into it and dedups candidates directly.
//   * OVER-CAP (sparse / huge S_k on fixed points / every order-2 type): C(alpha)
//     is held only as a generating set + BSGS (cgtEngine.h). coverGen() recomputes
//     the EXACT setwise stabilizer Stab_C(P) of the current partial cover per node
//     (full-strength dedup at every depth), and splitNode() feeds a shared
//     work-queue so the few root reps still saturate all worker threads.
//
// Parallelism: enumerable types fan their independent first-orbit reps across
// `kThreads` workers via an atomic index; over-cap types use a shared work-queue
// of partial covers (each worker pops one, expands it one deduped level, pushes
// the children back). The per-type/global dedup is merged after the workers join.
//
// Integration rules honored here:
//   * The ONLY entry point is the private member K20A2::runRepresentativeMethod,
//     called internally from K20A2::runExhaustiveSearch. Everything else lives in
//     the anonymous namespace below, so nothing new has external linkage.
//   * Worker count comes from K20A2::kThreads.
//   * Every print is gated by K20A2::m_bPrint (forwarded as g_bprint).
// =============================================================================

#include "k20a2.h"
#include "cgtEngine.h"   // shared computational-group-theory engine (BSGS, edgeStabGens, setwiseStab)
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <array>
#include <deque>
#include <set>
#include <map>
#include <unordered_map>
#include <string>
#include <algorithm>
#include <chrono>
#include <functional>
#include <thread>
#include <atomic>
#include <mutex>

namespace {  // everything here is file-local (no external linkage)

// ---- problem-size constants (single knob: N) --------------------------------
constexpr int N      = 20;             // players / graph vertices (K20)
constexpr int NM     = N - 1;          // 19 = factors in a 1-factorization
constexpr int NHALF  = N / 2;          // 10 = edges in one perfect matching (factor)
constexpr int NEDGES = N * (N - 1) / 2;// 190 = edges of K20
constexpr int INF    = 1 << 20;        // "no edge / not covered" sentinel (> any edge id)

typedef std::array<uint8_t, N> Match;  // a factor in adjacency form: m[u] = partner of u
typedef std::array<uint8_t, N> Perm;   // a vertex permutation: p[u] = image of u

// ---- edge-id tables (filled once by initEdges) ------------------------------
int     eidT[N][N];          // eidT[u][v] = canonical id (0..189) of edge {u,v}; symmetric
uint8_t edgeU[NEDGES];       // edgeU[e] = lower  endpoint of edge id e
uint8_t edgeV[NEDGES];       // edgeV[e] = higher endpoint of edge id e

void initEdges() {
    int c = 0;               // running edge id, assigned in (u<v) lexicographic order
    for (int i = 0; i < N; i++)
        for (int j = i + 1; j < N; j++) {
            eidT[i][j] = eidT[j][i] = c;   // both orientations map to the same id
            edgeU[c] = (uint8_t)i;
            edgeV[c] = (uint8_t)j;
            c++;
        }
}

// two factors are "perfect" iff their union is a single Hamiltonian N-cycle
inline bool is_perfect(const Match& a, const Match& b) {
    int u = 0;               // walk the alternating a/b path starting from vertex 0
    for (int i = 0; i < NHALF; i++) {
        u = a[u]; u = b[u];
        if (u == 0 && i < NHALF - 1) return false;  // closed early => more than one cycle
    }
    return u == 0;           // returned to 0 exactly on the last step => single cycle
}

// image of a factor M under permutation alpha: (alpha M)[alpha u] = alpha(M[u])
inline Match applyAlpha(const Match& m, const uint8_t* al) {
    Match r;                 // result factor (adjacency form)
    for (int u = 0; u < N; u++) r[al[u]] = al[m[u]];
    return r;
}

// permutation composition (a o b)(x) = a(b(x)), and inverse. Group elements act on a
// factor via applyAlpha; the actions compose as
//   applyAlpha(applyAlpha(m, b), a) == applyAlpha(m, composePerm(a, b)).
inline Perm composePerm(const Perm& a, const Perm& b) { Perm r; for (int x = 0; x < N; x++) r[x] = a[b[x]]; return r; }
inline Perm invPerm(const Perm& a) { Perm r; for (int x = 0; x < N; x++) r[a[x]] = (uint8_t)x; return r; }

// Canonical serialization of the alpha-orbit of factor m (the set {m, alpha m, ...}),
// as the sorted concatenation of the orbit's factors. Used to (a) group first-orbit
// tasks into C(alpha)-orbits and (b) test whether a centralizer element fixes an orbit.
// Note g commutes with alpha, so g(orbit of m) == orbit of g(m); hence applying g to m
// and re-taking the orbit gives g's image of the whole orbit.
std::string orbitKey(const Match& m, const uint8_t* alpha, int order) {
    std::vector<Match> orb;  // the alpha-orbit of m
    orb.push_back(m);
    Match cur = m;
    for (int s = 1; s < order; s++) { cur = applyAlpha(cur, alpha); if (cur == m) break; orb.push_back(cur); }
    std::sort(orb.begin(), orb.end());
    std::string s; s.reserve(orb.size() * N);
    for (auto& F : orb) for (int u = 0; u < N; u++) s.push_back((char)F[u]);
    return s;
}

// convert a factor from adjacency form to the project's "src" pair-list form
void adj_to_src(const uint8_t* adj, unsigned char* src) {
    bool visited[N] = { false };   // visited[u] = vertex u already emitted into src
    src[0] = 0; src[1] = adj[0];
    visited[0] = true; visited[adj[0]] = true;
    int idx = 2;                   // write cursor into src
    for (int u = 1; u < N; u++) {
        if (!visited[u]) {
            src[idx] = (uint8_t)u; src[idx + 1] = adj[u];
            visited[u] = true; visited[adj[u]] = true;
            idx += 2;
        }
    }
}

// ---- centralizer C(alpha) = { g : g*alpha == alpha*g } ----------------------
// Structure: g maps each alpha-cycle onto an alpha-cycle of the SAME length, with
// a rotation. |C(alpha)| = prod over cycle-lengths l of ( l^{m_l} * m_l! ). If
// that exceeds `cap` we do NOT enumerate (return false): such types have a huge
// centralizer but their search rejects almost instantly, so breaking is
// unnecessary. Returns true iff `out` was filled with all of C(alpha).
bool buildCentralizer(const uint8_t* alpha, std::vector<Perm>& out, long long cap) {
    std::vector<std::vector<int>> cyc;   // cyc[k] = vertices of cycle k, in alpha order
    bool vis[N] = { false };             // vis[u] = vertex u already assigned to a cycle
    for (int i = 0; i < N; i++) {
        if (vis[i]) continue;
        std::vector<int> c;              // the cycle currently being traced
        int x = i;
        do { vis[x] = true; c.push_back(x); x = alpha[x]; } while (x != i);
        cyc.push_back(c);
    }
    std::vector<int> lens;               // lens[k] = length of cycle k
    for (auto& c : cyc) lens.push_back((int)c.size());

    // size check against the cap (avoid building a multi-million-element group)
    {
        std::vector<int> uq;             // distinct cycle lengths present
        for (int L : lens) if (std::find(uq.begin(), uq.end(), L) == uq.end()) uq.push_back(L);
        long long sz = 1;                // running |C(alpha)| product
        for (int L : uq) {
            int m = 0; for (int x : lens) if (x == L) m++;   // m = #cycles of length L
            for (int i = 0; i < m; i++)  { sz *= L; if (sz > cap) { out.clear(); return false; } } // l^m
            for (int i = 2; i <= m; i++) { sz *= i; if (sz > cap) { out.clear(); return false; } } // m!
        }
    }

    std::vector<int> uniqLens;           // distinct cycle lengths (enumeration classes)
    for (int L : lens) if (std::find(uniqLens.begin(), uniqLens.end(), L) == uniqLens.end()) uniqLens.push_back(L);
    out.clear();
    Perm g; for (int i = 0; i < N; i++) g[i] = (uint8_t)i;   // permutation being built (start = identity)

    std::function<void(int)> recLen = [&](int li) {          // recurse over length classes
        if (li == (int)uniqLens.size()) { out.push_back(g); return; }
        int L = uniqLens[li];
        std::vector<int> idx;            // indices of the cycles having this length L
        for (int k = 0; k < (int)cyc.size(); k++) if (lens[k] == L) idx.push_back(k);
        int m = (int)idx.size();         // number of length-L cycles
        std::vector<int> perm(m);        // perm[s] = which length-L cycle source s maps to
        for (int i = 0; i < m; i++) perm[i] = i;
        std::vector<int> off(m, 0);      // off[s] = rotation offset applied to source cycle s

        std::function<void(int)> recOff = [&](int t) {       // choose offsets per source cycle
            if (t == m) {
                for (int s = 0; s < m; s++) {
                    const std::vector<int>& src = cyc[idx[s]];
                    const std::vector<int>& tgt = cyc[idx[perm[s]]];
                    for (int j = 0; j < L; j++) g[src[j]] = (uint8_t)tgt[(j + off[s]) % L];
                }
                recLen(li + 1);
                return;
            }
            for (int o = 0; o < L; o++) { off[t] = o; recOff(t + 1); }
        };
        std::function<void(int)> recPerm = [&](int t) {      // choose the cycle->cycle bijection
            if (t == m) { recOff(0); return; }
            for (int v = t; v < m; v++) { std::swap(perm[t], perm[v]); recPerm(t + 1); std::swap(perm[t], perm[v]); }
        };
        recPerm(0);
    };
    recLen(0);
    return true;
}

// A small GENERATING SET of C(alpha) (always cheap to build, even when C(alpha) itself is far
// too large to enumerate). Used to walk C(alpha)-orbits of tasks for the over-cap types whose
// centralizer is dominated by the S_k symmetry on alpha's fixed points (k! is astronomical but
// trivially generated). For each cycle-length class l with m cycles, the wreath C_l wr S_m is
// generated by: rotate one cycle (if l>1), swap two cycles, and cycle all m cycles.
void buildGenerators(const uint8_t* alpha, std::vector<Perm>& gens) {
    gens.clear();
    std::vector<std::vector<int>> cyc;   // cyc[k] = vertices of cycle k, in alpha order
    bool vis[N] = { false };
    for (int i = 0; i < N; i++) {
        if (vis[i]) continue;
        std::vector<int> c; int x = i;
        do { vis[x] = true; c.push_back(x); x = alpha[x]; } while (x != i);
        cyc.push_back(c);
    }
    std::vector<int> lens; for (auto& c : cyc) lens.push_back((int)c.size());
    std::vector<int> uniqLens;
    for (int L : lens) if (std::find(uniqLens.begin(), uniqLens.end(), L) == uniqLens.end()) uniqLens.push_back(L);
    auto identityPerm = []() { Perm g; for (int i = 0; i < N; i++) g[i] = (uint8_t)i; return g; };
    for (int L : uniqLens) {
        std::vector<int> idx;            // cycle indices of this length
        for (int k = 0; k < (int)cyc.size(); k++) if (lens[k] == L) idx.push_back(k);
        int m = (int)idx.size();
        if (L > 1) {                     // rotate cycle idx[0] by one step (base C_l generator)
            Perm g = identityPerm(); const std::vector<int>& K = cyc[idx[0]];
            for (int j = 0; j < L; j++) g[K[j]] = (uint8_t)K[(j + 1) % L];
            gens.push_back(g);
        }
        if (m > 1) {                     // swap cycles idx[0],idx[1] aligned (transposition in S_m)
            Perm g = identityPerm(); const std::vector<int>& A = cyc[idx[0]]; const std::vector<int>& B = cyc[idx[1]];
            for (int j = 0; j < L; j++) { g[A[j]] = (uint8_t)B[j]; g[B[j]] = (uint8_t)A[j]; }
            gens.push_back(g);
            if (m > 2) {                 // cycle all m cycles (m-cycle in S_m): idx[s] -> idx[s+1]
                Perm h = identityPerm();
                for (int s = 0; s < m; s++) { const std::vector<int>& P = cyc[idx[s]]; const std::vector<int>& Q = cyc[idx[(s + 1) % m]]; for (int j = 0; j < L; j++) h[P[j]] = (uint8_t)Q[j]; }
                gens.push_back(h);
            }
        }
    }
}


// Partition a set of candidate factor-orbits into orbits under the group G = <gens>, used
// for the OVER-CAP types whose centralizer C(alpha) cannot be enumerated (so the group is
// only available as a small generating set). Each candidate is a factor-orbit given by a
// representative matching reps[i], its serialization keys[i], and the lookup keyIdx (orbit
// key -> index). G acts on a factor-orbit O via applyAlpha (every g in C(alpha) commutes
// with alpha, hence maps alpha-orbits to alpha-orbits); the image orbit is identified by
// orbitKey. For each G-orbit we output one representative index and a generating set of (a
// subgroup of) its setwise stabilizer Stab_G(O), obtained as the Schreier generators along
// a candidate-restricted BFS tree.
//
// Soundness for class counting: every Schreier generator lies in Stab_G(O) (it maps O to
// itself), so the recursion is always handed a genuine subgroup of the true stabilizer;
// candidate marking only ever collapses two orbits that are genuine G-images of one another;
// and the global canonKey is the final dedup. The BFS is restricted to the candidate set, so
// the stabilizer it recovers may be an under-approximation of Stab_G(O) -- that can only cost
// extra search work, never a missed class.
void schreierDedup(const std::vector<Match>& reps, const std::vector<std::string>& keys,
                   const std::unordered_map<std::string, int>& keyIdx,
                   const std::vector<Perm>& gens, int order, const uint8_t* alpha,
                   std::vector<int>& outRepIdx, std::vector<std::vector<Perm>>& outStab) {
    (void)keys;                                          // (kept for signature symmetry/readability)
    int n = (int)reps.size();
    Perm idp; for (int x = 0; x < N; x++) idp[x] = (uint8_t)x;   // identity permutation
    std::vector<char> marked(n, 0);                      // marked[j]: j already assigned to some output orbit
    for (int i = 0; i < n; i++) {
        if (marked[i]) continue;
        std::vector<Perm> transv(n);     // transv[j] maps orbit_i -> orbit_j (valid where visited[j])
        std::vector<char> visited(n, 0);
        transv[i] = idp; visited[i] = 1; marked[i] = 1;
        std::vector<int> q; q.push_back(i);              // BFS queue of candidate indices in orbit_i's class
        std::set<Perm> stabSet;                          // distinct nontrivial Schreier gens of Stab_G(orbit_i)
        for (size_t qh = 0; qh < q.size(); qh++) {
            int x = q[qh];
            for (const auto& g : gens) {
                std::string ik = orbitKey(applyAlpha(reps[x], g.data()), alpha, order);
                auto it = keyIdx.find(ik);
                if (it == keyIdx.end()) continue;        // g moves this orbit out of the candidate set
                int t = it->second;
                Perm gx = composePerm(g, transv[x]);     // maps orbit_i -> orbit_t
                if (!visited[t]) { visited[t] = 1; transv[t] = gx; marked[t] = 1; q.push_back(t); }
                else { Perm s = composePerm(invPerm(transv[t]), gx); if (s != idp) stabSet.insert(s); }  // in Stab_G(orbit_i)
            }
        }
        outRepIdx.push_back(i);
        outStab.emplace_back(cgt::reduceGens<N>(std::vector<Perm>(stabSet.begin(), stabSet.end())));  // bound generator growth
    }
}

// ---- read-only data shared by all worker threads for one cycle type ---------
struct RepShared {
    uint8_t alpha[N];                 // representative automorphism (vertex perm)
    int order = 4;                    // target automorphism order
    std::vector<Perm> Calpha;         // all elements of C(alpha) (empty if over the cap)
    std::vector<Perm> Cinv;           // Cinv[i] = inverse permutation of Calpha[i]
    std::vector<int>  rootActive;     // {0,1,...,|Calpha|-1}: the root active-gamma set
    std::vector<Perm> rootCgens;      // generators of C(alpha) (over-cap path)
    cgt::BSGS<N> rootBSGS;            // BSGS of C(alpha), built once per type; reused by setwiseStab each node
};

// One parallel work item: a representative first-factor-orbit (the orbit of the factor
// containing edge (0,1)) plus the centralizer subgroup that fixes that orbit. The subtree
// for this orbit only needs `stab` for symmetry-breaking (Stab is far smaller than C(alpha),
// which is what makes the dense types fast).
struct Task {
    Match m0;                   // the factor containing edge (0,1) (defines the first orbit)
    std::vector<int> stab;      // ENUMERABLE path: indices into Calpha of the gammas fixing this orbit
    std::vector<Perm> stabGens; // OVER-CAP path: a generating set of Stab_{C(alpha)}(this orbit)
};

// ---- shared progress/print state (one classification run) -------------------
std::atomic<long long> g_nodes{ 0 };  // total exact-cover nodes across all workers
std::atomic<long long> g_emits{ 0 };  // total distinct classes found so far (across workers)
std::mutex g_print_mtx;               // serializes the periodic progress line
std::chrono::steady_clock::time_point g_t0;         // run start time
std::chrono::steady_clock::time_point g_last_print; // time of last progress line
bool g_bprint = false;                // mirror of K20A2::m_bPrint (gate for ALL prints)
long long g_last_print_nodes = 0;     // g_nodes at the previous progress line (for the rate)
int  g_order = 0, g_typeIdx = 0, g_numTypes = 0;    // current order, current cycle type 1/N
char g_typeStr[80] = "";              // current cycle type string (e.g. "2 2 2 2 2 2 2 2 2")
std::atomic<int> g_reps_total{ 0 };   // representative first-orbit tasks for the current type
std::atomic<int> g_reps_done{ 0 };    // of those, how many subtrees have finished (within-type progress)

// ---- harvest mode (per-order class-count target) ----------------------------
// When g_target>0, the run STOPS as soon as g_harvest holds g_target distinct classes,
// skipping the (intractable) proof-of-exhaustion tail. g_harvest is the labeling-
// independent distinct set across ALL workers (each emit() inserts under g_harvest_mtx),
// so the count is exact even mid-run. Result is COMPLETE only if g_target == the true
// count (e.g. a known/published number); otherwise it is an explicit LOWER BOUND.
int g_target = 0;                     // 0 = run to completion; >0 = stop after this many distinct classes
std::atomic<bool> g_stop{ false };    // set once g_harvest reaches g_target; all workers/recursion bail
std::mutex g_harvest_mtx;             // guards g_harvest (locked only on a genuine new class -> rare)
std::set<std::string> g_harvest;      // distinct canonical keys seen so far this order (harvest mode only)

// One periodic progress line (every 30s) shared by cover()/coverGen(): shows the current
// order + cycle type, how many of the type's representative subtrees are done (the best
// completion-trend signal), node throughput since the last line, and classes found so far.
inline void progressTick(long long& nodes_flush, std::atomic<long long>& g_nodesRef) {
    if ((nodes_flush & 0x3FF) != 0) return;   // flush every 1024 nodes/worker (the 30s gate paces the actual print)
    g_nodesRef.fetch_add(nodes_flush, std::memory_order_relaxed);
    nodes_flush = 0;
    if (!g_bprint) return;
    auto now = std::chrono::steady_clock::now();
    std::lock_guard<std::mutex> lk(g_print_mtx);
    double sinceLast = std::chrono::duration<double>(now - g_last_print).count();
    if (sinceLast < 30.0) return;
    g_last_print = now;
    long long cur = g_nodesRef.load();
    double el = std::chrono::duration<double>(now - g_t0).count();
    double rate = sinceLast > 0 ? (cur - g_last_print_nodes) / sinceLast : 0;
    g_last_print_nodes = cur;
    printf("    [rep] o%d type %d/%d \"%s\" reps %d/%d  elapsed=%.0fs nodes=%lld (%.0f/s) classes=%lld\n",
           g_order, g_typeIdx, g_numTypes, g_typeStr, g_reps_done.load(), g_reps_total.load(),
           el, cur, rate, g_emits.load());
    fflush(stdout);
}

// ---- per-thread search worker ----------------------------------------------
struct RepWorker {
    const RepShared* sh = nullptr;    // shared read-only data for the current type

    std::vector<Match> chosen;        // factors committed so far on this search path
    std::vector<std::array<uint8_t, 2 * NHALF>> facEdges; // facEdges[f] = edge ids of factor f
    std::vector<int> facMin;          // facMin[f] = smallest edge id in factor f ("MEC")
    bool covered[N][N];               // covered[u][v] = edge {u,v} used by a chosen factor
    int  edgeFactor[NEDGES];          // edgeFactor[e] = chosen factor index covering edge e, or -1
    uint8_t path_end[NM][N];          // path_end[k][x] = endpoint of the alternating path through x
                                      //   in (chosen factor k) U (partial matching); Hamiltonian prune
    int lastAut = 0;                  // |Aut|*2 of the P1F most recently canonicalized

    std::set<std::string>      canon;     // distinct canonical keys found by THIS worker
    std::map<std::string, int> autOf;     // canonical key -> 2*|Aut| (this worker)
    std::vector<Match>* collectInto = nullptr;  // non-null => collect first-orbit tasks, don't recurse

    long long nodes = 0;              // local node counter (flushed to g_nodes periodically)
    long long nodes_flush = 0;        // local nodes not yet added to g_nodes

    void clearState() {               // reset search state for a fresh task
        chosen.clear(); facEdges.clear(); facMin.clear();
        memset(covered, 0, sizeof(covered));
        for (int e = 0; e < NEDGES; e++) edgeFactor[e] = -1;
    }

    // every edge of an alpha-orbit member must be free, so the alpha-images of any
    // edge we place must currently be uncovered (a cheap, sound early reject)
    inline bool imgUncovered(int u, int v) {
        int au = u, av = v;
        for (int j = 1; j < sh->order; j++) {
            au = sh->alpha[au]; av = sh->alpha[av];
            int x = au < av ? au : av, y = au < av ? av : au;
            if (x == u && y == v) break;          // edge-orbit closed
            if (covered[x][y]) return false;
        }
        return true;
    }

    // build the alpha-orbit {M, alphaM, ...} and validate it: every member must be
    // perfect with all chosen factors and the members must be mutually perfect.
    bool buildAndValidateOrbit(const Match& M, std::vector<Match>& orbit) {
        orbit.clear(); orbit.push_back(M);
        Match cur = M;
        for (int s = 1; s < sh->order; s++) { cur = applyAlpha(cur, sh->alpha); if (cur == M) break; orbit.push_back(cur); }
        for (const auto& F : orbit) for (const auto& C : chosen) if (!is_perfect(F, C)) return false;
        for (size_t x = 0; x < orbit.size(); x++) for (size_t y = x + 1; y < orbit.size(); y++) if (!is_perfect(orbit[x], orbit[y])) return false;
        return true;
    }

    void commitOrbit(const std::vector<Match>& orbit) {   // append the orbit's factors
        for (const auto& F : orbit) {
            int fidx = (int)chosen.size();
            chosen.push_back(F);
            std::array<uint8_t, 2 * NHALF> el;   // this factor's edge ids
            int ec = 0;               // edge count written into el
            int mn = INF;             // min edge id of this factor
            for (int u = 0; u < N; u++) { int v = F[u]; if (u < v) { covered[u][v] = true; int e = eidT[u][v]; el[ec++] = (uint8_t)e; edgeFactor[e] = fidx; if (e < mn) mn = e; } }
            facEdges.push_back(el); facMin.push_back(mn);
        }
    }

    void rollbackTo(size_t base) {    // undo factors added after index `base`
        for (size_t k = base; k < chosen.size(); k++) {
            const Match& F = chosen[k];
            for (int u = 0; u < N; u++) { int v = F[u]; if (u < v) { covered[u][v] = false; edgeFactor[eidT[u][v]] = -1; } }
        }
        chosen.resize(base); facEdges.resize(base); facMin.resize(base);
    }

    // labeling-independent canonical form: lexicographically-smallest serialization
    // over all (ordered factor pair, direction, anchor) Hamiltonian-cycle frames.
    // Side effect: lastAut = number of frames achieving the minimum = 2*|Aut|.
    std::string canonKey(const std::vector<Match>& f) {
        std::string best;             // current lexicographically smallest serialization
        bool first = true;            // have we seen any valid frame yet?
        int bestCount = 0;            // how many frames tie for `best`
        for (int i = 0; i < NM; i++) for (int j = 0; j < NM; j++) {
            if (i == j) continue;
            const Match& A = f[i]; const Match& B = f[j];   // the factor pair forming the basis cycle
            for (int dir = 0; dir < 2; dir++) {
                const Match& X = dir ? B : A; const Match& Y = dir ? A : B;  // traversal order
                for (int a = 0; a < N; a++) {               // anchor (cycle start vertex)
                    uint8_t c[N];      // the Hamiltonian cycle as a vertex sequence
                    bool seen[N]; memset(seen, 0, sizeof(seen));  // seen[v] during cycle build
                    int pos = 0, cur = a; bool ok = true;
                    for (int t = 0; t < NHALF; t++) {
                        if (seen[cur]) { ok = false; break; } seen[cur] = true; c[pos++] = (uint8_t)cur; cur = X[cur];
                        if (seen[cur]) { ok = false; break; } seen[cur] = true; c[pos++] = (uint8_t)cur; cur = Y[cur];
                    }
                    if (!ok || cur != a || pos != N) continue;   // pair didn't form a single N-cycle here
                    uint8_t p[N]; for (int t = 0; t < N; t++) p[c[t]] = (uint8_t)t;  // relabel cycle -> 0..N-1
                    std::vector<Match> rel(NM);   // the relabeled factor set
                    for (int k = 0; k < NM; k++) { Match g; for (int u = 0; u < N; u++) g[p[u]] = p[f[k][u]]; rel[k] = g; }
                    std::sort(rel.begin(), rel.end());        // canonicalize factor order
                    std::string s; s.reserve(NM * N);         // serialized relabeled P1F
                    for (auto& g : rel) for (int u = 0; u < N; u++) s.push_back((char)g[u]);
                    if (first || s < best) { best = s; first = false; bestCount = 1; }
                    else if (s == best) bestCount++;
                }
            }
        }
        lastAut = bestCount;
        return best;
    }

    void emit() {
        std::string key = canonKey(chosen);          // labeling-independent class key
        if (canon.insert(key).second) {              // first time this worker sees the class
            autOf[key] = lastAut;
            g_emits.fetch_add(1, std::memory_order_relaxed);  // (approx; merged-dedup is exact)
            if (g_target > 0) {                      // harvest mode: track the EXACT cross-worker distinct count
                std::lock_guard<std::mutex> lk(g_harvest_mtx);
                g_harvest.insert(key);
                if ((int)g_harvest.size() >= g_target) g_stop.store(true, std::memory_order_relaxed);  // target hit -> bail everywhere
            }
        }
    }

    void processM(const Match& M, const std::vector<int>& active) {
        std::vector<Match> orbit;     // the alpha-orbit {M, alphaM, ...} (a block of factors)
        if (!buildAndValidateOrbit(M, orbit)) return;
        size_t base = chosen.size();  // roll-back point
        commitOrbit(orbit);
        cover(active);                // recurse on the rest of the cover
        rollbackTo(base);
    }

    // generate every perfect matching of the still-uncovered graph that contains
    // the forced first edge, pairing the lowest unused vertex with each legal
    // partner; prune partials that cannot stay Hamiltonian with a chosen factor.
    void genM(Match& M, bool* used, int edge_count) {
        int u = 0; while (u < N && used[u]) u++;   // lowest still-unmatched vertex
        if (u == N) {                              // matching complete -> record the orbit member
            std::vector<Match> orb;                // (cover always drives genM in collection mode)
            if (buildAndValidateOrbit(M, orb)) collectInto->push_back(M);
            return;
        }
        int nchosen = (int)chosen.size();
        used[u] = true;
        for (int v = u + 1; v < N; v++) {
            if (used[v] || covered[u][v]) continue;        // v taken, or edge already used
            if (!imgUncovered(u, v)) continue;             // alpha-images of {u,v} must be free
            bool ok = true;                                // Hamiltonian-with-chosen pre-check
            if (edge_count < NHALF - 1) for (int k = 0; k < nchosen; k++) if (path_end[k][u] == v) { ok = false; break; }
            if (!ok) continue;
            uint8_t epu[NM], epv[NM];   // path endpoints touched by adding edge {u,v}, per chosen factor
            for (int k = 0; k < nchosen; k++) { epu[k] = path_end[k][u]; epv[k] = path_end[k][v]; path_end[k][epu[k]] = epv[k]; path_end[k][epv[k]] = epu[k]; }
            M[u] = (uint8_t)v; M[v] = (uint8_t)u; used[v] = true;
            genM(M, used, edge_count + 1);
            used[v] = false;
            for (int k = 0; k < nchosen; k++) { path_end[k][epu[k]] = (uint8_t)u; path_end[k][epv[k]] = (uint8_t)v; }  // undo
        }
        used[u] = false;
    }

    // Orderly-generation cover step. `active` = the current stabilizer subgroup
    // G <= C(alpha) (as indices into sh->Calpha) that fixes the chosen factor set.
    // We collect every candidate factor-orbit covering the lowest uncovered edge,
    // keep ONE representative per G-orbit of candidates, and recurse each chosen
    // orbit O with its sub-stabilizer Stab_G(O). Sound for class counting: any
    // g in G maps {solutions extending P+O1} bijectively onto {solutions extending
    // P+O2} whenever g*O1 = O2, so exploring one representative finds every class
    // the others would (the global canonical key dedups any residual repeats).
    void cover(const std::vector<int>& active) {
        // local node counter; periodically flush to the shared atomic and (if enabled)
        // print one progress line per 30s. No ETA, no percentage — elapsed time only.
        nodes++; nodes_flush++;
        progressTick(nodes_flush, g_nodes);
        if (g_stop.load(std::memory_order_relaxed)) return;   // harvest target reached -> unwind the search

        int u0 = -1, v0 = -1;         // endpoints of the lowest uncovered edge (the cover anchor)
        for (int u = 0; u < N && u0 < 0; u++) for (int v = u + 1; v < N; v++) if (!covered[u][v]) { u0 = u; v0 = v; break; }
        if (u0 < 0) { if ((int)chosen.size() == NM) emit(); return; }   // all edges covered -> a P1F
        if (!imgUncovered(u0, v0)) return;            // no alpha-orbit can cover this edge

        // 1. collect every candidate matching whose orbit covers the forced anchor edge.
        int nchosen = (int)chosen.size();
        for (int k = 0; k < nchosen; k++) for (int x = 0; x < N; x++) path_end[k][x] = chosen[k][x];  // reset path_end
        Match M; for (int i = 0; i < N; i++) M[i] = 0xFF;     // matching being built (0xFF = unset)
        bool used[N]; memset(used, 0, sizeof(used));          // used[v] = v already matched in M
        M[u0] = (uint8_t)v0; M[v0] = (uint8_t)u0; used[u0] = used[v0] = true;   // force the anchor edge
        for (int k = 0; k < nchosen; k++) { uint8_t eu = path_end[k][u0], ev = path_end[k][v0]; path_end[k][eu] = ev; path_end[k][ev] = eu; }
        std::vector<Match> cands;
        collectInto = &cands; genM(M, used, 1); collectInto = nullptr;
        if (cands.empty()) return;

        // 2. collapse to one member per distinct alpha-orbit (genM may emit several
        //    matchings that lie in the same orbit), keyed by the orbit serialization.
        std::unordered_map<std::string, int> keyIdx;   // orbit key -> representative index
        std::vector<std::string> keys;                 // keys[i] = orbit key of reps[i]
        std::vector<Match> reps;                        // one matching per distinct orbit
        for (auto& C : cands) {
            std::string k = orbitKey(C, sh->alpha, sh->order);
            if (keyIdx.emplace(k, (int)reps.size()).second) { keys.push_back(std::move(k)); reps.push_back(C); }
        }

        // 3a. no usable symmetry (G empty: over-cap/sparse types) -> explore every orbit;
        //     the global canonical key removes any labeling duplicates downstream.
        if (active.empty()) {
            for (auto& R : reps) processM(R, active);
            return;
        }

        // 3b. dedup candidate orbits by G and recurse each representative with Stab_G(O).
        std::vector<char> marked(reps.size(), 0);       // marked[j] = orbit j is in an already-handled G-orbit
        for (size_t i = 0; i < reps.size(); i++) {
            if (marked[i]) continue;
            std::vector<int> stab;                      // gammas in G fixing this orbit = Stab_G(O)
            for (int gi : active) {
                std::string img = orbitKey(applyAlpha(reps[i], sh->Calpha[gi].data()), sh->alpha, sh->order);
                if (img == keys[i]) { stab.push_back(gi); continue; }   // gamma fixes the orbit
                auto it = keyIdx.find(img);             // gamma maps it onto another candidate -> same G-orbit
                if (it != keyIdx.end()) marked[it->second] = 1;
            }
            processM(reps[i], stab);                     // commit orbit, recurse with the sub-stabilizer
        }
    }

    // Cover step for OVER-CAP types, where C(alpha) is too large to enumerate (sh->Calpha
    // empty). The candidate-collection is identical to cover(); the symmetry reduction
    // recomputes the EXACT Stab_C(P) of the current partial cover from scratch each node
    // (setwiseStab over C(alpha)'s generators), then dedups candidates by the stabilizer of
    // the anchor edge within it. Computing Stab_C(P) exactly per node (rather than carrying an
    // under-approximating generator chain) is what makes the dense types -- order-2 2^9 in
    // particular -- terminate: the dedup is full strength at every depth. Sound for counting
    // (dedup is by a genuine subgroup; global canonKey is the final net).
    void coverGen() {
        nodes++; nodes_flush++;
        progressTick(nodes_flush, g_nodes);
        if (g_stop.load(std::memory_order_relaxed)) return;   // harvest target reached -> unwind the search

        int u0 = -1, v0 = -1;         // lowest uncovered edge (the cover anchor)
        for (int u = 0; u < N && u0 < 0; u++) for (int v = u + 1; v < N; v++) if (!covered[u][v]) { u0 = u; v0 = v; break; }
        if (u0 < 0) { if ((int)chosen.size() == NM) emit(); return; }
        if (!imgUncovered(u0, v0)) return;

        // 1. collect every candidate matching whose orbit covers the forced anchor edge
        int nchosen = (int)chosen.size();
        for (int k = 0; k < nchosen; k++) for (int x = 0; x < N; x++) path_end[k][x] = chosen[k][x];
        Match M; for (int i = 0; i < N; i++) M[i] = 0xFF;
        bool used[N]; memset(used, 0, sizeof(used));
        M[u0] = (uint8_t)v0; M[v0] = (uint8_t)u0; used[u0] = used[v0] = true;
        for (int k = 0; k < nchosen; k++) { uint8_t eu = path_end[k][u0], ev = path_end[k][v0]; path_end[k][eu] = ev; path_end[k][ev] = eu; }
        std::vector<Match> cands;
        collectInto = &cands; genM(M, used, 1); collectInto = nullptr;
        if (cands.empty()) return;

        // 2. collapse to one member per distinct alpha-orbit
        std::unordered_map<std::string, int> keyIdx;
        std::vector<std::string> keys;
        std::vector<Match> reps;
        for (auto& C : cands) {
            std::string k = orbitKey(C, sh->alpha, sh->order);
            if (keyIdx.emplace(k, (int)reps.size()).second) { keys.push_back(std::move(k)); reps.push_back(C); }
        }

        // 3. exact stabilizer of the current partial cover within C(alpha)
        std::vector<Perm> Gp = cgt::setwiseStab<N>(sh->rootBSGS, chosen, covered);

        // 3a. trivial stabilizer -> no usable symmetry anywhere below, explore every orbit
        if (Gp.empty()) { for (auto& R : reps) processMGen(R); return; }

        // 3b. dedup candidates by the stabilizer of the anchor edge within Stab_C(P). Its
        //     generators fix (u0,v0), so the candidate-restricted BFS stays inside the
        //     candidate set and recovers the full Stab_C(P)(edge)-orbit dedup at this node.
        std::vector<Perm> hgens = cgt::edgeStabGens<N>(Gp, u0, v0);
        std::vector<int> repIdx;
        std::vector<std::vector<Perm>> stabs;            // (Schreier stabs unused: Stab is recomputed per node)
        schreierDedup(reps, keys, keyIdx, hgens, sh->order, sh->alpha, repIdx, stabs);
        for (size_t r = 0; r < repIdx.size(); r++) processMGen(reps[repIdx[r]]);
    }

    void processMGen(const Match& M) {                    // over-cap sibling of processM
        std::vector<Match> orbit;
        if (!buildAndValidateOrbit(M, orbit)) return;
        size_t base = chosen.size();
        commitOrbit(orbit);
        coverGen();
        rollbackTo(base);
    }

    // One-level expansion of the current partial cover (state already committed): collect the
    // DEDUPED child partial covers (each = current chosen + a child orbit) into outChildren
    // instead of recursing; complete covers are emitted here. Mirrors coverGen's dedup so the
    // children are exactly the subtrees coverGen would explore. Used to build a fine task
    // frontier so the dense over-cap types (few root reps) can use all worker threads.
    void splitNode(std::vector<std::vector<Match>>& outChildren) {
        nodes++; nodes_flush++;
        progressTick(nodes_flush, g_nodes);
        if (g_stop.load(std::memory_order_relaxed)) return;   // harvest target reached -> unwind the search
        int u0 = -1, v0 = -1;
        for (int u = 0; u < N && u0 < 0; u++) for (int v = u + 1; v < N; v++) if (!covered[u][v]) { u0 = u; v0 = v; break; }
        if (u0 < 0) { if ((int)chosen.size() == NM) emit(); return; }     // leaf -> emit, no children
        if (!imgUncovered(u0, v0)) return;                                // dead end
        int nchosen = (int)chosen.size();
        for (int k = 0; k < nchosen; k++) for (int x = 0; x < N; x++) path_end[k][x] = chosen[k][x];
        Match M; for (int i = 0; i < N; i++) M[i] = 0xFF;
        bool used[N]; memset(used, 0, sizeof(used));
        M[u0] = (uint8_t)v0; M[v0] = (uint8_t)u0; used[u0] = used[v0] = true;
        for (int k = 0; k < nchosen; k++) { uint8_t eu = path_end[k][u0], ev = path_end[k][v0]; path_end[k][eu] = ev; path_end[k][ev] = eu; }
        std::vector<Match> cands;
        collectInto = &cands; genM(M, used, 1); collectInto = nullptr;
        if (cands.empty()) return;
        std::unordered_map<std::string, int> keyIdx; std::vector<std::string> keys; std::vector<Match> reps;
        for (auto& C : cands) { std::string k = orbitKey(C, sh->alpha, sh->order); if (keyIdx.emplace(k, (int)reps.size()).second) { keys.push_back(std::move(k)); reps.push_back(C); } }
        std::vector<Perm> Gp = cgt::setwiseStab<N>(sh->rootBSGS, chosen, covered);
        std::vector<int> repIdx;
        if (Gp.empty()) { for (size_t i = 0; i < reps.size(); i++) repIdx.push_back((int)i); }
        else { std::vector<Perm> hgens = cgt::edgeStabGens<N>(Gp, u0, v0); std::vector<std::vector<Perm>> stabs; schreierDedup(reps, keys, keyIdx, hgens, sh->order, sh->alpha, repIdx, stabs); }
        for (int ri : repIdx) {
            std::vector<Match> orbit;
            if (!buildAndValidateOrbit(reps[ri], orbit)) continue;
            std::vector<Match> child = chosen; for (auto& F : orbit) child.push_back(F);
            outChildren.push_back(std::move(child));
        }
    }


    // collect the valid first-factor-orbits (parallel tasks) into `tasks`, by running
    // the root-level matching generation once. Each task is the matching M0 that
    // contains the lowest edge (0,1); its orbit is the first block of the cover.
    void collectTasks(std::vector<Match>& tasks) {
        clearState();
        // anchor edge is (0,1) on the empty board
        Match M; for (int i = 0; i < N; i++) M[i] = 0xFF;
        bool used[N]; memset(used, 0, sizeof(used));
        M[0] = 1; M[1] = 0; used[0] = used[1] = true;
        collectInto = &tasks;
        genM(M, used, 1);
        collectInto = nullptr;
    }

    // run one task: commit its first-orbit, then exact-cover the rest. Dispatch on whether
    // C(alpha) was enumerated: enumerable types use the index-based cover() with the exact
    // stabilizer; over-cap types use the generator-based coverGen() with Schreier stabilizers.
    void runTask(const Task& tk) {
        clearState();
        std::vector<Match> orbit;
        if (!buildAndValidateOrbit(tk.m0, orbit)) return;   // safety (tasks are pre-validated)
        commitOrbit(orbit);
        if (!sh->Calpha.empty()) cover(tk.stab);
        else coverGen();
    }
};

// ---- cycle-type enumeration -------------------------------------------------
int gcd_(int a, int b) { while (b) { int t = a % b; a = b; b = t; } return a; }
int lcm_(int a, int b) { return a / gcd_(a, b) * b; }

void enumTypes(int order, std::vector<std::vector<int>>& out) {
    std::vector<int> divs;            // divisors of order that are >= 2 (cycle lengths > 1)
    for (int d = 2; d <= order; d++) if (order % d == 0) divs.push_back(d);
    std::vector<int> cur;             // parts chosen so far (the partial type)
    std::function<void(int, int)> rec = [&](int idx, int rem) {
        if (idx == (int)divs.size()) {
            int L = 1; for (int p : cur) L = lcm_(L, p);     // order of this type
            if (L == order) { std::vector<int> t = cur; for (int i = 0; i < rem; i++) t.push_back(1); out.push_back(t); }
            return;
        }
        int d = divs[idx]; int saved = (int)cur.size();
        for (int cnt = 0; cnt * d <= rem; cnt++) { rec(idx + 1, rem - cnt * d); cur.push_back(d); }
        cur.resize(saved);
    };
    rec(0, N);
}

}  // anonymous namespace

// =============================================================================
// Private K20A2 entry point (declared in k20a2.h, called only from
// runExhaustiveSearch). Enumerates all K20 P1Fs with an automorphism of `order`,
// deduplicates them, forwards each distinct class to resultCallback, and prints a
// summary (only when m_bPrint). Worker count = kThreads.
// =============================================================================
void K20A2::runRepresentativeMethod(int order, int target) {
    initEdges();
    // Centralizer-enumeration cap, overridable at runtime via env REP_SYMCAP (e.g. raise it
    // to pull a borderline dense type like order-3 3^5.1^5 (|C|=3.5M) onto the fast path).
    long long symCap = K20_REP_SYMCAP;
    if (const char* env = std::getenv("REP_SYMCAP")) symCap = atoll(env);
    int nThreads = (kThreads > 0) ? kThreads : 1;   // worker count from the standard KThreads knob

    g_bprint = m_bPrint;
    g_nodes.store(0); g_emits.store(0); g_last_print_nodes = 0;
    g_target = target; g_stop.store(false); g_harvest.clear();   // harvest mode (target>0 -> stop after `target` classes)
    g_t0 = std::chrono::steady_clock::now();
    g_last_print = g_t0;
    g_order = order;

    std::vector<std::vector<int>> types;            // all cycle types of this order
    enumTypes(order, types);
    g_numTypes = (int)types.size();
    int g_typeCounter = 0;                          // 1-based index of the type currently processed
    if (m_bPrint) {
        printf("[K20-REP] order=%d : %zu cycle types, kThreads=%d, symCap=%lld\n", order, types.size(), nThreads, symCap);
        if (target > 0) printf("[K20-REP] HARVEST MODE: stop after %d distinct classes (COMPLETE only if %d == true count, else a LOWER BOUND)\n", target, target);
        fflush(stdout);
    }

    std::set<std::string>      gcanon;   // global distinct canonical keys (across all types)
    std::map<std::string, int> gautOf;   // global canonical key -> 2*|Aut|

    static const bool only_fpf = (getenv("REP_FPF") != nullptr);   // [DIAG] isolate fixed-point-free type(s)
    for (auto& type : types) {
        if (g_stop.load(std::memory_order_relaxed)) break;   // harvest target already reached -> skip remaining types
        if (only_fpf) { bool hasFix = false; for (int p : type) if (p == 1) hasFix = true; if (hasFix) continue; }
        RepShared sh;                    // shared read-only data for this type
        sh.order = order;
        for (int i = 0; i < N; i++) sh.alpha[i] = (uint8_t)i;
        int base = 0;                    // first vertex of the next cycle block
        char tstr[64]; int tp = 0;       // human-readable type string, e.g. "4 4 4 4 2 "
        for (int p : type) {
            if (p > 1) { for (int i = 0; i < p; i++) sh.alpha[base + i] = (uint8_t)(base + (i + 1) % p); tp += snprintf(tstr + tp, sizeof(tstr) - tp, "%d ", p); }
            base += p;
        }
        buildCentralizer(sh.alpha, sh.Calpha, symCap);   // empty => run this type without symmetry-breaking
        sh.Cinv.resize(sh.Calpha.size());
        for (size_t i = 0; i < sh.Calpha.size(); i++) for (int x = 0; x < N; x++) sh.Cinv[i][sh.Calpha[i][x]] = (uint8_t)x;
        sh.rootActive.resize(sh.Calpha.size());
        for (size_t i = 0; i < sh.Calpha.size(); i++) sh.rootActive[i] = (int)i;

        // 1. collect the raw first-orbit tasks (single-threaded root enumeration)
        auto tcollect = std::chrono::steady_clock::now();   // [diag] start of collect+dedup phase
        std::vector<Match> raw;
        { RepWorker collector; collector.sh = &sh; collector.collectTasks(raw); }

        // 2. reduce tasks by C(alpha): keep ONE representative per C(alpha)-orbit. Each orbit is
        //    handled once, so the cost is O(#reps * <group walk>), not O(#rawTasks * |C|).
        // First collapse raw tasks to one matching per distinct alpha-orbit (collectTasks may emit
        // several matchings from the same orbit), keyed by the alpha-orbit serialization.
        std::unordered_map<std::string, int> orbIdx;   // alpha-orbit key -> unique-task index
        std::vector<Match> uniq;                        // one matching per distinct alpha-orbit
        std::vector<std::string> uniqKey;               // its alpha-orbit key
        for (auto& m : raw) {
            std::string k = orbitKey(m, sh.alpha, order);
            if (orbIdx.emplace(k, (int)uniq.size()).second) { uniq.push_back(m); uniqKey.push_back(std::move(k)); }
        }

        std::vector<Task> reps;
        std::vector<char> seen(uniq.size(), 0);         // seen[j] = unique-task j already covered
        if (!sh.Calpha.empty()) {
            // C(alpha) enumerated (dense types): mark each C(alpha)-orbit by applying every gamma,
            // and record the orbit's STABILIZER so its (large, productive) subtree is symmetry-broken.
            for (size_t i = 0; i < uniq.size(); i++) {
                if (seen[i]) continue;
                seen[i] = 1;
                std::vector<int> stab;                  // gammas fixing this orbit (setwise)
                for (int gi = 0; gi < (int)sh.Calpha.size(); gi++) {
                    std::string img = orbitKey(applyAlpha(uniq[i], sh.Calpha[gi].data()), sh.alpha, order);
                    if (img == uniqKey[i]) stab.push_back(gi);     // gamma fixes the orbit
                    auto it = orbIdx.find(img);
                    if (it != orbIdx.end()) seen[it->second] = 1;  // mark this C(alpha)-orbit member
                }
                reps.push_back(Task{ uniq[i], std::move(stab) });
            }
        } else {
            // C(alpha) too large to enumerate (sparse types, AND every order-2 type, whose fixed
            // points alone make |C| huge): represent the group by a GENERATING SET and, for each
            // first-orbit task, compute Schreier generators of its setwise stabilizer in C(alpha)
            // via a candidate-restricted BFS. These flow into coverGen so the over-cap subtrees get
            // the same per-level orbit dedup the enumerable types get from explicit elements.
            buildGenerators(sh.alpha, sh.rootCgens);              // seed for exact Stab_C(P) in coverGen
            sh.rootBSGS = cgt::buildBSGS<N>(sh.rootCgens);        // built once per type, reused per node
            std::vector<Perm> hgens = cgt::edgeStabGens<N>(sh.rootCgens, 0, 1);  // collectTasks forces the anchor edge (0,1)
            std::vector<int> repIdx; std::vector<std::vector<Perm>> stabs;
            schreierDedup(uniq, uniqKey, orbIdx, hgens, order, sh.alpha, repIdx, stabs);
            for (size_t r = 0; r < repIdx.size(); r++)
                reps.push_back(Task{ uniq[repIdx[r]], {}, std::move(stabs[r]) });
        }

        // 3. fan the representative tasks out across kThreads workers; each subtree breaks
        //    symmetry with only its orbit's stabilizer.
        size_t before = gcanon.size();   // classes known before this type
        auto ts = std::chrono::steady_clock::now();
        if (m_bPrint) {  // [diag] dedup finished; show task reduction and stabilizer size before the cover
            size_t stab0 = reps.empty() ? 0 : (sh.Calpha.empty() ? reps[0].stabGens.size() : reps[0].stab.size());
            printf("[K20-REP]  type %-14s rawTasks %zu uniq %zu reps %zu stab0 %zu (dedup %.1fs) -> cover...\n",
                   tstr, raw.size(), uniq.size(), reps.size(), stab0,
                   std::chrono::duration<double>(ts - tcollect).count());
            fflush(stdout);
        }
        // publish current-type progress state for the 30s [rep] line
        g_typeIdx = ++g_typeCounter; g_reps_done.store(0); g_reps_total.store(0);
        snprintf(g_typeStr, sizeof(g_typeStr), "%s", tstr);
        std::vector<RepWorker> workers((size_t)nThreads);
        std::vector<std::thread> pool;
        // Worker-shared state. MUST outlive the threads (joined after the if/else below), so it is
        // declared here in the same scope as `pool`/the join -- NOT inside the branches, where it
        // would be destroyed before the threads referencing it run.
        std::deque<std::vector<Match>> queue;             // OVER-CAP: shared work queue of partial covers
        std::mutex qmtx;                                  // OVER-CAP: guards `queue`
        std::atomic<int> active{ 0 };                     // OVER-CAP: tasks currently being expanded
        std::atomic<size_t> next{ 0 };                    // ENUMERABLE: next representative-task index
        if (sh.Calpha.empty()) {
            // OVER-CAP: a few root reps (order-2 2^9 has 5) would idle most cores with a per-rep
            // fan-out. Use a SHARED WORK QUEUE of partial covers: each worker pops one, expands
            // it ONE deduped level (splitNode -> children + leaf emits), and pushes the children
            // back. This spreads the per-node stabilizer work (shallow AND deep) across all
            // threads with load balancing. LIFO (push/pop back) keeps it DFS-like so the queue
            // stays bounded. g_reps_total/done track tasks-created vs tasks-completed.
            const size_t QCAP = 200000;                   // queue cap (~60MB); over it, workers drain subtrees locally
            { RepWorker seed; seed.sh = &sh; for (auto& tk : reps) { std::vector<Match> orb; if (seed.buildAndValidateOrbit(tk.m0, orb)) queue.push_back(std::move(orb)); } }
            g_reps_total.store((int)queue.size());
            for (int t = 0; t < nThreads; t++) {
                workers[t].sh = &sh;
                pool.emplace_back([&, t]() {
                    RepWorker& w = workers[t];
                    for (;;) {
                        if (g_stop.load(std::memory_order_relaxed)) break;   // harvest target reached
                        std::vector<Match> task; bool have = false;
                        {
                            std::unique_lock<std::mutex> lk(qmtx);
                            if (!queue.empty()) { task = std::move(queue.back()); queue.pop_back(); active.fetch_add(1); have = true; }
                            else if (active.load() == 0) break;     // no work and nobody producing -> done
                        }
                        if (!have) { std::this_thread::yield(); continue; }
                        w.clearState(); w.commitOrbit(task);
                        std::vector<std::vector<Match>> kids; w.splitNode(kids);
                        bool pushed = false;
                        {
                            std::unique_lock<std::mutex> lk(qmtx);
                            if (queue.size() < QCAP) { for (auto& k : kids) queue.push_back(std::move(k)); pushed = true; }
                        }
                        if (pushed) g_reps_total.fetch_add((int)kids.size(), std::memory_order_relaxed);
                        else for (auto& k : kids) { w.clearState(); w.commitOrbit(k); w.coverGen(); }  // queue full -> drain this subtree locally (bounds memory)
                        g_reps_done.fetch_add(1, std::memory_order_relaxed);
                        active.fetch_sub(1);
                    }
                });
            }
        } else {
            // ENUMERABLE: many root reps already -> simple per-rep fan-out via an atomic index.
            g_reps_total.store((int)reps.size());
            for (int t = 0; t < nThreads; t++) {
                workers[t].sh = &sh;
                pool.emplace_back([&, t]() {
                    RepWorker& w = workers[t];
                    for (;;) { if (g_stop.load(std::memory_order_relaxed)) break; size_t i = next.fetch_add(1); if (i >= reps.size()) break; w.runTask(reps[i]); g_reps_done.fetch_add(1, std::memory_order_relaxed); }
                });
            }
        }
        for (auto& th : pool) th.join();

        // 4. merge each worker's distinct classes into the global dedup
        for (auto& w : workers)
            for (auto& kv : w.autOf)
                if (gcanon.insert(kv.first).second) gautOf[kv.first] = kv.second;

        if (m_bPrint) {
            double tsec = std::chrono::duration<double>(std::chrono::steady_clock::now() - ts).count();
            printf("[K20-REP]  type %-14s |C|=%-8zu rawTasks %-7zu reps %-7zu new +%zu (cum %zu) %.1fs\n",
                   tstr, sh.Calpha.size(), raw.size(), reps.size(), gcanon.size() - before, gcanon.size(), tsec);
            fflush(stdout);
        }
    }

    // emit one representative P1F per distinct class, in "src" form, into the pipeline
    unsigned char src[NM * N];        // scratch buffer for one P1F in pair form
    for (const std::string& key : gcanon) {
        for (int k = 0; k < NM; k++) {
            uint8_t adj[N];           // factor k of this class, in adjacency form (decoded from key)
            for (int u = 0; u < N; u++) adj[u] = (uint8_t)key[k * N + u];
            adj_to_src(adj, src + k * N);
        }
        if (resultCallback) resultCallback(cbClass, src, 0, 1, 2);
    }

    if (m_bPrint) {
        std::map<int, int> hist;      // hist[|Aut|] = number of classes with that group order
        for (auto& kv : gautOf) hist[kv.second / 2]++;   // stored value is 2*|Aut|
        double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - g_t0).count();
        printf("[K20-REP] |Aut| histogram:");
        for (auto& h : hist) printf(" %d:%d", h.first, h.second);
        bool harvested = (g_target > 0 && g_stop.load(std::memory_order_relaxed));
        if (harvested)
            printf("\n[K20-REP] order=%d HARVESTED %zu classes (target %d reached -- NOT proven complete) (%.1fs)\n", order, gcanon.size(), g_target, secs);
        else
            printf("\n[K20-REP] order=%d TOTAL DISTINCT CLASSES = %zu  (%.1fs)\n", order, gcanon.size(), secs);
        fflush(stdout);
    }
}
