// =============================================================================
// cgtEngine.h -- small computational-group-theory engine on N points, shared by the
// automorphism-representative classifiers (k18a2rep / k16A2rep / k20a2rep). Templated on
// N so each classifier instantiates it at its own vertex count.
//
// Contents (namespace cgt):
//   * permutation helpers: idPerm, composePerm, invPerm, applyPerm
//   * orbitTransversal     : orbit of a point with coset transversal
//   * BSGS + buildBSGS      : Schreier-Sims base & strong generating set
//   * reduceGens            : bound a generating set to a strong generating set
//   * edgeStabGens          : pointwise stabilizer of an edge's two endpoints
//   * setwiseStab           : exact setwise stabilizer of a set of factors (perfect
//                             matchings), via property-pruned backtrack with K-minimality
//
// All validated standalone: BSGS orders match prod_l l^{m_l} m_l!; setwiseStab matches
// brute-force enumeration on enumerable groups (small stabilizers AND whole-group cases).
// =============================================================================
#pragma once
#include <cstdint>
#include <array>
#include <vector>
#include <set>
#include <functional>
#include <algorithm>

namespace cgt {

template<int N> using Perm = std::array<uint8_t, N>;

template<int N> inline Perm<N> idPerm() { Perm<N> p; for (int i = 0; i < N; i++) p[i] = (uint8_t)i; return p; }
template<int N> inline Perm<N> composePerm(const Perm<N>& a, const Perm<N>& b) { Perm<N> r; for (int x = 0; x < N; x++) r[x] = a[b[x]]; return r; } // (a o b)(x)=a(b(x))
template<int N> inline Perm<N> invPerm(const Perm<N>& a) { Perm<N> r; for (int x = 0; x < N; x++) r[a[x]] = (uint8_t)x; return r; }
// image of a factor m (adjacency form) under permutation g: (g m)[g u] = g(m[u])
template<int N> inline Perm<N> applyPerm(const Perm<N>& m, const uint8_t* g) { Perm<N> r; for (int u = 0; u < N; u++) r[g[u]] = g[m[u]]; return r; }

// orbit of point b under <gens>, with a transversal rep[p] (a perm mapping b -> p) and a
// visited mask. The orbit on N points has size <= N, so this is always cheap.
template<int N> inline void orbitTransversal(int b, const std::vector<Perm<N>>& gens,
                                             std::array<Perm<N>, N>& rep, std::array<char, N>& vis, std::vector<int>& orb) {
    for (int p = 0; p < N; p++) vis[p] = 0;
    orb.clear();
    rep[b] = idPerm<N>(); vis[b] = 1; orb.push_back(b);
    for (size_t h = 0; h < orb.size(); h++) {
        int x = orb[h];
        for (const auto& g : gens) { int y = g[x]; if (!vis[y]) { vis[y] = 1; rep[y] = composePerm<N>(g, rep[x]); orb.push_back(y); } }
    }
}

// Base and strong generating set for a permutation group on N points (deterministic
// Schreier-Sims). Used to bound generating sets (reduceGens) and as the search tree for
// setwiseStab. |G| = product of the per-level orbit sizes.
template<int N> struct BSGS {
    std::vector<int> base;                     // base points
    std::vector<std::vector<Perm<N>>> S;       // S[i] = strong gens fixing base[0..i-1]
    std::vector<std::vector<int>> orbit;       // orbit[i] = orbit of base[i] under S[i]
    std::vector<std::array<Perm<N>, N>> urep;  // urep[i][p] = transversal element base[i] -> p
    std::vector<std::array<char, N>> inorb;    // inorb[i][p] = (p in orbit[i])
};

template<int N> inline int bsgsStrip(const BSGS<N>& B, const Perm<N>& g, Perm<N>& residue) {
    residue = g;
    for (size_t i = 0; i < B.base.size(); i++) {
        int p = residue[B.base[i]];
        if (!B.inorb[i][p]) return (int)i;
        residue = composePerm<N>(invPerm<N>(B.urep[i][p]), residue);
    }
    return (int)B.base.size();
}

template<int N> inline BSGS<N> buildBSGS(const std::vector<Perm<N>>& gens0) {
    BSGS<N> B;
    auto recompute = [&](size_t i) {
        std::array<Perm<N>, N> rep; std::array<char, N> vis; std::vector<int> orb;
        orbitTransversal<N>(B.base[i], B.S[i], rep, vis, orb);
        B.orbit[i] = orb; B.urep[i] = rep; for (int p = 0; p < N; p++) B.inorb[i][p] = vis[p];
    };
    auto ensureLevel = [&](int mp) { B.base.push_back(mp); B.S.push_back({}); B.orbit.push_back({}); B.urep.push_back({}); B.inorb.push_back({}); };
    for (const auto& g : gens0) {
        if (B.base.empty()) { int mp = -1; for (int x = 0; x < N; x++) if (g[x] != x) { mp = x; break; } if (mp < 0) continue; ensureLevel(mp); }
        B.S[0].push_back(g);
    }
    if (B.base.empty()) return B;                 // trivial group
    for (size_t i = 0; i < B.base.size(); i++) recompute(i);

    size_t i = B.base.size();
    while (i >= 1) {
        size_t li = i - 1; bool extended = false;
        const std::vector<int> orb = B.orbit[li];
        for (size_t oi = 0; oi < orb.size() && !extended; oi++) {
            int p = orb[oi];
            for (const auto& g : B.S[li]) {
                Perm<N> sg = composePerm<N>(invPerm<N>(B.urep[li][g[p]]), composePerm<N>(g, B.urep[li][p]));  // Schreier gen
                Perm<N> res; int dropLevel = bsgsStrip<N>(B, sg, res);
                bool isId = true; for (int x = 0; x < N; x++) if (res[x] != x) { isId = false; break; }
                if (!isId) {
                    if (dropLevel == (int)B.base.size()) {
                        int mp = -1; for (int x = 0; x < N; x++) if (res[x] != x) { mp = x; break; }
                        ensureLevel(mp); recompute(B.base.size() - 1);
                    }
                    int addAt = dropLevel < (int)B.base.size() ? dropLevel : (int)B.base.size() - 1;
                    for (size_t lvl = li + 1; lvl <= (size_t)addAt; lvl++) B.S[lvl].push_back(res);
                    for (size_t lvl = li + 1; lvl < B.base.size(); lvl++) recompute(lvl);
                    i = B.base.size(); extended = true; break;
                }
            }
        }
        if (!extended) i--;
    }
    return B;
}

template<int N> inline long long bsgsOrder(const BSGS<N>& B) { long long o = 1; for (size_t i = 0; i < B.base.size(); i++) o *= (long long)B.orbit[i].size(); return o; }

// Reduce a (possibly large/redundant) generating set to a bounded strong generating set for
// the same group, so generator counts do not explode down a recursion.
template<int N> inline std::vector<Perm<N>> reduceGens(const std::vector<Perm<N>>& gens) {
    if (gens.size() <= 4) return gens;
    BSGS<N> B = buildBSGS<N>(gens);
    std::set<Perm<N>> uniq;
    for (auto& lvl : B.S) for (auto& g : lvl) uniq.insert(g);
    return std::vector<Perm<N>>(uniq.begin(), uniq.end());
}

// Generators of the POINTWISE stabilizer of vertices a then b within <gens> = the subgroup
// fixing both endpoints of the anchor edge {a,b}. These fix the edge, so a candidate-
// restricted BFS driven by them stays inside the candidate set. Pointwise (rather than
// setwise) keeps it simple/robust; it under-approximates the setwise edge stabilizer by at
// most a factor of two -- sound, negligible.
template<int N> inline std::vector<Perm<N>> edgeStabGens(const std::vector<Perm<N>>& gens, int a, int b) {
    auto pointStab = [](const std::vector<Perm<N>>& gg, int pt) {
        std::array<Perm<N>, N> rep; std::array<char, N> vis; std::vector<int> orb;
        orbitTransversal<N>(pt, gg, rep, vis, orb);
        std::set<Perm<N>> S;
        for (int x : orb) for (const auto& g : gg) {
            Perm<N> sg = composePerm<N>(invPerm<N>(rep[g[x]]), composePerm<N>(g, rep[x]));
            bool id = true; for (int t = 0; t < N; t++) if (sg[t] != t) { id = false; break; }
            if (!id) S.insert(sg);
        }
        return reduceGens<N>(std::vector<Perm<N>>(S.begin(), S.end()));
    };
    return pointStab(pointStab(gens, a), b);
}

// Exact setwise stabilizer of a set of `factors` (perfect matchings) within the group whose
// BSGS is BG: the subgroup of elements that permute the factor set among themselves (an
// automorphism of the partial 1-factorization). cov[a][b] (a<b) = "edge {a,b} is covered",
// used for cheap partial pruning. Property-pruned backtrack over BG, collecting a generating
// set with K-minimality pruning (one element per coset of the stabilizer-so-far). Validated
// against brute-force enumeration (small stabilizers AND whole-group cases).
template<int N> inline std::vector<Perm<N>> setwiseStab(const BSGS<N>& BG,
                                                        const std::vector<Perm<N>>& factors, const bool cov[][N]) {
    int m = (int)BG.base.size();
    if (m == 0) return {};                                   // trivial group
    std::set<Perm<N>> fset(factors.begin(), factors.end());
    auto isCov = [&](int a, int b) { return a < b ? cov[a][b] : cov[b][a]; };
    std::vector<Perm<N>> Kgens;                              // stabilizer generators found so far

    auto predFull = [&](const Perm<N>& g) -> bool {          // g permutes the factor set
        for (auto& F : factors) if (!fset.count(applyPerm<N>(F, g.data()))) return false;
        return true;
    };
    auto partialOK = [&](const Perm<N>& cur, int l) -> bool {// covered-edge consistency among assigned base points
        int bl = BG.base[l], cl = cur[bl];
        for (int j = 0; j < l; j++) { int bj = BG.base[j], cj = cur[bj]; if (isCov(bl, bj) != isCov(cl, cj)) return false; }
        return true;
    };
    auto kMinimal = [&](const Perm<N>& cur, int l, int cl) -> bool {  // cl minimal in its Stab_K(prefix)-orbit?
        std::vector<Perm<N>> sk = Kgens;
        for (int j = 0; j < l; j++) {
            int fp = cur[BG.base[j]];
            std::array<Perm<N>, N> rep; std::array<char, N> vis; std::vector<int> orb; orbitTransversal<N>(fp, sk, rep, vis, orb);
            std::set<Perm<N>> ns;
            for (int x : orb) for (auto& g : sk) { Perm<N> s = composePerm<N>(invPerm<N>(rep[g[x]]), composePerm<N>(g, rep[x])); bool id = true; for (int t = 0; t < N; t++) if (s[t] != t) { id = false; break; } if (!id) ns.insert(s); }
            sk.assign(ns.begin(), ns.end());
        }
        std::array<Perm<N>, N> rep; std::array<char, N> vis; std::vector<int> orb; orbitTransversal<N>(cl, sk, rep, vis, orb);
        int mn = cl; for (int x : orb) if (x < mn) mn = x;
        return mn == cl;
    };

    std::function<void(int, Perm<N>)> rec = [&](int l, Perm<N> pre) {
        if (l == m) {
            bool id = true; for (int x = 0; x < N; x++) if (pre[x] != x) { id = false; break; }
            if (!id && predFull(pre)) Kgens.push_back(pre);
            return;
        }
        std::vector<int> orb = BG.orbit[l];
        std::sort(orb.begin(), orb.end(), [&](int a, int b) { return pre[a] < pre[b]; });   // ascending image order
        for (int o : orb) {
            Perm<N> cur = composePerm<N>(pre, BG.urep[l][o]);
            if (!partialOK(cur, l)) continue;
            if (!kMinimal(cur, l, cur[BG.base[l]])) continue;
            rec(l + 1, cur);
        }
    };
    rec(0, idPerm<N>());
    return Kgens;
}

}  // namespace cgt
