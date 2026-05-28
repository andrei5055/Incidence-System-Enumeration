#include <iostream>
#include <vector>
#include <utility>
#include <iomanip>

using namespace std;

typedef unsigned char tchar;

void triangularGraph(tchar N, tchar *A, bool verifyMatrix)
{
    const int v = N * (N - 1) / 2;
    // ------------------------------------------------------------
    // Enumerate all 2-subsets (a,b), a < b
    // ------------------------------------------------------------

    vector<pair<tchar, tchar>> vertex_of;
    vertex_of.reserve(v);

    for (tchar a = 0; a < N; ++a)
    {
        for (auto b = a + 1; b < N; ++b)
            vertex_of.push_back({ a, b });
    }

    // ------------------------------------------------------------
    // Build adjacency matrix
    // ------------------------------------------------------------
    for (int i = 0; i < v; ++i) {
        const auto a = vertex_of[i].first;
        const auto b = vertex_of[i].second;

        for (int j = i + 1; j < v; ++j) {
            const auto c = vertex_of[j].first;
            const auto d = vertex_of[j].second;

            // Count common elements
            int common = 0;

            if (a == c || a == d) ++common;
            if (b == c || b == d) ++common;

            // Adjacent iff exactly one common element
            if (common == 1)
                A[j * v + i] = A[i * v + j] = 1;
        }
    }

    if (!verifyMatrix)
		return;

    // ------------------------------------------------------------
    // Print basic verification
    // ------------------------------------------------------------
    long long edges = 0;
    const auto* pRow = A;
    for (int i = 0; i < v; ++i, pRow += v) {
        int deg = 0;
        for (int j = 0; j < v; ++j)
            deg += pRow[j];
        
        cout << "deg(" << i << ") = " << deg << "\n";
        edges += deg;
    }

    edges /= 2;

    cout << "\n";
    cout << "Vertices : " << v << "\n";
    cout << "Edges    : " << edges << "\n";

    // ------------------------------------------------------------
    // Print adjacency matrix
    // ------------------------------------------------------------

    cout << "\nAdjacency matrix:\n\n";

    pRow = A;
    for (int i = 0; i < v; ++i) {
        for (int j = 0; j < v; ++j)
            cout << *pRow++;

        cout << "\n";
    }
}

