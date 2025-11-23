// generateLatinSquare.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
// converted from Java to c++ to calculate Latin Squares. 
// Java code was developed by Paul Hankin.
// According to Paul Hankin: "his Java implementation loosely based on the Java implementation described 
// by Ignacio Gallego Sagastume which implements the rather ingenious method of Jacobson and Matthews".
// For more details about Java implementation see https://blog.paulhankin.net/latinsquares/

#include <iostream>
#include <tuple>
#include <vector>
#include <random>
#include <ctime>
#include <array>
using namespace std;

std::mt19937 rng(42);
std::pair<int, int> rand2(int a, int b) {
    std::uniform_int_distribution<int> dist(0, 1);
    if (dist(rng) == 0) {
        return { a, b };
    }
    return { b, a };
}

std::vector<std::vector<int>> Latin(int n) {
    std::vector<std::vector<int>> xy(n, std::vector<int>(n));
    std::vector<std::vector<int>> xz(n, std::vector<int>(n));
    std::vector<std::vector<int>> yz(n, std::vector<int>(n));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int k = (i + j) % n;
            xy[i][j] = k;
            xz[i][k] = j;
            yz[j][k] = i;
        }
    }

    int mxy = 0, mxz = 0, myz = 0;
    array<int, 3> m;
    bool proper = true;
    int minIter = n * n * n;

    std::uniform_int_distribution<int> dist_n(0, n - 1);

    for (int iter = 0; iter < minIter || !proper; iter++) {
        int i, j, k, i2, j2, k2;
        int i2_, j2_, k2_;

        if (proper) {
            i = dist_n(rng);
            j = dist_n(rng);
            k = dist_n(rng);
            while (xy[i][j] == k) {
                i = dist_n(rng);
                j = dist_n(rng);
                k = dist_n(rng);
            }
            i2 = yz[j][k];
            j2 = xz[i][k];
            k2 = xy[i][j];
            i2_ = i;
            j2_ = j;
            k2_ = k;
        }
        else {
            i = m[0];
            j = m[1];
            k = m[2];
            std::tie(i2, i2_) = rand2(yz[j][k], myz);
            std::tie(j2, j2_) = rand2(xz[i][k], mxz);
            std::tie(k2, k2_) = rand2(xy[i][j], mxy);
        }

        proper = (xy[i2][j2] == k2);
        if (!proper) {
            m = { i2, j2, k2 };
            mxy = xy[i2][j2];
            myz = yz[j2][k2];
            mxz = xz[i2][k2];
        }

        xy[i][j] = k2_;
        xy[i][j2] = k2;
        xy[i2][j] = k2;
        xy[i2][j2] = k;

        yz[j][k] = i2_;
        yz[j][k2] = i2;
        yz[j2][k] = i2;
        yz[j2][k2] = i;

        xz[i][k] = j2_;
        xz[i][k2] = j2;
        xz[i2][k] = j2;
        xz[i2][k2] = j;
    }
    return xy;
}

int getLS(unsigned char * ls, const int n, const int iRandomStart) {
    rng.seed(iRandomStart);
    auto latin_square = Latin(n);
    for (const auto& row : latin_square) {
        for (int c : row) {
            //std::cout << c << "  ";
            *ls = c;
            ls++;
        }
        //std::cout << "\n";
    }
    return 0;
}
