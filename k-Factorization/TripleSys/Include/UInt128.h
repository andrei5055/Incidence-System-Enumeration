#pragma once

#include <iostream>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <cstdint>

class UInt128
{
public:

    uint64_t hi;
    uint64_t lo;

    UInt128(uint64_t low = 0) : hi(0), lo(low)               {}
    UInt128(uint64_t high, uint64_t low) : hi(high), lo(low) {}

    // comparison
    bool operator<(const UInt128& b) const
    {
        if (hi != b.hi)
            return hi < b.hi;

        return lo < b.lo;
    }

    bool operator==(const UInt128& b) const { return hi == b.hi && lo == b.lo; }
    bool operator!=(const UInt128& b) const { return !(*this == b); }
    bool operator>(const UInt128& b) const  { return b < *this; }
    bool operator<=(const UInt128& b) const { return !(*this > b); }
    bool operator>=(const UInt128& b) const { return !(*this < b); }

    // arithmetic operations:
    UInt128 operator+(const UInt128& b) const {
        UInt128 r(lo + b.lo);
        r.hi = hi + b.hi + (r.lo < lo);
        return r;
    }

    UInt128 operator-(const UInt128& b) const {
        UInt128 r(lo - b.lo);
        r.hi = hi - b.hi - (lo < b.lo);
        return r;
    }

    // shifts
    UInt128 operator<<(int s) const {
        if (s == 0)
            return *this;

        if (s >= 128)
            return UInt128(0, 0);

        if (s >= 64)
            return UInt128(lo << (s - 64), 0);

        return UInt128(
            (hi << s) | (lo >> (64 - s)),
            lo << s
        );
    }

    UInt128 operator>>(int s) const {
        if (s == 0)
            return *this;

        if (s >= 128)
            return UInt128(0, 0);

        if (s >= 64)
            return UInt128(0, hi >> (s - 64));

        return UInt128(
            hi >> s,
            (lo >> s) | (hi << (64 - s))
        );
    }

    // bit access
    bool bit(int i) const {
        if (i < 64)
            return (lo >> i) & 1ULL;

        return (hi >> (i - 64)) & 1ULL;
    }

    void setBit(int i) {
        if (i < 64)
            lo |= (1ULL << i);
        else
            hi |= (1ULL << (i - 64));
    }

    UInt128 operator*(const UInt128& b) const {
        // schoolbook multiplication
        // using 32-bit chunks
        const uint64_t MASK = 0xffffffffULL;

        uint64_t c[8] = {};
        uint64_t a[4] = { lo & MASK, lo >> 32, hi & MASK, hi >> 32 };
        uint64_t bb[4] = { b.lo & MASK, b.lo >> 32, b.hi & MASK, b.hi >> 32 };
        for (int i = 0; i < 4; ++i) {
            uint64_t carry = 0;
            for (int j = 0; j < 4; ++j) {
                uint64_t k = i + j;
                uint64_t t = c[k] + a[i] * bb[j] + carry;
                c[k] = t & MASK;
                carry = t >> 32;
            }

            c[i + 4] += carry;
        }

		return UInt128(c[2] | (c[3] << 32), c[0] | (c[1] << 32));
    }

    UInt128 operator*(uint64_t x) const { return (*this) * UInt128(x); }
    friend UInt128 operator*(uint64_t a, const UInt128& b) {
        return UInt128(a) * b;
    }

    // division
    UInt128 operator/(const UInt128& d) const {
        if (d == UInt128(0))
            throw std::runtime_error("division by zero");

        UInt128 q, r;
        for (int i = 127; i >= 0; --i) {
            r = r << 1;

            if (bit(i))
                r.lo |= 1;

            if (r >= d) {
                r = r - d;
                q.setBit(i);
            }
        }

        return q;
    }

    UInt128 operator%(const UInt128& d) const {
        if (d == UInt128(0))
            throw std::runtime_error("division by zero");

        UInt128 r;
        for (int i = 127; i >= 0; --i) {
            r = r << 1;
            if (bit(i))
                r.lo |= 1;

            if (r >= d)
                r = r - d;
        }

        return r;
    }

    // decimal conversion
    std::string toString() const {
        if (*this == UInt128(0))
            return "0";

        UInt128 x = *this;
        std::string s;
        while (x != UInt128(0))
        {
            UInt128 q = x / UInt128(10);
            UInt128 r = x % UInt128(10);

            s.push_back(char('0' + r.lo));
            x = q;
        }

        std::reverse(s.begin(), s.end());
        return s;
    }
};

std::ostream& operator<<(std::ostream& os, const UInt128& x);