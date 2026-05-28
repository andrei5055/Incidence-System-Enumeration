#pragma once
#include "KBase.h"

struct Mask256 { uint64_t m[3] = { 0, 0, 0 }; };

class KBipartBase : public KBase<Mask256> {
protected:
    KBipartBase(const FactorParams& factParam, bool bPrint) : KBase<Mask256>(factParam, bPrint) {}
    struct alignas(64) PackedAdj {
        uint64_t m[4] = { 0, 0, 0, 0 };
        __m128i v_e2o;
        __m128i v_o2e;
    };

    struct FastSortedFactor : FastSortedFactorBase<64> {};
    struct FastRowTriplet { FastSortedFactor r[3]; };

    virtual FastSortedFactor get_fast_sorted(const uint8_t* adj) = 0;

    void* cbClass = NULL;
    int kThreads;
    ResultCallback resultCallback;
    std::vector<PackedAdj> packed_pool;
    FastSortedFactor r1_can, r2_can;
};

inline auto lookup_v = [](__m256i table, __m256i idx) {
    __m256i table_swapped = _mm256_permute2x128_si256(table, table, 0x01);
    __m256i res_in_place = _mm256_shuffle_epi8(table, idx);
    __m256i res_swapped = _mm256_shuffle_epi8(table_swapped, idx);

    alignas(32) static const uint8_t dest_lane_data[32] = {
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0x10,0x10,0x10,0x10,0x10,0x10,0x10,0x10,
        0x10,0x10,0x10,0x10,0x10,0x10,0x10,0x10
    };
    __m256i dest_lane_bit = _mm256_load_si256((__m256i*)dest_lane_data);
    __m256i src_lane_bit = _mm256_and_si256(idx, _mm256_set1_epi8(0x10));
    __m256i match = _mm256_cmpeq_epi8(dest_lane_bit, src_lane_bit);
    return _mm256_blendv_epi8(res_swapped, res_in_place, match);
 };
