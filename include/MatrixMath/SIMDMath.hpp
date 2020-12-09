/**
 * @file SIMDMath.hpp
 * @author Yifan Zhang (scirocco_gti@yeah.net)
 * @brief simple math utils for simd
 * @date 2020-11-24 11:24:04
 * @modified: 2020-11-24 11:24:33
 */

#ifndef SIMDMATH_HPP
#define SIMDMATH_HPP

#include <immintrin.h>

using b_type = __m256i;
// using vector_type = std::vector<u_int64_t>;
const size_t inc = sizeof(b_type) / sizeof(u_int64_t);
/*
using b_type = mipp::Reg<int64_t>;
using vector_type = mipp::vector<int64_t>;
size_t inc = mipp::N<int64_t>();
*/
/*
using b_type = xsimd::simd_type<u_int64_t>;
using vector_type = std::vector<u_int64_t, XSIMD_DEFAULT_ALLOCATOR(u_int64_t)>;
size_t inc = b_type::size;
*/

u_int64_t v0[inc] = {0xfffffffffffffffe, 0xfffffffffffffffe, 0xfffffffffffffffe,
                     0xfffffffffffffffe};
u_int64_t v1[inc] = {0x5555555555555555, 0x5555555555555555, 0x5555555555555555,
                     0x5555555555555555};
u_int64_t v2[inc] = {0x3333333333333333, 0x3333333333333333, 0x3333333333333333,
                     0x3333333333333333};
u_int64_t v3[inc] = {0x0f0f0f0f0f0f0f0f, 0x0f0f0f0f0f0f0f0f, 0x0f0f0f0f0f0f0f0f,
                     0x0f0f0f0f0f0f0f0f};
u_int64_t v4[inc] = {0x00ff00ff00ff00ff, 0x00ff00ff00ff00ff, 0x00ff00ff00ff00ff,
                     0x00ff00ff00ff00ff};
u_int64_t v5[inc] = {0x0000ffff0000ffff, 0x0000ffff0000ffff, 0x0000ffff0000ffff,
                     0x0000ffff0000ffff};
u_int64_t v6[inc] = {0x00000000ffffffff, 0x00000000ffffffff, 0x00000000ffffffff,
                     0x00000000ffffffff};

/*
b_type b1 = xsimd::load_aligned(&v1[0]);
b_type b2 = xsimd::load_aligned(&v2[0]);
b_type b3 = xsimd::load_aligned(&v3[0]);
b_type b4 = xsimd::load_aligned(&v4[0]);
b_type b5 = xsimd::load_aligned(&v5[0]);
b_type b6 = xsimd::load_aligned(&v6[0]);
*/
/*
b_type b0 = &v0[0];
b_type b1 = &v1[0];
b_type b2 = &v2[0];
b_type b3 = &v3[0];
b_type b4 = &v4[0];
b_type b5 = &v5[0];
b_type b6 = &v6[0];
*/
b_type b0 = _mm256_load_si256((__m256i *)v0);
b_type b1 = _mm256_load_si256((__m256i *)v1);
b_type b2 = _mm256_load_si256((__m256i *)v2);
b_type b3 = _mm256_load_si256((__m256i *)v3);
b_type b4 = _mm256_load_si256((__m256i *)v4);
b_type b5 = _mm256_load_si256((__m256i *)v5);
b_type b6 = _mm256_load_si256((__m256i *)v6);

b_type hamming(b_type &n) {
    n = (n & b1) + ((n >> 1) & b1);
    n = (n & b2) + ((n >> 2) & b2);
    n = (n & b3) + ((n >> 4) & b3);
    n = (n & b4) + ((n >> 8) & b4);
    n = (n & b5) + ((n >> 16) & b5);
    n = (n & b6) + ((n >> 32) & b6);
    return n;
}

#endif