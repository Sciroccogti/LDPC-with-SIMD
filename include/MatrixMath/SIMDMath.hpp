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
const size_t inc = sizeof(b_type) / sizeof(double);
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

double v0[inc] = {(double)0xfffffffffffffffe, (double)0xfffffffffffffffe,
                  (double)0xfffffffffffffffe, (double)0xfffffffffffffffe};
double v1[inc] = {(double)0x5555555555555555, (double)0x5555555555555555,
                  (double)0x5555555555555555, (double)0x5555555555555555};
double v2[inc] = {(double)0x3333333333333333, (double)0x3333333333333333,
                  (double)0x3333333333333333, (double)0x3333333333333333};
double v3[inc] = {(double)0x0f0f0f0f0f0f0f0f, (double)0x0f0f0f0f0f0f0f0f,
                  (double)0x0f0f0f0f0f0f0f0f, (double)0x0f0f0f0f0f0f0f0f};
double v4[inc] = {(double)0x00ff00ff00ff00ff, (double)0x00ff00ff00ff00ff,
                  (double)0x00ff00ff00ff00ff, (double)0x00ff00ff00ff00ff};
double v5[inc] = {(double)0x0000ffff0000ffff, (double)0x0000ffff0000ffff,
                  (double)0x0000ffff0000ffff, (double)0x0000ffff0000ffff};
double v6[inc] = {(double)0x00000000ffffffff, (double)0x00000000ffffffff,
                  (double)0x00000000ffffffff, (double)0x00000000ffffffff};

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
b_type b0 = _mm256_cvtpd_epu64(_mm256_load_pd(v0));
b_type b1 = _mm256_cvtpd_epu64(_mm256_load_pd(v1));
b_type b2 = _mm256_cvtpd_epu64(_mm256_load_pd(v2));
b_type b3 = _mm256_cvtpd_epu64(_mm256_load_pd(v3));
b_type b4 = _mm256_cvtpd_epu64(_mm256_load_pd(v4));
b_type b5 = _mm256_cvtpd_epu64(_mm256_load_pd(v5));
b_type b6 = _mm256_cvtpd_epu64(_mm256_load_pd(v6));

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