/**
 * @file SIMDMath.hpp
 * @author Yifan Zhang (scirocco_gti@yeah.net)
 * @brief simple math utils for xsimd
 * @date 2020-11-24 11:24:04
 * @modified: 2020-11-24 11:24:33
 */

#ifndef SIMDMATH_HPP
#define SIMDMATH_HPP
#include "xsimd/xsimd.hpp"

using b_type = xsimd::simd_type<u_int64_t>;
using vector_type = std::vector<u_int64_t, XSIMD_DEFAULT_ALLOCATOR(u_int64_t)>;

size_t inc = b_type::size;
vector_type v1(inc, 0x5555555555555555);
vector_type v2(inc, 0x3333333333333333);
vector_type v3(inc, 0x0f0f0f0f0f0f0f0f);
vector_type v4(inc, 0x00ff00ff00ff00ff);
vector_type v5(inc, 0x0000ffff0000ffff);
vector_type v6(inc, 0x00000000ffffffff);
b_type b1 = xsimd::load_aligned(&v1[0]);
b_type b2 = xsimd::load_aligned(&v2[0]);
b_type b3 = xsimd::load_aligned(&v3[0]);
b_type b4 = xsimd::load_aligned(&v4[0]);
b_type b5 = xsimd::load_aligned(&v5[0]);
b_type b6 = xsimd::load_aligned(&v6[0]);

b_type hamming(b_type n) {
    n = (n & b1) + ((n >> 1) & b1);
    n = (n & b2) + ((n >> 2) & b2);
    n = (n & b3) + ((n >> 4) & b3);
    n = (n & b4) + ((n >> 8) & b4);
    n = (n & b5) + ((n >> 16) & b5);
    n = (n & b6) + ((n >> 32) & b6);
    return n;
}

#endif