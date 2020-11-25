/**
 * @file NoramlMath.hpp
 * @author Yifan Zhang (scirocco_gti@yeah.net)
 * @brief simple math utils for non-xsimd
 * @date 2020-11-24 12:08:33
 * @modified: 2020-11-24 12:09:17
 */

#ifndef NORMALMATH_HPP
#define NORMALMATH_HPP

unsigned long hamming(unsigned long n) {
    n = (n & 0x5555555555555555) + ((n >> 1) & 0x5555555555555555);
    n = (n & 0x3333333333333333) + ((n >> 2) & 0x3333333333333333);
    n = (n & 0x0f0f0f0f0f0f0f0f) + ((n >> 4) & 0x0f0f0f0f0f0f0f0f);
    n = (n & 0x00ff00ff00ff00ff) + ((n >> 8) & 0x00ff00ff00ff00ff);
    n = (n & 0x0000ffff0000ffff) + ((n >> 16) & 0x0000ffff0000ffff);
    n = (n & 0x00000000ffffffff) + ((n >> 32) & 0x00000000ffffffff);
    return n;
}

#endif