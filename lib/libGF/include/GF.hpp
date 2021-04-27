/**
 * @file GF.hpp
 * @author Yifan Zhang (scirocco_gti@yeah.net)
 * @brief 0 -> 0 -> 0, 1 -> alpha0 -> 1, 2 -> alpha1 -> 2, 9 -> alpha8 -> 29
 * @date 2021-03-27 20:09:21
 * @modified: 2021-04-27 15:11:05
 */

#ifndef GF_HPP
#define GF_HPP

#include <cstdint>

// template <int Q>
// class GF {
//   private:
//     /* data */
//   public:
//     GF(/* args */);
//     ~GF();
// };

uint8_t GF_plus(const uint8_t &a, const uint8_t &b, const int Q);

uint8_t GF_mul(const uint8_t &a, const uint8_t &b, const int Q);
uint8_t GF_div(const uint8_t &a, const uint8_t &b, const int Q);

uint8_t GF_p2v(const uint8_t &a, const int Q);
uint8_t GF_v2p(const uint8_t &a, const int Q);

#endif