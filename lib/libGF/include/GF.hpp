/**
 * @file GF.hpp
 * @author Yifan Zhang (scirocco_gti@yeah.net)
 * @brief 0 -> 0 -> 0, 1 -> alpha0 -> 1, 2 -> alpha1 -> 2, 9 -> alpha8 -> 29
 * @date 2021-03-27 20:09:21
 * @modified: 2021-04-01 21:54:51
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

short GF_plus(const short &a, const short &b, const int Q);

short GF_mul(const short &a, const short &b, const int Q);
short GF_div(const short &a, const short &b, const int Q);

short GF_p2v(const short &a, const int Q);
short GF_v2p(const short &a, const int Q);

#endif