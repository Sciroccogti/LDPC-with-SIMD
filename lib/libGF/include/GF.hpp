/**
 * @file GF.hpp
 * @author Yifan Zhang (scirocco_gti@yeah.net)
 * @brief
 * @date 2021-03-27 20:09:21
 * @modified: 2021-03-28 17:02:29
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

short GF_2pow(const short &a, const int Q);
short GF_log2(const short &a, const int Q);

#endif