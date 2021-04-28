/**
 * @file NBLLR.hpp
 * @author Yifan Zhang (scirocco_gti@yeah.net)
 * @brief
 * @date 2021-04-20 00:30:58
 * @modified: 2021-04-27 16:14:32
 */

#ifndef NBLLR_HPP
#define NBLLR_HPP

#include <cstdint>

/**
 * @brief class for NBLLR, has two value: 1. LLR; 2. corresponding Q
 * 
 */
class NBLLR {
  private:
    float data;
    uint8_t Q;  // should in [1, GF]
  public:
    NBLLR();
    NBLLR(float d, uint8_t q);
    ~NBLLR();
    float getLLR() const;
    uint8_t getQ();
};

bool operator<(const NBLLR& llr1, const NBLLR& llr2);
bool operator>(const NBLLR& llr1, const NBLLR& llr2);

#endif