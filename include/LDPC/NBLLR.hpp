/**
 * @file NBLLR.hpp
 * @author Yifan Zhang (scirocco_gti@yeah.net)
 * @brief
 * @date 2021-04-20 00:30:58
 * @modified: 2021-04-21 13:02:45
 */

#ifndef NBLLR_HPP
#define NBLLR_HPP

/**
 * @brief class for NBLLR, has two value: 1. LLR; 2. corresponding Q
 * 
 */
class NBLLR {
  private:
    double data;
    int Q;  // should in [1, GF]
  public:
    NBLLR(double d, int q);
    ~NBLLR();
    double getLLR() const;
    int getQ();
};

bool operator<(const NBLLR& llr1, const NBLLR& llr2);
bool operator>(const NBLLR& llr1, const NBLLR& llr2);

#endif