/**
 * @file Modem.hpp
 * @author Yifan Zhang (scirocco_gti@yeah.net)
 * @brief
 * @date 2020-12-11 13:52:36
 * @modified: 2020-12-11 14:03:59
 */

#ifndef MODEM_HPP
#define MODEM_HPP

#include "MatrixMath/MatrixMath.hpp"

#define MODEM_BPSK 0

class Modem {
  private:
    int fc;     // carrier frequency
    int fs;     // sampling frequency
    double Tb;  // secs of each symbol
    int type;   // choose from MODEM_BPSK
    int L;      // number of samples in the time of one symbol

  public:
    Modem(int freqc, int freqs, double Timeb, int t = MODEM_BPSK);
    ~Modem();
    Eigen::RowVectorXd modulate(const Eigen::RowVectorXi &c);
};

#endif