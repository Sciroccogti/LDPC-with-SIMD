/**
 * @file Modem.hpp
 * @author Yifan Zhang (scirocco_gti@yeah.net)
 * @brief
 * @date 2020-12-11 13:52:36
 * @modified: 2021-08-08 22:00:48
 */

#ifndef MODEM_HPP
#define MODEM_HPP

#include "MatrixMath/MatrixMath.hpp"

#define MODEM_BPSK 0

class Modem {
  private:
    int fc;     // carrier frequency
    int fs;     // sampling frequency
    float Tb;  // secs of each symbol
    int type;   // choose from MODEM_BPSK
    int L;      // number of samples in the time of one symbol

  public:
    Modem(int freqc, int freqs, float Timeb, int t = MODEM_BPSK);
    ~Modem();
    Eigen::RowVectorXf modulate(const Eigen::RowVectorXi &c);
    Eigen::RowVectorXf demodulate(const Eigen::RowVectorXf &s);
    int getL();
};

int Compare(const Eigen::MatrixXi &X, const Eigen::MatrixXi &Y);
int Comparef(const Eigen::MatrixXf &X, const Eigen::MatrixXf &Y);
int CompareBPSK(const Eigen::RowVectorXi &X, const Eigen::RowVectorXi &Y);
Eigen::MatrixXi TransBPSK(const Eigen::MatrixXi& X);
Eigen::MatrixXi RetransBPSK(const Eigen::MatrixXi &X);

#endif