/*
 * File: LDPC.hpp
 * File Created: Thursday, 29th October 2020 11:45:19
 * Author: Yifan Zhang (scirocco_gti@yeah.net)
 * Last Modified: Thursday, 29th October 2020 22:38:43
 */

#ifndef LDPC_HPP
#define LDPC_HPP

#include "alist/Alist.hpp"

class LDPC {
  private:
    alist_matrix H_mat;
    int K;  // length of message

  public:
    LDPC();
    LDPC(Alist<alist_matrix>);
    LDPC(const char* filename);
    ~LDPC();
    Alist<alist_matrix> Encoder();
    Alist<alist_matrix> Decoder();
    Alist<alist_matrix> getH();
    int getK();
};

LDPC::LDPC() {
    H_mat = Alist<alist_matrix>().getData();
    K = 0;
}

LDPC::LDPC(Alist<alist_matrix> d) {
    H_mat = d.getData();
    // G.nRow = N - M, G.nCol = N
    // N = n, M = n - k
    K = H_mat.N - H_mat.M;
}

LDPC::LDPC(const char* filename) {
    H_mat = Alist<alist_matrix>(filename).getData();
    K = H_mat.N - H_mat.M;
}

LDPC::~LDPC() {}

Alist<alist_matrix> LDPC::getH() {
    return H_mat;
}

int LDPC::getK() {
    return K;
}

#endif