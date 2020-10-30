/*
 * File: LDPC.hpp
 * File Created: Thursday, 29th October 2020 11:45:19
 * Author: Yifan Zhang (scirocco_gti@yeah.net)
 * Last Modified: Friday, 30th October 2020 18:55:56
 */

#ifndef LDPC_HPP
#define LDPC_HPP

#include "alist/Alist.hpp"

class LDPC {
  private:
    Eigen::SparseMatrix<int> H_mat;
    int K;  // length of message

  public:
    LDPC();
    LDPC(Alist<alist_matrix>);
    LDPC(const char* filename);
    ~LDPC();
    Eigen::RowVectorX<int> Encoder(Eigen::RowVectorX<int> m);
    Alist<alist_matrix> Decoder();
    Eigen::SparseMatrix<int> getG();
    Eigen::SparseMatrix<int> getH();
    int getK();
};

LDPC::LDPC() {
    H_mat = Eigen::SparseMatrix<int>();
    K = 0;
}

LDPC::LDPC(Alist<alist_matrix> A) {
    // G.nRow = N - M, G.nCol = N
    // N = n, M = n - k
    K = A.getnCol() - A.getnRow();
    H_mat = A.getMat();
}

LDPC::LDPC(const char* filename) {
    Alist<alist_matrix> A = Alist<alist_matrix>(filename);
    K = A.getnCol() - A.getnRow();
    H_mat = A.getMat();
}

LDPC::~LDPC() {}

// Eigen::RowVectorX<int> LDPC::Encoder(Eigen::RowVectorX<int> m) {}

Eigen::SparseMatrix<int> LDPC::getG() {
    return H_mat;
}

Eigen::SparseMatrix<int> LDPC::getH() {
    return H_mat;
}

int LDPC::getK() {
    return K;
}

#endif