/**
 * @file LDPC.hpp
 * @author Yifan Zhang (scirocco_gti@yeah.net)
 * @brief
 * @date 2020-10-30 22:52:43
 * @modified: 2020-10-31 12:04:23
 */

#ifndef LDPC_HPP
#define LDPC_HPP

#include "Alist/Alist.hpp"
#include "MatrixMath/MatrixMath.hpp"

class LDPC {
  private:
    Eigen::SparseMatrix<int> H_mat;
    Eigen::SparseMatrix<int> G_mat;
    int K;  // length of message

  public:
    LDPC();
    LDPC(Alist<alist_matrix>);
    LDPC(const char* filename);
    // TODO: add H check
    ~LDPC();
    Eigen::RowVectorX<int>& Encoder(Eigen::RowVectorX<int>& m);
    Alist<alist_matrix> Decoder();
    Eigen::SparseMatrix<int> getG();
    Eigen::SparseMatrix<int> getH();
    int getK();
};

LDPC::LDPC() {
    H_mat = Eigen::SparseMatrix<int>();
    G_mat = Eigen::SparseMatrix<int>();
    K = 0;
}

LDPC::LDPC(Alist<alist_matrix> A) {
    // G.nRow = N - M, G.nCol = N
    // N = n, M = n - k
    K = A.getnCol() - A.getnRow();
    H_mat = A.getMat();
    G_mat = Eigen::SparseMatrix<int>(transform_H_to_G(H_mat).sparseView());
}

LDPC::LDPC(const char* filename) {
    Alist<alist_matrix> A = Alist<alist_matrix>(filename);
    K = A.getnCol() - A.getnRow();
    H_mat = A.getMat();
    auto G = transform_H_to_G(H_mat);
    G_mat = Eigen::SparseMatrix<int>(G.sparseView());
}

LDPC::~LDPC() {}

// Eigen::RowVectorX<int> LDPC::Encoder(Eigen::RowVectorX<int> m) {}

Eigen::SparseMatrix<int> LDPC::getG() {
    return G_mat;
}

Eigen::SparseMatrix<int> LDPC::getH() {
    return H_mat;
}

int LDPC::getK() {
    return K;
}

#endif