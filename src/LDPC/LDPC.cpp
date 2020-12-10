#include "LDPC/LDPC.hpp"

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

// Eigen::RowVectorXi LDPC::Encoder(Eigen::RowVectorXi m) {}

Eigen::SparseMatrix<int> LDPC::getG() {
    return G_mat;
}

Eigen::SparseMatrix<int> LDPC::getH() {
    return H_mat;
}

int LDPC::getK() {
    return K;
}
