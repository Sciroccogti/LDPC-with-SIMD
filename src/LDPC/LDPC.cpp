#include "LDPC/LDPC.hpp"

LDPC::LDPC() {
    H_mat = Eigen::SparseMatrix<int>();
    G_mat = Eigen::SparseMatrix<int>();
    K = 0;
}

LDPC::LDPC(Eigen::SparseMatrix<int> H) {
    // G.nRow = N - M, G.nCol = N
    // N = n, M = n - k
    // K = A.getnCol() - A.getnRow();
    H_mat = H;
    G_mat = Eigen::SparseMatrix<int>(transform_H_to_G(H_mat).sparseView());
    K = G_mat.rows();
}

LDPC::LDPC(Alist<alist_matrix> A) {
    *this = LDPC(A.getMat());
}

LDPC::LDPC(const char* filename) {
    Alist<alist_matrix> A = Alist<alist_matrix>(filename);
    *this = LDPC(A.getMat());
}

LDPC::~LDPC() {}

Eigen::RowVectorXi LDPC::encode(Eigen::RowVectorXi& m) {
    return binaryproduct(m, G_mat.toDense());
}

Eigen::SparseMatrix<int> LDPC::getG() {
    return G_mat;
}

Eigen::SparseMatrix<int> LDPC::getH() {
    return H_mat;
}

int LDPC::getK() {
    return K;
}
