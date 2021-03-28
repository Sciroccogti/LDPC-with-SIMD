#include "LDPC/NBLDPC.hpp"
#include "MatrixMath/MatrixMath.hpp"

NBLDPC::NBLDPC() {
    H_mat = Eigen::SparseMatrix<int>();
    G_mat = Eigen::SparseMatrix<int>();
    K = 0;
    N = 0;
    GF = 0;
}

NBLDPC::NBLDPC(Eigen::SparseMatrix<int> H) {
    // G.nRow = N - M, G.nCol = N
    // N = n, M = n - k
    // K = A.getnCol() - A.getnRow();
    H_mat = H;
    G_mat = Eigen::SparseMatrix<int>(NBtransform_H_to_G(H_mat).sparseView());
    assert(
        !(binaryproduct(G_mat, H_mat.transpose()).any()));  // G*H should be 0
    K = G_mat.rows();
    N = G_mat.cols();
    Eigen::MatrixXi diff =
        G_mat.block(0, N - K, K, K) - Eigen::MatrixXi::Identity(K, K);
    if (diff.any()) {
        isSystematic = false;
    } else {
        isSystematic = true;
    }
}
