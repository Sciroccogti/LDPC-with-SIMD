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
    G_mat =
        Eigen::SparseMatrix<int>(NBtransform_H_to_G(H_mat, GF).sparseView());
    assert(
        !(NBproduct(G_mat, H_mat.transpose(), GF).any()));  // G*H should be 0
    K = G_mat.rows();
    N = G_mat.cols();
    Eigen::MatrixXi diff =
        G_mat.block(0, N - K, K, K).unaryExpr([](const int x) {
            return x == 0 ? 0 : 1;
        }) -
        Eigen::MatrixXi::Identity(K, K);
    if (diff.any()) {
        isSystematic = false;
    } else {
        isSystematic = true;
    }
}

NBLDPC::NBLDPC(Alist<nbalist_matrix> A) {
    GF = A.getGF();
    memset(num_mlist, 0, A.getnRow() * sizeof(int));
    memcpy(num_mlist, A.getData().num_mlist, A.getnRow() * sizeof(int));
    memset(num_nlist, 0, A.getnCol() * sizeof(int));
    memcpy(num_nlist, A.getData().num_nlist, A.getnCol() * sizeof(int));
    new (this) NBLDPC(A.getMat());  // no need to delete, cause memory is reused
}

NBLDPC::NBLDPC(const char* filename) {
    Alist<nbalist_matrix> A = Alist<nbalist_matrix>(filename);
    GF = A.getGF();
    num_mlist = (int*)malloc(A.getnRow() * sizeof(int));
    memcpy(num_mlist, A.getData().num_mlist, A.getnRow() * sizeof(int));
    num_nlist = (int*)malloc(A.getnCol() * sizeof(int));
    memcpy(num_nlist, A.getData().num_nlist, A.getnCol() * sizeof(int));
    new (this) NBLDPC(A.getMat());  // no need to delete, cause memory is reused
}

NBLDPC::~NBLDPC() {
    free(num_mlist);
    free(num_nlist);
}
