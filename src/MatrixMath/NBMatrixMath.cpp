#include <GF.hpp>

#include "MatrixMath/MatrixMath.hpp"

Eigen::MatrixXi NBtransform_H_to_G(const Eigen::MatrixXi& H, const int GF) {
    int n_code = H.cols();

    // DOUBLE GAUSS-JORDAN:
    Eigen::MatrixXi Href_colonnes = H.transpose();
    Eigen::MatrixXi Q_mat = gaussjordan(Href_colonnes).transpose();
    Eigen::MatrixXi Href_diag = Href_colonnes.transpose();
    gaussjordan(Href_diag);

    int n_bits = n_code - Href_diag.nonZeros(); //.sum();
    n_bits = 1;
    Eigen::MatrixXi Y(n_code, n_bits);
    Y.setZero();

    // form identity
    for (size_t i = 0; i < n_bits; i++) {
        Y(n_code - n_bits + i, i) = 1;
    }

    return NBproduct(Q_mat, Y, GF).transpose();
}

Eigen::MatrixXi NBtransform_H_to_G(const Eigen::SparseMatrix<int>& H,
                                   const int GF) {
    Eigen::MatrixXi Hdense(H);
    return NBtransform_H_to_G(Hdense, GF);
}

Eigen::MatrixXi NBproduct(const Eigen::MatrixXi& X, const Eigen::MatrixXi& Y,
                          const int GF) {
    assert(X.cols() == Y.rows());
    Eigen::MatrixXi ret(X.rows(), Y.cols());
    for (int i = 0; i < X.rows(); i++) {
        for (int j = 0; j < Y.cols(); j++) {
            int ret_ij = 0;
            for (int k = 0; k < X.cols(); k++) {
                ret_ij += GF_mul(X(i, k), Y(k, j), GF);
            }
            ret(i, j) = ret_ij;
        }
    }
    return ret;
}