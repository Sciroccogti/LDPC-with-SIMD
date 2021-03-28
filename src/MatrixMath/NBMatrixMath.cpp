#include "MatrixMath/MatrixMath.hpp"

Eigen::MatrixXi NBtransform_H_to_G(const Eigen::MatrixXi& H) {
    int n_code = H.cols();

    // DOUBLE GAUSS-JORDAN:
    Eigen::MatrixXi Href_colonnes = H.transpose();
    Eigen::MatrixXi Q = gaussjordan(Href_colonnes).transpose();
    Eigen::MatrixXi Href_diag = Href_colonnes.transpose();
    gaussjordan(Href_diag);

    int n_bits = n_code - Href_diag.sum();

    Eigen::MatrixXi Y(n_code, n_bits);
    Y.setZero();

    // form identity
    for (size_t i = 0; i < n_bits; i++) {
        Y(n_code - n_bits + i, i) = 1;
    }

    return binaryproduct(Q, Y).transpose();
}

Eigen::MatrixXi NBtransform_H_to_G(const Eigen::SparseMatrix<int>& H) {
    Eigen::MatrixXi Hdense(H);
    return NBtransform_H_to_G(Hdense);
}