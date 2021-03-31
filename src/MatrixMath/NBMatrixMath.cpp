#include <GF.hpp>
#include <iostream>

#include "MatrixMath/MatrixMath.hpp"

/**
 * @brief Returns the index of the maximum value
 *
 * @param a : the input matrix
 * @return int
 */
int NBargmax(const Eigen::MatrixXi a) {
    int tmp = a(0, 0);
    int ret = 0;
    for (size_t i = 0; i < a.rows(); i++) {
        for (size_t j = 0; j < a.cols(); j++) {
            if (a(i, j) > tmp) {
                tmp = a(i, j);
                ret = i * a.cols() + j;
            }
        }
    }
    return 0;
}

Eigen::MatrixXi NBcoefdot(const Eigen::MatrixXi& M, const short c,
                          const int GF) {
    Eigen::MatrixXi ret = M;
    for (size_t i = 0; i < ret.rows(); i++) {
        for (size_t j = 0; j < ret.cols(); j++) {
            ret(i, j) = GF_mul(M(i, j), c, GF);
        }
    }
    return ret;
}

Eigen::MatrixXi NBtransform_H_to_G(const Eigen::MatrixXi& H, const int GF) {
    int n_code = H.cols();

    // DOUBLE GAUSS-JORDAN:
    Eigen::MatrixXi Href_colonnes = H.transpose();
    Eigen::MatrixXi Q_mat = NBgaussjordan(Href_colonnes, GF).transpose();
    std::cout <<"221 * 1 = "<< GF_mul(221, 1, GF) << std::endl;
    std::cout << "Q * H:\n" << NBproduct(H, Q_mat.transpose(), GF) << std::endl;
    std::cout << "Q:\n" << Q_mat << std::endl;
    Eigen::MatrixXi Href_diag = Href_colonnes.transpose();
    NBgaussjordan(Href_diag, GF);
    int n_bits =
        n_code -
        Href_diag.unaryExpr([](const int x) { return x == 0 ? 0 : 1; }).sum();
    Eigen::MatrixXi Y(n_code, n_bits);
    Y.setZero();

    // form identity
    for (size_t i = 0; i < n_bits; i++) {
        Y(n_code - n_bits + i, i) = 1;
    }

    std::cout << "Y:\n" << Y << std::endl;
    Eigen::MatrixXi ret = NBproduct(Q_mat, Y, GF).transpose();
    // for (int i = 0; i < n_bits; i++) {
    //     int coef = GF_div(1, ret(i, ret.rows() + i), GF);
    //     ret.block(i, 0, 1, ret.cols()) =
    //         NBcoefdot(ret.block(i, 0, 1, ret.cols()), coef, GF);
    // }

    return ret;
}

Eigen::MatrixXi NBtransform_H_to_G(const Eigen::SparseMatrix<int>& H,
                                   const int GF) {
    Eigen::MatrixXi Hdense(H);
    return NBtransform_H_to_G(Hdense, GF);
}

/**
 * @brief Compute the binary row reduced echelon form of X.
 *        That is, turning X into up-right triangle.
 *        https://github.com/hichamjanati/pyldpc/blob/master/pyldpc/utils.py#L38
 * @param X : will be reused to return result!
 * @param GF : GF order
 * @return Eigen::MatrixXi: the inverse transform
 */
Eigen::MatrixXi NBgaussjordan(Eigen::MatrixXi& X, const int GF) {
    size_t m = X.rows(), n = X.cols();
    Eigen::MatrixXi P = Eigen::MatrixXi::Identity(m, m);

    int pivot_old = -1;
    for (size_t j = 0; j < n; j++) {
        // current column
        Eigen::MatrixXi filtre_down =
            X.block(pivot_old + 1, j, m - pivot_old - 1, 1);

        // The first 1 found in filtre_down.
        // `argmax()` returns the index of the first nonzero .
        int pivot = argmax(filtre_down) + pivot_old + 1;

        if (X(pivot, j)) {
            pivot_old++;
            if (pivot_old != pivot) {
                swap_rows(X, pivot, pivot_old);
                swap_rows(P, pivot, pivot_old);
            }

            for (size_t i = 0; i < m; i++) {
                // row(i) -= row(pivot_old), as X(i, j) == 1
                if (i != pivot_old && X(i, j)) {
                    int coef = GF_div(X(i, j), X(pivot_old, j), GF);
                    P.row(i) = NBplus(
                        P.row(i), NBcoefdot(P.row(pivot_old), coef, GF), GF);
                    X.row(i) = (NBplus(
                        X.row(i), NBcoefdot(X.row(pivot_old), coef, GF), GF));
                }
            }
        }

        if (pivot_old == m - 1) {
            break;
        }
    }

    return P;
}

Eigen::MatrixXi NBproduct(const Eigen::MatrixXi& X, const Eigen::MatrixXi& Y,
                          const int GF) {
    assert(X.cols() == Y.rows());
    Eigen::MatrixXi ret(X.rows(), Y.cols());
    for (int i = 0; i < X.rows(); i++) {
        for (int j = 0; j < Y.cols(); j++) {
            int ret_ij = 0;
            for (int k = 0; k < X.cols(); k++) {
                // printf("%d %d %d\n", X(i, k), Y(k, j), GF_mul(X(i, k), Y(k, j), GF));
                ret_ij = GF_plus(ret_ij, GF_mul(X(i, k), Y(k, j), GF), GF);
            }
            ret(i, j) = ret_ij;
        }
    }
    return ret;
}

Eigen::MatrixXi NBplus(const Eigen::MatrixXi& X, const Eigen::MatrixXi& Y,
                       const int GF) {
    assert(X.cols() == Y.cols() && X.rows() == Y.rows());
    Eigen::MatrixXi ret = X;
    for (int i = 0; i < Y.rows(); i++) {
        for (int j = 0; j < Y.cols(); j++) {
            ret(i, j) = GF_plus(X(i, j), Y(i, j), GF);
        }
    }
    return ret;
}