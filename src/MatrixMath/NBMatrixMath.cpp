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

    Eigen::MatrixXi ret = NBproduct(Q_mat, Y, GF).transpose();
    // Eigen::MatrixXi Q = H;
    // int n_bits = NBGauss(Q, GF);
    // Eigen::MatrixXi Q_cut =
    //     Q.transpose().block(n_bits, 0, Q.cols() - n_bits, Q.rows());
    // Eigen::MatrixXi I = Eigen::MatrixXi::Identity(Q_cut.rows(),
    // Q_cut.rows()); Eigen::MatrixXi ret(Q_cut.rows(), Q_cut.cols() +
    // Q_cut.rows()); ret << Q_cut, I;

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

/**
 * @brief
 * https://github.com/kir1994/LDPC/blob/686d33b9bc818e1d050e91ad1ecef56c24299b8b/LDPC/GFLinAlg.cpp#L36
 *
 * @param X
 * @param GF
 * @return int
 */
int NBGauss(Eigen::MatrixXi& X, const int GF) {
    int n_rows = X.rows();
    Eigen::RowVectorXi permute(X.cols());
    for (size_t i = 0; i < X.cols(); i++) {
        permute[i] = i;
    }

    for (size_t i = 0; i < n_rows; i++) {
        bool isChanged = false;
        for (size_t col_cur = i; col_cur < X.cols(); col_cur++) {
            int C = permute[col_cur];
            for (size_t j = i; j < n_rows; j++) {
                if (X(j, C)) {
                    isChanged = true;
                    int coef = GF_div(1, X(j, C), GF);
                    X.row(j) = NBcoefdot(X.row(j), coef, GF);
                    if (j > i) {
                        for (size_t k = 0; k < X.cols(); k++) {
                            X(i, k) ^= X(j, k);
                        }
                    }
                    break;
                }
            }

            if (isChanged && col_cur != i) {
                int tmp = permute[col_cur];
                permute[col_cur] = permute[i];
                permute[i] = tmp;
            }
            break;
        }
        if (!isChanged) {
            n_rows = i;
            break;
        }

        int C = permute[i];
        for (size_t j = 0; j < n_rows; j++) {
            if (j == i) {
                continue;
            }
            if (X(j, C)) {
                int ret_jC = X(j, C);
                for (size_t k = 0; k < X.cols(); k++) {
                    int tmp = GF_mul(X(i, k), ret_jC, GF);
                    X(j, k) ^= tmp;
                }
            }
        }
    }

    return n_rows;
}

// TODO: #10 really should use xor?!
Eigen::MatrixXi NBproduct(const Eigen::MatrixXi& X, const Eigen::MatrixXi& Y,
                          const int GF) {
    assert(X.cols() == Y.rows());
    Eigen::MatrixXi ret(X.rows(), Y.cols());
    for (int i = 0; i < X.rows(); i++) {
        for (int j = 0; j < Y.cols(); j++) {
            int ret_ij = 0;
            for (int k = 0; k < X.cols(); k++) {
                ret_ij = ret_ij ^ GF_mul(X(i, k), Y(k, j), GF);
            }
            ret(i, j) = ret_ij;
        }
    }
    return ret;
}

// TODO: really should use xor?!
Eigen::MatrixXi NBplus(const Eigen::MatrixXi& X, const Eigen::MatrixXi& Y,
                       const int GF) {
    assert(X.cols() == Y.cols() && X.rows() == Y.rows());
    Eigen::MatrixXi ret = X;
    for (int i = 0; i < Y.rows(); i++) {
        for (int j = 0; j < Y.cols(); j++) {
            ret(i, j) = X(i, j) ^ Y(i, j);
        }
    }
    return ret;
}

/**
 * @brief turn GF to Binary, expanding log2(GF) times in cols
 *          e.g. 6 -> 01010000
 * @param X
 * @param GF
 * @return Eigen::MatrixXi
 */
Eigen::MatrixXi NB2Bin(const Eigen::MatrixXi& X, const int GF) {
    int rate = log2(GF);  // ret will be <X.rows(), X.cols() * rate>
    assert(pow(2, rate) == GF);
    Eigen::MatrixXi ret(X.rows(), X.cols() * rate);

    for (size_t i = 0; i < X.rows(); i++) {
        for (size_t j = 0; j < X.cols(); j++) {
            for (size_t r = 0; r < rate; r++) {
                ret(i, j * rate + r) = (X(i, j) >> r) & 1;  // 6 -> 01010000
            }
        }
    }

    return ret;
}

/**
 * @brief turn Binary to GF, shrinking log2(GF) times in cols
 *          e.g. 01010000 -> 6
 * @param X
 * @param GF
 * @return Eigen::MatrixXi
 */
Eigen::MatrixXi Bin2GF(const Eigen::MatrixXi& X, const int GF) {
    int rate = log2(GF);  // ret will be <X.rows(), X.cols() / rate>
    assert(X.cols() % rate == 0);
    Eigen::MatrixXi ret(X.rows(), X.cols() / rate);

    for (int i = 0; i < ret.rows(); i++) {
        for (int j = 0; j < ret.cols(); j++) {
            int ret_ij = 0;
            for (int r = rate - 1; r >= 0; r--) {
                ret_ij =
                    (ret_ij << 1) + X(i, j * rate + r);  // 6 -> 01010000
            }
            ret(i, j) = ret_ij;
        }
    }

    return ret;
}