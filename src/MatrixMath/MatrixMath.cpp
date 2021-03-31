#include "MatrixMath/MatrixMath.hpp"

#include <iostream>
void swap_columns(Eigen::MatrixXi& mat, size_t idx1, size_t idx2) {
    auto n_row = mat.rows();
    std::vector<int> tmp(n_row);
    for (size_t l = 0; l < n_row; l++) tmp[l] = mat(l, idx1);
    for (size_t l = 0; l < n_row; l++) mat(l, idx1) = mat(l, idx2);
    for (size_t l = 0; l < n_row; l++) mat(l, idx2) = tmp[l];
}

void swap_rows(Eigen::MatrixXi& mat, size_t idx1, size_t idx2) {
    auto n_col = mat.cols();
    std::vector<int> tmp(n_col);
    for (size_t l = 0; l < n_col; l++) tmp[l] = mat(idx1, l);
    for (size_t l = 0; l < n_col; l++) mat(idx1, l) = mat(idx2, l);
    for (size_t l = 0; l < n_col; l++) mat(idx2, l) = tmp[l];
}

/**
 * @brief https://stackoverflow.com/a/21068014/12851294
 *
 * @param matrix
 * @param rowToRemove
 */
void removeRow(Eigen::MatrixXi& matrix, unsigned int rowToRemove) {
    unsigned int numRows = matrix.rows() - 1;
    unsigned int numCols = matrix.cols();

    if (rowToRemove < numRows)
        matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) =
            matrix.block(rowToRemove + 1, 0, numRows - rowToRemove, numCols);

    matrix.conservativeResize(numRows, numCols);
}

/**
 * @brief resize while keep the TOP RIGHT still A|B -> B/0
 *
 * @param src
 * @param n_rows
 * @param n_cols
 * @return Eigen::MatrixXi
 */
Eigen::MatrixXi resize_topright(const Eigen::MatrixXi& src, const size_t n_rows,
                                const size_t n_cols) {
    Eigen::MatrixXi mat(n_rows, n_cols);
    mat.fill(0);

    int diff = n_cols - src.cols();
    size_t min_rows = n_rows <= src.rows() ? n_rows : src.rows();
    size_t min_cols = n_cols <= src.cols() ? n_cols : src.cols();

    if (diff > 0) {
        mat.block(0, diff, min_rows, min_cols) =
            src.block(0, 0, min_rows, min_cols);
    } else {
        mat.block(0, 0, min_rows, min_cols) =
            src.block(0, -diff, min_rows, min_cols);
    }

    return mat;
}

/**
 * @brief Returns the index of the nonzero value
 *
 * @param a : the input matrix
 * @return int
 */
int argmax(const Eigen::MatrixXi a) {
    for (size_t i = 0; i < a.rows(); i++) {
        for (size_t j = 0; j < a.cols(); j++) {
            if (a(i, j) != 0) {
                return i * a.cols() + j;
            }
        }
    }
    return 0;
}

/**
 * @brief https://github.com/hichamjanati/pyldpc/blob/master/pyldpc/code.py#L58
 *
 * @param H
 * @return Eigen::MatrixXi
 */
Eigen::MatrixXi transform_H_to_G(const Eigen::MatrixXi& H) {
    int n_code = H.cols();

    // DOUBLE GAUSS-JORDAN:
    Eigen::MatrixXi Href_colonnes = H.transpose();
    Eigen::MatrixXi Q = gaussjordan(Href_colonnes).transpose();
    std::cout << "Q:\n" << Q << std::endl;
    std::cout << "Q * H:\n" << H * Q.transpose() << std::endl;
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

/**
 * @brief Compute a coding matrix G in systematic format with an identity block.
 *
 * @param H
 * @return Eigen::MatrixXi
 */
Eigen::MatrixXi transform_H_to_G_sys(const Eigen::MatrixXi& H) {
    int M = H.rows();
    int N = H.cols();
    int K = N - M;
    Eigen::MatrixXi Hp = resize_topright(H.transpose(), M, M);
    // form_diagonal(Hp, 1, 2);
    Eigen::MatrixXf Hpinvf = Hp.cast<float>().inverse();
    Eigen::MatrixXf Hpinv2 = Hpinvf / -Hpinvf.minCoeff();

    Eigen::MatrixXi Hpinv =
        Hpinv2.unaryExpr([](const float x) { return (int)(abs(x) + 0.5) % 2; });

    Eigen::MatrixXi Hs = resize_topright(H, M, K).transpose();
    Eigen::MatrixXi GH =
        resize_topright((Hs * Hpinv.cast<int>()).transpose(), N, K);
    for (size_t r = 0; r < K; r++) {
        GH(r + M, r) = 1;
    }

    return GH.transpose();
}

Eigen::MatrixXi transform_H_to_G_sys(const Eigen::SparseMatrix<int>& H) {
    Eigen::MatrixXi Hdense(H);
    return transform_H_to_G_sys(Hdense);
}

Eigen::MatrixXi transform_H_to_G(const Eigen::SparseMatrix<int>& H) {
    Eigen::MatrixXi Hdense(H);
    return transform_H_to_G(Hdense);
}

/**
 * @brief Compute the binary row reduced echelon form of X.
 *        That is, turning X into up-right triangle.
 *        https://github.com/hichamjanati/pyldpc/blob/master/pyldpc/utils.py#L38
 * @param X : will be reused to return result!
 * @return Eigen::MatrixXi: the inverse transform
 */
Eigen::MatrixXi gaussjordan(Eigen::MatrixXi& X) {
    size_t m = X.rows(), n = X.cols();
    Eigen::MatrixXi P = Eigen::MatrixXi::Identity(m, m);

    int pivot_old = -1;
    for (size_t j = 0; j < n; j++) {
        Eigen::MatrixXi filtre_down =
            X.block(pivot_old + 1, j, m - pivot_old - 1, 1);

        // The first 1 found in filtre_down.
        // `argmax()` returns the index of the first 1.
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
                    P.row(i) =
                        (P.row(i) - P.row(pivot_old))
                            .unaryExpr([](const int x) { return abs(x); });
                    X.row(i) =
                        (X.row(i) - X.row(pivot_old))
                            .unaryExpr([](const int x) { return abs(x); });
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
 * https://github.com/hichamjanati/pyldpc/blob/master/pyldpc/utils.py#L161
 *
 * @param X : will be reused to return result!
 * @param b : will be reused to return result!
 */
void gausselimination(Eigen::MatrixXi& X, Eigen::RowVectorXi& b) {
    int N = X.rows(), K = X.cols();
    for (int i = 0; i < K; i++) {
        int count = 0;
        int minPivot = 0;

        for (int j = i; j < N; j++) {
            if (X(j, i)) {
                count++;
                if (j > minPivot) {
                    minPivot = j;
                }
            }
        }
        if (count == 0) {
            continue;
        }

        if (minPivot != i) {
            swap_rows(X, i, minPivot);
            int tmp = b[i];
            b[i] = b[minPivot];
            b[minPivot] = tmp;
        }

        for (int j = i + 1; j < N; j++) {
            if (X(j, i)) {
                X.row(j) = abs(X.row(j) - X.row(i));
                b[j] = abs(b[j] - b[i]);
            }
        }
    }
    return;
}

Eigen::MatrixXi binaryproduct(const Eigen::MatrixXi& X,
                              const Eigen::MatrixXi& Y) {
    assert(X.cols() == Y.rows());

    Eigen::MatrixXi ret = X * Y;

    return ret.unaryExpr([](const int x) { return x % 2; });
}

Eigen::MatrixXd cos(const Eigen::MatrixXd& X) {
    Eigen::MatrixXd ret = X;

    return ret.unaryExpr([](const double x) { return cos(x); });
}

Eigen::MatrixXi abs(const Eigen::MatrixXi& X) {
    Eigen::MatrixXi ret = X;

    return ret.unaryExpr([](const int x) { return abs(x); });
}

/**
 * @brief Repeat the input vector for n times.
 * @example X = [1, 2, 3, 4];
 *          Y = repeat(X, 3); // [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4]
 * @param X
 * @param n : the times to repeat
 * @return Eigen::MatrixXi
 */
Eigen::RowVectorXi repeat(const Eigen::RowVectorXi& X, const int n) {
    Eigen::RowVectorXi ret(X.size() * n);
    for (size_t i = 0; i < X.size(); i++) {
        for (size_t j = 0; j < n; j++) {
            ret(i * n + j) = X(i);
        }
    }

    return ret;
}

Eigen::MatrixXd multiplyd(const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y) {
    assert(X.cols() == Y.cols() && X.rows() == Y.rows());

    Eigen::MatrixXd ret(X.rows(), X.cols());
    for (size_t i = 0; i < X.rows(); i++) {
        for (size_t j = 0; j < X.cols(); j++) {
            ret(i, j) = X(i, j) * Y(i, j);
        }
    }
    return ret;
}

Eigen::MatrixXi multiplyi(const Eigen::MatrixXi& X, const Eigen::MatrixXi& Y) {
    assert(X.cols() == Y.cols() && X.rows() == Y.rows());

    Eigen::MatrixXi ret(X.rows(), X.cols());
    for (size_t i = 0; i < X.rows(); i++) {
        for (size_t j = 0; j < X.cols(); j++) {
            ret(i, j) = X(i, j) * Y(i, j);
        }
    }
    return ret;
}

/**
 * @brief https://toto-share.com/2011/11/cc-convolution-source-code/
 *
 * @param X
 * @param Y
 * @return Eigen::RowVectorXd
 */
Eigen::RowVectorXd convolve(const Eigen::RowVectorXd& X,
                            const Eigen::RowVectorXd& Y) {
    int nconv = X.size() + Y.size() - 1;
    Eigen::RowVectorXd ret = Eigen::RowVectorXd::Ones(nconv);

    for (size_t i = 0; i < nconv; i++) {
        size_t current_i = i;
        double tmp = 0;
        for (size_t j = 0; j < Y.size(); j++) {
            if (current_i >= 0 && current_i < X.size()) {
                tmp += (X(current_i) * Y(j));
            }
            current_i--;
        }
        ret(i) = tmp;
    }

    return ret;
}

/**
 * @brief
 * 
 * @param X 
 * @param Y 
 * @return Eigen::MatrixXi 
 */
Eigen::MatrixXi xori(const Eigen::MatrixXi& X, const Eigen::MatrixXi& Y) {
    assert(X.cols() == Y.cols() && X.rows() == Y.rows());

    Eigen::MatrixXi ret(X.rows(), X.cols());
    for (size_t i = 0; i < X.rows(); i++) {
        for (size_t j = 0; j < X.cols(); j++) {
            // printf("%d, %d, %d\n", X(i, j), Y(i, j), X(i, j) ^ Y(i, j));
            ret(i, j) = X(i, j) ^ Y(i, j);
        }
    }
    return ret;
}
