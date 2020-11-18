#include "MatrixMath/MatrixMath.hpp"

void swap_columns(Eigen::MatrixXi& mat, size_t idx1, size_t idx2) {
    auto n_row = mat.rows();
    std::vector<int> tmp(n_row);
    for (size_t l = 0; l < n_row; l++) tmp[l] = mat(l, idx1);
    for (size_t l = 0; l < n_row; l++) mat(l, idx1) = mat(l, idx2);
    for (size_t l = 0; l < n_row; l++) mat(l, idx2) = tmp[l];
}

void swap_rows(Eigen::MatrixXi& mat, size_t idx1, size_t idx2) {
    Eigen::RowVectorX<int> tmp1 = mat.row(idx1), tmp2 = mat.row(idx2);
    mat.row(idx1) = tmp1;
    mat.row(idx2) = tmp2;
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
 * @brief for systematic code
 *
 * @param H
 * @return Eigen::MatrixXi
 */
Eigen::MatrixXi transform_H_to_G(const Eigen::MatrixXi& H) {
    Eigen::MatrixXi G = H;
    int M = H.rows();
    int N = H.cols();
    int K = N - M;
    Positions_pair_vector swapped_cols = form_diagonal(G, 0, 2);
    form_identity(G);

    // erase the just created M*M identity in the left part of H and add the K*K
    // G.block(0, N - )
    G = resize_topright(G, N, K);
    // identity below
    for (auto i = M; i < N; i++) {  // Add rising diagonal identity at the end
        G(i, i - M) = 1;
    }

    // G is now VERTICAL

    // Re-organization: get G
    for (size_t l = swapped_cols.size(); l > 0; l--) {
        // std::swap(G[swapped_cols[l - 1].first], G[swapped_cols[l -
        // 1].second]);
        swap_rows(G, swapped_cols[l - 1].first, swapped_cols[l - 1].second);
    }

    return G;
}

/**
 * @brief LU decomp
 *
 * @param H
 * @return Eigen::MatrixXi
 */
Eigen::MatrixXi transform_H_to_G_LU(const Eigen::MatrixXi& H) {
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

Eigen::MatrixXi transform_H_to_G_LU(const Eigen::SparseMatrix<int>& H) {
    Eigen::MatrixXi Hdense(H);
    return transform_H_to_G_LU(Hdense);
}

Eigen::MatrixXi transform_H_to_G(const Eigen::SparseMatrix<int>& H) {
    Eigen::MatrixXi Hdense(H);
    return transform_H_to_G(Hdense);
}

/**
 * @brief turn a matrix to Top Left diagonal
 *
 * @param mat input matrix
 * @param direction 0: TOP LEFT; 1: BOTTOM LEFT
 * @param GFq the Galois field the mat takes, default to be 2: "1" = 0 or 1
 * @return Positions_pair_vector swapped_cols
 */
Positions_pair_vector form_diagonal(Eigen::MatrixXi& mat, int direction = 0,
                                    int GFq = 2) {
    int n_row = mat.rows();
    int n_col = mat.cols();

    Positions_pair_vector swapped_cols;

    switch (direction) {
        case 0:
            for (size_t i = 0; i < n_row; i++) {
                bool found = mat(i, i);

                if (!found) {
                    // try to find an other row which as a 1 in column i
                    for (size_t j = i + 1; j < n_row; j++)
                        if (mat(j, i)) {
                            // std::swap(mat.row(i), mat.row(j));
                            swap_rows(mat, i, j);
                            found = true;
                            break;
                        }
                    // no other row after (i+1) of the same column i with a 1
                    if (!found) {
                        // find an other column which is good on row i
                        for (auto j = i + 1; j < n_col; j++) {
                            if (mat(i, j)) {
                                swapped_cols.push_back(std::make_pair(i, j));

                                swap_columns(mat, i, j);

                                found = true;
                                break;
                            }
                        }
                    }
                }

                if (found) {
                    // there is a 1 on row i of the column i
                    // then remove any 1 of the column i from the row (i+1)
                    for (auto j = i + 1; j < n_row; j++)
                        if (mat(j, i))
                            std::transform(mat.row(i).begin() + i,
                                           mat.row(i).end(),  // ref
                                           mat.row(j).begin() + i,
                                           mat.row(j).begin() + i,
                                           std::not_equal_to<int>());
                } else {
                    // the row is the null vector then delete it
                    removeRow(mat, i);
                    i--;
                    n_row--;
                }
            }
            break;
        case 1:

            // for (size_t i = 0; i < n_col / 2; i++) {
            //     swap_columns(mat, i, n_col - i - 1);
            // }
            for (size_t i = n_row; i > 0; i--) {
                auto ref_row = i - 1;
                auto ref_col = n_row - ref_row - 1;
                bool found = mat(ref_row, ref_col);
                // std::cout << mat << "\n\n";
                if (!found) {
                    // try to find an other row which as a 1 in column ref_col
                    for (auto j = ref_row; j > 0; j--) {
                        auto tested_row = j - 1;
                        if (mat(tested_row, ref_col)) {
                            // std::swap(mat.row(ref_row), mat.row(tested_row));
                            swap_rows(mat, ref_row, tested_row);
                            found = true;
                            break;
                        }
                    }

                    if (!found)  // no other row before (ref_row-1) of the same
                                 // column ref_col with a 1
                    {
                        for (auto j = ref_col + 1; j < n_col;
                             j++)  // find an other column which is good on row
                                   // ref_row
                        {
                            if (mat(ref_row, j)) {
                                swapped_cols.push_back(
                                    std::make_pair(ref_col, j));

                                swap_columns(mat, ref_col, j);

                                found = true;
                                break;
                            }
                        }
                    }
                }

                if (found) {
                    // there is a 1 on row ref_row of the column ref_col
                    // then remove any 1 of the column ref_col from the row
                    // (ref_row-1)
                    for (auto j = ref_row; j > 0; j--) {
                        auto tested_row = j - 1;
                        if (mat(tested_row, ref_col))
                            std::transform(
                                mat.row(ref_row).begin() + ref_col,
                                mat.row(ref_row).end(),  // ref
                                mat.row(tested_row).begin() + ref_col,
                                mat.row(tested_row).begin() + ref_col,
                                std::not_equal_to<int>());
                    }
                } else {
                    // the row is the null vector then delete it
                    // mat.erase_row(ref_row);
                    removeRow(mat, ref_row);
                    n_row--;
                }
            }
            for (size_t i = 0; i < n_col / 2; i++) {
                swap_columns(mat, i, n_col - i - 1);
            }

            // std::cout << mat << "\n\n";

            break;
    }
    return swapped_cols;
}

void form_identity(Eigen::MatrixXi& mat) {
    int n_row = mat.rows();
    int n_col = mat.cols();
    int diff = n_col - n_row;
    for (size_t c = n_row - 1; c > 0; c--) {
        size_t ref_row = c;
        for (size_t r = c; r > 0; r--)
            if (mat(r - 1, c))
                std::transform(mat.row(ref_row).begin() + c,
                               mat.row(ref_row).end(),  // ref
                               mat.row(r - 1).begin() + c,
                               mat.row(r - 1).begin() + c,
                               std::not_equal_to<int>());
    }
}
