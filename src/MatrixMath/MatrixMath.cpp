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
 * @brief resize while keep the top right still A|B -> B/0
 *
 * @param mat
 * @param n_rows
 * @param n_cols
 */
void self_resize(Eigen::MatrixXi& mat, const size_t n_rows,
                 const size_t n_cols) {
    // Container::resize(n_rows, std::vector<int>(n_cols, 0));
    //     auto n_erase = get_n_cols() - n_cols;
    //     for (size_t r = 0; r < n_rows; r++)
    //         (*this)[r].erase((*this)[r].begin(), (*this)[r].begin() +
    //         n_erase);
    auto diff = mat.cols() - n_cols;
    Eigen::MatrixXi src(mat);
    mat.resize(n_rows, n_cols);
    mat.block(0, 0, src.rows(), n_cols) =
        src.block(0, diff, src.rows(), n_cols);
    mat.block(src.rows(), 0, n_cols, n_cols).fill(0);
}

Eigen::MatrixXi transform_H_to_G(const Eigen::MatrixXi& H) {
    Eigen::MatrixXi G = H;
    int M = H.rows();
    int N = H.cols();
    int K = N - M;
    Positions_pair_vector swapped_cols = form_diagonal(G, 2);
    form_identity(G);

    // erase the just created M*M identity in the left part of H and add the K*K
    // G.block(0, N - )
    self_resize(G, N, K);
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

Eigen::MatrixXi transform_H_to_G(const Eigen::SparseMatrix<int>& H) {
    Eigen::MatrixXi Hdense(H);
    return transform_H_to_G(Hdense);
}

/**
 * @brief turn a matrix to Top Left diagonal
 *
 * @param mat input matrix
 * @param GFq the Galois field the mat takes, default to be 2: "1" = 0 or 1
 * @return Positions_pair_vector swapped_cols
 */
Positions_pair_vector form_diagonal(Eigen::MatrixXi mat, int GFq = 2) {
    int n_row = mat.rows();
    int n_col = mat.cols();

    Positions_pair_vector swapped_cols;

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

            if (!found)  // no other row after (i+1) of the same column i with a
                         // 1
            {
                for (auto j = i + 1; j < n_col;
                     j++)  // find an other column which is good on row i
                {
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
                    std::transform(
                        mat.row(i).begin() + i, mat.row(i).end(),  // ref
                        mat.row(j).begin() + i, mat.row(j).begin() + i,
                        std::not_equal_to<int>());
        } else {
            // the row is the null vector then delete it
            removeRow(mat, i);
            i--;
            n_row--;
        }
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
