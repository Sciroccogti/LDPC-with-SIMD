/**
 * @file MatrixMath.hpp
 * @author Yifan Zhang (scirocco_gti@yeah.net)
 * @brief transplanted from aff3ct
 * https://github.com/aff3ct/aff3ct/tree/master/src/Tools/Code/LDPC/Matrix_handler
 * @date 2020-10-31 10:57:44
 * @modified: 2020-10-31 14:02:56
 */

#ifndef MATRIXMATH_HPP
#define MATRIXMATH_HPP

#include <Eigen/Eigen>

using Positions_pair_vector = std::vector<std::pair<size_t, size_t>>;
using Container = std::vector<std::vector<int>>;

Eigen::MatrixXi transform_H_to_G(const Eigen::MatrixXi& H);
Eigen::MatrixXi transform_H_to_G(const Eigen::SparseMatrix<int>& H);
Eigen::MatrixXi transform_H_to_G_LU(const Eigen::MatrixXi& H);
Eigen::MatrixXi transform_H_to_G_LU(const Eigen::SparseMatrix<int>& H);
Positions_pair_vector form_diagonal(Eigen::MatrixXi&, int, int);
void form_identity(Eigen::MatrixXi&);
Eigen::MatrixXi bgemmt(const Eigen::MatrixXi& A, const Eigen::MatrixXi& tB);

#endif