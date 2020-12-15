/**
 * @file MatrixMath.hpp
 * @author Yifan Zhang (scirocco_gti@yeah.net)
 * @brief
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
Eigen::MatrixXi transform_H_to_G_sys(const Eigen::MatrixXi& H);
Eigen::MatrixXi transform_H_to_G_sys(const Eigen::SparseMatrix<int>& H);
Eigen::MatrixXi gaussjordan(Eigen::MatrixXi& X);
Eigen::MatrixXi binaryproduct(const Eigen::MatrixXi& X,
                              const Eigen::MatrixXi& Y);
Eigen::MatrixXd cos(const Eigen::MatrixXd& X);
Eigen::RowVectorXi repeat(const Eigen::RowVectorXi& X, const int n);
Eigen::MatrixXd multiply(const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y);

#endif