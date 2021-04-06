/**
 * @file MatrixMath.hpp
 * @author Yifan Zhang (scirocco_gti@yeah.net)
 * @brief
 * @date 2020-10-31 10:57:44
 * @modified: 2021-04-07 00:50:25
 */

#ifndef MATRIXMATH_HPP
#define MATRIXMATH_HPP

#include <Eigen/Eigen>

using Positions_pair_vector = std::vector<std::pair<size_t, size_t>>;
using Container = std::vector<std::vector<int>>;

void swap_rows(Eigen::MatrixXi& mat, size_t idx1, size_t idx2);
int argmax(const Eigen::MatrixXi a);
Eigen::MatrixXi transform_H_to_G(const Eigen::MatrixXi& H);
Eigen::MatrixXi transform_H_to_G(const Eigen::SparseMatrix<int>& H);
Eigen::MatrixXi transform_H_to_G_sys(const Eigen::MatrixXi& H);
Eigen::MatrixXi transform_H_to_G_sys(const Eigen::SparseMatrix<int>& H);
Eigen::MatrixXi gaussjordan(Eigen::MatrixXi& X);
void gausselimination(Eigen::MatrixXi& X, Eigen::RowVectorXi& b);
Eigen::MatrixXi binaryproduct(const Eigen::MatrixXi& X,
                              const Eigen::MatrixXi& Y);
Eigen::MatrixXd cos(const Eigen::MatrixXd& X);
Eigen::MatrixXi abs(const Eigen::MatrixXi& X);
Eigen::RowVectorXi repeat(const Eigen::RowVectorXi& X, const int n);
Eigen::MatrixXd multiplyd(const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y);
Eigen::MatrixXi multiplyi(const Eigen::MatrixXi& X, const Eigen::MatrixXi& Y);
Eigen::RowVectorXd convolve(const Eigen::RowVectorXd& X,
                            const Eigen::RowVectorXd& Y);
Eigen::MatrixXi xori(const Eigen::MatrixXi& X, const Eigen::MatrixXi& Y);

Eigen::MatrixXi NBtransform_H_to_G(const Eigen::MatrixXi& H, const int GF);
Eigen::MatrixXi NBtransform_H_to_G(const Eigen::SparseMatrix<int>& H,
                                   const int GF);
Eigen::MatrixXi NBgaussjordan(Eigen::MatrixXi& X, const int GF);
int NBGauss(Eigen::MatrixXi &X, const int GF);
Eigen::MatrixXi NBproduct(const Eigen::MatrixXi& X, const Eigen::MatrixXi& Y,
                          const int GF);
Eigen::MatrixXi NBplus(const Eigen::MatrixXi& X, const Eigen::MatrixXi& Y,
                          const int GF);
Eigen::MatrixXi NB2Bin(const Eigen::MatrixXi& X, const int GF);
Eigen::MatrixXi Bin2GF(const Eigen::MatrixXi& X, const int GF);

#endif