/**
 * @file Channel.hpp
 * @author Sciroccogti (scirocco_gti@yeah.net)
 * @brief
 * @date 2021-03-11 21:03:04
 * @modified: 2021-04-12 01:17:46
 */
#ifndef CHANNEL_HPP
#define CHANNEL_HPP

#include <random>

#include "MatrixMath/MatrixMath.hpp"

Eigen::RowVectorXd AWGN(const Eigen::RowVectorXd& origin, const double snr,
                        const int Q, std::default_random_engine engine);

double LLR_AWGN(const double x, const double snr);
Eigen::MatrixXd LLR_BinAWGN2GF(const Eigen::RowVectorXd& X, const int GF,
                               const double snr);

#endif