/**
 * @file Channel.hpp
 * @author Sciroccogti (scirocco_gti@yeah.net)
 * @brief
 * @date 2021-03-11 21:03:04
 * @modified: 2021-04-27 13:13:32
 */
#ifndef CHANNEL_HPP
#define CHANNEL_HPP

#include <random>

#include "MatrixMath/MatrixMath.hpp"

Eigen::RowVectorXf AWGN(const Eigen::RowVectorXf& origin, const float snr,
                        const int Q, std::default_random_engine engine);

float LLR_AWGN(const float x, const float snr);
Eigen::MatrixXf LLR_BinAWGN2GF(const Eigen::RowVectorXf& X, const int GF,
                               const float snr);

#endif