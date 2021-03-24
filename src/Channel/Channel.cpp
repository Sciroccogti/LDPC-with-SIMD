#include "Channel/Channel.hpp"

/**
 * @brief add AWGN noise to QPSK signal
 *
 * @param origin original signal
 * @param snr
 * @param Q Q of QPSK // TODO:?
 * @return Eigen::RowVectorXd noised signal
 */
Eigen::RowVectorXd AWGN(const Eigen::RowVectorXd &origin, const double snr,
                        const int Q, std::default_random_engine engine) {
    // https://github.com/aff3ct/aff3ct/blob/master/src/Factory/Tools/Noise/Noise.cpp#L150
    double sigma = sqrt(1 / (2 * snr * log2(Q)));
    // https://github.com/aff3ct/aff3ct/blob/master/src/Tools/Algo/Draw_generator/Gaussian_noise_generator/Standard/Gaussian_noise_generator_std.cpp#L28
    std::normal_distribution<double> normal(0, sigma);
    Eigen::RowVectorXd ret = origin;
    for (int i = 0; i < ret.size(); i++) {
        ret[i] += normal(engine);
    }
    return ret;
}