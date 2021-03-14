#include "Channel/Channel.hpp"

#include <random>

/**
 * @brief add AWGN noise to QPSK signal
 *
 * @param origin original signal
 * @param snr
 * @param Q Q of QPSK // TODO:?
 * @return Eigen::RowVectorXd noised signal
 */
Eigen::RowVectorXd AWGN(const Eigen::RowVectorXd &origin, const double snr,
                        const int Q) {
    // github.com/aff3ct/aff3ct/blob/master/src/Factory/Tools/Noise/Noise.cpp#L150
    double sigma = sqrt(1 / (2 * snr * log2(Q)));
    static std::default_random_engine engine(time(0));
    static std::normal_distribution<double> normal(0, sigma);
    Eigen::RowVectorXd ret = origin;
    for (int i = 0; i < ret.size(); i++) {
        ret[i] += normal(engine);
    }
    return ret;
}