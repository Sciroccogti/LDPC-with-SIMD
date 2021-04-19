#include "Channel/Channel.hpp"

/**
 * @brief add AWGN noise to QPSK signal
 *
 * @param origin original signal
 * @param snr
 * @param Q Q of QPSK // TODO:?
 * @return Eigen::RowVectorXd noised signal
 */
Eigen::RowVectorXd AWGN(const Eigen::RowVectorXd& origin, const double snr,
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

/**
 * @brief - 2x / sigma^2 TODO: seems to be - (2x-1) / sigma^2
 *
 * @param x
 * @param snr
 * @return double
 */
double LLR_AWGN(const double x, const double snr) {
    // 4 * snr = 2 / sigma^2
    return -4 * snr * x;
}

/**
 * @brief should be according to Bin2GF()
 *
 * @param X
 * @param GF
 * @param snr
 * @return Eigen::MatrixXd <GF - 1, X.cols() / rate>
 */
Eigen::MatrixXd LLR_BinAWGN2GF(const Eigen::RowVectorXd& X, const int GF,
                               const double snr) {
    int rate = log2(GF);  // ret will be <GF - 1, X.cols() / rate>
    assert(X.cols() % rate == 0);

    Eigen::MatrixXd ret(GF - 1, X.cols() / rate);
    for (int j = 0; j < ret.cols(); j++) {
        // q is 1 ~ GF -1, but stored at 0~GF-2, cause no LLR for 0
        for (int q = 1; q < GF; q++) {
            double ret_qj = 0;
            for (int r = rate - 1; r >= 0; r--) {
                // whether this bit has contribute to the q
                if (q & (1 << r)) {
                    ret_qj +=-X(j * rate + r);// LLR_AWGN(X(j * rate + r), snr);  // 01010000 -> 6
                }
            }
            ret(q - 1, j) = ret_qj;
        }
    }
    return ret;
}
