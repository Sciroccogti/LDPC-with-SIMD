#include "Channel/Channel.hpp"

/**
 * @brief add AWGN noise to QPSK signal
 *
 * @param origin original signal
 * @param snr
 * @param Q Q of QPSK // TODO:?
 * @return Eigen::RowVectorXf noised signal
 */
Eigen::RowVectorXf AWGN(const Eigen::RowVectorXf& origin, const float snr,
                        const int Q, std::default_random_engine engine) {
    // https://github.com/aff3ct/aff3ct/blob/master/src/Factory/Tools/Noise/Noise.cpp#L150
    float sigma = sqrt(1 / (2 * snr * log2(Q)));
    // https://github.com/aff3ct/aff3ct/blob/master/src/Tools/Algo/Draw_generator/Gaussian_noise_generator/Standard/Gaussian_noise_generator_std.cpp#L28
    std::normal_distribution<float> normal(0, sigma);
    Eigen::RowVectorXf ret = origin;
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
 * @return float
 */
float LLR_AWGN(const float x, const float snr) {
    // 4 * snr = 2 / sigma^2
    return -4 * snr * x;
}

/**
 * @brief should be according to Bin2GF()
 *
 * @param X
 * @param GF
 * @param snr
 * @return Eigen::MatrixXf <GF, X.cols() / rate>
 */
Eigen::MatrixXf LLR_BinAWGN2GF(const Eigen::RowVectorXf& X, const int GF,
                               const float snr) {
    int rate = log2(GF);  // ret will be <GF, X.cols() / rate>
    assert(X.cols() % rate == 0);

    Eigen::MatrixXf ret(GF, X.cols() / rate);
    for (int j = 0; j < ret.cols(); j++) {
        ret(0, j) = 0;  // LLR for 0 set to 0, 'cause log(1)=0
        for (int q = 1; q < GF; q++) {
            float ret_qj = 0;
            for (int r = rate - 1; r >= 0; r--) {
                // whether this bit has contribute to the q
                if (q & (1 << r)) {
                    ret_qj += LLR_AWGN(X(j * rate + r), snr);  // 01010000 -> 6
                }
                // ret_qj += - LLR_AWGN(X(j * rate + r) - (q & (1 << r)), snr) / 8;
            }
            ret(q, j) = ret_qj;
        }
    }
    return ret / (float)rate / 2.0;
}
