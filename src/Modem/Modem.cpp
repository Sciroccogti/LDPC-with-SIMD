#include "Modem/Modem.hpp"

/**
 * @brief Construct a new Modem:: Modem object
 *
 * @param freqc carrier frequency
 * @param freqs sampling frequency
 * @param Timeb secs of each symbol
 * @param t type of Modem, choose from MODEM_BPSK, default to be MODEM_BPSK
 */
Modem::Modem(int freqc, int freqs, float Timeb, int t) {
    fc = freqc;
    fs = freqs;
    Tb = Timeb;
    type = t;
    L = Tb * fs;
}

Modem::~Modem() {}

/**
 * @brief https://github.com/DenisMedeiros/BPSKModDemod
 *
 * @param c : encoded code word
 * @return Eigen::RowVectorXf
 */
Eigen::RowVectorXf Modem::modulate(const Eigen::RowVectorXi &c) {
    // turn 0 1 sequence into 1 -1 sequence
    Eigen::RowVectorXi c_symbol = RetransBPSK(c);
    int length = L * c.size();
    Eigen::RowVectorXf ret;

    switch (type) {
        case MODEM_BPSK:

            // Generate the square wave where each Tb lasts L samples
            Eigen::RowVectorXf wave = repeat(c_symbol, L).cast<float>();

            // Generate the carrier
            Eigen::RowVectorXf t =
                Eigen::RowVectorXf::LinSpaced(length, 0, length);
            Eigen::RowVectorXf carrier = cos(2 * M_PI * fc * t);

            ret = multiplyd(carrier, wave);
            break;
    }

    return ret;
}

// TODO: #8 Real BPSK will introduce error
Eigen::RowVectorXf Modem::demodulate(const Eigen::RowVectorXf &x) {
    Eigen::RowVectorXf ret;
    switch (type) {
        case MODEM_BPSK:

            // // Generate the square wave where each Tb lasts L samples
            // Eigen::RowVectorXf wave = repeat(c_symbol, L).cast<float>();

            // Generate the carrier
            Eigen::RowVectorXf t =
                Eigen::RowVectorXf::LinSpaced(x.size(), 0, x.size());
            Eigen::RowVectorXf carrier = cos(2 * M_PI * fc * t);

            Eigen::RowVectorXf signald = multiplyd(carrier, x);

            signald = convolve(signald, Eigen::RowVectorXf::Ones(L));
            signald = signald.block(0, L - 1, 1, signald.size() - L);

            // Eigen::RowVectorXi signali = signald.unaryExpr(
            //     [](const float x) { return x > 0 ? 1 : -1; });

            // take one in every L values
            Eigen::RowVectorXf tmpRet((int)ceil(signald.size() / (float)L));
            for (size_t i = 0; i < signald.size(); i += L) {
                tmpRet(i / L) = signald(i + int(L / 2));
            }

            ret = tmpRet;
            break;
    }
    return ret;
}

int Modem::getL() {
    return L;
}

/**
 * @brief Compare two Matrix
 *
 * @param X
 * @param Y
 * @return int : error count
 */
int Compare(const Eigen::MatrixXi &X, const Eigen::MatrixXi &Y) {
    assert(X.size() == Y.size());
    int diffCount = 0;
    for (size_t i = 0; i < X.size(); i++) {
        int Xi = X(i), Yi = Y(i);
        assert(Xi == 0 || Xi == 1);
        assert(Yi == 0 || Yi == 1);

        if (Xi != Yi) {
            diffCount++;
        }
    }
    return diffCount;
}

/**
 * @brief Compare two Matrix
 *
 * @param X
 * @param Y
 * @return int : error count
 */
int Compare(const Eigen::MatrixXf &X, const Eigen::MatrixXf &Y) {
    assert(X.size() == Y.size());
    int diffCount = 0;
    for (size_t i = 0; i < X.size(); i++) {
        float Xi = X(i), Yi = Y(i);

        if (Xi != Yi) {
            diffCount++;
        }
    }
    return diffCount;
}

/**
 * @brief Compare two sequence, one is 0,1 sequence, latter one is -1,1 sequence
 *
 * @param X  0, 1 sequence
 * @param Y -1, 1 sequence
 * @return int : error count
 */
int CompareBPSK(const Eigen::RowVectorXi &X, const Eigen::RowVectorXi &Y) {
    assert(X.size() == Y.size());

    int diffCount = 0;
    for (size_t i = 0; i < X.size(); i++) {
        int Xi = X(i), Yi = Y(i);
        assert(Xi == 0 || Xi == 1);
        assert(Yi == -1 || Yi == 1);

        if (2 * Xi - 1 != Yi) {
            diffCount++;
        }
    }
    return diffCount;
}

/**
 * @brief Trans 1,-1 sequence to 0,1 sequence. 1=>0, -1=>1
 *
 * @param X
 * @return Eigen::MatrixXi
 */
Eigen::MatrixXi TransBPSK(const Eigen::MatrixXi &X) {
    Eigen::MatrixXi ret = X;

    return ret.unaryExpr([](const int x) {
        assert(x == 1 || x == -1);
        return (1 - x) / 2;
    });
}

/**
 * @brief Trans 0,1 sequence to 1,-1 sequence. 0=>1, 1=> -1
 *
 * @param X
 * @return Eigen::MatrixXi
 */
Eigen::MatrixXi RetransBPSK(const Eigen::MatrixXi &X) {
    Eigen::MatrixXi ret = X;

    return ret.unaryExpr([](const int x) {
        assert(x == 1 || x == 0);
        return 1 - 2 * x;
    });
}