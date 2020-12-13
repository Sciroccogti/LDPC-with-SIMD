#include "Modem/Modem.hpp"

/**
 * @brief Construct a new Modem:: Modem object
 *
 * @param freqc carrier frequency
 * @param freqs sampling frequency
 * @param Timeb secs of each symbol
 * @param t type of Modem, choose from MODEM_BPSK, default to be MODEM_BPSK
 */
Modem::Modem(int freqc, int freqs, double Timeb, int t) {
    fc = freqc;
    fs = freqs;
    Tb = Timeb;
    type = t;
    L = Tb * fs;
}

Modem::~Modem() {}

Eigen::RowVectorXd Modem::modulate(const Eigen::RowVectorXi &c) {
    // turn 0 1 sequence into -1 1 sequence
    Eigen::RowVectorXi c_symbol =
        c.unaryExpr([](const int x) { return 2 * x - 1; });
    int length = L * c.size();
    Eigen::RowVectorXd ret;

    switch (type) {
        case MODEM_BPSK:

            // Generate the square wave where each Tb lasts L samples
            Eigen::RowVectorXd wave = repeat(c_symbol, L).cast<double>();

            // Generate the carrier
            Eigen::RowVectorXd t =
                Eigen::RowVectorXd::LinSpaced(length, 0, length);
            Eigen::RowVectorXd carrier = cos(2 * M_PI * fc * t);

            ret = multiply(carrier, wave);
            break;
    }

    return ret;
}