#include "LDPC/NBLLR.hpp"

/**
 * @brief Construct a new NBLLR::NBLLR object
 *
 * @param d data
 * @param q Q, should in [1, GF]
 */
NBLLR::NBLLR(float d, int q) {
    data = d;
    Q = q;
}

NBLLR::~NBLLR() {}

float NBLLR::getLLR() const {
    return data;
}

int NBLLR::getQ() {
    return Q;
}

bool operator<(const NBLLR& llr1, const NBLLR& llr2) {
    return llr1.getLLR() < llr2.getLLR();
}

bool operator>(const NBLLR& llr1, const NBLLR& llr2) {
    return llr1.getLLR() > llr2.getLLR();
}