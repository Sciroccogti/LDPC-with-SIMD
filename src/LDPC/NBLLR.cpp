#include "LDPC/NBLLR.hpp"

NBLLR::NBLLR() {
    data = 0;
    Q = 0;
}

/**
 * @brief Construct a new NBLLR::NBLLR object
 *
 * @param d data
 * @param q Q, should in [1, GF]
 */
NBLLR::NBLLR(float d, uint8_t q) {
    data = d;
    Q = q;
}

NBLLR::~NBLLR() {}

float NBLLR::getLLR() const {
    return data;
}

uint8_t NBLLR::getQ() {
    return Q;
}

bool operator<(const NBLLR& llr1, const NBLLR& llr2) {
    return llr1.getLLR() < llr2.getLLR();
}

bool operator>(const NBLLR& llr1, const NBLLR& llr2) {
    return llr1.getLLR() > llr2.getLLR();
}