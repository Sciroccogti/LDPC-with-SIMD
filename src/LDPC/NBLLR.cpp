#include "LDPC/NBLLR.hpp"

/**
 * @brief Construct a new NBLLR::NBLLR object
 *
 * @param d data
 * @param q Q, should in [1, GF]
 */
NBLLR::NBLLR(double d, int q) {
    data = d;
    Q = q;
}

NBLLR::~NBLLR() {}

double NBLLR::getData() const {
    return data;
}

int NBLLR::getQ() {
    return Q;
}

bool operator<(const NBLLR& llr1, const NBLLR& llr2) {
    return llr1.getData() < llr2.getData();
}

bool operator>(const NBLLR& llr1, const NBLLR& llr2) {
    return llr1.getData() > llr2.getData();
}