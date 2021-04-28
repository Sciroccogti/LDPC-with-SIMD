#include <iostream>
#include <queue>

#include "GF.hpp"
#include "LDPC/Tanner.hpp"
#include "MatrixMath/MatrixMath.hpp"

using namespace std;

/**
 * @brief Construct a new NBNode:: NBNode object
 *  inValues are initialized as 0
 * @param d degree
 */
NBNode::NBNode(int d, int gf, int nmax) {
    degree = d;
    GF = gf;
    n_max = nmax;
    assert(n_max <= GF);
    // Nodes_.reserve(d);

    inCount = 0;
    // init inValues
    for (int i = 0; i < degree; i++) {
        std::vector<float> dataQ(GF);
        inValuesQ_.push_back(dataQ);
    }

    n_maxCount = 0;
    // init inn_maxValue
    for (int i = 0; i < degree; i++) {
        n_maxValue_.push_back(vector<NBLLR>());
    }
}

NBNode::~NBNode() {}

// return true if size == degree
bool NBNode::isReady() {
    return NBNodes_.size() == degree;
}

void NBNode::setInValue(std::vector<float> dataQ) {
    inValuesQ_[inCount] = dataQ;
    inCount++;
    if (inCount >= degree) {
        inCount = 0;
    }
}

void NBNode::setinn_maxValue(std::vector<NBLLR> vn_max) {
    n_maxValue_[n_maxCount] = vn_max;
    n_maxCount++;
    if (n_maxCount >= degree) {
        n_maxCount = 0;
    }
}

/**
 * @brief Construct a new NBVNode::NBVNode object
 *
 * @param d degree
 * @param v_ GF initial values
 * @param gf
 */
NBVNode::NBVNode(int d, Eigen::RowVectorXf vQ, const int gf, const int nmax)
    : NBNode(d, gf, nmax) {
    vector<float> stdvQ(vQ.data(), vQ.data() + vQ.size());
    valueQ = stdvQ;
    LLRQ = stdvQ;
    assert(valueQ.size() == GF);
}

/**
 * @brief
 *
 * @param n NBCNode
 * @param h value in H_mat
 */
void NBVNode::Link(NBNode* n, int h) {
    assert(!n->isVN());  // assure n is CN
    assert(NBNodes_.size() + 1 <= degree);
    NBNodes_.push_back(n);
    n->Link(this, h);  // register self to n
}

void NBVNode::Update(int mode) {
    valueQ = LLRQ;
    // update own value
    for (int i = 0; i < degree; i++) {
        for (int q = 0; q < GF; q++) {
            valueQ[q] += inValuesQ_[i][q];
        }
    }

    // update outputs
    for (int i = 0; i < degree; i++) {
        NBNodes_[i]->setInValue(valueQ);

        // https://docs.microsoft.com/zh-cn/troubleshoot/cpp/stl-priority-queue-class-custom-type
        // less: top is maximum
        // greater: top is minimum
        priority_queue<NBLLR, vector<NBLLR>, greater<vector<NBLLR>::value_type>>
            pqLLR;  // size should be GF
        for (int q = 0; q < GF; q++) {
            float cur_LLR = valueQ[q] - inValuesQ_[i][q];

            if (pqLLR.size() < n_max) {
                // havn't get enough
                pqLLR.push(NBLLR(cur_LLR, q));
            } else if (cur_LLR > pqLLR.top().getLLR()) {
                // wont push unless cur_LLR is large enough
                pqLLR.pop();  // remove the minimum
                pqLLR.push(NBLLR(cur_LLR, q));
            }
        }

        std::vector<NBLLR> LLR_(n_max);  // size should be n_max
        // pqLLR is min first, so should be reversed
        for (int j = n_max - 1; j >= 0; j--) {
            LLR_[j] = pqLLR.top();
            pqLLR.pop();
        }

        NBNodes_[i]->setinn_maxValue(LLR_);
        // priority_queue<NBLLR, vector<NBLLR>, less<vector<NBLLR>::value_type>>
        //     pqcopy(pqLLR);
        // printf("pq:\n");
        // for (size_t j = 0; j < n_maxLLR_.size(); j++) {
        //     printf("%.2lf:%d   ", n_maxLLR_[j].getLLR(),
        //     n_maxLLR_[j].getQ());
        // }
        // printf("\n");
    }
}

/**
 * @brief directly fetch GF estimate
 *
 * @return int
 */
int NBVNode::getValue() {
    float max = valueQ[0];
    int ret = 0;
    for (int i = 1; i < GF; i++) {
        // printf("%.2lf ", valueQ[i]);
        if (max < valueQ[i]) {
            max = valueQ[i];
            ret = i;
        }
    }
    // printf(" -> %.0lf ", max);
    // if (max < 0) {  // output 0 as others are larger than 0's probability
    //     return 0;
    // }
    return ret;
}

bool NBVNode::isVN() {
    return true;
}

/**
 * @brief Construct a new NBCNode::NBCNode object
 *
 * @param d degree
 * @param f normalize factor
 * @param gf
 */
NBCNode::NBCNode(int d, float f, const int gf, const int nmax)
    : NBNode(d, gf, nmax) {
    assert(f <= 1);
    factor = f;
    Hrow_ = vector<uint8_t>(degree);
}

/**
 * @brief
 *
 * @param n NBVNode
 * @param h value in H_mat
 */
void NBCNode::Link(NBNode* n, int h) {
    Hrow_[NBNodes_.size()] = h;
    assert(NBNodes_.size() + 1 <= degree);
    NBNodes_.push_back(n);
}

void NBCNode::Update(int mode) {
    assert(mode >= BP_EMS && mode <= BP_QSPA);
    assert(inCount == 0);  // assure all inValue is changed

    // select current VN
    for (int cur_VN = 0; cur_VN < degree; cur_VN++) {
        switch (mode) {
            case BP_EMS: {
                // output to cur_VN
                vector<float> output(GF);
                vector<uint8_t> confset(degree - 1);
                int confsetCount = 0;  // current No. of confset
                int cur = 0;           // cursur for nconf_q_1
                // max LLR sum for 0 <= Q < GF
                vector<float> max(GF);

                // -1 if still in conf(q, 1), 1 if not reaches end, 0 if
                // reaches end
                int conf_q_1 = -1;
                do {
                    // for (auto&& i : confset) {
                    //     printf("%3d ", i);
                    // }
                    // printf("\n");
                    uint8_t prodsum = 0;  // calculate the Q of cur_VN
                    float sum = 0;        // sum of the LLR
                    int hasPassedcur_VN = 0;
                    // cursor among other VNs
                    for (size_t i = 0; i < degree; i++) {
                        if (i == cur_VN) {
                            hasPassedcur_VN = 1;
                            continue;
                        }

                        uint8_t Qi;

                        // conf_q_1 should choose among all Q, no only n_max
                        if (conf_q_1 == -1) {
                            // 0 should be the maxLLR
                            if (confset[i - hasPassedcur_VN] == 0) {
                                Qi = n_maxValue_[i][0].getQ();
                            } else {
                                Qi = confset[i - hasPassedcur_VN];
                            }

                            sum += inValuesQ_[i][Qi];
                        } else {  // conf(n_max, d-1)
                            Qi = n_maxValue_[i][confset[i - hasPassedcur_VN]]
                                     .getQ();
                            sum += n_maxValue_[i][confset[i - hasPassedcur_VN]]
                                       .getLLR();
                        }

                        uint8_t Hi = Hrow_[i];
                        // sum(Qi * Hi)
                        prodsum = GF_plus(GF_mul(Qi, Hi, GF), prodsum, GF);
                    }

                    uint8_t cur_VNQ = GF_div(prodsum, Hrow_[cur_VN], GF);

                    sum += inValuesQ_[cur_VN][cur_VNQ];
                    if (max[cur_VNQ] < sum) {
                        max[cur_VNQ] = sum;
                        output[cur_VNQ] = sum - inValuesQ_[cur_VN][cur_VNQ];
                    }
                    // until conf_q_1 = 0
                } while (conf_q_1 = getConfset(confset, confsetCount, cur));
                NBNodes_[cur_VN]->setInValue(output);
            } break;

            case BP_QSPA: {
                // Eigen::RowVectorXf prod = Eigen::RowVectorXf::Ones(GF);
                // Eigen::RowVectorXi sgn = Eigen::RowVectorXi::Ones(GF);
                // for (int j = 0; j < degree; j++) {
                //     if (j == cur_VN) {
                //         continue;  // skip current VN
                //     }
                //     for (int q = 0; q < GF; q++) {
                //         prod[q] *= tanh(fabs(inValuesQ_[j][q]) / 2);
                //         sgn[q] *= inValuesQ_[j][q] >= 0 ? 1 : -1;
                //     }
                // }

                // for (int q = 0; q < GF; q++) {
                //     prod[q] = prod[q] < 1
                //                   ? prod[q]
                //                   : 1.0 - numeric_limits<float>::epsilon();
                //     prod[q] = sgn[q] * 2 * fabs(atanh(prod[q]));
                // }

                // NBNodes_[cur_VN]->setInValue(prod);
            } break;
            default:
                break;
        }
    }
}

/**
 * @brief get next Confset.
 * We use conf(n_max, b) to present a Confset, where `n_max` means the n max LLR
 * to choose from, `b` means b VNs can have LLRs other than the top n_max LLR.
 * This version doesn't care about which is the cur_VN, 'cause confset has only
 * `degree - 1` elements, which has already skipped the cur_VN
 *
 * CAUTION: Q=0 means the maxLLR, not really Q = 0
 *
 * @param confset should have degree - 1 elements
 * @param confsetCount current No. of confset
 * @param cur cursur for nconf_q_1
 * @return -1 if still in conf(q, 1), 1 if not reaches end, 0 if reaches end
 */
int NBCNode::getConfset(vector<uint8_t>& confset, int& confsetCount, int& cur) {
    assert(confset.size() == degree - 1);
    // if confset is all zero, then reset
    bool isAllZero = true;
    for (uint8_t i : confset) {
        if (i) {
            isAllZero = false;
            break;
        }
    }
    if (isAllZero) {
        confsetCount = 0;
        cur = 0;
    }

    // number of conf(q, 1)
    const int nconf_q_1 = (degree - 1) * (GF - 1) + 1;

    // if still in conf(q, 1)
    if (confsetCount < nconf_q_1) {
        if (confset[cur] >= GF - 1) {
            confset[cur] = 0;  // reset confset[cur] to greatest one
            if (cur < degree - 2) {
                cur++;              // start to process next VN
                confset[cur] += 1;  // skip next 0 to avoid all zero confset
            } else {
                // turn to conf(n_max, d - 1)
                for (size_t i = 0; i < confset.size(); i++) {
                    confset[i] = 1;
                }
            }
        } else {
            confset[cur] += 1;  // move confset[cur] to smaller one
        }

    } else {
        // conf(n_max, d - 1)
        for (int i = 0; i < degree - 1; i++) {
            confset[i] += 1;  // move confset[i] to smaller one

            if (i == degree - 2 && confset[i] == n_max) {  // reaches end
                return 0;
            } else if (confset[i] >= n_max) {
                confset[i] = 0;
                // continue to modify next VN
            } else {
                break;  // don't modify next VN
            }
        }
    }

    confsetCount++;
    if (confsetCount >= nconf_q_1) {
        return 1;
    } else {
        return -1;
    }
}

bool NBCNode::isVN() {
    return false;
}
