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
        inValuesQ_.push_back(Eigen::RowVectorXd::Zero(GF));
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

void NBNode::setInValue(Eigen::RowVectorXd dataQ) {
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
NBVNode::NBVNode(int d, Eigen::RowVectorXd vQ, const int gf, const int nmax)
    : NBNode(d, gf, nmax) {
    valueQ = vQ;
    LLRQ = vQ;
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
        valueQ = inValuesQ_[i] + valueQ;
    }

    // update outputs
    for (int i = 0; i < degree; i++) {
        NBNodes_[i]->setInValue(valueQ);

        // https://docs.microsoft.com/zh-cn/troubleshoot/cpp/stl-priority-queue-class-custom-type
        // less: top is maximum
        priority_queue<NBLLR, vector<NBLLR>, less<vector<NBLLR>::value_type>>
            pqLLR;  // size should be GF
        for (int q = 0; q < GF; q++) {
            pqLLR.push(NBLLR(valueQ[q] - inValuesQ_[i][q], q));
        }

        std::vector<NBLLR> LLR_;  // size should be GF
        for (int j = 0; j < GF; j++) {
            LLR_.push_back(pqLLR.top());
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
    double max = valueQ[0];
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
NBCNode::NBCNode(int d, double f, const int gf, const int nmax)
    : NBNode(d, gf, nmax) {
    assert(f <= 1);
    factor = f;
    Hrow_ = Eigen::RowVectorXi::Zero(degree);
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

    for (int cur_VN = 0; cur_VN < degree; cur_VN++) {
        switch (mode) {
            case BP_EMS: { /*
                Eigen::RowVectorXd min = Eigen::RowVectorXd::Zero(GF);
                Eigen::RowVectorXi sgn = Eigen::RowVectorXi::Ones(GF);
                bool isFirst = true;
                for (int j = 0; j < degree; j++) {
                    if (j == i) {
                        continue;  // skip current VN
                    }
                    Eigen::RowVectorXd absJ(GF);
                    for (int q = 0; q < GF; q++) {
                        absJ[q] = fabs(inValuesQ_[j][q]);
                        if (!isFirst) {
                            min[q] = absJ[q] < min[q] ? absJ[q] : min[q];
                        }

                        sgn[q] *= inValuesQ_[j][q] >= 0 ? 1 : -1;
                    }

                    if (isFirst) {
                        isFirst = false;
                        min = absJ;
                    }
                }

                for (int q = 0; q < GF; q++) {
                    min[q] = sgn[q] * min[q] * factor;
                }

                NBNodes_[i]->setInValue(min);*/

                // select current VN
                // output to cur_VN
                Eigen::RowVectorXd output = Eigen::RowVectorXd::Zero(GF);
                Eigen::RowVectorXi confset =
                    Eigen::RowVectorXi::Zero(degree - 1);
                // max LLR sum for 0 <= Q < GF
                Eigen::RowVectorXd max = Eigen::RowVectorXd::Zero(GF);

                // -1 if still in conf(q, 1), 1 if not reaches end, 0 if
                // reaches end
                int conf_q_1 = -1;
                do {
                    // cout << confset << endl;
                    int prodsum = 0;  // calculate the Q of cur_VN
                    double sum = 0;   // sum of the LLR
                    int hasPassedcur_VN = 0;
                    // cursor among other VNs
                    for (size_t i = 0; i < degree; i++) {
                        if (i == cur_VN) {
                            hasPassedcur_VN = 1;
                            continue;
                        }

                        int Qi =
                            n_maxValue_[i][confset[i - hasPassedcur_VN]].getQ();
                        // conf_q_1 should choose among all Q, no only n_max
                        if (conf_q_1 == -1) {
                            sum += inValuesQ_[i][Qi];
                        } else {
                            sum += n_maxValue_[i][confset[i - hasPassedcur_VN]]
                                       .getLLR();
                        }

                        int Hi = Hrow_[i];
                        // sum(Qi * Hi)
                        prodsum = GF_plus(GF_mul(Qi, Hi, GF), prodsum, GF);
                    }

                    int cur_VNQ = GF_div(prodsum, Hrow_[cur_VN], GF);

                    sum += inValuesQ_[cur_VN][cur_VNQ];
                    if (max[cur_VNQ] < sum) {
                        max[cur_VNQ] = sum;
                        output[cur_VNQ] = sum - inValuesQ_[cur_VN][cur_VNQ];
                    }
                } while (conf_q_1 = getConfset(confset));  // until conf_q_1 = 0
                NBNodes_[cur_VN]->setInValue(output);

                int maxQ = 0;
                double maxoutput = output[0];
                for (size_t i = 1; i < output.size(); i++) {
                    if (output[i] > maxoutput) {
                        maxoutput = output[i];
                        maxQ = i;
                    }
                }
                // printf("%d: %d\n", cur_VN, maxQ);

            } break;
            case BP_QSPA: {
                Eigen::RowVectorXd prod = Eigen::RowVectorXd::Ones(GF);
                Eigen::RowVectorXi sgn = Eigen::RowVectorXi::Ones(GF);
                for (int j = 0; j < degree; j++) {
                    if (j == cur_VN) {
                        continue;  // skip current VN
                    }
                    for (int q = 0; q < GF; q++) {
                        prod[q] *= tanh(fabs(inValuesQ_[j][q]) / 2);
                        sgn[q] *= inValuesQ_[j][q] >= 0 ? 1 : -1;
                    }
                }

                for (int q = 0; q < GF; q++) {
                    prod[q] = prod[q] < 1
                                  ? prod[q]
                                  : 1.0 - numeric_limits<double>::epsilon();
                    prod[q] = sgn[q] * 2 * fabs(atanh(prod[q]));
                }

                NBNodes_[cur_VN]->setInValue(prod);
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
 * @param confset should have degree - 1 elements
 * @return -1 if still in conf(q, 1), 1 if not reaches end, 0 if reaches end
 */
int NBCNode::getConfset(Eigen::RowVectorXi& confset) {
    assert(confset.size() == degree - 1);
    static int confsetCount = 0;  // current No. of confset
    static int cur = 0;           // cursur for nconf_q_1
    // if confset is all zero, then reset
    if (!confset.any()) {
        confsetCount = 0;
        cur = 0;
    }
    // number of conf(q, 1)
    const int nconf_q_1 = (degree - 1) * (GF - 1) + 1;

    // if still in conf(q, 1)
    if (confsetCount < nconf_q_1) {
        confset[cur] += 1;  // move confset[cur] to smaller one

        if (confset[cur] >= GF && cur < degree - 2) {
            confset[cur] = 0;   // reset confset[cur] to greatest one
            cur++;              // start to process next VN
            confset[cur] += 1;  // skip next 0 to avoid all zero confset
        }
    }

    confsetCount++;
    // TODO
    if (confsetCount == nconf_q_1) {
        return 0;
    } else {
        return -1;
    }
}

bool NBCNode::isVN() {
    return false;
}
