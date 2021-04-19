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
NBNode::NBNode(int d, int gf) {
    degree = d;
    GF = gf;
    // Nodes_.reserve(d);

    inCount = 0;
    // init inValues
    for (int i = 0; i < degree; i++) {
        inValuesQ_.push_back(Eigen::RowVectorXd::Zero(GF - 1));
    }

    n_maxCount = 0;
    // init inn_maxValue
    for (int i = 0; i < degree; i++) {
        n_maxValue_.push_back(vector<NBLLR>());
    }
}

NBNode::~NBNode() {}

// Link node n to this node
void NBNode::Link(NBNode* n) {
    assert(NBNodes_.size() + 1 <= degree);
    NBNodes_.push_back(n);
}

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
NBVNode::NBVNode(int d, Eigen::RowVectorXd vQ, const int gf) : NBNode(d, gf) {
    valueQ = vQ;
    LLRQ = vQ;
    assert(valueQ.size() == GF - 1);
}

void NBVNode::Link(NBNode* n) {
    assert(!n->isVN());  // assure n is CN
    NBNode::Link(n);     // register n to self
    n->Link(this);       // register self to n
}

void NBVNode::Update(int n_max) {
    valueQ = LLRQ;
    // update own value
    for (int i = 0; i < degree; i++) {
        valueQ = inValuesQ_[i] + valueQ;
    }

    // update outputs
    for (int i = 0; i < degree; i++) {
        NBNodes_[i]->setInValue(valueQ - inValuesQ_[i]);

        // https://docs.microsoft.com/zh-cn/troubleshoot/cpp/stl-priority-queue-class-custom-type
        priority_queue<NBLLR, vector<NBLLR>, less<vector<NBLLR>::value_type>>
            pqLLR;  // size should be GF - 1
        for (int q = 0; q < GF - 1; q++) {
            pqLLR.push(NBLLR(valueQ[q], q + 1));
        }

        std::vector<NBLLR> n_maxLLR_;  // size should be n_max
        for (int j = 0; j < n_max; j++) {
            n_maxLLR_.push_back(pqLLR.top());
            pqLLR.pop();
        }

        NBNodes_[i]->setinn_maxValue(n_maxLLR_);
        // priority_queue<NBLLR, vector<NBLLR>, less<vector<NBLLR>::value_type>>
        //     pqcopy(pqLLR);
        printf("pq:\n");
        for (size_t j = 0; j < n_maxLLR_.size(); j++) {
            printf("%.2lf:%d   ", n_maxLLR_[j].getData(), n_maxLLR_[j].getQ());
        }
        printf("\n");
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
    for (int i = 1; i < GF - 1; i++) {
        // printf("%.2lf ", valueQ[i]);
        if (max < valueQ[i]) {
            max = valueQ[i];
            ret = i;
        }
    }
    // printf(" -> %.0lf ", max);
    if (max < 0) {  // output 0 as others are larger than 0's probability
        return 0;
    }
    return ret + 1;
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
NBCNode::NBCNode(int d, double f, const int gf) : NBNode(d, gf) {
    assert(f <= 1);
    factor = f;
}

void NBCNode::Link(NBNode* n) {
    assert(n->isVN());  // assure n is VN
    NBNode::Link(n);    // register n to self
}

void NBCNode::Update(int mode) {
    assert(mode >= BP_QNMS && mode <= BP_QSPA);
    assert(inCount == 0);  // assure all inValue is changed

    for (int i = 0; i < degree; i++) {
        switch (mode) {
            case BP_QNMS: {
                Eigen::RowVectorXd min = Eigen::RowVectorXd::Zero(GF - 1);
                Eigen::RowVectorXi sgn = Eigen::RowVectorXi::Ones(GF - 1);
                bool isFirst = true;
                for (int j = 0; j < degree; j++) {
                    if (j == i) {
                        continue;  // skip current VN
                    }
                    Eigen::RowVectorXd absJ(GF - 1);
                    for (int q = 0; q < GF - 1; q++) {
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

                for (int q = 0; q < GF - 1; q++) {
                    min[q] = sgn[q] * min[q] * factor;
                }

                NBNodes_[i]->setInValue(min);
            } break;
            case BP_QSPA: {
                Eigen::RowVectorXd prod = Eigen::RowVectorXd::Ones(GF - 1);
                Eigen::RowVectorXi sgn = Eigen::RowVectorXi::Ones(GF - 1);
                for (int j = 0; j < degree; j++) {
                    if (j == i) {
                        continue;  // skip current VN
                    }
                    for (int q = 0; q < GF - 1; q++) {
                        prod[q] *= tanh(fabs(inValuesQ_[j][q]) / 2);
                        sgn[q] *= inValuesQ_[j][q] >= 0 ? 1 : -1;
                    }
                }

                for (int q = 0; q < GF - 1; q++) {
                    prod[q] = prod[q] < 1
                                  ? prod[q]
                                  : 1.0 - numeric_limits<double>::epsilon();
                    prod[q] = sgn[q] * 2 * fabs(atanh(prod[q]));
                }

                NBNodes_[i]->setInValue(prod);
            } break;
            default:
                break;
        }
    }
}

bool NBCNode::isVN() {
    return false;
}
