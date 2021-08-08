#include "LDPC/Tanner.hpp"

const char* Modes_[4] = {"NMS", "SPA", "EMS", "QSPA"};

/**
 * @brief Construct a new Node:: Node object
 *  inValues are initialized as 0
 * @param d degree
 */
Node::Node(int d) {
    degree = d;
    inCount = 0;
    // Nodes_.reserve(d);

    // init inValues
    for (int i = 0; i < degree; i++) {
        inValues_.push_back(0);
    }
}

Node::~Node() {}

// Link node n to this node
void Node::Link(Node* n) {
    assert(Nodes_.size() + 1 <= degree);
    Nodes_.push_back(n);
}

// return true if size == degree
bool Node::isReady() {
    return Nodes_.size() == degree;
}

void Node::setInValue(float data) {
    inValues_[inCount] = data;
    inCount++;
    if (inCount >= degree) {
        inCount = 0;
    }
}

VNode::VNode(int d, float v) : Node(d) {
    value = v;
    LLR = v;
}

void VNode::Link(Node* n) {
    assert(!n->isVN());  // assure n is CN
    Node::Link(n);       // register n to self
    n->Link(this);       // register self to n
}

void VNode::Update(int mode) {
    // update own value
    // for (float i : inValues_) {
    //     value += i;
    // }
    value = LLR;
    for (int i = 0; i < degree; i++) {
        value += inValues_[i];
    }

    // update outputs
    for (int i = 0; i < degree; i++) {
        Nodes_[i]->setInValue(value - inValues_[i]);
    }
}

float VNode::getValue() {
    float ret = LLR;
    for (int i = 0; i < degree; i++) {
        ret += inValues_[i];
    }
    return ret;
}

bool VNode::isVN() {
    return true;
}

/**
 * @brief Construct a new CNode::CNode object
 *
 * @param d degree
 * @param f normalize factor
 */
CNode::CNode(int d, float f) : Node(d) {
    assert(f <= 1);
    factor = f;
}

void CNode::Link(Node* n) {
    assert(n->isVN());  // assure n is VN
    Node::Link(n);      // register n to self
}

void CNode::Update(int mode) {
    assert(mode >= 0 && mode <= BP_SPA);
    assert(inCount == 0);  // assure all inValue is changed

    switch (mode) {
        case BP_NMS: {
            float min = MAXFLOAT, min2 = MAXFLOAT;  // min2 store the second min
            int sgn = 1;
            for (int i = 0; i < degree; i++) {
                float absJ = fabs(inValues_[i]);
                if (absJ < min) {
                    min2 = min;
                    min = absJ;
                } else if (absJ < min2) {
                    min2 = absJ;
                }
                sgn *= inValues_[i] >= 0 ? 1 : -1;
            }

            for (int i = 0; i < degree; i++) {
                int sgnTmp = sgn / (inValues_[i] >= 0 ? 1 : -1);
                if (fabs(inValues_[i]) - min <= 0) {
                    Nodes_[i]->setInValue(sgnTmp * min2 * factor);
                } else {
                    Nodes_[i]->setInValue(sgnTmp * min * factor);
                }
            }
        } break;
        case BP_SPA: {
            float prod = 1;
            int sgn = 1;
            for (int i = 0; i < degree; i++) {
                prod *= tanh(fabs(inValues_[i]) / 2);
                sgn *= inValues_[i] >= 0 ? 1 : -1;
            }

            prod =
                prod < 1 ? prod : 1.0 - std::numeric_limits<float>::epsilon();
            for (int i = 0; i < degree; i++) {
                int sgnTmp = sgn / (inValues_[i] >= 0 ? 1 : -1);
                Nodes_[i]->setInValue(
                    sgnTmp * 2 * atanh(prod / tanh(fabs(inValues_[i]) / 2)));
            }
        } break;
        default:
            break;
    }
}

bool CNode::isVN() {
    return false;
}
