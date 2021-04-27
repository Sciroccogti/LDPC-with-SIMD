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
    return value;
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

    for (int i = 0; i < degree; i++) {
        switch (mode) {
            case BP_NMS: {
                float min = 0;
                int sgn = 1;
                bool isFirst = true;
                for (int j = 0; j < degree; j++) {
                    float absJ = fabs(inValues_[j]);
                    if (j == i) {
                        continue;  // skip current VN
                    }
                    if (isFirst) {
                        isFirst = false;
                        min = absJ;
                    } else {
                        min = absJ < min ? absJ : min;
                    }
                    sgn *= inValues_[j] >= 0 ? 1 : -1;
                }

                Nodes_[i]->setInValue(sgn * min * factor);
            } break;
            case BP_SPA: {
                float prod = 1;
                int sgn = 1;
                for (int j = 0; j < degree; j++) {
                    if (j == i) {
                        continue;  // skip current VN
                    }
                    prod *= tanh(fabs(inValues_[j]) / 2);
                    sgn *= inValues_[j] >= 0 ? 1 : -1;
                }

                prod = prod < 1 ? prod
                                : 1.0 - std::numeric_limits<float>::epsilon();
                Nodes_[i]->setInValue(sgn * 2 * atanh(prod));
            } break;
            default:
                break;
        }
    }
}

bool CNode::isVN() {
    return false;
}
