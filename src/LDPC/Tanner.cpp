#include "LDPC/Tanner.hpp"

Node::Node(int d) {
    degree = d;
    inCount = 0;
    // Nodes_.reserve(d);

    // init inValues
    for (int i = 0; i < degree; i++) {
        inValues_.push_back(0);
    }
    // init outValues
    for (int i = 0; i < degree; i++) {
        outValues_.push_back(0);
    }
}

Node::~Node() {}

// Link node n to this node
void Node::Link(Node* n) {
    assert(Nodes_.size() + 1 <= degree);
    Nodes_.push_back(n);
}

bool Node::isReady() {
    return Nodes_.size() == degree;
}

void Node::setInValue(double data) {
    inValues_[inCount] = data;
    inCount++;
    if (inCount >= degree) {
        inCount = 0;
    }
}

VNode::VNode(int d, double v) : Node(d) {
    value = v;
}

void VNode::Link(Node* n) {
    assert(!n->isVN());  // assure n is CN
    Node::Link(n);       // register n to self
    n->Link(this);       // register self to n
}

void VNode::Update() {
    // update own value
    // for (double i : inValues_) {
    //     value += i;
    // }
    for (int i = 0; i < degree; i++) {
        value += inValues_[i];
    }

    // update outputs
    for (int i = 0; i < degree; i++) {
        Nodes_[i]->setInValue(value - inValues_[i]);
    }
}

double VNode::getValue() {
    return value;
}

bool VNode::isVN() {
    return true;
}

CNode::CNode(int d) : Node(d) {}

void CNode::Link(Node* n) {
    assert(n->isVN());  // assure n is VN
    Node::Link(n);      // register n to self
}

void CNode::Update() {
    assert(inCount == 0);  // assure all inValue is changed

    for (int i = 0; i < degree; i++) {
        double min = 0;
        int sgn = 1;
        bool isFirst = true;
        for (int j = 0; j < degree; j++) {
            double absJ = fabs(inValues_[j]);
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

        Nodes_[i]->setInValue(sgn * min);
    }
}

bool CNode::isVN() {
    return false;
}