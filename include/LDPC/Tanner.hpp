/**
 * @file Tanner.hpp
 * @author Sciroccogti (scirocco_gti@yeah.net)
 * @brief
 * @date 2021-03-09 17:32:42
 * @modified: 2021-03-09 17:32:51
 */

#ifndef TANNER_HPP
#define TANNER_HPP

#include <cassert>
#include <cmath>
#include <vector>

#define BP_NMS 0
#define BP_SPA 1

extern const char * Modes_[2];

class Node {
  protected:
    int degree;
    int inCount;                // used to record No. of current inputting node
    std::vector<Node*> Nodes_;  // Linked Nodes
    std::vector<double> inValues_;

  public:
    Node(int d);
    ~Node();
    void Link(Node* n);
    virtual bool isVN() = 0;
    void setInValue(double data);
    bool isReady();
    virtual void Update(int mode) = 0;
};

class VNode : public Node {
  private:
    double value;

  public:
    VNode(int d, double v);
    void Link(Node* n);
    void Update(int mode);
    double getValue();
    bool isVN();
};

class CNode : public Node {
  private:
    double factor;

  public:
    CNode(int d, double f);
    void Link(Node* n);
    void Update(int mode);
    bool isVN();
};

#endif