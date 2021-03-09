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

class Node {
  protected:
    int degree;
    int inCount;
    std::vector<Node*> Nodes_;  // Linked Nodes
    std::vector<double> inValues_;
    std::vector<double> outValues_;

  public:
    Node(int d);
    ~Node();
    void Link(Node* n);
    virtual bool isVN() = 0;
    void setInValue(double data);
    bool isReady();
    virtual void Update() = 0;
};

class VNode : public Node {
  private:
    double value;

  public:
    VNode(int d, double v);
    void Link(Node* n);
    void Update();
    double getValue();
    bool isVN();
};

class CNode : public Node {
  public:
    CNode(int d);
    void Link(Node* n);
    void Update();
    bool isVN();
};

#endif