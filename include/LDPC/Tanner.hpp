/**
 * @file Tanner.hpp
 * @author Sciroccogti (scirocco_gti@yeah.net)
 * @brief
 * @date 2021-03-09 17:32:42
 * @modified: 2021-04-20 01:55:48
 */

#ifndef TANNER_HPP
#define TANNER_HPP

#include <cassert>
#include <cmath>
#include <vector>

#include "MatrixMath/MatrixMath.hpp"
#include "LDPC/NBLLR.hpp"

#define BP_NMS 0
#define BP_SPA 1
#define BP_QNMS 2
#define BP_QSPA 3

extern const char* Modes_[4];

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
    double LLR;

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

class NBNode {
  protected:
    int GF;
    int degree;
    int inCount;                // used to record No. of current inputting node
    int n_maxCount;
    std::vector<NBNode*> NBNodes_;  // Linked Nodes
    std::vector<Eigen::RowVectorXd> inValuesQ_;
    std::vector<std::vector<NBLLR>> n_maxValue_;

  public:
    NBNode(int d, int gf);
    ~NBNode();
    void Link(NBNode* n);
    virtual bool isVN() = 0;
    void setInValue(Eigen::RowVectorXd dataQ);
    void setinn_maxValue(std::vector<NBLLR> vn_max);
    bool isReady();
    virtual void Update(int mode) = 0;
};

class NBVNode : public NBNode {
  private:
    Eigen::RowVectorXd valueQ;
    Eigen::RowVectorXd LLRQ;

  public:
    NBVNode(int d, Eigen::RowVectorXd vQ, const int GF);
    void Link(NBNode* n);
    void Update(int mode);
    int getValue();
    bool isVN();
};

class NBCNode : public NBNode {
  private:
    double factor;

  public:
    NBCNode(int d, double f, const int gf);
    void Link(NBNode* n);
    void Update(int mode);
    bool isVN();
};

#endif