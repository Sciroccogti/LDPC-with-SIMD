/**
 * @file Tanner.hpp
 * @author Sciroccogti (scirocco_gti@yeah.net)
 * @brief
 * @date 2021-03-09 17:32:42
 * @modified: 2021-04-22 11:31:12
 */

#ifndef TANNER_HPP
#define TANNER_HPP

#include <cassert>
#include <cmath>
#include <vector>

#include "LDPC/NBLLR.hpp"
#include "MatrixMath/MatrixMath.hpp"

#define BP_NMS 0
#define BP_SPA 1
#define BP_EMS 2
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
    int n_max;
    int inCount;  // used to record No. of current inputting node
    int n_maxCount;
    std::vector<NBNode*> NBNodes_;  // Linked Nodes
    std::vector<Eigen::RowVectorXd> inValuesQ_;
    // every element's size is n_max, number of elements is `degree`
    std::vector<std::vector<NBLLR>> n_maxValue_;

  public:
    NBNode(int d, int gf, int nmax);
    ~NBNode();
    virtual void Link(NBNode* n, int h) = 0;
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
    NBVNode(int d, Eigen::RowVectorXd vQ, const int GF, const int nmax);
    void Link(NBNode* n, int h);
    void Update(int n_max);
    int getValue();
    bool isVN();
};

class NBCNode : public NBNode {
  private:
    double factor;
    Eigen::RowVectorXi Hrow_;

  public:
    NBCNode(int d, double f, const int gf, const int nmax);
    void Link(NBNode* n, int h);
    void Update(int mode);
    int getConfset(Eigen::RowVectorXi& confset, int& confsetCount, int& cur);
    bool isVN();
};

#endif