/**
 * @file Tanner.hpp
 * @author Sciroccogti (scirocco_gti@yeah.net)
 * @brief
 * @date 2021-03-09 17:32:42
 * @modified: 2021-04-27 13:13:37
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
    std::vector<float> inValues_;

  public:
    Node(int d);
    ~Node();
    void Link(Node* n);
    virtual bool isVN() = 0;
    void setInValue(float data);
    bool isReady();
    virtual void Update(int mode) = 0;
};

class VNode : public Node {
  private:
    float value;
    float LLR;

  public:
    VNode(int d, float v);
    void Link(Node* n);
    void Update(int mode);
    float getValue();
    bool isVN();
};

class CNode : public Node {
  private:
    float factor;

  public:
    CNode(int d, float f);
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
    std::vector<std::vector<float>> inValuesQ_;
    // every element's size is n_max, number of elements is `degree`
    std::vector<std::vector<NBLLR>> n_maxValue_;

  public:
    NBNode(int d, int gf, int nmax);
    ~NBNode();
    virtual void Link(NBNode* n, int h) = 0;
    virtual bool isVN() = 0;
    void setInValue(std::vector<float> dataQ);
    void setinn_maxValue(std::vector<NBLLR> vn_max);
    bool isReady();
    virtual void Update(int mode) = 0;
};

class NBVNode : public NBNode {
  private:
    std::vector<float> valueQ;
    std::vector<float> LLRQ;

  public:
    NBVNode(int d, Eigen::RowVectorXf vQ, const int GF, const int nmax);
    void Link(NBNode* n, int h);
    void Update(int n_max);
    int getValue();
    bool isVN();
};

class NBCNode : public NBNode {
  private:
    float factor;
    std::vector<int> Hrow_;

  public:
    NBCNode(int d, float f, const int gf, const int nmax);
    void Link(NBNode* n, int h);
    void Update(int mode);
    int getConfset(std::vector<int>& confset, int& confsetCount, int& cur);
    bool isVN();
};

#endif