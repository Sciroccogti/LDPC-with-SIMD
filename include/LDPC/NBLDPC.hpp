/**
 * @file NBLDPC.hpp
 * @author Yifan Zhang (scirocco_gti@yeah.net)
 * @brief
 * @date 2021-03-27 10:36:37
 * @modified: 2021-03-29 15:41:54
 */

#ifndef NBLDPC_HPP
#define NBLDPC_HPP

#include "Alist/Alist.hpp"

class NBLDPC {
  private:
    Eigen::SparseMatrix<int> H_mat;  // rows=N-K, cols=N
    Eigen::SparseMatrix<int> G_mat;  // rows=K, cols=N
    int* num_mlist;
    int* num_nlist;
    int K;  // length of message
    int N;
    int GF;  // GF order
    bool isSystematic;

  public:
    NBLDPC();
    NBLDPC(Eigen::SparseMatrix<int> H);
    NBLDPC(Alist<nbalist_matrix>);
    NBLDPC(const char* filename);
    ~NBLDPC();

    Eigen::RowVectorXi encode(Eigen::RowVectorXi& m) const;
    Eigen::SparseMatrix<int> getG() const;
    Eigen::SparseMatrix<int> getH() const;
    int getK() const;
    int getN() const;
    bool getIsSys() const;
};

#endif
