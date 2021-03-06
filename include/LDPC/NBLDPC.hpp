/**
 * @file NBLDPC.hpp
 * @author Yifan Zhang (scirocco_gti@yeah.net)
 * @brief
 * @date 2021-03-27 10:36:37
 * @modified: 2021-04-27 13:09:49
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
    Eigen::RowVectorXi decode(Eigen::MatrixXf& LLR, int iter_max, float factor,
                              float snr, int mode, int n_max) const;
    Eigen::RowVectorXi recoverMessage(Eigen::RowVectorXi& d) const;
    Eigen::SparseMatrix<int> getG() const;
    Eigen::SparseMatrix<int> getH() const;
    int getK() const;
    int getN() const;
    int getGF() const;
    bool getIsSys() const;
};

#endif
