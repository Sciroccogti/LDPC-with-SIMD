/**
 * @file LDPC.hpp
 * @author Yifan Zhang (scirocco_gti@yeah.net)
 * @brief
 * @date 2020-10-30 22:52:43
 * @modified: 2020-10-31 12:04:23
 */

#ifndef LDPC_HPP
#define LDPC_HPP

#include "Alist/Alist.hpp"
#include "LDPC/Tanner.hpp"
#include "MatrixMath/MatrixMath.hpp"

class LDPC {
  private:
    Eigen::SparseMatrix<int> H_mat; // rows=N-K, cols=N
    Eigen::SparseMatrix<int> G_mat; // rows=K, cols=N
    int* num_mlist;
    int* num_nlist;
    int K;  // length of message
    int N;
    bool isSystematic;

  public:
    LDPC();
    LDPC(Eigen::SparseMatrix<int> H);
    LDPC(Alist<alist_matrix>);
    LDPC(const char* filename);
    ~LDPC();

    Eigen::RowVectorXi encode(Eigen::RowVectorXi& m) const;
    Eigen::RowVectorXi decode(Eigen::RowVectorXd& r, int iter_max,
                              double factor, int mode) const;
    Eigen::RowVectorXi recoverMessage(Eigen::RowVectorXi& d) const;
    Eigen::SparseMatrix<int> getG() const;
    Eigen::SparseMatrix<int> getH() const;
    int getK() const;
    int getN() const;
    bool getIsSys() const;
};

#endif