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
#include "MatrixMath/MatrixMath.hpp"
#include "LDPC/Tanner.hpp"

class LDPC {
  private:
    Eigen::SparseMatrix<int> H_mat;
    Eigen::SparseMatrix<int> G_mat;
    int* num_mlist;
    int* num_nlist;
    int K;  // length of message

  public:
    LDPC();
    LDPC(Eigen::SparseMatrix<int> H);
    LDPC(Alist<alist_matrix>);
    LDPC(const char* filename);
    // TODO: add H check
    ~LDPC();

    Eigen::RowVectorXi encode(Eigen::RowVectorXi& m);
    Eigen::RowVectorXi decode(Eigen::RowVectorXi& r);
    Eigen::SparseMatrix<int> getG();
    Eigen::SparseMatrix<int> getH();
    int getK();
};

#endif