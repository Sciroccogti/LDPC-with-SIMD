#include "LDPC/LDPC.hpp"

#include <iostream>

LDPC::LDPC() {
    H_mat = Eigen::SparseMatrix<int>();
    G_mat = Eigen::SparseMatrix<int>();
    K = 0;
}

LDPC::LDPC(Eigen::SparseMatrix<int> H) {
    // G.nRow = N - M, G.nCol = N
    // N = n, M = n - k
    // K = A.getnCol() - A.getnRow();
    H_mat = H;
    G_mat = Eigen::SparseMatrix<int>(transform_H_to_G(H_mat).sparseView());
    K = G_mat.rows();
}

LDPC::LDPC(Alist<alist_matrix> A) {
    memset(num_mlist, 0, A.getnRow() * sizeof(int));
    memcpy(num_mlist, A.getData().num_mlist, A.getnRow() * sizeof(int));
    memset(num_nlist, 0, A.getnCol() * sizeof(int));
    memcpy(num_nlist, A.getData().num_nlist, A.getnCol() * sizeof(int));
    new (this) LDPC(A.getMat());
}

LDPC::LDPC(const char* filename) {
    Alist<alist_matrix> A = Alist<alist_matrix>(filename);
    num_mlist = (int*)malloc(A.getnRow() * sizeof(int));
    memcpy(num_mlist, A.getData().num_mlist, A.getnRow() * sizeof(int));
    num_nlist = (int*)malloc(A.getnCol() * sizeof(int));
    memcpy(num_nlist, A.getData().num_nlist, A.getnCol() * sizeof(int));
    new (this) LDPC(A.getMat());
}

LDPC::~LDPC() {
    free(num_mlist);
    free(num_nlist);
}

Eigen::RowVectorXi LDPC::encode(Eigen::RowVectorXi& m) {
    return binaryproduct(m, G_mat.toDense());
}

Eigen::RowVectorXi LDPC::decode(Eigen::RowVectorXd& r) {
    std::vector<VNode*> VNodes_;
    std::vector<CNode*> CNodes_;
    int N = H_mat.cols();
    int M = H_mat.rows();
    assert(r.size() == N);

    // init Nodes
    for (int i = 0; i < M; i++) {
        CNode* c = new CNode(num_mlist[i]);
        CNodes_.push_back(c);
    }

    // TODO: use feature of SparseMatrix
    for (int i = 0; i < N; i++) {
        // init LLR directly by r
        VNode* v = new VNode(num_nlist[i], r[i]);
        VNodes_.push_back(v);
        for (int j = 0; j < M; j++) {
            if (H_mat.coeffRef(j, i) == 1) {
                VNodes_[i]->Link(CNodes_[j]);
            }
        }
        assert(VNodes_[i]->isReady());
    }

    // bool isDone = false;
    Eigen::RowVectorXi ret(N);
    do {
        for (VNode* v : VNodes_) {
            v->Update();
        }
        for (CNode* c : CNodes_) {
            c->Update();
        }

        // update ret
        for (int i = 0; i < N; i++) {
            if (VNodes_[i]->getValue() >= 0) {
                ret[i] = 1;
            } else {
                ret[i] = -1;
            }
        }
    } while (binaryproduct(ret, H_mat.toDense().transpose()).any());
     
    return ret;
}

Eigen::SparseMatrix<int> LDPC::getG() {
    return G_mat;
}

Eigen::SparseMatrix<int> LDPC::getH() {
    return H_mat;
}

int LDPC::getK() {
    return K;
}
