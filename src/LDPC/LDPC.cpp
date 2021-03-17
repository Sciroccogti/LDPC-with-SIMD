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
    N = G_mat.cols();
    Eigen::MatrixXi diff =
        G_mat.block(0, N - K, K, K) - Eigen::MatrixXi::Identity(K, K);
    if (diff.any()) {
        isSystematic = false;
    } else {
        isSystematic = true;
    }
}

LDPC::LDPC(Alist<alist_matrix> A) {
    memset(num_mlist, 0, A.getnRow() * sizeof(int));
    memcpy(num_mlist, A.getData().num_mlist, A.getnRow() * sizeof(int));
    memset(num_nlist, 0, A.getnCol() * sizeof(int));
    memcpy(num_nlist, A.getData().num_nlist, A.getnCol() * sizeof(int));
    new (this) LDPC(A.getMat());  // no need to delete, cause memory is reused
}

LDPC::LDPC(const char* filename) {
    Alist<alist_matrix> A = Alist<alist_matrix>(filename);
    num_mlist = (int*)malloc(A.getnRow() * sizeof(int));
    memcpy(num_mlist, A.getData().num_mlist, A.getnRow() * sizeof(int));
    num_nlist = (int*)malloc(A.getnCol() * sizeof(int));
    memcpy(num_nlist, A.getData().num_nlist, A.getnCol() * sizeof(int));
    new (this) LDPC(A.getMat());  // no need to delete, cause memory is reused
}

LDPC::~LDPC() {
    free(num_mlist);
    free(num_nlist);
}

Eigen::RowVectorXi LDPC::encode(Eigen::RowVectorXi& m) const {
    return binaryproduct(m, G_mat.toDense());
}

/**
 * @brief decode the received sequence
 *
 * @param r received sequence
 * @param iter_max stop early criterion
 * @param factor normalize factor, should not larger than 1
 * @return Eigen::RowVectorXi
 */
Eigen::RowVectorXi LDPC::decode(Eigen::RowVectorXd& r, int iter_max,
                                double factor) const {
    std::vector<VNode*> VNodes_;
    std::vector<CNode*> CNodes_;
    int M = H_mat.rows();
    Eigen::MatrixXi Hdense = H_mat.toDense();
    assert(r.size() == N);

    // init Nodes
    for (int i = 0; i < M; i++) {
        CNode* c = new CNode(num_mlist[i], factor);
        CNodes_.push_back(c);
    }

    // TODO: use feature of SparseMatrix
    for (int i = 0; i < N; i++) {
        // init LLR directly by r
        VNode* v = new VNode(num_nlist[i], r[i]);
        VNodes_.push_back(v);
        for (int j = 0; j < M; j++) {
            if (Hdense(j, i) == 1) {
                VNodes_[i]->Link(CNodes_[j]);
            }
        }
        assert(VNodes_[i]->isReady());
    }

    Eigen::RowVectorXi ret(N);
    int count = 0;
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
                ret[i] = 0;
            } else {
                ret[i] = 1;
            }
        }
        count++;
    } while (binaryproduct(ret, H_mat.toDense().transpose()).any() &&
             count <= iter_max);  // stop criterion
    // std::cout << "迭代次数： " << count << std::endl;

    for (CNode* c : CNodes_) {
        delete c;
    }

    for (VNode* d : VNodes_) {
        delete d;
    }

    return ret;
}

/**
 * @brief recover message from decoded sequence
 *     https://github.com/hichamjanati/pyldpc/blob/master/pyldpc/decoder.py#L186
 * @param d decoded sequence, should be 0,1 sequence!
 * @return Eigen::RowVectorXi
 */
Eigen::RowVectorXi LDPC::recoverMessage(Eigen::RowVectorXi& d) const {
    Eigen::RowVectorXi ret(K);

    if (isSystematic) {
        ret = d.rightCols(K);
    } else {
        Eigen::MatrixXi tG = G_mat.toDense().transpose();
        Eigen::RowVectorXi d_ = d.unaryExpr([](const int x) {
            assert(x == 1 || x == 0);
            return 2 * x - 1;
        });
        gausselimination(tG, d_);
        ret.setZero();
        ret[K - 1] = d[K - 1];
        for (int i = K - 2; i >= 0; i--) {
            ret[i] = d[i];
            Eigen::RowVectorXi tmpGrow = tG.block(i, i + 1, 1, K - i - 1);
            ret[i] -= tmpGrow.dot(ret.block(0, i + 1, 1, K - i - 1)) % 2;
        }
        // mod 2 and transbpsk
        ret = ret.unaryExpr([](const int x) { return (x % 2 + 1) / 2; });
    }
    return ret;
}

Eigen::SparseMatrix<int> LDPC::getG() const {
    return G_mat;
}

Eigen::SparseMatrix<int> LDPC::getH() const {
    return H_mat;
}

int LDPC::getK() const {
    return K;
}

int LDPC::getN() const {
    return N;
}

bool LDPC::getIsSys() const {
    return isSystematic;
}
