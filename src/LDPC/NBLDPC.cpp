#include "LDPC/NBLDPC.hpp"

#include <iostream>

#include "LDPC/Tanner.hpp"
#include "MatrixMath/MatrixMath.hpp"

NBLDPC::NBLDPC() {
    H_mat = Eigen::SparseMatrix<int>();
    G_mat = Eigen::SparseMatrix<int>();
    K = 0;
    N = 0;
    GF = 0;
}

NBLDPC::NBLDPC(Eigen::SparseMatrix<int> H) {
    // G.nRow = N - M, G.nCol = N
    // N = n, M = n - k
    // K = A.getnCol() - A.getnRow();
    H_mat = H;
    // G_mat = Eigen::SparseMatrix<int>(H_mat.cols() - H_mat.rows(),
    // H_mat.cols()); G_mat.setZero();
    G_mat =
        Eigen::SparseMatrix<int>(NBtransform_H_to_G(H_mat, GF).sparseView());
    assert(
        !(NBproduct(G_mat, H_mat.transpose(), GF).any()));  // G*H should be 0
    K = G_mat.rows();
    N = G_mat.cols();
    Eigen::MatrixXi diff =
        G_mat.block(0, N - K, K, K).unaryExpr([](const int x) {
            return x == 0 ? 0 : 1;
        }) -
        Eigen::MatrixXi::Identity(K, K);
    if (diff.any()) {
        isSystematic = false;
    } else {
        isSystematic = true;
    }
}

NBLDPC::NBLDPC(Alist<nbalist_matrix> A) {
    GF = A.getGF();
    memset(num_mlist, 0, A.getnRow() * sizeof(int));
    memcpy(num_mlist, A.getData().num_mlist, A.getnRow() * sizeof(int));
    memset(num_nlist, 0, A.getnCol() * sizeof(int));
    memcpy(num_nlist, A.getData().num_nlist, A.getnCol() * sizeof(int));
    new (this) NBLDPC(A.getMat());  // no need to delete, cause memory is reused
}

NBLDPC::NBLDPC(const char* filename) {
    Alist<nbalist_matrix> A = Alist<nbalist_matrix>(filename);
    GF = A.getGF();
    num_mlist = (int*)malloc(A.getnRow() * sizeof(int));
    memcpy(num_mlist, A.getData().num_mlist, A.getnRow() * sizeof(int));
    num_nlist = (int*)malloc(A.getnCol() * sizeof(int));
    memcpy(num_nlist, A.getData().num_nlist, A.getnCol() * sizeof(int));
    new (this) NBLDPC(A.getMat());  // no need to delete, cause memory is reused
}

NBLDPC::~NBLDPC() {
    free(num_mlist);
    free(num_nlist);
}

Eigen::RowVectorXi NBLDPC::encode(Eigen::RowVectorXi& m) const {
    return NBproduct(m, G_mat, GF);
}

/**
 * @brief decode the received LLR
 *
 * @param LLR received LLR, <GF, N>
 * @param iter_max stop early criterion
 * @param factor normalize factor, should not larger than 1
 * @param snr Channel snr
 * @param mode decode mode, 0: NMS, 1: SPA
 * @param n_max
 * @return Eigen::RowVectorXi
 */
Eigen::RowVectorXi NBLDPC::decode(Eigen::MatrixXd& LLR, int iter_max,
                                  double factor, double snr, int mode,
                                  int n_max) const {
    std::vector<NBVNode*> VNodes_;  // size: N
    std::vector<NBCNode*> CNodes_;  // size: M
    int M = H_mat.rows();
    Eigen::MatrixXi Hdense = H_mat.toDense();
    assert(LLR.cols() == N && LLR.rows() == GF);

    // init Nodes
    for (int i = 0; i < M; i++) {
        NBCNode* c = new NBCNode(num_mlist[i], factor, GF, n_max);
        CNodes_.push_back(c);
    }

    // TODO: use feature of SparseMatrix
    for (int i = 0; i < N; i++) {
        // init LLR directly by r
        NBVNode* v = new NBVNode(num_nlist[i], LLR.col(i), GF, n_max);
        VNodes_.push_back(v);
        for (int j = 0; j < M; j++) {
            if (Hdense(j, i)) {
                VNodes_[i]->Link(CNodes_[j], Hdense(j, i));
            }
        }
        assert(VNodes_[i]->isReady());
    }

    Eigen::RowVectorXi ret(N);
    int count = 0;
    do {
        for (NBVNode* v : VNodes_) {
            v->Update(3);
        }
        for (NBCNode* c : CNodes_) {
            c->Update(mode);
        }

        // update ret
        for (int i = 0; i < N; i++) {
            // printf("%d: ", i);
            ret[i] = VNodes_[i]->getValue();
            // printf("\n");
        }
        std::cout << ret << std::endl;
        count++;
    } while (NBproduct(H_mat.toDense(), ret.transpose(), GF).any() &&
             count < iter_max);  // stop criterion
    // std::cout << "迭代次数： " << count << std::endl;

    for (NBCNode* c : CNodes_) {
        delete c;
    }

    for (NBVNode* d : VNodes_) {
        delete d;
    }

    return ret;
}

Eigen::SparseMatrix<int> NBLDPC::getG() const {
    return G_mat;
}

Eigen::SparseMatrix<int> NBLDPC::getH() const {
    return H_mat;
}

int NBLDPC::getK() const {
    return K;
}

int NBLDPC::getN() const {
    return N;
}

int NBLDPC::getGF() const {
    return GF;
}

bool NBLDPC::getIsSys() const {
    return isSystematic;
}