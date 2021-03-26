#include <Alist/Alist.hpp>

template <>
Alist<alist_matrix>::Alist(const char* filename) {
    FILE* fp = fopen(filename, "r");
    read_alist(fp, &data);
    isInited = true;
    GF = 2;

    fclose(fp);
}

template <>
Alist<nbalist_matrix>::Alist(const char* filename) {
    FILE* fp = fopen(filename, "r");
    read_nbalist(fp, &data);
    isInited = true;
    GF =data.GF;

    fclose(fp);
}

template <>
Eigen::SparseMatrix<int> Alist<alist_matrix>::getMat() {
    Eigen::SparseMatrix<int> ret(data.M, data.N);
    ret.reserve(Eigen::VectorXi::Constant(data.N, data.biggest_num_n));
    int i, j;
    for (i = 0; i < data.N; i++) {
        for (j = 0; j < data.num_nlist[i]; j++) {
            ret.insert(data.nlist[i][j] - 1, i) = 1;
        }
    }
    ret.makeCompressed();
    return ret;
}

template <>
Eigen::SparseMatrix<int> Alist<nbalist_matrix>::getMat() {
    Eigen::SparseMatrix<int> ret(data.M, data.N);
    ret.reserve(Eigen::VectorXi::Constant(data.N, data.biggest_num_n));
    int i, j;
    for (i = 0; i < data.N; i++) {
        for (j = 0; j < data.num_nlist[i]; j++) {
            ret.insert(data.nlist[i][j] - 1, i) = data.nGFlist[i][j];
        }
    }
    ret.makeCompressed();
    return ret;
}
