/*
 * File: Alist.hpp
 * File Created: Thursday, 29th October 2020 15:11:30
 * Author: Yifan Zhang (scirocco_gti@yeah.net)
 * Last Modified: Friday, 30th October 2020 18:02:46
 */

#ifndef ALIST_HPP
#define ALIST_HPP

#include <Eigen/Eigen>

#include "alist/alist_matrix.h"

// a class for .alist IO
template <class T>
class Alist {
  private:
    T data;
    bool isInited;

  public:
    Alist();
    Alist(T d);
    Alist(const char* filename);
    ~Alist();
    int load(const char* filename);
    void save(const char* filename);
    T getData();
    Eigen::SparseMatrix<int> getMat();
};

template <class T>
Alist<T>::Alist() {
    isInited = false;
}

template <class T>
Alist<T>::Alist(T d) {
    data = d;
}

template <class T>
Alist<T>::~Alist() {
    if (isInited) {
        int i;
        for (i = 0; i < data.M; i++) {
            free(data.mlist[i]);
        }
        for (i = 0; i < data.N; i++) {
            free(data.nlist[i]);
        }

        free(data.mlist);
        free(data.nlist);
        free(data.num_mlist);
        free(data.num_nlist);
    }
}

template <class T>
Alist<T>::Alist(const char* filename) {
    FILE* fp = fopen(filename, "r");
    read_alist(fp, &data);
    isInited = true;
    fclose(fp);
}

template <class T>
int Alist<T>::load(const char* filename) {
    FILE* fp = fopen(filename, "r");
    int ret = read_alist(fp, &data);
    isInited = true;
    fclose(fp);
    return ret;
}

template <class T>
void Alist<T>::save(const char* filename) {
    FILE* fp = fopen(filename, "w+");
    write_alist(fp, &data);
    fclose(fp);
}

template <class T>
T Alist<T>::getData() {
    return data;
}

template <class T>
Eigen::SparseMatrix<int> Alist<T>::getMat() {
    Eigen::SparseMatrix<int> ret(data.N, data.M);
    ret.reserve(Eigen::VectorXi::Constant(data.M, data.biggest_num_m));
    int i, j;
    for (i = 0; i < data.M; i++) {
        for (j = 0; j < data.num_mlist[i]; j++) {
            ret.insert(data.mlist[i][j] - 1, i) = 1;
        }
    }
    ret.makeCompressed();
    return ret;
}

#endif