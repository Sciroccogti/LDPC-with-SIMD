/**
 * @file Alist.hpp
 * @author Yifan Zhang (scirocco_gti@yeah.net)
 * @brief
 * @date 2020-10-29 15:11:30
 * @modified: 2021-03-26 21:49:42
 */

#ifndef ALIST_HPP
#define ALIST_HPP

#include <Eigen/Eigen>

#include "Alist/alist_matrix.h"
#include "Alist/nbalist_matrix.h"

// a class for .alist IO
template <class T>
class Alist {
  private:
    int GF;
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
    int getnRow();  // M
    int getnCol();  // N
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
void Alist<T>::save(const char* filename) {
    FILE* fp = fopen(filename, "w+");
    write_alist(fp, &data);
    fclose(fp);
}

template <class T>
T Alist<T>::getData() {
    return data;
}

// M
template <class T>
int Alist<T>::getnRow() {
    return data.M;
}

// N
template <class T>
int Alist<T>::getnCol() {
    return data.N;
}

#endif