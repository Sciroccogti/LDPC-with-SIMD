/*
 * File: Alist.hpp
 * File Created: Thursday, 29th October 2020 15:11:30
 * Author: Yifan Zhang (scirocco_gti@yeah.net)
 * Last Modified: Thursday, 29th October 2020 19:43:47
 */

#ifndef ALIST_H
#define ALIST_H

#include "alist/alist_matrix.h"

template <class T>
class Alist {
  private:
    T data;
    bool isInited;

  public:
    Alist();
    Alist(T d);
    ~Alist();
    Alist(const char* filename);
    int load(const char* filename);
    void save(const char* filename);
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
    read_alist(fp, &data);
    isInited = true;
    fclose(fp);
}

template <class T>
void Alist<T>::save(const char* filename) {
    FILE* fp = fopen(filename, "w+");
    write_alist(fp, &data);
    fclose(fp);
}

#endif