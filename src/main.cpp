#include "main.hpp"

#include <omp.h>

#include <bitset>
#include <iostream>
#include <vector>

#include "Alist/Alist.hpp"
#include "MatrixMath/NoramlMath.hpp"
#include "MatrixMath/SIMDMath.hpp"
// #include "LDPC/LDPC.hpp"

int main(int argc, char* argv[]) {
    // Eigen::initParallel();
    char* alist_path = NULL;
    bool enable_SIMD = true;

    if (opt(argc, argv, alist_path, enable_SIMD)) {
        return -1;
    }
#ifdef MIPP_ALIGNED_LOADS
    printf("aligned!\n");
#endif

    // Eigen::SparseMatrix<int> H = ldpc.getH(), G = ldpc.getG();
    // std::cout << G << H << std::endl;
    // std::cout << G * H.transpose() << std::endl;

    Alist<alist_matrix> G_alist(alist_path);
    const size_t M = 64, N = 128;
    Eigen::SparseMatrix<int> G = G_alist.getMat();

    int count[23] = {0};

    double start = omp_get_wtime();
    if (enable_SIMD) {
        // store a 64 bit with int64
        vector_type G_int;
        // size for which the vectorization is possible
        size_t vec_size = N - N % inc;

        // store a row as a uint64, and every 4 uint64 as a batch
        // outerSize: M
        // https://eigen.tuxfamily.org/dox/group__TutorialSparse.html
        for (int k = 0; k < G.outerSize(); k++) {
            // bitset to turn bits into uint64
            // https://www.cplusplus.com/reference/bitset/bitset/
            std::bitset<M> col(0);
            for (Eigen::SparseMatrix<int>::InnerIterator it(G, k); it; ++it) {
                col[it.row()] = it.value();
            }
            G_int.push_back(col.to_ullong());
        }

        vector_type two(inc, 2);
        b_type two_vec = &two[0];

        start = omp_get_wtime();
#pragma omp parallel for reduction(+ : count[:23])
        for (size_t i = 0; i < 4096000; i++) {
            int weight = 0;
            bool is_wasted = false;
            vector_type I(inc, i);

#pragma omp parallel for
            for (int j = 0; j < vec_size; j += inc) {
                if (is_wasted) {
                    continue;
                }
                b_type G_vec = &G_int[j];
                b_type I_vec = &I[0];
                b_type tmp_vec = G_vec & I_vec;
                tmp_vec = hamming(tmp_vec);
                b_type tmp2_vec = (tmp_vec & b0);
                for (size_t k = 0; k < inc; k++) {
                    if (tmp_vec[k] != tmp2_vec[k]) {
                        if (weight >= 22) {
                            is_wasted = true;
                            break;
                        }
                        weight++;
                    }
                }
            }
            if (is_wasted) {
                continue;
            }
            count[weight]++;
        }

    } else {  // Non SIMD

        // store a 64 bit with int64
        u_int64_t G_int[N];
#pragma omp parallel for
        for (size_t i = 0; i < N; i++) {
            std::bitset<M> col(0);
            for (Eigen::SparseMatrix<int>::InnerIterator it(G, i); it; ++it) {
                col[it.row()] = it.value();
            }
            G_int[i] = col.to_ullong();
        }
#pragma omp parallel for reduction(+ : count[:23])
        for (size_t i = 0; i < 4096000; i++) {
            int weight = 0;
            bool is_wasted = false;
            for (size_t j = 0; j < N; j++) {
                u_int64_t G_tmp = G_int[j] & i;
                G_tmp = hamming(G_tmp);
                if (G_tmp % 2) {
                    if (weight >= 22) {
                        is_wasted = true;
                        break;
                    }
                    weight++;
                }
            }
            if (is_wasted) {
                continue;
            }
            count[weight]++;
        }
    }
    double end = omp_get_wtime();

    printf("Use Time:%f\n", end - start);
    for (int i = 0; i < 23; i++) {
        if (count[i]) {
            printf("%d: %d\n", i, count[i]);
        }
    }

    //     printf("Threads: %d\n", Eigen::nbThreads());
}