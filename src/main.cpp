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

    // Eigen::SparseMatrix<int> H = ldpc.getH(), G = ldpc.getG();
    // std::cout << G << H << std::endl;
    // std::cout << G * H.transpose() << std::endl;

    Alist<alist_matrix> G_alist(alist_path);
    const size_t M = 64, N = 128;
    Eigen::SparseMatrix<int> G = G_alist.getMat();

    int count[129] = {0};
    double start = omp_get_wtime();
    if (enable_SIMD) {
        // store a 64 bit with int64
        std::vector<xsimd::batch<u_int64_t, 4>> G_int;

        // store a row as a uint64, and every 4 uint64 as a batch
        // outerSize: M
        // https://eigen.tuxfamily.org/dox/group__TutorialSparse.html
        for (int k = 0; k < G.outerSize(); k += 4) {
            u_int64_t col_int[4];
            for (int i = 0; i < 4; i++) {
                // bitset to turn bits into uint64
                // https://www.cplusplus.com/reference/bitset/bitset/
                std::bitset<M> col(0);
                for (Eigen::SparseMatrix<int>::InnerIterator it(G, k + i); it;
                     ++it) {
                    // sequence of bitset and Matrix is different
                    col[M - it.row()] = it.value();
                }
                col_int[i] = col.to_ullong();
            }
            G_int.push_back(xsimd::batch<u_int64_t, 4>(col_int));
        }

#pragma omp parallel for
        for (size_t i = 0; i < 4096000; i++) {
            xsimd::batch<u_int64_t, 4> I(i, i, i, i);
            xsimd::batch<u_int64_t, 4> two(2, 2, 2, 2);
            int weight = 0;

            for (int j = 0; j < N / 4; j++) {
                xsimd::batch<u_int64_t, 4> G_tmp = G_int[j] & I;
                G_tmp = hamming(G_tmp);
                xsimd::batch<u_int64_t, 4> bit = G_tmp % two;
#pragma omp parrallel for
                for (size_t k = 0; k < 4; k++) {
                    if (bit[k] != 0) {
                        weight++;
                    }
                }
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
                // sequence of bitset and Matrix is different
                col[M - it.row()] = it.value();
            }
            G_int[i] = col.to_ullong();
        }

#pragma omp parallel for
        for (size_t i = 0; i < 4096000; i++) {
            int weight = 0;
#pragma omp parallel for
            for (size_t j = 0; j < N; j++) {
                u_int64_t G_tmp = G_int[j] & i;
                G_tmp = hamming(G_tmp);
                if (G_tmp % 2) {
                    weight++;
                }
            }
            count[weight]++;
        }
    }
    double end = omp_get_wtime();

    printf("Use Time:%f\n", end - start);
    for (int i = 0; i < 129; i++) {
        printf("%d: %d\n", i, count[i]);
    }

    //     printf("Threads: %d\n", Eigen::nbThreads());
}