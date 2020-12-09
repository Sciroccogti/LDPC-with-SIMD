#include "main.hpp"

#include <omp.h>

#include <iostream>

#include "Alist/Alist.hpp"
#include "Combine/Combine.hpp"
#include "MatrixMath/NoramlMath.hpp"
#include "MatrixMath/SIMDMath.hpp"
// #include "LDPC/LDPC.hpp"

int main(int argc, char* argv[]) {
    // Eigen::initParallel();
    char *alist_path = NULL, *output_path = NULL;
    bool enable_SIMD = true;

    if (opt(argc, argv, alist_path, enable_SIMD, output_path)) {
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

    int count[129] = {0};
    readOutput(output_path, count);

    double start = omp_get_wtime();
    if (enable_SIMD) {
        // store a 64 bit with double
        u_int64_t G_int[N];
        b_type G_vec[N];
        // size for which the vectorization is possible
        size_t vec_size = N - N % inc;

        // store a row as a uint64, and every 4 uint64 as a batch
        // outerSize: N
        // https://eigen.tuxfamily.org/dox/group__TutorialSparse.html
        for (int k = 0; k < G.outerSize(); k++) {
            // bitset to turn bits into uint64
            // https://www.cplusplus.com/reference/bitset/bitset/
            std::bitset<M> col(0);
            for (Eigen::SparseMatrix<int>::InnerIterator it(G, k); it; ++it) {
                col[it.row()] = it.value();
            }
            G_int[k] = col.to_ullong();
        }

#pragma omp parallel for
        for (int i = 0; i < vec_size; i += inc) {
            G_vec[i / inc] = _mm256_load_si256((__m256i*)&G_int[i]);
        }

        // vector_type two(inc, 2);
        // b_type two_vec = &two[0];

        start = omp_get_wtime();
#pragma omp parallel for reduction(+ : count[:129])
        for (size_t i = 0; i < 4096000; i++) {
            int weight = 0;
            __attribute__((aligned(32)))
            u_int64_t v0[8] = {i, i, i, i, i, i, i, i};
            b_type I_vec = _mm256_load_si256((__m256i*)&v0[0]);

#pragma omp parallel for
            for (int j = 0; j < vec_size; j += inc) {
                b_type tmp_vec = G_vec[j / inc] & I_vec;
                hamming(tmp_vec);
                for (size_t k = 0; k < inc; k++) {
                    if (tmp_vec[k] % 2) {
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
                col[it.row()] = it.value();
            }
            G_int[i] = col.to_ullong();
        }

        start = omp_get_wtime();
#pragma omp parallel for reduction(+ : count[:129])
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
    writeOutput(output_path, count);

    //     printf("Threads: %d\n", Eigen::nbThreads());
}