#include <getopt.h>
#include <omp.h>

#include <bitset>
#include <iostream>
#include <vector>

#include "Alist/Alist.hpp"
#include "xsimd/xsimd.hpp"
// #include "LDPC/LDPC.hpp"
xsimd::batch<u_int64_t, 4> hamming(xsimd::batch<u_int64_t, 4> n) {
    n = (n & 0x5555555555555555) + ((n >> 1) & 0x5555555555555555);
    n = (n & 0x3333333333333333) + ((n >> 2) & 0x3333333333333333);
    n = (n & 0x0f0f0f0f0f0f0f0f) + ((n >> 4) & 0x0f0f0f0f0f0f0f0f);
    n = (n & 0x00ff00ff00ff00ff) + ((n >> 8) & 0x00ff00ff00ff00ff);
    n = (n & 0x0000ffff0000ffff) + ((n >> 16) & 0x0000ffff0000ffff);
    n = (n & 0x00000000ffffffff) + ((n >> 32) & 0x00000000ffffffff);
    return n;
}
int main(int argc, char* argv[]) {
    // Eigen::initParallel();
#ifdef __AVX512F__
    printf("AVX!\n");
#endif
    int opt;
    static struct option long_options[] = {
        {"help", no_argument, NULL, 'h'},
        {"file", required_argument, NULL, 'H'},
        {"output", optional_argument, NULL, 'o'},
        {0, 0, 0, 0}};
    static char* const short_options = (char*)"hH:o::";

    int option_index = 0;
    int ret = 0;

    char* alist_path = NULL;

    while ((ret = getopt_long(argc, argv, short_options, long_options,
                              &option_index)) != -1) {
        switch (ret) {
            case 'h':
                printf("usage: LPDC-with-SIMD [-h] [-H DEC_H_PATH]\n");
                printf("\nrequired arguments:\n");
                printf("  -H DEC_H_PATH, --dec-h-path DEC_H_PATH\n");
                printf("\t\t\tpath to H matrix in .alist\n");
                printf("\noptional arguments:\n");
                printf("  -h, --help\t\tshow this help message and exit\n");
                return 0;
            case 'H':
                alist_path = optarg;
                printf("parity matrix: %s\n", alist_path);
                break;
            case 'o':
                printf("HAVE option: -c\n");
                printf("The argument of -c is %s\n\n", optarg);
                break;
            default:
                break;
        }
    }

    if (alist_path == NULL) {
        printf("Please specify H matrix with argument -H!\n");
        return 1;
    }
    // xsimd::batch<u_int64_t, 4> a(1, 1, 1, 1);
    // xsimd::batch<u_int64_t, 4> b(-74, 2, 2, 2);
    // std::cout<<b<<std::endl;
    // a^=b;
    // std::cout<<a<<std::endl;
    // a = hamming(a);
    // std::cout<<a<<std::endl;
    // return 0;
    // LDPC ldpc(alist_path);

    // Eigen::SparseMatrix<int> H = ldpc.getH(), G = ldpc.getG();
    // std::cout << G << H << std::endl;
    // std::cout << G * H.transpose() << std::endl;
    Alist<alist_matrix> G_alist(alist_path);
    const size_t M = 64, N = 128;
    Eigen::SparseMatrix<int> G = G_alist.getMat();
    std::vector<xsimd::batch<u_int64_t, 4>> G_int;  // store a 64 bit with 2 int

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

    int count[129] = {0};
    double start = omp_get_wtime();

#pragma omp parallel for
    for (size_t i = 0; i < 4096000; i++) {
        u_int64_t I0 = i;

        xsimd::batch<u_int64_t, 4> I(I0, I0, I0, I0);
        xsimd::batch<u_int64_t, 4> two(2, 2, 2, 2);
        int weight = 0;

        for (int j = 0; j < N / 4; j++) {
            xsimd::batch<u_int64_t, 4> G_tmp = G_int[j];
            G_tmp &= I;
            G_tmp = hamming(G_tmp);
            xsimd::batch<u_int64_t, 4> bit = G_tmp % two;
            for (size_t k = 0; k < 4; k++) {
                if (bit[k] != 0) {
                    weight++;
                }
            }
        }
        count[weight]++;
    }
    double end = omp_get_wtime();

    printf("Use Time:%f\n", end - start);
    for (int i = 0; i < 129; i++)
    {
        printf("%d: %d\n", i, count[i]);
    }
    
    //     printf("Threads: %d\n", Eigen::nbThreads());
}