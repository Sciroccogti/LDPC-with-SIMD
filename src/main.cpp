#include <getopt.h>
#include <omp.h>

#include <bitset>
#include <iostream>

#include "Alist/Alist.hpp"
// #include "LDPC/LDPC.hpp"
using namespace std;

// https://stackoverflow.com/q/18838959
struct logical_xor {
    bool operator()(bool a, bool b) const {
        return a != b;
    }
};

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

    // LDPC ldpc(alist_path);

    // Eigen::SparseMatrix<int> H = ldpc.getH(), G = ldpc.getG();
    // std::cout << G << H << std::endl;
    // std::cout << G * H.transpose() << std::endl;
    Alist<alist_matrix> G_alist(alist_path);
    const size_t M = 64;
    Eigen::Array<bool, M, 128> G = G_alist.getMat().cast<bool>();
    Eigen::Array<bool, M, 1> I;
    std::cout << I.transpose() << std::endl;

    double start = omp_get_wtime();
    int count[129] = {0};
#pragma omp parallel for
    for (size_t i = 0; i < 4096000; i++) {
        I.setRandom(M);
        // I.Zero(M);
        int weight = (G.colwise() * I).colwise().redux(logical_xor()).count();
        count[weight]++;
        if (weight < 22) {
            std::cout << I.transpose() << std::endl;
        }
    }

    for (size_t i = 0; i < 129; i++) {
        printf("%ld:\t\t%d\n", i, count[i]);
    }

    // #pragma omp parallel for
    //     for (size_t i = 0; i < M; i++) {
    //         for (size_t j = 0; i < 2; i++) {
    //             ;
    //         }
    //     }
    double end = omp_get_wtime();

    printf("Use Time:%f\n", end - start);
    //     printf("Threads: %d\n", Eigen::nbThreads());
}