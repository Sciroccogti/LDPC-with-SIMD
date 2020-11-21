#include <getopt.h>
#include <omp.h>

#include <iostream>

#include "Alist/Alist.hpp"
// #include "LDPC/LDPC.hpp"
using namespace std;
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
    Eigen::SparseMatrix<int> G = G_alist.getMat();
    double start = omp_get_wtime();
#pragma omp parallel for
    for (size_t i = 0; i < 1000000000; i++) {
        G* G.transpose();
    }
    double end = omp_get_wtime();

    printf("hhUse Time:%f\n", end - start);
    printf("%d\n", Eigen::nbThreads());
}