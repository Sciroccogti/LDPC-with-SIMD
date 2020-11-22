#include <getopt.h>
#include <omp.h>

#include <bitset>
#include <iostream>

#include "Alist/Alist.hpp"
#include "xsimd/xsimd.hpp"
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
    const size_t M = 64;
    Eigen::SparseMatrix<int> G = G_alist.getMat();
    // Eigen::Vector<bool, M> I();
    xsimd::batch<double, 4> a(1.5, 2.5, 3.5, 4.5);
    xsimd::batch<double, 4> b(2.5, 3.5, 4.5, 5.5);
    auto mean = (a + b) / 2;
    std::cout << mean << std::endl;
    double start = omp_get_wtime();
    // u_long I_ulong = 1001;
    // printf("%ld\n", I_ulong);
    // std::bitset<M> I(I_ulong);
    // printf("%s\n", I.to_string().c_str());

    // #pragma omp parallel for
    //     for (size_t i = 0; i < M; i++) {
    //         for (size_t j = 0; i < 2; i++) {
    //             ;
    //         }
    //     }
    //     double end = omp_get_wtime();

    //     printf("Use Time:%f\n", end - start);
    //     printf("Threads: %d\n", Eigen::nbThreads());
}