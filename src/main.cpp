#include "main.hpp"

#include <omp.h>

#include <iostream>
#include <vector>

#include "Alist/Alist.hpp"
#include "Combine/Combine.hpp"
#include "LDPC/LDPC.hpp"
#include "MatrixMath/NoramlMath.hpp"
#include "MatrixMath/SIMDMath.hpp"

int main(int argc, char* argv[]) {
    Config conf = {.alist_path = NULL,
                   .enable_SIMD = true,
                   .output_path = NULL,
                   .enable_MIPP = false};

    if (opt(argc, argv, conf)) {
        return -1;
    }
#ifdef MIPP_ALIGNED_LOADS
    printf("aligned!\n");
#endif

    LDPC ldpc(conf.alist_path);
    Eigen::SparseMatrix<int> H = ldpc.getH(), G = ldpc.getG();
    std::cout << H << std::endl;
    printf("\n");
    std::cout << G << std::endl;
    // std::cout << G * H.transpose() << std::endl;
}