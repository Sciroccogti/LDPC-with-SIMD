#include "main.hpp"

#include <omp.h>

#include <algorithm>
#include <iostream>
#include <random>
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
    int n_info = ldpc.getK();
    int sigma = 2;
    static std::default_random_engine e(time(0));
    static std::normal_distribution<double> normal(0, sigma);

    std::srand(time(nullptr));
    Eigen::RowVectorXi m = Eigen::RowVectorXf::Random(n_info).unaryExpr(
        [](const float x) { return (int)ceil(x); });
    // .unaryExpr(
    // [](double dummy) { return (int)normal(e); });
    std::cout << m << std::endl;
    std::cout << ldpc.encode(m) << std::endl;
}