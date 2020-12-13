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
#include "Modem/Modem.hpp"
#include "UI/pyplot.hpp"

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
    std::cout << "message  : " << m << std::endl;

    Eigen::RowVectorXi c = ldpc.encode(m);
    std::cout << "code word: " << c << std::endl;

    Modem modem(100, 4 * 100, 0.1);
    int length = modem.getL() * c.size();
    Eigen::RowVectorXd x = Eigen::RowVectorXd::LinSpaced(length, 0, length);
    Eigen::RowVectorXd y = modem.modulate(c);
    double *x_array = (double *)malloc(length * sizeof(double)),
           *y_array = (double *)malloc(length * sizeof(double));
    Eigen::RowVectorXd::Map(x_array, x.cols()) = x;
    Eigen::RowVectorXd::Map(y_array, y.cols()) = y;
    if (pyplot(x_array, y_array, length)) {
        printf("ERROR during pyplot!\n");
    }
    free(x_array);
    free(y_array);
}