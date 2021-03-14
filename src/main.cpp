#include "main.hpp"

#include <omp.h>

#include <algorithm>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>

#include "Alist/Alist.hpp"
#include "Channel/Channel.hpp"
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
    int K = ldpc.getK();  // length of message
    int N = ldpc.getN();  // length of message
    int FEcount = 0;
    int BEcount = 0;
    int count = 0;
    while (FEcount <= 50) {
        std::srand(time(nullptr));
        Eigen::RowVectorXi m = Eigen::RowVectorXf::Random(K).unaryExpr(
            [](const float x) { return (int)ceil(x); });  // message
        // std::cout << "message  : " << m << std::endl;

        Eigen::RowVectorXi c = ldpc.encode(m);  // code word encoded from m
        // std::cout << "code word: " << c << std::endl;

        // Modem modem(100, 4 * 100, 0.1);
        // Eigen::RowVectorXd y = modem.modulate(c);
        // Real BPSK will introduce error
        Eigen::RowVectorXd y = RetransBPSK(c).cast<double>();

        double snrdB = 0.0;                        // log Eb/N0
        double snr = pow(10, snrdB / 10) * K / N;  // linear Es/N0
        // std::cout << snr << std::endl;
        // std::cout << sqrt(1 / (snr))<<std::endl;
        Eigen::RowVectorXd y_noised = AWGN(y, snr, 2);
        // std::cout << y - y_noised << std::endl;

        // std::stringstream ss;
        // ss << c;
        // int length = modem.getL() * c.size();
        // Eigen::RowVectorXd x = Eigen::RowVectorXd::LinSpaced(length, 0,
        // length); double *x_array = (double *)malloc(length * sizeof(double)),
        //        *y_array = (double *)malloc(length * sizeof(double));
        // Eigen::RowVectorXd::Map(x_array, x.cols()) = x;
        // // Eigen::RowVectorXd::Map(y_array, y.cols()) = y;
        // // if (pyplot(x_array, y_array, length, ss.str().c_str())) {
        // //     printf("ERROR during pyplot!\n");
        // // }

        // Eigen::RowVectorXd::Map(y_array, y_noised.cols()) = y_noised;
        // if (pyplot(x_array, y_array, length, ss.str().c_str())) {
        //     printf("ERROR during pyplot!\n");
        // }

        // free(x_array);
        // free(y_array);

        // std::cout << y << std::endl;

        // Eigen::RowVectorXd r = modem.demodulate(y);
        Eigen::RowVectorXd r = y_noised;
        // std::cout << "received: " << r << std::endl;

        Eigen::RowVectorXi d = ldpc.decode(r, 30);
        // std::cout << "decoded: " << d << std::endl;

        // std::cout << Compare(c, d) << std::endl;

        Eigen::RowVectorXi m_ = ldpc.recoverMessage(d);
        // std::cout << "message  : " << m_ << std::endl;
        // std::cout << Compare(m, m_) << std::endl;
        int BE = Compare(m, m_);
        if (BE) {
            FEcount++;
        }
        BEcount += BE;
        count++;
    }
    printf("BER: %.2e\n", (double)BEcount / (count * K));
    printf("FER: %.2e\n", (double)FEcount / count);
}