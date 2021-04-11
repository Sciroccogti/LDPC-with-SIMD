#include <omp.h>

#include <GF.hpp>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <mutex>
#include <random>
#include <sstream>
#include <thread>
#include <vector>

#include "Alist/Alist.hpp"
#include "Channel/Channel.hpp"
#include "Combine/Combine.hpp"
#include "LDPC/LDPC.hpp"
#include "LDPC/NBLDPC.hpp"
#include "MatrixMath/NoramlMath.hpp"
#include "MatrixMath/SIMDMath.hpp"
#include "Modem/Modem.hpp"
#include "UI/pyplot.hpp"
#include "utils/utils.hpp"

std::mutex mtx;

void decode(const LDPC *ldpc, const Config *conf, const double SNR, int *count,
            int *BEcount, int *FEcount);

void NBdecode(const NBLDPC *NBldpc, const Config *conf, const double SNR,
              int *count, int *BEcount, int *FEcount);

int main(int argc, char *argv[]) {
    Config conf = {.alist_path = NULL,
                   .enable_SIMD = true,
                   .output_path = NULL,
                   .threads = std::thread::hardware_concurrency() - 1,
                   .enable_MIPP = false,
                   .factor = 1.0,
                   .SNRmin = 0.0,
                   .SNRmax = 3.0,
                   .SNRstep = 0.5,
                   .iter_max = 30,
                   .FEcount = 100,
                   .mode = BP_NMS};

    if (opt(argc, argv, conf)) {
        return -1;
    }
#ifdef MIPP_ALIGNED_LOADS
    printf("aligned!\n");
#endif

    // Alist<nbalist_matrix> A = Alist<nbalist_matrix>(conf.alist_path);
    // std::cout << A.getMat() << std::endl;

    // LDPC ldpc(conf.alist_path);
    NBLDPC NBldpc(conf.alist_path);
    // std::cout << "G:\n" << NBldpc.getG() << std::endl;
    // std::cout << "H:\n" << NBldpc.getH() << std::endl;
    printf("%s read in successfully.\n", conf.alist_path);
    int K = NBldpc.getK();  // length of message
    // printf("K = %d\n", K);
    // int N = NBldpc.getN();  // length of message
    // printf("N = %d\n", N);

    // int FEcount = 0;
    // int BEcount = 0;
    // int count = 0;
    // NBdecode(&NBldpc, &conf, 0, &count, &BEcount, &FEcount);

    for (double SNR = conf.SNRmin; SNR <= conf.SNRmax; SNR += conf.SNRstep) {
        int FEcount = 0;
        int BEcount = 0;
        int count = 0;

        std::vector<std::thread> threads_(conf.threads);
        for (int i = 0; i < conf.threads; i++) {
            threads_[i] = std::thread(NBdecode, &NBldpc, &conf, SNR, &count,
                                      &BEcount, &FEcount);
        }

        auto start = std::chrono::steady_clock::now();
        for (int i = 0; i < conf.threads; i++) {
            threads_[i].join();
        }
        auto end = std::chrono::steady_clock::now();
        double BER = (double)BEcount / (count * K);
        double FER = (double)FEcount / count;
        // count is ns
        double duration = (end - start).count() / 1000000000.0;
        printf("\nBER: %.2e\n", BER);
        printf("FER: %.2e\n", FER);
        printf("Time: %.2f sec\n", duration);
        writeResult(conf.output_path, SNR, BER, FER, duration);
        for (int i = 0; i < conf.threads; i++) {
            threads_[i].~thread();
        }
    }
    return 0;
}

void decode(const LDPC *ldpc, const Config *conf, const double SNR, int *count,
            int *BEcount, int *FEcount) {
    int K = ldpc->getK();  // length of message
    int N = ldpc->getN();  // length of message

    // printf("%lf\n", SNR);
    double snr = pow(10, SNR / 10) * K / N;  // linear Es/N0

    std::srand(time(nullptr) +
               std::hash<std::thread::id>()(std::this_thread::get_id()));
    std::default_random_engine engine(
        time(nullptr) +
        std::hash<std::thread::id>()(std::this_thread::get_id()));

    mtx.lock();
    while (*FEcount <= conf->FEcount) {
        mtx.unlock();
        Eigen::RowVectorXi m = Eigen::RowVectorXf::Random(K).unaryExpr(
            [](const float x) { return (int)ceil(x); });  // message
        // std::cout << "message  : " << m << std::endl;

        Eigen::RowVectorXi c = ldpc->encode(m);  // code word encoded from m
        // std::cout << "code word: " << c << std::endl;

        // Modem modem(100, 4 * 100, 0.1);
        // Eigen::RowVectorXd y = modem.modulate(c);
        // Real BPSK will introduce error
        Eigen::RowVectorXd y = RetransBPSK(c).cast<double>();

        Eigen::RowVectorXd y_noised = AWGN(y, snr, 2, engine);
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

        Eigen::RowVectorXi d =
            ldpc->decode(r, conf->iter_max, conf->factor, snr, conf->mode);
        // std::cout << "decoded: " << d << std::endl;

        Eigen::RowVectorXi m_ = ldpc->recoverMessage(d);
        // std::cout << "message  : " << m_ << std::endl;
        // std::cout << Compare(m, m_) << std::endl;
        int BE = Compare(m, m_);

        mtx.lock();
        if (*FEcount >= conf->FEcount) {
            break;
        }
        if (BE) {
            (*FEcount)++;
            // printf("\rFEcount = %3d", *FEcount);
        }
        *BEcount += BE;
        (*count)++;
    }
    mtx.unlock();
}

void NBdecode(const NBLDPC *NBldpc, const Config *conf, const double SNR,
              int *count, int *BEcount, int *FEcount) {
    int K = NBldpc->getK();  // length of message
    int N = NBldpc->getN();  // length of message
    int GF = NBldpc->getGF();

    // printf("%lf\n", SNR);
    double snr = pow(10, SNR / 10) * K / N;  // linear Es/N0

    std::srand(time(nullptr) +
               std::hash<std::thread::id>()(std::this_thread::get_id()));
    std::default_random_engine engine(
        time(nullptr) +
        std::hash<std::thread::id>()(std::this_thread::get_id()));

    mtx.lock();
    while (*FEcount <= conf->FEcount) {
        mtx.unlock();
        // 0 < randf < 2
        Eigen::RowVectorXf randf =
            Eigen::RowVectorXf::Random(K) + Eigen::RowVectorXf::Ones(K);
        randf *= (GF * 0.5);  // 0 < randf < GF
        // message
        Eigen::RowVectorXi m =
            randf.unaryExpr([](const float x) { return (int)floor(x); });
        std::cout << "message  : " << m << std::endl;

        Eigen::RowVectorXi c = NBldpc->encode(m);  // code word encoded from m
        std::cout << "code word: " << c << std::endl;

        // // Modem modem(100, 4 * 100, 0.1);
        // // Eigen::RowVectorXd y = modem.modulate(c);
        // // Real BPSK will introduce error
        // Eigen::RowVectorXd y = RetransBPSK(c).cast<double>();
        Eigen::RowVectorXd y = RetransBPSK(NB2Bin(c, GF)).cast<double>();
        std::cout << y << std::endl;
        Eigen::RowVectorXi test = Bin2GF(TransBPSK(y.cast<int>()), GF);
        std::cout << test << std::endl;

        Eigen::RowVectorXd y_noised = AWGN(y, snr, 2, engine);
        // std::cout << y - y_noised << std::endl;

        // // std::stringstream ss;
        // // ss << c;
        // // int length = modem.getL() * c.size();
        // // Eigen::RowVectorXd x = Eigen::RowVectorXd::LinSpaced(length, 0,
        // // length); double *x_array = (double *)malloc(length *
        // sizeof(double)),
        // //        *y_array = (double *)malloc(length * sizeof(double));
        // // Eigen::RowVectorXd::Map(x_array, x.cols()) = x;
        // // // Eigen::RowVectorXd::Map(y_array, y.cols()) = y;
        // // // if (pyplot(x_array, y_array, length, ss.str().c_str())) {
        // // //     printf("ERROR during pyplot!\n");
        // // // }

        // // Eigen::RowVectorXd::Map(y_array, y_noised.cols()) = y_noised;
        // // if (pyplot(x_array, y_array, length, ss.str().c_str())) {
        // //     printf("ERROR during pyplot!\n");
        // // }

        // // free(x_array);
        // // free(y_array);

        // // std::cout << y << std::endl;

        // // Eigen::RowVectorXd r = modem.demodulate(y);
        Eigen::RowVectorXd r = y_noised;
        std::cout << "received: " << r << std::endl;

        Eigen::MatrixXd LLR = LLR_BinAWGN2GF(r, GF, snr);

        Eigen::RowVectorXi d =
            NBldpc->decode(LLR, conf->iter_max, conf->factor, snr, conf->mode);
        std::cout << "decoded: " << d << std::endl;

        // // Eigen::RowVectorXi m_ = NBldpc->recoverMessage(d);
        // // // std::cout << "message  : " << m_ << std::endl;
        // // // std::cout << Compare(m, m_) << std::endl;
        // // int BE = Compare(m, m_);

        // mtx.lock();
        // if (*FEcount >= conf->FEcount) {
        //     break;
        // }
        // // if (BE) {
        // (*FEcount)++;
        // //     // printf("\rFEcount = %3d", *FEcount);
        // // }
        // // *BEcount += BE;
        // (*count)++;
        break;
    }
    mtx.unlock();
}
