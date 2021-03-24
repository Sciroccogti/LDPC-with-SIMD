/**
 * @file utils.hpp
 * @author Sciroccogti (scirocco_gti@yeah.net)
 * @brief simple utils for main.cpp
 * @date 2020-11-24 11:25:54
 * @modified: 2021-03-16 11:30:19
 */

#ifndef UTILS_HPP
#define UTILS_HPP

#include <getopt.h>
#include <stdio.h>
#include <string.h>
#include <yaml-cpp/yaml.h>

#include <fstream>

struct Config {
    char* alist_path;      // path to the .alist which stores the matrix
    bool enable_SIMD;      // whether to enable SIMD
    char* output_path;     // path to the .yaml which stores the output
    unsigned int threads;  // threads to use
    bool enable_MIPP;      // whether to enable MIPP for SIMD

    // B-LDPC
    double factor;  // normalize factor for NMS
    double SNRmin;     // Eb/N0 of AWGN
    double SNRmax;
    double SNRstep;
    int iter_max;   // stop criterion
    int FEcount;    // Frame error count
    int mode;       // mode of decoding algorithm
};

int opt(int argc, char* argv[], Config& conf);

void readOutput(const char* filename, int* count);
void writeConf(const char* filename, Config& conf);

void writeResult(const char* filename, double SNR, double BER, double FER,
                 double duration);

#endif