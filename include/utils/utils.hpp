/**
 * @file utils.hpp
 * @author Sciroccogti (scirocco_gti@yeah.net)
 * @brief simple utils for main.cpp
 * @date 2020-11-24 11:25:54
 * @modified: 2021-04-27 13:11:00
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
    float factor;  // normalize factor for NMS
    // NB-LDPC
    int n_max;      // n_max for EMS
    float SNRmin;  // Eb/N0 of AWGN
    float SNRmax;
    float SNRstep;
    int iter_max;  // stop criterion
    int FEcount;   // Frame error count
    int mode;      // mode of decoding algorithm
};

int opt(int argc, char* argv[], Config& conf);

void readOutput(const char* filename, int* count);
void writeConf(const char* filename, Config& conf);

void writeResult(const char* filename, float SNR, float BER, float FER,
                 float duration);

#endif