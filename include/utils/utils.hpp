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
    char* alist_path;   // path to the .alist which stores the matrix
    bool enable_SIMD;   // whether to enable SIMD
    char* output_path;  // path to the .yaml which stores the output
    bool enable_MIPP;   // whether to enable MIPP for SIMD

    // B-LDPC
    double factor;  // normalize factor for NMS
    double SNR;     // Eb/N0 of AWGN
    int iter_max;   // stop criterion
    int FEcount;    // Frame error count
};

int opt(int argc, char* argv[], Config& conf);

void readOutput(const char* filename, int* count);

void writeOutput(const char* filename, Config& conf);

#endif