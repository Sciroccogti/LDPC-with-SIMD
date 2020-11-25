/**
 * @file main.hpp
 * @author Yifan Zhang (scirocco_gti@yeah.net)
 * @brief simple utils for main.cpp
 * @date 2020-11-24 11:25:54
 * @modified: 2020-11-24 12:11:30
 */

#ifndef MAIN_HPP
#define MAIN_HPP

#include <getopt.h>
#include <stdio.h>
#include <string.h>

// get options
int opt(int argc, char* argv[], char*& alist_path, bool& enable_SIMD) {
    int opt;
    int option_index = 0;
    int ret = 0;

    static struct option long_options[] = {
        {"help", no_argument, NULL, 'h'},
        {"dec-h-path", required_argument, NULL, 'H'},
        {"simd", required_argument, NULL, 's'},
        {0, 0, 0, 0}};
    static char* const short_options = (char*)"hH:s:";

    while ((ret = getopt_long(argc, argv, short_options, long_options,
                              &option_index)) != -1) {
        switch (ret) {
            case 'h':
                printf("usage: LPDC-with-SIMD [-h] [-H DEC_H_PATH]\n");
                printf("\nrequired arguments:\n");
                printf("  -H DEC_H_PATH, --dec-h-path DEC_H_PATH\n");
                printf("\t\t\tpath to H matrix in .alist\n");
                printf("\noptional arguments:\n");
                printf("  -s ON/1/OFF/0, --simd ON/1/OFF/0\n");
                printf("\t\t\twhether to enable SIMD");
                printf("  -h, --help\t\tshow this help message and exit\n");
                return -1;
            case 'H':
                alist_path = optarg;
                printf("parity matrix: %s\n", alist_path);
                break;
            case 's':
                if (strcasecmp(optarg, "OFF") || strcasecmp(optarg, "0")) {
                    enable_SIMD = false;
                }
                break;
            default:
                break;
        }
    }

    if (alist_path == NULL) {
        printf("Please specify H matrix with argument -H!\n");
        return -1;
    }
    if (enable_SIMD) {
#ifdef __AVX512F__
        printf("using AVX512\n");
#else
#ifdef __AVX2__
        printf("using AVX2\n");
#else
        printf("No AVX found!\n");
        return -1;
#endif
#endif
    }
    return 0;
}

#endif