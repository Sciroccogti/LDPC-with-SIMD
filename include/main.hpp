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
#include <yaml-cpp/yaml.h>

#include <fstream>

struct Config {
    char* alist_path;   // path to the .alist which stores the matrix
    bool enable_SIMD;   // whether to enable SIMD
    char* output_path;  // path to the .yaml which stores the output
    bool enable_MIPP;   // whether to enable MIPP for SIMD
};

/**
 * @brief get options
 *
 * @param argc
 * @param argv
 * @param alist_path
 * @param enable_SIMD
 * @param output_path
 * @return int : return 0 if all right, -1 if something wrong
 */
int opt(int argc, char* argv[], Config& conf) {
    int opt;
    int option_index = 0;
    int ret = 0;
    bool mipp;

    static struct option long_options[] = {
        {"help", no_argument, NULL, 'h'},
        {"dec-h-path", required_argument, NULL, 'H'},
        {"simd", required_argument, NULL, 's'},
        {"output", required_argument, NULL, 'o'},
        {"MIPP", required_argument, NULL, 'M'},
        {0, 0, 0, 0}};
    static char* const short_options = (char*)"hH:s:o:M:";

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
                printf(
                    "\t\t\twhether to enable SIMD, default to be "
                    "ON\n");
                printf("  -o OUTPUT_PATH, --output OUTPUT_PATH\n");
                printf(
                    "\t\t\tpath to output in .yaml,"
                    " default to be output.yaml\n");
                printf("  -M ON/1/OFF/0, --MIPP ON/1/OFF/0\n");
                printf(
                    "\t\t\twhether to enable MIPP for SIMD,"
                    " default to be OFF\n");
                printf(
                    "  -h, --help\t\tshow this help message and "
                    "exit\n");
                return -1;
            case 'H':
                conf.alist_path = optarg;
                printf("parity matrix: %s\n", conf.alist_path);
                break;
            case 's':
                if (strcasecmp(optarg, "OFF") || strcasecmp(optarg, "0")) {
                    conf.enable_SIMD = false;
                }
                break;
            case 'o':
                conf.output_path = optarg;
                printf("save output to: %s\n", conf.output_path);
                break;
            case 'M':
                if (strcasecmp(optarg, "ON") || strcasecmp(optarg, "1")) {
                    conf.enable_MIPP = true;
                }
                mipp = conf.enable_MIPP;
#if mipp
#define ENABLE_MIPP 1
#else
#define ENABLE_MIPP 0
#endif
                break;
            default:
                break;
        }
    }

    if (conf.alist_path == NULL) {
        printf("Please specify H matrix with argument -H!\n");
        return -1;
    }

    if (conf.enable_SIMD) {
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

    if (conf.output_path == NULL) {
        conf.output_path = (char*)"output.yaml";
        printf("save output default to: %s\n", conf.output_path);
    }

    return 0;
}

/**
 * @brief read output from savefile
 *
 * @param filename path to savefile
 * @param count[22] should be initialized to all zero beforehead
 */
void readOutput(const char* filename, int* count) {
    try {
        YAML::Node yaml = YAML::LoadFile(filename);
        for (int i = 0; i < 23; i++) {
            count[i] = yaml["count"][std::to_string(i)].as<int>();
        }
    } catch (const std::exception& e) {
        YAML::Node yaml;
        YAML::Node countNode;
        for (int i = 0; i < 23; i++) {
            countNode[std::to_string(i)] = 0;
        }
        yaml["count"] = countNode;

        std::ofstream fout(filename);
        fout << yaml;
        fout.close();
    }
}

void writeOutput(const char* filename, int* count) {
    YAML::Node yaml;
    YAML::Node countNode;
    for (int i = 0; i < 23; i++) {
        countNode[std::to_string(i)] = count[i];
    }
    yaml["count"] = countNode;

    std::ofstream fout(filename);
    fout << yaml;
    fout.close();
}

#endif