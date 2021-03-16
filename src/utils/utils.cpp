#include "utils/utils.hpp"

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
        {"factor", required_argument, NULL, 'f' + 128},  // +128 to avoid ascii
        {"SNR", required_argument, NULL, 'S' + 128},
        {"iter", required_argument, NULL, 'i'},
        {"FE", required_argument, NULL, 'e'},
        {0, 0, 0, 0}};
    static char* const short_options = (char*)"hH:s:o:M:i:e:";

    while ((ret = getopt_long(argc, argv, short_options, long_options,
                              &option_index)) != -1) {
        switch (ret) {
            case 'h':
                printf("usage: LPDC-with-SIMD [-h] [-H DEC_H_PATH]\n");
                printf("\nrequired arguments:\n");
                printf("  --dec-h-path, -H <file [read only]>\n");
                printf("\t\t\tpath to H matrix in .alist\n");
                printf("\noptional arguments:\n");
                printf("  --simd, -s <bool: ON/1/OFF/0>\n");
                printf(
                    "\t\t\twhether to enable SIMD, default to be "
                    "ON\n");
                printf("  --output, -o <file [write only]>\n");
                printf(
                    "\t\t\tpath to output in .yaml,"
                    " default to be output.yaml\n");
                printf("  --MIPP, -M <bool: ON/1/OFF/0>\n");
                printf(
                    "\t\t\twhether to enable MIPP for SIMD,"
                    " default to be OFF\n");
                printf(
                    "  --factor <double>"
                    "\tthe normalize factor for NMS,"
                    " should not greater than 1,"
                    " default to be 1.0\n");
                printf(
                    "  --SNR <double>"
                    "\tthe Eb/N0 of AWGN channel,"
                    " default to be 0.0\n");
                printf(
                    "  --iter, -i <int>"
                    "\tthe Num. of iterations,"
                    " default to be 30\n");
                printf(
                    "  --FE, -e <int>"
                    "\tFrame error count,"
                    " default to be 100\n");
                printf(
                    "  --help, -h\t\tshow this help message and "
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
            case 'f' + 128:
                conf.factor = atof(optarg);
                break;
            case 'S' + 128:
                conf.SNR = atof(optarg);
                break;
            case 'i':
                conf.iter_max = atoi(optarg);
                break;
            case 'e':
                conf.FEcount = atoi(optarg);
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

    writeOutput(conf.output_path, conf);

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

void writeOutput(const char* filename, Config& conf) {
    YAML::Node yaml;
    YAML::Node basic;
    basic["alist_path"] = conf.alist_path;
    basic["enable_SIMD"] = conf.enable_SIMD;
    basic["output_path"] = conf.output_path;
    basic["enable_MIPP"] = conf.enable_MIPP;
    yaml["basic"] = basic;

    YAML::Node B_LDPC;
    B_LDPC["factor"] = conf.factor;
    B_LDPC["SNR"] = conf.SNR;
    B_LDPC["iter_max"] = conf.iter_max;
    B_LDPC["FEcount"] = conf.FEcount;
    yaml["B_LDPC"] = B_LDPC;

    std::ofstream fout(filename);
    fout << yaml;
    fout.close();
}