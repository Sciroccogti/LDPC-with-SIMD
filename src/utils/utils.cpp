#include "utils/utils.hpp"

#include <ctime>

#include "LDPC/Tanner.hpp"

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
        {"simd", required_argument, NULL, 's' + 128},
        {"output", required_argument, NULL, 'o'},
        {"thread", required_argument, NULL, 't'},
        {"MIPP", required_argument, NULL, 'M' + 128},
        {"factor", required_argument, NULL, 'f' + 128},  // +128 to avoid ascii
        {"SNRmin", required_argument, NULL, 'm'},
        {"SNRmax", required_argument, NULL, 'M'},
        {"SNRstep", required_argument, NULL, 's'},
        {"iter", required_argument, NULL, 'i'},
        {"FE", required_argument, NULL, 'e'},
        {"dec-implem", required_argument, NULL, 'd' + 128},
        {0, 0, 0, 0}};
    static char* const short_options = (char*)"hH:o:t:m:M:s:i:e:";

    while ((ret = getopt_long(argc, argv, short_options, long_options,
                              &option_index)) != -1) {
        switch (ret) {
            case 'h':
                printf("usage: LPDC-with-SIMD [-h] [-H DEC_H_PATH]\n");
                printf("\nrequired arguments:\n");
                printf("  --dec-h-path, -H <file [read only]>\n");
                printf("\t\t\tpath to H matrix in .alist\n");
                printf("\noptional arguments:\n");
                printf("  --simd <bool: ON/1/OFF/0>\n");
                printf(
                    "\t\t\twhether to enable SIMD, default to be "
                    "ON\n");
                printf("  --output, -o <file [write only]>\n");
                printf(
                    "\t\t\tpath to output in .yaml,"
                    " default to be output.yaml\n");
                printf("  --thread, -t <int>\n");
                printf(
                    "\t\t\tthreads to use,"
                    " default to be the maximum of CPU - 1\n");
                printf("  --MIPP <bool: ON/1/OFF/0>\n");
                printf(
                    "\t\t\twhether to enable MIPP for SIMD,"
                    " default to be OFF\n");
                printf(
                    "  --factor <double>"
                    "\tthe normalize factor for NMS,"
                    " should not greater than 1,"
                    " default to be 1.0\n");
                printf(
                    "  --SNRmin, -m <double>\n"
                    "\t\t\tthe start Eb/N0 of AWGN channel,"
                    " default to be 0.0\n");
                printf(
                    "  --SNRmax, -M <double>\n"
                    "\t\t\tthe stop Eb/N0 of AWGN channel,"
                    " default to be 3.0\n");
                printf(
                    "  --SNRstep, -s <double>\n"
                    "\t\t\tthe step of Eb/N0 increment of AWGN channel,"
                    " default to be 0.5\n");
                printf(
                    "  --iter, -i <int>"
                    "\tthe Num. of iterations,"
                    " default to be 30\n");
                printf(
                    "  --FE, -e <int>"
                    "\tFrame error count,"
                    " default to be 100\n");
                printf(
                    "  --dec-implem <string>\n"
                    "\t\t\tdecoding algorithm,"
                    " choose from \"NMS\", \"SPA\", \"QSPA\","
                    " default to be \"NMS\"\n");
                printf(
                    "\n  --help, -h\t\tshow this help message and "
                    "exit\n");
                return -1;
            case 'H':
                conf.alist_path = optarg;
                printf("parity matrix: %s\n", conf.alist_path);
                break;
            case 's' + 128:
                if (strcasecmp(optarg, "OFF") || strcasecmp(optarg, "0")) {
                    conf.enable_SIMD = false;
                }
                break;
            case 'o':
                conf.output_path = optarg;
                printf("save output to: %s\n", conf.output_path);
                break;
            case 't':
                conf.threads = atoi(optarg);
                break;
            case 'M' + 128:
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
            case 'm':
                conf.SNRmin = atof(optarg);
                break;
            case 'M':
                conf.SNRmax = atof(optarg);
                break;
            case 's':
                conf.SNRstep = atof(optarg);
                break;
            case 'i':
                conf.iter_max = atoi(optarg);
                break;
            case 'e':
                conf.FEcount = atoi(optarg);
                break;
            case 'd' + 128: {
                int i = 0;
                bool isInited = false;
                for (const char* m : Modes_) {
                    if (!strcasecmp(m, optarg)) {
                        conf.mode = i;
                        isInited = true;
                        break;
                    }
                    i++;
                }
            } break;
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

    writeConf(conf.output_path, conf);

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

void writeConf(const char* filename, Config& conf) {
    YAML::Node yaml;
    YAML::Node basic;
    basic["alist_path"] = conf.alist_path;
    basic["enable_SIMD"] = conf.enable_SIMD;
    basic["output_path"] = conf.output_path;
    basic["threads"] = conf.threads;
    basic["enable_MIPP"] = conf.enable_MIPP;
    yaml["basic"] = basic;

    YAML::Node B_LDPC;
    B_LDPC["factor"] = conf.factor;
    B_LDPC["iter_max"] = conf.iter_max;
    B_LDPC["FEcount"] = conf.FEcount;
    B_LDPC["mode"] = Modes_[conf.mode];
    yaml["B_LDPC"] = B_LDPC;

    std::ofstream fout(filename);
    fout << yaml << std::endl;
    fout.close();
}

void writeResult(const char* filename, double SNR, double BER, double FER,
                 double duration) {
    YAML::Node yaml;
    YAML::Node result;
    char tmp[20], snr[20];
    sprintf(snr, "%.2lf", SNR);
    result["SNR"] = snr;
    sprintf(tmp, "%.2e", BER);
    result["BER"] = tmp;
    sprintf(tmp, "%.2e", FER);
    result["FER"] = tmp;
    sprintf(tmp, "%.2f", duration);
    result["sec"] = tmp;
    yaml[snr] = result;

    std::ofstream fout(filename, std::ios_base::app);
    fout << std::scientific << yaml << std::endl;
    fout.close();
}