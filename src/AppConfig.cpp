#include "AppConfig.h"
#include <cstdlib>
#include <getopt.h>
#include <iostream>

AppConfig::AppConfig()
    : inputFileName(""), outputFilePath(""),
      inputSPECalibrationPath("/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/"
                              "yaml/bnl1t_spe_fit_results_241111.csv"),
      eventNumber(10000), triggerType(0), sampleSize(1920) // default triggerType: top paddle
{}

AppConfig AppConfig::parseArgs(int argc, char **argv) {
    AppConfig config;
    int c;
    while (true) {
        static struct option long_options[] = {
            {"inputFileName", required_argument, 0, 'i'},
            {"outputFilePath", required_argument, 0, 'o'},
            {"eventNumber", required_argument, 0, 'n'},
            {"sampleSize", required_argument, 0, 's'},
            {"triggerType", required_argument, 0, 't'},
            {"inputSPECalibrationPath", required_argument, 0, 'c'},
            {0, 0, 0, 0}};

        int option_index = 0;
        c = getopt_long(argc, argv, "i:o:n:t:c:s:", long_options, &option_index);
        if (c == -1)
            break;

        switch (c) {
            case 'i':
                config.inputFileName = optarg;
                break;
            case 'o':
                config.outputFilePath = optarg;
                break;
            case 'n':
                config.eventNumber = std::stoi(optarg);
                break;
            case 's':
                config.sampleSize = std::stoi(optarg);
                break;
            case 't':
                config.triggerType = std::stoi(optarg);
                break;
            case 'c':
                config.inputSPECalibrationPath = optarg;
                break;
            default:
                std::cout << "Usage: " << argv[0] << " [OPTIONS]" << std::endl;
                std::cout << "Options:" << std::endl;
                std::cout << "  -i <inputFilePath>          Path to the input file."
                    << std::endl;
                std::cout
                    << "  -o <outputFilePath>         Path to the output file."
                    << std::endl;
                std::cout
                    << "  -n <eventNumber>            Event number to process."
                    << std::endl;
                std::cout
                    << "  -n <sampleSize>            Number of samples (b1 ~ b4, time window / 2)."
                    << std::endl;
                std::cout << "  -t <triggerType>            Trigger type:"
                    << std::endl;
                std::cout << "                              0: tp" << std::endl;
                std::cout << "                              1: alpha" << std::endl;
                std::cout << "                              2: majority" << std::endl;
                std::cout << "                              3: crossing muon" << std::endl;
                std::cout << "                              4: crossing muon (before majority setup)" << std::endl;
                std::cout << "  -c <inputSPECalibrationPath>  Path to the SPE "
                    "calibration file."
                    << std::endl;
                exit(EXIT_FAILURE);
        }
    }
    if (config.inputFileName.empty() || config.outputFilePath.empty()) {
        std::cout << "Usage: " << argv[0] << " [OPTIONS]" << std::endl;
        std::cout << "Options:" << std::endl;
        std::cout << "  -i <inputFilePath>          Path to the input file."
            << std::endl;
        std::cout
            << "  -o <outputFilePath>         Path to the output file."
            << std::endl;
        std::cout
            << "  -n <eventNumber>            Event number to process."
            << std::endl;
        std::cout
            << "  -n <sampleSize>            Number of samples (b1 ~ b4, time window / 2)."
            << std::endl;
        std::cout << "  -t <triggerType>            Trigger type:"
            << std::endl;
        std::cout << "                              0: tp" << std::endl;
        std::cout << "                              1: alpha" << std::endl;
        std::cout << "                              2: majority" << std::endl;
        std::cout << "                              3: crossing muon" << std::endl;
        std::cout << "                              4: crossing muon (before majority setup)" << std::endl;
        std::cout << "  -c <inputSPECalibrationPath>  Path to the SPE "
            "calibration file."
            << std::endl;
        exit(EXIT_FAILURE);
    }
    return config;
}
