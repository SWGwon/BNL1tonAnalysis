#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <mutex>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include <cstdlib>  // for exit()
#include <getopt.h> // for getopt_long
#include <unistd.h> // for getopt_long

#include <TFile.h>
#include <TKey.h>
#include <TROOT.h>
#include <TTree.h>

namespace fs = std::filesystem;

// Function to collect data from text files for each channel
std::map<std::string, std::vector<double>>
collectChannelData(const std::vector<std::string> &fileNames,
                   const std::vector<std::string> &channel_names,
                   const std::string &basePath,
                   const std::string &analysisFolder) {

    std::map<std::string, std::vector<double>> channelData;

    // Initialize the map with empty vectors for each channel
    for (const auto &channel : channel_names) {
        channelData[channel] = std::vector<double>();
    }

    // Loop over all file names
    for (const auto &name : fileNames) {
        // Build the path to the analysis folder for this file
        std::string analysisPath = basePath + analysisFolder + name + "/";

        // Loop over each channel
        for (const auto &channel : channel_names) {
            std::string txtFileName =
                analysisPath + "npe_channel_" + channel + ".txt";

            std::ifstream infile(txtFileName);
            if (!infile.is_open()) {
                std::cerr << "Failed to open file " << txtFileName << std::endl;
                continue; // Skip if the file cannot be opened
            }

            double value;
            while (infile >> value) {
                channelData[channel].push_back(value);
            }
            infile.close();
        }
    }

    return channelData;
}

// 디렉토리를 스캔하여 "muon"으로 시작하고 ".root"로 끝나는 파일 목록을 반환하는
// 함수
std::vector<std::string> getFileList(const std::string &directoryPath) {
    std::vector<std::string> fileList;

    try {
        // 디렉토리 스캔
        for (const auto &entry : fs::directory_iterator(directoryPath)) {
            if (entry.is_regular_file()) {
                std::string fileName = entry.path().filename().string();
                // "muon"으로 시작하고 ".root"로 끝나는 파일만 포함
                if (fileName.find("muon") == 0 &&
                    entry.path().extension() == ".root") {
                    fileList.push_back(
                        entry.path()
                            .stem()
                            .string()); // 확장자 제외한 파일명 저장
                }
            }
        }
    } catch (const fs::filesystem_error &e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return fileList;
}

int main(int argc, char *argv[]) {
    std::string liquid_config =
        "";         // Initial value set to -1 to check if it's provided
    int config = 0; // Default value 0
    std::string timeStamp = "";

    // 옵션 파싱을 위한 설정
    struct option long_options[] = {
        {"liquid_config", required_argument, 0, 'L'},
        {"config", optional_argument, 0, 'C'},
        {"timestamp", optional_argument, 0, 'T'},
        {0, 0, 0, 0} // 옵션 목록 종료 표시
    };

    int opt;
    int option_index = 0;
    while ((opt = getopt_long(argc, argv, "L:C:T:", long_options,
                              &option_index)) != -1) {
        switch (opt) {
        case 'L': // --waterLevel 또는 -W
            liquid_config = optarg;
            break;
        case 'C': // --config 또는 -C
            config = std::stoi(optarg);
            break;
        case 'T':
            timeStamp = optarg;
            break;
        default:
            std::cerr << "Usage: " << argv[0]
                      << " --liquid_config <value> [--config <value>]\n";
            exit(EXIT_FAILURE);
        }
    }

    // waterLevel이 전달되지 않았을 경우 프로그램 종료
    if (liquid_config == "") {
        std::cerr << "Error: --liquid_config (-L) must be provided\n";
        exit(EXIT_FAILURE);
    }

    std::string basePath =
        "/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/"; // 기본 경로
    std::string analysisFolder = "analysis/";

    std::string directoryPath =
        basePath + "data/" + "crossing_muon_" + liquid_config + "/";

    // 외부 파일에서 파일 이름 목록을 읽어오기
    std::vector<std::string> fileNames = getFileList(directoryPath);

    std::vector<std::string> bot_and_side_channels = {
        "adc_b1_ch1",  "adc_b1_ch2",  "adc_b1_ch3",  "adc_b1_ch4",
        "adc_b1_ch5",  "adc_b1_ch6",  "adc_b1_ch7",  "adc_b1_ch8",
        "adc_b1_ch9",  "adc_b1_ch10", "adc_b1_ch11", "adc_b1_ch12",
        "adc_b1_ch13", "adc_b1_ch14", "adc_b1_ch15", "adc_b2_ch0",
        "adc_b2_ch1",  "adc_b2_ch2",  "adc_b2_ch3",  "adc_b2_ch4",
        "adc_b2_ch5",  "adc_b2_ch6",  "adc_b2_ch7",  "adc_b2_ch8",
        "adc_b2_ch9",  "adc_b2_ch10", "adc_b2_ch11", "adc_b2_ch12",
        "adc_b2_ch13", "adc_b2_ch14", "adc_b3_ch0",  "adc_b3_ch1",
        "adc_b3_ch2",  "adc_b3_ch3",  "adc_b3_ch4",  "adc_b3_ch5",
        "adc_b3_ch6",  "adc_b3_ch7",  "adc_b3_ch8",  "adc_b3_ch9",
        "adc_b3_ch10", "adc_b3_ch11", "adc_b3_ch12", "adc_b3_ch13",
        "adc_b3_ch14", "adc_b3_ch15", "adc_b4_ch0",  "adc_b4_ch1",
        "adc_b4_ch2",  "adc_b4_ch3",  "adc_b4_ch4",  "adc_b4_ch5",
        "adc_b4_ch6",  "adc_b4_ch7",  "adc_b4_ch8",  "adc_b4_ch9",
        "adc_b4_ch10", "adc_b4_ch11",
    };

    // Collect data for each channel
    std::map<std::string, std::vector<double>> channelData = collectChannelData(
        fileNames, bot_and_side_channels, basePath, analysisFolder);

    for (const auto &channel : bot_and_side_channels) {
        std::cout << "Channel: " << channel
                  << ", Data Points: " << channelData[channel].size()
                  << std::endl;
    }

    // After collecting the channelData
    // Create a ROOT file
    std::string outputFileName = "combined_data_" + liquid_config + ".root";
    TFile *outFile = new TFile(outputFileName.c_str(), "RECREATE");
    TTree *tree = new TTree("tree", "Combined Data Tree");

    // Map to hold the branch variables
    std::map<std::string, double> branchVariables;

    // Create branches for each channel
    for (const auto &channel : bot_and_side_channels) {
        branchVariables[channel] = 0.0;
        tree->Branch(channel.c_str(), &branchVariables[channel]);
    }

    // Determine the number of entries (assuming all channels have the same
    // number of entries)
    size_t nEntries = 0;
    if (!bot_and_side_channels.empty()) {
        nEntries = channelData[bot_and_side_channels[0]].size();
    }

    // Fill the tree
    for (size_t i = 0; i < nEntries; ++i) {
        for (const auto &channel : bot_and_side_channels) {
            if (channelData[channel].size() > i) {
                branchVariables[channel] = channelData[channel][i];
            } else {
                branchVariables[channel] = 0.0; // Handle missing data
            }
        }
        tree->Fill();
    }

    // Write and close the ROOT file
    tree->Write();
    outFile->Close();
    delete outFile;

    return 0;
}
