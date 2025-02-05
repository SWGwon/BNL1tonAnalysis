#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TROOT.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <getopt.h> // for getopt_long
#include <unistd.h> // for getopt_long

#include <algorithm>
#include <cmath>
#include <cstdlib> // for exit()
#include <filesystem>
#include <fstream>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace fs = std::filesystem;
const double ADC_TO_MV = 2000.0 / (std::pow(2, 14) - 1);
const int SAMPLE_TO_NS = 2;
const int cut_pC = 15;
const int cut_pC_tp = 15;

TH1D *hist_bp1_area = nullptr;

// 디렉토리를 스캔하여 "muon"으로 시작하고 ".root"로 끝나는 파일 목록을 반환하는
// 함수
std::vector<std::string> getFileList(const std::string &directoryPath,
                                     const std::string &searchString = "") {
    std::vector<std::string> fileList;

    try {
        for (const auto &entry : fs::directory_iterator(directoryPath)) {
            if (entry.is_regular_file()) {
                std::string fileName = entry.path().filename().string();

                // 파일명이 "muon"으로 시작하고, 추가 문자열을 포함하며,
                // 확장자가
                // ".root"인지 확인
                if ((fileName.find("injection") == 0 ||
                     fileName.find("muon") == 0) &&
                    fileName.find(searchString) != std::string::npos &&
                    entry.path().extension() == ".root") {
                    fileList.push_back(entry.path().stem().string());
                }
            }
        }
    } catch (const fs::filesystem_error &e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return fileList;
}

// median 계산 함수
double getMedian(std::vector<unsigned short> &arr) {
    size_t size = arr.size();
    if (size == 0) {
        throw std::domain_error("Array is empty");
    }
    std::sort(arr.begin(), arr.end());
    return (size % 2 == 0) ? (arr[size / 2 - 1] + arr[size / 2]) / 2.0
                           : arr[size / 2];
}

// 함수: 데이터 로드
std::vector<int> GetCrossingMuonEventIds(const std::string &fpath, int config) {
    hist_bp1_area = new TH1D("hist_bp1_area", "hist_bp1_area", 100, 0, 100);
    std::vector<int> crossingMuonEventIds = {};
    TFile *file = TFile::Open(fpath.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "파일을 찾을 수 없습니다: " << fpath << std::endl;
        return {};
    }

    // List the contents of the file for debugging
    std::cout << "Listing contents of file: " << fpath << std::endl;
    file->ls();

    TTree *tree = (TTree *)file->Get("daq");
    if (!tree) {
        std::cerr << "Failed to get TTree 'daq' from file: " << fpath
                  << std::endl;
        file->Close();
        return {};
    }
    UInt_t event_id;
    tree->SetBranchAddress("event_id", &event_id);
    unsigned short adc_b1_ch0[2000] = {};
    tree->SetBranchAddress("adc_b1_ch0", adc_b1_ch0);
    // unsigned short adc_b5_ch33[2000] = {};
    // tree->SetBranchAddress("adc_b5_ch33", adc_b5_ch33);

    int nentries = tree->GetEntries();
    for (int i = 0; i < nentries; i++) {
        tree->GetEntry(i);
        if (tree->GetEntry(i) < 0) {
            std::cout << fpath << ", is broken" << std::endl;
            break;
        }

        // 그냥 지역 벡터 변수를 만든 뒤 배열의 내용을 복사
        std::vector<unsigned short> adc_b1_ch0_vec(adc_b1_ch0,
                                                   adc_b1_ch0 + 2000);
        // std::vector<unsigned short> adc_b5_ch33_vec(adc_b5_ch33, adc_b5_ch33
        // + 2000);

        double median = getMedian(adc_b1_ch0_vec);
        // double median_tp = getMedian(adc_b5_ch33_vec);

        double area_bp1 = 0;
        for (auto adc_value : adc_b1_ch0_vec) {
            area_bp1 += (median - adc_value) * ADC_TO_MV * 2 / 50;
        }
        hist_bp1_area->Fill(area_bp1);

        if (area_bp1 > cut_pC) {
            std::cout << "event_id: " << event_id << std::endl;
            std::cout << "area_bp1: " << area_bp1 << std::endl;
            // std::cout << "area_tp: " << area_tp << std::endl;
            crossingMuonEventIds.push_back(event_id);
        }
    }
    file->Close();
    return crossingMuonEventIds;
}

// 각 파일을 처리하는 함수 (스레드 또는 싱글 스레드용)
void processFiles(const std::string &directory,
                  const std::string &outputDirectory,
                  const std::vector<std::string> &fileNames) {
    for (const auto &fileName : fileNames) {
        try {
            std::string fpath = directory + fileName + ".root";
            std::cout << "file: " << fpath << std::endl;
            std::vector<int> crossingMuonEventIds =
                GetCrossingMuonEventIds(fpath, 0);

            if (crossingMuonEventIds.size() == 0)
                continue;

            std::ofstream outputFile(outputDirectory + fileName +
                                     "_good_eventId_with_paddle_CXX.txt");
            for (const auto &id : crossingMuonEventIds) {
                outputFile << id << std::endl;
            }
            outputFile.close();
            TCanvas can;
            can.SetLogy();
            hist_bp1_area->Draw();
            std::string a = outputDirectory + fileName + "_bp1_area.pdf";
            can.SaveAs(a.c_str());
            delete hist_bp1_area;

        } catch (...) {
            continue;
        }
    }
}

void createDirectoryIfNotExists(const std::string &outputDirectory) {
    // Check if the directory exists
    if (!std::filesystem::exists(outputDirectory)) {
        // Create the directory
        try {
            std::filesystem::create_directories(outputDirectory);
            std::cout << "Directory created: " << outputDirectory << std::endl;
        } catch (const std::filesystem::filesystem_error &e) {
            std::cerr << "Error creating directory: " << e.what() << std::endl;
        }
    } else {
        std::cout << "Directory already exists: " << outputDirectory
                  << std::endl;
    }
}

// 메인 함수
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

    std::string directoryPath =
        "/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/data/RAW_" +
        liquid_config + "/";
    std::string outputDirectory =
        "/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/data/crossing_muon_" +
        liquid_config + "/";
    createDirectoryIfNotExists(outputDirectory);

    std::vector<std::string> fileNames = getFileList(directoryPath, timeStamp);

    processFiles(directoryPath, outputDirectory, fileNames);

    std::cout << "All done." << std::endl;
    return 0;
}
