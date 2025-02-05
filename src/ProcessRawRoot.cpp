//#include <Eigen/Dense>
//#include <getopt.h>
//#include <iomanip>
//#include <iostream>
//#include <string>
//#include <vector>
//#include "TApplication.h"
//#include "TCanvas.h"
//#include "TH2D.h"
//#include "TQObject.h"
//#include "TStyle.h"
//#include <TFile.h>
//#include <TTree.h>
//#include <TGraph.h>
//#include <iostream>
//#include <fstream>
//#include <sstream>
//#include <map>
//#include <string>
//#include <TBranch.h>
//#include <TLeaf.h>
//#include <chrono>
//
//#include "Waveform.h"
//
//const int MAX_SAMPLE_SIZE = 3000;
//int event_number_to_process = 10000;
//int TRRIGER = 0; //0: top paddle, 1: alpha, 2: majority
//int INTEGRATE_RANGE_START = 0;
//int INTEGRATE_RANGE_END = 3000;
//
//std::map<std::string, float> spe_mean; // spe_mean 저장
//std::map<std::string, std::map<std::string, float>> spe_fit_results; // ch_name 기준으로 결과 저장
//std::map<std::string, std::vector<double>> amp_pe;
//std::map<std::string, std::vector<double>> amp_mV;
//std::map<std::string, std::vector<double>> pe;
//
//std::string trigtag_tpch = "adc_b5_ch33";
//std::string trigtag_alphach = "adc_b5_ch34";
//std::string trigtag_majoritych = "adc_b5_ch35";
//std::vector<std::string> triggers = {
//    trigtag_tpch, trigtag_alphach, trigtag_majoritych
//};
//
//std::vector<std::string> pmts = {
//    "adc_b1_ch1", "adc_b1_ch2", "adc_b1_ch3", "adc_b1_ch4",
//    "adc_b1_ch5", "adc_b1_ch6", "adc_b1_ch7", "adc_b1_ch8",
//    "adc_b1_ch9", "adc_b1_ch10", "adc_b1_ch11", "adc_b1_ch12",
//    "adc_b1_ch13", "adc_b1_ch14", "adc_b1_ch15", "adc_b2_ch0",
//    "adc_b2_ch1", "adc_b2_ch2", "adc_b2_ch3", "adc_b2_ch4",
//    "adc_b2_ch5", "adc_b2_ch6", "adc_b2_ch7", "adc_b2_ch8",
//    "adc_b2_ch9", "adc_b2_ch10", "adc_b2_ch11", "adc_b2_ch12",
//    "adc_b2_ch13", "adc_b2_ch14",
//    "adc_b3_ch0", "adc_b3_ch1", "adc_b3_ch2", "adc_b3_ch3", "adc_b3_ch4", "adc_b3_ch5",
//    "adc_b3_ch6", "adc_b3_ch7", "adc_b3_ch8", "adc_b3_ch9", "adc_b3_ch10", "adc_b3_ch11",
//    "adc_b3_ch12", "adc_b3_ch13", "adc_b3_ch14", "adc_b3_ch15",
//    "adc_b4_ch0", "adc_b4_ch1", "adc_b4_ch2", "adc_b4_ch3", "adc_b4_ch4", "adc_b4_ch5", "adc_b4_ch6", "adc_b4_ch7","adc_b4_ch8", "adc_b4_ch9", "adc_b4_ch10", "adc_b4_ch11",
//};
//
//std::vector<std::string> bot_pmt_channels = {
//    "adc_b1_ch1", "adc_b1_ch2", "adc_b1_ch3", "adc_b1_ch4",
//    "adc_b1_ch5", "adc_b1_ch6", "adc_b1_ch7", "adc_b1_ch8",
//    "adc_b1_ch9", "adc_b1_ch10", "adc_b1_ch11", "adc_b1_ch12",
//    "adc_b1_ch13", "adc_b1_ch14", "adc_b1_ch15", "adc_b2_ch0",
//    "adc_b2_ch1", "adc_b2_ch2", "adc_b2_ch3", "adc_b2_ch4",
//    "adc_b2_ch5", "adc_b2_ch6", "adc_b2_ch7", "adc_b2_ch8",
//    "adc_b2_ch9", "adc_b2_ch10", "adc_b2_ch11", "adc_b2_ch12",
//    "adc_b2_ch13", "adc_b2_ch14"
//};
//
//void setSPEResult(std::string inputSPECalibrationPath) {
//    std::ifstream file(inputSPECalibrationPath);
//    if (!file.is_open()) {
//        std::cerr << "Your spe_fit_results_file cannot be loaded properly!" << std::endl;
//        std::exit(EXIT_FAILURE);
//    }
//
//    std::string line, header;
//    std::getline(file, header); // 첫 줄 읽기 (헤더)
//
//    while (std::getline(file, line)) {
//        std::istringstream iss(line);
//        std::string ch_id, ch_name, pmt, spe_mean_str;
//        float spe_mean_value;
//
//        // ','로 나누기 (CSV 형식)
//        std::getline(iss, ch_id, ',');
//        std::getline(iss, ch_name, ',');
//        std::getline(iss, pmt, ',');
//        std::getline(iss, spe_mean_str, ',');
//        //std::cout << spe_mean_str << std::endl;
//
//        try {
//            spe_mean_value = std::stof(spe_mean_str);
//        } catch (const std::exception& e) {
//            std::cerr << "Error parsing float value: " << e.what() << std::endl;
//            std::exit(EXIT_FAILURE);
//        }
//
//        // spe_mean과 spe_fit_results에 값 저장
//        spe_mean[ch_name] = spe_mean_value;
//        spe_fit_results[ch_name]["spe_mean"] = spe_mean_value;
//    }
//
//    file.close();
//}
//
//
//void setBranches(TTree* inputTree, std::vector<std::string>& pmts, std::vector<UShort_t*>& dataStorage) {
//    for (size_t i = 0; i < pmts.size(); ++i) {
//        UShort_t* tempArray = new UShort_t[MAX_SAMPLE_SIZE]();
//        dataStorage.push_back(tempArray); // 포인터를 저장
//        inputTree->SetBranchAddress(pmts[i].c_str(), tempArray); // 브랜치와 연결
//    }
//}
//
//void processFile(std::string inputFileName) {
//    TFile* inputFile = new TFile(inputFileName.c_str());
//    TTree* inputTree = (TTree*)inputFile->Get("daq");
//
//    std::cout << "entry number: " << inputTree->GetEntries() << std::endl;
//
//    UInt_t event_id; inputTree->SetBranchAddress("event_id", &event_id);
//
//    std::vector<UShort_t*> dataStorage;
//    std::vector<UShort_t*> trigDataStorage;
//
//    setBranches(inputTree, pmts, dataStorage);
//    setBranches(inputTree, triggers, trigDataStorage);
//
//
//    int trigger_type = -1;
//
//    for (int ievt = 0; ievt < event_number_to_process; ++ievt) {
//        inputTree->GetEntry(ievt);
//        if (event_id % 1000 == 0)
//            std::cout << "event_id: " << event_id << std::endl;
//
//        Waveform topPaddleWaveform(trigDataStorage[0]);
//        Waveform alphaWaveform(trigDataStorage[1]);
//        Waveform majorityWaveform(trigDataStorage[2]);
//        if (Waveform::hasValueLessThan(topPaddleWaveform.getSamples(), 3000))
//            trigger_type = 0;
//        if (Waveform::hasValueLessThan(alphaWaveform.getSamples(), 3000))
//            trigger_type = 1;
//        if (Waveform::hasValueLessThan(majorityWaveform.getSamples(), 3000))
//            trigger_type = 2;
//
//        if (TRRIGER != trigger_type) {
//            std::cout << "skip " << event_id << " not desired trigger" << std::endl;
//            continue;
//        }
//
//        double totalPE = 0;
//
//        //std::cout << "event_id: " << event_id << std::endl;
//        for (size_t i = 0; i < pmts.size(); ++i) {
//            std::string ch_name = pmts[i];
//
//            Waveform wf(dataStorage[i]);
//
//            wf.subtractFlatBaseline(0, 100);
//
//            wf.setAmpPE(spe_mean[ch_name]);
//
//            wf.correctDaisyChainTrgDelay(ch_name);
//
//            double pe_value = wf.getPE(INTEGRATE_RANGE_START, INTEGRATE_RANGE_END);
//            pe[ch_name].push_back(pe_value);
//        }
//        //std::cout << "total PE: " << totalPE << std::endl;
//    }
//}
//
//void saveBinary(const std::string& filename, const std::vector<double>& data) {
//    std::ofstream file(filename, std::ios::binary);
//    std::cout << "saving " << filename << std::endl;
//    if (!file.is_open()) {
//        std::cerr << "Error opening file for writing!" << std::endl;
//        return;
//    }
//
//    // 데이터 크기 저장
//    size_t size = data.size();
//    file.write(reinterpret_cast<const char*>(&size), sizeof(size));
//
//    // 데이터 저장
//    file.write(reinterpret_cast<const char*>(data.data()), size * sizeof(double));
//    file.close();
//}
//
//// 파일 이름에서 확장자를 제거하고 식별자를 추출
//std::string extractFileID(const std::string& filePath) {
//    // 파일 이름 추출
//    size_t lastSlash = filePath.find_last_of('/');
//    std::string fileName = (lastSlash == std::string::npos) ? filePath : filePath.substr(lastSlash + 1);
//
//    // 확장자 제거 (.root)
//    size_t lastDot = fileName.find_last_of('.');
//    if (lastDot != std::string::npos) {
//        fileName = fileName.substr(0, lastDot);
//    }
//
//    return fileName; // 파일 이름 반환
//}
//
//int main(int argc, char **argv) {
//    auto start = std::chrono::high_resolution_clock::now();
//    std::string inputFileName;
//    std::string outputFilePath;
//    std::string inputSPECalibrationPath = "/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/yaml/bnl1t_spe_fit_results_241111.csv";
//
//    int c;
//    while (1) {
//        static struct option long_options[] = {
//            {"inputFileName", required_argument, 0, 'i'},
//            {"outputFilePath", required_argument, 0, 'o'},
//            {"eventNumber", required_argument, 0, 'n'},
//            {"trigger type", required_argument, 0, 't'},
//            {"inputSPECalibrationPath", required_argument, 0, 'c'},
//            {0, 0, 0, 0}};
//
//        int option_index = 0;
//        c = getopt_long(argc, argv, "i:o:n:t:c:", long_options, &option_index);
//
//        if (c == -1)
//            break;
//
//        switch (c) {
//        case 'i':
//            inputFileName = optarg;
//            break;
//        case 'o':
//            outputFilePath = optarg;
//            break;
//        case 'n':
//            event_number_to_process = std::stoi(optarg);
//            break;
//        case 't':
//            TRRIGER = std::stoi(optarg);
//            break;
//        case 'c':
//            inputSPECalibrationPath = optarg;
//            break;
//        default:
//            std::cerr << "Usage: " << argv[0]
//                      << " -i inputFileName -o outputFilePath -n number of events to process -t trigger type(0: tp, 1: alpha, 2: majority) -c inputSPECalibrationPath (optional)" << std::endl;
//            exit(1);
//        }
//    }
//
//    if (inputFileName.empty()) {
//        std::cerr << "Usage: " << argv[0]
//                      << " -i inputFileName -o outputFilePath -n number of events to process -t trigger type(0: tp, 1: alpha, 2: majority) -c inputSPECalibrationPath (optional)" << std::endl;
//        return 1;
//    }
//    if (TRRIGER == 0) {
//        INTEGRATE_RANGE_START = 380;
//        INTEGRATE_RANGE_END = 450;
//    } else if (TRRIGER == 1) {
//        INTEGRATE_RANGE_START = 320;
//        INTEGRATE_RANGE_END = 390;
//    } else if (TRRIGER ==2) {
//        INTEGRATE_RANGE_START = 160;
//        INTEGRATE_RANGE_END = 230;
//    }
//
//    setSPEResult(inputSPECalibrationPath);
//
//    processFile(inputFileName);
//
//    std::string fileID = extractFileID(inputFileName);
//    for (auto& pmt : pmts) {
//        saveBinary(outputFilePath + "npe_channel_" + fileID + "_" + pmt, pe[pmt]);
//    }
//        auto end = std::chrono::high_resolution_clock::now();
//
//    std::chrono::duration<double> elapsed = end - start;
//    std::cout << "Execution time: " << elapsed.count() << " seconds" << std::endl;
//
//    return 0;
//}
#include "AppConfig.h"
#include "DataProcessor.h"
#include <chrono>
#include <iostream>

int main(int argc, char** argv) {
    auto start = std::chrono::high_resolution_clock::now();

    // 커맨드라인 인자 파싱 후 설정 객체 생성
    AppConfig config = AppConfig::parseArgs(argc, argv);

    // DataProcessor 객체 생성 및 실행
    DataProcessor processor(config);
    processor.run();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Execution time: " << elapsed.count() << " seconds" << std::endl;
    return 0;
}
