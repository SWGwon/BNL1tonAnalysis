#include "DataProcessor1ton.h"
#include <numeric>

#include "Waveform.h"
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

namespace {
    constexpr int MAX_SAMPLE_SIZE = 3000;
    constexpr int ADC_SATURATION_VALUE = 15400;
    constexpr int BASELINE_END_BIN = 100;
    constexpr int PE_INTEGRATION_PRE_BINS = 10;
    constexpr int PE_INTEGRATION_POST_BINS = 20;
    constexpr double ALPHA_PEAK_THRESHOLD_MV = 16.0 / 50.0;
}

DataProcessor1ton::DataProcessor1ton(const AppConfig &config) : config_(config) {
    for (const auto &pmt : pmts_) {
        pe_[pmt] = std::vector<double>();
    }
}

void DataProcessor1ton::run() {
    setSPEResult(config_.inputSPECalibrationPath);
    processFile();
    saveRootOutput();
}

void DataProcessor1ton::processFile() {
    TFile *inputFile = new TFile(config_.inputFileName.c_str());
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open input file " << config_.inputFileName << std::endl;
        return;
    }
    TTree *inputTree = (TTree *)inputFile->Get("daq");
    if (!inputTree) {
        std::cerr << "Error: TTree 'daq' not found in file " << config_.inputFileName << std::endl;
        return;
    }

    PMTData1ton pmtData = setupAllBranches(inputTree);
    UInt_t event_id;
    inputTree->SetBranchAddress("event_id", &event_id);

    const int event_number = (config_.eventNumber > 0 && config_.eventNumber < inputTree-> GetEntries())
        ? config_.eventNumber
        : inputTree->GetEntries();


    int start = 0, end = 0;
    if (config_.triggerType == 0 || config_.triggerType == 3) {
        start = 380;
        end = 450;
    } else if (config_.triggerType == 1) {
        start = 320;
        end = 360;
    } else if (config_.triggerType == 2) {
        start = 160;
        end = 230;
    } else if (config_.triggerType == 4) {
        start = 150;
        end = 210;
    }

    for (int ievt = 0; ievt < event_number; ++ievt) {
        inputTree->GetEntry(ievt);
        if (ievt % 1000 == 0)
            std::cout << event_id << std::endl;

        EventSelectionResult1ton selection = selectEvent(pmtData, event_id);
        if (!selection.isSelected) {
            continue;
        }

        auto processedWaveforms = processWaveforms(pmtData.pmts);
        if (processedWaveforms.empty()) {
            std::cout << "Skipping event " << event_id << " due to invalid waveform size." << std::endl;
            continue;
        }

        const int peakTime = findPeakTime(processedWaveforms);
        if (peakTime < 0) continue;

        if (peakTime < start || peakTime > start + 40) {
            continue;
        }

        auto peValues = calculatePEs(processedWaveforms, peakTime);

        for (size_t i = 0; i < pmts_.size(); ++i) {
            pe_[pmts_[i]].push_back(peValues[i]);
        }
        double totalPE = std::accumulate(peValues.begin(), peValues.end(), 0.0);
        std::cout << "Event " << event_id << " processed. Peak time: " << peakTime << ", Total PE: " << totalPE << std::endl;
    }
    inputFile->Close();
}

void DataProcessor1ton::saveRootOutput() {
    // Extract the file ID from the input file name.
    std::string fileID = extractFileID(config_.inputFileName);

    TFile* outputFile = nullptr;
    // Create the output ROOT file.
    if (config_.triggerType == -1)
        outputFile = new TFile((config_.outputFilePath + "/output_all_event_" + fileID + ".root").c_str(), "RECREATE");
    if (config_.triggerType == 0)
        outputFile = new TFile((config_.outputFilePath + "/output_top_paddle_" + fileID + ".root").c_str(), "RECREATE");
    if (config_.triggerType == 1)
        outputFile = new TFile((config_.outputFilePath + "/output_alpha_" + fileID + ".root").c_str(), "RECREATE");
    if (config_.triggerType == 2)
        outputFile = new TFile((config_.outputFilePath + "/output_majority_" + fileID + ".root").c_str(), "RECREATE");
    if (config_.triggerType == 3 || config_.triggerType == 4)
        outputFile = new TFile((config_.outputFilePath + "/output_crossing_muon_" + fileID + ".root").c_str(), "RECREATE");

    // Create a new TTree.
    TTree* outputTree = new TTree("tree", "tree");

    // Create a vector to hold the output PE values for each PMT.
    // Using a std::vector is preferred over a variable-length array.
    std::vector<double> outputPE(pmts_.size(), 0.0);
    // Create a branch for each PMT.
    int branchIndex = 0;
    for (const auto &pmt : pmts_) {
        // The branch name is the PMT name and the branch variable is the corresponding element in outputPE.
        outputTree->Branch(pmt.c_str(), &outputPE[branchIndex]);
        branchIndex++;
    }
    double totalPE = 0;
    outputTree->Branch("totalPE", &totalPE);

    // Check that we have at least one PMT.
    if (pmts_.empty()) {
        std::cerr << "No PMTs available for output." << std::endl;
        outputFile->Close();
        return;
    }

    // Determine the number of events from the first PMT's vector.
    // Assumes that all PMTs have the same number of events.
    size_t numEvents = pe_[pmts_.front()].size();

    // Loop over all events.
    for (size_t ievt = 0; ievt < numEvents; ++ievt) {
        totalPE = 0.0;
        branchIndex = 0;
        // Loop over all PMTs to fill the output array.
        for (const auto &pmt : pmts_) {
            // Check that the current event index is within bounds for the PMT.
            if (ievt < pe_[pmt].size()) {
                outputPE[branchIndex] = pe_[pmt][ievt];
            } else {
                // In case of missing data, set a default value.
                outputPE[branchIndex] = 0.0;
            }
            totalPE += outputPE[branchIndex];
            branchIndex++;
        }
        // Fill the tree with the event's data.
        outputTree->Fill();
    }

    // Write the tree to the file and close it.
    outputTree->Write();
    outputFile->Close();
}

void DataProcessor1ton::setSPEResult(const std::string &calibrationPath) {
    std::ifstream file(calibrationPath);
    if (!file.is_open()) {
        std::cerr << "Your spe_fit_results_file cannot be loaded properly!"
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }

    std::string line, header;
    std::getline(file, header);

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string ch_id, ch_name, pmt, spe_mean_str;
        double spe_mean_value;

        std::getline(iss, ch_id, ',');
        std::getline(iss, ch_name, ',');
        std::getline(iss, pmt, ',');
        std::getline(iss, spe_mean_str, ',');

        try {
            spe_mean_value = std::stod(spe_mean_str);
        } catch (const std::exception &e) {
            std::cerr << "Error parsing float value: " << e.what() << std::endl;
            std::exit(EXIT_FAILURE);
        }

        spe_mean_[ch_name] = spe_mean_value;
        spe_fit_results_[ch_name]["spe_mean"] = spe_mean_value;
        std::cout << "set spe mean for " << ch_name << " : "
                  << spe_mean_[ch_name] << std::endl;
    }

    file.close();
}

std::string DataProcessor1ton::extractFileID(const std::string &filePath) {
    size_t lastSlash = filePath.find_last_of('/');
    std::string fileName = (lastSlash == std::string::npos)
                               ? filePath
                               : filePath.substr(lastSlash + 1);
    size_t lastDot = fileName.find_last_of('.');
    if (lastDot != std::string::npos) {
        fileName = fileName.substr(0, lastDot);
    }
    return fileName;
}

void DataProcessor1ton::setupBranches(TTree* tree, const std::vector<std::string>& branchNames,
                   std::vector<std::unique_ptr<UShort_t[]>>& storage,
                   int sampleSize) {
    for (const auto& name : branchNames) {
        auto tempArray = std::make_unique<UShort_t[]>(sampleSize);
        tree->SetBranchStatus(name.c_str(), 1);
        tree->SetBranchAddress(name.c_str(), tempArray.get());
        tree->AddBranchToCache(name.c_str());
        storage.push_back(std::move(tempArray));
    }
}

PMTData1ton DataProcessor1ton::setupAllBranches(TTree* tree) {
    PMTData1ton pmtData;

    setupBranches(tree, pmts_, pmtData.pmts, MAX_SAMPLE_SIZE);
    setupBranches(tree, triggers_, pmtData.triggers, MAX_SAMPLE_SIZE);
    setupBranches(tree, {"adc_b4_ch12"}, pmtData.alphaPMT, MAX_SAMPLE_SIZE);
    setupBranches(tree, {"adc_b4_ch13", "adc_b4_ch14"}, pmtData.topPaddles, MAX_SAMPLE_SIZE);
    setupBranches(tree, {"adc_b1_ch0"}, pmtData.botPaddle, MAX_SAMPLE_SIZE);

    return pmtData;
}

EventSelectionResult1ton DataProcessor1ton::selectEvent(const PMTData1ton& pmtData, int
event_id) {
    EventSelectionResult1ton result;

    if (config_.triggerType == 3 || config_.triggerType == 4) {
        Waveform botPaddle(pmtData.botPaddle[0].get());
        Waveform tp1(pmtData.topPaddles[0].get());
        Waveform tp2(pmtData.topPaddles[1].get());
        if (Waveform::hasValueLessThan(botPaddle.getSamples(), ADC_SATURATION_VALUE) &&
            Waveform::hasValueLessThan(tp1.getSamples(), ADC_SATURATION_VALUE) &&
            Waveform::hasValueLessThan(tp2.getSamples(), ADC_SATURATION_VALUE)) {
            result.isCrossingMuon = true;
        } else {
            return result;
        }
    }

    if (config_.triggerType == 1) {
        Waveform alpha(pmtData.alphaPMT[0].get());
        alpha.subtractFlatBaseline(0, BASELINE_END_BIN);
        if (!alpha.hasPeakAboveThreshold(ALPHA_PEAK_THRESHOLD_MV)) {
            return result;
        }
    }

    result.isSelected = true;
    return result;
}

std::vector<Waveform> DataProcessor1ton::processWaveforms(const std::vector<std::unique_ptr<UShort_t[]>>& pmtRawData) {
    std::vector<Waveform> processedWaveforms;
    processedWaveforms.reserve(pmts_.size());

    for (size_t i = 0; i < pmts_.size(); ++i) {
        Waveform wf(pmtRawData[i].get());
        std::string ch_name = pmts_[i];

        if (wf.getSamples().size() != config_.sampleSize) {
            std::cout << "configed sample size: " << config_.sampleSize << ", but " << ch_name << " sample size: " << wf.getSamples().size() << std::endl;
            return {};
        }

        wf.subtractFlatBaseline(0, BASELINE_END_BIN);

        // spe_mean 값이 없는 경우에 대한 예외 처리
        auto it = spe_mean_.find(pmts_[i]);
        if (it == spe_mean_.end()) {
            std::cerr << "Warning: SPE value not found for " << pmts_[i] << ". Skipping PE conversion." << std::endl;
            // PE 변환 없이 추가하거나, 이벤트를 건너뛸 수 있음
        } else {
            wf.setAmpPE(it->second);
        }

        wf.correctDaisyChainTrgDelay(pmts_[i]);
        processedWaveforms.push_back(std::move(wf));
    }
    return processedWaveforms;
}

int DataProcessor1ton::findPeakTime(const std::vector<Waveform>& waveforms) {
    if (waveforms.empty() || waveforms[0].getSamples().empty()) {
        return -1;
    }

    std::vector<double> summedWaveform(waveforms[0].getSamples().size(), 0.0);
    for (const auto& wf : waveforms) {
        // getAmpPE 대신 mV 단위 파형을 합산 (PE는 spe_mean에 따라 스케일이 다르므로)
        const auto& samples = wf.getAmpMV();
        for (size_t i = 0; i < samples.size(); ++i) {
            summedWaveform[i] += samples[i];
        }
    }

    auto maxIt = std::max_element(summedWaveform.begin(), summedWaveform.end());
    return std::distance(summedWaveform.begin(), maxIt);
}

std::vector<double> DataProcessor1ton::calculatePEs(const std::vector<Waveform>& waveforms, int peakTime) {
    std::vector<double> peValues;
    peValues.reserve(waveforms.size());

    int lowlim = peakTime - PE_INTEGRATION_PRE_BINS;
    if (lowlim < 0) lowlim = 0;

    for (size_t i = 0; i < waveforms.size(); ++i) {
        const auto& wf = waveforms[i];
        const std::string& ch_name = pmts_[i];

        int upperlim = peakTime + PE_INTEGRATION_POST_BINS;
        if (upperlim >= wf.getSamples().size()) {
            upperlim = wf.getSamples().size() - 1;
        }

        double pe_value = 0;
        // 특정 트리거 타입에 따라 다른 적분 구간을 사용하는 로직
        if (config_.triggerType == 1) { // Alpha event
            pe_value = wf.getPE(310, 350, spe_mean_.at(ch_name));
        } else {
            pe_value = wf.getPE(lowlim, upperlim, spe_mean_.at(ch_name));
        }
        peValues.push_back(pe_value);
    }
    return peValues;
}
