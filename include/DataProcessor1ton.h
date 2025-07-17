#ifndef DATAPROCESSOR1TON_H
#define DATAPROCESSOR1TON_H

#include "AppConfig.h"
#include <string>
#include <vector>
#include <map>
#include <memory>
#include "Waveform.h"
#include <TTree.h>

struct PMTData1ton {
    std::vector<std::unique_ptr<UShort_t[]>> pmts;
    std::vector<std::unique_ptr<UShort_t[]>> triggers;
    std::vector<std::unique_ptr<UShort_t[]>> alphaPMT;
    std::vector<std::unique_ptr<UShort_t[]>> topPaddles;
    std::vector<std::unique_ptr<UShort_t[]>> botPaddle;
};

struct EventSelectionResult1ton {
    bool isSelected = false;
    bool isCrossingMuon = false;
};

class DataProcessor1ton {
public:
    explicit DataProcessor1ton(const AppConfig& config);
    void dailyCheck();

    void setPMTs(std::vector<std::string> pmts) {this->pmts_ = pmts;};
    void setTriggers(std::vector<std::string> triggers) {this->triggers_ = triggers;};
    void run();

    void DrawWaveforms(int event_id, std::vector<Waveform> waveforms);
private:
    AppConfig config_;
    std::vector<std::string> pmts_; //data analysis
    std::vector<std::string> pmtsAll_; //daily check
    std::vector<std::string> triggers_;
    std::map<std::string, std::vector<double>> pe_; 
    std::map<std::string, double> spe_mean_;
    std::map<std::string, std::map<std::string, double>> spe_fit_results_; // ch_name 기준으로 결과 저장

    void setSPEResult(const std::string& calibrationPath);
    void processFile();
    void setupBranches(TTree* tree, const std::vector<std::string>& branchNames,
                   std::vector<std::unique_ptr<UShort_t[]>>& storage,
                   int sampleSize);

    void saveRootOutput();
    void saveBinaryOutput();
    std::string extractFileID(const std::string& filePath);

    std::vector<Waveform> processWaveforms(const std::vector<std::unique_ptr<UShort_t[]>>& pmtRawData);
    EventSelectionResult1ton selectEvent(const PMTData1ton& pmtData, int event_id);
    PMTData1ton setupAllBranches(TTree* tree);
    int findPeakTime(const std::vector<Waveform>& waveforms);
    std::vector<double> calculatePEs(const std::vector<Waveform>& waveforms, int peakTime);
};

#endif // DATAPROCESSOR1TON_H
