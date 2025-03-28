#ifndef DATAPROCESSOR_H
#define DATAPROCESSOR_H

#include "AppConfig.h"
#include <string>
#include <vector>
#include <map>
#include "Waveform.h"

class DataProcessor {
public:
    explicit DataProcessor(const AppConfig& config);
    void dailyCheck();
    void dailyCheck30t();
    void run();
    void saveBinaryOutput();
    void saveRootOutput();
    void setPMTs(std::vector<std::string> pmts) {this->pmts_ = pmts;};
    void setTriggers(std::vector<std::string> triggers) {this->triggers_ = triggers;};
    void DrawWaveforms(int event_id, std::vector<Waveform> waveforms);
private:
    AppConfig config_;
    std::vector<std::string> pmts_; //data analysis
    std::vector<std::string> pmtsAll_; //daily check
    std::vector<std::string> pmts30t; //daily check 30t
    std::vector<std::string> triggers_;
    std::vector<std::string> triggers_30t_;
    std::map<std::string, std::vector<double>> pe_; 
    std::map<std::string, double> spe_mean_;
    std::map<std::string, std::map<std::string, double>> spe_fit_results_; // ch_name 기준으로 결과 저장

    void setSPEResult(const std::string& calibrationPath);
    void processFile();

    std::string extractFileID(const std::string& filePath);
};

#endif // DATAPROCESSOR_H
