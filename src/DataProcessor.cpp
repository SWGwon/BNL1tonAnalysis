#include "DataProcessor.h"
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

DataProcessor::DataProcessor(const AppConfig &config) : config_(config) {
    pmts_ = {"adc_b1_ch1",  "adc_b1_ch2",  "adc_b1_ch3",  "adc_b1_ch4",
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
             "adc_b4_ch10", "adc_b4_ch11"};
    pmtsAll_ = {"adc_b1_ch0",  "adc_b1_ch1",  "adc_b1_ch2",  "adc_b1_ch3",
                "adc_b1_ch4",  "adc_b1_ch5",  "adc_b1_ch6",  "adc_b1_ch7",
                "adc_b1_ch8",  "adc_b1_ch9",  "adc_b1_ch10", "adc_b1_ch11",
                "adc_b1_ch12", "adc_b1_ch13", "adc_b1_ch14", "adc_b1_ch15",
                "adc_b2_ch0",  "adc_b2_ch1",  "adc_b2_ch2",  "adc_b2_ch3",
                "adc_b2_ch4",  "adc_b2_ch5",  "adc_b2_ch6",  "adc_b2_ch7",
                "adc_b2_ch8",  "adc_b2_ch9",  "adc_b2_ch10", "adc_b2_ch11",
                "adc_b2_ch12", "adc_b2_ch13", "adc_b2_ch14", "adc_b2_ch15",
                "adc_b3_ch0",  "adc_b3_ch1",  "adc_b3_ch2",  "adc_b3_ch3",
                "adc_b3_ch4",  "adc_b3_ch5",  "adc_b3_ch6",  "adc_b3_ch7",
                "adc_b3_ch8",  "adc_b3_ch9",  "adc_b3_ch10", "adc_b3_ch11",
                "adc_b3_ch12", "adc_b3_ch13", "adc_b3_ch14", "adc_b3_ch15",
                "adc_b4_ch0",  "adc_b4_ch1",  "adc_b4_ch2",  "adc_b4_ch3",
                "adc_b4_ch4",  "adc_b4_ch5",  "adc_b4_ch6",  "adc_b4_ch7",
                "adc_b4_ch8",  "adc_b4_ch9",  "adc_b4_ch10", "adc_b4_ch11",
                "adc_b4_ch12", "adc_b4_ch13", "adc_b4_ch14", "adc_b4_ch15",
                "adc_b5_ch0",  "adc_b5_ch1",  "adc_b5_ch2",  "adc_b5_ch3",
                "adc_b5_ch4",  "adc_b5_ch5",  "adc_b5_ch6",  "adc_b5_ch7",
                "adc_b5_ch8",  "adc_b5_ch9",  "adc_b5_ch10", "adc_b5_ch11",
                "adc_b5_ch12", "adc_b5_ch13", "adc_b5_ch14", "adc_b5_ch15",
                "adc_b5_ch16", "adc_b5_ch17", "adc_b5_ch18", "adc_b5_ch19",
                "adc_b5_ch20", "adc_b5_ch21", "adc_b5_ch22", "adc_b5_ch23",
                "adc_b5_ch24", "adc_b5_ch25", "adc_b5_ch26", "adc_b5_ch27",
                "adc_b5_ch28", "adc_b5_ch29", "adc_b5_ch30", "adc_b5_ch31",
                "adc_b5_ch32"};
    pmts30t = {"adc_b1_ch0",  "adc_b1_ch1",  "adc_b1_ch2",  "adc_b1_ch3",
               "adc_b1_ch4",  "adc_b1_ch5",  "adc_b1_ch6",  "adc_b1_ch7",
               "adc_b1_ch8",  "adc_b1_ch9",  "adc_b1_ch10", "adc_b1_ch11",
               "adc_b1_ch12", "adc_b1_ch13", "adc_b1_ch14", "adc_b1_ch15",
               "adc_b2_ch0",  "adc_b2_ch1",  "adc_b2_ch2",  "adc_b2_ch3",
               "adc_b2_ch4",  "adc_b2_ch5",  "adc_b2_ch6",  "adc_b2_ch7",
               "adc_b2_ch8",  "adc_b2_ch9",  "adc_b2_ch10", "adc_b2_ch11",
               "adc_b2_ch12", "adc_b2_ch13", "adc_b2_ch14", "adc_b2_ch15",
               "adc_b3_ch0",  "adc_b3_ch1",  "adc_b3_ch2",  "adc_b3_ch3"};
    triggers_ = {"adc_b5_ch33", "adc_b5_ch34", "adc_b5_ch35"};

    for (const auto &pmt : pmts_) {
        pe_[pmt] = std::vector<double>();
    }
}

void DataProcessor::setSPEResult(const std::string &calibrationPath) {
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

std::string DataProcessor::extractFileID(const std::string &filePath) {
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

void DataProcessor::processFile() {
    const int MAX_SAMPLE_SIZE = 3000;
    TFile *inputFile = new TFile(config_.inputFileName.c_str());
    TTree *inputTree = (TTree *)inputFile->Get("daq");
    UInt_t event_id;
    inputTree->SetBranchAddress("event_id", &event_id);

    std::vector<UShort_t *> dataStorage;
    std::vector<UShort_t *> trigDataStorage;

    for (auto &pmt : pmts_) {
        UShort_t *tempArray = new UShort_t[MAX_SAMPLE_SIZE]();
        dataStorage.push_back(tempArray);
        inputTree->SetBranchAddress(pmt.c_str(), tempArray);
    }

    for (auto &trig : triggers_) {
        UShort_t *tempArray = new UShort_t[MAX_SAMPLE_SIZE]();
        trigDataStorage.push_back(tempArray);
        inputTree->SetBranchAddress(trig.c_str(), tempArray);
    }

    int trigger_type = -1;
    int event_number = 0;
    if (config_.eventNumber > inputTree->GetEntries())
        event_number = inputTree->GetEntries();
    else
        event_number = config_.eventNumber;
    for (int ievt = 0; ievt < event_number; ++ievt) {
        inputTree->GetEntry(ievt);
        // checking trigger
        Waveform topPaddleWaveform(trigDataStorage[0]);
        Waveform alphaWaveform(trigDataStorage[1]);
        Waveform majorityWaveform(trigDataStorage[2]);
        if (Waveform::hasValueLessThan(topPaddleWaveform.getSamples(), 3000))
            trigger_type = 0;
        if (Waveform::hasValueLessThan(alphaWaveform.getSamples(), 3000))
            trigger_type = 1;
        if (Waveform::hasValueLessThan(majorityWaveform.getSamples(), 3000))
            trigger_type = 2;

        if (config_.triggerType != trigger_type) {
            std::cout << "skip " << event_id << " not desired trigger"
                      << std::endl;
            if (trigger_type == 0)
                std::cout << "this event is top paddle triggered" << std::endl;
            if (trigger_type == 1)
                std::cout << "this event is alpha triggered" << std::endl;
            if (trigger_type == 2)
                std::cout << "this event is majority triggered" << std::endl;
            continue;
        }

        for (size_t i = 0; i < pmts_.size(); ++i) {
            std::string ch_name = pmts_[i];
            Waveform wf(dataStorage[i]);
            wf.subtractFlatBaseline(0, 100);
            wf.setAmpPE(spe_mean_[ch_name]);
            wf.correctDaisyChainTrgDelay(ch_name);

            int start = 0, end = 0;
            if (config_.triggerType == 0) {
                start = 380;
                end = 450;
            } else if (config_.triggerType == 1) {
                start = 320;
                end = 390;
            } else if (config_.triggerType == 2) {
                start = 160;
                end = 230;
            }
            double pe_value = wf.getPE(start, end);
            pe_[ch_name].push_back(pe_value);
        }
    }
    inputFile->Close();
}

void DataProcessor::saveOutput() {
    std::string fileID = extractFileID(config_.inputFileName);
    for (auto &pmt : pmts_) {
        std::string filename =
            config_.outputFilePath + "npe_channel_" + fileID + "_" + pmt;
        std::cout << "saving " << filename << std::endl;
        std::ofstream file(filename, std::ios::binary);
        if (!file.is_open()) {
            std::cerr << "Error opening file for writing: " << filename
                      << std::endl;
            continue;
        }
        size_t size = pe_[pmt].size();
        file.write(reinterpret_cast<const char *>(&size), sizeof(size));
        file.write(reinterpret_cast<const char *>(pe_[pmt].data()),
                   size * sizeof(double));
        file.close();
    }
}
void DataProcessor::saveRootOutput() {
    std::string fileID = extractFileID(config_.inputFileName);
}

void DataProcessor::run() {
    setSPEResult(config_.inputSPECalibrationPath);
    processFile();
    saveOutput();
}

void DataProcessor::dailyCheck() {
    // Initialize file, tree, and histograms
    std::string fileID = extractFileID(config_.inputFileName);
    const int MAX_SAMPLE_SIZE = 3000;
    TFile *inputFile = new TFile(config_.inputFileName.c_str());
    TTree *inputTree = (TTree *)inputFile->Get("daq");

    // Set branch for event_id
    UInt_t event_id;
    inputTree->SetBranchAddress("event_id", &event_id);

    // Allocate storage for each PMT channel waveform data and create histograms
    std::vector<UShort_t *> dataStorage;
    std::map<std::string, TH1D *> histPMTPE;
    TH1D *histMaxWaveformIndex = new TH1D(
        "histMaxWaveformIndex",
        "histMaxWaveformIndex;relative time x 2ns;counts", 70, 60, 500);
    TH2D *hist_top_fire =
        new TH2D("hist_top_fire", "top hodoscope", 10, 0, 10, 10, 0, 10);
    TH2D *hist_bot_fire =
        new TH2D("hist_bot_fire", "bot hodoscope", 10, 0, 10, 10, 0, 10);

    for (const auto &pmt : pmtsAll_) {
        UShort_t *tempArray =
            new UShort_t[MAX_SAMPLE_SIZE](); // zero-initialize the array
        dataStorage.push_back(tempArray);
        inputTree->SetBranchAddress(pmt.c_str(), tempArray);
        histPMTPE[pmt] =
            new TH1D(pmt.c_str(), (pmt + ";mV ns;counts").c_str(), 50, 0, 0);
    }

    // Utility lambda to update hodoscope fire counts
    auto updateFire = [](const std::vector<int> &vec1,
                         const std::vector<int> &vec2, int fire[10][10]) {
        if (vec1.empty() && !vec2.empty()) {
            // If vec1 is empty, use the loop index for the column
            for (size_t jj = 0; jj < vec2.size(); ++jj)
                fire[9][jj]++;
        } else if (vec2.empty() && !vec1.empty()) {
            // If vec2 is empty, use the loop index for the row
            for (size_t ii = 0; ii < vec1.size(); ++ii)
                fire[ii][9]++;
        } else if (!vec1.empty() && !vec2.empty()) {
            // If both vectors are non-empty, use the stored indices in the
            // vectors
            for (size_t i = 0; i < vec1.size(); ++i)
                for (size_t j = 0; j < vec2.size(); ++j)
                    fire[vec1[i]][vec2[j]]++;
        }
    };

    // Lambda to save a canvas with a given filename suffix
    auto saveCanvas = [&](const std::string &suffix, TCanvas *canvas) {
        std::string outName = config_.outputFilePath + fileID + suffix;
        canvas->SaveAs(outName.c_str());
    };

    // Event loop
    int event_number = 0;
    if (config_.eventNumber > inputTree->GetEntries())
        event_number = inputTree->GetEntries();
    else
        event_number = config_.eventNumber;
    for (int ievt = 0; ievt < event_number; ++ievt) {
        inputTree->GetEntry(ievt);
        if (ievt % 1000 == 0)
            std::cout << event_id << std::endl;

        // Vectors for hodoscope channels (to store indices where the threshold
        // is exceeded)
        std::vector<int> top1, top2, bot1, bot2;
        std::vector<double> summedWaveform;
        double schn[97] = {0};

        // Process waveform for each PMT channel
        for (size_t i = 0; i < pmtsAll_.size(); ++i) {
            const std::string &ch_name = pmtsAll_[i];
            Waveform wf(dataStorage[i]);
            wf.subtractFlatBaseline(0, 100);
            double spe = (spe_mean_.count(ch_name) > 0) ? spe_mean_[ch_name]
                                                        : 0.267534455;
            wf.setAmpPE(spe);
            wf.correctDaisyChainTrgDelay(ch_name);

            // Sum waveforms for channels of interest (channels 1 to 30)
            if (i == 1) {
                summedWaveform = wf.getAmpPE();
            } else if (i > 1 && i < 31) {
                std::vector<double> tempWf = wf.getAmpPE();
                std::transform(summedWaveform.begin(), summedWaveform.end(),
                               tempWf.begin(), summedWaveform.begin(),
                               std::plus<double>());
            }

            // Calculate PE value for the channel and fill its histogram
            double pe_value = wf.getPE(0, wf.getAmpPE().size() - 1) / 2.0;
            pe_[ch_name].push_back(pe_value);
            schn[i] = pe_value;
            histPMTPE[ch_name]->Fill(pe_value);
        }

        // Fill the histogram for the maximum waveform index based on the summed
        // waveform
        if (!summedWaveform.empty()) {
            auto maxIt =
                std::max_element(summedWaveform.begin(), summedWaveform.end());
            int maxIndex = std::distance(summedWaveform.begin(), maxIt);
            histMaxWaveformIndex->Fill(maxIndex);
        }

        // Update hodoscope vectors based on a threshold (> 0.5)
        // Note: The indices (64 + ih, etc.) should match the hodoscope channels
        // in your data
        for (int ih = 0; ih < 8; ++ih) {
            if (schn[64 + ih] > 0.5)
                top1.push_back(ih);
            if (schn[64 + 8 + ih] > 0.5)
                top2.push_back(ih);
            if (schn[64 + 16 + ih] > 0.5)
                bot1.push_back(ih);
            if (schn[64 + 24 + ih] > 0.5)
                bot2.push_back(ih);
        }

        // Local fire arrays for this event
        int top_fire[10][10] = {0};
        int bot_fire[10][10] = {0};
        // updateFire(top1, top2, top_fire);
        // updateFire(bot1, bot2, bot_fire);
        if (!top1.empty() || !top2.empty()) {
            if (top1.empty() && !top2.empty()) {
                for (size_t jj = 0; jj < top2.size(); ++jj) {
                    top_fire[9][jj]++;
                }
            }
            if (top2.empty() && !top1.empty()) {
                for (size_t ii = 0; ii < top1.size(); ++ii) {
                    top_fire[ii][9]++;
                }
            }
            if (!top1.empty() && !top2.empty()) {
                for (size_t i = 0; i < top1.size(); ++i) {
                    for (size_t j = 0; j < top2.size(); ++j) {
                        top_fire[top1[i]][top2[j]]++;
                    }
                }
            }
        }
        if (!bot1.empty() || !bot2.empty()) {
            if (bot1.empty() && !bot2.empty()) {
                for (size_t jj = 0; jj < bot2.size(); ++jj) {
                    bot_fire[9][jj]++;
                }
            }
            if (bot2.empty() && !bot1.empty()) {
                for (size_t ii = 0; ii < bot1.size(); ++ii) {
                    bot_fire[ii][9]++;
                }
            }
            if (!bot1.empty() && !bot2.empty()) {
                for (size_t i = 0; i < bot1.size(); ++i) {
                    for (size_t j = 0; j < bot2.size(); ++j) {
                        bot_fire[bot1[i]][bot2[j]]++;
                    }
                }
            }
        }
        for (int iii = 0; iii < 10; ++iii) {
            for (int jjj = 0; jjj < 10; ++jjj) {
                int temp = hist_top_fire->GetBinContent(iii + 1, jjj + 1);
                hist_top_fire->SetBinContent(iii + 1, jjj + 1,
                                             temp + top_fire[iii][jjj]);

                temp = hist_bot_fire->GetBinContent(iii + 1, jjj + 1);
                hist_bot_fire->SetBinContent(iii + 1, jjj + 1,
                                             temp + bot_fire[iii][jjj]);
            }
        }

        // Add the event's fire counts to the overall hodoscope histograms
        for (int irow = 0; irow < 10; ++irow) {
            for (int jcol = 0; jcol < 10; ++jcol) {
                int tempTopBin = hist_top_fire->GetBin(irow + 1, jcol + 1);
                hist_top_fire->AddBinContent(tempTopBin, top_fire[irow][jcol]);
                int tempBotBin = hist_bot_fire->GetBin(irow + 1, jcol + 1);
                hist_bot_fire->AddBinContent(tempBotBin, bot_fire[irow][jcol]);
                //hist_top_fire->AddBinContent(irow + 1, jcol + 1,
                                             //top_fire[irow][jcol]);
                //hist_bot_fire->AddBinContent(irow + 1, jcol + 1,
                                             //bot_fire[irow][jcol]);
            }
        }
    } // End of event loop


    TCanvas *c0 = new TCanvas("", "", 1280, 960);
    c0->SetLogy();
    histMaxWaveformIndex->Draw();
    saveCanvas("_peak.png", c0);

    // Save PMT histograms grouped by canvas layout
    struct CanvasSetup {
        int numPads;
        int rows;
        int cols;
        std::string suffix;
    };
    std::vector<CanvasSetup> canvasSetups = {
        {16, 4, 4, "_b1.png"}, {16, 4, 4, "_b2.png"}, {16, 4, 4, "_b3.png"},
        {16, 4, 4, "_b4.png"}, {8, 2, 4, "_b5.png"},  {8, 2, 4, "_b6.png"},
        {8, 2, 4, "_b7.png"},  {8, 2, 4, "_b8.png"}};

    int tempPMTCount = 0;
    for (const auto &setup : canvasSetups) {
        TCanvas *canvas = new TCanvas("", "", 1280, 960);
        canvas->Divide(setup.rows, setup.cols);
        for (int i = 0; i < setup.numPads; ++i) {
            canvas->cd(i + 1);
            gPad->SetLogy();
            const std::string &ch_name = pmtsAll_[tempPMTCount];
            histPMTPE[ch_name]->Draw();
            ++tempPMTCount;
        }
        saveCanvas(setup.suffix, canvas);
    }

    // Save hodoscope histograms
    TCanvas *c_top = new TCanvas("", "", 1280, 960);
    c_top->SetLogz();
    hist_top_fire->SetStats(0);
    hist_top_fire->Draw("colz");
    saveCanvas("_top_hodoscope.png", c_top);

    TCanvas *c_bot = new TCanvas("", "", 1280, 960);
    c_bot->SetLogz();
    hist_bot_fire->SetStats(0);
    hist_bot_fire->Draw("colz");
    saveCanvas("_bot_hodoscope.png", c_bot);

    inputFile->Close();
}

void DataProcessor::dailyCheck30t() {
    std::string fileID = extractFileID(config_.inputFileName);
    const int MAX_SAMPLE_SIZE = 3000;
    TFile* inputFile = new TFile(config_.inputFileName.c_str());
    TTree* inputTree = dynamic_cast<TTree*>(inputFile->Get("daq"));
    if (!inputTree) {
        std::cerr << "Error: TTree 'daq' not found in file " << config_.inputFileName << std::endl;
        exit(EXIT_FAILURE);
    }
    UInt_t event_id;
    inputTree->SetBranchAddress("event_id", &event_id);

    std::vector<UShort_t*> dataStorage;
    std::vector<UShort_t*> trigDataStorage;  // (Unused in this function)

    std::map<std::string, TH1D*> histPMTPE;
    TH1D* histMaxWaveformIndex = new TH1D("histMaxWaveformIndex",
                                          "histMaxWaveformIndex;relative time x 2ns;counts",
                                          70, 60, 500);
    TH1D* histBotMV = new TH1D("histBotMV", "histBotMV;mV ns;counts", 50, 0, 0);
    TH1D* histSideMV = new TH1D("histSideMV", "histSideMV;mV ns;counts", 50, 0, 0);

    // Set branch addresses for each PMT channel and create corresponding histograms
    for (auto& pmt : pmts30t) {
        UShort_t* tempArray = new UShort_t[MAX_SAMPLE_SIZE]();  // zero-initialized array
        dataStorage.push_back(tempArray);
        inputTree->SetBranchAddress(pmt.c_str(), tempArray);
        histPMTPE[pmt] = new TH1D(pmt.c_str(), (pmt + ";mV ns;counts").c_str(), 50, 0, 0);
    }

    int event_number = 0;
    if (config_.eventNumber > inputTree->GetEntries())
        event_number = inputTree->GetEntries();
    else
        event_number = config_.eventNumber;
    for (int ievt = 0; ievt < event_number; ++ievt) {
        inputTree->GetEntry(ievt);
        if (ievt % 1000 == 0)
            std::cout << event_id << std::endl;

        int tempIndex = 0;
        double tempBotMV = 0;
        double tempSideMV = 0;
        std::vector<double> summedWaveform;

        // Process waveform for each PMT channel
        for (size_t i = 0; i < pmts30t.size(); ++i) {
            std::string ch_name = pmts30t[i];
            Waveform wf(dataStorage[i]);
            wf.subtractFlatBaseline(0, 100);
            wf.setAmpPE(1.0, 1.0);
            wf.correctDaisyChainTrgDelay(ch_name);

            // Sum waveforms: use the first channel as the base and add subsequent channels
            if (i == 0) {
                summedWaveform = wf.getAmpPE();
            } else {
                std::vector<double> tempWf = wf.getAmpPE();
                std::transform(summedWaveform.begin(), summedWaveform.end(),
                               tempWf.begin(),
                               summedWaveform.begin(), // store result back in summedWaveform
                               std::plus<double>());
            }

            // Set integration range (here using entire sample length)
            int start = 0;
            int end = MAX_SAMPLE_SIZE;
            double pe_value = wf.getPE(start, wf.getAmpPE().size() - 1);
            pe_[ch_name].push_back(pe_value / 2.0);
            histPMTPE[ch_name]->Fill(pe_value / 2.0);

            // Accumulate integrated values for bottom and side channels
            if (tempIndex < 12) {
                tempBotMV += pe_value / 2.0;
            } else {
                tempSideMV += pe_value / 2.0;
            }
            tempIndex++;
        } // end PMT channel loop

        // Fill histograms for bottom and side integrated values
        histBotMV->Fill(tempBotMV);
        histSideMV->Fill(tempSideMV);

        // Determine the index of the maximum value in the summed waveform
        auto maxIt = std::max_element(summedWaveform.begin(), summedWaveform.end());
        int maxIndex = std::distance(summedWaveform.begin(), maxIt);
        histMaxWaveformIndex->Fill(maxIndex);
    } // end event loop

    // Generate output filename using output file path and extracted fileID
    std::string filename = config_.outputFilePath + fileID;

    int tempPMTCount = 0;

    // Plot histogram for maximum waveform index
    TCanvas* c0 = new TCanvas();
    c0->SetLogy();
    histMaxWaveformIndex->Draw();
    std::cout << "saving " << filename + "_peak.pdf" << std::endl;
    c0->SaveAs((filename + "_peak.pdf").c_str());

    // Lambda to draw a set of PMT histograms on a divided canvas
    auto drawCanvas = [&](const std::string& canvasSuffix, int nRows, int nCols, int numHistograms) {
        TCanvas* c = new TCanvas();
        c->Divide(nCols, nRows);
        for (int i = 0; i < numHistograms; ++i) {
            c->cd(i + 1);
            c->cd(i + 1)->SetLogy();
            std::string ch_name = pmts30t[tempPMTCount];
            histPMTPE[ch_name]->Draw();
            tempPMTCount++;
        }
        std::cout << "saving " << filename + "_" + canvasSuffix + ".pdf" << std::endl;
        c->SaveAs((filename + "_" + canvasSuffix + ".pdf").c_str());
    };

    // Draw PMT histogram canvases
    drawCanvas("b1", 4, 4, 16);
    drawCanvas("b2", 4, 4, 16);
    drawCanvas("b3", 2, 2, 4);

    // Draw canvases for bottom and side integrated values
    TCanvas* c4 = new TCanvas();
    c4->SetLogy();
    histBotMV->Draw();
    std::cout << "saving " << filename + "_BotMV.pdf" << std::endl;
    c4->SaveAs((filename + "_BotMV.pdf").c_str());

    TCanvas* c5 = new TCanvas();
    c5->SetLogy();
    histSideMV->Draw();
    std::cout << "saving " << filename + "_SideMV.pdf" << std::endl;
    c5->SaveAs((filename + "_SideMV.pdf").c_str());

    inputFile->Close();
}
