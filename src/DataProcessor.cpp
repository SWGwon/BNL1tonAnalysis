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
    pmts30t = {"adc_b1_ch0",  "adc_b1_ch1",  "adc_b1_ch2",  "adc_b1_ch3",
               "adc_b1_ch4",  "adc_b1_ch5",  "adc_b1_ch6",  "adc_b1_ch7",
               "adc_b1_ch8",  "adc_b1_ch9",  "adc_b1_ch10", "adc_b1_ch11",
               "adc_b1_ch12", "adc_b1_ch13", "adc_b1_ch14", "adc_b1_ch15",
               "adc_b2_ch0",  "adc_b2_ch1",  "adc_b2_ch2",  "adc_b2_ch3",
               "adc_b2_ch4",  "adc_b2_ch5",  "adc_b2_ch6",  "adc_b2_ch7",
               "adc_b2_ch8",  "adc_b2_ch9",  "adc_b2_ch10", "adc_b2_ch11",
               "adc_b2_ch12", "adc_b2_ch13", "adc_b2_ch14", "adc_b2_ch15",
               "adc_b3_ch0",  "adc_b3_ch1",  "adc_b3_ch2",  "adc_b3_ch3"};
    triggers_30t_ = {"adc_b4_ch13", "adc_b4_ch22", "adc_b4_ch23"};

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
    std::vector<UShort_t *> botPaddleDataStorage;

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

    UShort_t *tempArray = new UShort_t[MAX_SAMPLE_SIZE]();
    botPaddleDataStorage.push_back(tempArray);
    inputTree->SetBranchAddress("adc_b1_ch0", tempArray);


    //unsigned short adc_b1_ch0[2000] = {};
    //inputTree->SetBranchAddress("adc_b1_ch0", adc_b1_ch0);

    TCanvas* c = new TCanvas();
    TH1D* hist_area_bp = new TH1D("hist_area_bp", "hist_area_bp", 100,0,100);


    int trigger_type = -1;
    int event_number = 0;
    if (config_.eventNumber > inputTree->GetEntries())
        event_number = inputTree->GetEntries();
    else
        event_number = config_.eventNumber;
    for (int ievt = 0; ievt < event_number; ++ievt) {
        inputTree->GetEntry(ievt);
        if (ievt % 1000 == 0)
            std::cout << event_id << std::endl;

        // checking trigger
        if (config_.triggerType != -1) {
            Waveform topPaddleWaveform(trigDataStorage[0]);
            Waveform alphaWaveform(trigDataStorage[1]);
            Waveform majorityWaveform(trigDataStorage[2]);
            if (Waveform::hasValueLessThan(topPaddleWaveform.getSamples(), 3000))
                trigger_type = 0;
            if (Waveform::hasValueLessThan(alphaWaveform.getSamples(), 3000))
                trigger_type = 1;
            if (Waveform::hasValueLessThan(majorityWaveform.getSamples(), 3000))
                trigger_type = 2;



            if (config_.triggerType == 3) {
                //if (trigger_type != 2)
                    //continue;
                //else {
                    double area_bp1 = 0;
                    const double ADC_TO_MV = 2000.0 / (std::pow(2, 14) - 1);
                    Waveform botPaddleWaveform(botPaddleDataStorage[0]);
                    //botPaddleWaveform.subtractFlatBaseline(0,100);
                    //botPaddleWaveform.setAmpPE(1.0);
                    //area_bp1 = botPaddleWaveform.getPE(0, botPaddleWaveform.getSamples().size() - 1);
                    //hist_area_bp->Fill(area_bp1);
                    if (Waveform::hasValueLessThan(botPaddleWaveform.getSamples(), 15400)) {
                        TCanvas c;
                        c.Divide(2,1);
                        TGraph* gr = botPaddleWaveform.drawMVAsGraph("");
                        TGraph* gr2 = topPaddleWaveform.drawMVAsGraph("");
                        //TGraph* gr = new TGraph(2000, x.data(), y.data());
                        c.cd(1);
                        gr->Draw();
                        c.cd(2);
                        gr2->Draw();

                        c.SaveAs(Form("%d.pdf",ievt));
                    } else {
                        continue;
                    }
                //}
            } else if (config_.triggerType != trigger_type) {
                //std::cout << "skip " << event_id << " not desired trigger"
                //    << std::endl;
                //if (trigger_type == 0)
                //    std::cout << "this event is top paddle triggered" << std::endl;
                //if (trigger_type == 1)
                //    std::cout << "this event is alpha triggered" << std::endl;
                //if (trigger_type == 2)
                //    std::cout << "this event is majority triggered" << std::endl;
                continue;
            }
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
    c->SetLogy();
    hist_area_bp->Draw();
    c->SaveAs("hist_area_bp.pdf");

    inputFile->Close();
}

void DataProcessor::saveBinaryOutput() {
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
    if (config_.triggerType == 3)
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
            branchIndex++;
        }
        // Fill the tree with the event's data.
        outputTree->Fill();
    }

    // Write the tree to the file and close it.
    outputTree->Write();
    outputFile->Close();
}

void DataProcessor::run() {
    setSPEResult(config_.inputSPECalibrationPath);
    processFile();
    saveBinaryOutput();
    saveRootOutput();
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

    for (const auto &pmt : pmts_) {
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
        for (size_t i = 0; i < pmts_.size(); ++i) {
            const std::string &ch_name = pmts_[i];
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
            const std::string &ch_name = pmts_[tempPMTCount];
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

    TFile* outputFile = new TFile((config_.outputFilePath + "/output_30t_" + fileID + ".root").c_str(), "RECREATE");
    TTree* outputTree = new TTree("tree", "tree");

    std::vector<UShort_t*> dataStorage;
    std::vector<UShort_t*> trigDataStorage;  // (Unused in this function)

    std::map<std::string, TH1D*> histPMTPE;
    TH1D* histMaxWaveformIndex = new TH1D("histMaxWaveformIndex",
                                          "histMaxWaveformIndex;relative time x 2ns;counts",
                                          200, 0, 1000);
    TH1D* histBotMV = new TH1D("histBotMV", "histBotMV;mV ns;counts", 100, 0, 60000);
    TH1D* histSideMV = new TH1D("histSideMV", "histSideMV;mV ns;counts", 100, 0, 60000);
    TH2D* histBotSideMV = new TH2D("histBotSideMV", (fileID + ";Bot mV;Side mV").c_str(), 100,0,60000, 100, 0, 60000);

    int branchIndex = 0;
    std::vector<double> outputPE(pmts30t.size(), 0.0);

    // Set branch addresses for each PMT channel and create corresponding histograms
    for (auto& pmt : pmts30t) {
        UShort_t* tempArray = new UShort_t[MAX_SAMPLE_SIZE]();  // zero-initialized array
        dataStorage.push_back(tempArray);
        inputTree->SetBranchAddress(pmt.c_str(), tempArray);
        histPMTPE[pmt] = new TH1D(pmt.c_str(), (pmt + ";mV ns;counts").c_str(), 50, 0, 0);
        outputTree->Branch(pmt.c_str(), &outputPE[branchIndex]);
        branchIndex++;
    }

    double bottomPE = 0;
    outputTree->Branch("bottomPE", &bottomPE);
    double sidePE = 0;
    outputTree->Branch("sidePE", &sidePE);

    for (auto &trig : triggers_30t_) {
        UShort_t *tempArray = new UShort_t[MAX_SAMPLE_SIZE]();
        trigDataStorage.push_back(tempArray);
        inputTree->SetBranchAddress(trig.c_str(), tempArray);
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

        bool alpha_triggered = false;
        Waveform alpha(trigDataStorage[0]);
        bool majority_triggered = false;
        Waveform majority(trigDataStorage[1]);
        bool adc_b4_ch23_fired = false;
        Waveform adc_b4_ch23(trigDataStorage[2]);
        if (Waveform::hasValueLessThan(alpha.getSamples(), 3000)) {
            alpha_triggered = true;
        }
        if (Waveform::hasValueLessThan(majority.getSamples(), 3000))
            majority_triggered = true;
        if (Waveform::hasValueLessThan(adc_b4_ch23.getSamples(), 3000))
            adc_b4_ch23_fired = true;

        if (alpha_triggered) {
            if (majority_triggered) {
                std::cout << "alpha majority together" << std::endl;
                alpha.subtractFlatBaseline(0, 100);
                alpha.setAmpPE(1.0, 1.0); //spe_mean, factor = 1.0,1.0
                alpha.correctDaisyChainTrgDelay("adc_b4_ch13");
                majority.subtractFlatBaseline(0, 100);
                majority.setAmpPE(1.0, 1.0); //spe_mean, factor = 1.0,1.0
                majority.correctDaisyChainTrgDelay("adc_b4_ch22");
                TCanvas* can = new TCanvas;
                can->Divide(2,1);
                TGraph* gr2 = majority.drawMVAsGraph(Form("majority_waveform_event_%d", event_id));
                can->cd(1);
                gr2->Draw();
                TGraph* gr = alpha.drawMVAsGraph(Form("alpha_waveform_event_%d", event_id));
                can->cd(2);
                gr->Draw();
                can->SaveAs(Form("adc_b4_ch13_%d.pdf", ievt));
            }
            continue;
        }

        if (adc_b4_ch23_fired) continue;
        if (!majority_triggered) continue;

        int branchIndex = 0;
        // Process waveform for each PMT channel
        for (size_t i = 0; i < pmts30t.size(); ++i) {
            std::string ch_name = pmts30t[i];
            Waveform wf(dataStorage[i]);
            if (wf.getSamples().size() < 100) {
                std::cout << wf.getSamples().size() << std::endl;
                std::cout << ievt << std::endl;
                std::cout << ch_name << std::endl;
                break;
            }
            wf.subtractFlatBaseline(0, 100);
            wf.setAmpPE(1.0, 1.0); //spe_mean, factor = 1.0,1.0
            wf.correctDaisyChainTrgDelay(ch_name);

            // Sum waveforms: use the first channel as the base and add subsequent channels
            if (i == 0) {
                summedWaveform = wf.getAmpMV();
            } else {
                std::vector<double> tempWf = wf.getAmpMV();
                std::transform(summedWaveform.begin(), summedWaveform.end(),
                               tempWf.begin(),
                               summedWaveform.begin(), // store result back in summedWaveform
                               std::plus<double>());
            }

            // Set integration range (here using entire sample length)
            int start = 0;
            int end = MAX_SAMPLE_SIZE;
            double mV_value = 0;
            std::vector<double> ampMV = wf.getAmpMV();
            for (int i = start; i <= ampMV.size() - 1; ++i) {
                mV_value += ampMV[i];
            }
            double pe_value = wf.getPE(start, wf.getAmpPE().size() - 1);
            pe_[ch_name].push_back(pe_value / 2.0);
            histPMTPE[ch_name]->Fill(pe_value / 2.0);
            
            outputPE[branchIndex] = pe_value / 2.0;

            // Accumulate integrated values for bottom and side channels
            if (tempIndex < 12) {
                tempBotMV += mV_value;
            } else {
                tempSideMV += mV_value;
            }
            tempIndex++;
            branchIndex++;
        } // end PMT channel loop

        // Fill histograms for bottom and side integrated values
        histBotMV->Fill(tempBotMV);
        bottomPE = tempBotMV;
        histSideMV->Fill(tempSideMV);
        sidePE = tempSideMV;
        histBotSideMV->Fill(tempBotMV,tempSideMV);

        outputTree->Fill();

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

    TCanvas* c6 = new TCanvas();
    c6->SetLogz();
    histBotSideMV->Draw();
    c6->SaveAs((filename + "_2dMV.pdf").c_str());

    inputFile->Close();

    histBotMV->Write();
    histSideMV->Write();
    histBotSideMV->Write();
    outputTree->Write();
    outputFile->Close();
}
