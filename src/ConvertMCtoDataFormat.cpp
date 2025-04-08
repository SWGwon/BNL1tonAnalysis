#include "AppConfig.h"
#include <iostream>

#include <TFile.h>
#include <TTree.h>

int main(int argc, char** argv) {
    AppConfig config = AppConfig::parseArgs(argc, argv);

    TFile *inputMCFile = new TFile(config.inputFileName.c_str());
    TTree *mcTree = (TTree*)inputMCFile->Get("output");
    std::vector<int>* mcPMTNPE = nullptr;
    mcTree->SetBranchAddress("mcPMTNPE", &mcPMTNPE);
    std::vector<int>* mcPMTID = nullptr;
    mcTree->SetBranchAddress("mcPMTID", &mcPMTID);

    TFile* outputFile = new TFile((config.outputFilePath + "/converted_MC.root").c_str(), "RECREATE");
    TTree* outputTree = new TTree("tree", "tree");
    const int numbin = 80;
    const int pelim = numbin * 20;


    double adc_b1_ch1; outputTree->Branch("adc_b1_ch1", &adc_b1_ch1, "adc_b1_ch1/D");
    double adc_b1_ch2; outputTree->Branch("adc_b1_ch2", &adc_b1_ch2, "adc_b1_ch2/D");
    double adc_b1_ch3; outputTree->Branch("adc_b1_ch3", &adc_b1_ch3, "adc_b1_ch3/D");
    double adc_b1_ch4; outputTree->Branch("adc_b1_ch4", &adc_b1_ch4, "adc_b1_ch4/D");
    double adc_b1_ch5; outputTree->Branch("adc_b1_ch5", &adc_b1_ch5, "adc_b1_ch5/D");
    double adc_b1_ch6; outputTree->Branch("adc_b1_ch6", &adc_b1_ch6, "adc_b1_ch6/D");
    double adc_b1_ch7; outputTree->Branch("adc_b1_ch7", &adc_b1_ch7, "adc_b1_ch7/D");
    double adc_b1_ch8; outputTree->Branch("adc_b1_ch8", &adc_b1_ch8, "adc_b1_ch8/D");
    double adc_b1_ch9; outputTree->Branch("adc_b1_ch9", &adc_b1_ch9, "adc_b1_ch9/D");
    double adc_b1_ch10; outputTree->Branch("adc_b1_ch10", &adc_b1_ch10, "adc_b1_ch10/D");
    double adc_b1_ch11; outputTree->Branch("adc_b1_ch11", &adc_b1_ch11, "adc_b1_ch11/D");
    double adc_b1_ch12; outputTree->Branch("adc_b1_ch12", &adc_b1_ch12, "adc_b1_ch12/D");
    double adc_b1_ch13; outputTree->Branch("adc_b1_ch13", &adc_b1_ch13, "adc_b1_ch13/D");
    double adc_b1_ch14; outputTree->Branch("adc_b1_ch14", &adc_b1_ch14, "adc_b1_ch14/D");
    double adc_b1_ch15; outputTree->Branch("adc_b1_ch15", &adc_b1_ch15, "adc_b1_ch15/D");
    double adc_b2_ch0; outputTree->Branch("adc_b2_ch0", &adc_b2_ch0, "adc_b2_ch0/D");
    double adc_b2_ch1; outputTree->Branch("adc_b2_ch1", &adc_b2_ch1, "adc_b2_ch1/D");
    double adc_b2_ch2; outputTree->Branch("adc_b2_ch2", &adc_b2_ch2, "adc_b2_ch2/D");
    double adc_b2_ch3; outputTree->Branch("adc_b2_ch3", &adc_b2_ch3, "adc_b2_ch3/D");
    double adc_b2_ch4; outputTree->Branch("adc_b2_ch4", &adc_b2_ch4, "adc_b2_ch4/D");
    double adc_b2_ch5; outputTree->Branch("adc_b2_ch5", &adc_b2_ch5, "adc_b2_ch5/D");
    double adc_b2_ch6; outputTree->Branch("adc_b2_ch6", &adc_b2_ch6, "adc_b2_ch6/D");
    double adc_b2_ch7; outputTree->Branch("adc_b2_ch7", &adc_b2_ch7, "adc_b2_ch7/D");
    double adc_b2_ch8; outputTree->Branch("adc_b2_ch8", &adc_b2_ch8, "adc_b2_ch8/D");
    double adc_b2_ch9; outputTree->Branch("adc_b2_ch9", &adc_b2_ch9, "adc_b2_ch9/D");
    double adc_b2_ch10; outputTree->Branch("adc_b2_ch10", &adc_b2_ch10, "adc_b2_ch10/D");
    double adc_b2_ch11; outputTree->Branch("adc_b2_ch11", &adc_b2_ch11, "adc_b2_ch11/D");
    double adc_b2_ch12; outputTree->Branch("adc_b2_ch12", &adc_b2_ch12, "adc_b2_ch12/D");
    double adc_b2_ch13; outputTree->Branch("adc_b2_ch13", &adc_b2_ch13, "adc_b2_ch13/D");
    double adc_b2_ch14; outputTree->Branch("adc_b2_ch14", &adc_b2_ch14, "adc_b2_ch14/D");
    double adc_b3_ch0; outputTree->Branch("adc_b3_ch0", &adc_b3_ch0, "adc_b3_ch0/D");
    double adc_b3_ch1; outputTree->Branch("adc_b3_ch1", &adc_b3_ch1, "adc_b3_ch1/D");
    double adc_b3_ch2; outputTree->Branch("adc_b3_ch2", &adc_b3_ch2, "adc_b3_ch2/D");
    double adc_b3_ch3; outputTree->Branch("adc_b3_ch3", &adc_b3_ch3, "adc_b3_ch3/D");
    double adc_b3_ch4; outputTree->Branch("adc_b3_ch4", &adc_b3_ch4, "adc_b3_ch4/D");
    double adc_b3_ch5; outputTree->Branch("adc_b3_ch5", &adc_b3_ch5, "adc_b3_ch5/D");
    double adc_b3_ch6; outputTree->Branch("adc_b3_ch6", &adc_b3_ch6, "adc_b3_ch6/D");
    double adc_b3_ch7; outputTree->Branch("adc_b3_ch7", &adc_b3_ch7, "adc_b3_ch7/D");
    double adc_b3_ch8; outputTree->Branch("adc_b3_ch8", &adc_b3_ch8, "adc_b3_ch8/D");
    double adc_b3_ch9; outputTree->Branch("adc_b3_ch9", &adc_b3_ch9, "adc_b3_ch9/D");
    double adc_b3_ch10; outputTree->Branch("adc_b3_ch10", &adc_b3_ch10, "adc_b3_ch10/D");
    double adc_b3_ch11; outputTree->Branch("adc_b3_ch11", &adc_b3_ch11, "adc_b3_ch11/D");
    double adc_b3_ch12; outputTree->Branch("adc_b3_ch12", &adc_b3_ch12, "adc_b3_ch12/D");
    double adc_b3_ch13; outputTree->Branch("adc_b3_ch13", &adc_b3_ch13, "adc_b3_ch13/D");
    double adc_b3_ch14; outputTree->Branch("adc_b3_ch14", &adc_b3_ch14, "adc_b3_ch14/D");
    double adc_b3_ch15; outputTree->Branch("adc_b3_ch15", &adc_b3_ch15, "adc_b3_ch15/D");
    double adc_b4_ch0; outputTree->Branch("adc_b4_ch0", &adc_b4_ch0, "adc_b4_ch0/D");
    double adc_b4_ch1; outputTree->Branch("adc_b4_ch1", &adc_b4_ch1, "adc_b4_ch1/D");
    double adc_b4_ch2; outputTree->Branch("adc_b4_ch2", &adc_b4_ch2, "adc_b4_ch2/D");
    double adc_b4_ch3; outputTree->Branch("adc_b4_ch3", &adc_b4_ch3, "adc_b4_ch3/D");
    double adc_b4_ch4; outputTree->Branch("adc_b4_ch4", &adc_b4_ch4, "adc_b4_ch4/D");
    double adc_b4_ch5; outputTree->Branch("adc_b4_ch5", &adc_b4_ch5, "adc_b4_ch5/D");
    double adc_b4_ch6; outputTree->Branch("adc_b4_ch6", &adc_b4_ch6, "adc_b4_ch6/D");
    double adc_b4_ch7; outputTree->Branch("adc_b4_ch7", &adc_b4_ch7, "adc_b4_ch7/D");
    double adc_b4_ch8; outputTree->Branch("adc_b4_ch8", &adc_b4_ch8, "adc_b4_ch8/D");
    double adc_b4_ch9; outputTree->Branch("adc_b4_ch9", &adc_b4_ch9, "adc_b4_ch9/D");
    double adc_b4_ch10; outputTree->Branch("adc_b4_ch10", &adc_b4_ch10, "adc_b4_ch10/D");
    double adc_b4_ch11; outputTree->Branch("adc_b4_ch11", &adc_b4_ch11, "adc_b4_ch11/D");


    int event_number = 0;
    if (config.eventNumber > mcTree->GetEntries())
        event_number = mcTree->GetEntries();
    else
        event_number = config.eventNumber;

    for (int ievt = 0; ievt < event_number; ++ievt) {
        mcTree->GetEntry(ievt);

        adc_b1_ch1 = 0;
        adc_b1_ch2 = 0;
        adc_b1_ch3 = 0;
        adc_b1_ch4 = 0;
        adc_b1_ch5 = 0;
        adc_b1_ch6 = 0;
        adc_b1_ch7 = 0;
        adc_b1_ch8 = 0;
        adc_b1_ch9 = 0;
        adc_b1_ch10 = 0;
        adc_b1_ch11 = 0;
        adc_b1_ch12 = 0;
        adc_b1_ch13 = 0;
        adc_b1_ch14 = 0;
        adc_b1_ch15 = 0;
        adc_b2_ch0 = 0;
        adc_b2_ch1 = 0;
        adc_b2_ch2 = 0;
        adc_b2_ch3 = 0;
        adc_b2_ch4 = 0;
        adc_b2_ch5 = 0;
        adc_b2_ch6 = 0;
        adc_b2_ch7 = 0;
        adc_b2_ch8 = 0;
        adc_b2_ch9 = 0;
        adc_b2_ch10 = 0;
        adc_b2_ch11 = 0;
        adc_b2_ch12 = 0;
        adc_b2_ch13 = 0;
        adc_b2_ch14 = 0;
        adc_b3_ch0 = 0;
        adc_b3_ch1 = 0;
        adc_b3_ch2 = 0;
        adc_b3_ch3 = 0;
        adc_b3_ch4 = 0;
        adc_b3_ch5 = 0;
        adc_b3_ch6 = 0;
        adc_b3_ch7 = 0;
        adc_b3_ch8 = 0;
        adc_b3_ch9 = 0;
        adc_b3_ch10 = 0;
        adc_b3_ch11 = 0;
        adc_b3_ch12 = 0;
        adc_b3_ch13 = 0;
        adc_b3_ch14 = 0;
        adc_b3_ch15 = 0;
        adc_b4_ch0 = 0;
        adc_b4_ch1 = 0;
        adc_b4_ch2 = 0;
        adc_b4_ch3 = 0;
        adc_b4_ch4 = 0;
        adc_b4_ch5 = 0;
        adc_b4_ch6 = 0;
        adc_b4_ch7 = 0;
        adc_b4_ch8 = 0;
        adc_b4_ch9 = 0;
        adc_b4_ch10 = 0;
        adc_b4_ch11 = 0;

        try {
            adc_b1_ch1 = mcPMTNPE->at(0);
        } catch (...) {
            adc_b1_ch1 = 0;
        }
        try {
            adc_b1_ch2 = mcPMTNPE->at(1);
        } catch (...) {
            adc_b1_ch2 = 0;
        }
        try {
            adc_b1_ch3 = mcPMTNPE->at(2);
        } catch (...) {
            adc_b1_ch3 = 0;
        }
        try {
            adc_b1_ch4 = mcPMTNPE->at(3);
        } catch (...) {
            adc_b1_ch4 = 0;
        }
        try {
            adc_b1_ch5 = mcPMTNPE->at(4);
        } catch (...) {
            adc_b1_ch5 = 0;
        }
        try {
            adc_b1_ch6 = mcPMTNPE->at(5);
        } catch (...) {
            adc_b1_ch6 = 0;
        }
        try {
            adc_b1_ch7 = mcPMTNPE->at(6);
        } catch (...) {
            adc_b1_ch7 = 0;
        }
        try {
            adc_b1_ch8 = mcPMTNPE->at(8);
        } catch (...) {
            adc_b1_ch8 = 0;
        }
        try {
            adc_b1_ch9 = mcPMTNPE->at(9);
        } catch (...) {
            adc_b1_ch9 = 0;
        }
        try {
            adc_b1_ch10 = mcPMTNPE->at(9);
        } catch (...) {
            adc_b1_ch10 = 0;
        }
        try {
            adc_b1_ch11 = mcPMTNPE->at(10);
        } catch (...) {
            adc_b1_ch11 = 0;
        }
        try {
            adc_b1_ch12 = mcPMTNPE->at(11);
        } catch (...) {
            adc_b1_ch12 = 0;
        }
        try {
            adc_b1_ch13 = mcPMTNPE->at(12);
        } catch (...) {
            adc_b1_ch13 = 0;
        }
        try {
            adc_b1_ch14 = mcPMTNPE->at(13);
        } catch (...) {
            adc_b1_ch14 = 0;
        }
        try {
            adc_b1_ch15 = mcPMTNPE->at(14);
        } catch (...) {
            adc_b1_ch15 = 0;
        }

        try {
            adc_b2_ch0 = mcPMTNPE->at(15);
        } catch (...) {
            adc_b2_ch0 = 0;
        }
        try {
            adc_b2_ch1 = mcPMTNPE->at(16);
        } catch (...) {
            adc_b2_ch1 = 0;
        }
        try {
            adc_b2_ch2 = mcPMTNPE->at(17);
        } catch (...) {
            adc_b2_ch2 = 0;
        }
        try {
            adc_b2_ch3 = mcPMTNPE->at(18);
        } catch (...) {
            adc_b2_ch3 = 0;
        }
        try {
            adc_b2_ch4 = mcPMTNPE->at(19);
        } catch (...) {
            adc_b2_ch4 = 0;
        }
        try {
            adc_b2_ch5 = mcPMTNPE->at(20);
        } catch (...) {
            adc_b2_ch5 = 0;
        }
        try {
            adc_b2_ch6 = mcPMTNPE->at(21);
        } catch (...) {
            adc_b2_ch6 = 0;
        }
        try {
            adc_b2_ch7 = mcPMTNPE->at(22);
        } catch (...) {
            adc_b2_ch7 = 0;
        }
        try {
            adc_b2_ch8 = mcPMTNPE->at(23);
        } catch (...) {
            adc_b2_ch8 = 0;
        }
        try {
            adc_b2_ch9 = mcPMTNPE->at(24);
        } catch (...) {
            adc_b2_ch9 = 0;
        }
        try {
            adc_b2_ch10 = mcPMTNPE->at(25);
        } catch (...) {
            adc_b2_ch10 = 0;
        }
        try {
            adc_b2_ch11 = mcPMTNPE->at(26);
        } catch (...) {
            adc_b2_ch11 = 0;
        }
        try {
            adc_b2_ch12 = mcPMTNPE->at(27);
        } catch (...) {
            adc_b2_ch12 = 0;
        }
        try {
            adc_b2_ch13 = mcPMTNPE->at(28);
        } catch (...) {
            adc_b2_ch13 = 0;
        }
        try {
            adc_b2_ch14 = mcPMTNPE->at(29);
        } catch (...) {
            adc_b2_ch14 = 0;
        }

        try {
            adc_b3_ch0 = mcPMTNPE->at(34);
        } catch (...) {
            adc_b3_ch0 = 0;
        }
        try {
            adc_b3_ch1 = mcPMTNPE->at(35);
        } catch (...) {
            adc_b3_ch1 = 0;
        }
        try {
            adc_b3_ch2 = mcPMTNPE->at(36);
        } catch (...) {
            adc_b3_ch2 = 0;
        }
        try {
            adc_b3_ch3 = mcPMTNPE->at(37);
        } catch (...) {
            adc_b3_ch3 = 0;
        }

        try {
            adc_b3_ch4 = mcPMTNPE->at(38);
        } catch (...) {
            adc_b3_ch4 = 0;
        }
        try {
            adc_b3_ch5 = mcPMTNPE->at(39);
        } catch (...) {
            adc_b3_ch5 = 0;
        }
        try {
            adc_b3_ch6 = mcPMTNPE->at(40);
        } catch (...) {
            adc_b3_ch6 = 0;
        }
        try {
            adc_b3_ch7 = mcPMTNPE->at(41);
        } catch (...) {
            adc_b3_ch7 = 0;
        }

        try {
            adc_b3_ch8 = mcPMTNPE->at(42);
        } catch (...) {
            adc_b3_ch8 = 0;
        }
        try {
            adc_b3_ch9 = mcPMTNPE->at(43);
        } catch (...) {
            adc_b3_ch9 = 0;
        }
        try {
            adc_b3_ch10 = mcPMTNPE->at(44);
        } catch (...) {
            adc_b3_ch10 = 0;
        }
        try {
            adc_b3_ch11 = mcPMTNPE->at(45);
        } catch (...) {
            adc_b3_ch11 = 0;
        }

        try {
            adc_b3_ch12 = mcPMTNPE->at(30);
        } catch (...) {
            adc_b3_ch12 = 0;
        }
        try {
            adc_b3_ch13 = mcPMTNPE->at(31);
        } catch (...) {
            adc_b3_ch13 = 0;
        }
        try {
            adc_b3_ch14 = mcPMTNPE->at(32);
        } catch (...) {
            adc_b3_ch14 = 0;
        }
        try {
            adc_b3_ch15 = mcPMTNPE->at(33);
        } catch (...) {
            adc_b3_ch15 = 0;
        }
        
        //adc_b4_ch0 = mcPMTNPE->at(0);
        //adc_b4_ch1 = mcPMTNPE->at(0);
        //adc_b4_ch2 = mcPMTNPE->at(0);
        //adc_b4_ch3 = mcPMTNPE->at(0);
        //adc_b4_ch4 = mcPMTNPE->at(0);
        //adc_b4_ch5 = mcPMTNPE->at(0);
        //adc_b4_ch6 = mcPMTNPE->at(0);
        //adc_b4_ch7 = mcPMTNPE->at(0);
        //adc_b4_ch8 = mcPMTNPE->at(0);
        //adc_b4_ch9 = mcPMTNPE->at(0);
        //adc_b4_ch10 = mcPMTNPE->at(0);
        //adc_b4_ch11 = mcPMTNPE->at(0);

        outputTree->Fill();
    }

    outputTree->Write();
    outputFile->Close();
}
