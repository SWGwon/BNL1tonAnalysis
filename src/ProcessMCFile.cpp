#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "TApplication.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TQObject.h"
#include "TStyle.h"
#include <TFile.h>
#include <TTree.h>

void ProcessFile(std::string inputFileName, std::string outputFileName) {
    TFile *inputFile = new TFile(inputFileName.c_str());
    TTree *inputTree = (TTree *)inputFile->Get("output");

    int mcpdg;
    inputTree->SetBranchAddress("mcpdg", &mcpdg);
    double mcx;
    inputTree->SetBranchAddress("mcx", &mcx);
    double mcy;
    inputTree->SetBranchAddress("mcy", &mcy);
    double mcz;
    inputTree->SetBranchAddress("mcz", &mcz);
    double mcu;
    inputTree->SetBranchAddress("mcu", &mcu);
    double mcv;
    inputTree->SetBranchAddress("mcv", &mcv);
    double mcw;
    inputTree->SetBranchAddress("mcw", &mcw);
    double mcke;
    inputTree->SetBranchAddress("mcke", &mcke);
    double mct;
    inputTree->SetBranchAddress("mct", &mct);
    int evid;
    inputTree->SetBranchAddress("evid", &evid);
    int subev;
    inputTree->SetBranchAddress("subev", &subev);
    int nhits;
    inputTree->SetBranchAddress("nhits", &nhits);
    double triggerTime;
    inputTree->SetBranchAddress("triggerTime", &triggerTime);
    int mcparticlecount;
    inputTree->SetBranchAddress("mcparticlecount", &mcparticlecount);
    int mcpecount;
    inputTree->SetBranchAddress("mcpecount", &mcpecount);
    int mcnhits;
    inputTree->SetBranchAddress("mcnhits", &mcnhits);
    double scintEdep;
    inputTree->SetBranchAddress("scintEdep", &scintEdep);
    double scintEdepQuenched;
    inputTree->SetBranchAddress("scintEdepQuenched", &scintEdepQuenched);
    double scintPhotons;
    inputTree->SetBranchAddress("scintPhotons", &scintPhotons);
    double remPhotons;
    inputTree->SetBranchAddress("remPhotons", &remPhotons);
    double cherPhotons;
    inputTree->SetBranchAddress("cherPhotons", &cherPhotons);
    std::vector<int> *mcpdgs = {};
    inputTree->SetBranchAddress("mcpdgs", &mcpdgs);
    std::vector<double> *mcxs = {};
    inputTree->SetBranchAddress("mcxs", &mcxs);
    std::vector<double> *mcys = {};
    inputTree->SetBranchAddress("mcys", &mcys);
    std::vector<double> *mczs = {};
    inputTree->SetBranchAddress("mczs", &mczs);
    std::vector<double> *mcus = {};
    inputTree->SetBranchAddress("mcus", &mcus);
    std::vector<double> *mcvs = {};
    inputTree->SetBranchAddress("mcvs", &mcvs);
    std::vector<double> *mcws = {};
    inputTree->SetBranchAddress("mcws", &mcws);
    std::vector<double> *mckes = {};
    inputTree->SetBranchAddress("mckes", &mckes);
    std::vector<double> *mcts = {};
    inputTree->SetBranchAddress("mcts", &mcts);
    std::vector<double> *hitPMTTime = {};
    inputTree->SetBranchAddress("hitPMTTime", &hitPMTTime);
    std::vector<double> *hitPMTCharge = {};
    inputTree->SetBranchAddress("hitPMTCharge", &hitPMTCharge);
    std::vector<double> *hitPMTDigitizedTime = {};
    inputTree->SetBranchAddress("hitPMTDigitizedTime", &hitPMTDigitizedTime);
    std::vector<double> *hitPMTDigitizedCharge = {};
    inputTree->SetBranchAddress("hitPMTDigitizedCharge",
                                &hitPMTDigitizedCharge);
    std::vector<int> *hitPMTNCrossings = {};
    inputTree->SetBranchAddress("hitPMTNCrossings", &hitPMTNCrossings);
    std::vector<int> *mcPMTID = {};
    inputTree->SetBranchAddress("mcPMTID", &mcPMTID);
    std::vector<int> *mcPMTNPE = {};
    inputTree->SetBranchAddress("mcPMTNPE", &mcPMTNPE);
    std::vector<double> *mcPETime = {};
    inputTree->SetBranchAddress("mcPETime", &mcPETime);
    std::vector<int> *mcPEProcess = {};
    inputTree->SetBranchAddress("mcPEProcess", &mcPEProcess);
    std::vector<double> *mcPEWavelength = {};
    inputTree->SetBranchAddress("mcPEWavelength", &mcPEWavelength);
    std::vector<double> *mcPEx = {};
    inputTree->SetBranchAddress("mcPEx", &mcPEx);
    std::vector<double> *mcPEy = {};
    inputTree->SetBranchAddress("mcPEy", &mcPEy);
    std::vector<double> *mcPEz = {};
    inputTree->SetBranchAddress("mcPEz", &mcPEz);
    std::vector<int> *hitPMTID = {};
    inputTree->SetBranchAddress("hitPMTID", &hitPMTID);

    TFile *outputFile = new TFile(outputFileName.c_str(), "RECREATE");
    TTree *tree = new TTree("tree", "tree");

    float outputAngle;
    tree->Branch("outputAngle", &outputAngle, "outputAngle/F");
    float E;
    tree->Branch("E", &E, "E/F");
    int PMTPE[46] = {};
    tree->Branch("PMTPE", &PMTPE, "PMTPE[46]/I");

    for (int i = 0; i < inputTree->GetEntries(); ++i) {
        inputTree->GetEntry(i);
        outputAngle = -1;
        E = mcke;
        for (int j = 0; j < 46; ++j) {
            PMTPE[j] = -1;
        }

        for (int j = 0; j < hitPMTID->size(); ++j) {
            int tempPMTId = hitPMTID->at(j);
            int tempPMTPE = mcPMTNPE->at(j);

            PMTPE[tempPMTId] = tempPMTPE;
        }

        tree->Fill();
    }

    tree->Write();
    outputFile->Close();
}

int main(int argc, char **argv) {
    std::string inputFileName;
    std::string outputFileName;

    int c;
    while (1) {
        static struct option long_options[] = {
            {"inputFileName", required_argument, 0, 'i'},
            {"outputFileName", required_argument, 0, 'o'},
            {0, 0, 0, 0}};

        int option_index = 0;
        c = getopt_long(argc, argv, "i:o:", long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {
        case 'i':
            inputFileName = optarg;
            break;
        case 'o':
            outputFileName = optarg;
            break;
        default:
            std::cerr << "Usage: " << argv[0]
                      << " -i inputFileName -o outputFileName" << std::endl;
            exit(1);
        }
    }

    if (inputFileName.empty() || outputFileName.empty()) {
        std::cerr << "Usage: " << argv[0]
                  << " -i inputFileName -o outputFileName" << std::endl;
        return 1;
    }

    ProcessFile(inputFileName, outputFileName);
    return 0;
}
