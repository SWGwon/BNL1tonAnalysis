#include <filesystem>
#include <fstream>
#include <iostream>
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

std::vector<std::string> getFileList(const std::string &directoryPath,
                                     const std::string &searchString = "") {
    std::vector<std::string> fileList;

    try {
        for (const auto &entry : fs::directory_iterator(directoryPath)) {
            if (entry.is_regular_file()) {
                std::string fileName = entry.path().filename().string();

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

std::string cleanString(const std::string &str) {
    std::string cleaned;
    for (char ch : str) {
        if (std::isdigit(ch) || ch == '-') {
            cleaned += ch;
        }
    }
    return cleaned;
}

std::vector<int> getArrayFromFile(const std::string &filename) {
    std::vector<int> array;
    std::ifstream file(filename);

    if (file.is_open()) {
        std::string line;
        std::cout << "Reading event IDs from: " << filename << std::endl;
        while (std::getline(file, line)) {
            std::stringstream ss(line);
            std::string token;
            while (std::getline(ss, token, ',')) {
                std::string cleanedToken = cleanString(token);
                if (!cleanedToken.empty()) {
                    try {
                        int num = std::stoi(cleanedToken);
                        array.push_back(num);
                        std::cout << "Found event ID: " << num
                                  << std::endl; // 디버깅 정보 출력
                    } catch (const std::invalid_argument &e) {
                        std::cerr << "Invalid number in file: " << cleanedToken
                                  << std::endl;
                    }
                }
            }
        }
        file.close();
    } else {
        std::cerr << "Error opening file: " << filename << std::endl;
    }

    return array;
}

std::vector<int> getArrayFromFileCXX(const std::string &filename) {
    std::vector<int> numbers;
    std::ifstream inputFile(filename);

    if (!inputFile) {
        std::cerr << "can not open file: " << filename << std::endl;
        return numbers; 
    }

    int number;
    while (inputFile >> number) {  
        numbers.push_back(number);
    }

    inputFile.close();
    return numbers;
}

void MakeRootFile(const std::string &basePath, std::string File,
                  const std::vector<int> &event_ids, std::string waterLevel) {
    std::cout << "Opening ROOT file: " << File << std::endl;
    TFile *file = TFile::Open(File.c_str(), "READ");
    if (!file || file->IsZombie()) {
        // std::cerr << "Error opening file: " << File << std::endl;
        std::cout << "Attempting to recover the file..." << std::endl;
        file = TFile::Open(File.c_str(), "READ");
        if (file && file->Recover()) {
            std::cout << "File recovered successfully." << std::endl;
        } else {
            std::cerr << "Failed to recover the file." << std::endl;
            if (file)
                file->Close();
            return;
        }
    }

    // Get the list of keys from the input file before creating the output file
    TList *keyList = file->GetListOfKeys();
    if (!keyList) {
        std::cerr << "Error getting list of keys from file: " << File
                  << std::endl;
        file->Close();
        return;
    }

    // Create the output file
    std::string outputFileName = basePath + "data/crossing_muon_" + waterLevel +
                                 File.substr(File.find_last_of("/") + 1);
    // std::cout << "Creating output ROOT file: " << outputFileName <<
    // std::endl;
    TFile *outputFile = TFile::Open(outputFileName.c_str(), "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error creating output file: " << outputFileName
                  << std::endl;
        file->Close();
        return;
    }

    TTree *tree;
    TIter next(keyList);
    TKey *key;

    // Iterate over all trees in the input ROOT file
    while ((key = (TKey *)next())) {
        std::string className = key->GetClassName();
        std::string treeName = key->GetName();

        // Check if the object is a TTree
        if (className == "TTree") {
            tree = (TTree *)file->Get(treeName.c_str());
            tree->SetMaxVirtualSize(0);
            // std::cout << "Processing TTree: " << treeName << std::endl;

            // If the tree is "daq", filter entries based on event IDs
            if (treeName == "daq") {
                TTree *newTree = tree->CloneTree(0); // Clone an empty tree
                tree->SetMaxVirtualSize(0);
                unsigned int event_id;
                tree->SetBranchAddress("event_id", &event_id);

                // Convert event_ids vector to a set for faster lookup
                std::set<int> eventIdSet(event_ids.begin(), event_ids.end());

                // Loop over all entries in the tree
                for (int i = 0; i < tree->GetEntries(); ++i) {
                    // tree->GetEntry(i);
                    int readStatus = tree->GetEntry(i);
                    if (readStatus <= 0) {
                        std::cerr << "Error reading entry " << i << " in tree "
                                  << treeName << std::endl;
                        continue;
                    }

                    // If the event_id matches, add it to the new tree
                    if (eventIdSet.find(event_id) != eventIdSet.end()) {
                        newTree->Fill();
                        // std::cout << "Added event ID: " << event_id << " to
                        // new tree." << std::endl;
                    }
                }

                // Save the new tree to the output file
                outputFile->cd();
                newTree->Write();
                // Do not delete newTree; let ROOT handle it
            } else {
                // For other trees like "run_info", clone them directly
                outputFile->cd();
                TTree *clonedTree = tree->CloneTree();
                clonedTree->Write();
                // std::cout << "Cloned TTree: " << treeName << " as is." <<
                // std::endl;
                //  Do not delete clonedTree; let ROOT handle it
            }
        }
    }

    std::cout << "New tree written to: " << outputFileName << std::endl;
    outputFile->Close();
    file->Close();
}

int main(int argc, char *argv[]) {
    std::string liquid_config =
        ""; // Initial value set to -1 to check if it's provided
    std::string timeStamp = "";

    struct option long_options[] = {
        {"liquid_config", required_argument, 0, 'L'},
        {"timestamp", optional_argument, 0, 'T'},
        {0, 0, 0, 0}
    };

    int opt;
    int option_index = 0;
    while ((opt = getopt_long(argc, argv, "L:C:T:", long_options,
                              &option_index)) != -1) {
        switch (opt) {
        case 'L':
            liquid_config = optarg;
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

    if (liquid_config == "") {
        std::cerr << "Error: --liquid_config (-L) must be provided\n";
        exit(EXIT_FAILURE);
    }

    std::string basePath =
        "/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/";

    std::string directoryPath =
        "/Users/gwon/WbLS/1ton_analysis/drop_apr17_2024/data/RAW_" +
        liquid_config + "/";

    std::vector<std::string> fileNames = getFileList(directoryPath, timeStamp);

    for (auto &f : fileNames) {
        std::cout << f << std::endl;
    }
    for (const std::string &fileName : fileNames) {
        std::cout << "Processing file: " << fileName << std::endl;

        std::string goodEventPath = basePath + "data/crossing_muon_" +
                                    liquid_config + "/" + fileName +
                                    "_good_eventId_with_paddle_CXX.txt";
        std::vector<int> event_ids = getArrayFromFileCXX(goodEventPath);

        std::cout << "event_ids.size(): " << event_ids.size() << std::endl;

        if (event_ids.empty()) {
            continue;
        }

        std::string rootFilePath = directoryPath + fileName + ".root";

        MakeRootFile(basePath, rootFilePath, event_ids, liquid_config + "/");
    }

    return 0;
}
