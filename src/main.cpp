#include <Eigen/Dense>
#include <OsqpEigen/OsqpEigen.h>
#include <iomanip>
#include <iostream>
#include <vector>

#include "TApplication.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TLine.h"
#include "TQObject.h"
#include "TStyle.h"
#include <Eigen/Dense>
#include <TFile.h>
#include <TTree.h>
#include <qpOASES.hpp>

USING_NAMESPACE_QPOASES

int TOTAL_PMT_COUNT = 46;
int TOTAL_NPE_UPPER_CUT = 100000;
int TOTAL_NPE_LOWER_CUT = -100;
TH1D *histTotalNPEData;
TH1D *histTotalNPEMC;

int PE_CUT_EACH_PMT = 50;

std::vector<double> solveForF(const std::vector<std::vector<double>> &M,
                              const std::vector<double> &D,
                              const std::vector<double> &F_min,
                              const std::vector<double> &F_max) {
    int n = M.size();    // Number of rows in M
    int m = M[0].size(); // Number of columns in M

    if (D.size() != n) {
        throw std::invalid_argument(
            "Size of vector D does not match number of rows in matrix M.");
    }

    if (F_min.size() != m || F_max.size() != m) {
        throw std::invalid_argument(
            "Size of F_min and F_max must match number of variables.");
    }

    // Convert std::vector to Eigen matrices
    Eigen::MatrixXd matrixM(n, m);
    Eigen::VectorXd vectorD(n);
    for (int i = 0; i < n; ++i) {
        vectorD(i) = D[i];
        for (int j = 0; j < m; ++j) {
            matrixM(i, j) = M[i][j];
        }
    }

    // Set up the QP problem:
    // Minimize (1/2) * F^T * H * F + q^T * F
    // Subject to lb <= F <= ub

    // Compute H = 2 * M^T * M
    Eigen::MatrixXd H = 2.0 * matrixM.transpose() * matrixM;

    // Compute q = -2 * M^T * D
    Eigen::VectorXd q = -2.0 * matrixM.transpose() * vectorD;

    // Set up the bounds
    Eigen::VectorXd lb =
        Eigen::Map<const Eigen::VectorXd>(F_min.data(), F_min.size());
    Eigen::VectorXd ub =
        Eigen::Map<const Eigen::VectorXd>(F_max.data(), F_max.size());

    // Set up the OSQP solver
    OsqpEigen::Solver solver;
    solver.settings()->setVerbosity(false);
    solver.settings()->setWarmStart(true);

    // Since we have only variable bounds, we set constraints matrix A to
    // identity
    Eigen::SparseMatrix<double> A(m, m);
    A.setIdentity();

    // Convert H to sparse matrix
    Eigen::SparseMatrix<double> H_sparse = H.sparseView();

    // Initialize solver data
    solver.data()->setNumberOfVariables(m);
    solver.data()->setNumberOfConstraints(m);
    solver.data()->setHessianMatrix(H_sparse);
    solver.data()->setGradient(q);
    solver.data()->setLinearConstraintsMatrix(A);
    solver.data()->setLowerBound(lb);
    solver.data()->setUpperBound(ub);

    // Check if data is ready
    // if (!solver.isDataReady()) {
    //    throw std::runtime_error("Solver data is not ready.");
    //}

    // Initialize solver
    if (!solver.initSolver()) {
        throw std::runtime_error("Failed to initialize solver.");
    }

    // Solve the problem
    if (solver.solveProblem() != OsqpEigen::ErrorExitFlag::NoError) {
        throw std::runtime_error("Solver failed to find a solution.");
    }

    // Retrieve the solution
    Eigen::VectorXd vectorF = solver.getSolution();

    // Convert the solution to std::vector
    std::vector<double> resultF(m);
    for (int i = 0; i < m; ++i) {
        resultF[i] = vectorF(i);
    }

    return resultF;
}
// std::vector<double> solveForF(const std::vector<std::vector<double>>& M,
//                               const std::vector<double>& D) {
//     int n = M.size();  // 행렬 M의 행 크기 (n x m에서 n)
//     int m = M[0].size();  // 행렬 M의 열 크기 (n x m에서 m)
//
//     if (D.size() != n) {
//         std::cout << "D.size(): " << D.size() << std::endl;
//         std::cout << "n: " << n << std::endl;
//         throw std::invalid_argument("D 벡터의 크기가 행렬 M의 행 크기와 맞지
//         않습니다.");
//     }
//
//     // Eigen 행렬과 벡터로 변환
//     Eigen::MatrixXd matrixM(n, m);
//     Eigen::VectorXd vectorD(n);  // D는 n차원의 열 벡터
//     for (int i = 0; i < n; ++i) {
//         for (int j = 0; j < m; ++j) {
//             matrixM(i, j) = M[i][j];
//         }
//         vectorD(i) = D[i];
//     }
//
//     // QR 분해를 사용하여 선형 방정식 풀기 (M * F = D)
//     Eigen::VectorXd vectorF = matrixM.colPivHouseholderQr().solve(vectorD);
//
//     // 결과를 std::vector로 변환
//     std::vector<double> resultF(m);
//     for (int j = 0; j < m; ++j) {
//         resultF[j] = vectorF(j);
//     }
//
//     return resultF;
// }

void printMatrix(const std::vector<std::vector<double>> &matrix) {
    // 행렬의 크기를 알아내기 위해 행과 열의 수 계산
    size_t rows = matrix.size();
    size_t cols = matrix[0].size();

    // 행렬 출력
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            // setw(10)으로 칸을 맞추고 setprecision(4)로 소수점 4자리까지 출력
            // std::cout << std::setw(1) << std::setprecision(0) << matrix[i][j]
            // << " ";
            std::cout << std::setw(2) << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void printVector(const std::vector<double> &vector) {
    // 행렬의 크기를 알아내기 위해 행과 열의 수 계산
    size_t cols = vector.size();

    // 행렬 출력
    for (size_t j = 0; j < cols; ++j) {
        // setw(10)으로 칸을 맞추고 setprecision(4)로 소수점 4자리까지 출력
        // std::cout << std::setw(1) << std::setprecision(0) << matrix[i][j] <<
        // " ";
        std::cout << std::setw(2) << vector[j] << " ";
    }
    std::cout << std::endl;
}

std::vector<std::vector<double>> MakeResponseMatrixMC(TFile *inputFile,
                                                      int event_num) {
    std::vector<std::vector<double>> Response_matrix = {};

    TTree *tree = (TTree *)inputFile->Get("output");

    int mcpdg;
    tree->SetBranchAddress("mcpdg", &mcpdg);
    double mcx;
    tree->SetBranchAddress("mcx", &mcx);
    double mcy;
    tree->SetBranchAddress("mcy", &mcy);
    double mcz;
    tree->SetBranchAddress("mcz", &mcz);
    double mcu;
    tree->SetBranchAddress("mcu", &mcu);
    double mcv;
    tree->SetBranchAddress("mcv", &mcv);
    double mcw;
    tree->SetBranchAddress("mcw", &mcw);
    double mcke;
    tree->SetBranchAddress("mcke", &mcke);
    double mct;
    tree->SetBranchAddress("mct", &mct);
    int evid;
    tree->SetBranchAddress("evid", &evid);
    int subev;
    tree->SetBranchAddress("subev", &subev);
    int nhits;
    tree->SetBranchAddress("nhits", &nhits);
    double triggerTime;
    tree->SetBranchAddress("triggerTime", &triggerTime);
    int mcparticlecount;
    tree->SetBranchAddress("mcparticlecount", &mcparticlecount);
    int mcpecount;
    tree->SetBranchAddress("mcpecount", &mcpecount);
    int mcnhits;
    tree->SetBranchAddress("mcnhits", &mcnhits);
    double scintEdep;
    tree->SetBranchAddress("scintEdep", &scintEdep);
    double scintEdepQuenched;
    tree->SetBranchAddress("scintEdepQuenched", &scintEdepQuenched);
    double scintPhotons;
    tree->SetBranchAddress("scintPhotons", &scintPhotons);
    double remPhotons;
    tree->SetBranchAddress("remPhotons", &remPhotons);
    double cherPhotons;
    tree->SetBranchAddress("cherPhotons", &cherPhotons);
    std::vector<int> *mcpdgs = {};
    tree->SetBranchAddress("mcpdgs", &mcpdgs);
    std::vector<double> *mcxs = {};
    tree->SetBranchAddress("mcxs", &mcxs);
    std::vector<double> *mcys = {};
    tree->SetBranchAddress("mcys", &mcys);
    std::vector<double> *mczs = {};
    tree->SetBranchAddress("mczs", &mczs);
    std::vector<double> *mcus = {};
    tree->SetBranchAddress("mcus", &mcus);
    std::vector<double> *mcvs = {};
    tree->SetBranchAddress("mcvs", &mcvs);
    std::vector<double> *mcws = {};
    tree->SetBranchAddress("mcws", &mcws);
    std::vector<double> *mckes = {};
    tree->SetBranchAddress("mckes", &mckes);
    std::vector<double> *mcts = {};
    tree->SetBranchAddress("mcts", &mcts);
    std::vector<double> *hitPMTTime = {};
    tree->SetBranchAddress("hitPMTTime", &hitPMTTime);
    std::vector<double> *hitPMTCharge = {};
    tree->SetBranchAddress("hitPMTCharge", &hitPMTCharge);
    std::vector<double> *hitPMTDigitizedTime = {};
    tree->SetBranchAddress("hitPMTDigitizedTime", &hitPMTDigitizedTime);
    std::vector<double> *hitPMTDigitizedCharge = {};
    tree->SetBranchAddress("hitPMTDigitizedCharge", &hitPMTDigitizedCharge);
    std::vector<int> *hitPMTNCrossings = {};
    tree->SetBranchAddress("hitPMTNCrossings", &hitPMTNCrossings);
    std::vector<int> *mcPMTID = {};
    tree->SetBranchAddress("mcPMTID", &mcPMTID);
    std::vector<int> *mcPMTNPE = {};
    tree->SetBranchAddress("mcPMTNPE", &mcPMTNPE);
    std::vector<double> *mcPETime = {};
    tree->SetBranchAddress("mcPETime", &mcPETime);
    std::vector<int> *mcPEProcess = {};
    tree->SetBranchAddress("mcPEProcess", &mcPEProcess);
    std::vector<double> *mcPEWavelength = {};
    tree->SetBranchAddress("mcPEWavelength", &mcPEWavelength);
    std::vector<double> *mcPEx = {};
    tree->SetBranchAddress("mcPEx", &mcPEx);
    std::vector<double> *mcPEy = {};
    tree->SetBranchAddress("mcPEy", &mcPEy);
    std::vector<double> *mcPEz = {};
    tree->SetBranchAddress("mcPEz", &mcPEz);
    std::vector<int> *hitPMTID = {};
    tree->SetBranchAddress("hitPMTID", &hitPMTID);

    for (int ievt = 0; ievt < event_num; ++ievt) {
        tree->GetEntry(ievt);

        int tempTotalNPE = 0;
        std::vector<double> tempMC_nPE_PMT(TOTAL_PMT_COUNT);
        for (int j = 0; j < hitPMTID->size(); ++j) {
            if (hitPMTID->at(j) > TOTAL_PMT_COUNT)
                continue;
            tempTotalNPE += mcPMTNPE->at(j);
            tempMC_nPE_PMT.at(hitPMTID->at(j)) = mcPMTNPE->at(j);
        }
        Response_matrix.push_back(tempMC_nPE_PMT);
    }

    return Response_matrix;
}

void MakeResponseMatrixMC2(TFile *inputFile, int event_num,
                           std::vector<std::vector<double>> &Response_matrix) {
    TTree *tree = (TTree *)inputFile->Get("tree");
    if (!tree) {
        std::cerr << "Error: TTree 'tree' not found in the input file."
                  << std::endl;
        return;
    } else {
        std::cout << "Successfully retrieved TTree 'tree'." << std::endl;
    }

    int PMTPE[46];
    if (tree->SetBranchAddress("PMTPE", PMTPE) != 0) {
        std::cerr << "Error: Could not set branch address for 'PMTPE'."
                  << std::endl;
        return;
    } else {
        std::cout << "Successfully set branch address for 'PMTPE'."
                  << std::endl;
    }

    int tempEventNum = 0;
    Long64_t totalEntries = tree->GetEntries();
    std::cout << "Total entries in tree: " << totalEntries << std::endl;

    for (Long64_t ievt = 0; ievt < totalEntries; ++ievt) {
        std::cout << "Processing event " << ievt << std::endl;
        Long64_t nb = tree->GetEntry(ievt);
        if (nb <= 0) {
            std::cerr << "Warning: No bytes read for entry " << ievt
                      << std::endl;
            continue;
        }

        int tempTotalNPE = 0;
        std::vector<double> tempMC_nPE_PMT(TOTAL_PMT_COUNT);
        for (int j = 0; j < TOTAL_PMT_COUNT; ++j) {
            tempTotalNPE += PMTPE[j];
            tempMC_nPE_PMT[j] = PMTPE[j];
        }

        std::cout << "Total NPE for event " << ievt << ": " << tempTotalNPE
                  << std::endl;

        // if (tempTotalNPE < TOTAL_NPE_LOWER_CUT) {
        //     std::cout << "Total NPE below lower cut. Skipping event." <<
        //     std::endl; continue;
        // }
        // if (tempTotalNPE > TOTAL_NPE_UPPER_CUT) {
        //     std::cout << "Total NPE above upper cut. Skipping event." <<
        //     std::endl; continue;
        // }

        Response_matrix.push_back(tempMC_nPE_PMT);
        tempEventNum++;
        if (tempEventNum == event_num)
            break;
    }

    std::cout << "Processed " << tempEventNum << " events." << std::endl;

    return;
}

void MakeResponseMatrixMC3(TFile *inputFile, int event_num,
                           std::vector<std::vector<double>> &Response_matrix) {
    TTree *tree = (TTree *)inputFile->Get("output");
    if (!tree) {
        std::cerr << "Error: TTree 'output' not found in the input file."
                  << std::endl;
        return;
    } else {
        std::cout << "Successfully retrieved TTree 'output'." << std::endl;
    }

    std::vector<int> *mcPMTNPE = {};
    if (tree->SetBranchAddress("mcPMTNPE", &mcPMTNPE) != 0) {
        std::cerr << "Error: Could not set branch address for 'mcPMTNPE'."
                  << std::endl;
        return;
    } else {
        std::cout << "Successfully set branch address for 'mcPMTNPE'."
                  << std::endl;
    }
    std::vector<int> *mcPMTID = {};
    if (tree->SetBranchAddress("mcPMTID", &mcPMTID) != 0) {
        std::cerr << "Error: Could not set branch address for 'mcPMTID'."
                  << std::endl;
        return;
    } else {
        std::cout << "Successfully set branch address for 'mcPMTID'."
                  << std::endl;
    }

    int tempEventNum = 0;
    Long64_t totalEntries = tree->GetEntries();
    std::cout << "Total entries in tree: " << totalEntries << std::endl;

    for (Long64_t ievt = 0; ievt < totalEntries; ++ievt) {
        std::cout << "Processing event " << ievt << std::endl;
        Long64_t nb = tree->GetEntry(ievt);
        if (nb <= 0) {
            std::cerr << "Warning: No bytes read for entry " << ievt
                      << std::endl;
            continue;
        }

        int tempTotalNPE = 0;
        std::vector<double> tempMC_nPE_PMT(TOTAL_PMT_COUNT);
        for (int j = 0; j < TOTAL_PMT_COUNT; ++j) {
            tempTotalNPE += mcPMTNPE->at(j);
            tempMC_nPE_PMT[j] = mcPMTNPE->at(j);
        }

        std::cout << "Total NPE for event " << ievt << ": " << tempTotalNPE
                  << std::endl;

        // if (tempTotalNPE < TOTAL_NPE_LOWER_CUT) {
        //     std::cout << "Total NPE below lower cut. Skipping event." <<
        //     std::endl; continue;
        // }
        // if (tempTotalNPE > TOTAL_NPE_UPPER_CUT) {
        //     std::cout << "Total NPE above upper cut. Skipping event." <<
        //     std::endl; continue;
        // }

        Response_matrix.push_back(tempMC_nPE_PMT);
        tempEventNum++;
        if (tempEventNum == event_num)
            break;
    }

    std::cout << "Processed " << tempEventNum << " events." << std::endl;

    return;
}

std::vector<std::vector<double>> MakeResponseMatrixData(TFile *inputFile,
                                                        int event_num) {
    std::vector<std::vector<double>> Response_matrix = {};

    TTree *tree = (TTree *)inputFile->Get("tree");
    double adc_b1_ch1;
    tree->SetBranchAddress("adc_b1_ch1", &adc_b1_ch1);
    double adc_b1_ch2;
    tree->SetBranchAddress("adc_b1_ch2", &adc_b1_ch2);
    double adc_b1_ch3;
    tree->SetBranchAddress("adc_b1_ch3", &adc_b1_ch3);
    double adc_b1_ch4;
    tree->SetBranchAddress("adc_b1_ch4", &adc_b1_ch4);
    double adc_b1_ch5;
    tree->SetBranchAddress("adc_b1_ch5", &adc_b1_ch5);
    double adc_b1_ch6;
    tree->SetBranchAddress("adc_b1_ch6", &adc_b1_ch6);
    double adc_b1_ch7;
    tree->SetBranchAddress("adc_b1_ch7", &adc_b1_ch7);
    double adc_b1_ch8;
    tree->SetBranchAddress("adc_b1_ch8", &adc_b1_ch8);
    double adc_b1_ch9;
    tree->SetBranchAddress("adc_b1_ch9", &adc_b1_ch9);
    double adc_b1_ch10;
    tree->SetBranchAddress("adc_b1_ch10", &adc_b1_ch10);
    double adc_b1_ch11;
    tree->SetBranchAddress("adc_b1_ch11", &adc_b1_ch11);
    double adc_b1_ch12;
    tree->SetBranchAddress("adc_b1_ch12", &adc_b1_ch12);
    double adc_b1_ch13;
    tree->SetBranchAddress("adc_b1_ch13", &adc_b1_ch13);
    double adc_b1_ch14;
    tree->SetBranchAddress("adc_b1_ch14", &adc_b1_ch14);
    double adc_b1_ch15;
    tree->SetBranchAddress("adc_b1_ch15", &adc_b1_ch15);
    double adc_b2_ch0;
    tree->SetBranchAddress("adc_b2_ch0", &adc_b2_ch0);
    double adc_b2_ch1;
    tree->SetBranchAddress("adc_b2_ch1", &adc_b2_ch1);
    double adc_b2_ch2;
    tree->SetBranchAddress("adc_b2_ch2", &adc_b2_ch2);
    double adc_b2_ch3;
    tree->SetBranchAddress("adc_b2_ch3", &adc_b2_ch3);
    double adc_b2_ch4;
    tree->SetBranchAddress("adc_b2_ch4", &adc_b2_ch4);
    double adc_b2_ch5;
    tree->SetBranchAddress("adc_b2_ch5", &adc_b2_ch5);
    double adc_b2_ch6;
    tree->SetBranchAddress("adc_b2_ch6", &adc_b2_ch6);
    double adc_b2_ch7;
    tree->SetBranchAddress("adc_b2_ch7", &adc_b2_ch7);
    double adc_b2_ch8;
    tree->SetBranchAddress("adc_b2_ch8", &adc_b2_ch8);
    double adc_b2_ch9;
    tree->SetBranchAddress("adc_b2_ch9", &adc_b2_ch9);
    double adc_b2_ch10;
    tree->SetBranchAddress("adc_b2_ch10", &adc_b2_ch10);
    double adc_b2_ch11;
    tree->SetBranchAddress("adc_b2_ch11", &adc_b2_ch11);
    double adc_b2_ch12;
    tree->SetBranchAddress("adc_b2_ch12", &adc_b2_ch12);
    double adc_b2_ch13;
    tree->SetBranchAddress("adc_b2_ch13", &adc_b2_ch13);
    double adc_b2_ch14;
    tree->SetBranchAddress("adc_b2_ch14", &adc_b2_ch14);
    double adc_b3_ch0;
    tree->SetBranchAddress("adc_b3_ch0", &adc_b3_ch0);
    double adc_b3_ch1;
    tree->SetBranchAddress("adc_b3_ch1", &adc_b3_ch1);
    double adc_b3_ch2;
    tree->SetBranchAddress("adc_b3_ch2", &adc_b3_ch2);
    double adc_b3_ch3;
    tree->SetBranchAddress("adc_b3_ch3", &adc_b3_ch3);
    double adc_b3_ch4;
    tree->SetBranchAddress("adc_b3_ch4", &adc_b3_ch4);
    double adc_b3_ch5;
    tree->SetBranchAddress("adc_b3_ch5", &adc_b3_ch5);
    double adc_b3_ch6;
    tree->SetBranchAddress("adc_b3_ch6", &adc_b3_ch6);
    double adc_b3_ch7;
    tree->SetBranchAddress("adc_b3_ch7", &adc_b3_ch7);
    double adc_b3_ch8;
    tree->SetBranchAddress("adc_b3_ch8", &adc_b3_ch8);
    double adc_b3_ch9;
    tree->SetBranchAddress("adc_b3_ch9", &adc_b3_ch9);
    double adc_b3_ch10;
    tree->SetBranchAddress("adc_b3_ch10", &adc_b3_ch10);
    double adc_b3_ch11;
    tree->SetBranchAddress("adc_b3_ch11", &adc_b3_ch11);
    double adc_b3_ch12;
    tree->SetBranchAddress("adc_b3_ch12", &adc_b3_ch12);
    double adc_b3_ch13;
    tree->SetBranchAddress("adc_b3_ch13", &adc_b3_ch13);
    double adc_b3_ch14;
    tree->SetBranchAddress("adc_b3_ch14", &adc_b3_ch14);
    double adc_b3_ch15;
    tree->SetBranchAddress("adc_b3_ch15", &adc_b3_ch15);
    double adc_b4_ch0;
    tree->SetBranchAddress("adc_b4_ch0", &adc_b4_ch0);
    double adc_b4_ch1;
    tree->SetBranchAddress("adc_b4_ch1", &adc_b4_ch1);
    double adc_b4_ch2;
    tree->SetBranchAddress("adc_b4_ch2", &adc_b4_ch2);
    double adc_b4_ch3;
    tree->SetBranchAddress("adc_b4_ch3", &adc_b4_ch3);
    double adc_b4_ch4;
    tree->SetBranchAddress("adc_b4_ch4", &adc_b4_ch4);
    double adc_b4_ch5;
    tree->SetBranchAddress("adc_b4_ch5", &adc_b4_ch5);
    double adc_b4_ch6;
    tree->SetBranchAddress("adc_b4_ch6", &adc_b4_ch6);
    double adc_b4_ch7;
    tree->SetBranchAddress("adc_b4_ch7", &adc_b4_ch7);
    double adc_b4_ch8;
    tree->SetBranchAddress("adc_b4_ch8", &adc_b4_ch8);
    double adc_b4_ch9;
    tree->SetBranchAddress("adc_b4_ch9", &adc_b4_ch9);
    double adc_b4_ch10;
    tree->SetBranchAddress("adc_b4_ch10", &adc_b4_ch10);
    double adc_b4_ch11;
    tree->SetBranchAddress("adc_b4_ch11", &adc_b4_ch11);
    int tempEventNum = 0;

    for (int ievt = 0; ievt < tree->GetEntries(); ++ievt) {
        tree->GetEntry(ievt);

        std::vector<double> tempMC_nPE_PMT = {};
        // First, collect all the ADC channel variables into an array or vector
        std::vector<double *> adcChannels;
        if (adc_b1_ch1 > PE_CUT_EACH_PMT || adc_b1_ch1 < 0)
            continue;
        if (adc_b1_ch2 > PE_CUT_EACH_PMT || adc_b1_ch2 < 0)
            continue;
        if (adc_b1_ch3 > PE_CUT_EACH_PMT || adc_b1_ch3 < 0)
            continue;
        if (adc_b1_ch4 > PE_CUT_EACH_PMT || adc_b1_ch4 < 0)
            continue;
        if (adc_b1_ch5 > PE_CUT_EACH_PMT || adc_b1_ch5 < 0)
            continue;
        if (adc_b1_ch6 > PE_CUT_EACH_PMT || adc_b1_ch6 < 0)
            continue;
        if (adc_b1_ch7 > PE_CUT_EACH_PMT || adc_b1_ch7 < 0)
            continue;
        if (adc_b1_ch8 > PE_CUT_EACH_PMT || adc_b1_ch8 < 0)
            continue;
        if (adc_b1_ch9 > PE_CUT_EACH_PMT || adc_b1_ch9 < 0)
            continue;
        if (adc_b1_ch10 > PE_CUT_EACH_PMT || adc_b1_ch10 < 0)
            continue;
        if (adc_b1_ch11 > PE_CUT_EACH_PMT || adc_b1_ch11 < 0)
            continue;
        if (adc_b1_ch12 > PE_CUT_EACH_PMT || adc_b1_ch12 < 0)
            continue;
        if (adc_b1_ch13 > PE_CUT_EACH_PMT || adc_b1_ch13 < 0)
            continue;
        if (adc_b1_ch14 > PE_CUT_EACH_PMT || adc_b1_ch14 < 0)
            continue;
        if (adc_b1_ch15 > PE_CUT_EACH_PMT || adc_b1_ch15 < 0)
            continue;
        if (adc_b2_ch0 > PE_CUT_EACH_PMT || adc_b2_ch0 < 0)
            continue;
        if (adc_b2_ch1 > PE_CUT_EACH_PMT || adc_b2_ch1 < 0)
            continue;
        if (adc_b2_ch2 > PE_CUT_EACH_PMT || adc_b2_ch2 < 0)
            continue;
        if (adc_b2_ch3 > PE_CUT_EACH_PMT || adc_b2_ch3 < 0)
            continue;
        if (adc_b2_ch4 > PE_CUT_EACH_PMT || adc_b2_ch4 < 0)
            continue;
        if (adc_b2_ch5 > PE_CUT_EACH_PMT || adc_b2_ch5 < 0)
            continue;
        if (adc_b2_ch6 > PE_CUT_EACH_PMT || adc_b2_ch6 < 0)
            continue;
        if (adc_b2_ch7 > PE_CUT_EACH_PMT || adc_b2_ch7 < 0)
            continue;
        if (adc_b2_ch8 > PE_CUT_EACH_PMT || adc_b2_ch8 < 0)
            continue;
        if (adc_b2_ch9 > PE_CUT_EACH_PMT || adc_b2_ch9 < 0)
            continue;
        if (adc_b2_ch10 > PE_CUT_EACH_PMT || adc_b2_ch10 < 0)
            continue;
        if (adc_b2_ch11 > PE_CUT_EACH_PMT || adc_b2_ch11 < 0)
            continue;
        if (adc_b2_ch12 > PE_CUT_EACH_PMT || adc_b2_ch12 < 0)
            continue;
        if (adc_b2_ch13 > PE_CUT_EACH_PMT || adc_b2_ch13 < 0)
            continue;
        if (adc_b2_ch14 > PE_CUT_EACH_PMT || adc_b2_ch14 < 0)
            continue;
        if (adc_b3_ch0 > PE_CUT_EACH_PMT || adc_b3_ch0 < 0)
            continue;
        if (adc_b3_ch1 > PE_CUT_EACH_PMT || adc_b3_ch1 < 0)
            continue;
        if (adc_b3_ch2 > PE_CUT_EACH_PMT || adc_b3_ch2 < 0)
            continue;
        if (adc_b3_ch3 > PE_CUT_EACH_PMT || adc_b3_ch3 < 0)
            continue;
        if (adc_b3_ch4 > PE_CUT_EACH_PMT || adc_b3_ch4 < 0)
            continue;
        if (adc_b3_ch5 > PE_CUT_EACH_PMT || adc_b3_ch5 < 0)
            continue;
        if (adc_b3_ch6 > PE_CUT_EACH_PMT || adc_b3_ch6 < 0)
            continue;
        if (adc_b3_ch7 > PE_CUT_EACH_PMT || adc_b3_ch7 < 0)
            continue;
        if (adc_b3_ch8 > PE_CUT_EACH_PMT || adc_b3_ch8 < 0)
            continue;
        if (adc_b3_ch9 > PE_CUT_EACH_PMT || adc_b3_ch9 < 0)
            continue;
        if (adc_b3_ch10 > PE_CUT_EACH_PMT || adc_b3_ch10 < 0)
            continue;
        if (adc_b3_ch11 > PE_CUT_EACH_PMT || adc_b3_ch11 < 0)
            continue;
        if (adc_b3_ch12 > PE_CUT_EACH_PMT || adc_b3_ch12 < 0)
            continue;
        if (adc_b3_ch13 > PE_CUT_EACH_PMT || adc_b3_ch13 < 0)
            continue;
        if (adc_b3_ch14 > PE_CUT_EACH_PMT || adc_b3_ch14 < 0)
            continue;
        if (adc_b3_ch15 > PE_CUT_EACH_PMT || adc_b3_ch15 < 0)
            continue;
        if (adc_b4_ch0 > PE_CUT_EACH_PMT || adc_b4_ch0 < 0)
            continue;
        if (adc_b4_ch1 > PE_CUT_EACH_PMT || adc_b4_ch1 < 0)
            continue;
        if (adc_b4_ch2 > PE_CUT_EACH_PMT || adc_b4_ch2 < 0)
            continue;
        if (adc_b4_ch3 > PE_CUT_EACH_PMT || adc_b4_ch3 < 0)
            continue;
        if (adc_b4_ch4 > PE_CUT_EACH_PMT || adc_b4_ch4 < 0)
            continue;
        if (adc_b4_ch5 > PE_CUT_EACH_PMT || adc_b4_ch5 < 0)
            continue;
        if (adc_b4_ch6 > PE_CUT_EACH_PMT || adc_b4_ch6 < 0)
            continue;
        if (adc_b4_ch7 > PE_CUT_EACH_PMT || adc_b4_ch7 < 0)
            continue;
        if (adc_b4_ch8 > PE_CUT_EACH_PMT || adc_b4_ch8 < 0)
            continue;
        if (adc_b4_ch9 > PE_CUT_EACH_PMT || adc_b4_ch9 < 0)
            continue;
        if (adc_b4_ch10 > PE_CUT_EACH_PMT || adc_b4_ch10 < 0)
            continue;
        if (adc_b4_ch11 > PE_CUT_EACH_PMT || adc_b4_ch11 < 0)
            continue;
        adcChannels.push_back(&adc_b1_ch1);
        adcChannels.push_back(&adc_b1_ch2);
        adcChannels.push_back(&adc_b1_ch3);
        adcChannels.push_back(&adc_b1_ch4);
        adcChannels.push_back(&adc_b1_ch5);
        adcChannels.push_back(&adc_b1_ch6);
        adcChannels.push_back(&adc_b1_ch7);
        adcChannels.push_back(&adc_b1_ch8);
        adcChannels.push_back(&adc_b1_ch9);
        adcChannels.push_back(&adc_b1_ch10);
        adcChannels.push_back(&adc_b1_ch11);
        adcChannels.push_back(&adc_b1_ch12);
        adcChannels.push_back(&adc_b1_ch13);
        adcChannels.push_back(&adc_b1_ch14);
        adcChannels.push_back(&adc_b1_ch15);
        adcChannels.push_back(&adc_b2_ch0);
        adcChannels.push_back(&adc_b2_ch1);
        adcChannels.push_back(&adc_b2_ch2);
        adcChannels.push_back(&adc_b2_ch3);
        adcChannels.push_back(&adc_b2_ch4);
        adcChannels.push_back(&adc_b2_ch5);
        adcChannels.push_back(&adc_b2_ch6);
        adcChannels.push_back(&adc_b2_ch7);
        adcChannels.push_back(&adc_b2_ch8);
        adcChannels.push_back(&adc_b2_ch9);
        adcChannels.push_back(&adc_b2_ch10);
        adcChannels.push_back(&adc_b2_ch11);
        adcChannels.push_back(&adc_b2_ch12);
        adcChannels.push_back(&adc_b2_ch13);
        adcChannels.push_back(&adc_b2_ch14);
        adcChannels.push_back(&adc_b3_ch0);
        adcChannels.push_back(&adc_b3_ch1);
        adcChannels.push_back(&adc_b3_ch2);
        adcChannels.push_back(&adc_b3_ch3);
        adcChannels.push_back(&adc_b3_ch4);
        adcChannels.push_back(&adc_b3_ch5);
        adcChannels.push_back(&adc_b3_ch6);
        adcChannels.push_back(&adc_b3_ch7);
        adcChannels.push_back(&adc_b3_ch8);
        adcChannels.push_back(&adc_b3_ch9);
        adcChannels.push_back(&adc_b3_ch10);
        adcChannels.push_back(&adc_b3_ch11);
        adcChannels.push_back(&adc_b3_ch12);
        adcChannels.push_back(&adc_b3_ch13);
        adcChannels.push_back(&adc_b3_ch14);
        adcChannels.push_back(&adc_b3_ch15);

        // Now, use a loop to process up to TOTAL_PMT_COUNT channels
        double tempTotalNPE = 0;
        for (int i = 0; i < TOTAL_PMT_COUNT && i < adcChannels.size(); ++i) {
            tempMC_nPE_PMT.push_back(*adcChannels[i]);
            tempTotalNPE += *adcChannels[i];
        }

        if (tempTotalNPE < TOTAL_NPE_LOWER_CUT)
            continue;
        if (tempTotalNPE > TOTAL_NPE_UPPER_CUT)
            continue;
        Response_matrix.push_back(tempMC_nPE_PMT);
        tempEventNum++;
        if (tempEventNum == event_num)
            break;
    }

    return Response_matrix;
}
std::vector<double> GetTotalNPEVectorMC(TFile *inputFile, int event_num) {
    std::vector<double> total_npe = {};

    TTree *tree = (TTree *)inputFile->Get("output");

    int mcpdg;
    tree->SetBranchAddress("mcpdg", &mcpdg);
    double mcx;
    tree->SetBranchAddress("mcx", &mcx);
    double mcy;
    tree->SetBranchAddress("mcy", &mcy);
    double mcz;
    tree->SetBranchAddress("mcz", &mcz);
    double mcu;
    tree->SetBranchAddress("mcu", &mcu);
    double mcv;
    tree->SetBranchAddress("mcv", &mcv);
    double mcw;
    tree->SetBranchAddress("mcw", &mcw);
    double mcke;
    tree->SetBranchAddress("mcke", &mcke);
    double mct;
    tree->SetBranchAddress("mct", &mct);
    int evid;
    tree->SetBranchAddress("evid", &evid);
    int subev;
    tree->SetBranchAddress("subev", &subev);
    int nhits;
    tree->SetBranchAddress("nhits", &nhits);
    double triggerTime;
    tree->SetBranchAddress("triggerTime", &triggerTime);
    int mcparticlecount;
    tree->SetBranchAddress("mcparticlecount", &mcparticlecount);
    int mcpecount;
    tree->SetBranchAddress("mcpecount", &mcpecount);
    int mcnhits;
    tree->SetBranchAddress("mcnhits", &mcnhits);
    double scintEdep;
    tree->SetBranchAddress("scintEdep", &scintEdep);
    double scintEdepQuenched;
    tree->SetBranchAddress("scintEdepQuenched", &scintEdepQuenched);
    double scintPhotons;
    tree->SetBranchAddress("scintPhotons", &scintPhotons);
    double remPhotons;
    tree->SetBranchAddress("remPhotons", &remPhotons);
    double cherPhotons;
    tree->SetBranchAddress("cherPhotons", &cherPhotons);
    std::vector<int> *mcpdgs = {};
    tree->SetBranchAddress("mcpdgs", &mcpdgs);
    std::vector<double> *mcxs = {};
    tree->SetBranchAddress("mcxs", &mcxs);
    std::vector<double> *mcys = {};
    tree->SetBranchAddress("mcys", &mcys);
    std::vector<double> *mczs = {};
    tree->SetBranchAddress("mczs", &mczs);
    std::vector<double> *mcus = {};
    tree->SetBranchAddress("mcus", &mcus);
    std::vector<double> *mcvs = {};
    tree->SetBranchAddress("mcvs", &mcvs);
    std::vector<double> *mcws = {};
    tree->SetBranchAddress("mcws", &mcws);
    std::vector<double> *mckes = {};
    tree->SetBranchAddress("mckes", &mckes);
    std::vector<double> *mcts = {};
    tree->SetBranchAddress("mcts", &mcts);
    std::vector<double> *hitPMTTime = {};
    tree->SetBranchAddress("hitPMTTime", &hitPMTTime);
    std::vector<double> *hitPMTCharge = {};
    tree->SetBranchAddress("hitPMTCharge", &hitPMTCharge);
    std::vector<double> *hitPMTDigitizedTime = {};
    tree->SetBranchAddress("hitPMTDigitizedTime", &hitPMTDigitizedTime);
    std::vector<double> *hitPMTDigitizedCharge = {};
    tree->SetBranchAddress("hitPMTDigitizedCharge", &hitPMTDigitizedCharge);
    std::vector<int> *hitPMTNCrossings = {};
    tree->SetBranchAddress("hitPMTNCrossings", &hitPMTNCrossings);
    std::vector<int> *mcPMTID = {};
    tree->SetBranchAddress("mcPMTID", &mcPMTID);
    std::vector<int> *mcPMTNPE = {};
    tree->SetBranchAddress("mcPMTNPE", &mcPMTNPE);
    std::vector<double> *mcPETime = {};
    tree->SetBranchAddress("mcPETime", &mcPETime);
    std::vector<int> *mcPEProcess = {};
    tree->SetBranchAddress("mcPEProcess", &mcPEProcess);
    std::vector<double> *mcPEWavelength = {};
    tree->SetBranchAddress("mcPEWavelength", &mcPEWavelength);
    std::vector<double> *mcPEx = {};
    tree->SetBranchAddress("mcPEx", &mcPEx);
    std::vector<double> *mcPEy = {};
    tree->SetBranchAddress("mcPEy", &mcPEy);
    std::vector<double> *mcPEz = {};
    tree->SetBranchAddress("mcPEz", &mcPEz);
    std::vector<int> *hitPMTID = {};
    tree->SetBranchAddress("hitPMTID", &hitPMTID);

    for (int ievt = 0; ievt < event_num; ++ievt) {
        tree->GetEntry(ievt);

        int tempTotalNPE = 0;
        for (int j = 0; j < TOTAL_PMT_COUNT; ++j) {
            tempTotalNPE += mcPMTNPE->at(j);
        }
        total_npe.push_back(tempTotalNPE);
    }
    return total_npe;
}

std::vector<double> GetTotalNPEVectorMC2(TFile *inputFile, int event_num) {
    std::vector<double> total_npe = {};

    TTree *tree = (TTree *)inputFile->Get("tree");

    int PMTPE[46];
    tree->SetBranchAddress("PMTPE", &PMTPE);

    int tempEventNum = 0;
    for (int ievt = 0; ievt < tree->GetEntries(); ++ievt) {
        tree->GetEntry(ievt);

        int tempTotalNPE = 0;
        for (int j = 0; j < TOTAL_PMT_COUNT; ++j) {
            tempTotalNPE += PMTPE[j];
        }
        histTotalNPEMC->Fill(tempTotalNPE);
        // if (tempTotalNPE < TOTAL_NPE_LOWER_CUT+100) continue;
        // if (tempTotalNPE > TOTAL_NPE_UPPER_CUT+100) continue;
        total_npe.push_back(tempTotalNPE);
        tempEventNum++;
        if (tempEventNum == event_num)
            break;
    }
    return total_npe;
}

std::vector<double> GetTotalNPEVectorMC3(TFile *inputFile, int event_num) {
    std::vector<double> total_npe = {};

    TTree *tree = (TTree *)inputFile->Get("output");
    if (!tree) {
        std::cerr << "Error: TTree 'output' not found in the input file."
                  << std::endl;
    } else {
        std::cout << "Successfully retrieved TTree 'output'." << std::endl;
    }

    std::vector<int> *mcPMTNPE = {};
    if (tree->SetBranchAddress("mcPMTNPE", &mcPMTNPE) != 0) {
        std::cerr << "Error: Could not set branch address for 'mcPMTNPE'."
                  << std::endl;
    } else {
        std::cout << "Successfully set branch address for 'mcPMTNPE'."
                  << std::endl;
    }
    std::vector<int> *mcPMTID = {};
    if (tree->SetBranchAddress("mcPMTID", &mcPMTID) != 0) {
        std::cerr << "Error: Could not set branch address for 'mcPMTID'."
                  << std::endl;
    } else {
        std::cout << "Successfully set branch address for 'mcPMTID'."
                  << std::endl;
    }

    int tempEventNum = 0;
    Long64_t totalEntries = tree->GetEntries();
    std::cout << "Total entries in tree: " << totalEntries << std::endl;

    for (Long64_t ievt = 0; ievt < totalEntries; ++ievt) {
        std::cout << "Processing event " << ievt << std::endl;
        Long64_t nb = tree->GetEntry(ievt);
        if (nb <= 0) {
            std::cerr << "Warning: No bytes read for entry " << ievt
                      << std::endl;
            continue;
        }

        int tempTotalNPE = 0;
        std::vector<double> tempMC_nPE_PMT(TOTAL_PMT_COUNT);
        for (int j = 0; j < TOTAL_PMT_COUNT; ++j) {
            tempTotalNPE += mcPMTNPE->at(j);
        }
        total_npe.push_back(tempTotalNPE);
        histTotalNPEMC->Fill(tempTotalNPE);

        if (tempEventNum == event_num)
            break;
    }
    return total_npe;
}

std::vector<double> GetTotalNPEVectorData(TFile *inputFile, int event_num) {
    std::vector<double> total_npe = {};
    int tempEventNum = 0;

    TTree *tree = (TTree *)inputFile->Get("tree");
    double adc_b1_ch1;
    tree->SetBranchAddress("adc_b1_ch1", &adc_b1_ch1);
    double adc_b1_ch2;
    tree->SetBranchAddress("adc_b1_ch2", &adc_b1_ch2);
    double adc_b1_ch3;
    tree->SetBranchAddress("adc_b1_ch3", &adc_b1_ch3);
    double adc_b1_ch4;
    tree->SetBranchAddress("adc_b1_ch4", &adc_b1_ch4);
    double adc_b1_ch5;
    tree->SetBranchAddress("adc_b1_ch5", &adc_b1_ch5);
    double adc_b1_ch6;
    tree->SetBranchAddress("adc_b1_ch6", &adc_b1_ch6);
    double adc_b1_ch7;
    tree->SetBranchAddress("adc_b1_ch7", &adc_b1_ch7);
    double adc_b1_ch8;
    tree->SetBranchAddress("adc_b1_ch8", &adc_b1_ch8);
    double adc_b1_ch9;
    tree->SetBranchAddress("adc_b1_ch9", &adc_b1_ch9);
    double adc_b1_ch10;
    tree->SetBranchAddress("adc_b1_ch10", &adc_b1_ch10);
    double adc_b1_ch11;
    tree->SetBranchAddress("adc_b1_ch11", &adc_b1_ch11);
    double adc_b1_ch12;
    tree->SetBranchAddress("adc_b1_ch12", &adc_b1_ch12);
    double adc_b1_ch13;
    tree->SetBranchAddress("adc_b1_ch13", &adc_b1_ch13);
    double adc_b1_ch14;
    tree->SetBranchAddress("adc_b1_ch14", &adc_b1_ch14);
    double adc_b1_ch15;
    tree->SetBranchAddress("adc_b1_ch15", &adc_b1_ch15);
    double adc_b2_ch0;
    tree->SetBranchAddress("adc_b2_ch0", &adc_b2_ch0);
    double adc_b2_ch1;
    tree->SetBranchAddress("adc_b2_ch1", &adc_b2_ch1);
    double adc_b2_ch2;
    tree->SetBranchAddress("adc_b2_ch2", &adc_b2_ch2);
    double adc_b2_ch3;
    tree->SetBranchAddress("adc_b2_ch3", &adc_b2_ch3);
    double adc_b2_ch4;
    tree->SetBranchAddress("adc_b2_ch4", &adc_b2_ch4);
    double adc_b2_ch5;
    tree->SetBranchAddress("adc_b2_ch5", &adc_b2_ch5);
    double adc_b2_ch6;
    tree->SetBranchAddress("adc_b2_ch6", &adc_b2_ch6);
    double adc_b2_ch7;
    tree->SetBranchAddress("adc_b2_ch7", &adc_b2_ch7);
    double adc_b2_ch8;
    tree->SetBranchAddress("adc_b2_ch8", &adc_b2_ch8);
    double adc_b2_ch9;
    tree->SetBranchAddress("adc_b2_ch9", &adc_b2_ch9);
    double adc_b2_ch10;
    tree->SetBranchAddress("adc_b2_ch10", &adc_b2_ch10);
    double adc_b2_ch11;
    tree->SetBranchAddress("adc_b2_ch11", &adc_b2_ch11);
    double adc_b2_ch12;
    tree->SetBranchAddress("adc_b2_ch12", &adc_b2_ch12);
    double adc_b2_ch13;
    tree->SetBranchAddress("adc_b2_ch13", &adc_b2_ch13);
    double adc_b2_ch14;
    tree->SetBranchAddress("adc_b2_ch14", &adc_b2_ch14);
    double adc_b3_ch0;
    tree->SetBranchAddress("adc_b3_ch0", &adc_b3_ch0);
    double adc_b3_ch1;
    tree->SetBranchAddress("adc_b3_ch1", &adc_b3_ch1);
    double adc_b3_ch2;
    tree->SetBranchAddress("adc_b3_ch2", &adc_b3_ch2);
    double adc_b3_ch3;
    tree->SetBranchAddress("adc_b3_ch3", &adc_b3_ch3);
    double adc_b3_ch4;
    tree->SetBranchAddress("adc_b3_ch4", &adc_b3_ch4);
    double adc_b3_ch5;
    tree->SetBranchAddress("adc_b3_ch5", &adc_b3_ch5);
    double adc_b3_ch6;
    tree->SetBranchAddress("adc_b3_ch6", &adc_b3_ch6);
    double adc_b3_ch7;
    tree->SetBranchAddress("adc_b3_ch7", &adc_b3_ch7);
    double adc_b3_ch8;
    tree->SetBranchAddress("adc_b3_ch8", &adc_b3_ch8);
    double adc_b3_ch9;
    tree->SetBranchAddress("adc_b3_ch9", &adc_b3_ch9);
    double adc_b3_ch10;
    tree->SetBranchAddress("adc_b3_ch10", &adc_b3_ch10);
    double adc_b3_ch11;
    tree->SetBranchAddress("adc_b3_ch11", &adc_b3_ch11);
    double adc_b3_ch12;
    tree->SetBranchAddress("adc_b3_ch12", &adc_b3_ch12);
    double adc_b3_ch13;
    tree->SetBranchAddress("adc_b3_ch13", &adc_b3_ch13);
    double adc_b3_ch14;
    tree->SetBranchAddress("adc_b3_ch14", &adc_b3_ch14);
    double adc_b3_ch15;
    tree->SetBranchAddress("adc_b3_ch15", &adc_b3_ch15);
    double adc_b4_ch0;
    tree->SetBranchAddress("adc_b4_ch0", &adc_b4_ch0);
    double adc_b4_ch1;
    tree->SetBranchAddress("adc_b4_ch1", &adc_b4_ch1);
    double adc_b4_ch2;
    tree->SetBranchAddress("adc_b4_ch2", &adc_b4_ch2);
    double adc_b4_ch3;
    tree->SetBranchAddress("adc_b4_ch3", &adc_b4_ch3);
    double adc_b4_ch4;
    tree->SetBranchAddress("adc_b4_ch4", &adc_b4_ch4);
    double adc_b4_ch5;
    tree->SetBranchAddress("adc_b4_ch5", &adc_b4_ch5);
    double adc_b4_ch6;
    tree->SetBranchAddress("adc_b4_ch6", &adc_b4_ch6);
    double adc_b4_ch7;
    tree->SetBranchAddress("adc_b4_ch7", &adc_b4_ch7);
    double adc_b4_ch8;
    tree->SetBranchAddress("adc_b4_ch8", &adc_b4_ch8);
    double adc_b4_ch9;
    tree->SetBranchAddress("adc_b4_ch9", &adc_b4_ch9);
    double adc_b4_ch10;
    tree->SetBranchAddress("adc_b4_ch10", &adc_b4_ch10);
    double adc_b4_ch11;
    tree->SetBranchAddress("adc_b4_ch11", &adc_b4_ch11);

    for (int ievt = 0; ievt < tree->GetEntries(); ++ievt) {
        tree->GetEntry(ievt);

        // First, collect all the ADC channel variables into an array or vector
        if (adc_b1_ch1 > PE_CUT_EACH_PMT || adc_b1_ch1 < 0)
            continue;
        if (adc_b1_ch2 > PE_CUT_EACH_PMT || adc_b1_ch2 < 0)
            continue;
        if (adc_b1_ch3 > PE_CUT_EACH_PMT || adc_b1_ch3 < 0)
            continue;
        if (adc_b1_ch4 > PE_CUT_EACH_PMT || adc_b1_ch4 < 0)
            continue;
        if (adc_b1_ch5 > PE_CUT_EACH_PMT || adc_b1_ch5 < 0)
            continue;
        if (adc_b1_ch6 > PE_CUT_EACH_PMT || adc_b1_ch6 < 0)
            continue;
        if (adc_b1_ch7 > PE_CUT_EACH_PMT || adc_b1_ch7 < 0)
            continue;
        if (adc_b1_ch8 > PE_CUT_EACH_PMT || adc_b1_ch8 < 0)
            continue;
        if (adc_b1_ch9 > PE_CUT_EACH_PMT || adc_b1_ch9 < 0)
            continue;
        if (adc_b1_ch10 > PE_CUT_EACH_PMT || adc_b1_ch10 < 0)
            continue;
        if (adc_b1_ch11 > PE_CUT_EACH_PMT || adc_b1_ch11 < 0)
            continue;
        if (adc_b1_ch12 > PE_CUT_EACH_PMT || adc_b1_ch12 < 0)
            continue;
        if (adc_b1_ch13 > PE_CUT_EACH_PMT || adc_b1_ch13 < 0)
            continue;
        if (adc_b1_ch14 > PE_CUT_EACH_PMT || adc_b1_ch14 < 0)
            continue;
        if (adc_b1_ch15 > PE_CUT_EACH_PMT || adc_b1_ch15 < 0)
            continue;
        if (adc_b2_ch0 > PE_CUT_EACH_PMT || adc_b2_ch0 < 0)
            continue;
        if (adc_b2_ch1 > PE_CUT_EACH_PMT || adc_b2_ch1 < 0)
            continue;
        if (adc_b2_ch2 > PE_CUT_EACH_PMT || adc_b2_ch2 < 0)
            continue;
        if (adc_b2_ch3 > PE_CUT_EACH_PMT || adc_b2_ch3 < 0)
            continue;
        if (adc_b2_ch4 > PE_CUT_EACH_PMT || adc_b2_ch4 < 0)
            continue;
        if (adc_b2_ch5 > PE_CUT_EACH_PMT || adc_b2_ch5 < 0)
            continue;
        if (adc_b2_ch6 > PE_CUT_EACH_PMT || adc_b2_ch6 < 0)
            continue;
        if (adc_b2_ch7 > PE_CUT_EACH_PMT || adc_b2_ch7 < 0)
            continue;
        if (adc_b2_ch8 > PE_CUT_EACH_PMT || adc_b2_ch8 < 0)
            continue;
        if (adc_b2_ch9 > PE_CUT_EACH_PMT || adc_b2_ch9 < 0)
            continue;
        if (adc_b2_ch10 > PE_CUT_EACH_PMT || adc_b2_ch10 < 0)
            continue;
        if (adc_b2_ch11 > PE_CUT_EACH_PMT || adc_b2_ch11 < 0)
            continue;
        if (adc_b2_ch12 > PE_CUT_EACH_PMT || adc_b2_ch12 < 0)
            continue;
        if (adc_b2_ch13 > PE_CUT_EACH_PMT || adc_b2_ch13 < 0)
            continue;
        if (adc_b2_ch14 > PE_CUT_EACH_PMT || adc_b2_ch14 < 0)
            continue;
        if (adc_b3_ch0 > PE_CUT_EACH_PMT || adc_b3_ch0 < 0)
            continue;
        if (adc_b3_ch1 > PE_CUT_EACH_PMT || adc_b3_ch1 < 0)
            continue;
        if (adc_b3_ch2 > PE_CUT_EACH_PMT || adc_b3_ch2 < 0)
            continue;
        if (adc_b3_ch3 > PE_CUT_EACH_PMT || adc_b3_ch3 < 0)
            continue;
        if (adc_b3_ch4 > PE_CUT_EACH_PMT || adc_b3_ch4 < 0)
            continue;
        if (adc_b3_ch5 > PE_CUT_EACH_PMT || adc_b3_ch5 < 0)
            continue;
        if (adc_b3_ch6 > PE_CUT_EACH_PMT || adc_b3_ch6 < 0)
            continue;
        if (adc_b3_ch7 > PE_CUT_EACH_PMT || adc_b3_ch7 < 0)
            continue;
        if (adc_b3_ch8 > PE_CUT_EACH_PMT || adc_b3_ch8 < 0)
            continue;
        if (adc_b3_ch9 > PE_CUT_EACH_PMT || adc_b3_ch9 < 0)
            continue;
        if (adc_b3_ch10 > PE_CUT_EACH_PMT || adc_b3_ch10 < 0)
            continue;
        if (adc_b3_ch11 > PE_CUT_EACH_PMT || adc_b3_ch11 < 0)
            continue;
        if (adc_b3_ch12 > PE_CUT_EACH_PMT || adc_b3_ch12 < 0)
            continue;
        if (adc_b3_ch13 > PE_CUT_EACH_PMT || adc_b3_ch13 < 0)
            continue;
        if (adc_b3_ch14 > PE_CUT_EACH_PMT || adc_b3_ch14 < 0)
            continue;
        if (adc_b3_ch15 > PE_CUT_EACH_PMT || adc_b3_ch15 < 0)
            continue;
        if (adc_b4_ch0 > PE_CUT_EACH_PMT || adc_b4_ch0 < 0)
            continue;
        if (adc_b4_ch1 > PE_CUT_EACH_PMT || adc_b4_ch1 < 0)
            continue;
        if (adc_b4_ch2 > PE_CUT_EACH_PMT || adc_b4_ch2 < 0)
            continue;
        if (adc_b4_ch3 > PE_CUT_EACH_PMT || adc_b4_ch3 < 0)
            continue;
        if (adc_b4_ch4 > PE_CUT_EACH_PMT || adc_b4_ch4 < 0)
            continue;
        if (adc_b4_ch5 > PE_CUT_EACH_PMT || adc_b4_ch5 < 0)
            continue;
        if (adc_b4_ch6 > PE_CUT_EACH_PMT || adc_b4_ch6 < 0)
            continue;
        if (adc_b4_ch7 > PE_CUT_EACH_PMT || adc_b4_ch7 < 0)
            continue;
        if (adc_b4_ch8 > PE_CUT_EACH_PMT || adc_b4_ch8 < 0)
            continue;
        if (adc_b4_ch9 > PE_CUT_EACH_PMT || adc_b4_ch9 < 0)
            continue;
        if (adc_b4_ch10 > PE_CUT_EACH_PMT || adc_b4_ch10 < 0)
            continue;
        if (adc_b4_ch11 > PE_CUT_EACH_PMT || adc_b4_ch11 < 0)
            continue;
        std::vector<double *> adcChannels;
        adcChannels.push_back(&adc_b1_ch1);
        adcChannels.push_back(&adc_b1_ch2);
        adcChannels.push_back(&adc_b1_ch3);
        adcChannels.push_back(&adc_b1_ch4);
        adcChannels.push_back(&adc_b1_ch5);
        adcChannels.push_back(&adc_b1_ch6);
        adcChannels.push_back(&adc_b1_ch7);
        adcChannels.push_back(&adc_b1_ch8);
        adcChannels.push_back(&adc_b1_ch9);
        adcChannels.push_back(&adc_b1_ch10);
        adcChannels.push_back(&adc_b1_ch11);
        adcChannels.push_back(&adc_b1_ch12);
        adcChannels.push_back(&adc_b1_ch13);
        adcChannels.push_back(&adc_b1_ch14);
        adcChannels.push_back(&adc_b1_ch15);
        adcChannels.push_back(&adc_b2_ch0);
        adcChannels.push_back(&adc_b2_ch1);
        adcChannels.push_back(&adc_b2_ch2);
        adcChannels.push_back(&adc_b2_ch3);
        adcChannels.push_back(&adc_b2_ch4);
        adcChannels.push_back(&adc_b2_ch5);
        adcChannels.push_back(&adc_b2_ch6);
        adcChannels.push_back(&adc_b2_ch7);
        adcChannels.push_back(&adc_b2_ch8);
        adcChannels.push_back(&adc_b2_ch9);
        adcChannels.push_back(&adc_b2_ch10);
        adcChannels.push_back(&adc_b2_ch11);
        adcChannels.push_back(&adc_b2_ch12);
        adcChannels.push_back(&adc_b2_ch13);
        adcChannels.push_back(&adc_b2_ch14);
        adcChannels.push_back(&adc_b3_ch0);
        adcChannels.push_back(&adc_b3_ch1);
        adcChannels.push_back(&adc_b3_ch2);
        adcChannels.push_back(&adc_b3_ch3);
        adcChannels.push_back(&adc_b3_ch4);
        adcChannels.push_back(&adc_b3_ch5);
        adcChannels.push_back(&adc_b3_ch6);
        adcChannels.push_back(&adc_b3_ch7);
        adcChannels.push_back(&adc_b3_ch8);
        adcChannels.push_back(&adc_b3_ch9);
        adcChannels.push_back(&adc_b3_ch10);
        adcChannels.push_back(&adc_b3_ch11);
        adcChannels.push_back(&adc_b3_ch12);
        adcChannels.push_back(&adc_b3_ch13);
        adcChannels.push_back(&adc_b3_ch14);
        adcChannels.push_back(&adc_b3_ch15);

        // Now, use a loop to process up to TOTAL_PMT_COUNT channels
        double tempTotalNPE = 0;
        for (int i = 0; i < TOTAL_PMT_COUNT && i < adcChannels.size(); ++i) {
            tempTotalNPE += *adcChannels[i];
        }

        histTotalNPEData->Fill(tempTotalNPE);
        if (tempTotalNPE < TOTAL_NPE_LOWER_CUT)
            continue;
        if (tempTotalNPE > TOTAL_NPE_UPPER_CUT)
            continue;
        total_npe.push_back(tempTotalNPE);
        tempEventNum++;
        if (tempEventNum == event_num)
            break;
    }

    return total_npe;
}

void visualizeMatrix(const std::vector<std::vector<double>> &matrix,
                     const std::string &title = "Matrix Visualization") {
    // 행렬의 크기 확인
    int nRows = matrix.size();
    if (nRows == 0) {
        std::cerr << "Error: Matrix has no rows." << std::endl;
        return;
    }
    int nCols = matrix[0].size();
    for (const auto &row : matrix) {
        if (row.size() != nCols) {
            std::cerr << "Error: All rows must have the same number of columns."
                      << std::endl;
            return;
        }
    }

    // ROOT 캔버스 생성
    TCanvas *c = new TCanvas("c", title.c_str(), 800, 600);

    // 2D 히스토그램 생성
    TH2D *h = new TH2D("h", title.c_str(), nCols, 0, nCols, nRows, 0, nRows);

    // 히스토그램에 데이터 채우기
    for (int i = 0; i < nRows; ++i) {
        for (int j = 0; j < nCols; ++j) {
            h->SetBinContent(j + 1, nRows - i, matrix[i][j]);
        }
    }

    // 스타일 설정
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRainBow);
    gStyle->SetPaintTextFormat("4.1f");

    // 히스토그램 그리기
    h->GetXaxis()->SetTitle("Columns");
    h->GetYaxis()->SetTitle("Rows");
    h->Draw("COLZ TEXT");

    // 캔버스 업데이트
    c->Update();
}

int main(int argc, char **argv) {
    TApplication app("app", &argc, argv);
    // Check if enough arguments are provided
    if (argc >= 3)
        TOTAL_PMT_COUNT = std::atoi(argv[2]);
    if (argc >= 4)
        TOTAL_NPE_LOWER_CUT = std::atoi(argv[3]);
    if (argc >= 5)
        TOTAL_NPE_UPPER_CUT = std::atoi(argv[4]);
    // TFile* fileMC = new TFile("MC_1ton_1GeV_downward.root");
    // TFile* fileMC = new TFile("result_water.root");
    // TFile* fileData = new TFile("combined_data_full_1.root");
    TFile *fileMC = new TFile("mc_result_wbls03.root");
    TFile *fileData = new TFile("combined_data_025_WbLS.root");

    histTotalNPEData =
        new TH1D("histTotalNPEData", "histTotalNPEData", 50, 0, 1000);
    histTotalNPEMC = new TH1D("histTotalNPEMC", "histTotalNPEMC", 50, 0, 1000);

    int event_num = std::atoi(argv[1]);

    std::cout << "event_num: " << event_num << std::endl;

    std::vector<std::vector<double>> MC_Response_matrix;
    try {
        // MakeResponseMatrixMC2(fileMC, event_num, MC_Response_matrix);
        MakeResponseMatrixMC3(fileMC, event_num, MC_Response_matrix);
    } catch (const std::exception &e) {
        std::cerr << "Exception caught: " << e.what() << std::endl;
        return 1;
    }
    std::cout << "MC_Response_matrix ready: " << std::endl;

    std::vector<std::vector<double>> Data_Response_matrix =
        MakeResponseMatrixData(fileData, event_num);
    std::cout << "Data_Response_matrix ready: " << std::endl;
    // std::vector<double> MC_nPE_event = GetTotalNPEVectorMC(fileMC,
    // event_num);
    std::vector<double> MC_nPE_event = GetTotalNPEVectorMC2(fileMC, event_num);
    // std::vector<double> MC_nPE_event = GetTotalNPEVectorMC3(fileMC,
    // event_num); std::vector<double> MC_nPE_event =
    // GetTotalNPEVectorData(fileData, event_num);
    std::cout << "MC_nPE_event ready: " << std::endl;
    std::vector<double> Data_nPE_event =
        GetTotalNPEVectorData(fileData, event_num);
    std::cout << "Data_nPE_event ready: " << std::endl;

    std::vector<double> scale_event = {};

    std::cout << "MC_Response_matrix.size(): " << MC_Response_matrix.size()
              << std::endl;
    printMatrix(Data_Response_matrix);
    std::cout << "Data_Response_matrix.size(): " << Data_Response_matrix.size()
              << std::endl;
    std::cout << "Data_nPE_event.size(): " << Data_nPE_event.size()
              << std::endl;
    std::cout << "MC_nPE_event.size(): " << MC_nPE_event.size() << std::endl;
    printVector(MC_nPE_event);
    // visualizeMatrix(Data_Response_matrix);
    // visualizeMatrix(MC_Response_matrix);
    // app.Run(kTRUE);

    // --- 이차 계획법 문제 풀이 시작 ---

    // Eigen 타입으로 변환
    int m = Data_Response_matrix.size();
    int n = Data_Response_matrix[0].size();

    Eigen::MatrixXd Data(m, n);
    Eigen::VectorXd MC(m);

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            Data(i, j) = Data_Response_matrix[i][j];
        }
    }

    for (int i = 0; i < m; ++i) {
        MC(i) = MC_nPE_event[i];
    }

    // QP 문제 설정 및 풀이
    Eigen::MatrixXd Q = Data.transpose() * Data;
    Eigen::VectorXd c = -Data.transpose() * MC;

    real_t *H = Q.data();
    real_t *g = c.data();

    Eigen::VectorXd lb = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd ub = Eigen::VectorXd::Constant(n, 1e20);

    real_t *lb_ptr = lb.data();
    real_t *ub_ptr = ub.data();

    QProblem qp(n, 0);

    Options options;
    options.printLevel = PL_HIGH;
    qp.setOptions(options);

    int nWSR = 100;

    int ret = qp.init(H, g, nullptr, lb_ptr, ub_ptr, nullptr, nullptr, nWSR);

    Eigen::VectorXd V(n);
    if (ret == qpOASES::SUCCESSFUL_RETURN) {
        qp.getPrimalSolution(V.data());

        std::cout << "Optimal V:\n" << V << std::endl;
    } else {
        std::cout << "QP failed to solve. Return code: " << ret << std::endl;
    }
    std::vector<double> F(n);
    for (int i = 0; i < n; ++i) {
        F[i] = V(i);
    }

    /*
    int PMTID_to_fix_1 = -1;
    PMTID_to_fix_1 = std::stoi(argv[1]);

    if (PMTID_to_fix_1 < 0 || PMTID_to_fix_1 > 45) {
        std::cout << "invalid PMT ID" << std::endl;
        return 1;
    }

    for (int i = 0; i < event_num; ++i) {
        double tempScale =
    MC_Response_matrix.at(i).at(PMTID_to_fix_1)/Data_Response_matrix.at(i).at(PMTID_to_fix_1);
        scale_event.push_back(tempScale);
    };

    std::vector<std::vector<double>> Scaled_Data_Response_matrix = {};
    std::vector<double> Scaled_Data_nPE_event;
    int event_index = 0;

    for (const auto& v : Data_Response_matrix) {
        std::vector<double> scaled_v;
        // Loop over each element in the vector
        for (const auto& e : v) {
            // Scale the element and add it to the new vector
            double scaled_e = e * scale_event.at(event_index);
            scaled_v.push_back(scaled_e);
        }
        // Add the scaled vector to the new matrix
        Scaled_Data_Response_matrix.push_back(scaled_v);
        event_index++;
    }
    //std::vector<double> F = solveForF(Data_Response_matrix, MC_nPE_event);

    std::vector<double> vec_to_substract = {};
    for (auto& vec : Scaled_Data_Response_matrix) {
        if (PMTID_to_fix_1 >= 0 && PMTID_to_fix_1 < vec.size()) {
            vec_to_substract.push_back(vec.at(PMTID_to_fix_1));
            vec.erase(vec.begin() + PMTID_to_fix_1);
        } else {
            std::cerr << "인덱스가 범위를 벗어났습니다." << std::endl;
        }
    }

    std::vector<double> substracted_NPE_vector;
    for (size_t i = 0; i < MC_nPE_event.size(); ++i) {
        double tempNPE = MC_nPE_event.at(i) - vec_to_substract.at(i);
        substracted_NPE_vector.push_back(tempNPE);
    }


    //std::vector<double> F = solveForF(Data_Response_matrix, Data_nPE_event,
    F_min, F_max);
    //

    //std::vector<double> F = solveForF(Scaled_Data_Response_matrix,
    substracted_NPE_vector);
    // Define the bounds for F
    //int m = Scaled_Data_Response_matrix[0].size();
    int m = Data_Response_matrix[0].size();
    std::vector<double> F_min(m, 0.0);    // Lower bounds (e.g., non-negative)
    std::vector<double> F_max(m, 1000.0);    // Upper bounds

    std::vector<double> F = {};
    try {
        //F = solveForF(Scaled_Data_Response_matrix, substracted_NPE_vector,
    F_min, F_max); F = solveForF(Data_Response_matrix, Data_nPE_event, F_min,
    F_max);

        // Output the results
        for (size_t i = 0; i < F.size(); ++i) {
            std::cout << "F[" << i << "] = " << F[i] << std::endl;
        }
    } catch (const std::exception& ex) {
        std::cerr << "Error solving for F: " << ex.what() << std::endl;
    }
    F.insert(F.begin() + PMTID_to_fix_1, 1);

    std::vector<std::vector<double>> vec_F;
    vec_F.push_back(F);
    visualizeMatrix(vec_F);
    */
    std::vector<std::vector<double>> vec_F;
    vec_F.push_back(F);
    visualizeMatrix(vec_F);

    TCanvas *c1 = new TCanvas;
    histTotalNPEData->Draw();
    // Get the y-axis range to draw the lines from the bottom to the top of the
    // histogram
    Double_t ymin = histTotalNPEData->GetMinimum();
    Double_t ymax = histTotalNPEData->GetMaximum();

    // Create a dashed line at TOTAL_NPE_LOWER_CUT
    TLine *lineLower =
        new TLine(TOTAL_NPE_LOWER_CUT, ymin, TOTAL_NPE_LOWER_CUT, ymax);
    lineLower->SetLineColor(kRed); // Set line color (you can choose any color)
    lineLower->SetLineStyle(2);    // Set line style to dashed
    lineLower->SetLineWidth(2);    // Set line style to dashed
    lineLower->Draw("same");       // Draw on the same canvas

    // Create a dashed line at TOTAL_NPE_UPPER_CUT
    TLine *lineUpper =
        new TLine(TOTAL_NPE_UPPER_CUT, ymin, TOTAL_NPE_UPPER_CUT, ymax);
    lineUpper->SetLineColor(kRed); // Set line color (you can choose any color)
    lineUpper->SetLineStyle(2);    // Set line style to dashed
    lineUpper->SetLineWidth(2);    // Set line style to dashed
    lineUpper->Draw("same");
    c1->Update();

    TCanvas *c2 = new TCanvas;
    histTotalNPEMC->Draw();
    c2->Update();

    // F_j 값 출력
    std::cout << "F vector: " << std::endl;
    for (double f : F) {
        std::cout << f << " ";
    }
    std::cout << std::endl;

    app.Run();
    return 0;
}
