void CreateLegend(TLegend *leg, TObject *obj, const char *label,
                  const char *option = "l") {
    leg->AddEntry(obj, label, option);
}
double GetTruncatedMean(TH1D* hist, int xmin, int xmax) {
    double sum = 0;
    double weightedSum = 0;

    int binMin = hist->FindBin(xmin);
    int binMax = hist->FindBin(xmax);

    for (int i = binMin; i <= binMax; ++i) {
        double binContent = hist->GetBinContent(i);
        double binCenter = hist->GetBinCenter(i);

        sum += binContent;            // 총 bin 내용 합
        weightedSum += binContent * binCenter; // bin 값에 가중치를 적용한 합
    }

    // 평균 계산
    double truncatedMean = sum > 0 ? weightedSum / sum : 0;

    std::cout << "Truncated Mean: " << truncatedMean << std::endl;

    return truncatedMean;
}
void DrawPEDistribution() {
    // TFile* data = new
    // TFile("/Users/gwon/WbLS/1ton_analysis/BNL1tonAnalysis/build/combined_data_025_WbLS.root");
    TFile *data = new TFile("/Users/gwon/WbLS/1ton_analysis/BNL1tonAnalysis/"
                            "build/combined_data_060_WbLS.root");
     //TFile* data = new
     //TFile("/Users/gwon/WbLS/1ton_analysis/BNL1tonAnalysis/build/combined_data_full_1.root");
    // TFile* mc = new
    // TFile("/Users/gwon/WbLS/1ton_analysis/BNL1tonAnalysis/build/mc_03_wbls.root");
    TFile *mc = new TFile(
        "/Users/gwon/WbLS/1ton_analysis/BNL1tonAnalysis/build/mc_05_wbls.root");
    TFile *mc_noLY = new TFile("/Users/gwon/WbLS/1ton_analysis/BNL1tonAnalysis/"
                               "build/mc_05_wbls_noLY.root");
    // TFile* mc_noLY = new
    // TFile("/Users/gwon/WbLS/1ton_analysis/BNL1tonAnalysis/build/mc_04_wbls_noLY.root");
    TFile *mc_1000 = new TFile("/Users/gwon/WbLS/1ton_analysis/BNL1tonAnalysis/"
                               "build/mc_03_wbls_1000.root");
    // TFile* mc_1000 = new
    // TFile("/Users/gwon/WbLS/1ton_analysis/BNL1tonAnalysis/build/result_water.root");
    // TFile* mc = new
    // TFile("/Users/gwon/WbLS/1ton_analysis/BNL1tonAnalysis/build/test_1000.root");
    // TFile* mc = new
    // TFile("/Users/gwon/WbLS/1ton_analysis/BNL1tonAnalysis/build/test.root");
    // TFile* mc = new
    // TFile("/Users/gwon/WbLS/1ton_analysis/BNL1tonAnalysis/build/test_03wbls.root");
    // TFile* mc_noLY = new
    // TFile("/Users/gwon/WbLS/1ton_analysis/BNL1tonAnalysis/build/test_03wbls_noLY.root");

    const int numbin = 40;
    const int pelim = numbin * 20;
    TH1D *hist_data_bottomPE =
        new TH1D("hist_data_bottomPE", "hist_data_bottomPE", numbin, 0, pelim);
    TH1D *hist_data_sidePE =
        new TH1D("hist_data_sidePE", "hist_data_sidePE", numbin, 0, pelim);
    // TH1D* hist_data_sidePE = new TH1D("hist_data_sidePE", "hist_data_sidePE",
    // 20,0,200);
    TH1D *hist_data_totalPE =
        new TH1D("hist_data_totalPE", "hist_data_totalPE", numbin, 0, pelim);

    TH1D *hist_mc_bottomPE =
        new TH1D("hist_mc_bottomPE", "hist_mc_bottomPE", numbin, 0, pelim);
    TH1D *hist_mc_sidePE =
        new TH1D("hist_mc_sidePE", "hist_mc_sidePE", numbin, 0, pelim);
    // TH1D* hist_mc_sidePE = new TH1D("hist_mc_sidePE", "hist_mc_sidePE",
    // 20,0,200);
    TH1D *hist_mc_totalPE =
        new TH1D("hist_mc_totalPE", "hist_mc_totalPE", numbin, 0, pelim);

    TH1D *hist_mc_bottomPE_noLY = new TH1D(
        "hist_mc_bottomPE_noLY", "hist_mc_bottomPE_noLY", numbin, 0, pelim);
    TH1D *hist_mc_sidePE_noLY = new TH1D(
        "hist_mc_sidePE_noLY", "hist_mc_sidePE_noLY", numbin, 0, pelim);
    TH1D *hist_mc_totalPE_noLY = new TH1D(
        "hist_mc_totalPE_noLY", "hist_mc_totalPE_noLY", numbin, 0, pelim);

    TH1D *hist_mc_bottomPE_1000 = new TH1D(
        "hist_mc_bottomPE_1000", "hist_mc_bottomPE_1000", numbin, 0, pelim);
    TH1D *hist_mc_sidePE_1000 = new TH1D(
        "hist_mc_sidePE_1000", "hist_mc_sidePE_1000", numbin, 0, pelim);
    TH1D *hist_mc_totalPE_1000 = new TH1D(
        "hist_mc_totalPE_1000", "hist_mc_totalPE_1000", numbin, 0, pelim);

    TTree *dataTree = (TTree *)data->Get("tree");
    double adc_b1_ch1;
    dataTree->SetBranchAddress("adc_b1_ch1", &adc_b1_ch1);
    double adc_b1_ch2;
    dataTree->SetBranchAddress("adc_b1_ch2", &adc_b1_ch2);
    double adc_b1_ch3;
    dataTree->SetBranchAddress("adc_b1_ch3", &adc_b1_ch3);
    double adc_b1_ch4;
    dataTree->SetBranchAddress("adc_b1_ch4", &adc_b1_ch4);
    double adc_b1_ch5;
    dataTree->SetBranchAddress("adc_b1_ch5", &adc_b1_ch5);
    double adc_b1_ch6;
    dataTree->SetBranchAddress("adc_b1_ch6", &adc_b1_ch6);
    double adc_b1_ch7;
    dataTree->SetBranchAddress("adc_b1_ch7", &adc_b1_ch7);
    double adc_b1_ch8;
    dataTree->SetBranchAddress("adc_b1_ch8", &adc_b1_ch8);
    double adc_b1_ch9;
    dataTree->SetBranchAddress("adc_b1_ch9", &adc_b1_ch9);
    double adc_b1_ch10;
    dataTree->SetBranchAddress("adc_b1_ch10", &adc_b1_ch10);
    double adc_b1_ch11;
    dataTree->SetBranchAddress("adc_b1_ch11", &adc_b1_ch11);
    double adc_b1_ch12;
    dataTree->SetBranchAddress("adc_b1_ch12", &adc_b1_ch12);
    double adc_b1_ch13;
    dataTree->SetBranchAddress("adc_b1_ch13", &adc_b1_ch13);
    double adc_b1_ch14;
    dataTree->SetBranchAddress("adc_b1_ch14", &adc_b1_ch14);
    double adc_b1_ch15;
    dataTree->SetBranchAddress("adc_b1_ch15", &adc_b1_ch15);
    double adc_b2_ch0;
    dataTree->SetBranchAddress("adc_b2_ch0", &adc_b2_ch0);
    double adc_b2_ch1;
    dataTree->SetBranchAddress("adc_b2_ch1", &adc_b2_ch1);
    double adc_b2_ch2;
    dataTree->SetBranchAddress("adc_b2_ch2", &adc_b2_ch2);
    double adc_b2_ch3;
    dataTree->SetBranchAddress("adc_b2_ch3", &adc_b2_ch3);
    double adc_b2_ch4;
    dataTree->SetBranchAddress("adc_b2_ch4", &adc_b2_ch4);
    double adc_b2_ch5;
    dataTree->SetBranchAddress("adc_b2_ch5", &adc_b2_ch5);
    double adc_b2_ch6;
    dataTree->SetBranchAddress("adc_b2_ch6", &adc_b2_ch6);
    double adc_b2_ch7;
    dataTree->SetBranchAddress("adc_b2_ch7", &adc_b2_ch7);
    double adc_b2_ch8;
    dataTree->SetBranchAddress("adc_b2_ch8", &adc_b2_ch8);
    double adc_b2_ch9;
    dataTree->SetBranchAddress("adc_b2_ch9", &adc_b2_ch9);
    double adc_b2_ch10;
    dataTree->SetBranchAddress("adc_b2_ch10", &adc_b2_ch10);
    double adc_b2_ch11;
    dataTree->SetBranchAddress("adc_b2_ch11", &adc_b2_ch11);
    double adc_b2_ch12;
    dataTree->SetBranchAddress("adc_b2_ch12", &adc_b2_ch12);
    double adc_b2_ch13;
    dataTree->SetBranchAddress("adc_b2_ch13", &adc_b2_ch13);
    double adc_b2_ch14;
    dataTree->SetBranchAddress("adc_b2_ch14", &adc_b2_ch14);
    double adc_b3_ch0;
    dataTree->SetBranchAddress("adc_b3_ch0", &adc_b3_ch0);
    double adc_b3_ch1;
    dataTree->SetBranchAddress("adc_b3_ch1", &adc_b3_ch1);
    double adc_b3_ch2;
    dataTree->SetBranchAddress("adc_b3_ch2", &adc_b3_ch2);
    double adc_b3_ch3;
    dataTree->SetBranchAddress("adc_b3_ch3", &adc_b3_ch3);
    double adc_b3_ch4;
    dataTree->SetBranchAddress("adc_b3_ch4", &adc_b3_ch4);
    double adc_b3_ch5;
    dataTree->SetBranchAddress("adc_b3_ch5", &adc_b3_ch5);
    double adc_b3_ch6;
    dataTree->SetBranchAddress("adc_b3_ch6", &adc_b3_ch6);
    double adc_b3_ch7;
    dataTree->SetBranchAddress("adc_b3_ch7", &adc_b3_ch7);
    double adc_b3_ch8;
    dataTree->SetBranchAddress("adc_b3_ch8", &adc_b3_ch8);
    double adc_b3_ch9;
    dataTree->SetBranchAddress("adc_b3_ch9", &adc_b3_ch9);
    double adc_b3_ch10;
    dataTree->SetBranchAddress("adc_b3_ch10", &adc_b3_ch10);
    double adc_b3_ch11;
    dataTree->SetBranchAddress("adc_b3_ch11", &adc_b3_ch11);
    double adc_b3_ch12;
    dataTree->SetBranchAddress("adc_b3_ch12", &adc_b3_ch12);
    double adc_b3_ch13;
    dataTree->SetBranchAddress("adc_b3_ch13", &adc_b3_ch13);
    double adc_b3_ch14;
    dataTree->SetBranchAddress("adc_b3_ch14", &adc_b3_ch14);
    double adc_b3_ch15;
    dataTree->SetBranchAddress("adc_b3_ch15", &adc_b3_ch15);
    double adc_b4_ch0;
    dataTree->SetBranchAddress("adc_b4_ch0", &adc_b4_ch0);
    double adc_b4_ch1;
    dataTree->SetBranchAddress("adc_b4_ch1", &adc_b4_ch1);
    double adc_b4_ch2;
    dataTree->SetBranchAddress("adc_b4_ch2", &adc_b4_ch2);
    double adc_b4_ch3;
    dataTree->SetBranchAddress("adc_b4_ch3", &adc_b4_ch3);
    double adc_b4_ch4;
    dataTree->SetBranchAddress("adc_b4_ch4", &adc_b4_ch4);
    double adc_b4_ch5;
    dataTree->SetBranchAddress("adc_b4_ch5", &adc_b4_ch5);
    double adc_b4_ch6;
    dataTree->SetBranchAddress("adc_b4_ch6", &adc_b4_ch6);
    double adc_b4_ch7;
    dataTree->SetBranchAddress("adc_b4_ch7", &adc_b4_ch7);
    double adc_b4_ch8;
    dataTree->SetBranchAddress("adc_b4_ch8", &adc_b4_ch8);
    double adc_b4_ch9;
    dataTree->SetBranchAddress("adc_b4_ch9", &adc_b4_ch9);
    double adc_b4_ch10;
    dataTree->SetBranchAddress("adc_b4_ch10", &adc_b4_ch10);
    double adc_b4_ch11;
    dataTree->SetBranchAddress("adc_b4_ch11", &adc_b4_ch11);
    for (int ievt = 0; ievt < dataTree->GetEntries(); ++ievt) {
        dataTree->GetEntry(ievt);
        int tempBottomPE = 0;
        int tempSidePE = 0;
        int tempTotalPE = 0;
        tempBottomPE += adc_b1_ch1;
        tempBottomPE += adc_b1_ch2;
        tempBottomPE += adc_b1_ch3;
        tempBottomPE += adc_b1_ch4;
        tempBottomPE += adc_b1_ch5;
        tempBottomPE += adc_b1_ch6;
        tempBottomPE += adc_b1_ch7;
        tempBottomPE += adc_b1_ch8;
        tempBottomPE += adc_b1_ch9;
        tempBottomPE += adc_b1_ch10;
        tempBottomPE += adc_b1_ch11;
        tempBottomPE += adc_b1_ch12;
        tempBottomPE += adc_b1_ch13;
        tempBottomPE += adc_b1_ch14;
        tempBottomPE += adc_b1_ch15;
        tempBottomPE += adc_b2_ch0;
        tempBottomPE += adc_b2_ch1;
        tempBottomPE += adc_b2_ch2;
        tempBottomPE += adc_b2_ch3;
        tempBottomPE += adc_b2_ch4;
        tempBottomPE += adc_b2_ch5;
        tempBottomPE += adc_b2_ch6;
        tempBottomPE += adc_b2_ch7;
        tempBottomPE += adc_b2_ch8;
        tempBottomPE += adc_b2_ch9;
        tempBottomPE += adc_b2_ch10;
        tempBottomPE += adc_b2_ch11;
        tempBottomPE += adc_b2_ch12;
        tempBottomPE += adc_b2_ch13;
        tempBottomPE += adc_b2_ch14;

        tempSidePE += adc_b3_ch0;
        tempSidePE += adc_b3_ch1;
        tempSidePE += adc_b3_ch2;
        tempSidePE += adc_b3_ch3;
        tempSidePE += adc_b3_ch4;
        tempSidePE += adc_b3_ch5;
        tempSidePE += adc_b3_ch6;
        tempSidePE += adc_b3_ch7;
        tempSidePE += adc_b3_ch8;
        tempSidePE += adc_b3_ch9;
        tempSidePE += adc_b3_ch10;
        tempSidePE += adc_b3_ch11;
        tempSidePE += adc_b3_ch12;
        tempSidePE += adc_b3_ch13;
        tempSidePE += adc_b3_ch14;
        tempSidePE += adc_b3_ch15;

        tempTotalPE += adc_b1_ch1;
        tempTotalPE += adc_b1_ch2;
        tempTotalPE += adc_b1_ch3;
        tempTotalPE += adc_b1_ch4;
        tempTotalPE += adc_b1_ch5;
        tempTotalPE += adc_b1_ch6;
        tempTotalPE += adc_b1_ch7;
        tempTotalPE += adc_b1_ch8;
        tempTotalPE += adc_b1_ch9;
        tempTotalPE += adc_b1_ch10;
        tempTotalPE += adc_b1_ch11;
        tempTotalPE += adc_b1_ch12;
        tempTotalPE += adc_b1_ch13;
        tempTotalPE += adc_b1_ch14;
        tempTotalPE += adc_b1_ch15;
        tempTotalPE += adc_b2_ch0;
        tempTotalPE += adc_b2_ch1;
        tempTotalPE += adc_b2_ch2;
        tempTotalPE += adc_b2_ch3;
        tempTotalPE += adc_b2_ch4;
        tempTotalPE += adc_b2_ch5;
        tempTotalPE += adc_b2_ch6;
        tempTotalPE += adc_b2_ch7;
        tempTotalPE += adc_b2_ch8;
        tempTotalPE += adc_b2_ch9;
        tempTotalPE += adc_b2_ch10;
        tempTotalPE += adc_b2_ch11;
        tempTotalPE += adc_b2_ch12;
        tempTotalPE += adc_b2_ch13;
        tempTotalPE += adc_b2_ch14;

        tempTotalPE += adc_b3_ch0;
        tempTotalPE += adc_b3_ch1;
        tempTotalPE += adc_b3_ch2;
        tempTotalPE += adc_b3_ch3;
        tempTotalPE += adc_b3_ch4;
        tempTotalPE += adc_b3_ch5;
        tempTotalPE += adc_b3_ch6;
        tempTotalPE += adc_b3_ch7;
        tempTotalPE += adc_b3_ch8;
        tempTotalPE += adc_b3_ch9;
        tempTotalPE += adc_b3_ch10;
        tempTotalPE += adc_b3_ch11;
        tempTotalPE += adc_b3_ch12;
        tempTotalPE += adc_b3_ch13;
        tempTotalPE += adc_b3_ch14;
        tempTotalPE += adc_b3_ch15;

        hist_data_bottomPE->Fill(tempBottomPE);
        hist_data_sidePE->Fill(tempSidePE);
        hist_data_totalPE->Fill(tempTotalPE);
    }

    TTree *mcTree = (TTree *)mc->Get("tree");
    int PMTPE[46];
    if (mcTree->SetBranchAddress("PMTPE", PMTPE) != 0) {
        std::cerr << "Error: Could not set branch address for 'PMTPE'."
                  << std::endl;
        return;
    } else {
        std::cout << "Successfully set branch address for 'PMTPE'."
                  << std::endl;
    }

    for (int ievt = 0; ievt < mcTree->GetEntries(); ++ievt) {
        mcTree->GetEntry(ievt);
        int tempBottomPE = 0;
        int tempSidePE = 0;
        int tempTotalPE = 0;
        for (int pmtid = 0; pmtid < 46; ++pmtid) {
            tempTotalPE += PMTPE[pmtid];
            if (pmtid < 30)
                tempBottomPE += PMTPE[pmtid];
            else
                tempSidePE += PMTPE[pmtid];
        }
        hist_mc_bottomPE->Fill(tempBottomPE);
        hist_mc_sidePE->Fill(tempSidePE);
        hist_mc_totalPE->Fill(tempTotalPE);
    }

    TTree *mcTree_noLY = (TTree *)mc_noLY->Get("tree");
    int PMTPE_noLY[46];
    if (mcTree_noLY->SetBranchAddress("PMTPE", PMTPE_noLY) != 0) {
        std::cerr << "Error: Could not set branch address for 'PMTPE'."
                  << std::endl;
        return;
    } else {
        std::cout << "Successfully set branch address for 'PMTPE'."
                  << std::endl;
    }
    for (int ievt = 0; ievt < mcTree_noLY->GetEntries(); ++ievt) {
        mcTree_noLY->GetEntry(ievt);
        int tempBottomPE = 0;
        int tempSidePE = 0;
        int tempTotalPE = 0;
        for (int pmtid = 0; pmtid < 46; ++pmtid) {
            tempTotalPE += PMTPE_noLY[pmtid];
            if (pmtid < 30)
                tempBottomPE += PMTPE_noLY[pmtid];
            else
                tempSidePE += PMTPE_noLY[pmtid];
        }
        hist_mc_bottomPE_noLY->Fill(tempBottomPE);
        hist_mc_sidePE_noLY->Fill(tempSidePE);
        hist_mc_totalPE_noLY->Fill(tempTotalPE);
    }

    TTree *mcTree_1000 = (TTree *)mc_1000->Get("tree");
    int PMTPE_1000[46];
    if (mcTree_1000->SetBranchAddress("PMTPE", PMTPE_1000) != 0) {
        std::cerr << "Error: Could not set branch address for 'PMTPE'."
                  << std::endl;
        return;
    } else {
        std::cout << "Successfully set branch address for 'PMTPE'."
                  << std::endl;
    }
    for (int ievt = 0; ievt < mcTree_1000->GetEntries(); ++ievt) {
        mcTree_1000->GetEntry(ievt);
        int tempBottomPE = 0;
        int tempSidePE = 0;
        int tempTotalPE = 0;
        for (int pmtid = 0; pmtid < 46; ++pmtid) {
            tempTotalPE += PMTPE_1000[pmtid];
            if (pmtid < 30)
                tempBottomPE += PMTPE_1000[pmtid];
            else
                tempSidePE += PMTPE_1000[pmtid];
        }
        hist_mc_bottomPE_1000->Fill(tempBottomPE);
        hist_mc_sidePE_1000->Fill(tempSidePE);
        hist_mc_totalPE_1000->Fill(tempTotalPE);
    }

    // gStyle->SetOptStat(0);
    TCanvas *can = new TCanvas;
    can->Divide(2, 2);

    // Create a legend
    TLegend *leg1 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    TLegend *leg2 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    TLegend *leg3 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg3->SetBorderSize(0);
    leg3->SetFillStyle(0);

    // First pad
    can->cd(1);
    hist_mc_bottomPE->SetLineColor(2);
    hist_mc_bottomPE->DrawNormalized();
    CreateLegend(leg1, hist_mc_bottomPE, "MC Bottom PE");
    hist_mc_bottomPE_noLY->SetLineColor(1);
    hist_mc_bottomPE_noLY->DrawNormalized("same");
    CreateLegend(leg1, hist_mc_bottomPE_noLY, "MC Bottom PE no LY");
    // hist_mc_bottomPE_1000->SetLineColor(1);
    // hist_mc_bottomPE_1000->DrawNormalized("same");
    // CreateLegend(leg1, hist_mc_bottomPE_1000, "MC Bottom PE LY 1000");
    hist_data_bottomPE->DrawNormalized("same");
    CreateLegend(leg1, hist_data_bottomPE, "Data Bottom PE");
    // Draw the legend
    leg1->Draw();

    // Second pad
    can->cd(2);
    hist_mc_sidePE->SetLineColor(2);
    hist_mc_sidePE->DrawNormalized();
    CreateLegend(leg2, hist_mc_sidePE, "MC Side PE");
    hist_mc_sidePE_noLY->SetLineColor(1);
    // hist_mc_sidePE_noLY->DrawNormalized("same");
    // CreateLegend(leg2, hist_mc_sidePE_noLY, "MC Side PE no LY");
    hist_data_sidePE->DrawNormalized("same");
    CreateLegend(leg2, hist_data_sidePE, "Data Side PE");
    // Draw the legend
    leg2->Draw();

    // Third pad
    can->cd(3);
    hist_mc_totalPE->SetLineColor(2);
    hist_mc_totalPE->DrawNormalized();
    CreateLegend(leg3, hist_mc_totalPE, "MC Total PE");
    hist_mc_totalPE_noLY->SetLineColor(1);
    hist_mc_totalPE_noLY->DrawNormalized("same");
    CreateLegend(leg3, hist_mc_totalPE_noLY, "MC Total PE no LY");
    hist_data_totalPE->DrawNormalized("same");
    CreateLegend(leg3, hist_data_totalPE, "Data Total PE");
    // Draw the legend
    leg3->Draw();

    TCanvas *can2 = new TCanvas;
    can2->Divide(2, 2);
    can2->cd(1);
    hist_mc_bottomPE->Draw();
    can2->cd(2);
    hist_mc_sidePE->Draw();
    can2->cd(3);
    hist_mc_bottomPE_noLY->Draw();
    can2->cd(4);
    hist_mc_sidePE_noLY->Draw();

    TCanvas *can3 = new TCanvas;
    can3->Divide(2, 1);
    can3->cd(1);
    hist_data_bottomPE->Draw();
    std::cout << "bottome Truncated Mean: " << GetTruncatedMean(hist_data_bottomPE, 200, 400) << std::endl;
    can3->cd(2);
    hist_data_sidePE->Draw();
    std::cout << "side Truncated Mean: " << GetTruncatedMean(hist_data_sidePE, 100, 200) << std::endl;
}

