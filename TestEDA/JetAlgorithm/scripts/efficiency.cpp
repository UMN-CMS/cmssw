//Standard Library
#include <string>
#include <iostream>

//Root Includes
#include <TH1D.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>

int efficiency()
{
    std::string inFile = "/home/aevans/CMS/CMSSW_7_1_0_pre5/src/TestEDA/JetAlgorithm/plugins/jets.root";
    std::string outFile = "eff.root";

    TFile* inTfile = new TFile(inFile.c_str());

    const int numberOfHists = 8;
    TH1D* histograms[numberOfHists];
    inTfile->GetObject("/jetFinder/EtMatched", histograms[0]);
    inTfile->GetObject("/jetFinder/Et", histograms[1]);
    inTfile->GetObject("/jetFinder/EtaMatched", histograms[2]);
    inTfile->GetObject("/jetFinder/Eta", histograms[3]);
    inTfile->GetObject("/jetFinder/PhiMatched", histograms[4]);
    inTfile->GetObject("/jetFinder/Phi", histograms[5]);
    inTfile->GetObject("/jetFinder/absEtaMatched", histograms[6]);
    inTfile->GetObject("/jetFinder/absEta", histograms[7]);
    for (int i = 0; i < numberOfHists; ++i)
    {
        if (!histograms[i])
        {
            std::cout << "Histogram " << i << " failed to load!" << std::endl;
            return 1;
        }
        else
        {
            histograms[i]->Sumw2();
        }
    }

    TGraphAsymmErrors* etEffPlot = new TGraphAsymmErrors(histograms[0], histograms[1], "-b");
    TGraphAsymmErrors* etaEffPlot = new TGraphAsymmErrors(histograms[2], histograms[3], "-b");
    TGraphAsymmErrors* phiEffPlot = new TGraphAsymmErrors(histograms[4], histograms[5], "-b");
    TGraphAsymmErrors* absEtaEffPlot = new TGraphAsymmErrors(histograms[6], histograms[7], "-b");

    TFile f(outFile.c_str(), "recreate");
    etEffPlot->Write();
    etaEffPlot->Write();
    phiEffPlot->Write();
    absEtaEffPlot->Write();
    f.Close();







}
