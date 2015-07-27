//Standard Library
#include <string>
#include <iostream>

//Root Includes
#include <TH1D.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>

int DataJetEtPlotter()
{
    std::string inFile = "/home/aevans/CMS/slc6test/CMSSW_7_1_0_pre5/src/TestEDA/JetAlgorithm/plugins/newDataJetTest.root";
    std::string outFile = "/home/aevans/CMS/slc6test/CMSSW_7_1_0_pre5/src/TestEDA/JetAlgorithm/plugins/newDataJetTestPlots.root";

    TFile* inTfile = new TFile(inFile.c_str());

    TStyle* style = new TStyle("Style", "Style for P-TDR");  //roughly taken from Alex's ToyMC style.py in ZFinder
    style->SetCanvasBorderMode(0);
    style->SetCanvasColor(kWhite);
    style->SetCanvasDefX(0);
    style->SetCanvasDefY(0);
    // For the Pad:
    style->SetPadBorderMode(0);
    style->SetPadColor(kWhite);
    style->SetPadGridX(false);
    style->SetPadGridY(false);
    style->SetGridColor(kBlack);
    style->SetGridStyle(3);
    style->SetGridWidth(1);
    // For the frame:
    style->SetFrameBorderMode(0);
    style->SetFrameBorderSize(1);
    style->SetFrameFillColor(kWhite);
    style->SetFrameFillStyle(0);
    style->SetFrameLineColor(kBlack);
    style->SetFrameLineStyle(1);
    style->SetFrameLineWidth(1);


    const int numberOfHists = 4;
    TH1D* histograms[numberOfHists];
    TH1D* intHistograms[numberOfHists];
    char histName[1024]; 
    
    inTfile->GetObject("/jetFinder/trigJetEtwHEcoin", histograms[3]);
    inTfile->GetObject("/jetFinder/trigJetEt", histograms[0]);
    inTfile->GetObject("/jetFinder/L1ForJetPt", histograms[1]);
    inTfile->GetObject("/jetFinder/numTrigJets", histograms[2]);
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

    TCanvas* canvas1 = new TCanvas;
    canvas1->Divide(1,2);
    p2 = canvas1->cd(2);
    p1 = canvas1->cd(1);
    p1->SetLogy();
    int numEvents = histograms[2]->GetEntries();
    cout << numEvents << endl;
    histograms[1]->Scale(1.0/numEvents);
    histograms[0]->Scale(1.0/numEvents);
    histograms[0]->Rebin(5);
    
    histograms[1]->Clone("diffRatio");
    diffRatio->Divide(histograms[0]);
    
    histograms[1]->Draw("histe");
    histograms[0]->Draw("sameshiste");
    histograms[1]->SetLineColor(2);
    histograms[0]->SetLineColor(4);

    canvas1->cd(2);
    diffRatio->Draw("histe");
    
    TCanvas* canvas2 = new TCanvas;
    canvas2->Divide(1,2);
    p4 = canvas2->cd(2);
    p3 = canvas2->cd(1);
    p3->SetLogy();
    
    for(int j = numberOfHists-3; j >= 0; j--)
    {
        
        sprintf(histName, "hist_%d",j);
        
        TH1D* histo = histograms[j]->Clone(histName);
        for(int ibin = 1; ibin <= histograms[j]->GetNbinsX(); ibin++)
        {
            double integral = histo->GetIntegral();
            histo->SetBinContent(ibin, histo->Integral(ibin,histograms[j]->GetNbinsX()));
        }
        histo->Draw("sameshiste");
        histo->SetLineColor(1+j*2);
        intHistograms[j]=histo;
    }
    intHistograms[1]->Clone("intRatio");
    intRatio->Divide(intHistograms[0]);
    canvas2->cd(2);
    intRatio->Draw("histe");
    





    TFile f(outFile.c_str(), "recreate");
    canvas1->Write();
    canvas2->Write();
    
    
    f.Close();







}
