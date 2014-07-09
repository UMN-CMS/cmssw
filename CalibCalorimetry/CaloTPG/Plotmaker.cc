#include <TLegend.h>
#include <TDirectory.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TProfile.h>
#include <math.h> 
#include <map>
#include "TH1.h"
#include "TH2.h"
#include <string>
#include <TFile.h>
#include <iostream>
#include<TPaveStats.h>

//////////
#include <TGraphErrors.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TF1.h>
#include <iostream>
///////
#include "TObject.h"
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <TMath.h>
#include <TTree.h>
#include <TString.h>
#include <TSystem.h>
#include "TLorentzVector.h"
#include <TDirectory.h>


#include <cmath> 
#include <map>
#include "Math/GenVector/VectorUtil.h"


#include <TMultiGraph.h>
#include <TGraphErrors.h>

#include <cmath> 
#define PI 3.14159265
using namespace std;

void DataHistoGram(string FileOrgin)
{
    string type = ".png";
    string plotname;

    FileOrgin = "fdata/" + FileOrgin;
    string hold = FileOrgin + ".root";

    TFile* theFile = new TFile(hold.c_str());
    if(!theFile)
    {
        cout << "Error. Can not open file. try again" << endl;
        return;
    }

    TH1 *PsdEt, *PsdRtoEt, *PsdHit; //passed Et and Passed ratio of a 3x3 square to the location
    //TH2 *PsdEtvrsiEta, *PsdEtvrsiPhi;
    //We will figure this out later terrible way of doing this for now
    string PsdEtLoc = "/DataAna/finebit/Et"; // The location of both psdEt and PsdRtoEt
    string PsdRtoEtLoc = "/DataAna/finebit/EtRatio";
    string PsdhitLoc = "/DataAna/finebit/Hits";
    string PsdEtvrsiEtaLoc;
    string PsdEtvrsiPhiLoc;

    PsdEt = (TH1F*) theFile->Get(PsdEtLoc.c_str());
    cout << " Done opening " << PsdEtLoc << endl;
    PsdRtoEt = (TH1F*) theFile->Get(PsdRtoEtLoc.c_str());
    cout << " Done opening " << PsdRtoEtLoc << endl;
    PsdHit = (TH1I*) theFile->Get(PsdhitLoc.c_str());
    cout << "Done opening " << PsdhitLoc << endl;

    gSystem->mkdir(("PDFFiles/" + FileOrgin).c_str());

    TCanvas* c1 = new TCanvas("c1", "", 800, 700);
    c1->SetFillColor(10);
    c1->SetFillColor(10);
    c1->cd();

    
    gStyle->SetOptStat("eoumi");

    gStyle->SetStatY(0.9);
    // Set y-position (fraction of pad size)
    gStyle->SetStatX(0.9);
    // Set x-position (fraction of pad size)
    gStyle->SetStatW(0.2);
    // Set width of stat-box (fraction of pad size)
    gStyle->SetStatH(0.1);
    // Set height of stat-box (fraction of pad size)
    PsdEt->GetXaxis()->SetTitle("Et(GeV)");
    PsdEt->Draw();


    //gStyle->SetOptStat("");
    c1->SetLogy();

    c1->Update();
    plotname = "PDFFiles/" + FileOrgin + "/Passed_Et" + type;

    c1->Print(plotname.c_str());
    delete c1;

    TCanvas* c2 = new TCanvas("c2", "", 800, 700);


    c2->SetFillColor(10);
    c2->SetFillColor(10);
    c2->cd();
    gStyle->SetOptStat("eoumi");
    gStyle->SetStatY(0.9);
    // Set y-position (fraction of pad size)
    gStyle->SetStatX(0.33);
    // Set x-position (fraction of pad size)
    gStyle->SetStatW(0.2);
    // Set width of stat-box (fraction of pad size)
    gStyle->SetStatH(0.1);
    // Set height of stat-box (fraction of pad size)
    PsdRtoEt->GetXaxis()->SetTitle("Number of HF Trigger candidates");
    PsdRtoEt->Draw();


    c2->SetLogy();
    /////////////////////

    c2->Update();
    plotname = "PDFFiles/" + FileOrgin + "/Passed_Et_ratio" + type;
    c2->Print(plotname.c_str());
    delete c2;

    TCanvas* c3 = new TCanvas("c3", "", 800, 700);

    c3->SetFillColor(10);
    c3->SetFillColor(10);
    c3->cd();

    gStyle->SetOptStat("eoumi");

    gStyle->SetStatY(0.9);
    // Set y-position (fraction of pad size)
    gStyle->SetStatX(0.9);
    // Set x-position (fraction of pad size)
    gStyle->SetStatW(0.2);
    // Set width of stat-box (fraction of pad size)
    gStyle->SetStatH(0.1);
    PsdHit->Draw();



    c3->SetLogy();
    c3->Update();
    plotname = "PDFFiles/" + FileOrgin + "/Hits" + type;

    c3->Print(plotname.c_str());
    delete c3;
}

void DataCompareHistoGram(string FileOrgin)
{
    string type = ".png";
    string plotname;
    string saveloc = "fdataComp/" + FileOrgin;

    FileOrgin = "fdata/" + FileOrgin;
    string hold = FileOrgin + ".root";

    TFile* theFile = new TFile(hold.c_str());
    if(!theFile)
    {
        cout << "Error. Can not open file. try again" << endl;
        return;
    }

    TH1 *PsdEt, *PsdRtoEt, *PsdHit; //passed Et and Passed ratio of a 3x3 square to the location
    TH1 *PsdAllEt, *PsdAllRtoEt, *PsdAllHit;
    //TH2 *PsdEtvrsiEta, *PsdEtvrsiPhi;

    //We will figure this out later terrible way of doing this for now
    string PsdEtLoc = "/DataAna/finebit/Et"; // The location of both psdEt and PsdRtoEt
    string PsdRtoEtLoc = "/DataAna/finebit/EtRatio";
    string PsdhitLoc = "/DataAna/finebit/Hits";
    string PsdAllEtLoc = "/DataAna/no_fine/Et"; // The location of both psdEt and PsdRtoEt
    string PsdAllRtoEtLoc = "/DataAna/no_fine/EtRatio";
    string PsdAllhitLoc = "/DataAna/no_fine/Hits";
    string PsdEtvrsiEtaLoc;
    string PsdEtvrsiPhiLoc;

    PsdEt = (TH1F*) theFile->Get(PsdEtLoc.c_str());
    cout << " Done opening " << PsdEtLoc << endl;
    PsdRtoEt = (TH1F*) theFile->Get(PsdRtoEtLoc.c_str());
    cout << " Done opening " << PsdRtoEtLoc << endl;
    PsdHit = (TH1I*) theFile->Get(PsdhitLoc.c_str());
    cout << "Done opening " << PsdhitLoc << endl;

    PsdAllEt = (TH1F*) theFile->Get(PsdAllEtLoc.c_str());
    cout << " Done opening " << PsdAllEtLoc << endl;
    PsdAllRtoEt = (TH1F*) theFile->Get(PsdAllRtoEtLoc.c_str());
    cout << " Done opening " << PsdAllRtoEtLoc << endl;
    PsdAllHit = (TH1I*) theFile->Get(PsdAllhitLoc.c_str());
    cout << "Done opening " << PsdAllhitLoc << endl;
    gSystem->mkdir(("PDFFiles/" + saveloc).c_str());
    TCanvas* c1 = new TCanvas("c1", "", 800, 700);
    c1->SetFillColor(10);
    c1->SetFillColor(10);

    c1->cd();


    PsdEt->SetLineColor(8);

    PsdAllEt->SetLineColor(2);
    PsdAllEt->Draw();
    PsdEt->Draw("same");
    gStyle->SetOptStat("");
    c1->SetLogy();
    /////////////////////
    TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg->SetHeader("Et");
    leg-> SetFillColor(kWhite);
    leg-> SetShadowColor(kWhite);
    leg->AddEntry(PsdEt, "Et using trigprimcuts");
    leg->AddEntry(PsdAllEt, "Et no trigprimcuts");
    leg->Draw();
    c1->Update();
    plotname = "PDFFiles/" + saveloc + "/Passed_Et" + type;

    c1->Print(plotname.c_str());
    delete c1;

    TCanvas* c2 = new TCanvas("c2", "", 800, 700);


    c2->SetFillColor(10);
    c2->SetFillColor(10);
    c2->cd();
    PsdRtoEt->Draw();
    PsdAllRtoEt->SetLineColor(8);
    PsdAllRtoEt->Draw("same");
    //gStyle->SetOptStat("eoumi");
    gStyle->SetOptStat("");
    leg = new TLegend(0.15, 0.7, 0.4, 0.9);
    leg->SetHeader("Et/3X3");
    leg-> SetFillColor(kWhite);
    leg-> SetShadowColor(kWhite);
    leg->AddEntry(PsdRtoEt, "Ratio using trigprimcuts");
    leg->AddEntry(PsdAllRtoEt, "Ratio no trigprimcuts");
    leg->Draw();
    // gStyle->SetOptStat("");

    c2->SetLogy();
    /////////////////////

    c2->Update();
    plotname = "PDFFiles/" + saveloc + "/Passed_Et_ratio" + type;
    c2->Print(plotname.c_str());
    delete c2;

    TCanvas* c3 = new TCanvas("c3", "", 800, 700);

    c3->SetFillColor(10);
    c3->SetFillColor(10);
    c3->cd();
    //PsdHit->SetMaximum(1.1 * PsdAllHit->GetBin(PsdAllHit->GetMaximumBin()));
    PsdAllHit->Draw();
    PsdHit->SetLineColor(8);
    PsdHit->Draw("same");

    //gStyle->SetOptStat("eoumi");
    gStyle->SetOptStat("");

    c3->SetLogy();
    leg = new TLegend(0.25, 0.7, 0.55, 0.9);
    leg->SetHeader("Hits for one event");
    leg-> SetFillColor(kWhite);
    leg-> SetShadowColor(kWhite);
    leg->AddEntry(PsdHit, "Hits trig cuts");
    leg->AddEntry(PsdAllHit, "Hits no trigprimcuts");
    leg->Draw();
    c3->Update();
    plotname = "PDFFiles/" + saveloc + "/Hits" + type;

    c3->Print(plotname.c_str());
    delete c3;

}

void MCHisto(string FileOrgin)
{
    string type = ".png";
    string plotname;
    string hold = "fmc/" + FileOrgin + ".root";

    TFile* theFile = new TFile(hold.c_str());
    if(!theFile)
    {
        cout << "Error. Can not open file. try again" << endl;
        return;
    }

    TH1 *PsdEt, *PsdRtoEt;
    TH2 *CompVrsGen, *UnCompVrsGen;
    string PsdEtLoc = "/CutAs/finebit/Et_of_Electron"; // The location of both psdEt and PsdRtoEt
    string PsdRtoEtLoc = "/CutAs/finebit/Et_Ratio_of_Electron";
    //string PsdhitLoc = "/CusAs/finebit/Hits";
    string ComELoc = "/CutAs/compressed_VS_Generated";
    string UnComELoc = "/CutAs/ETuncomp_VS_Generated";
    PsdEt = (TH1F*) theFile->Get(PsdEtLoc.c_str());
    //cout << " Done opening " << PsdEtLoc << endl;
    PsdRtoEt = (TH1F*) theFile->Get(PsdRtoEtLoc.c_str());
    // cout << " Done opening " << PsdRtoEtLoc << endl;
    CompVrsGen = (TH2F*) theFile->Get(ComELoc.c_str());
    cout << " Done opening " << ComELoc << endl;
    UnCompVrsGen = (TH2F*) theFile->Get(UnComELoc.c_str());
    cout << " Done opening " << UnComELoc << endl;


    cout << "attempted directory name :" << "PDFFiles/fmc/" + FileOrgin << endl;
    gSystem->mkdir(("PDFFiles/fmc/" + FileOrgin).c_str());
    TCanvas* c1 = new TCanvas("c1", "", 800, 700);
    c1->SetFillColor(10);
    c1->SetFillColor(10);

    c1->cd();
    PsdEt->GetXaxis()->SetTitle("Et(Gev)");
    gStyle->SetOptStat("eoumi");

    gStyle->SetStatY(0.9);
    // Set y-position (fraction of pad size)
    gStyle->SetStatX(0.9);
    // Set x-position (fraction of pad size)
    gStyle->SetStatW(0.2);
    // Set width of stat-box (fraction of pad size)
    gStyle->SetStatH(0.1);
    PsdEt->Draw();

    c1->SetLogy();
    /////////////////////

    c1->Update();
    plotname = "PDFFiles/fmc/" + FileOrgin + "/Passed_Et" + type;

    c1->Print(plotname.c_str());
    delete c1;


    TCanvas* c2 = new TCanvas("c1", "", 800, 700);
    c2->SetFillColor(10);
    c2->SetFillColor(10);

    c2->cd();

    PsdRtoEt->GetXaxis()->SetTitle("Et/3X3");
    gStyle->SetOptStat("eoumi");

    gStyle->SetStatY(0.9);
    // Set y-position (fraction of pad size)
    gStyle->SetStatX(0.33);
    // Set x-position (fraction of pad size)
    gStyle->SetStatW(0.2);
    // Set width of stat-box (fraction of pad size)
    gStyle->SetStatH(0.1);
    // Set height of stat-box (fraction of pad size)
    PsdRtoEt->Draw();

    c2->SetLogy();
    /////////////////////

    c2->Update();
    plotname = "PDFFiles/fmc/" + FileOrgin + "/Passed_EtRatio" + type;

    c2->Print(plotname.c_str());
    delete c2;

    TCanvas* c3 = new TCanvas("c3", "", 900, 700);
    c3->SetFillColor(10);
    c3->SetFillColor(10);
    c3->cd();
    gStyle->SetOptStat("");


    //int q;


    TF1 * Fit2 = new TF1("linefit", "x*[0]", 0, 120);
    gStyle->SetOptStat("");
    gStyle->SetOptFit(1011);
    CompVrsGen->Fit(Fit2, "RQ");


    //////////
    gPad->Update();


    TPaveStats *stats3 = (TPaveStats *) CompVrsGen->GetListOfFunctions()->FindObject("stats");
    stats3->SetX1NDC(.1);
    stats3->SetX2NDC(.3);

    cout << "test 3" << endl;
    CompVrsGen->Draw("COLZ");

    stats3->Draw();
    Fit2->Draw("sames");
    /////////////

    plotname = "PDFFiles/fmc/" + FileOrgin + "/CompressedVsGenEt" + type;

    c3->Print(plotname.c_str());
    delete c3;

    TCanvas* c6 = new TCanvas("c6", "", 900, 700);
    c6->SetFillColor(10);
    c6->SetFillColor(10);
    c6->cd();
    //gStyle->SetOptStat("");
    TF1 * Fit3 = new TF1("linefit", "x*[0]", 0, 120);
    gStyle->SetOptStat("");
    gStyle->SetOptFit(1011);
    UnCompVrsGen->Fit(Fit3, "RQ");
    cout << "test 2" << endl;

    //////////
    gPad->Update();
    cout << "test 2.1" << endl;

    TPaveStats *stats4 = (TPaveStats *) UnCompVrsGen->GetListOfFunctions()->FindObject("stats");
    stats4->SetX1NDC(.1);
    stats4->SetX2NDC(.3);

    cout << "test 3" << endl;
    UnCompVrsGen->Draw("COLZ");
    stats4->Draw();
    Fit3->Draw("sames");


    plotname = "PDFFiles/fmc/" + FileOrgin + "/UnCompressedVsGenEt" + type;

    c6->Print(plotname.c_str());
    delete c6;

    TCanvas* c4 = new TCanvas("c4", "", 900, 700);
    c4->SetFillColor(10);
    c4->SetFillColor(10);
    c4->cd();
    cout << "test 1" << endl;
    TProfile *ConGen = new TProfile("compressed_vrs_GenEt", "compressed_vrs_GenEt", 120, -0.001, 120, 0, 100);
    for(int i = 1; i < 120; i++)
    {
        for(int k = 1; k < 140; k++)
        {
            ConGen->Fill(i - 1, k - 1, CompVrsGen->GetBinContent(i, k));
            //cout<<" okay this is odd :"<<CompVrsGen->GetBinContent(i,k)<<endl;
        }
    }

    ConGen->GetXaxis()->SetTitle("Compressed Et");
    ConGen->GetYaxis()->SetTitle("Generated Et");
    TF1 * Fit0 = new TF1("linefit", "[0]+x*[1]", 0, 120);
    gStyle->SetOptStat("");
    gStyle->SetOptFit(1011);
    ConGen->Fit(Fit0, "RQ");
    ConGen->Draw();
    cout << "test 2" << endl;

    //////////
    gPad->Update();
    cout << "test 2.1" << endl;
    TPaveStats *stats = (TPaveStats *) ConGen->GetListOfFunctions()->FindObject("stats");
    stats->SetX1NDC(.1);
    stats->SetX2NDC(.3);

    cout << "test 3" << endl;

    stats->Draw();
    /////////////
    cout << "test 4" << endl;
    Fit0->Draw("sames");
    plotname = "PDFFiles/fmc/" + FileOrgin + "/CompVsGenEtD" + type;
    c4->Update();
    cout << "test 5" << endl;
    c4->Print(plotname.c_str());
    delete c4;
    //////////////////////////////////////////////
    cout << "test 1" << endl;
    TCanvas* c5 = new TCanvas("c5", "", 900, 700);
    c5->SetFillColor(10);
    c5->SetFillColor(10);
    c5->cd();


    TProfile *UnConGen = new TProfile("Uncompressed_vrs_GenEt", "Uncompressed_vrs_GenEt", 120, 0, 120, 0, 100);
    for(int i = 0; i < 120; i++)
    {
        for(int k = 0; k < 140; k++)
        {
            UnConGen->Fill(i, k, UnCompVrsGen->GetBinContent(i, k));
            //cout<<" okay this is odd :"<<CompVrsGen->GetBinContent(i,k)<<endl;
        }
    }

    UnCompVrsGen->GetXaxis()->SetTitle("Compressed Et");
    UnCompVrsGen->GetYaxis()->SetTitle("Generated Et");
    cout << "test 3" << endl;

    //TF1 * Fit1 = new TF1("linefit", "[0]*x", 0, 120);
    TF1 * Fit1 = new TF1("linefit", "[0]+x*[1]", 0, 120);
    UnConGen->Fit(Fit1, "RQ");
    //Fit0->(0, Fit0->GetParameter(0));
    UnConGen->Draw();
    gStyle->SetOptStat("");
    gStyle->SetOptFit(1011);

    //////////
    TPaveStats *stats2 = (TPaveStats *) UnConGen->GetListOfFunctions()->FindObject("stats");
    stats2->SetX1NDC(.1);
    stats2->SetX2NDC(.3);

    gPad->Update();
    stats2->Draw();
    /////////////

    cout << "test 4" << endl;
    plotname = "PDFFiles/fmc/" + FileOrgin + "/UnCompVsGenEtD" + type;
    c5->Update();
    c5->Print(plotname.c_str());
    delete c5;

}

void MCHistoComp(string FileOrgin)
{
    string type = ".png";
    string plotname;
    string saveloc = "fmcComp/" + FileOrgin;

    string hold = "fmc/" + FileOrgin + ".root";

    TFile* theFile = new TFile(hold.c_str());
    if(!theFile)
    {
        cout << "Error. Can not open file. try again" << endl;
        return;
    }
    TH1 *PsdEt, *PsdRtoEt;
    TH1 *PsdAllEt, *PsdAllRtoEt;
    string PsdEtLoc = "/CutAs/finebit/Et_of_Electron"; // The location of both psdEt and PsdRtoEt
    string PsdRtoEtLoc = "/CutAs/finebit/Et_Ratio_of_Electron";
    string PsdhitLoc = "/CusAs/finebit/Hits";
    string PsdAllEtLoc = "/CutAs/no_fine/Et_of_Electron"; // The location of both psdEt and PsdRtoEt
    string PsdAllRtoEtLoc = "/CutAs/no_fine/Et_Ratio_of_Electron";
    PsdEt = (TH1F*) theFile->Get(PsdEtLoc.c_str());
    cout << " Done opening " << PsdEtLoc << endl;
    PsdRtoEt = (TH1F*) theFile->Get(PsdRtoEtLoc.c_str());
    cout << " Done opening " << PsdRtoEtLoc << endl;

    PsdAllEt = (TH1F*) theFile->Get(PsdAllEtLoc.c_str());
    cout << " Done opening " << PsdAllEtLoc << endl;
    PsdAllRtoEt = (TH1F*) theFile->Get(PsdAllRtoEtLoc.c_str());
    cout << " Done opening " << PsdAllRtoEtLoc << endl;


    cout << "attempted directory name :" << "PDFFiles/fmc/" + FileOrgin << endl;
    gSystem->mkdir(("PDFFiles/" + saveloc).c_str());
    TCanvas* c1 = new TCanvas("c1", "", 800, 700);
    c1->SetFillColor(10);
    c1->SetFillColor(10);

    c1->cd();

    PsdEt->GetXaxis()->SetTitle("Et(Gev)");
    PsdEt->Draw();

    PsdAllEt->SetLineColor(8);
    PsdAllEt->Draw("same");
    //gStyle->SetOptStat("eoumi");
    gStyle->SetOptStat("");
    c1->SetLogy();
    /////////////////////
    TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg->SetHeader("Et");
    leg-> SetFillColor(kWhite);
    leg-> SetShadowColor(kWhite);
    leg->AddEntry(PsdEt, "Et using trigprimcuts");
    leg->AddEntry(PsdAllEt, "Et no trigprimcuts");
    leg->Draw();
    c1->Update();
    plotname = "PDFFiles/" + saveloc + "/Passed_Et" + type;

    c1->Print(plotname.c_str());
    delete c1;


    TCanvas* c2 = new TCanvas("c2", "", 800, 700);
    c2->SetFillColor(10);
    c2->SetFillColor(10);
    c2->cd();
    PsdRtoEt->GetXaxis()->SetTitle("Et/3X3");
    PsdRtoEt->Draw();
    PsdAllRtoEt->SetLineColor(8);
    PsdAllRtoEt->Draw("same");
    //gStyle->SetOptStat("eoumi");
    gStyle->SetOptStat("");
    c2->SetLogy();
    /////////////////////
    leg = new TLegend(0.25, 0.7, 0.55, 0.9);
    //leg->SetHeader("Et/3X3");
    leg-> SetFillColor(kWhite);
    leg-> SetShadowColor(kWhite);
    leg->AddEntry(PsdRtoEt, "Ratio using trigprimcuts");
    leg->AddEntry(PsdAllRtoEt, "Ratio no trigprimcuts");
    leg->Draw();
    c2->Update();
    plotname = "PDFFiles/" + saveloc + "/Passed_EtRatio" + type;

    c2->Print(plotname.c_str());
    delete c2;

}

void DataLumiComp(string FileOrgin)
{

    string type = ".png";
    string plotname;
    string saveloc = "fdata/" + FileOrgin + "/LumiTrig";
    cout << "test 0" << endl;
    FileOrgin = "fdata/" + FileOrgin;
    string hold = FileOrgin + ".root";
    cout << "test 1" << endl;
    TFile* theFile = new TFile(hold.c_str());
    if(!theFile)
    {
        cout << "Error. Can not open file. try again" << endl;
        return;
    }

    TH2 *PrimNoFineTrigCanLum, *PrimFineTrigCanLum;
    //TH2 *NormPrimNoFineTrigCanLum, *NormPrimFineTrigCanLum;
    TH2F NormPrimFineTrigCanLum("Luminosity_Vrs_Trig_Candidates", "Luminosity_Vrs_Trig_Candidates; Luminosity ;Trigger_Candidates", 30, 1900, 7900, 10, 0, 10);
    TH2F NormPrimNoFineTrigCanLum("Luminosity_Vrs_Trig_Candidates", "Luminosity_Vrs_Trig_Candidates; Luminosity ;Trigger_Candidates", 30, 1900, 7900, 10, 0, 10);
    string PrimNoFineTrigCanLumloc = "/DataAna/no_fine/Luminosity_Vrs_Trig_Candidates";
    string PrimFineTrigCanLumloc = "/DataAna/finebit/Luminosity_Vrs_Trig_Candidates";

    cout << "test 2" << endl;
    PrimNoFineTrigCanLum = (TH2F*) theFile->Get(PrimNoFineTrigCanLumloc.c_str());
    PrimFineTrigCanLum = (TH2F*) theFile->Get(PrimFineTrigCanLumloc.c_str());

    gSystem->mkdir(("PDFFiles/" + saveloc).c_str());


    TProfile *ProfNoFineTrigCalcLum = new TProfile("luminosity_vs_Trigger_Candidates", "luminosity_vs_Trigger_Candidates", 30, 1900, 7900, 0, 10);
    for(int i = 1; i < 30; i++)
    {
        for(int k = 0; k < 10; k++)
        {
            ProfNoFineTrigCalcLum->Fill(1700 + 200 * i, k - 1, PrimNoFineTrigCanLum->GetBinContent(i, k));
            //cout<<" okay this is odd :"<<PrimNoFineTrigCanLum->GetBinContent(i, k)<<endl;
        }
    }

    TCanvas* c1 = new TCanvas("c1", "", 900, 700);
    c1->SetFillColor(10);
    c1->SetFillColor(10);

    c1->cd();

    c1->SetLogy();
    ProfNoFineTrigCalcLum->GetYaxis()->SetTitle("luminosity Hz/Ub");
    ProfNoFineTrigCalcLum->GetYaxis()->SetTitleOffset(1.5);
    ProfNoFineTrigCalcLum->GetXaxis()->SetTitle("Number of Trigger Candidates");
    ProfNoFineTrigCalcLum->Draw();
    //PrimFineTrigCanLum->Draw("COLZ");
    cout << "test 4" << endl;
    plotname = "PDFFiles/" + saveloc + "/nofinetrig" + type;
    gStyle->SetOptStat("");

    c1->Update();
    c1->Print(plotname.c_str());
    delete c1;


    TCanvas* c2 = new TCanvas("c2", "", 900, 700);
    c2->SetFillColor(10);
    c2->SetFillColor(10);
    c2->cd();

    TProfile *ProfFineTrigCalcLum = new TProfile("luminosity_vs_Trigger_Candidates", "luminosity_vs_Trigger_Candidates", 30, 1900, 7900, -.01, 10.01);
    for(int i = 1; i < 30; i++)
    {
        for(int k = 0; k < 10; k++)
        {
            ProfFineTrigCalcLum->Fill(1700 + 200 * i, k - 1, PrimFineTrigCanLum->GetBinContent(i, k));
            //if(k==0)cout<<" okay this is odd :"<<PrimFineTrigCanLum->GetBinContent(i,k)<<endl;

        }
    }

    c2->SetLogy();
    ProfFineTrigCalcLum->GetYaxis()->SetTitle("luminosity Hz/Ub");
    ProfFineTrigCalcLum->GetYaxis()->SetTitleOffset(1.5);
    ProfFineTrigCalcLum->GetXaxis()->SetTitle("Number of Trigger Candidates");
    ProfFineTrigCalcLum->Draw();
    cout << "test 6" << endl;
    plotname = "PDFFiles/" + saveloc + "/finetrig" + type;
    //gStyle->SetOptStat("");
    c2->Update();
    c2->Print(plotname.c_str());
    delete c2;

    for(int i = 1; i <= 30; i++)
    {
        float lumeventf = 0;
        float lumeventnf = 0;
        for(int k = 1; k <= 10; k++)
        {
            lumeventnf += PrimNoFineTrigCanLum->GetBinContent(i, k);
            lumeventf += PrimFineTrigCanLum->GetBinContent(i, k);
        }
        //cout<<"lumeventf :"<<lumeventf<<endl;


        for(int k = 1; k <= 10; k++)
        {
            if(lumeventnf != 0)
            {
                if(k < 10 && i == 5 && false)
                {
                    cout << "lumeventf :" << lumeventnf << endl;
                    cout << "bin content :" << PrimNoFineTrigCanLum->GetBinContent(i, k) << endl;
                    cout << "and number that is going into it" << PrimNoFineTrigCanLum->GetBinContent(i, k) / lumeventnf << endl;
                }
                NormPrimNoFineTrigCanLum.SetBinContent(i, k, (PrimNoFineTrigCanLum->GetBinContent(i, k) / lumeventnf)*50000);

            } else NormPrimNoFineTrigCanLum.SetBinContent(i, k, 0);
            if(lumeventf != 0) NormPrimFineTrigCanLum.SetBinContent(i, k, (PrimFineTrigCanLum->GetBinContent(i, k) / lumeventf)*50000);
            else NormPrimFineTrigCanLum.SetBinContent(i, k, 0);
        }
    }

    cout << "test 7" << endl;
    TCanvas* c3 = new TCanvas("c3", "", 900, 700);
    c3->SetFillColor(10);
    c3->SetFillColor(10);
    c3->cd();
    c3->SetLogz();
    NormPrimFineTrigCanLum.GetYaxis()->SetTitle("luminosity Hz/Ub");
    NormPrimFineTrigCanLum.GetYaxis()->SetTitleOffset(1.5);
    NormPrimFineTrigCanLum.GetXaxis()->SetTitle("Number of Trigger Candidates");
    NormPrimFineTrigCanLum.Draw("COLZ");
    plotname = "PDFFiles/" + saveloc + "/finetrignorm" + type;
    gStyle->SetOptStat("");
    c3->Update();
    c3->Print(plotname.c_str());
    delete c3;



    TCanvas* c4 = new TCanvas("c3", "", 900, 700);
    c4->SetFillColor(10);
    c4->SetFillColor(10);
    c4->cd();
    c4->SetLogz();
    NormPrimNoFineTrigCanLum.GetYaxis()->SetTitle("luminosity Hz/Ub");
    NormPrimNoFineTrigCanLum.GetYaxis()->SetTitleOffset(1.5);
    NormPrimNoFineTrigCanLum.GetXaxis()->SetTitle("Number of Trigger Candidates");
    //NormPrimNoFineTrigCanLum.GetXaxis()->SetTimeOffset(.15);
    NormPrimNoFineTrigCanLum.Draw("COLZ");
    plotname = "PDFFiles/" + saveloc + "/nofinetrignorm" + type;
    gStyle->SetOptStat("");
    c4->Update();
    c4->Print(plotname.c_str());
    delete c4;


    TCanvas* c5 = new TCanvas("c5", "", 900, 700);
    c5->SetFillColor(10);
    c5->SetFillColor(10);
    c5->cd();
    c5->SetLogz();
    

    c5->SetLogy();
    ProfFineTrigCalcLum->GetYaxis()->SetTitle("luminosity Hz/Ub");
    ProfFineTrigCalcLum->GetYaxis()->SetTitleOffset(1.5);
    ProfFineTrigCalcLum->GetXaxis()->SetTitle("Number of Trigger Candidates");
    ProfFineTrigCalcLum->SetMarkerColor(kBlue+1);
    gStyle->SetOptStat("");
    ProfFineTrigCalcLum->SetMarkerStyle(kFullSquare);
    ProfFineTrigCalcLum->Draw();
    ProfNoFineTrigCalcLum->SetMarkerStyle(kStar);
    ProfNoFineTrigCalcLum->SetMarkerColor(kGreen +2);
    ProfNoFineTrigCalcLum->Draw("same");
    
    
    TLegend* leg = new TLegend(0.15, 0.8, 0.3, 0.9);
    //leg->SetHeader("Et");
    leg-> SetFillColor(kWhite);
    leg-> SetShadowColor(kWhite);
    leg->AddEntry(ProfFineTrigCalcLum, "Fine Grain");
    leg->AddEntry(ProfNoFineTrigCalcLum, "No Fine Gran");
    leg->Draw();
    c5->Update();
    plotname = "PDFFiles/" + saveloc + "/CompLumfine" + type;
    c5->Print(plotname.c_str());
    delete c5;
    
    
}

void MCVrsData(string FileOrgin)
{
    MCHisto(FileOrgin);
    DataHistoGram(FileOrgin);
    DataLumiComp(FileOrgin);
    MCHistoComp(FileOrgin);
    DataCompareHistoGram(FileOrgin);


    string type = ".png";
    string plotname;

    string BackLoc = "fdata/" + FileOrgin;
    string ElectLoc = "fmc/" + FileOrgin;
    string hold = BackLoc + ".root";
    string saveloc = "Comp/" + FileOrgin;
    TFile* theFile = new TFile(hold.c_str());
    if(!theFile)
    {
        cout << "Error. Can not open file. try again" << endl;
        return;
    }

    TH1 *PsdEtNon, *PsdRtoEtNon; //passed Et and Passed ratio of a 3x3 square to the location
    TH1 *PsdEtEle, *PsdRtoEtEle;
    //TH2 *PsdEtvrsiEta, *PsdEtvrsiPhi;
    //We will figure this out later terrible way of doing this for now
    string PsdEtNonLoc = "/DataAna/finebit/Et"; // The location of both psdEt and PsdRtoEt
    string PsdRtoEtNonLoc = "/DataAna/finebit/EtRatio";

    PsdEtNon = (TH1F*) theFile->Get(PsdEtNonLoc.c_str());
    cout << " Done opening " << PsdEtNonLoc << endl;
    PsdRtoEtNon = (TH1F*) theFile->Get(PsdRtoEtNonLoc.c_str());
    cout << " Done opening " << PsdRtoEtNonLoc << endl;
    gSystem->mkdir(("PDFFiles/" + saveloc).c_str());


    theFile = new TFile((ElectLoc + ".root").c_str());
    string PsdEtEleLoc = "/CutAs/finebit/Et_of_Electron"; // The location of both psdEt and PsdRtoEt
    string PsdRtoEtEleLoc = "/CutAs/finebit/Et_Ratio_of_Electron";

    PsdEtEle = (TH1F*) theFile->Get(PsdEtEleLoc.c_str());
    cout << " Done opening " << PsdEtEleLoc << endl;
    PsdRtoEtEle = (TH1F*) theFile->Get(PsdRtoEtEleLoc.c_str());
    cout << " Done opening " << PsdRtoEtEleLoc << endl;

    PsdEtNon->Scale(1 / PsdEtNon->GetEntries());

    PsdRtoEtNon->Scale(1 / PsdRtoEtNon->GetEntries());

    PsdEtEle->Scale(1 / PsdEtEle->GetEntries());

    PsdRtoEtEle->Scale(1 / PsdEtEle->GetEntries());

    TCanvas* c1 = new TCanvas("c1", "", 800, 700);
    c1->SetFillColor(10);
    c1->SetFillColor(10);

    c1->cd();
    PsdEtEle->GetXaxis()->SetTitle("Et(Gev)");

    float max1 = PsdEtNon->GetBin(PsdEtNon->GetMaximumBin());

    float max2 = PsdEtEle->GetBin(PsdEtEle->GetMaximumBin());
    //PsdEtNon->SetMaximum(1.1 * max(max1, max2));
    PsdEtNon->SetLineColor(kGreen + 1);
    PsdEtNon->Draw();
    PsdEtEle->SetLineColor(kBlue + 2);
    PsdEtEle->Draw("same");

    gStyle->SetOptStat("");
    //c1->SetLogy();
    //b///////////////////
    TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg->SetHeader("Et");
    leg-> SetFillColor(kWhite);
    leg-> SetShadowColor(kWhite);
    leg->AddEntry(PsdEtNon, "Et of Backround");
    leg->AddEntry(PsdEtEle, "Et of electrons");
    leg->Draw();
    c1->Update();
    plotname = "PDFFiles/" + saveloc + "/Passed_Et" + type;

    c1->Print(plotname.c_str());
    delete c1;



    TCanvas* c2 = new TCanvas("c2", "", 800, 700);
    c2->SetFillColor(10);
    c2->SetFillColor(10);

    c2->cd();

    max1 = PsdRtoEtNon->GetBin(PsdRtoEtNon->GetMaximumBin());

    max2 = PsdRtoEtEle->GetBin(PsdRtoEtEle->GetMaximumBin());
    PsdRtoEtNon->SetMaximum(1.1 * max(max1, max2));
    PsdRtoEtNon->SetLineColor(kGreen + 1);
    PsdRtoEtNon->GetXaxis()->SetTitle("Et/3X3");
    PsdRtoEtNon->Draw();
    PsdRtoEtEle->SetLineColor(kBlue + 2);
    PsdRtoEtEle->Draw("same");

    //gStyle->SetOptStat("eoumi");
    gStyle->SetOptStat("");
    c2->SetLogy();
    /////////////////////
    leg = new TLegend(0.25, 0.7, 0.55, 0.9);
    leg->SetHeader("Et/3X3");
    leg-> SetFillColor(kWhite);
    leg-> SetShadowColor(kWhite);
    leg->AddEntry(PsdRtoEtEle, "Ratio Electrons");
    leg->AddEntry(PsdRtoEtNon, "Ratio Background");
    leg->Draw();
    c2->Update();
    plotname = "PDFFiles/" + saveloc + "/EtRatio" + type;

    c2->Print(plotname.c_str());
    delete c2;



}



