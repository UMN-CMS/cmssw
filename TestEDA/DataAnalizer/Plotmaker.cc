#include <TLegend.h>
#include <TDirectory.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <math.h> 
#include <map>
#include "TH1.h"
#include "TH2.h"
#include <string>
#include <TFile.h>
#include <iostream>



#include <TMultiGraph.h>
#include <TGraphErrors.h>

#include <cmath> 
#define PI 3.14159265
using namespace std;

void MakeHistoGram(string FileOrgin)
{
    string type = ".pdf";
    string plotname;
    cout<<"test 0"<<endl;
    string hold=FileOrgin+".root";
    cout<<hold<<" please"<<endl;
    TFile* theFile = new TFile(hold.c_str());
    if(!theFile)
    {
        cout << "Error. Can not open file. try again" << endl;
        return;
    }
    cout<<"test 1"<<endl;
    TH1 *PsdEt, *PsdRtoEt, *PsdHit; //passed Et and Passed ratio of a 3x3 square to the location
    //TH2 *PsdEtvrsiEta, *PsdEtvrsiPhi;
    //We will figure this out later terrible way of doing this for now
    string PsdEtLoc="fdata/"+FileOrgin+"/DataAna/Et"; // The location of both psdEt and PsdRtoEt
    string PsdRtoEtLoc="fdata/"+FileOrgin+"/DataAna/EtRatio";
    string PsdhitLoc="fdata/"+FileOrgin+Data+"/DataAna/Hits";
    string PsdEtvrsiEtaLoc;
    string PsdEtvrsiPhiLoc;

    PsdEt = (TH1F*) theFile->Get(PsdEtLoc.c_str());
    cout << " Done opening " << PsdEtLoc << endl;
    PsdRtoEt = (TH1F*) theFile->Get(PsdRtoEtLoc.c_str());
    cout << " Done opening " << PsdRtoEtLoc << endl;
    PsdHit =(TH1I*) theFile->Get(PsdhitLoc.c_str());
    cout<<"Done opening "<<PsdhitLoc<<endl;
    gSystem->mkdir("PDFFiles"+FileOrgin)
    TCanvas* c1 = new TCanvas("c1", "", 800, 700);
    c1->SetFillColor(10);
    c1->SetFillColor(10);

    c1->cd();

    PsdEt->Draw();
    gStyle->SetOptStat("eoumi");
    gStyle->SetOptStat("");
    c1->SetLogy();
/////////////////////

    c1->Update();
    plotname ="PDFFiles/" +FileOrgin+"/Passed_Et" type;
    
    c1->Print(plotname.c_str());
    delete c1;
    
    TCanvas* c2 = new TCanvas("c2", "", 800, 700);
    
    
    c2->SetFillColor(10);
    c2->SetFillColor(10);
    c2->cd();
    PsdRtoEt->Draw();

    gStyle->SetOptStat("eoumi");
    gStyle->SetOptStat("");

    c2->SetLogy();
/////////////////////

    c2->Update();
    plotname = "PDFFiles" +FileOrgin+"/Passed_Et_ratio"+ type;
    c2->Print(plotname.c_str());
    delete c2;
    
     TCanvas* c3 = new TCanvas("c3", "", 800, 700);
     
    c3->SetFillColor(10);
    c3->SetFillColor(10);
    c3->cd();
    PsdHit->Draw();
    
    gStyle->SetOptStat("eoumi");
    gStyle->SetOptStat("");

    c3->SetLogy();
    
    c3->Update();
    plotname ="PDFFiles" +FileOrgin+"/Hits" type;
    
    c3->Print(plotname.c_str());
    delete c3;
}

void MCHisto()
{
    
}