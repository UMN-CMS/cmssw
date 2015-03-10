//C++ include
#include <iostream>
//ROOT include
#include "TH2D.h"
#include "TH1D.h"
//CMSSW include
#include "TestEDA/CutAssessment/interface/Genconverter.h"
#include "DataFormats/Math/interface/deltaR.h"

//NOTE: I don't know how to run this code :( a parasitic version can be found in JetAlgorithm.cc:q

using namespace std;

void deltaRCheck()
{

    int cand1 [16] = {1, 1, 1, 1, 2, 2, 2, 2, 72, 72, 72, 72, 71, 71, 71, 71};
    int cand2 [16] = {1, 2, 72, 71, 1, 2, 72, 71, 1, 2, 72, 71, 1, 2, 72, 71};
    double dR [16];
    double etaVal = 3.6;
    
    for(unsigned int cand = 0; cand < 16; ++cand)
    {
        dR [cand] = deltaR(etaVal, Genconverter().IPhi2Phi(cand1[cand]), etaVal, Genconverter().IPhi2Phi(cand2[cand]));
        std::cout << dR [cand] << std::endl;

    }


























}
