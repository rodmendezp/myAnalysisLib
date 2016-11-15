#ifndef PCORRFULL_H
#define PCORRFULL_H

#include <iostream>
#include <fstream>
#include <string>

#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TClasTool.h"
#include "TEVNTClass.h"
#include "TIdentificator.h"

class PCorrFull{
public:
    PCorrFull(Float_t eBeam, TClasTool *ct, TString fName, Float_t thCut = 16);
    ~PCorrFull();
    void fillF1(Bool_t printGFit);
    void fitF1();
    void fillF2(Bool_t printGFit);
    void fitF2();
    void writeParams(TString txtfName);
    void printPlots();
private:
    Int_t nEntries;
    Int_t nRowsEVNT;
    Int_t nSect;
    Int_t phiBins;
    Int_t thBins;
    Int_t f1Bins;
    Int_t f2Bins;
    Float_t eBeam;
    Float_t kMassP;
    Float_t thCut;
    Float_t thMin;
    Float_t thMax;
    Float_t f1Min;
    Float_t f1Max;
    Float_t f2Min;
    Float_t f2Max;
    Float_t f1Params[6][3]; // Parameters f1(x) = p1+p2*x+p2*x*x
    Float_t f2Params[6][4]; // Parameters f2(x) = p1+(p2+p3*x+p4*x*x)*exp(-x)
    Float_t phiMin[6];
    Float_t phiMax[6];

    TString rootfName;
    TH1F *hW;
    TH1F *hWextra;
    TH2F *hF1[6];
    TH2F *hF2[6];
    TH1F *hF1Mean[6];
    TH1F *hF2Mean[6];
    TEVNTClass *fEVNT;
    TClasTool *fCT;
    TIdentificator *fId;

    Bool_t isInitElec();
    Bool_t isPartProton(Int_t j);
};

#endif // FITPCORR_H
