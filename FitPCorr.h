#ifndef FITPCORR_H
#define FITPCORR_H

#include <string>
#include <iostream>
#include <fstream>

#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TClasTool.h"
#include "TEVNTClass.h"
#include "TIdentificator.h"

/*
 * Momentum correction is calculated based on the document
 * "Electron Momentum Correction for CLAS at 4.4 GeV" from Stepanyan
 * It selects events with the following conditions:
 * - 0.8 <= W <= 1.05
 * - theta >= 16
 * - One good electron and a proton
*/

class FitPCorr
{
public:
    FitPCorr(Float_t eBeam, TClasTool *ct);
    ~FitPCorr();
    void fillHists(TString rootfName, Bool_t printGFit, Float_t thCut);
    void fitHists(Bool_t setInitParams, Bool_t tryF2e);
    void writeParams(TString txtfName);
    void printPlots();
    void writeHistsBins(TString txtfName);
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
    Float_t phiMin[6];
    Float_t phiMax[6];
    Float_t thetaMin[4];
    Float_t thetaMax[4];
    Float_t f1Params[6][3]; // Parameters f1(x) = p1+p2*x+p2*x*x
    Float_t f2Params[6][4]; // Parameters f2(x) = p1+(p2+p3*x+p4*x*x)*exp(-x)
    Float_t f2eParams[6][3]; // Parameters f2e(x) = p1*x^p2+p3

    TString rootfName;
	TH1F *hW;
    TH1F *hWextra;
	TH1F *hZ;
    TH1F *hF1Mean[6];
    TH1F *hF2Mean[6];
    TH2F *hF1[6];
    TH2F *hF2[6];
    TH2F *hF1Extra[6][4];

	TEVNTClass *fEVNT;
    TClasTool *fCT;
    TIdentificator *fId;

    Bool_t isInitElec();
    Bool_t isPartProton(Int_t j);
    Float_t ratioF1(Float_t p, Float_t theta);
    Float_t pCalc(Float_t theta);
    void initHists();
    void writeHists(TFile *f);
};

#endif // FITPCORR_H
