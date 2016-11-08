#ifndef PCORRSIMPLE_H
#define PCORRSIMPLE_H

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

class PCorrSimple{
public:
    PCorrSimple(Float_t eBeam, TClasTool *ct);
    ~PCorrSimple();
    void mainFun(TString rootfName, Float_t thCut);
private:
    Int_t nEntries;
    Int_t nRowsEVNT;
    Int_t nSect;
    Float_t eBeam;
    Float_t thCut;
    Float_t kMassP;
    Float_t phiMin[6];
    Float_t phiMax[6];
    Float_t f1Params[6][3];
    Float_t f2Params[6][4];
    Float_t f2eParams[6][3];
    TString rootfName;

    TEVNTClass *fEVNT;
    TClasTool *fCT;
    TIdentificator *fId;

    Bool_t isInitElec();
    Bool_t isPartProton(Int_t j);
    Float_t f1(Float_t phi, Int_t s);
    Float_t f2(Float_t th, Int_t s);
    Float_t f2e(Float_t th, Int_t s);
    Float_t newQ2(Float_t p);
    Float_t newNu(Float_t p);
    Float_t newW(Float_t q2, Float_t nu);
    void LoadParams(TString fName);
    void LoadParamsE(TString fName);
};

#endif // PCORRSIMPLE_H
