#ifndef FITMCORRSTEPAN_H
#define FITMCORRSTEPAN_H

#include <string>

#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
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
 * 
*/

class FitMCorrStepan
{
public:
    FitMCorrStepan(Float_t eBeam, TClasTool *ct);
    ~FitMCorrStepan();
	void fillHists(TIdentificator *t);
    void fillExtraHists(TIdentificator *t);
	void fitHists();
	void writeParams();
private:
	Float_t eBeam;
	Float_t kMassP = 0.938;
	
    Int_t nEntries;
    Int_t nRowsEVNT;
    Int_t nSect = 6;
	Int_t thetaMin[4] = {0, 16, 20, 25};
	Int_t thetaMax[4] = {16, 20, 25, 90};
	Int_t phiMin[6] = {-30, 30, 90, 150, 210, 270};
	Int_t phiMax[6] = {30, 90, 150, 210, 270, 330};

    TString rootfName;
	
	TH1F *hW;
	TH1F *hZ;
    TH1F *hF1Mean[6];
    TH1F *hF2Mean[6];
    TH2F *hF1[6];
    TH2F *hF2[6];

    TH2F *hF1Extra[6][4];

	TEVNTClass *fEVNT;
    TClasTool *fCT;
    TIdentificator *fId;

    Bool_t isInitElec(TIdentificator *t);
    Bool_t isPartProton(TIdentificator *t, Int_t j);
    Float_t ratioF1(Float_t p, Float_t theta);
    Float_t pCalc(Float_t theta);
}

#endif // FITMCORRSTEPAN_H

