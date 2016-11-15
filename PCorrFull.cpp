#include "PCorrFull.h"

using namespace std;

PCorrFull::PCorrFull(Float_t eBeam, TClasTool *ct, TString fName, Float_t thCut)
{
    nSect = 6;
    kMassP = 0.938;
    thetaMin[0] = 0;
    thetaMin[1] = 16;
    thetaMin[2] = 20;
    thetaMin[3] = 25;
    thetaMax[0] = 16;
    thetaMax[1] = 20;
    thetaMax[2] = 25;
    thetaMax[3] = 90;
    thMin = 16;
    thMax = 34;
    phiMin[0] = -30;
    phiMax[0] = 30;
    for(Int_t i = 1; i < nSect; i++){
        phiMin[i] = phiMin[i-1] + 60;
        phiMax[i] = phiMax[i-1] + 60;
    }

    phiBins = 60;
    thBins = 18;
    f1Bins = 20;
    f2Bins = 20;
    f1Min = 0.95;
    f1Max = 1.05;
    f2Min = 0.95;
    f2Max = 1.05;

    this->eBeam = eBeam;
    this->fCT = ct;
    this->thCut = thCut;
    this->nEntries = (Int_t) fCT->GetEntries();
    this->rootfName = fName;
    fId = new TIdentificator(fCT);
}

PCorrFull::~PCorrFull()
{

}

void PCorrFull::fillF1(Bool_t printGFit)
{
    Bool_t hasProton;
    TFile *f = new TFile(rootfName, "RECREATE");

    hW = new TH1F("hW", "W no DIS", 100, 0.7, 1.2);
    hWextra = new TH1F("hWextra", "W no DIS", 100, 0.7, 1.2);

    TString sHF1Title;
    TString sHF1MeanTitle;

    for(Int_t i = 0; i < nSect; i++){
        sHF1Title = Form("#phi vs f_{1} (Sector %d)", i+1);
        hF1[i] = new TH2F(Form("hF1_phi%d", i+1), sHF1Title, phiBins, phiMin[i], phiMax[i], f1Bins, f1Min, f1Max);
        sHF1MeanTitle = Form("#phi vs mean f_{1} (Sector %d)", i+1);
        hF1Mean[i] = new TH1F(Form("hF1_phi%d_m", i+1), sHF1MeanTitle, phiBins, phiMin[i], phiMax[i]);
    }

    cout << "Number of Events: " << nEntries << endl;

    for(Int_t i = 0; i < nEntries; i++){
        fCT->Next();
        nRowsEVNT = fCT->GetNRows("EVNT");
        hasProton = false;
        if(nRowsEVNT > 1 && isInitElec()){
            if(fId->ThetaLab(0) >= 16){
                for(Int_t j = 0; j < nRowsEVNT; j++){
                    hasProton = isPartProton(j);
                    if(hasProton) break;
                }
            }
            fEVNT = (TEVNTClass*) fCT->GetBankRow("EVNT", 0);
            if(hasProton && fEVNT->GetZ() < -10){
                hWextra->Fill(fId->W());
                if(fId->W() >= 0.8 && fId->W() <= 1.05){

                }
            }
        }
    }
}

void PCorrFull::fitF1()
{
    TFile *f = new TFile(rootfName, "UPDATE");
    TH1F *hF1r = NULL;
    TF1 *f1[nSect];

    for(Int_t i = 0; i < nSect; i++)
        f1[i] = new TF1(Form("f1_phi%d", i+1), "[0]+[1]*x+[2]*x*x", phiMin[i], phiMax[i]);

    for(Int_t i = 0; i < nSect; i++){
        hF1r = (TH1F*) f->Get(Form("hF1_phi%d_m", i+1));
        hF1r->Fit(Form("f1_phi%", i+1), "REQ");
        f->cd();
        hF1r->Write(Form("%_FIT", hF1r->GetName()));
        for(Int_t j = 0; j < 3; j++)
            f1Params[i][j] = f1[i]->GetParameter(j);
    }

    for(Int_t i = 0; i < nSect; i++)
        f1[i]->Delete();
    hF1r->Delete();
    f->Close();
    delete f;
}

void PCorrFull::fillF2(Bool_t printGFit)
{
    Bool_t hasProton;
    TFile *f = new TFile(rootfName, "UPDATE");
    hW = new TH1F("hW_f1", "W no DIS", 100, 0.7, 1.2);
    TH1F *hF1Perc = new TH1F("hF1Perc", "f1 Correction", 60, 0.97, 1.03);

    TString sHF2Title;
    TString sHF2MeanTitle;
    for(Int_t i = 0; i < nSect; i++){

    }
}

