#include "fitPCorrStepan.h"

FitMCorrStepan::FitMCorrStepan(Float_t eBeam, TClasTool *ct)
{
    this->fCT = ct;
    this->eBeam = eBeam;
    this->nEntries = (Int_t) fCT->GetEntries();

    fId = new TIdentificator(fCT);
}

FitMCorrStepan::~FitMCorrStepan()
{
    delete hW;
    delete hZ;

    for(Int_t i = 0; i < nSect; i++){
        delete hF1[i];
        delete hF2[i];
    }
}

void FitMCorrStepan::fillHists(TString rootfName)
{
    this->rootfName = rootfName;
    TFile *f = new TFile(rootfName, "RECREATE");

    Bool_t hasProton;

    Int_t phiBins = 60;
    Int_t thBins = 20;
    Int_t f1Bins = 20;
    Int_t f2Bins = 20;

    Float_t f1Min = 0.9;
    Float_t f1Max = 1.1;
    Float_t f2Min = 0.9;
    Float_t f2Max = 1.1;
    Float_t thMin = 15;
    Float_t thMax = 35;

    hW = new TH1F("hW", "W no DIS", 50, 0.7, 1.2);
    hZ = new TH1F("hZ", "Z Coordinate", 100, -80, 20);

    TString sHF1Title;
    TString sHF2Title;
    TString sHF1MeanTitle;
    TString sHF2MeanTitle;

    for(Int_t i = 0; i < nSect; i++){
        sHF1Title = Form("F1 vs Phi (Sector %d)", i+1);
        sHF2Title = Form("F2 vs Theta (Sector %d)", i+1);
        hF1[i] = new TH2F(Form("hF1_phi%d", i+1), sHF1Title, phiBins, phiMin[i], phiMax[i], f1Bins, f1Min, f1Max);
        hF2[i] = new TH2F(Form("hF2_phi%d", i+1), sHF2Title, thBins, thMin, thMax, f2Bins, f2Min, f2Max);
        hF1Mean[i] = new TH1F(Form("hF1_phi%d_m", i+1), sHF1MeanTitle, phiBins, phiMin[i], phiMax[i]);
        hF2Mean[i] = new TH1F(Form("hF2_phi%d_m", i+1), sHF2MeanTitle, phiBins, phiMin[i], phiMax[i]);
    }

    for(Int_t i = 0; i < nEntries; i++){
        fCT->Next();
        nRowsEVNT = fCT->GetNRows("EVNT");
        hasProton = false;
        if(nRowsEVNT > 1 && isInitElec(t) && fId->W() >= 0.8 && fId->W() <= 1.05 && fId->ThetaLab(0) >= 16){
            for(Int_t j = 0; j < nRowsEVNT; j++){
                if(hasProton = isPartProton(fId, j)) break;
            }
        }
        if(hasProton){
            hZ->Fill(fEVNT->GetZ());
            hW->Fill(fId->W());
            Int_t sec = fId->Sector(0);
            hF1[sec]->Fill(fId->PhiLab(0), ratioF1(fId->Momentum(0), fId->ThetaLab(0)));
            hF2[sec]->Fill(fId->ThetaLab(0), ratioF1(fId->Momentum(0), fId->ThetaLab(0)));
        }
    }

    Int_t binContent;
    Int_t totBin;

    TH1F *hGauss;
    TF1 *fGauss;

    for(Int_t i = 0; i < nSect; i++){
        for(Int_t j = 0; j < phiBins; j++){
            totBin = 0;
            for(Int_t k = 0; k < f1Bins; k++)
                totBin = totBin + hF1[i]->GetBinContent(j+1, k+1);
            if(totBin = 0) continue;
            hGauss = new TH1F("hGauss", "hGauss", f1Bins, f1Min, f1Max);
            for(Int_t k = 0; k < f1Bins; k++){
                binContent = hF1[i]->GetBinContent(j+1, k+1);
                if(binContent != 0){
                    hGauss->Fill(hF1[i]->GetXaxis()->GetBinCenter(k+1), binContent);
                }
            }
            hGauss->Fit("gaus");
            fGauss = hGauss->GetFunction("gaus");

        }
    }

    f->cd();
    hW->Write();
    hZ->Write();
    for(Int_t i = 0; i < nSect; i++){
        hF1[i]->Write();
        hF2[i]->Write();
    }
    f->Close();
}

void FitMCorrStepan::fillExtraHists(TString rootfName){

}

void FitMCorrStepan::fitHists()
{

}

void FitMCorrStepan::writeParams()
{

}

Bool_t FitMCorrStepan::isInitElec()
{
    Bool_t result;
    result = fId->Status(0) > 0 && fId->Status(0) < 100 && fId->Charge(0) == -1 && fId->StatCC(0) > 0 && fId->StatSC(0) > 0;
    result = result && fId->StatDC(0) > 0 && fId->StatEC(0) > 0 && fId->DCStatus(0) > 0 && fId->SCStatus(0) == 33 && fId->Nphe(0) > 25;
    result = result && (t->Etot(0) / 0.27 / 1.15 + 0.4) > fId->Momentum(0) && (t->Etot(0) / 0.27 / 1.15 - 0.2 < fId->Momentum(0));
    result = result && (t->Ein(0) + fId->Eout(0) > 0.8 * 0.27 * fId->Momentum(0)) && (t->Ein(0) + fId->Eout(0) < 1.2 * 0.27 * fId->Momentum(0));
    result = result && fId->Eout(0) != 0; // && fId->FidCheckCut() == 1;
    return result;
}

Float_t FitMCorrStepan::ratioF1(Float_t p, Float_t theta)
{
    Float_t ratio;
    ratio = pCalc(theta)/p;
    return ratio;
}

Float_t FitMCorrStepan::pCalc(Float_t theta){
    Float_t p;
    theta = TMath::Pi()*theta/180;
    p = eBeam/(1 + E_beam * (1 - TMath::Cos(theta))/kMassP);
    return p;
}
