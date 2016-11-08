#include "PCorrSimple.h"
#include <cmath>

using namespace std;

PCorrSimple::PCorrSimple(Float_t eBeam, TClasTool *ct)
{
    nSect = 6;
    kMassP = 0.938;
    phiMin[0] = -30;
    phiMax[0] = 30;
    for(Int_t i = 1; i < 6; i++){
        phiMin[i] = phiMin[i-1] + 60;
        phiMax[i] = phiMax[i-1] + 60;
    }

    this->fCT = ct;
    this->eBeam = eBeam;
    this->nEntries = (Int_t) fCT->GetEntries();

    fId = new TIdentificator(fCT);
}

PCorrSimple::~PCorrSimple()
{
    delete fCT;
    delete fId;
}

void PCorrSimple::mainFun(TString rootfName, Float_t thCut)
{
    Bool_t hasProton;
    this->thCut = thCut;
    this->rootfName = rootfName;
    TFile *f = new TFile(rootfName, "RECREATE");

    TH1F *hPorc = new TH1F("hPorc", "", 100, 0.9, 1.1);
    TH1F *hW = new TH1F("hW", "", 100, 0.7, 1.2);
    TH1F *hnewPe = new TH1F("hnewPe", "", 500, 0, eBeam);
    TH1F *hPe = new TH1F("hPe", "", 500, 0, eBeam);

    Int_t sec;
    Float_t newP;

    LoadParams("fitParams.txt");
    LoadParamsE("fitPower2Params.txt");

    cout << "Number of Events: " << nEntries << endl;
    for(Int_t i = 0; i < nEntries; i++){
        fCT->Next();
        nRowsEVNT = fCT->GetNRows("EVNT");
        hasProton = false;
        if(nRowsEVNT > 1 && isInitElec() && fId->ThetaLab(0) >= thCut && fEVNT->GetZ() < -10 && fId->W() >= 0.8 && fId->W() <= 1.05){
            fEVNT = (TEVNTClass*) fCT->GetBankRow("EVNT", 0);
            for(Int_t j = 0; j < nRowsEVNT; j++){
                hasProton = isPartProton(j);
                if(hasProton) break;
            }
            if(hasProton){
                sec = fId->Sector(0);
                hPorc->Fill(f1(fId->PhiLab(0), sec)*f2e(fId->ThetaLab(0), sec));
                hPe->Fill(fId->Momentum(0));
                newP = fId->Momentum(0)*f1(fId->PhiLab(0), sec)*f2e(fId->ThetaLab(0), sec);
                hnewPe->Fill(newP);
                hW->Fill(newW(newQ2(newP), newNu(newP)));
            }
        }
    }

    f->cd();
    hPorc->Write();
    hW->Write();
    hnewPe->Write();
    hPe->Write();

    hPorc->Delete();
    hW->Delete();
    hnewPe->Delete();
    hPe->Delete();

    f->Close();
    delete f;
}

Bool_t PCorrSimple::isInitElec()
{
    Bool_t result;
    result = fId->Status(0) > 0 && fId->Status(0) < 100 && fId->Charge(0) == -1 && fId->StatCC(0) > 0 && fId->StatSC(0) > 0;
    result = result && fId->StatDC(0) > 0 && fId->StatEC(0) > 0 && fId->DCStatus(0) > 0 && fId->SCStatus(0) == 33 && fId->Nphe(0) > 25;
    result = result && (fId->Etot(0) / 0.27 / 1.15 + 0.4) > fId->Momentum(0) && (fId->Etot(0) / 0.27 / 1.15 - 0.2 < fId->Momentum(0));
    result = result && (fId->Ein(0) + fId->Eout(0) > 0.8 * 0.27 * fId->Momentum(0)) && (fId->Ein(0) + fId->Eout(0) < 1.2 * 0.27 * fId->Momentum(0));
    result = result && fId->Eout(0) != 0 && fId->FidCheckCut() == 1;
    return result;
}

Bool_t PCorrSimple::isPartProton(Int_t j)
{
    Bool_t result;
    result = j > 0 && fId->Charge(j) == 1 && fId->Status(j) > 0 && fId->Status(j) < 100 && fId->StatDC(j) > 0 && fId->DCStatus(j) > 0;
    result = result && fId->StatSC(j) > 0;
    if(fId->Momentum(j) > 1){
        result = result && fId->TimeCorr4(kMassP,j) >= -0.69 && fId->TimeCorr4(kMassP,j) <= 1.38;
    }
    else if(fId->Momentum(j) < 1){
        result = result && fId->TimeCorr4(kMassP,j) >= -3.78 && fId->TimeCorr4(kMassP,j) <= 6.75;
    }
    return result;
}

Float_t PCorrSimple::f1(Float_t phi, Int_t s)
{
    return f1Params[s][0]+f1Params[s][1]*phi+f1Params[s][2]*phi*phi;
}

Float_t PCorrSimple::f2(Float_t th, Int_t s)
{
    return f2Params[s][0]+(f2Params[s][1]+f2Params[s][1]*th+f2Params[s][2]*th*th)*TMath::Exp(-th);
}

Float_t PCorrSimple::f2e(Float_t th, Int_t s)
{
    return f2eParams[s][0]*pow(th, f2eParams[s][1])+f2eParams[s][2];
}

void PCorrSimple::LoadParams(TString fName)
{
    ifstream txtFile;
    TString aux;
    cout << "Reading " << fName << endl;
    txtFile.open(fName);
    for(Int_t i = 0; i < nSect; i++){
        txtFile >> aux >> aux >> aux >> f1Params[i][0] >> f1Params[i][1] >> f1Params[i][2];
        txtFile >> aux >> aux >> aux >> f2Params[i][0] >> f2Params[i][1] >> f2Params[i][2] >> f2Params[i][3];
    }
    txtFile.close();

    return;
}

void PCorrSimple::LoadParamsE(TString fName)
{
    ifstream txtFile;
    cout << "Reading " << fName << endl;
    txtFile.open(fName);
    for(Int_t i = 0; i < nSect; i++){
        txtFile >> f2eParams[i][0] >> f2eParams[i][1] >> f2eParams[i][2];
    }
    txtFile.close();

    return;
}

Float_t PCorrSimple::newNu(Float_t p)
{
    return eBeam - p;
}

Float_t PCorrSimple::newQ2(Float_t p)
{
    return 4*eBeam*p*sin(fId->ThetaLab(0)*TMath::Pi()/180./2)*sin(fId->ThetaLab(0)*TMath::Pi()/180./2);
}

Float_t PCorrSimple::newW(Float_t q2, Float_t nu)
{
    return TMath::Sqrt(kMassP * kMassP + 2. * kMassP * nu - q2);
}
