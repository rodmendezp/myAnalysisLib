#include "pCorrStepan.h"
#include "myROOTUtils.h"
#include "TTree.h"
#include "TMath.h"
#include "TClonesArray.h"

#define __COMPRESS__ 9
#define __BUFSIZE__  16384
#define PHYSICEVNT 2
#define __SCALER_BUFSIZE__ 1024
#define __COMPRESS__ 9
#define __MAXFILESIZE__ 2000000000

using namespace std;

pCorrStepan::pCorrStepan(string txtfName){
    TString aux;
    nSect = 6;
    ifstream txtFile;
    if(fexists(txtfName)){
        txtFile.open(txtfName.c_str());
        for(Int_t i = 0; i < nSect; i++){
            txtFile >> aux >> aux >> aux >> f1Params[i][0] >> f1Params[i][1] >> f1Params[i][2];
            txtFile >> aux >> aux >> aux >> f2Params[i][0] >> f2Params[i][1] >> f2Params[i][2] >> f2Params[i][3];
        }
        txtFile.close();
    }
}

pCorrStepan::~pCorrStepan(){

}

void pCorrStepan::pCorrFile(string rootfName){
    Int_t comp = __COMPRESS__;
    Int_t nEntries;
    Int_t sect;

    Float_t nMom;
    Float_t nPx;
    Float_t nPy;
    Float_t nPxy;
    Float_t nPz;

    string nrootfName = Form("%s.pCorr", rootfName.c_str());

    Last_Evt_Num = 0;
    Last_Key_Num = 0;

    TFile *f = new TFile(nrootfName.c_str(), "RECREATE", "Title", comp);
    f->SetCompressionLevel(comp);

    hHEADER = new THEADERClass();

    initClones();

    tree = new TTree("CLASEVENT","CLAS Event Tree");
    tree->SetAutoSave(1024*1024*1024);

    initBranches(tree);

    fCT = new TClasTool();
    fCT->InitDSTReader("ROOTDSRT");
    if(fexists(rootfName))
        fCT->Add(rootfName.c_str());

    fID = new TIdentificator(fCT);
    nEntries = (Int_t) fCT->GetEntries();

    int EventType = PHYSICEVNT;

    cout << "Number of Events = " << nEntries << endl;

    for(Int_t i = 0; i < nEntries; i++){
        clearClones();
        nEVNT_Store = 0;
        nECPB_Store = 0;
        nSCPB_Store = 0;
        nDCPB_Store = 0;
        nCCPB_Store = 0;
        nLCPB_Store = 0;
        nSTPB_Store = 0;
        nTGPB_Store = 0;
        nTAGR_Store = 0;
        if((i+1)%100 == 0) cout << "Filling Event " << i + 1 << endl;
        if(i == 10) break;
        fCT->Next();
        hHEADER = (THEADERClass*) fCT->GetHEADER();
        Last_Evt_Num = hHEADER->NEvent;
        for(Int_t j = 0; j < (Int_t) fCT->GetNRows("EVNT"); j++){
            fEVNT = (TEVNTClass*) fCT->GetBankRow("EVNT", j);
            nfEVNT = new TEVNTClass(fEVNT);
            sect = fID->Sector(j);
            nMom = fID->Momentum(j)*calcF1(sect, fID->PhiLab(j))*calcF2(sect, fID->ThetaLab(j));
            cout << "Phi = " << fID->PhiLab(j) << " Theta = " << fID->ThetaLab(j) << endl;
            cout << "Old Momentum = " << fID->Momentum(j) << " New Momentum " << nMom << endl;
            cout << "Ratio = " << nMom/fID->Momentum(j) << endl;
            nPz = nMom*TMath::Cos(fID->ThetaLab(j)*TMath::Pi()/180);
            nPxy = nMom*TMath::Sin(fID->ThetaLab(j)*TMath::Pi()/180);
            nPx = nPxy*TMath::Cos(fID->PhiLab(j)*TMath::Pi()/180);
            nPy = nPxy*TMath::Sin(fID->PhiLab(j)*TMath::Pi()/180);
            cout << "Old Px = " << fEVNT->Px << " New Px = " << nPx << endl;
            cout << "Old Py = " << fEVNT->Py << " New Py = " << nPy << endl;
            cout << "Old Pz = " << fEVNT->Pz << " New Pz = " << nPz << endl;
            nfEVNT->Px = nPx;
            nfEVNT->Py = nPy;
            nfEVNT->Pz = nPz;
            TClonesArray &tEVNTbank = *EVNTStore;
            new(tEVNTbank[nEVNT_Store++]) TEVNTClass(nfEVNT);
        }
        for(Int_t j = 0; j < (Int_t) fCT->GetNRows("ECPB"); j++){
            fECPB = (TECPBClass*) fCT->GetBankRow("ECPB", j);
            TClonesArray &tECPBbank = *ECPBStore;
            new(tECPBbank[nECPB_Store++]) TECPBClass(fECPB);
        }
        for(Int_t j = 0; j < (Int_t) fCT->GetNRows("SCPB"); j++){
            fSCPB = (TSCPBClass*) fCT->GetBankRow("SCPB", j);
            TClonesArray &tSCPBbank = *SCPBStore;
            new(tSCPBbank[nSCPB_Store++]) TSCPBClass(fSCPB);
        }
        for(Int_t j = 0; j < (Int_t) fCT->GetNRows("DCPB"); j++){
            fDCPB = (TDCPBClass*) fCT->GetBankRow("DCPB", j);
            TClonesArray &tDCPBbank = *DCPBStore;
            new(tDCPBbank[nDCPB_Store++]) TDCPBClass(fDCPB);
        }
        for(Int_t j = 0; j < (Int_t) fCT->GetNRows("CCPB"); j++){
            fCCPB = (TCCPBClass*) fCT->GetBankRow("CCPB", j);
            TClonesArray &tCCPBbank = *CCPBStore;
            new(tCCPBbank[nCCPB_Store++]) TCCPBClass(fCCPB);
        }
        for(Int_t j = 0; j < (Int_t) fCT->GetNRows("LCPB"); j++){
            fLCPB = (TLCPBClass*) fCT->GetBankRow("LCPB", j);
            TClonesArray &tLCPBbank = *LCPBStore;
            new(tLCPBbank[nLCPB_Store++]) TLCPBClass(fLCPB);
        }
        for(Int_t j = 0; j < (Int_t) fCT->GetNRows("STPB"); j++){
            fSTPB = (TSTPBClass*) fCT->GetBankRow("STPB", j);
            TClonesArray &tSTPBbank = *STPBStore;
            new(tSTPBbank[nSTPB_Store++]) TSTPBClass(fSTPB);
        }
        for(Int_t j = 0; j < (Int_t) fCT->GetNRows("TGPB"); j++){
            fTGPB = (TTGPBClass*) fCT->GetBankRow("TGPB", j);
            TClonesArray &tTGPBbank = *TGPBStore;
            new(tTGPBbank[nTGPB_Store++]) TTGPBClass(fTGPB);
        }
        for(Int_t j = 0; j < (Int_t) fCT->GetNRows("TAGR"); j++){
            fTAGR = (TTAGRClass*) fCT->GetBankRow("TAGR", j);
            TClonesArray &tTAGRbank = *TAGRStore;
            new(tTAGRbank[nTAGR_Store++]) TTAGRClass(fTAGR);
        }

        tree->Fill();
        Last_Key_Num++;  // Each Fill call increased Key by one.
    }

    f->cd();
    tree->Write();
    tree->Delete();
    f->Close();
}

void pCorrStepan::initClones(){
    fgEVNTStore = new TClonesArray("TEVNTClass",1);
    EVNTStore = fgEVNTStore;
    fgECPBStore = new TClonesArray("TECPBClass",1);
    ECPBStore  = fgECPBStore;
    fgSCPBStore = new TClonesArray("TSCPBClass",1);
    SCPBStore  = fgSCPBStore;
    fgDCPBStore = new TClonesArray("TDCPBClass",1);
    DCPBStore  = fgDCPBStore;
    fgCCPBStore = new TClonesArray("TCCPBClass",1);
    CCPBStore  = fgCCPBStore;
    fgLCPBStore = new TClonesArray("TLCPBClass",1);
    LCPBStore  = fgLCPBStore;
    fgSTPBStore = new TClonesArray("TSTPBClass",1);
    STPBStore  = fgSTPBStore;
    fgTGPBStore = new TClonesArray("TTGPBClass",1);
    TGPBStore  = fgTGPBStore;
    fgTAGRStore = new TClonesArray("TTAGRClass",1);
    TAGRStore  = fgTAGRStore;
}

void pCorrStepan::clearClones()
{
    if(fgEVNTStore) EVNTStore->Clear();
    if(fgECPBStore) ECPBStore->Clear();
    if(fgSCPBStore) SCPBStore->Clear();
    if(fgDCPBStore) DCPBStore->Clear();
    if(fgCCPBStore) CCPBStore->Clear();
    if(fgLCPBStore) LCPBStore->Clear();
    if(fgSTPBStore) STPBStore->Clear();
    if(fgTGPBStore) TGPBStore->Clear();
    if(fgTAGRStore) TAGRStore->Clear();
}

void pCorrStepan::initBranches(TTree *tree)
{
    int bufsize = __BUFSIZE__;
    int split   = 1;

    HEADERBranch = tree->Branch("HEADER","THEADERClass",&hHEADER,bufsize,split);
    tree->SetBranchStatus("HEADER",1);
    EVNTBranch = tree->Branch("EVNT",&fgEVNTStore,bufsize,split);
    ECPBBranch = tree->Branch("ECPB",&fgECPBStore,bufsize,split);
    SCPBBranch = tree->Branch("SCPB",&fgSCPBStore,bufsize,split);
    DCPBBranch = tree->Branch("DCPB",&fgDCPBStore,bufsize,split);
    CCPBBranch = tree->Branch("CCPB",&fgCCPBStore,bufsize,split);
    LCPBBranch = tree->Branch("LCPB",&fgLCPBStore,bufsize,split);
    STPBBranch = tree->Branch("STPB",&fgSTPBStore,bufsize,split);
    TGPBBranch = tree->Branch("TGPB",&fgTGPBStore,bufsize,split);
    TAGRBranch = tree->Branch("TAGR",&fgTAGRStore,bufsize,split);
    tree->SetBranchStatus("*",1);
}

Float_t pCorrStepan::calcF1(Int_t sect, Float_t phi)
{
    Float_t f1;
    f1 = f1Params[sect][0]+f1Params[sect][1]*phi+f1Params[sect][2]*phi*phi;
    return f1;
}

Float_t pCorrStepan::calcF2(Int_t sect, Float_t theta)
{
    Float_t f2;
    f2 = f2Params[sect][0]+(f2Params[sect][1]+f2Params[sect][2]*theta+f2Params[sect][3]*theta*theta)*TMath::Exp(theta);
    return f2;
}
