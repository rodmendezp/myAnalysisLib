#ifndef PCORRSTEPAN_H
#define PCORRSTEPAN_H

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
#include "TIdentificator.h"

#include "THEADERClass.h"
#include "TEVNTClass.h"
#include "TECPBClass.h"
#include "TSCPBClass.h"
#include "TDCPBClass.h"
#include "TCCPBClass.h"
#include "TLCPBClass.h"
#include "TSTPBClass.h"
#include "TTGPBClass.h"
#include "TTAGRClass.h"

class pCorrStepan
{
public:
    pCorrStepan(string txtfName);
    ~pCorrStepan();
    void loadF2eParams(string txtfName);
    void pCorrFile(string rootfName, Bool_t useF2e);
    Float_t calcF1(Int_t sect, Float_t phi);
    Float_t calcF2(Int_t sect, Float_t theta);
    Float_t calcF2e(Int_t sect, Float_t theta);
private:
    Int_t nSect;
    Float_t f1Params[6][3];
    Float_t f2Params[6][4];
    Float_t f2eParams[6][3];

    Int_t Last_Evt_Num;
    Int_t Last_Key_Num;

    TTree *tree;
    TClasTool *fCT;
    THEADERClass *fHEADER;
    TEVNTClass *fEVNT;
    TEVNTClass *nfEVNT;
    TECPBClass *fECPB;
    TSCPBClass *fSCPB;
    TDCPBClass *fDCPB;
    TCCPBClass *fCCPB;
    TLCPBClass *fLCPB;
    TSTPBClass *fSTPB;
    TTGPBClass *fTGPB;
    TTAGRClass *fTAGR;
    TVERTClass *fVERT;
    TMVRTClass *fMVRT;
    TTBERClass *fTBER;
    TIdentificator *fID;

    void initClones();
    void clearClones();
    void initBranches(TTree *tree);

    TClonesArray *EVNTStore;
    TClonesArray *fgEVNTStore;
    TClonesArray *ECPBStore;
    TClonesArray *fgECPBStore;
    TClonesArray *SCPBStore;
    TClonesArray *fgSCPBStore;
    TClonesArray *DCPBStore;
    TClonesArray *fgDCPBStore;
    TClonesArray *CCPBStore;
    TClonesArray *fgCCPBStore;
    TClonesArray *LCPBStore;
    TClonesArray *fgLCPBStore;
    TClonesArray *STPBStore;
    TClonesArray *fgSTPBStore;
    TClonesArray *TGPBStore;
    TClonesArray *fgTGPBStore;
    TClonesArray *TAGRStore;
    TClonesArray *fgTAGRStore;

    THEADERClass *hHEADER;
    TBranch *HEADERBranch;
    TBranch *EVNTBranch;
    Int_t nEVNT_Store;
    TBranch *ECPBBranch;
    Int_t nECPB_Store;
    TBranch *SCPBBranch;
    Int_t nSCPB_Store;
    TBranch *DCPBBranch;
    Int_t nDCPB_Store;
    TBranch *CCPBBranch;
    Int_t nCCPB_Store;
    TBranch *LCPBBranch;
    Int_t nLCPB_Store;
    TBranch *STPBBranch;
    Int_t nSTPB_Store;
    TBranch *TGPBBranch;
    Int_t nTGPB_Store;
    TBranch *TAGRBranch;
    Int_t nTAGR_Store;
};

#endif
