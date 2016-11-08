#ifndef MYROOTUTILS_H
#define MYROOTUTILS_H

#include <fstream>
#include <iostream>
#include <string>

#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TDirectory.h"

using namespace std;

Bool_t fexists(string fName);
void drawNicePlot(TH1 *h, TCanvas *ccNP, TString fName, TString title, TString xLabel, TString yLabel, Bool_t drawLegend);

#endif // MYROOTUTILS_H
