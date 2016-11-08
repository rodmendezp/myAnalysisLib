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
void drawNicePlot(string title, string xLabel, string yLabel, Bool_t drawLegend, TH1F *h);

#endif // MYROOTUTILS_H
