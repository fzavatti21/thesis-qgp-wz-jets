#ifndef MYLIB_H
#define MYLIB_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>
#include <random>
#include <tuple>
#include <cstdlib> 

#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TLatex.h"
#include "TProfile.h"
#include "TObject.h"
#include "TString.h"
#include "TGaxis.h"
#include "TError.h"
#include "TSystem.h"
#include "TColor.h"
// #include "TRandom3.h"

using namespace std;

#define MAXJETS 1000
#define MAX_XB  20000

// Functions
void mkDir(TString dir);
void CMS_latex();
void DrawLatex(Float_t x, Float_t y, Int_t color, const char* text, Float_t textSize = 0.06);
void canvasSetup(TCanvas* c);
void histSetup(TH1* histo, const TString& x_label, const TString& y_label, Color_t color, bool fill);
void RatioPlot(TCanvas* c, TH1* h1, TH1* h2, TString h1_leg, TString h2_leg, Double_t ratio_y_inf, Double_t ratio_y_sup, TString ratio_leg);
void RatioPlot_severalHists(TCanvas* c, vector<TH1*> hists, vector<TString> legends, vector <double> leg_coord, 
 Double_t ratio_y_inf, Double_t ratio_y_sup, TString ratio_leg, Double_t xMin, Double_t xMax, bool logy_index);
void acusticEnd();
void PlotSeveralHists(TCanvas* c, vector<TH1*> hists, vector<TString> legends);
void test_IntegratedOverPt   (TString Djet, TCanvas* canvas, TH2D* hNum2D[], TH2D* hDenom2D, Int_t nIter, TString algo, TString outDir);
void test_IntegratedOverTheta(TString Djet, TCanvas* canvas, TH2D* hNum2D[], TH2D* hDenom2D, Int_t nIter, TString algo, TString outDir);
void test_InBinsOfTheta      (TString Djet, TCanvas* canvas, TH2D* hNum2D[], TH2D* hDenom2D, Int_t nIter, TString algo, TString outDir);
void test_InBinsOfPt         (TString Djet, TCanvas* canvas, TH2D* hNum2D[], TH2D* hDenom2D, Int_t nIter, TString algo, TString outDir);
Int_t FindSelectedIterFromHist(TString Djet, TH1* pValue_h);
void normalizeTH2Column(TH2* h);
void normalizeTH2Row(TH2* h);

// Utilities    
    // Colors for cout
TString cReset     = "\033[0m";
TString cRed       = "\033[31m";
TString cGreen     = "\033[32m";
TString cYellow    = "\033[33m";
TString cBlue      = "\033[34m";
TString cMagenta   = "\033[35m";
TString cCyan      = "\033[36m";

    // My own colors for Plotting
Int_t rgbNorm = 255.;
Int_t kMyCyan = TColor::GetFreeColorIndex();
auto myCyan = new TColor(kMyCyan, 10./rgbNorm, 119./rgbNorm, 93./rgbNorm);
Int_t kMyRed  = TColor::GetFreeColorIndex();
auto myRed  = new TColor(kMyRed, 153./rgbNorm, 0./rgbNorm, 0./rgbNorm);

vector<TString> algos{"SD", "latekt"};
TString algo; 

Double_t theta_bins[] = {-1., 0., 1.6, 2., 2.3, 2.55, 3.};
// Double_t theta_bins[] = {-1., 0., 1.6, 2., 2.3, 2.55, 2.8, 3.};

// Double_t jetpt_bins[] = {60., 80., 100., 120., 180.};
Double_t jetpt_bins[] = {80., 100., 120., 180.};
// Double_t DCAS_bins[]  = {0., 1., 2., 4., 10., 18.};
Double_t DCAS_bins[]  = {0., 1., 2., 5., 8.};

Int_t    theta_nbins  = sizeof(theta_bins) / sizeof(theta_bins[0]) - 1;
Int_t    jetpt_nbins  = sizeof(jetpt_bins) / sizeof(jetpt_bins[0]) - 1;
Int_t    DCAS_nbins   = sizeof(DCAS_bins)  / sizeof(DCAS_bins[0])  - 1;


Double_t avoid_errors[8] = {0.0000001,0.0000001,0.0000001,0.0000001,0.0000001,0.0000001,0.0000001,0.0000001};

#endif 
