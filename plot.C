#include "myLib.h"
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>
#include <vector>
#include <fstream>
#include <memory>

using namespace std;

void CMS_latex(){
    TLatex* cmsColl = new TLatex();
    cmsColl->SetTextSize(0.04);
    cmsColl->SetTextFont(61);
    cmsColl->SetLineWidth(2);
    cmsColl->SetTextColor(kBlack);
    cmsColl->SetNDC();
    cmsColl->SetTextAlign(11);

    TLatex* cmsPrel = new TLatex();
    cmsPrel->SetTextSize(0.04);
    cmsPrel->SetTextFont(52);
    cmsPrel->SetLineWidth(2);
    cmsPrel->SetTextColor(kBlack);
    cmsPrel->SetNDC();
    cmsPrel->SetTextAlign(11);

    TLatex* lumiText = new TLatex();
    lumiText->SetTextSize(0.03);
    lumiText->SetTextFont(42);
    lumiText->SetLineWidth(2);
    lumiText->SetTextColor(kBlack);
    lumiText->SetNDC();
    lumiText->SetTextAlign(31);

    TString extraText = "Sim.";
    TString lumi = "Pythia Tune 4C, #sqrt{s} = 5.02 TeV, #hat{p}_{T}^{min}=150 GeV";

    cmsColl->DrawLatex(0.15-0.005, 0.942, "CMS");
    cmsPrel->DrawLatex(0.22-0.005, 0.942, extraText);
    lumiText->DrawLatex(0.97-0.005, 0.942, lumi);
}

//  struct Sample
struct Sample {
    TString name;
    TString fileList;
    Color_t color;

    // ---- spettri pT ----
    TH1D* h_ptWjet           = nullptr;
    TH1D* h_ptZjet           = nullptr;
    TH1D* h_ptWjet_pt150_250 = nullptr;
    TH1D* h_ptWjet_pt250_350 = nullptr;
    TH1D* h_ptWjet_pt350_500 = nullptr;
    TH1D* h_ptZjet_pt150_250 = nullptr;
    TH1D* h_ptZjet_pt250_350 = nullptr;
    TH1D* h_ptZjet_pt350_500 = nullptr;

    // ---- N2 ----
    TH1D* h_n2Fat     = nullptr;
    TH1D* h_n2Fat_pure     = nullptr;

    TH1D* h_n2Fat_WMatched = nullptr;
    TH1D* h_n2Fat_ZMatched = nullptr;
    TH1D* h_n2Fat_WZMatched = nullptr;
    TH1D* h_n2Fat_bkg = nullptr;

    TH1D* h_n2Fat_WMatched_pure = nullptr;
    TH1D* h_n2Fat_ZMatched_pure = nullptr;
    TH1D* h_n2Fat_WZMatched_pure = nullptr;
    TH1D* h_n2Fat_bkg_pure = nullptr;

   

    TH2D* h2D_n2_n2pure_WZ    = nullptr;



    // ----  n2Fat vs massa ----    
    TH2D* h2D_mFat_n2_W = nullptr;
    TH2D* h2D_mFat_n2_Z = nullptr;
    TH2D* h2D_mFat_n2_bkg = nullptr;

    //n2Fat_pure vs massa
    TH2D* h2D_mFat_n2pure_W = nullptr;
    TH2D* h2D_mFat_n2pure_Z = nullptr;
    TH2D* h2D_mFat_n2pure_bkg = nullptr;


    // ---- masse n2 ----
    TH1D* h_mFat_All_all             = nullptr;
    TH1D* h_mFat_All_n2           = nullptr;
    TH1D* h_mFat_WMatched_all        = nullptr;
    TH1D* h_mFat_WMatched_n2      = nullptr;
    TH1D* h_mFat_ZMatched_all        = nullptr;
    TH1D* h_mFat_ZMatched_n2      = nullptr;
    TH1D* h_mFat_All_n2_pt150_250 = nullptr;
    TH1D* h_mFat_All_n2_pt250_350 = nullptr;
    TH1D* h_mFat_All_n2_pt350_500 = nullptr;

    


//massa pure
    TH1D* h_mFat_All_n2pure           = nullptr;
    TH1D* h_mFat_All_n2pure_pt150_250 = nullptr;
    TH1D* h_mFat_All_n2pure_pt250_350 = nullptr;
    TH1D* h_mFat_All_n2pure_pt350_500 = nullptr;

    TH1D* h_mFat_All_n2pure_mix = nullptr;  // mFat embedded con n2pure < 0.18
    TH1D* h_mFat_All_n2pure_mix_pt150_250 = nullptr;
    TH1D* h_mFat_All_n2pure_mix_pt250_350 = nullptr;
    TH1D* h_mFat_All_n2pure_mix_pt350_500 = nullptr;



  

    
    

   



    // ---- tempi W: totali ----
    TH1D* h_WdecayTime_all  = nullptr;
    TH1D* h_WcoherTime_all  = nullptr;
    TH1D* h_WtotalTime_all  = nullptr;

    // ---- tempi W: per bin di pT (bin1:150-300, bin2:300-400, bin3:>400) ----
    TH1D* h_WtotalTime_bin1 = nullptr;
    TH1D* h_WtotalTime_bin2 = nullptr;
    TH1D* h_WtotalTime_bin3 = nullptr;

    // ---- tempi Z: totali ----
    TH1D* h_ZdecayTime_all  = nullptr;
    TH1D* h_ZcoherTime_all  = nullptr;
    TH1D* h_ZtotalTime_all  = nullptr;

    // ---- tempi Z: per bin di pT ----
    TH1D* h_ZtotalTime_bin1 = nullptr;
    TH1D* h_ZtotalTime_bin2 = nullptr;
    TH1D* h_ZtotalTime_bin3 = nullptr;

    // ---- groomed mass con taglio n2Fat < 0.18 ----
TH1D* h_mFat_groomed_n2           = nullptr;
TH1D* h_mFat_groomed_n2_pt150_250 = nullptr;
TH1D* h_mFat_groomed_n2_pt250_350 = nullptr;
TH1D* h_mFat_groomed_n2_pt350_500 = nullptr;
TH1D* h_mFat_WMatched_groomed_n2  = nullptr;
TH1D* h_mFat_ZMatched_groomed_n2  = nullptr;

// ---- groomed mass con taglio n2Fat_pure < 0.18 ----
TH1D* h_mFat_groomed_n2pure           = nullptr;
TH1D* h_mFat_groomed_n2pure_pt150_250 = nullptr;
TH1D* h_mFat_groomed_n2pure_pt250_350 = nullptr;
TH1D* h_mFat_groomed_n2pure_pt350_500 = nullptr;

// ---- groomed pure mass con taglio n2Fat_pure < 0.18 ----
TH1D* h_mFatpure_groomed_n2pure           = nullptr;
TH1D* h_mFatpure_groomed_n2pure_pt150_250 = nullptr;
TH1D* h_mFatpure_groomed_n2pure_pt250_350 = nullptr;
TH1D* h_mFatpure_groomed_n2pure_pt350_500 = nullptr;




TH1D* h_mFat_WMatched_groomed  = nullptr;
TH1D* h_mFat_ZMatched_groomed  = nullptr;

    // ---- accumulatori ----
    double sumW_gen = 0, sumZ_gen = 0;
    double sumW_filled = 0, sumZ_filled = 0;
    int nTotal = 0, nAll_n2cut = 0;
    int nAll_n2cut_150_250 = 0, nAll_n2cut_250_350 = 0, nAll_n2cut_350_500 = 0;
    int nWall = 0, nZall = 0;
    int nWMatched_all = 0, nZMatched_all = 0;
    int nWMatched_all_150_250 = 0, nZMatched_all_150_250 = 0;
    int nWMatched_all_250_350 = 0, nZMatched_all_250_350 = 0;
    int nWMatched_all_350_500 = 0, nZMatched_all_350_500 = 0;
    int nWMatched_n2cut = 0, nZMatched_n2cut = 0;

    Sample(const TString& n, const TString& fl, Color_t c)
        : name(n), fileList(fl), color(c)
    {
        TString pf = n + "_";

        h_mFat_All_all              = new TH1D(pf+"mFat_All_all",              "", 150, 0, 150);
        h_mFat_All_n2            = new TH1D(pf+"mFat_All_n2",            "", 150, 0, 150);
        h_mFat_WMatched_all         = new TH1D(pf+"mFat_WMatched_all",         "", 150, 0, 150);
        h_mFat_WMatched_n2       = new TH1D(pf+"mFat_WMatched_n2",       "", 150, 0, 150);
        h_mFat_ZMatched_all         = new TH1D(pf+"mFat_ZMatched_all",         "", 150, 0, 150);
        h_mFat_ZMatched_n2       = new TH1D(pf+"mFat_ZMatched_n2",       "", 150, 0, 150);
        h_mFat_All_n2_pt150_250  = new TH1D(pf+"mFat_All_n2_pt150_250",  "", 150, 0, 150);
        h_mFat_All_n2_pt250_350  = new TH1D(pf+"mFat_All_n2_pt250_350",  "", 150, 0, 150);
        h_mFat_All_n2_pt350_500  = new TH1D(pf+"mFat_All_n2_pt350_500",  "", 150, 0, 150);

        h_mFat_All_n2pure_mix = new TH1D(pf+"mFat_All_n2pure_mix", "", 150, 0, 150);
        h_mFat_All_n2pure_mix_pt150_250  = new TH1D(pf+"mFat_All_n2pure_mix_pt150_250",  "", 150, 0, 150);
        h_mFat_All_n2pure_mix_pt250_350  = new TH1D(pf+"mFat_All_n2pure_mix_pt250_350",  "", 150, 0, 150);
        h_mFat_All_n2pure_mix_pt350_500  = new TH1D(pf+"mFat_All_n2pure_mix_pt350_500",  "", 150, 0, 150);


       


        h_ptWjet           = new TH1D(pf+"ptWjet",           "", 160, 0, 800);
        h_ptZjet           = new TH1D(pf+"ptZjet",           "", 160, 0, 800);
        h_ptWjet_pt150_250 = new TH1D(pf+"ptWjet_pt150_250", "", 160, 0, 800);
        h_ptWjet_pt250_350 = new TH1D(pf+"ptWjet_pt250_350", "", 160, 0, 800);
        h_ptWjet_pt350_500 = new TH1D(pf+"ptWjet_pt350_500", "", 160, 0, 800);
        h_ptZjet_pt150_250 = new TH1D(pf+"ptZjet_pt150_250", "", 160, 0, 800);
        h_ptZjet_pt250_350 = new TH1D(pf+"ptZjet_pt250_350", "", 160, 0, 800);
        h_ptZjet_pt350_500 = new TH1D(pf+"ptZjet_pt350_500", "", 160, 0, 800);

        h_n2Fat     = new TH1D(pf+"n2Fat",       "", 100, 0, 1);
        h_n2Fat_pure     = new TH1D(pf+"n2Fat probe",       "", 100, 0, 1);
        h_n2Fat_WMatched  = new TH1D(pf+"n2Fat_WMatched",  "", 100, 0, 1);
        h_n2Fat_ZMatched  = new TH1D(pf+"n2Fat_ZMatched",  "", 100, 0, 1);
        h_n2Fat_WZMatched = new TH1D(pf+"n2Fat_WZMatched", "", 100, 0, 1);
        h_n2Fat_bkg  = new TH1D(pf+"n2Fat_bkg", "", 100, 0, 1);
        h_n2Fat_WMatched_pure  = new TH1D(pf+"n2Fat_WMatched probe",  "", 100, 0, 1);
        h_n2Fat_ZMatched_pure  = new TH1D(pf+"n2Fat_ZMatched probe",  "", 100, 0, 1);
        h_n2Fat_WZMatched_pure = new TH1D(pf+"n2Fat_WZMatched probe", "", 100, 0, 1);
        h_n2Fat_bkg_pure  = new TH1D(pf+"n2Fat_bkg_pure", "", 100, 0, 1);

       
        h2D_n2_n2pure_WZ = new TH2D(pf+"h2D_n2_n2pure_WZ",     "", 100, 0, 0.6, 100, 0, 0.6);

        h2D_n2_n2pure_WZ->GetXaxis()->SetTitle("N_{2} embedded");
        h2D_n2_n2pure_WZ->GetYaxis()->SetTitle("N_{2} pure");
       

         
         h_mFat_All_n2pure = new TH1D(pf+"mFat_All_n2pure", "", 150, 0, 150);
         h_mFat_All_n2pure_pt150_250 = new TH1D(pf+"mFat_All_n2pure_pt150_250", "", 150, 0, 150);
         h_mFat_All_n2pure_pt250_350 = new TH1D(pf+"mFat_All_n2pure_pt250_350", "", 150, 0, 150);
         h_mFat_All_n2pure_pt350_500 = new TH1D(pf+"mFat_All_n2pure_pt350_500", "", 150, 0, 150);

        
        

        h2D_mFat_n2pure_W      = new TH2D(pf+"h2D_mFat_n2pure_W",    "", 150, 0, 150, 100, 0, 0.6);
        h2D_mFat_n2pure_Z      = new TH2D(pf+"h2D_mFat_n2pure_Z",    "", 150, 0, 150, 100, 0, 0.6);
        h2D_mFat_n2pure_bkg    = new TH2D(pf+"h2D_mFat_n2pure_bkg",  "", 150, 0, 150, 100, 0, 0.6);


        h2D_mFat_n2_W = new TH2D(pf+"h2D_mFat_n2_W", "", 150, 0, 150, 100, 0, 0.6);
        h2D_mFat_n2_Z = new TH2D(pf+"h2D_mFat_n2_Z", "", 150, 0, 150, 100, 0, 0.6);
        h2D_mFat_n2_bkg = new TH2D(pf+"h2D_mFat_n2_bkg", "", 150, 0, 150, 100, 0, 0.6);

        h2D_mFat_n2_W->GetXaxis()->SetTitle("Inv. Mass [GeV]");
        h2D_mFat_n2pure_W->GetYaxis()->SetTitle("N_{2}");

        h2D_mFat_n2_Z->GetXaxis()->SetTitle("Inv. Mass [GeV]");
        h2D_mFat_n2pure_Z->GetYaxis()->SetTitle("N_{2}");

        h2D_mFat_n2_bkg->GetXaxis()->SetTitle("Inv. Mass [GeV]");
        h2D_mFat_n2pure_bkg->GetYaxis()->SetTitle("N_{2}");



        h_WdecayTime_all  = new TH1D(pf+"WdecayTime_all",  "", 100, 0, 2);
        h_WcoherTime_all  = new TH1D(pf+"WcoherTime_all",  "", 100, 0, 2);
        h_WtotalTime_all  = new TH1D(pf+"WtotalTime_all",  "", 100, 0, 5);
        h_WtotalTime_bin1 = new TH1D(pf+"WtotalTime_bin1", "", 100, 0, 5);
        h_WtotalTime_bin2 = new TH1D(pf+"WtotalTime_bin2", "", 100, 0, 5);
        h_WtotalTime_bin3 = new TH1D(pf+"WtotalTime_bin3", "", 100, 0, 5);

        h_ZdecayTime_all  = new TH1D(pf+"ZdecayTime_all",  "", 100, 0, 2);
        h_ZcoherTime_all  = new TH1D(pf+"ZcoherTime_all",  "", 100, 0, 2);
        h_ZtotalTime_all  = new TH1D(pf+"ZtotalTime_all",  "", 100, 0, 5);
        h_ZtotalTime_bin1 = new TH1D(pf+"ZtotalTime_bin1", "", 100, 0, 5);
        h_ZtotalTime_bin2 = new TH1D(pf+"ZtotalTime_bin2", "", 100, 0, 5);
        h_ZtotalTime_bin3 = new TH1D(pf+"ZtotalTime_bin3", "", 100, 0, 5);

        h_mFat_groomed_n2           = new TH1D(pf+"mFat_groomed_n2",           "", 150, 60, 150);
        h_mFat_groomed_n2_pt150_250 = new TH1D(pf+"mFat_groomed_n2_pt150_250", "", 150, 60, 150);
        h_mFat_groomed_n2_pt250_350 = new TH1D(pf+"mFat_groomed_n2_pt250_350", "", 150, 60, 150);
        h_mFat_groomed_n2_pt350_500 = new TH1D(pf+"mFat_groomed_n2_pt350_500", "", 150, 60, 150);
        h_mFat_WMatched_groomed_n2  = new TH1D(pf+"mFat_WMatched_groomed_n2",  "", 150, 60, 150);
        h_mFat_ZMatched_groomed_n2  = new TH1D(pf+"mFat_ZMatched_groomed_n2",  "", 150, 60, 150);


        h_mFat_groomed_n2pure           = new TH1D(pf+"mFat_groomed_n2pure",           "", 150, 60, 150);
        h_mFat_groomed_n2pure_pt150_250 = new TH1D(pf+"mFat_groomed_n2pure_pt150_250", "", 150, 60, 150);
        h_mFat_groomed_n2pure_pt250_350 = new TH1D(pf+"mFat_groomed_n2pure_pt250_350", "", 150, 60, 150);
        h_mFat_groomed_n2pure_pt350_500 = new TH1D(pf+"mFat_groomed_n2pure_pt350_500", "", 150, 60, 150);

        h_mFat_WMatched_groomed  = new TH1D(pf+"mFat_WMatched_groomed",  "", 150, 60, 150);
        h_mFat_ZMatched_groomed  = new TH1D(pf+"mFat_ZMatched_groomed",  "", 150, 60, 150);

        h_mFatpure_groomed_n2pure           = new TH1D(pf+"mFatpure_groomed_n2pure",           "", 150, 60, 150);
        h_mFatpure_groomed_n2pure_pt150_250 = new TH1D(pf+"mFatpure_groomed_n2pure_pt150_250", "", 150, 60, 150);
        h_mFatpure_groomed_n2pure_pt250_350 = new TH1D(pf+"mFatpure_groomed_n2pure_pt250_350", "", 150, 60, 150);
        h_mFatpure_groomed_n2pure_pt350_500 = new TH1D(pf+"mFatpure_groomed_n2pure_pt350_500", "", 150, 60, 150);



        for (TH1D* h : {
            h_mFat_All_all, h_mFat_All_n2,
            h_mFat_WMatched_all, h_mFat_WMatched_n2,
            h_mFat_ZMatched_all, h_mFat_ZMatched_n2,
            h_mFat_All_n2_pt150_250, h_mFat_All_n2_pt250_350, h_mFat_All_n2_pt350_500,
           
            h_ptWjet, h_ptZjet,
            h_ptWjet_pt150_250, h_ptWjet_pt250_350, h_ptWjet_pt350_500,
            h_ptZjet_pt150_250, h_ptZjet_pt250_350, h_ptZjet_pt350_500,
            h_n2Fat,h_n2Fat_pure,h_n2Fat_WMatched, h_n2Fat_ZMatched, h_n2Fat_WZMatched,
            h_n2Fat_WMatched_pure, h_n2Fat_ZMatched_pure, h_n2Fat_WZMatched_pure,
            h_mFat_All_n2pure,h_mFat_All_n2pure_pt150_250,h_mFat_All_n2pure_pt250_350,h_mFat_All_n2pure_pt350_500,
            h_n2Fat_bkg,h_n2Fat_bkg_pure,
            h_WdecayTime_all, h_WcoherTime_all, h_WtotalTime_all,
            h_WtotalTime_bin1, h_WtotalTime_bin2, h_WtotalTime_bin3,
            h_ZdecayTime_all, h_ZcoherTime_all, h_ZtotalTime_all,
            h_ZtotalTime_bin1, h_ZtotalTime_bin2, h_ZtotalTime_bin3,
            h_mFat_All_n2pure_mix,h_mFat_All_n2pure_mix_pt150_250,h_mFat_All_n2pure_mix_pt250_350,h_mFat_All_n2pure_mix_pt350_500,

            h_mFat_groomed_n2, h_mFat_groomed_n2_pt150_250,
            h_mFat_groomed_n2_pt250_350, h_mFat_groomed_n2_pt350_500,
            h_mFat_WMatched_groomed_n2, h_mFat_ZMatched_groomed_n2,
            h_mFat_WMatched_groomed, h_mFat_ZMatched_groomed,h_mFat_groomed_n2pure,h_mFat_groomed_n2pure_pt150_250,
            h_mFat_groomed_n2pure_pt250_350, h_mFat_groomed_n2pure_pt350_500,h_mFatpure_groomed_n2pure,h_mFatpure_groomed_n2pure_pt150_250,
            h_mFatpure_groomed_n2pure_pt250_350, h_mFatpure_groomed_n2pure_pt350_500
        }) h->Sumw2();

        h2D_mFat_n2_W->Sumw2();
        h2D_mFat_n2_Z->Sumw2();
        h2D_mFat_n2_bkg->Sumw2();
        for (TH2D* h : { h2D_mFat_n2pure_W,  h2D_mFat_n2pure_Z,  h2D_mFat_n2pure_bkg})
        h->Sumw2();
         h2D_n2_n2pure_WZ->Sumw2();
        }
};

 
//  fillSample
 
void fillSample(Sample& s) {
    ifstream list(s.fileList.Data());
    if (!list.is_open()) { cout << "ERRORE: " << s.fileList << endl; return; }
    TString fname;
    while (list >> fname) {
        unique_ptr<TFile> f(TFile::Open(fname));
        if (!f || f->IsZombie()) { cout << "Errore apertura: " << fname << endl; continue; }
        TTree* tEvent = (TTree*)f->Get("EventTree");
        if (!tEvent) { cout << "Errore: EventTree non trovato in " << fname << endl; continue; }

        Float_t mFat, n2Fat,n2Fat_pure, evweight, pFat;
        Float_t ev_isW, ev_isZ, ev_WMatchedFatJet, ev_ZMatchedFatJet;
        Float_t ptW, etaW, phiW, ptZ, etaZ, phiZ;
        Float_t hofatjet;
        Float_t WdecayTime, WcoherenceTime, WtotalTime;
        Float_t ZdecayTime, ZcoherenceTime, ZtotalTime;
        Float_t mFat_pure,ev_WMatchedprobe, ev_ZMatchedprobe, hoprobe;
        Float_t mFat_groomed,mFatpure_groomed;


        tEvent->SetBranchAddress("mFat",            &mFat);
        tEvent->SetBranchAddress("n2Fat",            &n2Fat);
        tEvent->SetBranchAddress("n2Fat_pure",            &n2Fat_pure);
        tEvent->SetBranchAddress("evweight",         &evweight);
        tEvent->SetBranchAddress("isWEvent",         &ev_isW);
        tEvent->SetBranchAddress("ptW",              &ptW);
        tEvent->SetBranchAddress("etaW",             &etaW);
        tEvent->SetBranchAddress("phiW",             &phiW);
        tEvent->SetBranchAddress("pFat",             &pFat);
        tEvent->SetBranchAddress("WMatchedFatJet",   &ev_WMatchedFatJet);
        tEvent->SetBranchAddress("isZEvent",         &ev_isZ);
        tEvent->SetBranchAddress("ZMatchedFatJet",   &ev_ZMatchedFatJet);
        tEvent->SetBranchAddress("ptZ",              &ptZ);
        tEvent->SetBranchAddress("etaZ",             &etaZ);
        tEvent->SetBranchAddress("phiZ",             &phiZ);
        tEvent->SetBranchAddress("hofatjet",         &hofatjet);
        tEvent->SetBranchAddress("WdecayTime",       &WdecayTime);
        tEvent->SetBranchAddress("WcoherenceTime",   &WcoherenceTime);
        tEvent->SetBranchAddress("WtotalTime",       &WtotalTime);
        tEvent->SetBranchAddress("ZdecayTime",       &ZdecayTime);
        tEvent->SetBranchAddress("ZcoherenceTime",   &ZcoherenceTime);
        tEvent->SetBranchAddress("ZtotalTime",       &ZtotalTime);
        tEvent->SetBranchAddress("mFat_pure", &mFat_pure);
        tEvent->SetBranchAddress("hoprobe", &hoprobe);
        tEvent->SetBranchAddress("WMatchedprobe", &ev_WMatchedprobe);
        tEvent->SetBranchAddress("ZMatchedprobe", &ev_ZMatchedprobe);
        tEvent->SetBranchAddress("mFat_groomed", &mFat_groomed);
        tEvent->SetBranchAddress("mFatpure_groomed", &mFatpure_groomed);


        for (Long64_t ie = 0; ie < tEvent->GetEntries(); ++ie) {
            tEvent->GetEntry(ie);

            if (hofatjet == 1) s.nTotal++;
            if (ev_isW == 1)   s.sumW_gen += evweight;
            if (ev_isZ == 1)   s.sumZ_gen += evweight;

              if (mFat > 0 && hofatjet == 1) s.h_mFat_All_all->Fill(mFat, evweight);

              if (mFat > 0 && n2Fat_pure >= 0 && n2Fat_pure < 0.18 && hofatjet == 1){  
                  s.h_mFat_All_n2pure_mix->Fill(mFat, evweight);
                   if      (pFat >= 150 && pFat < 250) s.h_mFat_All_n2pure_mix_pt150_250->Fill(mFat, evweight);
                   else if (pFat >= 250 && pFat < 350) s.h_mFat_All_n2pure_mix_pt250_350->Fill(mFat, evweight);
                    else if (pFat >= 350 && pFat < 500) s.h_mFat_All_n2pure_mix_pt350_500->Fill(mFat, evweight);
              }
            //massa con n2Fat<0.18
            if (mFat > 0 && n2Fat >= 0 && n2Fat < 0.18 && hofatjet == 1) {
                s.h_mFat_All_n2->Fill(mFat, evweight);
                s.nAll_n2cut++;
                if      (pFat >= 150 && pFat < 250) { s.h_mFat_All_n2_pt150_250->Fill(mFat, evweight); s.nAll_n2cut_150_250++; }
                else if (pFat >= 250 && pFat < 350) { s.h_mFat_All_n2_pt250_350->Fill(mFat, evweight); s.nAll_n2cut_250_350++; }
                else if (pFat >= 350 && pFat < 500) { s.h_mFat_All_n2_pt350_500->Fill(mFat, evweight); s.nAll_n2cut_350_500++; }
            }
            

             // massa con n2Fat_pure < 0.18
            if(mFat_pure > 0 && n2Fat_pure >= 0 && n2Fat_pure < 0.18 && hoprobe == 1){
             s.h_mFat_All_n2pure->Fill(mFat_pure, evweight);
            if      (pFat >= 150 && pFat < 250) s.h_mFat_All_n2pure_pt150_250->Fill(mFat_pure, evweight);
            else if (pFat >= 250 && pFat < 350) s.h_mFat_All_n2pure_pt250_350->Fill(mFat_pure, evweight);
            else if (pFat >= 350 && pFat < 500) s.h_mFat_All_n2pure_pt350_500->Fill(mFat_pure, evweight);
            }

            
            

            //N2 1 d all
            if (n2Fat >= 0 && hofatjet == 1) {
                s.h_n2Fat->Fill(n2Fat, evweight);
            }
           
          
            if (ev_isW == 1 && hofatjet == 1) s.nWall++;
            if (ev_isZ == 1 && hofatjet == 1) s.nZall++;

            if (ev_isW == 1 && ev_WMatchedFatJet == 1 && mFat > 0 && hofatjet == 1) {
                s.h_mFat_WMatched_all->Fill(mFat, evweight);
                s.nWMatched_all++;
                if (n2Fat >= 0) {
                s.h_n2Fat_WMatched->Fill(n2Fat, evweight);
                s.h_n2Fat_WZMatched->Fill(n2Fat, evweight);
                s.h2D_mFat_n2_W->Fill(mFat, n2Fat, evweight);
                 if (n2Fat_pure >= 0) s.h2D_n2_n2pure_WZ->Fill(n2Fat,n2Fat_pure, evweight);
                }
                
               s.h_ptWjet->Fill(pFat, evweight);
                s.sumW_filled += evweight;
               
                if      (pFat >= 150 && pFat < 250) { s.h_ptWjet_pt150_250->Fill(pFat, evweight); s.nWMatched_all_150_250++; }
                else if (pFat >= 250 && pFat < 350) { s.h_ptWjet_pt250_350->Fill(pFat, evweight); s.nWMatched_all_250_350++; }
                else if (pFat >= 350 && pFat < 500) { s.h_ptWjet_pt350_500->Fill(pFat, evweight); s.nWMatched_all_350_500++; }
                if (n2Fat >= 0 && n2Fat <= 0.18)    { s.h_mFat_WMatched_n2->Fill(mFat, evweight); s.nWMatched_n2cut++; }
            }

            if (ev_isW == 0 && ev_isZ == 0  && ev_WMatchedFatJet == 0 && ev_ZMatchedFatJet == 0 && n2Fat >= 0 && mFat > 0 && hofatjet == 1) {
                s.h_n2Fat_bkg->Fill(n2Fat, evweight);
            }
            
            if (ev_isW == 0 && ev_isZ == 0  && ev_WMatchedprobe == 0 && ev_ZMatchedprobe == 0 && n2Fat_pure >= 0 && mFat_pure > 0 && hoprobe == 1) {
            s.h_n2Fat_bkg_pure->Fill(n2Fat_pure, evweight);
            }




            if (ev_isZ == 1 && ev_ZMatchedFatJet == 1 && mFat > 0 && hofatjet == 1) {
                s.h_mFat_ZMatched_all->Fill(mFat, evweight);
                s.nZMatched_all++;
                if (n2Fat >= 0) {
        s.h_n2Fat_ZMatched->Fill(n2Fat, evweight);
        s.h_n2Fat_WZMatched->Fill(n2Fat, evweight);
        s.h2D_mFat_n2_Z->Fill(mFat, n2Fat, evweight);
        if (n2Fat_pure >= 0) s.h2D_n2_n2pure_WZ->Fill(n2Fat,n2Fat_pure, evweight);
    }
     

    // W matched su entrambi
    if (ev_isW == 1 && ev_WMatchedFatJet == 1 && ev_WMatchedprobe == 1 
    && hofatjet == 1 && hoprobe == 1 
    && n2Fat >= 0 && n2Fat_pure >= 0) {
    s.h2D_n2_n2pure_WZ->Fill(n2Fat, n2Fat_pure, evweight);
    }
    if (ev_isZ == 1 && ev_ZMatchedFatJet == 1 && ev_ZMatchedprobe == 1 
    && hofatjet == 1 && hoprobe == 1 
    && n2Fat >= 0 && n2Fat_pure >= 0) {
    s.h2D_n2_n2pure_WZ->Fill(n2Fat, n2Fat_pure, evweight);
    }

    
               s.h_ptZjet->Fill(pFat, evweight);
                s.sumZ_filled += evweight;
                if      (pFat >= 150 && pFat < 250) { s.h_ptZjet_pt150_250->Fill(pFat, evweight); s.nZMatched_all_150_250++; }
                else if (pFat >= 250 && pFat < 350) { s.h_ptZjet_pt250_350->Fill(pFat, evweight); s.nZMatched_all_250_350++; }
                else if (pFat >= 350 && pFat < 500) { s.h_ptZjet_pt350_500->Fill(pFat, evweight); s.nZMatched_all_350_500++; }
                if (n2Fat >= 0 && n2Fat <= 0.18)    { s.h_mFat_ZMatched_n2->Fill(mFat, evweight); s.nZMatched_n2cut++; }
            }
            
            if ( ev_ZMatchedFatJet == 0 && mFat > 0 && hofatjet == 1 && ev_WMatchedFatJet == 0 && n2Fat>=0){
                                s.h2D_mFat_n2_bkg->Fill(mFat, n2Fat, evweight);
            }

           
            
            if (ev_isW == 1 && ev_WMatchedprobe == 1 && mFat_pure > 0 && hoprobe == 1) {
            if (n2Fat_pure >= 0) {
                s.h_n2Fat_WMatched_pure->Fill(n2Fat_pure, evweight);
                s.h_n2Fat_WZMatched_pure->Fill(n2Fat_pure, evweight);
             s.h2D_mFat_n2pure_W->Fill(mFat_pure, n2Fat_pure, evweight);
             }
            }

            if (ev_isZ == 1 && ev_ZMatchedprobe == 1 && mFat_pure > 0 && hoprobe == 1) {
            if (n2Fat_pure >= 0) {
                s.h_n2Fat_ZMatched_pure->Fill(n2Fat_pure, evweight);
                s.h_n2Fat_WZMatched_pure->Fill(n2Fat_pure, evweight);
             s.h2D_mFat_n2pure_Z->Fill(mFat_pure, n2Fat_pure, evweight);
             }

            }
            if ( ev_ZMatchedprobe == 0 && mFat_pure > 0 && hoprobe == 1 && ev_WMatchedprobe == 0 && n2Fat_pure>=0){
                                s.h2D_mFat_n2pure_bkg->Fill(mFat, n2Fat, evweight);
            }

            


            // tempi W
            if (ev_isW == 1 && hofatjet == 1) {
                s.h_WdecayTime_all->Fill(WdecayTime,     evweight);
                s.h_WcoherTime_all->Fill(WcoherenceTime, evweight);
                s.h_WtotalTime_all->Fill(WtotalTime,     evweight);
                int bin = (ptW >= 150 && ptW < 300) ? 1 :
                          (ptW >= 300 && ptW < 400) ? 2 :
                          (ptW >= 400)              ? 3 : 0;
                if      (bin==1) s.h_WtotalTime_bin1->Fill(WtotalTime, evweight);
                else if (bin==2) s.h_WtotalTime_bin2->Fill(WtotalTime, evweight);
                else if (bin==3) s.h_WtotalTime_bin3->Fill(WtotalTime, evweight);
            }

            // tempi Z
            if (ev_isZ == 1 && hofatjet == 1) {
                s.h_ZdecayTime_all->Fill(ZdecayTime,     evweight);
                s.h_ZcoherTime_all->Fill(ZcoherenceTime, evweight);
                s.h_ZtotalTime_all->Fill(ZtotalTime,     evweight);
                int bin = (ptZ >= 150 && ptZ < 300) ? 1 :
                          (ptZ >= 300 && ptZ < 400) ? 2 :
                          (ptZ >= 400)              ? 3 : 0;
                if      (bin==1) s.h_ZtotalTime_bin1->Fill(ZtotalTime, evweight);
                else if (bin==2) s.h_ZtotalTime_bin2->Fill(ZtotalTime, evweight);
                else if (bin==3) s.h_ZtotalTime_bin3->Fill(ZtotalTime, evweight);
            }

             // groomed mass con taglio n2Fat < 0.18
          if (mFat_groomed >= 0 && n2Fat >= 0 && n2Fat < 0.18 && hofatjet == 1) {
              s.h_mFat_groomed_n2->Fill(mFat_groomed, evweight);
              if      (pFat >= 150 && pFat < 250) s.h_mFat_groomed_n2_pt150_250->Fill(mFat_groomed, evweight);
              else if (pFat >= 250 && pFat < 350) s.h_mFat_groomed_n2_pt250_350->Fill(mFat_groomed, evweight);
              else if (pFat >= 350 && pFat < 500) s.h_mFat_groomed_n2_pt350_500->Fill(mFat_groomed, evweight);
              if (ev_isW == 1 && ev_WMatchedFatJet == 1) s.h_mFat_WMatched_groomed_n2->Fill(mFat_groomed, evweight);
              if (ev_isZ == 1 && ev_ZMatchedFatJet == 1) s.h_mFat_ZMatched_groomed_n2->Fill(mFat_groomed, evweight);
                    }
                    if (mFat_groomed >= 0 && n2Fat_pure >= 0 && n2Fat_pure < 0.18 && hofatjet == 1) {
                        s.h_mFat_groomed_n2pure->Fill(mFat_groomed, evweight);
                      if      (pFat >= 150 && pFat < 250) s.h_mFat_groomed_n2pure_pt150_250->Fill(mFat_groomed, evweight);
                      else if (pFat >= 250 && pFat < 350) s.h_mFat_groomed_n2pure_pt250_350->Fill(mFat_groomed, evweight);
                      else if (pFat >= 350 && pFat < 500) s.h_mFat_groomed_n2pure_pt350_500->Fill(mFat_groomed, evweight);
                      }

              if (ev_isW == 1 && ev_WMatchedFatJet == 1) s.h_mFat_WMatched_groomed->Fill(mFat_groomed, evweight);
              if (ev_isZ == 1 && ev_ZMatchedFatJet == 1) s.h_mFat_ZMatched_groomed->Fill(mFat_groomed, evweight);

              if (mFatpure_groomed >= 0 && n2Fat_pure >= 0 && n2Fat_pure < 0.18 && hofatjet == 1) {
                        s.h_mFatpure_groomed_n2pure->Fill(mFatpure_groomed, evweight);
                      if      (pFat >= 150 && pFat < 250) s.h_mFatpure_groomed_n2pure_pt150_250->Fill(mFatpure_groomed, evweight);
                      else if (pFat >= 250 && pFat < 350) s.h_mFatpure_groomed_n2pure_pt250_350->Fill(mFatpure_groomed, evweight);
                      else if (pFat >= 350 && pFat < 500) s.h_mFatpure_groomed_n2pure_pt350_500->Fill(mFatpure_groomed, evweight);
                      }

        }
    }
}

 
//  printStats / printPtCounts
 
void printStats(const Sample& s) {
    cout << "\n=== STATISTICHE: " << s.name << " ===" << endl;
    cout << "eventi totali (hofatjet==1):  " << s.nTotal               << endl;
    cout << "All con N2 < 0.18:            " << s.nAll_n2cut           << endl;
    cout << "  di cui 150-250 GeV:         " << s.nAll_n2cut_150_250   << endl;
    cout << "  di cui 250-350 GeV:         " << s.nAll_n2cut_250_350   << endl;
    cout << "  di cui 350-500 GeV:         " << s.nAll_n2cut_350_500   << endl;
    cout << "W events:                     " << s.nWall                << endl;
    cout << "W matched (all):              " << s.nWMatched_all         << endl;
    cout << "  di cui 150-250 GeV:         " << s.nWMatched_all_150_250 << endl;
    cout << "  di cui 250-350 GeV:         " << s.nWMatched_all_250_350 << endl;
    cout << "  di cui 350-500 GeV:         " << s.nWMatched_all_350_500 << endl;
    cout << "W matched, N2 < 0.18:         " << s.nWMatched_n2cut       << endl;
    cout << "Z events:                     " << s.nZall                << endl;
    cout << "Z matched (all):              " << s.nZMatched_all         << endl;
    cout << "  di cui 150-250 GeV:         " << s.nZMatched_all_150_250 << endl;
    cout << "  di cui 250-350 GeV:         " << s.nZMatched_all_250_350 << endl;
    cout << "  di cui 350-500 GeV:         " << s.nZMatched_all_350_500 << endl;
    cout << "Z matched, N2 < 0.18:         " << s.nZMatched_n2cut       << endl;
    cout << "W+Z matched 150-250:          " << s.nWMatched_all_150_250 + s.nZMatched_all_150_250 << endl;
    cout << "W+Z matched 250-350:          " << s.nWMatched_all_250_350 + s.nZMatched_all_250_350 << endl;
    cout << "W+Z matched 350-500:          " << s.nWMatched_all_350_500 + s.nZMatched_all_350_500 << endl;
    cout << "=================================" << endl;
}

void printPtCounts(const Sample& s) {
    cout << "\n=== CONTEGGI JETS normalizzati: " << s.name << " ===" << endl;
    cout << "W jets totali:             " << s.h_ptWjet->Integral("width") << endl;
    cout << "W jets 150-250 GeV:        " << s.h_ptWjet_pt150_250->Integral("width") << endl;
    cout << "W jets 250-350 GeV:        " << s.h_ptWjet_pt250_350->Integral("width") << endl;
    cout << "W jets 350-500 GeV:        " << s.h_ptWjet_pt350_500->Integral("width") << endl;
    cout << "Z jets totali:             " << s.h_ptZjet->Integral("width") << endl;
    cout << "Z jets 150-250 GeV:        " << s.h_ptZjet_pt150_250->Integral("width") << endl;
    cout << "Z jets 250-350 GeV:        " << s.h_ptZjet_pt250_350->Integral("width") << endl;
    cout << "Z jets 350-500 GeV:        " << s.h_ptZjet_pt350_500->Integral("width") << endl;
    double wzTot = s.h_ptWjet->Integral("width") + s.h_ptZjet->Integral("width");
    double wz150 = s.h_ptWjet_pt150_250->Integral("width") + s.h_ptZjet_pt150_250->Integral("width");
    double wz250 = s.h_ptWjet_pt250_350->Integral("width") + s.h_ptZjet_pt250_350->Integral("width");
    double wz350 = s.h_ptWjet_pt350_500->Integral("width") + s.h_ptZjet_pt350_500->Integral("width");
    cout << "W+Z jets totali:           " << wzTot << endl;
    cout << "W+Z jets 150-250 GeV:      " << wz150 << endl;
    cout << "W+Z jets 250-350 GeV:      " << wz250 << endl;
    cout << "W+Z jets 350-500 GeV:      " << wz350 << endl;
    cout << "=================================" << endl;
}

 
//  Helpers di plot
 
void normToUnity(TH1* h) { if (h->Integral() > 0) h->Scale(1.0 / h->Integral()); }

Color_t sampleColor(const TString& n) {
    if (n == "quenched")   return kRed+1;
    if (n == "unquenched") return kBlue+1;
    return kGreen+2;
}
Color_t sampleFill(const TString& n) {
    if (n == "quenched")   return kRed-9;
    if (n == "unquenched") return kBlue-9;
    return kGreen-9;
}

void addLabels(TLatex* l, const TString& top) {
    l->DrawLatex(0.18, 0.88, top);
    l->DrawLatex(0.18, 0.84, "|#eta^{jet}| < 1.6");
    l->DrawLatex(0.18, 0.80, "R = 0.8");
}

// ---- 3 campioni sovrapposti ----
void draw3Overlay(TH1D* h1, TH1D* h2, TH1D* h3,
                  const TString& xTitle, const TString& tag, const TString& outDir,
                  const TString& l1, const TString& l2, const TString& l3,
                  const TString& topLabel, bool logy = false, bool isMass = false)
{
    TCanvas* c = new TCanvas("c3_"+tag, "", 900, 800);
    canvasSetup(c);
    double ymax = TMath::Max(h1->GetMaximum(), TMath::Max(h2->GetMaximum(), h3->GetMaximum()));
    if (logy) { gPad->SetLogy(); h1->SetMinimum(1e-3); h1->SetMaximum(ymax*100.); }
    else       { h1->SetMinimum(0);                     h1->SetMaximum(ymax*1.35); }
    TString xlab = isMass ? "Inv. Mass [GeV]" : xTitle;
    h1->SetLineColor(kRed+1);  h2->SetLineColor(kBlue+1);  h3->SetLineColor(kGreen+2);
    histSetup(h1, xlab, "", kRed+1,   false);
    histSetup(h2, xlab, "", kBlue+1,  false);
    histSetup(h3, xlab, "", kGreen+2, false);
    h1->SetLineWidth(2); h2->SetLineWidth(2); h3->SetLineWidth(2);
    TLegend* leg = new TLegend(0.55, 0.72, 0.90, 0.90);
    leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextFont(42); leg->SetTextSize(0.035);
    if (isMass) {
        h1->SetFillColor(kRed-9); h1->SetFillStyle(1001);
        h1->Draw("HIST"); h2->Draw("E1 PE SAME"); h3->Draw("E1 PE SAME");
        leg->AddEntry(h1,l1,"f"); leg->AddEntry(h2,l2,"lep"); leg->AddEntry(h3,l3,"lep");
    } else {
        h1->Draw("E1 PE"); h2->Draw("E1 PE SAME"); h3->Draw("E1 PE SAME");
        leg->AddEntry(h1,l1,"l"); leg->AddEntry(h2,l2,"l"); leg->AddEntry(h3,l3,"l");
    }
    leg->Draw();
    TLatex* latex = new TLatex(); latex->SetNDC(); latex->SetTextFont(42); latex->SetTextSize(0.038);
    addLabels(latex, topLabel);
    CMS_latex();
    c->SaveAs(outDir+"/"+tag+"_overlay_realOx.pdf");
    delete c;
}

// ---- 1 istogramma singolo shaded ----
void drawSingle(TH1D* h, const TString& tag, const TString& outDir,
                const TString& sampleName, const TString& xTitle,
                const TString& topLabel, bool logy = false, bool isMass = false)
{
    TCanvas* c = new TCanvas("cs_"+tag+"_"+sampleName, "", 900, 800);
    canvasSetup(c);
    Color_t col  = sampleColor(sampleName);
    Color_t fill = sampleFill(sampleName);
    h->SetLineColor(col); h->SetFillColor(fill); h->SetFillStyle(1001); h->SetLineWidth(2);
    TString xlab = isMass ? "Inv. Mass [GeV]" : xTitle;
    histSetup(h, xlab, "", col, false);
    if (logy) { gPad->SetLogy(); h->SetMinimum(1e-3); h->SetMaximum(h->GetMaximum()*100.); }
    else       { h->SetMinimum(0);                     h->SetMaximum(h->GetMaximum()*1.35); }
    h->Draw("HIST");
    TLatex* latex = new TLatex(); latex->SetNDC(); latex->SetTextFont(42); latex->SetTextSize(0.038);
    addLabels(latex, topLabel);
    CMS_latex();
    c->SaveAs(outDir+"/"+tag+"_"+sampleName+"_realOx.pdf");
    delete c;
}
// ---- 2 istogrammi sovrapposti per stesso campione ----
void drawOverlap2(TH1D* h1, TH1D* h2,
                  const TString& tag, const TString& outDir,
                  const TString& sampleName,
                  const TString& label1, const TString& label2,
                  const TString& xTitle,const TString& yTitle, const TString& topLabel)
{
    TCanvas* c = new TCanvas("cov2_"+tag+"_"+sampleName, "", 900, 800);
    canvasSetup(c);
    h1->SetLineColor(kBlue+1); h1->SetLineWidth(2);
    h1->SetFillColorAlpha(kBlue+1, 0.5); h1->SetFillStyle(1001);
    h2->SetLineColor(kRed+1);  h2->SetLineWidth(2);
    h2->SetFillColorAlpha(kRed+1, 0.5);  h2->SetFillStyle(1001);
    histSetup(h1, xTitle, yTitle, kBlue+1, false);
    histSetup(h2, xTitle, yTitle, kRed+1,  false);
    double m1 = h1->GetMaximum();
    double m2 = h2->GetMaximum();
    double ymax = TMath::Max(h1->GetMaximum(), h2->GetMaximum());
    h1->SetMinimum(0); h1->SetMaximum(ymax * 2);
    h2->SetMinimum(0); h2->SetMaximum(ymax * 2);

    if (m1 >= m2) {
    h1->Draw("HIST"); h2->Draw("HIST SAME");
} else {
    h2->Draw("HIST"); h1->Draw("HIST SAME");
}
    TLegend* leg = new TLegend(0.15, 0.65, 0.60, 0.80);

    leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextFont(42); leg->SetTextSize(0.035);
    leg->AddEntry(h1, label1, "f");
    leg->AddEntry(h2, label2, "f");
    leg->Draw();
    TLatex* latex = new TLatex(); latex->SetNDC(); latex->SetTextFont(42); latex->SetTextSize(0.038);
    addLabels(latex, topLabel);
    CMS_latex();
    c->SaveAs(outDir+"/"+tag+"_"+sampleName+"_realOx.pdf");
    delete c;
}

void drawOverlap3(TH1D* h1, TH1D* h2, TH1D* h3,
                  const TString& tag, const TString& outDir,
                  const TString& sampleName,
                  const TString& label1, const TString& label2, const TString& label3,
                  const TString& xTitle,const TString& yTitle, const TString& topLabel)
{
    TCanvas* c = new TCanvas("cov3_"+tag+"_"+sampleName, "", 900, 800);
    canvasSetup(c);

    h1->SetLineColor(kBlue+1);   h1->SetLineWidth(2);
    h1->SetFillColorAlpha(kBlue+1,   0.5); h1->SetFillStyle(1001);
    h2->SetLineColor(kRed+1);    h2->SetLineWidth(2);
    h2->SetFillColorAlpha(kRed+1,    0.5); h2->SetFillStyle(1001);
    h3->SetLineColor(kGreen+2);  h3->SetLineWidth(2);
    h3->SetFillColorAlpha(kGreen+2,  0.5); h3->SetFillStyle(1001);

    histSetup(h1, xTitle, yTitle, kBlue+1,  false);
    histSetup(h2, xTitle, yTitle, kRed+1,   false);
    histSetup(h3, xTitle, yTitle, kGreen+2, false);
    
    double m1 = h1->GetMaximum();
    double m2 = h2->GetMaximum();
    double m3 = h3->GetMaximum();

    double ymax = TMath::Max(m1, TMath::Max(m2, m3));

    h1->SetMaximum(ymax*2);
    h2->SetMaximum(ymax*2);
    h3->SetMaximum(ymax*2);
    

    // disegna prima il più alto così gli altri non vengono tagliati
    if (m1 >= m2 && m1 >= m3) {
        h1->Draw("HIST"); h2->Draw("HIST SAME"); h3->Draw("HIST SAME");
    } else if (m2 >= m1 && m2 >= m3 ) {
        h2->Draw("HIST"); h1->Draw("HIST SAME"); h3->Draw("HIST SAME");
    } else {
        h3->Draw("HIST"); h1->Draw("HIST SAME"); h2->Draw("HIST SAME");
    }

    TLegend* leg = new TLegend(0.15, 0.65, 0.60, 0.80);
    leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextFont(42); leg->SetTextSize(0.033);
    leg->AddEntry(h1, label1, "f");
    leg->AddEntry(h2, label2, "f");
    leg->AddEntry(h3, label3, "f");
    leg->Draw();

    TLatex* latex = new TLatex(); latex->SetNDC(); latex->SetTextFont(42); latex->SetTextSize(0.038);
    addLabels(latex, topLabel);
    CMS_latex();
    c->SaveAs(outDir+"/"+tag+"_"+sampleName+"_idealOss.pdf");
    delete c;
}

// ---- 2D per un campione ----
void drawSingle2D(TH2D* h, const TString& tag, const TString& outDir,
                  const TString& sampleName, const TString& topLabel)
{
    TCanvas* c = new TCanvas("c2d_"+tag+"_"+sampleName, "", 1000, 800);
    canvasSetup(c); gPad->SetRightMargin(0.15);
    h->Draw("COLZ");
    TLatex* latex = new TLatex(); latex->SetNDC(); latex->SetTextFont(42); latex->SetTextSize(0.038);
    latex->DrawLatex(0.18, 0.89, topLabel);
    latex->SetTextSize(0.032);
    latex->DrawLatex(0.18, 0.85, "p_{T}^{jet} #geq 150 GeV");
    latex->DrawLatex(0.18, 0.81, "|#eta^{jet}| < 1.6, R = 0.8");
    CMS_latex();
    c->SaveAs(outDir+"/"+tag+"_"+sampleName+"_realOx.pdf");
    delete c;
}

// ---- W+Z matched massa sovrapposti shaded per un campione ----
void drawWZmassSingle(TH1D* hW, TH1D* hZ, const TString& tag, const TString& outDir,
                      const TString& sampleName, const TString& topLabel)
{
    TCanvas* c = new TCanvas("cwz_"+tag+"_"+sampleName, "", 900, 800);
    canvasSetup(c);
    hW->SetLineColor(kBlue+1);  hW->SetFillColorAlpha(kBlue+1,  0.5); hW->SetFillStyle(1001); hW->SetLineWidth(2);
    hZ->SetLineColor(kGreen+2); hZ->SetFillColorAlpha(kGreen+2, 0.5); hZ->SetFillStyle(1001); hZ->SetLineWidth(2);
    histSetup(hW, "Inv. Mass [GeV]", "", kBlue+1,  false);
    histSetup(hZ, "Inv. Mass [GeV]", "", kGreen+2, false);
    double ymax = TMath::Max(hW->GetMaximum(), hZ->GetMaximum());
    hW->SetMinimum(0); hW->SetMaximum(ymax*1.35);
    hW->Draw("HIST"); hZ->Draw("HIST SAME");
    TLegend* leg = new TLegend(0.55, 0.78, 0.90, 0.90);
    leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextFont(42); leg->SetTextSize(0.035);
    leg->AddEntry(hW, "W matched", "f");
    leg->AddEntry(hZ, "Z matched", "f");
    leg->Draw();
    TLatex* latex = new TLatex(); latex->SetNDC(); latex->SetTextFont(42); latex->SetTextSize(0.038);
    addLabels(latex, topLabel);
    CMS_latex();
    c->SaveAs(outDir+"/"+tag+"_"+sampleName+"_realOx.pdf");
    delete c;
}

// ---- 3 curve massa (All, N2cut, WZ matched) shaded per un campione ----
void drawMass3curves(TH1D* hAll, TH1D* hN2, TH1D* hWZ,
                     const TString& tag, const TString& outDir,
                     const TString& sampleName, const TString& topLabel)
{
    TCanvas* c = new TCanvas("cm3_"+tag+"_"+sampleName, "", 900, 800);
    canvasSetup(c);
    hAll->SetLineColor(kRed+1);  hAll->SetFillColorAlpha(kRed+1,   0.5); hAll->SetFillStyle(1001); hAll->SetLineWidth(2);
    hN2->SetLineColor(kBlue+1);  hN2->SetFillColorAlpha(kBlue+1,   0.5); hN2->SetFillStyle(1001);  hN2->SetLineWidth(2);
    hWZ->SetLineColor(kGreen+2); hWZ->SetFillColorAlpha(kGreen+2,  0.5); hWZ->SetFillStyle(1001);  hWZ->SetLineWidth(2);
    histSetup(hAll, "Inv. Mass [GeV]", "", kRed+1,   false);
    histSetup(hN2,  "Inv. Mass [GeV]", "", kBlue+1,  false);
    histSetup(hWZ,  "Inv. Mass [GeV]", "", kGreen+2, false);
    double ymax = TMath::Max(hAll->GetMaximum(), TMath::Max(hN2->GetMaximum(), hWZ->GetMaximum()));
    hAll->SetMinimum(0); hAll->SetMaximum(ymax*1.35);
    hAll->Draw("HIST"); hN2->Draw("HIST SAME"); hWZ->Draw("HIST SAME");
    TLegend* leg = new TLegend(0.42, 0.72, 0.90, 0.90);
    leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextFont(42); leg->SetTextSize(0.033);
    leg->AddEntry(hAll, "All events (no N_{2} cut)",  "f");
    leg->AddEntry(hN2,  "All events, N_{2} < 0.18",   "f");
    leg->AddEntry(hWZ,  "W+Z matched (no N_{2} cut)", "f");
    leg->Draw();
    TLatex* latex = new TLatex(); latex->SetNDC(); latex->SetTextFont(42); latex->SetTextSize(0.038);
    addLabels(latex, topLabel);
    CMS_latex();
    c->SaveAs(outDir+"/"+tag+"_"+sampleName+"_realOx.pdf");
    delete c;
}

// ---- tempo singolo shaded, canvas QUADRATO ----
void drawTimeSingle(TH1D* h, const TString& tag, const TString& outDir,
                    const TString& sampleName, const TString& xTitle,
                    const TString& topLabel)
{
    normToUnity(h);
    TCanvas* c = new TCanvas("ctS_"+tag+"_"+sampleName, "", 800, 800);
    canvasSetup(c);
    Color_t col  = sampleColor(sampleName);
    Color_t fill = sampleFill(sampleName);
    h->SetLineColor(col); h->SetFillColor(fill); h->SetFillStyle(1001); h->SetLineWidth(2);
    histSetup(h, xTitle, "", col, false);
    h->SetMinimum(0); h->SetMaximum(h->GetMaximum()*1.35);
    h->Draw("HIST");
    TLatex* latex = new TLatex(); latex->SetNDC(); latex->SetTextFont(42); latex->SetTextSize(0.038);
    addLabels(latex, topLabel);
    CMS_latex();
    c->SaveAs(outDir+"/"+tag+"_"+sampleName+"_realOx.pdf");
    delete c;
}

// ---- total time con 3 bin di pT sovrapposti shaded, canvas QUADRATO ----
void drawTotalTimeBins(TH1D* hb1, TH1D* hb2, TH1D* hb3,
                       const TString& xTitle, const TString& tag,
                       const TString& outDir, const TString& sampleName,
                       const TString& topLabel)
{
    normToUnity(hb1); normToUnity(hb2); normToUnity(hb3);
    TCanvas* c = new TCanvas("ctB_"+tag+"_"+sampleName, "", 800, 800);
    canvasSetup(c);
    hb1->SetLineColor(kBlue+1);   hb1->SetFillColorAlpha(kBlue+1,   0.5); hb1->SetFillStyle(1001); hb1->SetLineWidth(2);
    hb2->SetLineColor(kGreen+2);  hb2->SetFillColorAlpha(kGreen+2,  0.5); hb2->SetFillStyle(1001); hb2->SetLineWidth(2);
    hb3->SetLineColor(kOrange+1); hb3->SetFillColorAlpha(kOrange+1, 0.5); hb3->SetFillStyle(1001); hb3->SetLineWidth(2);
    histSetup(hb1, xTitle, "", kBlue+1,   false);
    histSetup(hb2, xTitle, "", kGreen+2,  false);
    histSetup(hb3, xTitle, "", kOrange+1, false);
    double ymax = TMath::Max(hb1->GetMaximum(), TMath::Max(hb2->GetMaximum(), hb3->GetMaximum()));
    hb1->SetMinimum(0); hb1->SetMaximum(ymax*1.35);
    hb1->Draw("HIST"); hb2->Draw("HIST SAME"); hb3->Draw("HIST SAME");
    TLegend* leg = new TLegend(0.55, 0.65, 0.90, 0.88);
    leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextFont(42); leg->SetTextSize(0.032);
    leg->AddEntry(hb1, "150 #leq p_{T} < 300 GeV", "f");
    leg->AddEntry(hb2, "300 #leq p_{T} < 400 GeV", "f");
    leg->AddEntry(hb3, "p_{T} #geq 400 GeV",        "f");
    leg->Draw();
    TLatex* latex = new TLatex(); latex->SetNDC(); latex->SetTextFont(42); latex->SetTextSize(0.038);
    latex->DrawLatex(0.18, 0.88, topLabel);
    latex->DrawLatex(0.18, 0.84, sampleName);
    latex->DrawLatex(0.18, 0.80, "|#eta^{jet}| < 1.6, R = 0.8");
    CMS_latex();
    c->SaveAs(outDir+"/"+tag+"_"+sampleName+"_realOx.pdf");
    delete c;
}

 
//  MAIN
 
void plot() {

    gStyle->SetOptStat(0);
    gROOT->SetBatch(kTRUE);


     const TString outDir        = "/eos/home-f/fzavatti/bkgCONSub400SOFT_VB_1TRANCHE_beta034PROVACLEAR";
    gSystem->Exec("mkdir -p " + outDir);

    const TString fl_quenched   = "/eos/home-f/fzavatti/eventi6_R0.2_bkgCONsub400pclesSOFTz02_dq/cucu.txt";
    const TString fl_unquenched = "/eos/home-f/fzavatti/eventi6_R0.2_bkgCONsub400pclesSOFTz02_dq/cucu.txt";
    const TString fl_total      = "/eos/home-f/fzavatti/eventi6_R0.2_bkgCONsub400pclesSOFTz02_dq/cucu.txt";


    const double sigma_W_nb  = 5.392e-2;
    const double sigma_Z_nb  = 2.1611e-2;
    const double lumi_nb_inv = 1000.; //change Lumi for different cases
    const double A           = 16.;  //change to 208 for PbPb
    const double inel_frac   = 0.42;
    const double factorW = sigma_W_nb * lumi_nb_inv * A * A * inel_frac;
    const double factorZ = sigma_Z_nb * lumi_nb_inv * A * A * inel_frac;

    Sample quenched  ("quenched",   fl_quenched,   kRed+1  );
    Sample unquenched("unquenched", fl_unquenched, kBlue+1 );
    Sample total     ("total",      fl_total,      kGreen+2);

    auto normMass = [](Sample& s) {
     
        normToUnity(s.h_n2Fat);
        normToUnity(s.h_n2Fat_pure);

        // --- AGGIUNGI ---
normToUnity(s.h_n2Fat_WMatched);
normToUnity(s.h_n2Fat_ZMatched);
normToUnity(s.h_n2Fat_WZMatched);
normToUnity(s.h_n2Fat_bkg);

// --- AGGIUNGI ---
normToUnity(s.h_n2Fat_WMatched_pure);
normToUnity(s.h_n2Fat_ZMatched_pure);
normToUnity(s.h_n2Fat_WZMatched_pure);
normToUnity(s.h_n2Fat_bkg_pure);








normToUnity(s.h_mFat_WMatched_groomed_n2);
normToUnity(s.h_mFat_ZMatched_groomed_n2);

normToUnity(s.h_mFat_WMatched_groomed);
normToUnity(s.h_mFat_ZMatched_groomed);

double norm_probe = s.h_mFat_All_n2pure->Integral();
double norm_probe_groomed = s.h_mFatpure_groomed_n2pure->Integral();

if (norm_probe > 0) {
    s.h_mFat_All_n2->Scale(1.0 / norm_probe);
    s.h_mFat_All_n2pure->Scale(1.0 / norm_probe);
    s.h_mFat_All_n2pure_mix->Scale(1.0 / norm_probe);

    s.h_mFat_All_n2_pt150_250->Scale(1.0 / norm_probe);
    s.h_mFat_All_n2pure_pt150_250->Scale(1.0 / norm_probe);
    s.h_mFat_All_n2pure_mix_pt150_250->Scale(1.0 / norm_probe);
    s.h_mFat_All_n2_pt250_350->Scale(1.0 / norm_probe);
    s.h_mFat_All_n2pure_pt250_350->Scale(1.0 / norm_probe);
    s.h_mFat_All_n2pure_mix_pt250_350->Scale(1.0 / norm_probe);
    s.h_mFat_All_n2_pt350_500->Scale(1.0 / norm_probe);
    s.h_mFat_All_n2pure_pt350_500->Scale(1.0 / norm_probe);
    s.h_mFat_All_n2pure_mix_pt350_500->Scale(1.0 / norm_probe);
    
    s.h_mFat_groomed_n2->Scale(1.0 / norm_probe_groomed);
    s.h_mFatpure_groomed_n2pure->Scale(1.0 / norm_probe_groomed);
    s.h_mFat_groomed_n2pure->Scale(1.0 / norm_probe_groomed);

    s.h_mFat_groomed_n2_pt150_250->Scale(1.0 / norm_probe_groomed);
    s.h_mFat_groomed_n2_pt250_350->Scale(1.0 / norm_probe_groomed);
    s.h_mFat_groomed_n2_pt350_500->Scale(1.0 / norm_probe_groomed);

    s.h_mFat_groomed_n2pure_pt150_250->Scale(1.0 / norm_probe_groomed);
    s.h_mFat_groomed_n2pure_pt250_350->Scale(1.0 / norm_probe_groomed);
    s.h_mFat_groomed_n2pure_pt350_500->Scale(1.0 / norm_probe_groomed);

    s.h_mFatpure_groomed_n2pure_pt150_250->Scale(1.0 / norm_probe_groomed);
    s.h_mFatpure_groomed_n2pure_pt250_350->Scale(1.0 / norm_probe_groomed);
    s.h_mFatpure_groomed_n2pure_pt350_500->Scale(1.0 / norm_probe_groomed);


}

    };

    auto normPt = [&](Sample& s) {
        if (s.sumW_gen > 0) {
            s.h_ptWjet->Scale(1.0/s.sumW_gen,"width");           s.h_ptWjet->Scale(factorW);
            s.h_ptWjet_pt150_250->Scale(1.0/s.sumW_gen,"width"); s.h_ptWjet_pt150_250->Scale(factorW);
            s.h_ptWjet_pt250_350->Scale(1.0/s.sumW_gen,"width"); s.h_ptWjet_pt250_350->Scale(factorW);
            s.h_ptWjet_pt350_500->Scale(1.0/s.sumW_gen,"width"); s.h_ptWjet_pt350_500->Scale(factorW);
        }
        if (s.sumZ_gen > 0) {
            s.h_ptZjet->Scale(1.0/s.sumZ_gen,"width");           s.h_ptZjet->Scale(factorZ);
            s.h_ptZjet_pt150_250->Scale(1.0/s.sumZ_gen,"width"); s.h_ptZjet_pt150_250->Scale(factorZ);
            s.h_ptZjet_pt250_350->Scale(1.0/s.sumZ_gen,"width"); s.h_ptZjet_pt250_350->Scale(factorZ);
            s.h_ptZjet_pt350_500->Scale(1.0/s.sumZ_gen,"width"); s.h_ptZjet_pt350_500->Scale(factorZ);
        }
    };

    for (Sample* s : {&unquenched, &quenched, &total}) {
        cout << "\n -----------------------------------------" << endl;
        cout << ">>>  " << s->name << endl;
        cout << " -----------------------------------------" << endl;
        fillSample(*s);
        printStats(*s);
        normMass(*s);
        normPt(*s);
        printPtCounts(*s);
    }

    auto makeWZ = [](Sample& s) -> tuple<TH1D*,TH1D*,TH1D*> {
        TH1D* w150 = (TH1D*)s.h_ptWjet_pt150_250->Clone(s.name+"_ptWZjet_150"); w150->Add(s.h_ptZjet_pt150_250);
        TH1D* w250 = (TH1D*)s.h_ptWjet_pt250_350->Clone(s.name+"_ptWZjet_250"); w250->Add(s.h_ptZjet_pt250_350);
        TH1D* w350 = (TH1D*)s.h_ptWjet_pt350_500->Clone(s.name+"_ptWZjet_350"); w350->Add(s.h_ptZjet_pt350_500);
        return {w150, w250, w350};
    };
    auto [qWZ150, qWZ250, qWZ350] = makeWZ(quenched);
    auto [uWZ150, uWZ250, uWZ350] = makeWZ(unquenched);
    auto [tWZ150, tWZ250, tWZ350] = makeWZ(total);

    TH1D* h_ptWZ_q = (TH1D*)quenched.h_ptWjet->Clone("ptWZ_q");   h_ptWZ_q->Add(quenched.h_ptZjet);
    TH1D* h_ptWZ_u = (TH1D*)unquenched.h_ptWjet->Clone("ptWZ_u"); h_ptWZ_u->Add(unquenched.h_ptZjet);
    TH1D* h_ptWZ_t = (TH1D*)total.h_ptWjet->Clone("ptWZ_t");      h_ptWZ_t->Add(total.h_ptZjet);

    auto makeWZmass = [](Sample& s) -> TH1D* {
        TH1D* wz = (TH1D*)s.h_mFat_WMatched_all->Clone(s.name+"_mFat_WZMatched");
        wz->Add(s.h_mFat_ZMatched_all);
        normToUnity(wz);
        return wz;
    };
    TH1D* h_mWZ_q = makeWZmass(quenched);
    TH1D* h_mWZ_u = makeWZmass(unquenched);
    TH1D* h_mWZ_t = makeWZmass(total);

    const TString lQ = "Quenched";
    const TString lU = "Unquenched";
    const TString lT = "Total quenched";

    // ===========================================================
    //  PLOTS SOVRAPPOSTI (3 campioni)
    // ===========================================================
    draw3Overlay(quenched.h_ptWjet,       unquenched.h_ptWjet,       total.h_ptWjet,
                 "p_{T} [GeV]", "ptWjet_total",    outDir, lQ,lU,lT, "W jets",    true);
    draw3Overlay(quenched.h_ptZjet,       unquenched.h_ptZjet,       total.h_ptZjet,
                 "p_{T} [GeV]", "ptZjet_total",    outDir, lQ,lU,lT, "Z jets",    true);
    draw3Overlay(h_ptWZ_q, h_ptWZ_u, h_ptWZ_t,
                 "p_{T} [GeV]", "ptWZjet_total",   outDir, lQ,lU,lT, "W+Z jets", true);
    draw3Overlay(quenched.h_ptWjet_pt150_250, unquenched.h_ptWjet_pt150_250, total.h_ptWjet_pt150_250,
                 "p_{T} [GeV]", "ptWjet_150_250",  outDir, lQ,lU,lT, "W jets, 150<p_{T}<250 GeV", true);
    draw3Overlay(quenched.h_ptWjet_pt250_350, unquenched.h_ptWjet_pt250_350, total.h_ptWjet_pt250_350,
                 "p_{T} [GeV]", "ptWjet_250_350",  outDir, lQ,lU,lT, "W jets, 250<p_{T}<350 GeV", true);
    draw3Overlay(quenched.h_ptWjet_pt350_500, unquenched.h_ptWjet_pt350_500, total.h_ptWjet_pt350_500,
                 "p_{T} [GeV]", "ptWjet_350_500",  outDir, lQ,lU,lT, "W jets, 350<p_{T}<500 GeV", true);
    draw3Overlay(quenched.h_ptZjet_pt150_250, unquenched.h_ptZjet_pt150_250, total.h_ptZjet_pt150_250,
                 "p_{T} [GeV]", "ptZjet_150_250",  outDir, lQ,lU,lT, "Z jets, 150<p_{T}<250 GeV", true);
    draw3Overlay(quenched.h_ptZjet_pt250_350, unquenched.h_ptZjet_pt250_350, total.h_ptZjet_pt250_350,
                 "p_{T} [GeV]", "ptZjet_250_350",  outDir, lQ,lU,lT, "Z jets, 250<p_{T}<350 GeV", true);
    draw3Overlay(quenched.h_ptZjet_pt350_500, unquenched.h_ptZjet_pt350_500, total.h_ptZjet_pt350_500,
                 "p_{T} [GeV]", "ptZjet_350_500",  outDir, lQ,lU,lT, "Z jets, 350<p_{T}<500 GeV", true);
    draw3Overlay(qWZ150, uWZ150, tWZ150,
                 "p_{T} [GeV]", "ptWZjet_150_250", outDir, lQ,lU,lT, "W+Z jets, 150<p_{T}<250 GeV", true);
    draw3Overlay(qWZ250, uWZ250, tWZ250,
                 "p_{T} [GeV]", "ptWZjet_250_350", outDir, lQ,lU,lT, "W+Z jets, 250<p_{T}<350 GeV", true);
    draw3Overlay(qWZ350, uWZ350, tWZ350,
                 "p_{T} [GeV]", "ptWZjet_350_500", outDir, lQ,lU,lT, "W+Z jets, 350<p_{T}<500 GeV", true);

    // ===========================================================
    //  PLOTS SEPARATI per campione
    // ===========================================================
    for (Sample* s : {&quenched, &unquenched, &total}) {
        const TString& nm = s->name;
        TString caseLabel = (nm=="quenched")   ? "" :
                            (nm=="unquenched") ? "" : " ";


        drawSingle(s->h_mFat_All_all,             "mFat_All_all",         outDir, nm, "Inv. Mass [GeV]",
                   "All events (no N_{2} cut)  "+caseLabel,          false, true);
        drawSingle(s->h_mFat_All_n2,           "mFat_All_n2",       outDir, nm, "Inv. Mass [GeV]",
                   "All events, N_{2} < 0.18  "+caseLabel,           false, true);
       
        drawSingle(s->h_mFat_All_n2pure,           "mFat_All_n2pure",       outDir, nm, "Inv. Mass [GeV]",
                   "All events, N_{2}^{probe} < 0.18  "+caseLabel,           false, true);

        drawSingle(s->h_mFat_WMatched_all,        "mFat_WMatched_all",    outDir, nm, "Inv. Mass [GeV]",
                   "W matched (no N_{2} cut)  "+caseLabel,           false, true);
        drawSingle(s->h_mFat_WMatched_n2,      "mFat_WMatched_n2",  outDir, nm, "Inv. Mass [GeV]",
                   "W matched, N_{2} < 0.18  "+caseLabel,            false, true);
        drawSingle(s->h_mFat_ZMatched_all,        "mFat_ZMatched_all",    outDir, nm, "Inv. Mass [GeV]",
                   "Z matched (no N_{2} cut) "+caseLabel,           false, true);
        drawSingle(s->h_mFat_ZMatched_n2,      "mFat_ZMatched_n2",  outDir, nm, "Inv. Mass [GeV]",
                   "Z matched, N_{2} < 0.18  "+caseLabel,            false, true);
        drawSingle(s->h_mFat_All_n2_pt150_250, "mFat_n2_pt150_250", outDir, nm, "Inv. Mass [GeV]",
                   "All, N_{2}<0.18, 150<p_{T}<250 GeV "+caseLabel, false, true);
        drawSingle(s->h_mFat_All_n2_pt250_350, "mFat_n2_pt250_350", outDir, nm, "Inv. Mass [GeV]",
                   "All, N_{2}<0.18, 250<p_{T}<350 GeV  "+caseLabel, false, true);
        drawSingle(s->h_mFat_All_n2_pt350_500, "mFat_n2_pt350_500", outDir, nm, "Inv. Mass [GeV]",
                   "All, N_{2}<0.18, 350<p_{T}<500 GeV  "+caseLabel, false, true);



    
        drawSingle(s->h_mFat_All_n2pure_pt150_250, "mFat_n2pure_pt150_250", outDir, nm, "Inv. Mass [GeV]",
                   "All, N_{2}^{probe}<0.18, 150<p_{T}<250 GeV "+caseLabel, false, true);
        drawSingle(s->h_mFat_All_n2pure_pt250_350, "mFat_n2pure_pt250_350", outDir, nm, "Inv. Mass [GeV]",
                   "All, N_{2}^{probe}<0.18, 250<p_{T}<350 GeV  "+caseLabel, false, true);
        drawSingle(s->h_mFat_All_n2pure_pt350_500, "mFat_n2pure_pt350_500", outDir, nm, "Inv. Mass [GeV]",
                   "All, N_{2}^{probe}<0.18, 350<p_{T}<500 GeV  "+caseLabel, false, true);
           

        drawWZmassSingle(s->h_mFat_WMatched_all,   s->h_mFat_ZMatched_all,
                         "mFat_WZ_matched", outDir, nm,
                         "W and Z matched (no N_{2} cut)  "+caseLabel);
        drawWZmassSingle(s->h_mFat_WMatched_n2, s->h_mFat_ZMatched_n2,
                         "mFat_WZ_n2",   outDir, nm,
                         "W and Z matched, N_{2} < 0.18  "+caseLabel);

       

        drawSingle(s->h_n2Fat, "n2Fat", outDir, nm, "N_{2}",
                   "N_{2} (all events) "+caseLabel, false, false);
        
                          
        drawSingle(s->h_n2Fat_WMatched,  "n2Fat_WMatched",  outDir, nm, "N_{2}",
           "N_{2} W matched jets "+caseLabel,  false, false);
        drawSingle(s->h_n2Fat_ZMatched,  "n2Fat_ZMatched",  outDir, nm, "N_{2}",
           "N_{2} Z matched jets "+caseLabel,  false, false);
        drawSingle(s->h_n2Fat_WZMatched, "n2Fat_WZMatched", outDir, nm, "N_{2}",
           "N_{2} W+Z matched jets "+caseLabel, false, false);
          

           drawOverlap2(s->h_n2Fat, s->h_n2Fat_pure,
             "n2Fat_vs_n2pure", outDir, nm,
             "N_{2} (UE subtracted)", "N_{2} (pure Pythia)",
             "N_{2}", "1/NdN/dN_{2}", "N_{2} comparison " + caseLabel);

        drawOverlap2(s->h_n2Fat_WMatched, s->h_n2Fat_WMatched_pure,
             "n2Fat_vs_n2pure_WMatched", outDir, nm,
             "N_{2} W (UE subtracted)", "N_{2} W (pure Pythia)",
             "N_{2}", "1/NdN/dN_{2}", "N_{2} W matched comparison " + caseLabel);

        drawOverlap2(s->h_n2Fat_ZMatched, s->h_n2Fat_ZMatched_pure,
             "n2Fat_vs_n2pure_ZMatched", outDir, nm,
             "N_{2} Z (UE subtracted)", "N_{2} Z (pure Pythia)",
             "N_{2}", "1/NdN/dN_{2}", "N_{2} Z matched comparison " + caseLabel);

        drawOverlap2(s->h_n2Fat_WZMatched, s->h_n2Fat_WZMatched_pure,
             "n2Fat_vs_n2pure_WZMatched", outDir, nm,
             "N_{2} W+Z (UE subtracted)", "N_{2} W+Z (pure Pythia)",
             "N_{2}", "1/NdN/dN_{2}", "N_{2} W+Z matched comparison " + caseLabel);

       
             
    
        drawOverlap2(s->h_n2Fat_bkg, s->h_n2Fat_bkg_pure,
             "n2Fat_vs_n2pure_bkg", outDir, nm,
             "N_{2} bkg", "N_{2}^{probe} bkg",
             "N_{2}", "1/NdN/dN_{2}", "N_{2} vs N_{2}^{probe} bkg "+caseLabel);

        drawSingle2D(s->h2D_n2_n2pure_WZ, "n2_vs_n2pure_WZ", outDir, nm,"N_{2} vs N_{2}pure   "+caseLabel);
        
        /*drawSingle2D(s->h2D_mFat_n2_W, "mFat_vs_n2_W", outDir, nm,"N_{2} vs W Mass   "+caseLabel);
        drawSingle2D(s->h2D_mFat_n2_Z, "mFat_vs_n2_Z", outDir, nm,"N_{2} vs Z Mass  "+caseLabel);
        drawSingle2D(s->h2D_mFat_n2_bkg, "mFat_vs_n2_bkg", outDir, nm,"N_{2} vs bkg Mass "+caseLabel);

        // n2Fat_pure 2D
        drawSingle2D(s->h2D_mFat_n2pure_W,   "mFat_vs_n2pure_W",   outDir, nm, "N_{2}^{probe} vs W Mass "+caseLabel);
        drawSingle2D(s->h2D_mFat_n2pure_Z,   "mFat_vs_n2pure_Z",   outDir, nm, "N_{2}^{probe} vs Z Mass "+caseLabel);
        drawSingle2D(s->h2D_mFat_n2pure_bkg, "mFat_vs_n2pure_bkg", outDir, nm, "N_{2}^{probe} vs bkg Mass "+caseLabel);*/


drawSingle(s->h_n2Fat_WZMatched_pure, "n2pure_WZmatched", outDir, nm, "N2 pure",
           "W and Z matched "+caseLabel, false, false);



        drawOverlap3(s->h_mFat_All_n2, s->h_mFat_All_n2pure,s->h_mFat_All_n2pure_mix,
             "mFat_3massCompare", outDir, nm,
             "Mass embedded (N_{2}<0.18)",
             "Mass probe (N_{2}^{probe}<0.18)",
             "Mass embedded (N_{2}^{probe}<0.18)",
             "Inv. Mass [GeV]", "1/N_{p}dN/dM",
             "Masses "+caseLabel);


drawOverlap3(s->h_mFat_All_n2_pt150_250, s->h_mFat_All_n2pure_pt150_250, s->h_mFat_All_n2pure_mix_pt150_250,
    "mFat_3massCompare_pt150_250", outDir, nm,
    "Mass embedded (N_{2}<0.18)",
    "Mass probe (N_{2}^{probe}<0.18)",
    "Mass embedded (N_{2}^{probe}<0.18)",
    "Inv. Mass [GeV]", "1/N_{p}dN/dM",
    "Masses 150<p_{T}<250 GeV "+caseLabel);

drawOverlap3(s->h_mFat_All_n2_pt250_350, s->h_mFat_All_n2pure_pt250_350, s->h_mFat_All_n2pure_mix_pt250_350,
    "mFat_3massCompare_pt250_350", outDir, nm,
    "Mass embedded (N_{2}<0.18)",
    "Mass probe (N_{2}^{probe}<0.18)",
    "Mass embedded (N_{2}^{probe}<0.18)",
    "Inv. Mass [GeV]", "1/N_{p}dN/dM",
    "Masses 250<p_{T}<350 GeV "+caseLabel);

drawOverlap3(s->h_mFat_All_n2_pt350_500, s->h_mFat_All_n2pure_pt350_500, s->h_mFat_All_n2pure_mix_pt350_500,
    "mFat_3massCompare_pt350_500", outDir, nm,
    "Mass embedded (N_{2}<0.18)",
    "Mass probe (N_{2}^{probe}<0.18)",
    "Mass embedded (N_{2}^{probe}<0.18)",
    "Inv. Mass [GeV]", "1/N_{p}dN/dM",
    "Masses 350<p_{T}<500 GeV "+caseLabel);

    drawOverlap2(s->h_mFat_All_n2, s->h_mFat_All_n2pure,
             "mFat_vs_mFatpure", outDir, nm,
             "Mass embedded", "Mass probe",
             "Inv. Mass [GeV]", "1/N_{p}dN/dM", "Mass embedded vs probe "+caseLabel);

    drawOverlap2(s->h_mFat_All_n2_pt150_250, s->h_mFat_All_n2pure_pt150_250,
    "mFat_vs_mFatpure_pt150_250", outDir, nm,
    "Mass embedded (N_{2}<0.18)", "Mass probe (N_{2}^{probe}<0.18)",
    "Inv. Mass [GeV]", "1/N_{p}dN/dM", "Mass embedded vs probe, 150<p_{T}<250 GeV "+caseLabel);

drawOverlap2(s->h_mFat_All_n2_pt250_350, s->h_mFat_All_n2pure_pt250_350,
    "mFat_vs_mFatpure_pt250_350", outDir, nm,
    "Mass embedded (N_{2}<0.18)", "Mass probe (N_{2}^{probe}<0.18)",
    "Inv. Mass [GeV]", "1/N_{p}dN/dM", "Mass embedded vs probe, 250<p_{T}<350 GeV "+caseLabel);

drawOverlap2(s->h_mFat_All_n2_pt350_500, s->h_mFat_All_n2pure_pt350_500,
    "mFat_vs_mFatpure_pt350_500", outDir, nm,
    "Mass embedded (N_{2}<0.18)", "Mass probe (N_{2}^{probe}<0.18)",
    "Inv. Mass [GeV]",  "1/N_{p}dN/dM","Mass embedded vs probe, 350<p_{T}<500 GeV "+caseLabel);

// ---- confronto a 2: groomed(n2) vs groomed pure(n2pure) ----
//  //ex: mFat_groomed_vs_pure
drawOverlap2(s->h_mFat_groomed_n2, s->h_mFatpure_groomed_n2pure,
             "mFat_vs_pure_GROOMED", outDir, nm,
             "Groomed embedded mass (N_{2}<0.18)", "Groomed Probe mass (N_{2}^{probe}<0.18)",
             "Inv. Mass [GeV]", "1/N_{p}dN/dM", "Groomed vs pure "+caseLabel);
             


// ---- confronto a 3: groomed(n2) vs groomed pure(n2pure) vs groomed(n2pure) ----
drawOverlap3(s->h_mFat_groomed_n2, s->h_mFatpure_groomed_n2pure, s->h_mFat_groomed_n2pure,
             "mFat_3compareGROOMED", outDir, nm,
             "Groomed embedded mass (N_{2}<0.18)",
             "Groomed probe mass (N_{2}^{probe}<0.18)",
             "Groomed embedded mass (N_{2}^{probe}<0.18)",
             "Inv. Mass [GeV]", "1/N_{p}dN/dM",
             "Masses"+caseLabel);

drawOverlap3(s->h_mFat_groomed_n2_pt150_250, s->h_mFatpure_groomed_n2pure_pt150_250, s->h_mFat_groomed_n2pure_pt150_250,
    "mFat_3massCompareGROOMED_pt150_250", outDir, nm,
     "Groomed embedded mass (N_{2}<0.18)",
     "Groomed probe mass (N_{2}^{probe}<0.18)",
    "Groomed embedded mass (N_{2}^{probe}<0.18)",
    "Inv. Mass [GeV]", "1/N_{p}dN/dM",
    "Masses 150<p_{T}<250 GeV "+caseLabel);

drawOverlap3(s->h_mFat_groomed_n2_pt250_350, s->h_mFatpure_groomed_n2pure_pt250_350, s->h_mFat_groomed_n2pure_pt250_350,
    "mFat_3massCompareGROOMED_pt250_350", outDir, nm,
    "Groomed embedded mass (N_{2}<0.18)",
    "Groomed probe mass (N_{2}^{probe}<0.18)",
    "Groomed embedded mass (N_{2}^{probe}<0.18)",
    "Inv. Mass [GeV]", "1/N_{p}dN/dM",
    "Masses 250<p_{T}<350 GeV "+caseLabel);

drawOverlap3(s->h_mFat_groomed_n2_pt350_500, s->h_mFatpure_groomed_n2pure_pt350_500, s->h_mFat_groomed_n2pure_pt350_500,
    "mFat_3massCompareGROOMED_pt350_500", outDir, nm,
    "Groomed embedded mass (N_{2}<0.18)",
    "Groomed probe mass (N_{2}^{probe}<0.18)",
    "Groomed embedded mass (N_{2}^{probe}<0.18)",
    "Inv. Mass [GeV]", "1/N_{p}dN/dM",
    "Masses 350<p_{T}<500 GeV "+caseLabel);

    drawOverlap2(s->h_mFat_groomed_n2_pt150_250, s->h_mFatpure_groomed_n2pure_pt150_250,
    "mFat_vs_mFatpureGROOMED_pt150_250", outDir, nm,
    "Groomed embedded mass (N_{2}<0.18)", "Groomed probe mass (N_{2}^{probe}<0.18)",
    "Inv. Mass [GeV]", "1/N_{p}dN/dM", "Mass embedded vs pure, 150<p_{T}<250 GeV "+caseLabel);

drawOverlap2(s->h_mFat_groomed_n2_pt250_350, s->h_mFatpure_groomed_n2pure_pt250_350,
    "mFat_vs_mFatpureGROOMED_pt250_350", outDir, nm,
    "Groomed embedded mass (N_{2}<0.18)", "Groomed probe mass (N_{2}^{probe}<0.18)",
    "Inv. Mass [GeV]", "1/N_{p}dN/dM", "Mass embedded vs pure, 250<p_{T}<350 GeV "+caseLabel);

drawOverlap2(s->h_mFat_groomed_n2_pt350_500, s->h_mFatpure_groomed_n2pure_pt350_500,
    "mFat_vs_mFatpureGROOMED_pt350_500", outDir, nm,
    "Groomed embedded mass (N_{2}<0.18)", "Groomed probe mass (N_{2}^{probe}<0.18)",
    "Inv. Mass [GeV]",  "1/N_{p}dN/dM","Mass embedded vs pure, 350<p_{T}<500 GeV "+caseLabel);



       

             // ---- singoli groomed con n2 < 0.18 ----
drawSingle(s->h_mFat_groomed_n2, "mFat_groomed_n2", outDir, nm, "Inv. Mass [GeV]",
           "Groomed mass, N_{2}<0.18 "+caseLabel, false, true);
drawSingle(s->h_mFat_groomed_n2_pt150_250, "mFat_groomed_n2_pt150_250", outDir, nm, "Inv. Mass [GeV]",
           "Groomed mass, N_{2}<0.18, 150<p_{T}<250 GeV "+caseLabel, false, true);
drawSingle(s->h_mFat_groomed_n2_pt250_350, "mFat_groomed_n2_pt250_350", outDir, nm, "Inv. Mass [GeV]",
           "Groomed mass, N_{2}<0.18, 250<p_{T}<350 GeV "+caseLabel, false, true);
drawSingle(s->h_mFat_groomed_n2_pt350_500, "mFat_groomed_n2_pt350_500", outDir, nm, "Inv. Mass [GeV]",
           "Groomed mass, N_{2}<0.18, 350<p_{T}<500 GeV "+caseLabel, false, true);

// ---- W+Z matched ----
drawWZmassSingle(s->h_mFat_WMatched_groomed_n2, s->h_mFat_ZMatched_groomed_n2,
                 "mFat_WZ_groomed_n2", outDir, nm,
                 "Groomed mass W+Z matched, N_{2}<0.18 "+caseLabel);
drawWZmassSingle(s->h_mFat_WMatched_groomed, s->h_mFat_ZMatched_groomed,
                 "mFat_WZ_groomed", outDir, nm,
                 "Groomed mass W+Z matched "+caseLabel);
           

           

                        }

    cout << "\nDone! Plots saved in: " << outDir << endl;
}