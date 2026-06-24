// fit2.C
// PHASE 2: W signal fit only (f=1, g=1).
// Fits the W mass peak with a Crystal Ball function with left-side
// exponential tail, on the W-matched MC truth sample.
// Background parameters are fixed from Phase 1.
// Results (W Crystal Ball parameters) are saved and used in Phase 4.
//
// pT bins: pt150_350, pt350_500, pt500_700, pt700_1000
// To add a new bin: add an entry in GetConfig_Phase2() and in creating.C
//
// Usage:
//   EOS=/your/output/path
//   root -l -q "fit2.C(\"$EOS/WZ_histograms4bins.root\",\"$EOS/figuresB/WZfit_phase1_pt150_350/phase1_results.root\",\"WZfit_phase2\",\"pt150_350\")"
//   (repeat for each pT bin)
//
// NOTE: update EOS_BASE below to match your output directory.

#include <TFile.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TMath.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TROOT.h>
#include <TSystem.h>

using namespace std;

// NOTE: change EOS_BASE to match your output directory
#define EOS_BASE "/eos/home-f/fzavatti/eventi6_R0.2_bkgCONsub400pclesSOFTz02_tot"

struct FitConfig_Phase2 {
    Double_t A_init;
    Double_t fgauss_W_init;
    Double_t mean_W_init;
    Double_t sigma1_W_init;
    Double_t sigma2_W_init;
    Double_t lambda_W_init;
    Double_t trans_W_init;

    Double_t A_min,         A_max;
    Double_t fgauss_W_min,  fgauss_W_max;
    Double_t mean_W_min,    mean_W_max;
    Double_t sigma1_W_min,  sigma1_W_max;
    Double_t sigma2_W_min,  sigma2_W_max;
    Double_t lambda_W_min,  lambda_W_max;
    Double_t trans_W_min,   trans_W_max;

    Int_t rebin;
    TString pt_label;
};

FitConfig_Phase2 GetConfig_Phase2(TString pt_bin) {
    FitConfig_Phase2 cfg;

    if (pt_bin == "all") {
        cfg.A_init          = 9.37291e-05;
        cfg.fgauss_W_init   = 0.5;
        cfg.mean_W_init     = 73.4;
        cfg.sigma1_W_init   = 2.0;
        cfg.sigma2_W_init   = 2.0;
        cfg.lambda_W_init   = 0.4;
        cfg.trans_W_init    = 68.0;
        cfg.A_min = 0;           cfg.A_max = 1.0e-02;
        cfg.fgauss_W_min = 0;    cfg.fgauss_W_max = 1;
        cfg.mean_W_min = 65;     cfg.mean_W_max = 85;
        cfg.sigma1_W_min = 1;    cfg.sigma1_W_max = 8;
        cfg.sigma2_W_min = 1;    cfg.sigma2_W_max = 8;
        cfg.lambda_W_min = 0;    cfg.lambda_W_max = 1;
        cfg.trans_W_min  = 60;   cfg.trans_W_max  = 76;
        cfg.rebin    = 1;
        cfg.pt_label = "all p_{T}";
    }
    else if (pt_bin == "pt150_350") {
        cfg.A_init          = 3.37291e-03;
        cfg.fgauss_W_init   = 0.5;
        cfg.mean_W_init     = 73.0;
        cfg.sigma1_W_init   = 2.5;
        cfg.sigma2_W_init   = 2.0;
        cfg.lambda_W_init   = 0.4;
        cfg.trans_W_init    = 68.0;
        cfg.A_min = 0;           cfg.A_max = 1.0e-02;
        cfg.fgauss_W_min = 0;    cfg.fgauss_W_max = 1;
        cfg.mean_W_min = 60;     cfg.mean_W_max = 90;
        cfg.sigma1_W_min = 1;    cfg.sigma1_W_max = 3.5;
        cfg.sigma2_W_min = 1;    cfg.sigma2_W_max = 3.5;
        cfg.lambda_W_min = 0;    cfg.lambda_W_max = 1;
        cfg.trans_W_min  = 60;   cfg.trans_W_max  = 76;
        cfg.rebin    = 1;
        cfg.pt_label = "150 #leq p_{T} < 350 GeV";
    }
    else if (pt_bin == "pt350_500") {
        cfg.A_init          = 2.37291e-04;
        cfg.fgauss_W_init   = 0.5;
        cfg.mean_W_init     = 79.0;
        cfg.sigma1_W_init   = 1.5;
        cfg.sigma2_W_init   = 4.0;
        cfg.lambda_W_init   = 0.4;
        cfg.trans_W_init    = 68.0;
        cfg.A_min = 0;           cfg.A_max = 1.0e-02;
        cfg.fgauss_W_min = 0;    cfg.fgauss_W_max = 1;
        cfg.mean_W_min = 65;     cfg.mean_W_max = 82;
        cfg.sigma1_W_min = 1;    cfg.sigma1_W_max = 4;
        cfg.sigma2_W_min = 1;    cfg.sigma2_W_max = 4;
        cfg.lambda_W_min = 0;    cfg.lambda_W_max = 1;
        cfg.trans_W_min  = 60;   cfg.trans_W_max  = 76;
        cfg.rebin    = 1;
        cfg.pt_label = "350 #leq p_{T} < 500 GeV";
    }
    
    else if (pt_bin == "pt500_700") {
        cfg.A_init          = 1.5e-06;
        cfg.fgauss_W_init   = 0.5;
        cfg.mean_W_init     = 79.5;
        cfg.sigma1_W_init   = 2.5;
        cfg.sigma2_W_init   = 2.5;
        cfg.lambda_W_init   = 0.4;
        cfg.trans_W_init    = 74.0;
        cfg.A_min = 0;           cfg.A_max = 1.0e-02;
        cfg.fgauss_W_min = 0;    cfg.fgauss_W_max = 1;
        cfg.mean_W_min = 75;     cfg.mean_W_max = 85;
        cfg.sigma1_W_min = 0.5;  cfg.sigma1_W_max =3.;
        cfg.sigma2_W_min = 0.5;  cfg.sigma2_W_max =3;
        cfg.lambda_W_min = 0;    cfg.lambda_W_max = 1;
        cfg.trans_W_min  = 70.5;   cfg.trans_W_max  = 75;
        cfg.rebin    = 1;
        cfg.pt_label = "500 #leq p_{T} < 700 GeV";
    }
    else if (pt_bin == "pt700_1000") {
        cfg.A_init          = 1.0e-06;
        cfg.fgauss_W_init   = 0.5;
        cfg.mean_W_init     = 73.5;
        cfg.sigma1_W_init   = 1.5;
        cfg.sigma2_W_init   = 2.0;
        cfg.lambda_W_init   = 0.4;
        cfg.trans_W_init    = 68.0;
        cfg.A_min = 0;           cfg.A_max = 1.0e-02;
        cfg.fgauss_W_min = 0;    cfg.fgauss_W_max = 1;
        cfg.mean_W_min = 65;     cfg.mean_W_max = 82;
        cfg.sigma1_W_min = 0.5;  cfg.sigma1_W_max = 2.5;
        cfg.sigma2_W_min = 0.5;  cfg.sigma2_W_max = 2.5;
        cfg.lambda_W_min = 0;    cfg.lambda_W_max = 1;
        cfg.trans_W_min  = 60;   cfg.trans_W_max  = 70;
        cfg.rebin    = 1;
        cfg.pt_label = "700 #leq p_{T} < 1000 GeV";
    }
    else {
        cout << "ERROR: unrecognised pt_bin: " << pt_bin << endl;
    }
    return cfg;
}

// BACKGROUND: pol. 4 grade

Double_t background_func(Double_t *x, Double_t *p) {
    Double_t m = x[0];
    return (1.0 + p[0]*m + p[1]*m*m + p[2]*m*m*m + p[3]*m*m*m*m);
}

// W SIGNAL: Crystal Ball with left tail
//   p[0]=fgauss_W, p[1]=mean_W, p[2]=sigma1_W, p[3]=sigma2_W,
//   p[4]=lambda_W, p[5]=trans_W

Double_t W_signal_func(Double_t *x, Double_t *p) {
    Double_t G_trans = p[0] * TMath::Gaus(p[5], p[1], p[2], kFALSE) +
                       (1.0 - p[0]) * TMath::Gaus(p[5], p[1], p[3], kFALSE);
    if (x[0] > p[5]) return p[0] * TMath::Gaus(x[0], p[1], p[2], kFALSE) +
                            (1.0 - p[0]) * TMath::Gaus(x[0], p[1], p[3], kFALSE);
    else              return G_trans * TMath::Exp(p[4] * (x[0] - p[5]));
}

// Z SIGNAL: Crystal Ball with left tail

Double_t Z_signal_func(Double_t *x, Double_t *p) {
    Double_t G_trans = p[0] * TMath::Gaus(p[5], p[1], p[2], kFALSE) +
                       (1.0 - p[0]) * TMath::Gaus(p[5], p[1], p[3], kFALSE);
    if (x[0] > p[5]) return p[0] * TMath::Gaus(x[0], p[1], p[2], kFALSE) +
                            (1.0 - p[0]) * TMath::Gaus(x[0], p[1], p[3], kFALSE);
    else              return G_trans * TMath::Exp(p[4] * (x[0] - p[5]));
}

// TOT FUNCTION: 19 parameters
//   p[0]       = A            (amplitude )
//   p[1]       = f            (sig fraction, fisso = 0 in fase 1)
//   p[2]       = g            (fraction W vs Z,  fisso = 0 in fase 1)
//   p[3..8]    = W signal     (fgauss_W, mean_W, sigma1_W, sigma2_W, lambda_W, trans_W)
//   p[9..14]   = Z signal     (fgauss_Z, mean_Z, sigma1_Z, sigma2_Z, lambda_Z, trans_Z)
//   p[15..18]  = background   (p1, p2, p3, p4)

Double_t total_func(Double_t *x, Double_t *p) {
    Double_t W_pars[6]   = {p[3],  p[4],  p[5],  p[6],  p[7],  p[8]};
    Double_t Z_pars[6]   = {p[9],  p[10], p[11], p[12], p[13], p[14]};
    Double_t bkg_pars[4] = {p[15], p[16], p[17], p[18]};
    Double_t W   = W_signal_func(x, W_pars);
    Double_t Z   = Z_signal_func(x, Z_pars);
    Double_t bkg = background_func(x, bkg_pars);
    return p[0] * (p[1] * (p[2]*W + (1.0-p[2])*Z) + (1.0-p[1])*bkg);
}

void CMS_latex() {
    TLatex *cms = new TLatex();
    cms->SetTextSize(0.04); cms->SetTextFont(61); cms->SetNDC();
    cms->DrawLatex(0.22, 0.92, "CMS");
    TLatex *sim = new TLatex();
    sim->SetTextSize(0.035); sim->SetTextFont(52); sim->SetNDC();
    sim->DrawLatex(0.30, 0.92, "Simulation Preliminary");
    TLatex *lumi = new TLatex();
    lumi->SetTextSize(0.03); lumi->SetTextFont(42); lumi->SetNDC();
    lumi->SetTextAlign(31);
    lumi->DrawLatex(0.90, 0.92, "Pythia 8");
}

void fit2(
    TString input_file  = EOS_BASE "/WZ_histograms4bins.root",
    TString phase1_file = "auto",
    TString output_tag  = "WZfit_phase2",
    TString pt_bin      = "all",
    Double_t fit_min    = 60.0,
    Double_t fit_max    = 83.5

) {
    

    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TString full_tag = (pt_bin == "all") ? output_tag
                                         : Form("%s_%s", output_tag.Data(), pt_bin.Data());
    TString out_dir = Form(EOS_BASE "/figuresB/%s", full_tag.Data());
    gSystem->Exec(Form("mkdir -p %s", out_dir.Data()));

    TString ph1_tag  = (pt_bin == "all") ? "WZfit_phase1"
                                         : Form("WZfit_phase1_%s", pt_bin.Data());
    TString ph1_path = Form(EOS_BASE "/figuresB/%s/phase1_results.root", ph1_tag.Data());
    if (phase1_file != "auto") ph1_path = phase1_file;

    cout << "Reading phase 1 from: " << ph1_path << endl;

    TFile* f_phase1 = TFile::Open(ph1_path);
    if (!f_phase1 || f_phase1->IsZombie()) {
        cout <<"ERROR: Phase 1 file not found: " << ph1_path << endl;
        return;
    }
    TH1F* h_phase1 = (TH1F*)f_phase1->Get("results");
    if (!h_phase1) { cout <<"ERROR: Phase 1 results not found!" << endl; return; }

    Double_t p1_bkg = h_phase1->GetBinContent(1);
    Double_t p2_bkg = h_phase1->GetBinContent(2);
    Double_t p3_bkg = h_phase1->GetBinContent(3);
    Double_t p4_bkg = h_phase1->GetBinContent(4);

    
    f_phase1->Close();

    TFile* f = TFile::Open(input_file);
    if (!f || f->IsZombie()) { cout << "ERROR: " << input_file << endl; return; }

    TString hname = (pt_bin == "all") ? "mFat_n2_2DW"
                                      : Form("mFat_n2_2DW_%s", pt_bin.Data());
    TH2F* h2_sig = (TH2F*)f->Get(hname);
    if (!h2_sig) { cout << "ERROR: Histogram not found:" << hname << endl; f->ls(); return; }

    Double_t N2_cut  = 0.18;
    int binN2_min    = h2_sig->GetYaxis()->FindBin(0.0);
    int binN2_max    = h2_sig->GetYaxis()->FindBin(N2_cut);
    TH1D* hist_sig   = h2_sig->ProjectionX("hist_sig", binN2_min, binN2_max);

    FitConfig_Phase2 cfg = GetConfig_Phase2(pt_bin);
    if (cfg.rebin > 1) hist_sig->Rebin(cfg.rebin);
    Double_t binwidth = hist_sig->GetBinWidth(1);
    hist_sig->SetYTitle(Form("Events / %.1f GeV", binwidth));

   

    TF1* fit_phase2 = new TF1("fit_phase2", total_func, fit_min, fit_max, 19);
    fit_phase2->SetParName(0,  "A");
    fit_phase2->SetParName(1,  "f");
    fit_phase2->SetParName(2,  "g");
    fit_phase2->SetParName(3,  "fgauss_W");
    fit_phase2->SetParName(4,  "mean_W");
    fit_phase2->SetParName(5,  "sigma1_W");
    fit_phase2->SetParName(6,  "sigma2_W");
    fit_phase2->SetParName(7,  "lambda_W");
    fit_phase2->SetParName(8,  "trans_W");
    fit_phase2->SetParName(9,  "fgauss_Z");
    fit_phase2->SetParName(10, "mean_Z");
    fit_phase2->SetParName(11, "sigma1_Z");
    fit_phase2->SetParName(12, "sigma2_Z");
    fit_phase2->SetParName(13, "lambda_Z");
    fit_phase2->SetParName(14, "trans_Z");
    fit_phase2->SetParName(15, "p1_bkg");
    fit_phase2->SetParName(16, "p2_bkg");
    fit_phase2->SetParName(17, "p3_bkg");
    fit_phase2->SetParName(18, "p4_bkg");

    fit_phase2->SetParameter(0, cfg.A_init);
    fit_phase2->SetParLimits(0, cfg.A_min, cfg.A_max);
    fit_phase2->FixParameter(1, 1.0);   // f=1: solo segnale
    fit_phase2->FixParameter(2, 1.0);   // g=1: solo W

    // W (to fit)
    fit_phase2->SetParameter(3, cfg.fgauss_W_init);
    fit_phase2->SetParameter(4, cfg.mean_W_init);
    fit_phase2->SetParameter(5, cfg.sigma1_W_init);
    fit_phase2->SetParameter(6, cfg.sigma2_W_init);
    fit_phase2->SetParameter(7, cfg.lambda_W_init);
    fit_phase2->SetParameter(8, cfg.trans_W_init);
    fit_phase2->SetParLimits(3, cfg.fgauss_W_min,  cfg.fgauss_W_max);
    fit_phase2->SetParLimits(4, cfg.mean_W_min,    cfg.mean_W_max);
    fit_phase2->SetParLimits(5, cfg.sigma1_W_min,  cfg.sigma1_W_max);
    fit_phase2->SetParLimits(6, cfg.sigma2_W_min,  cfg.sigma2_W_max);
    fit_phase2->SetParLimits(7, cfg.lambda_W_min,  cfg.lambda_W_max);
    fit_phase2->SetParLimits(8, cfg.trans_W_min,   cfg.trans_W_max);

    // Z dummy (fixed, g=1)
    fit_phase2->FixParameter(9,  0.5);
    fit_phase2->FixParameter(10, 91.0);
    fit_phase2->FixParameter(11, 2.0);
    fit_phase2->FixParameter(12, 2.0);
    fit_phase2->FixParameter(13, 0.4);
    fit_phase2->FixParameter(14, 85.0);

    // Background pol4 fixed with phase 1
    fit_phase2->FixParameter(15, p1_bkg);
    fit_phase2->FixParameter(16, p2_bkg);
    fit_phase2->FixParameter(17, p3_bkg);
    fit_phase2->FixParameter(18, p4_bkg);

    fit_phase2->SetNpx(1000);

    cout << "\n RUNNING W FIT (Crystal Ball) " << endl;
    TFitResultPtr result = hist_sig->Fit(fit_phase2, "RSWL", "", fit_min, fit_max);

    Double_t A            = fit_phase2->GetParameter(0);
    Double_t err_A        = fit_phase2->GetParError(0);
    Double_t fgauss       = fit_phase2->GetParameter(3);
    Double_t err_fgauss   = fit_phase2->GetParError(3);
    Double_t mass_sig     = fit_phase2->GetParameter(4);
    Double_t err_mass_sig = fit_phase2->GetParError(4);
    Double_t sigma1_sig   = fit_phase2->GetParameter(5);
    Double_t err_sigma1   = fit_phase2->GetParError(5);
    Double_t sigma2_sig   = fit_phase2->GetParameter(6);
    Double_t err_sigma2   = fit_phase2->GetParError(6);
    Double_t lambda_sig   = fit_phase2->GetParameter(7);
    Double_t err_lambda   = fit_phase2->GetParError(7);
    Double_t trans_sig    = fit_phase2->GetParameter(8);
    Double_t err_trans    = fit_phase2->GetParError(8);

    // Chi2 manuale
    Double_t chi2_manual = 0.0;
    Int_t    npoints     = 0;
    for (int i = 1; i <= hist_sig->GetNbinsX(); i++) {
        Double_t x   = hist_sig->GetBinCenter(i);
        if (x < fit_min || x > fit_max) continue;
        Double_t obs = hist_sig->GetBinContent(i);
        Double_t err = hist_sig->GetBinError(i);
        Double_t exp = fit_phase2->Eval(x);
        if (err > 0) {
            chi2_manual += (obs - exp)*(obs - exp) / (err*err);
            npoints++;
        }
    }
    Int_t    ndf_manual     = npoints - fit_phase2->GetNumberFreeParameters();
    Double_t chi2ndf_manual = (ndf_manual > 0) ? chi2_manual / ndf_manual : -1;

    cout << "\n========================================" << endl;
    cout << "PHASE 2 RESULTS (W Crystal Ball) [" << pt_bin << "]" << endl;
    cout << "========================================" << endl;
    cout << "Chi2 manual / ndf = " << chi2_manual << " / " << ndf_manual
         << " = " << chi2ndf_manual << endl;
    cout << "A        = " << A          << " +/- " << err_A        << endl;
    cout << "f_gauss1 = " << fgauss     << " +/- " << err_fgauss   << endl;
    cout << "Mean_W   = " << mass_sig   << " +/- " << err_mass_sig << " GeV" << endl;
    cout << "Sigma1_W = " << sigma1_sig << " +/- " << err_sigma1   << " GeV" << endl;
    cout << "Sigma2_W = " << sigma2_sig << " +/- " << err_sigma2   << " GeV" << endl;
    cout << "Lambda_W = " << lambda_sig << " +/- " << err_lambda   << endl;
    cout << "Trans_W  = " << trans_sig  << " +/- " << err_trans    << " GeV" << endl;
    cout << "========================================\n" << endl;

    TCanvas* c = new TCanvas("c", "Phase 2, W Signal Fit (Crystal Ball)", 900, 800);
    c->SetTicks();
    c->SetLeftMargin(0.14); c->SetRightMargin(0.05);
    c->SetTopMargin(0.08);  c->SetBottomMargin(0.12);

    hist_sig->SetMarkerStyle(20); hist_sig->SetMarkerSize(0.9);
    hist_sig->SetMarkerColor(kBlack); hist_sig->SetLineColor(kBlack);
    hist_sig->SetXTitle("Inv. Mass [GeV]");
    hist_sig->GetXaxis()->SetRangeUser(65, 100);
    hist_sig->GetXaxis()->SetTitleSize(0.045); hist_sig->GetXaxis()->SetLabelSize(0.04);
    hist_sig->GetYaxis()->SetTitleSize(0.045); hist_sig->GetYaxis()->SetLabelSize(0.04);
    hist_sig->SetMinimum(0);
    hist_sig->SetMaximum(hist_sig->GetMaximum() * 1.4);
    hist_sig->Draw("PE");

    fit_phase2->SetLineColor(kRed+1); fit_phase2->SetLineWidth(3);
    fit_phase2->Draw("SAME");

    TLegend* leg1 = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg1->SetBorderSize(0); leg1->SetFillStyle(0);
    leg1->SetTextSize(0.032); leg1->SetTextFont(42);
    leg1->AddEntry(hist_sig,   "W Signal only",    "PE");
    leg1->AddEntry(fit_phase2, "Crystal Ball Fit", "L");
    leg1->Draw();

    TLegend* leg2 = new TLegend(0.52, 0.40, 0.92, 0.73);
    leg2->SetBorderSize(0); leg2->SetFillStyle(0);
    leg2->SetTextSize(0.028); leg2->SetTextFont(42); leg2->SetTextAlign(12);
    leg2->AddEntry((TObject*)0, Form("#chi^{2}/ndf = %.1f / %d = %.2f", chi2_manual, ndf_manual, chi2ndf_manual), "");
    leg2->AddEntry((TObject*)0, Form("A = %.2e #pm %.2e", A, err_A), "");
    leg2->AddEntry((TObject*)0, Form("f_{gauss1} = %.3f #pm %.3f", fgauss, err_fgauss), "");
    leg2->AddEntry((TObject*)0, Form("m_{W} = %.2f #pm %.2f GeV", mass_sig, err_mass_sig), "");
    leg2->AddEntry((TObject*)0, Form("#sigma_{1} = %.2f #pm %.2f GeV", sigma1_sig, err_sigma1), "");
    leg2->AddEntry((TObject*)0, Form("#sigma_{2} = %.2f #pm %.2f GeV", sigma2_sig, err_sigma2), "");
    leg2->AddEntry((TObject*)0, Form("#lambda_{W} = %.3f #pm %.3f", lambda_sig, err_lambda), "");
    leg2->AddEntry((TObject*)0, Form("x_{0}^{W} = %.2f #pm %.2f GeV", trans_sig, err_trans), "");
    leg2->Draw();

    CMS_latex();

    TLatex* info = new TLatex();
    info->SetNDC(); info->SetTextSize(0.030); info->SetTextFont(42);
    info->DrawLatex(0.17, 0.87, "|#eta^{jet}| < 1.6");
    info->DrawLatex(0.17, 0.83, "R = 0.8");
    info->DrawLatex(0.17, 0.79, "N_{2} < 0.18");
    info->DrawLatex(0.17, 0.75, cfg.pt_label);

    c->SaveAs(Form("%s/phase2_W_fit.pdf", out_dir.Data()));
    cout << "Plot saved: " << out_dir << "/phase2_W_fit.pdf" << endl;

    // Saving results
    // bin1=A, bin2=f_signal(1), bin3=fgauss_W, bin4=Mean_W, bin5=Sigma1_W,
    // bin6=Sigma2_W, bin7=Lambda_W, bin8=Trans_W, bin9=Chi2, bin10=NDF

    TFile* output = new TFile(Form("%s/phase2_results.root", out_dir.Data()), "RECREATE");
    hist_sig->Write("hist_sig");
    fit_phase2->Write("fit_phase2");

    TH1F* h_results = new TH1F("results", "Phase 2 Results", 14, 0, 14);
    h_results->SetBinContent(1, A);          h_results->SetBinError(1, err_A);
    h_results->SetBinContent(2, 1.0);
    h_results->SetBinContent(3, fgauss);     h_results->SetBinError(3, err_fgauss);
    h_results->SetBinContent(4, mass_sig);   h_results->SetBinError(4, err_mass_sig);
    h_results->SetBinContent(5, sigma1_sig); h_results->SetBinError(5, err_sigma1);
    h_results->SetBinContent(6, sigma2_sig); h_results->SetBinError(6, err_sigma2);
    h_results->SetBinContent(7, lambda_sig); h_results->SetBinError(7, err_lambda);
    h_results->SetBinContent(8, trans_sig);  h_results->SetBinError(8, err_trans);
    h_results->SetBinContent(9, chi2_manual);
    h_results->SetBinContent(10, ndf_manual);
    h_results->GetXaxis()->SetBinLabel(1,  "A");
    h_results->GetXaxis()->SetBinLabel(2,  "f_signal");
    h_results->GetXaxis()->SetBinLabel(3,  "f_gauss1");
    h_results->GetXaxis()->SetBinLabel(4,  "Mean_W");
    h_results->GetXaxis()->SetBinLabel(5,  "Sigma1_W");
    h_results->GetXaxis()->SetBinLabel(6,  "Sigma2_W");
    h_results->GetXaxis()->SetBinLabel(7,  "Lambda_W");
    h_results->GetXaxis()->SetBinLabel(8,  "Trans_W");
    h_results->GetXaxis()->SetBinLabel(9,  "Chi2");
    h_results->GetXaxis()->SetBinLabel(10, "NDF");
    h_results->Write();

    output->Close();
    f->Close();

    cout << "\n========================================" << endl;
    cout << "PHASE 2 DONE! [" << pt_bin << "]" << endl;
    cout << "Results saved in: " << out_dir << "/phase2_results.root" << endl;
    cout << "========================================\n" << endl;
}
