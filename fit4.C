// fit4.C
// PHASE 4: Full fit W + Z + background.
// All shape parameters are fixed from Phases 1, 2 and 3.
// Only the global amplitude A and the signal fraction f are free.
// The W fraction within the signal (g) is computed from the
// amplitudes fitted in Phases 2 and 3 as g = A_W/(A_W+A_Z).
// The key physics output is the reconstructed W and Z mass peaks
// as a function of pT bin, used to study time-dependent quenching.
//
// pT bins: pt150_350, pt350_500, pt500_700, pt700_1000
// To add a new bin: add an entry in GetConfig_Phase4() and in creating.C
//
// Usage:
//   EOS=/your/output/path
//   root -l -q "fit4.C(\"$EOS/WZ_histograms4bins.root\",\"auto\",\"auto\",\"auto\",\"WZfit_phase4\",\"pt150_350\")"
//   (repeat for each pT bin; "auto" reconstructs paths from EOS_BASE)
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
#define EOS_BASE "/eos/home-f/fzavatti/eventi6_R0.2_bkgCONsub400pclesSOFTz02_unq"

struct FitConfig_Phase4 {
    Double_t A_init;
    Double_t f_init;
    Double_t A_min, A_max;
    Double_t f_min, f_max;
    Double_t g_fixed;
    Int_t rebin;
    TString pt_label;
};

FitConfig_Phase4 GetConfig_Phase4(TString pt_bin) {
    FitConfig_Phase4 cfg;

    if (pt_bin == "all") {
        cfg.A_init  = 0.3;   cfg.f_init  = 0.040;
        cfg.A_min   = 0;     cfg.A_max   = 1;
        cfg.f_min   = 0.0;   cfg.f_max   = 1.0;
        cfg.g_fixed = -1;
        cfg.rebin    = 1;
        cfg.pt_label = "all p_{T}";
    }
    else if (pt_bin == "pt150_350") {
        cfg.A_init  = 0.4;   cfg.f_init  = 0.040;
        cfg.A_min   = 0.0;   cfg.A_max   = 1;
        cfg.f_min   = 0.0;   cfg.f_max   = 1.0;
        cfg.g_fixed = -1;
        cfg.rebin    = 1;
        cfg.pt_label = "150 #leq p_{T} < 350 GeV";
    }
    else if (pt_bin == "pt350_500") {
        cfg.A_init  = 0.0076; cfg.f_init  = 0.040;
        cfg.A_min   = 0;      cfg.A_max   = 1;
        cfg.f_min   = 0.0;    cfg.f_max   = 1.0;
        cfg.g_fixed = -1;
        cfg.rebin    = 1;
        cfg.pt_label = "350 #leq p_{T} < 500 GeV";
    }
    else if (pt_bin == "pt500_700") {
        cfg.A_init  = 0.0010; cfg.f_init  = 0.040;
        cfg.A_min   = 0;      cfg.A_max   = 1;
        cfg.f_min   = 0.0;    cfg.f_max   = 1.0;
        cfg.g_fixed = -1;
        cfg.rebin    = 1;
        cfg.pt_label = "500 #leq p_{T} < 700 GeV";
    }
    else if (pt_bin == "pt700_1000") {
        cfg.A_init  = 0.0001; cfg.f_init  = 0.040;
        cfg.A_min   = 0;      cfg.A_max   = 1;
        cfg.f_min   = 0.0;    cfg.f_max   = 1.0;
        cfg.g_fixed = -1;
        cfg.rebin    = 1;
        cfg.pt_label = "700 #leq p_{T} < 1000 GeV";
    }
    else {
        cout << "ERROR: unrecognised pt_bin: " << pt_bin << endl;
    }
    return cfg;
}

// BACKGROUND: polinomio 4 grade

Double_t background_func(Double_t *x, Double_t *p) {
    Double_t m = x[0];
    return (1.0 + p[0]*m + p[1]*m*m + p[2]*m*m*m + p[3]*m*m*m*m);
}

// W SIGNAL: Crystal Ball left tail
//   p[0]=fgauss_W, p[1]=mean_W, p[2]=sigma1_W, p[3]=sigma2_W,
//   p[4]=lambda_W, p[5]=trans_W

Double_t W_signal_func(Double_t *x, Double_t *p) {
    Double_t G_trans = p[0] * TMath::Gaus(p[5], p[1], p[2], kFALSE) +
                       (1.0 - p[0]) * TMath::Gaus(p[5], p[1], p[3], kFALSE);
    if (x[0] > p[5]) return p[0] * TMath::Gaus(x[0], p[1], p[2], kFALSE) +
                            (1.0 - p[0]) * TMath::Gaus(x[0], p[1], p[3], kFALSE);
    else              return G_trans * TMath::Exp(p[4] * (x[0] - p[5]));
}

// Z SIGNAL: Crystal Ball left tail

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

TString BuildPhasePath(TString base_tag, int phase_num, TString pt_bin) {
    TString tag = (pt_bin == "all") ? base_tag
                                    : Form("%s_%s", base_tag.Data(), pt_bin.Data());
    return Form(EOS_BASE "/figuresB/%s/phase%d_results.root", tag.Data(), phase_num);
}

void fit4(
    TString input_file  = EOS_BASE "/WZ_histograms4bins.root",
    TString phase1_file = "auto",
    TString phase2_file = "auto",
    TString phase3_file = "auto",
    TString output_tag  = "WZfit_phase4",
    TString pt_bin      = "all",
    Double_t fit_min    = 65.0,
    Double_t fit_max    = 100.0
) {
    

    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TString full_tag = (pt_bin == "all") ? output_tag
                                         : Form("%s_%s", output_tag.Data(), pt_bin.Data());
    TString out_dir = Form(EOS_BASE "/figuresB/%s", full_tag.Data());
    gSystem->Exec(Form("mkdir -p %s", out_dir.Data()));

    if (phase1_file == "auto") phase1_file = BuildPhasePath("WZfit_phase1", 1, pt_bin);
    if (phase2_file == "auto") phase2_file = BuildPhasePath("WZfit_phase2", 2, pt_bin);
    if (phase3_file == "auto") phase3_file = BuildPhasePath("WZfit_phase3", 3, pt_bin);

    cout << "Reading phase 1 from: " << phase1_file << endl;
    cout << "Reading phase 2 from: "<< phase2_file << endl;
    cout << "Reading phase 3 from: " << phase3_file << endl;

    //  Fase 1: parameters pol4 background 
    TFile* f1 = TFile::Open(phase1_file);
    if (!f1 || f1->IsZombie()) { cout << "ERROR: phase 1 file not found!" << endl; return; }
    TH1F* h1 = (TH1F*)f1->Get("results");
    if (!h1) { cout << "ERROR: phase 1 results not found!" << endl; return; }
    Double_t p1_bkg = h1->GetBinContent(1);
    Double_t p2_bkg = h1->GetBinContent(2);
    Double_t p3_bkg = h1->GetBinContent(3);
    Double_t p4_bkg = h1->GetBinContent(4);
    cout << "\nPHASE 1 (Background pol4):" << endl;
    cout << "  p1=" << p1_bkg << "  p2=" << p2_bkg
         << "  p3=" << p3_bkg << "  p4=" << p4_bkg << endl;
    f1->Close();

    //  Phase 2: parameters signal W 
    // bin1=A, bin2=f_signal, bin3=fgauss_W, bin4=Mean_W, bin5=Sigma1_W,
    // bin6=Sigma2_W, bin7=Lambda_W, bin8=Trans_W
    TFile* f2 = TFile::Open(phase2_file);
    if (!f2 || f2->IsZombie()) { cout << "ERROR: phase 2 file not found!" << endl; return; }
    TH1F* h2 = (TH1F*)f2->Get("results");
    if (!h2) { cout << "ERROR: phase 2 results not found!" << endl; return; }
    Double_t A_W        = h2->GetBinContent(1);
    Double_t fgauss_W   = h2->GetBinContent(3);
    Double_t mean_W     = h2->GetBinContent(4);
    Double_t err_mean_W = h2->GetBinError(4);
    Double_t sigma1_W   = h2->GetBinContent(5);
    Double_t sigma2_W   = h2->GetBinContent(6);
    Double_t lambda_W   = h2->GetBinContent(7);
    Double_t trans_W    = h2->GetBinContent(8);
    cout << "\nPHASE 2  (W Crystal Ball):" << endl;
    cout << "  A_W=" << A_W << "  fgauss=" << fgauss_W
         << "  mean=" << mean_W << " +/- " << err_mean_W << endl;
    cout << "  s1=" << sigma1_W << "  s2=" << sigma2_W
         << "  lambda=" << lambda_W << "  trans=" << trans_W << endl;
    f2->Close();

    //  Phase 3: parameters signal Z 
    // bin1=A, bin2=f_signal, bin3=fgauss_Z, bin4=Mean_Z, bin5=Sigma1_Z,
    // bin6=Sigma2_Z, bin7=Lambda_Z, bin8=Trans_Z
    TFile* f3 = TFile::Open(phase3_file);
    if (!f3 || f3->IsZombie()) { cout << "ERROR: phase 3 file not found!" << endl; return; }
    TH1F* h3 = (TH1F*)f3->Get("results");
    if (!h3) { cout << "ERROR: phase 3 results not found!" << endl; return; }
    Double_t A_Z        = h3->GetBinContent(1);
    Double_t fgauss_Z   = h3->GetBinContent(3);
    Double_t mean_Z     = h3->GetBinContent(4);
    Double_t err_mean_Z = h3->GetBinError(4);
    Double_t sigma1_Z   = h3->GetBinContent(5);
    Double_t sigma2_Z   = h3->GetBinContent(6);
    Double_t lambda_Z   = h3->GetBinContent(7);
    Double_t trans_Z    = h3->GetBinContent(8);
    cout << "\nPHASE 3 (Z Crystal Ball):" << endl;
    cout << "  A_Z=" << A_Z << "  fgauss=" << fgauss_Z
         << "  mean=" << mean_Z << " +/- " << err_mean_Z << endl;
    f3->Close();

    // Mean W and Z masses
    Double_t mean_avg     = 0.5 * (mean_W + mean_Z);
    Double_t err_mean_avg = 0.5 * TMath::Sqrt(err_mean_W*err_mean_W + err_mean_Z*err_mean_Z);
    cout << "\nAverage of W and Z masses:" << endl;
    cout << "  mean_avg = " << mean_avg << " +/- " << err_mean_avg << endl;

    FitConfig_Phase4 cfg = GetConfig_Phase4(pt_bin);
    Double_t g_param;
    if (cfg.g_fixed < 0) {
        g_param = (A_W + A_Z > 0) ? A_W / (A_W + A_Z) : 0.78;
        cout << "\ng computed from A_W/(A_W+A_Z) = " << g_param << endl;
    } else {
        g_param = cfg.g_fixed;
        cout << "\ng fixed: " << g_param << endl;
    }

    
    TFile* file_input = TFile::Open(input_file);
    if (!file_input || file_input->IsZombie()) {
        cout << "ERROR: " << input_file << endl; return;
    }

    TString hname = (pt_bin == "all") ? "mFat_n2_2D"
                                      : Form("mFat_n2_2D_%s", pt_bin.Data());
    TH2F* h2_total = (TH2F*)file_input->Get(hname);
    if (!h2_total) { cout << "ERROR: Hist not found: " << hname << endl; return; }

    Double_t N2_cut  = 0.18;
    int binN2_min    = h2_total->GetYaxis()->FindBin(0.0);
    int binN2_max    = h2_total->GetYaxis()->FindBin(N2_cut);
    TH1D* hist_total = h2_total->ProjectionX("hist_total", binN2_min, binN2_max);

 
    if (cfg.rebin > 1) hist_total->Rebin(cfg.rebin);
    Double_t binwidth = hist_total->GetBinWidth(1);
    hist_total->SetYTitle(Form("Events / %.1f GeV", binwidth));

    // fit function
    TF1* fit_phase4 = new TF1("fit_phase4", total_func, fit_min, fit_max, 19);
    fit_phase4->SetParName(0,  "A");
    fit_phase4->SetParName(1,  "f");
    fit_phase4->SetParName(2,  "g");
    fit_phase4->SetParName(3,  "fgauss_W");
    fit_phase4->SetParName(4,  "mean_W");
    fit_phase4->SetParName(5,  "sigma1_W");
    fit_phase4->SetParName(6,  "sigma2_W");
    fit_phase4->SetParName(7,  "lambda_W");
    fit_phase4->SetParName(8,  "trans_W");
    fit_phase4->SetParName(9,  "fgauss_Z");
    fit_phase4->SetParName(10, "mean_Z");
    fit_phase4->SetParName(11, "sigma1_Z");
    fit_phase4->SetParName(12, "sigma2_Z");
    fit_phase4->SetParName(13, "lambda_Z");
    fit_phase4->SetParName(14, "trans_Z");
    fit_phase4->SetParName(15, "p1_bkg");
    fit_phase4->SetParName(16, "p2_bkg");
    fit_phase4->SetParName(17, "p3_bkg");
    fit_phase4->SetParName(18, "p4_bkg");

    // Only A e f free
    fit_phase4->SetParameter(0, cfg.A_init);
    fit_phase4->SetParameter(1, cfg.f_init);
    fit_phase4->SetParLimits(0, cfg.A_min, cfg.A_max);
    fit_phase4->SetParLimits(1, cfg.f_min, cfg.f_max);

    //everything else fixed by previous phases
    fit_phase4->FixParameter(2,  g_param);
    fit_phase4->FixParameter(3,  fgauss_W);
    fit_phase4->FixParameter(4,  mean_W);
    fit_phase4->FixParameter(5,  sigma1_W);
    fit_phase4->FixParameter(6,  sigma2_W);
    fit_phase4->FixParameter(7,  lambda_W);
    fit_phase4->FixParameter(8,  trans_W);
    fit_phase4->FixParameter(9,  fgauss_Z);
    fit_phase4->FixParameter(10, mean_Z);
    fit_phase4->FixParameter(11, sigma1_Z);
    fit_phase4->FixParameter(12, sigma2_Z);
    fit_phase4->FixParameter(13, lambda_Z);
    fit_phase4->FixParameter(14, trans_Z);
    fit_phase4->FixParameter(15, p1_bkg);
    fit_phase4->FixParameter(16, p2_bkg);
    fit_phase4->FixParameter(17, p3_bkg);
    fit_phase4->FixParameter(18, p4_bkg);
    fit_phase4->SetNpx(1200);

    cout << "\n RUNNING COMPLETED FIT " << endl;
    TFitResultPtr result = hist_total->Fit(fit_phase4, "RWLS", "", fit_min, fit_max);

    Double_t A_final  = fit_phase4->GetParameter(0);
    Double_t err_A    = fit_phase4->GetParError(0);
    Double_t f_final  = fit_phase4->GetParameter(1);
    Double_t err_f    = fit_phase4->GetParError(1);

    // Chi2 
    Double_t chi2_manual = 0.0;
    Int_t    npoints     = 0;
    for (int i = 1; i <= hist_total->GetNbinsX(); i++) {
        Double_t x   = hist_total->GetBinCenter(i);
        if (x < fit_min || x > fit_max) continue;
        Double_t obs = hist_total->GetBinContent(i);
        Double_t err = hist_total->GetBinError(i);
        Double_t exp = fit_phase4->Eval(x);
        if (err > 0) {
            chi2_manual += (obs - exp)*(obs - exp) / (err*err);
            npoints++;
        }
    }
    Int_t    ndf_manual     = npoints - fit_phase4->GetNumberFreeParameters();
    Double_t chi2ndf_manual = (ndf_manual > 0) ? chi2_manual / ndf_manual : -1;

    cout << "\n========================================" << endl;
    cout << "RESULTS PHASE 4 [" << pt_bin << "]" << endl;
    cout << "========================================" << endl;
    cout << "Chi2 / ndf = " << chi2_manual << " / " << ndf_manual
         << " = " << chi2ndf_manual << endl;
    cout << "A       = " << A_final << " +/- " << err_A   << endl;
    cout << "f       = " << f_final << " +/- " << err_f   << endl;
    cout << "g       = " << g_param << " (fixed)"         << endl;
    cout << "mean_W  = " << mean_W  << " +/- " << err_mean_W  << " (fixed from phase2)" << endl;
    cout << "mean_Z  = " << mean_Z  << " +/- " << err_mean_Z  << " (fixed from phase3)" << endl;
    cout << "mean_avg= " << mean_avg << " +/- " << err_mean_avg << " (simple average)"  << endl;
    cout << "========================================\n" << endl;

    // plot separated components
    TF1* W_component = new TF1("W_component",
        [A_final, f_final, g_param, fgauss_W, mean_W, sigma1_W, sigma2_W, lambda_W, trans_W]
        (Double_t *x, Double_t*) {
            Double_t pars[6] = {fgauss_W, mean_W, sigma1_W, sigma2_W, lambda_W, trans_W};
            return A_final * f_final * g_param * W_signal_func(x, pars);
        }, fit_min, fit_max, 0);

    TF1* Z_component = new TF1("Z_component",
        [A_final, f_final, g_param, fgauss_Z, mean_Z, sigma1_Z, sigma2_Z, lambda_Z, trans_Z]
        (Double_t *x, Double_t*) {
            Double_t pars[6] = {fgauss_Z, mean_Z, sigma1_Z, sigma2_Z, lambda_Z, trans_Z};
            return A_final * f_final * (1.0 - g_param) * Z_signal_func(x, pars);
        }, fit_min, fit_max, 0);

    TF1* bkg_component = new TF1("bkg_component",
        [A_final, f_final, p1_bkg, p2_bkg, p3_bkg, p4_bkg](Double_t *x, Double_t*) {
            Double_t pars[4] = {p1_bkg, p2_bkg, p3_bkg, p4_bkg};
            return A_final * (1.0 - f_final) * background_func(x, pars);
        }, fit_min, fit_max, 0);

    W_component->SetLineColor(kRed+1);      W_component->SetLineWidth(2);   W_component->SetLineStyle(2);
    Z_component->SetLineColor(kBlue+1);     Z_component->SetLineWidth(2);   Z_component->SetLineStyle(2);
    bkg_component->SetLineColor(kOrange+1); bkg_component->SetLineWidth(2); bkg_component->SetLineStyle(2);

    TCanvas* c = new TCanvas("c", "Phase 4 - Complete Fit", 900, 800);
    c->SetTicks();
    c->SetLeftMargin(0.14); c->SetRightMargin(0.05);
    c->SetTopMargin(0.08);  c->SetBottomMargin(0.12);

    hist_total->SetMarkerStyle(20); hist_total->SetMarkerSize(0.9);
    hist_total->SetMarkerColor(kBlack); hist_total->SetLineColor(kBlack);
    hist_total->SetXTitle("Inv. Mass [GeV]");
    hist_total->GetXaxis()->SetRangeUser(65, 100);
    hist_total->GetXaxis()->SetTitleSize(0.045); hist_total->GetXaxis()->SetLabelSize(0.04);
    hist_total->GetYaxis()->SetTitleSize(0.045); hist_total->GetYaxis()->SetLabelSize(0.04);
    hist_total->SetMinimum(0);
    hist_total->SetMaximum(hist_total->GetMaximum() * 2.0);
    hist_total->Draw("PE");

    fit_phase4->SetLineColor(kGreen+2); fit_phase4->SetLineWidth(3);
    fit_phase4->Draw("SAME");
    W_component->Draw("SAME");
    Z_component->Draw("SAME");
    bkg_component->Draw("SAME");

    TLegend* leg1 = new TLegend(0.50, 0.72, 0.88, 0.90);
    leg1->SetBorderSize(0); leg1->SetFillStyle(0);
    leg1->SetTextSize(0.026); leg1->SetTextFont(42);
    leg1->AddEntry(hist_total,    "Data (W + Z + Bkg)",                        "PE");
    leg1->AddEntry(fit_phase4,    "Complete Fit",                               "L");
    leg1->AddEntry(W_component,   Form("W Component (g=%.2f)",   g_param),     "L");
    leg1->AddEntry(Z_component,   Form("Z Component (1-g=%.2f)", 1-g_param),   "L");
    leg1->AddEntry(bkg_component, "Background Component (Pol4)",                "L");
    leg1->Draw();

    TLegend* leg2 = new TLegend(0.54, 0.51, 0.88, 0.71);
    leg2->SetBorderSize(0); leg2->SetFillStyle(0);
    leg2->SetTextSize(0.026); leg2->SetTextFont(42); leg2->SetTextAlign(12);
    leg2->AddEntry((TObject*)0, Form("#chi^{2}/ndf = %.1f / %d = %.2f", chi2_manual, ndf_manual, chi2ndf_manual), "");
    leg2->AddEntry((TObject*)0, Form("A = %.2e #pm %.2e", A_final, err_A), "");
    leg2->AddEntry((TObject*)0, Form("f = %.3f #pm %.3f", f_final, err_f), "");
    leg2->AddEntry((TObject*)0, Form("g = %.3f (fixed)", g_param), "");
    leg2->AddEntry((TObject*)0, Form("M_{W} = %.2f #pm %.2f GeV", mean_W, err_mean_W), "");
    leg2->AddEntry((TObject*)0, Form("M_{Z} = %.2f #pm %.2f GeV", mean_Z, err_mean_Z), "");
    leg2->AddEntry((TObject*)0, Form("#bar{M} = %.2f #pm %.2f GeV", mean_avg, err_mean_avg), "");
    leg2->Draw();

    CMS_latex();

    TLatex* info = new TLatex();
    info->SetNDC(); info->SetTextSize(0.030); info->SetTextFont(42);
    info->DrawLatex(0.17, 0.87, "|#eta^{jet}| < 1.6");
    info->DrawLatex(0.17, 0.83, "R = 0.8");
    info->DrawLatex(0.17, 0.79, "N_{2} < 0.18");
    info->DrawLatex(0.17, 0.75, cfg.pt_label);

    c->SaveAs(Form("%s/phase4_complete_fit.pdf", out_dir.Data()));
    c->SaveAs(Form("%s/phase4_complete_fit.root", out_dir.Data()));
    cout << "Plot saved: " << out_dir << "/phase4_complete_fit.pdf" << endl;

    TFile* output = new TFile(Form("%s/phase4_results.root", out_dir.Data()), "RECREATE");
    hist_total->Write("hist_total");
    fit_phase4->Write("fit_phase4");
    W_component->Write("W_component");
    Z_component->Write("Z_component");
    bkg_component->Write("bkg_component");

    TH1F* h_results = new TH1F("results", "Phase 4 Results", 8, 0, 8);
    h_results->SetBinContent(1, A_final);     h_results->SetBinError(1, err_A);
    h_results->SetBinContent(2, f_final);     h_results->SetBinError(2, err_f);
    h_results->SetBinContent(3, g_param);
    h_results->SetBinContent(4, chi2_manual);
    h_results->SetBinContent(5, ndf_manual);
    h_results->GetXaxis()->SetBinLabel(1, "A");
    h_results->GetXaxis()->SetBinLabel(2, "f");
    h_results->GetXaxis()->SetBinLabel(3, "g");
    h_results->GetXaxis()->SetBinLabel(4, "Chi2");
    h_results->GetXaxis()->SetBinLabel(5, "NDF");
    h_results->Write();

    output->Close();
    file_input->Close();

    cout << "\n========================================" << endl;
    cout << "PHASE 4 COMPLETED! [" << pt_bin << "]" << endl;
    cout << "Results in: " << out_dir << "/phase4_results.root" << endl;
    cout << "========================================\n" << endl;
}
