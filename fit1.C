// fit1.C
// PHASE 1: Background-only fit (signal fraction f fixed to 0).
// Fits the QCD background with a 4th-degree polynomial in the
// mass sidebands, leaving the signal region free.
// Results (pol4 coefficients) are saved to a ROOT file and
// used as fixed input for Phases 2, 3 and 4.
//
// pT bins used in this analysis: pt150_350, pt350_500, pt500_700, pt700_1000
// These can be modified by adding/editing entries in GetConfig_Phase1()
// and adding the corresponding histograms in creating.C
//
// Usage:
//   EOS=/your/output/path
//   root -l -q "fit1.C(\"$EOS/WZ_histograms4bins.root\",\"WZfit_phase1\",\"pt150_350\")"
//   root -l -q "fit1.C(\"$EOS/WZ_histograms4bins.root\",\"WZfit_phase1\",\"pt350_500\")"
//   root -l -q "fit1.C(\"$EOS/WZ_histograms4bins.root\",\"WZfit_phase1\",\"pt500_700\")"
//   root -l -q "fit1.C(\"$EOS/WZ_histograms4bins.root\",\"WZfit_phase1\",\"pt700_1000\")"
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

//  Fit configuration per pT bin 
// Each pT bin has its own initial values and parameter limits for the fit.
// To add a new pT bin: add a new else-if block in GetConfig_Phase1()
// with the corresponding initial values and limits, and make sure the
// matching 2D histogram exists in the input ROOT file (see creating.C).

struct FitConfig_Phase1 {
    Double_t A_init;
    Double_t p1_bkg_init;
    Double_t p2_bkg_init;
    Double_t p3_bkg_init;
    Double_t p4_bkg_init;

    Double_t A_min,    A_max;
    Double_t p1_min,   p1_max;
    Double_t p2_min,   p2_max;
    Double_t p3_min,   p3_max;
    Double_t p4_min,   p4_max;

    Int_t rebin;
    TString pt_label;
};

// Returns fit configuration (initial values and parameter limits) for the
// requested pT bin. Initial values were tuned on the unquenched Pythia sample
// and may need retuning if the binning or dataset changes significantly.

FitConfig_Phase1 GetConfig_Phase1(TString pt_bin) {
    FitConfig_Phase1 cfg;

    if (pt_bin == "all") {
        cfg.A_init      = 1.0;
        cfg.p1_bkg_init = -0.0376651;
        cfg.p2_bkg_init =  0.000543263;
        cfg.p3_bkg_init = -3.547521e-6;
        cfg.p4_bkg_init =  8.82114e-9;

        cfg.A_min  = 0.0;    cfg.A_max  = 10.0;
        cfg.p1_min = -0.5;   cfg.p1_max =  0.5;
        cfg.p2_min = -0.01;  cfg.p2_max =  0.01;
        cfg.p3_min = -1e-4;  cfg.p3_max =  1e-4;
        cfg.p4_min = -1e-6;  cfg.p4_max =  1e-6;

        cfg.rebin    = 2;
        cfg.pt_label = "all p_{T}";
    }
    else if (pt_bin == "pt150_350") {
        cfg.A_init      = 0.05;
        cfg.p1_bkg_init = -0.0376651;
        cfg.p2_bkg_init =  0.000543263;
        cfg.p3_bkg_init = -3.547521e-6;
        cfg.p4_bkg_init =  8.82114e-9;

        cfg.A_min  = 0.001;  cfg.A_max  = 0.5;
        cfg.p1_min = -0.5;   cfg.p1_max =  0.5;
        cfg.p2_min = -0.01;  cfg.p2_max =  0.01;
        cfg.p3_min = -1e-4;  cfg.p3_max =  1e-4;
        cfg.p4_min = -1e-6;  cfg.p4_max =  1e-6;

        cfg.rebin    = 2;
        cfg.pt_label = "150 #leq p_{T} < 350 GeV";
    }
    else if (pt_bin == "pt350_500") {
        cfg.A_init      = 1.737013e-04;
        cfg.p1_bkg_init = -6.705490e-02;
        cfg.p2_bkg_init =  1.567275e-03;
        cfg.p3_bkg_init = -1.441343e-05;
        cfg.p4_bkg_init =  3.597460e-08;

        cfg.A_min  = 0.0001;  cfg.A_max  = 0.01;
        cfg.p1_min = -0.5;    cfg.p1_max =  0.5;
        cfg.p2_min = -0.01;   cfg.p2_max =  0.01;
        cfg.p3_min = -1e-4;   cfg.p3_max =  1e-4;
        cfg.p4_min = -1e-6;   cfg.p4_max =  1e-6;

        cfg.rebin    = 2;
        cfg.pt_label = "350 #leq p_{T} < 500 GeV";
    }
    else if (pt_bin == "pt500_700") {
        cfg.A_init      = 2.2e-05;
        cfg.p1_bkg_init = -0.0517166;
        cfg.p2_bkg_init =  0.00095632;
        cfg.p3_bkg_init = -7.01515e-06;
        cfg.p4_bkg_init =  1.76154e-08;

        cfg.A_min  = 1e-6;    cfg.A_max  = 0.001;
        cfg.p1_min = -0.5;    cfg.p1_max =  0.5;
        cfg.p2_min = -0.01;   cfg.p2_max =  0.01;
        cfg.p3_min = -1e-4;   cfg.p3_max =  1e-4;
        cfg.p4_min = -1e-6;   cfg.p4_max =  1e-6;

        cfg.rebin    = 2;
        cfg.pt_label = "500 #leq p_{T} < 700 GeV";
    }
    else if (pt_bin == "pt700_1000") {
        cfg.A_init      = 5.0e-06;
        cfg.p1_bkg_init = -0.0517166;
        cfg.p2_bkg_init =  0.00095632;
        cfg.p3_bkg_init = -7.01515e-06;
        cfg.p4_bkg_init =  1.76154e-08;

        cfg.A_min  = 1e-8;    cfg.A_max  = 0.0001;
        cfg.p1_min = -0.5;    cfg.p1_max =  0.5;
        cfg.p2_min = -0.01;   cfg.p2_max =  0.01;
        cfg.p3_min = -1e-4;   cfg.p3_max =  1e-4;
        cfg.p4_min = -1e-6;   cfg.p4_max =  1e-6;

        cfg.rebin    = 2;
        cfg.pt_label = "700 #leq p_{T} < 1000 GeV";
    }
    else {
     cout << "ERROR: unrecognised pt_bin: " << pt_bin << endl;
    }
    return cfg;
}

// BACKGROUND: pol. 4 grade
//   bkg(x) = (1 + p[0]*x + p[1]*x^2 + p[2]*x^3 + p[3]*x^4)

Double_t background_func(Double_t *x, Double_t *p) {
    Double_t m = x[0];
    return (1.0 + p[0]*m + p[1]*m*m + p[2]*m*m*m + p[3]*m*m*m*m);
}

// W SIGNAL: Crystal Ball con left tail 
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
//   p[0]=fgauss_Z, p[1]=mean_Z, p[2]=sigma1_Z, p[3]=sigma2_Z,
//   p[4]=lambda_Z, p[5]=trans_Z

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

    Double_t signal = p[2] * W + (1.0 - p[2]) * Z;
    return p[0] * (p[1] * signal + (1.0 - p[1]) * bkg);
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

void fit1(
    TString input_file = EOS_BASE "/WZ_histograms4bins.root",
    TString output_tag = "WZfit_phase1",
    TString pt_bin     = "all",
    Double_t fit_min   = 65.0,
    Double_t fit_max   = 100.0
) {
   
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TString full_tag = (pt_bin == "all") ? output_tag
                                         : Form("%s_%s", output_tag.Data(), pt_bin.Data());
    TString out_dir = Form(EOS_BASE "/figuresB/%s", full_tag.Data());
    gSystem->Exec(Form("mkdir -p %s", out_dir.Data()));

    FitConfig_Phase1 cfg = GetConfig_Phase1(pt_bin);

    TFile* f = TFile::Open(input_file);
    if (!f || f->IsZombie()) {
        cout << "ERROR: Cannot open " << input_file << endl;
        return;
    }

    TString hname = (pt_bin == "all") ? "mFat_n2_2Dbkg"
                                      : Form("mFat_n2_2Dbkg_%s", pt_bin.Data());
    TH2F* h2_bkg = (TH2F*)f->Get(hname);
    if (!h2_bkg) {
        cout << "ERROR: Histogram not found: " << hname << endl;
        f->ls(); return;
    }
    cout << "2D histogram loaded: " << hname << endl;

    Double_t N2_cut  = 0.18;
    int binN2_min    = h2_bkg->GetYaxis()->FindBin(0.0);
    int binN2_max    = h2_bkg->GetYaxis()->FindBin(N2_cut);
    TH1D* hist_bkg   = h2_bkg->ProjectionX("hist_bkg", binN2_min, binN2_max);

    if (cfg.rebin > 1) hist_bkg->Rebin(cfg.rebin);
    Double_t binwidth = hist_bkg->GetBinWidth(1);
    hist_bkg->SetYTitle(Form("Events / %.1f GeV", binwidth));

   

    // FIT function: 19 parameters
    // p[0]       = A
    // p[1]       = f           (fisso = 0)
    // p[2]       = g           (fisso = 0)
    // p[3..8]    = W signal    (fissi, dummy)
    // p[9..14]   = Z signal    (fissi, dummy)
    // p[15..18]  = bkg pol4    (p1, p2, p3, p4) -- to fit

    TF1* fit_phase1 = new TF1("fit_phase1", total_func, fit_min, fit_max, 19);

    fit_phase1->SetParName(0,  "A");
    fit_phase1->SetParName(1,  "f");
    fit_phase1->SetParName(2,  "g");
    fit_phase1->SetParName(3,  "fgauss_W");
    fit_phase1->SetParName(4,  "mean_W");
    fit_phase1->SetParName(5,  "sigma1_W");
    fit_phase1->SetParName(6,  "sigma2_W");
    fit_phase1->SetParName(7,  "lambda_W");
    fit_phase1->SetParName(8,  "trans_W");
    fit_phase1->SetParName(9,  "fgauss_Z");
    fit_phase1->SetParName(10, "mean_Z");
    fit_phase1->SetParName(11, "sigma1_Z");
    fit_phase1->SetParName(12, "sigma2_Z");
    fit_phase1->SetParName(13, "lambda_Z");
    fit_phase1->SetParName(14, "trans_Z");
    fit_phase1->SetParName(15, "p1_bkg");
    fit_phase1->SetParName(16, "p2_bkg");
    fit_phase1->SetParName(17, "p3_bkg");
    fit_phase1->SetParName(18, "p4_bkg");

    // Amplitude
    fit_phase1->SetParameter(0, cfg.A_init);
    fit_phase1->SetParLimits(0, cfg.A_min, cfg.A_max);

    // f=0 e g=0: only bkg in phase 1
    fit_phase1->FixParameter(1, 0.0);
    fit_phase1->FixParameter(2, 0.0);

    // W dummy (fixed)
    fit_phase1->FixParameter(3,  0.5);
    fit_phase1->FixParameter(4,  80.0);
    fit_phase1->FixParameter(5,  2.0);
    fit_phase1->FixParameter(6,  5.0);
    fit_phase1->FixParameter(7,  0.4);   
    fit_phase1->FixParameter(8,  74.0); 

    //  Z dummy (fixed)
    fit_phase1->FixParameter(9,  0.5);
    fit_phase1->FixParameter(10, 91.0);
    fit_phase1->FixParameter(11, 2.0);
    fit_phase1->FixParameter(12, 5.0);
    fit_phase1->FixParameter(13, 0.4);
    fit_phase1->FixParameter(14, 85.0);

    // coeff bkg pol (to fit)
    fit_phase1->SetParameter(15, cfg.p1_bkg_init);
    fit_phase1->SetParameter(16, cfg.p2_bkg_init);
    fit_phase1->SetParameter(17, cfg.p3_bkg_init);
    fit_phase1->SetParameter(18, cfg.p4_bkg_init);

    fit_phase1->SetParLimits(15, cfg.p1_min, cfg.p1_max);
    fit_phase1->SetParLimits(16, cfg.p2_min, cfg.p2_max);
    fit_phase1->SetParLimits(17, cfg.p3_min, cfg.p3_max);
    fit_phase1->SetParLimits(18, cfg.p4_min, cfg.p4_max);

    fit_phase1->SetNpx(1000);

    cout << "\n RUNNING FIT " << endl;
    TFitResultPtr result = hist_bkg->Fit(fit_phase1, "RWLS", "", fit_min, fit_max);

    Double_t chi2    = fit_phase1->GetChisquare();
    Int_t    ndf     = fit_phase1->GetNDF();
    Double_t chi2ndf = (ndf > 0) ? chi2/ndf : -1;

    Double_t A       = fit_phase1->GetParameter(0);
    Double_t err_A   = fit_phase1->GetParError(0);
    Double_t p1_bkg  = fit_phase1->GetParameter(15);
    Double_t ep1     = fit_phase1->GetParError(15);
    Double_t p2_bkg  = fit_phase1->GetParameter(16);
    Double_t ep2     = fit_phase1->GetParError(16);
    Double_t p3_bkg  = fit_phase1->GetParameter(17);
    Double_t ep3     = fit_phase1->GetParError(17);
    Double_t p4_bkg  = fit_phase1->GetParameter(18);
    Double_t ep4     = fit_phase1->GetParError(18);

    cout << "\n========================================" << endl;
    cout << "PHASE 1 RESULTS [" << pt_bin << "]" << endl;
    cout << "========================================" << endl;
    cout << "Chi2/ndf = " << chi2 << " / " << ndf << " = " << chi2ndf << endl;
    cout << "A  = " << A      << " +/- " << err_A << endl;
    cout << "p1 = " << p1_bkg << " +/- " << ep1   << " GeV^{-1}" << endl;
    cout << "p2 = " << p2_bkg << " +/- " << ep2   << " GeV^{-2}" << endl;
    cout << "p3 = " << p3_bkg << " +/- " << ep3   << " GeV^{-3}" << endl;
    cout << "p4 = " << p4_bkg << " +/- " << ep4   << " GeV^{-4}" << endl;
    cout << "========================================\n" << endl;

    // Chi2 
    Double_t chi2_manual = 0.0;
    Int_t    npoints     = 0;
    for (int i = 1; i <= hist_bkg->GetNbinsX(); i++) {
        Double_t x   = hist_bkg->GetBinCenter(i);
        if (x < fit_min || x > fit_max) continue;
        Double_t obs = hist_bkg->GetBinContent(i);
        Double_t err = hist_bkg->GetBinError(i);
        Double_t exp = fit_phase1->Eval(x);
        if (err > 0) {
            chi2_manual += (obs - exp)*(obs - exp) / (err*err);
            npoints++;
        }
    }
    Int_t    ndf_manual     = npoints - fit_phase1->GetNumberFreeParameters();
    Double_t chi2ndf_manual = (ndf_manual > 0) ? chi2_manual / ndf_manual : -1;
    cout << "Manual chi2 = " << chi2_manual
         << " / " << ndf_manual
         << " = "  << chi2ndf_manual << endl;

    // PLOT
    TCanvas* c = new TCanvas("c", "Phase 1 - Background Fit", 900, 800);
    c->SetTicks();
    c->SetLeftMargin(0.14); c->SetRightMargin(0.05);
    c->SetTopMargin(0.08);  c->SetBottomMargin(0.12);

    hist_bkg->SetMarkerStyle(20); hist_bkg->SetMarkerSize(0.9);
    hist_bkg->SetMarkerColor(kBlack); hist_bkg->SetLineColor(kBlack);
    hist_bkg->SetXTitle("Inv. Mass [GeV]");
    hist_bkg->GetXaxis()->SetRangeUser(65, 100);
    hist_bkg->GetXaxis()->SetTitleSize(0.045); hist_bkg->GetXaxis()->SetLabelSize(0.04);
    hist_bkg->GetYaxis()->SetTitleSize(0.045); hist_bkg->GetYaxis()->SetLabelSize(0.04);
    hist_bkg->SetMinimum(0);
    hist_bkg->SetMaximum(hist_bkg->GetMaximum() * 1.4);
    hist_bkg->Draw("PE");

    fit_phase1->SetLineColor(kGreen+2); fit_phase1->SetLineWidth(3);
    fit_phase1->Draw("SAME");

    TLegend* leg1 = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg1->SetBorderSize(0); leg1->SetFillStyle(0);
    leg1->SetTextSize(0.032); leg1->SetTextFont(42);
    leg1->AddEntry(hist_bkg,   "Background only", "PE");
    leg1->AddEntry(fit_phase1, "Pol4 Fit",        "L");
    leg1->Draw();

    TLegend* leg2 = new TLegend(0.52, 0.52, 0.92, 0.73);
    leg2->SetBorderSize(0); leg2->SetFillStyle(0);
    leg2->SetTextSize(0.028); leg2->SetTextFont(42); leg2->SetTextAlign(12);
    leg2->AddEntry((TObject*)0, Form("#chi^{2}/ndf = %.1f / %d = %.2f", chi2_manual, ndf_manual, chi2ndf_manual), "");
    leg2->AddEntry((TObject*)0, Form("A  = %.3e #pm %.3e", A, err_A), "");
    leg2->AddEntry((TObject*)0, "f_{signal} = 0.0 (fixed)", "");
    leg2->AddEntry((TObject*)0, Form("p_{1} = %.4e #pm %.4e", p1_bkg, ep1), "");
    leg2->AddEntry((TObject*)0, Form("p_{2} = %.4e #pm %.4e", p2_bkg, ep2), "");
    leg2->AddEntry((TObject*)0, Form("p_{3} = %.4e #pm %.4e", p3_bkg, ep3), "");
    leg2->AddEntry((TObject*)0, Form("p_{4} = %.4e #pm %.4e", p4_bkg, ep4), "");
    leg2->Draw();

    CMS_latex();

    TLatex* info = new TLatex();
    info->SetNDC(); info->SetTextSize(0.030); info->SetTextFont(42);
    info->DrawLatex(0.17, 0.87, "|#eta^{jet}| < 1.6");
    info->DrawLatex(0.17, 0.83, "R = 0.8");
    info->DrawLatex(0.17, 0.79, "N_{2} < 0.18");
    info->DrawLatex(0.17, 0.75, cfg.pt_label);

    c->SaveAs(Form("%s/phase1_background_fit.pdf", out_dir.Data()));
    cout << "Plot saved: " << out_dir << "/phase1_background_fit.pdf" << endl;

    //Saving results
    // Bin Map -> parameter (not varying, next phases read bin 1-4):
    // bin1=p1_bkg, bin2=p2_bkg, bin3=p3_bkg, bin4=p4_bkg, bin5=A, bin6=f_signal(0)

    TFile* output = new TFile(Form("%s/phase1_results.root", out_dir.Data()), "RECREATE");
    hist_bkg->Write("hist_bkg");
    fit_phase1->Write("fit_phase1");

    TH1F* h_results = new TH1F("results", "Phase 1 Results", 10, 0, 10);
    h_results->SetBinContent(1, p1_bkg); h_results->SetBinError(1, ep1);
    h_results->SetBinContent(2, p2_bkg); h_results->SetBinError(2, ep2);
    h_results->SetBinContent(3, p3_bkg); h_results->SetBinError(3, ep3);
    h_results->SetBinContent(4, p4_bkg); h_results->SetBinError(4, ep4);
    h_results->SetBinContent(5, A);      h_results->SetBinError(5, err_A);
    h_results->SetBinContent(6, 0.0);
    h_results->SetBinContent(8, chi2_manual);
    h_results->SetBinContent(9, ndf_manual);

    h_results->GetXaxis()->SetBinLabel(1, "p1_bkg");
    h_results->GetXaxis()->SetBinLabel(2, "p2_bkg");
    h_results->GetXaxis()->SetBinLabel(3, "p3_bkg");
    h_results->GetXaxis()->SetBinLabel(4, "p4_bkg");
    h_results->GetXaxis()->SetBinLabel(5, "A");
    h_results->GetXaxis()->SetBinLabel(6, "f_signal");
    h_results->GetXaxis()->SetBinLabel(8, "Chi2");
    h_results->GetXaxis()->SetBinLabel(9, "NDF");
    h_results->Write();

    output->Close();
    f->Close();

    cout << "\n========================================" << endl;
    cout << "PHASE 1 DONE [" << pt_bin << "]" << endl;
    cout << "Results saved in: " << out_dir << "/phase1_results.root" << endl;
    cout << "========================================\n" << endl;
}
