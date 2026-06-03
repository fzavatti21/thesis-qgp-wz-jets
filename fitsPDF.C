// fitsPDF.C
// Toy MC study to evaluate the statistical uncertainty on the
// reconstructed W and Z mass peaks as a function of pT bin.
// Generates pseudo-experiments by sampling from the fitted total
// function (W + Z + background), fits each toy with the same
// procedure used on data, and collects the distribution of
// fitted mean_W, mean_Z and their average.
//
// All parameters (true_mean_W, true_mean_Z, bkg_p1..p4, etc.)
// must be replaced with the results of the four-phase fit
// (fit1.C -> fit4.C) for your specific sample and scenario.
// nEventsList must be set to the number of events observed in
// each pT bin after the N2 < 0.18 cut.
//
// pT bins: currently 3 (pt150_250, pt250_350, pt350_500).
// To add a new bin: add one element to each array and
// extend the loop to iTest < 4 (or however many bins you have).
//
// Usage:
//   root -l fitsPDF.C
//
// NOTE: update the output path in gSystem->Exec and c->SaveAs
//       to match your directory.

// Background: polinomio di 4° grado
Double_t background_func(Double_t *x, Double_t *p) {
    return (1.0 + p[0]*x[0] + p[1]*x[0]*x[0] + p[2]*x[0]*x[0]*x[0] + p[3]*x[0]*x[0]*x[0]*x[0]);
}

// W SIGNAL: Crystal Ball con coda sinistra
// p[0]=fgauss_W, p[1]=mean_W, p[2]=sigma1_W, p[3]=sigma2_W, p[4]=lambda_W, p[5]=trans_W
Double_t W_signal_func(Double_t *x, Double_t *p) {
    Double_t G_trans = p[0] * TMath::Gaus(p[5], p[1], p[2], kTRUE) +
                       (1.0 - p[0]) * TMath::Gaus(p[5], p[1], p[3], kTRUE);
    if (x[0] > p[5]) {
        return p[0] * TMath::Gaus(x[0], p[1], p[2], kTRUE) +
               (1.0 - p[0]) * TMath::Gaus(x[0], p[1], p[3], kTRUE);
    } else {
        return G_trans * TMath::Exp(p[4] * (x[0] - p[5]));
    }
}

Double_t Z_signal_func(Double_t *x, Double_t *p) {
    Double_t G_trans = p[0] * TMath::Gaus(p[5], p[1], p[2], kTRUE) +
                       (1.0 - p[0]) * TMath::Gaus(p[5], p[1], p[3], kTRUE);
    if (x[0] > p[5]) {
        return p[0] * TMath::Gaus(x[0], p[1], p[2], kTRUE) +
               (1.0 - p[0]) * TMath::Gaus(x[0], p[1], p[3], kTRUE);
    } else {
        return G_trans * TMath::Exp(p[4] * (x[0] - p[5]));
    }
}

// p[0]=A, p[1]=f, p[2]=g
// p[3]=fgauss_W, p[4]=mean_W, p[5]=sigma1_W, p[6]=sigma2_W, p[7]=lambda_W, p[8]=trans_W
// p[9]=fgauss_Z, p[10]=mean_Z, p[11]=sigma1_Z, p[12]=sigma2_Z, p[13]=lambda_Z, p[14]=trans_Z
// p[15]=p1_bkg, p[16]=p2_bkg, p[17]=p3_bkg, p[18]=p4_bkg
Double_t total_func(Double_t *x, Double_t *p) {
    Double_t W_pars[6]   = {p[3],  p[4],  p[5],  p[6],  p[7],  p[8]};
    Double_t Z_pars[6]   = {p[9],  p[10], p[11], p[12], p[13], p[14]};
    Double_t bkg_pars[4] = {p[15], p[16], p[17], p[18]};
    Double_t W   = W_signal_func(x, W_pars);
    Double_t Z   = Z_signal_func(x, Z_pars);
    Double_t bkg = background_func(x, bkg_pars);
    return p[0] * (p[1] * (p[2]*W + (1.0-p[2])*Z) + (1.0-p[1])*bkg);
}


void fitsPDF() {
    gROOT->SetBatch(kTRUE);
    
    // NOTE: change this path to match your output directory
    gSystem->Exec("mkdir -p /eos/home-f/fzavatti/eventi6_R0.2_unq6/figPDFB");
   
   //  Input parameters from the four-phase fit 
// Replace all values below with the results from fit1.C-fit4.C
// for your specific sample (quenching scenario and collision system).
// nEventsList: number of events in each pT bin after N2 < 0.18 cut.
Int_t nEventsList[4] = {182229, 32198, 7441, 0};

TString plotTitles[4] = {
    "Jet p_{T}: 150-250 GeV",
    "Jet p_{T}: 250-350 GeV",
    "Jet p_{T}: 350-500 GeV",
    "Jet p_{T}: 500-700 GeV"   
};

Double_t A_init[4]        = {0.0237312,   0.000558138,  8.79538e-05,  0.0};
Double_t f_init[4]        = {0.0100841,   0.305174,     0.42358,      0.0};
Double_t true_g[4]        = {0.809137,    0.777015,     0.763333,     0.0};

// W Crystal Ball parameters from Phase 2 fit. Replace with your fitted values.
Double_t true_fgauss_W[4] = {0.254672,    0.774626,     0.737615,     0.0};
Double_t true_mean_W[4]   = {80.2995,     80.5468,      80.642,       0.0};
Double_t true_sigma1_W[4] = {3.17504,     0.996092,     0.882094,     0.0};
Double_t true_sigma2_W[4] = {1.1127,      3.52214,      3.46784,      0.0};
Double_t true_lambda_W[4] = {0.254968,    0.4,          0.498973,     0.0};
Double_t true_trans_W[4]  = {76.9153,     69.0,         78.363,       0.0};

// Z Crystal Ball parameters from Phase 3 fit. Replace with your fitted values.
Double_t true_fgauss_Z[4] = {0.381951,    0.229219,     0.193266,     0.0};
Double_t true_mean_Z[4]   = {91.0153,     91.278,       91.4182,      0.0};
Double_t true_sigma1_Z[4] = {2.79684,     3.82164,      4.52343,      0.0};
Double_t true_sigma2_Z[4] = {1.22219,     1.19811,      1.26534,      0.0};
Double_t true_lambda_Z[4] = {0.140374,    0.206877,     0.394484,     0.0};
Double_t true_trans_Z[4]  = {87.5,        83.962,       87.9399,      0.0};

// Background pol4 parameters from Phase 1 fit. Replace with your fitted values.
Double_t bkg_p1[4]        = {-0.0374666,  -0.0670179,   -0.0525983,   0.0};
Double_t bkg_p2[4]        = {0.000541508,  0.00156434,   0.000977926,  0.0};
Double_t bkg_p3[4]        = {-3.56788e-06, -1.42413e-05, -6.97326e-06, 0.0};
Double_t bkg_p4[4]        = {9.0035e-09,   4.46927e-08,  1.68456e-08,  0.0};

    // DEFINE TOYS NUMBER
    Int_t nToys = 1000; 

    for (Int_t iTest = 0; iTest < 4; iTest++) {
        Int_t nEvents = nEventsList[iTest];
        TString title = plotTitles[iTest];

        cout << "STARTING TOY MC FOR BIN " << iTest+1 << ": " << title << endl;

        // Istogrammi di controllo per salvare i risultati di OGNI fit
        TH1F *h_results_mean_W = new TH1F(Form("h_res_mean_W_%d", iTest), "Distribution of fitted mean_W", 100, true_mean_W[iTest]-5, true_mean_W[iTest]+5);
        TH1F *h_results_mean_Z = new TH1F(Form("h_res_mean_Z_%d", iTest), "Distribution of mean_Z", 100, true_mean_Z[iTest]-5, true_mean_Z[iTest]+5);
        
        // histogram of the mean
        TH1F *h_results_M_mean = new TH1F(Form("h_res_M_mean_%d", iTest), "Average fitted (M_W + M_Z)/2 ", 100, 78, 98);

        // generating function
        TF1 *f_gen = new TF1(Form("f_gen_%d", iTest), total_func, 65, 100, 19);
        f_gen->SetParameter(0,  1.0);
        f_gen->SetParameter(1,  f_init[iTest]);
        f_gen->SetParameter(2,  true_g[iTest]);
        f_gen->SetParameter(3,  true_fgauss_W[iTest]);
        f_gen->SetParameter(4,  true_mean_W[iTest]);
        f_gen->SetParameter(5,  true_sigma1_W[iTest]);
        f_gen->SetParameter(6,  true_sigma2_W[iTest]);
        f_gen->SetParameter(7,  true_lambda_W[iTest]);
        f_gen->SetParameter(8,  true_trans_W[iTest]);
        f_gen->SetParameter(9,  true_fgauss_Z[iTest]);
        f_gen->SetParameter(10, true_mean_Z[iTest]);
        f_gen->SetParameter(11, true_sigma1_Z[iTest]);
        f_gen->SetParameter(12, true_sigma2_Z[iTest]);
        f_gen->SetParameter(13, true_lambda_Z[iTest]);
        f_gen->SetParameter(14, true_trans_Z[iTest]);
        f_gen->SetParameter(15, bkg_p1[iTest]);
        f_gen->SetParameter(16, bkg_p2[iTest]);
        f_gen->SetParameter(17, bkg_p3[iTest]);
        f_gen->SetParameter(18, bkg_p4[iTest]);

        // LOOP TOY MC 
        for (Int_t iToy = 0; iToy < nToys; iToy++) {

            gRandom->SetSeed(33 + iToy); 

            TH1F *h = new TH1F("h_toy", title, 70, 65, 100);
            
            //Poisson fluctuation on the number of total events
            Int_t nEventsFluctuated = gRandom->Poisson(nEvents);
            for (Int_t i = 0; i < nEventsFluctuated; i++) {
                h->Fill(f_gen->GetRandom());
            }
            h->Rebin(3);

            // Setup fit function
            TF1 *f_fit = new TF1("f_fit_toy", total_func, 65, 100, 19);
            f_fit->SetNpx(1000);
            f_fit->SetParameter(0, A_init[iTest]);
            f_fit->SetParameter(1, f_init[iTest]);
            f_fit->FixParameter(2,  true_g[iTest]);
            f_fit->FixParameter(3,  true_fgauss_W[iTest]);
            f_fit->FixParameter(4,  true_mean_W[iTest]);
            f_fit->FixParameter(5,  true_sigma1_W[iTest]);
            f_fit->FixParameter(6,  true_sigma2_W[iTest]);
            f_fit->FixParameter(7,  true_lambda_W[iTest]);
            f_fit->FixParameter(8,  true_trans_W[iTest]);
            f_fit->FixParameter(9,  true_fgauss_Z[iTest]);
            f_fit->FixParameter(10, true_mean_Z[iTest]);
            f_fit->FixParameter(11, true_sigma1_Z[iTest]);
            f_fit->FixParameter(12, true_sigma2_Z[iTest]);
            f_fit->FixParameter(13, true_lambda_Z[iTest]);
            f_fit->FixParameter(14, true_trans_Z[iTest]);
            f_fit->FixParameter(15, bkg_p1[iTest]);
            f_fit->FixParameter(16, bkg_p2[iTest]);
            f_fit->FixParameter(17, bkg_p3[iTest]);
            f_fit->FixParameter(18, bkg_p4[iTest]);

            f_fit->FixParameter(1, f_init[iTest]);
             h->Fit(f_fit, "RQWLN");

            // Phase 2: release parameters
            f_fit->ReleaseParameter(1);  f_fit->SetParLimits(1, 0.0, 1.0);
            f_fit->ReleaseParameter(4);  f_fit->SetParLimits(4, true_mean_W[iTest]-3, true_mean_W[iTest]+3);
            f_fit->ReleaseParameter(10); f_fit->SetParLimits(10, true_mean_Z[iTest]-3, true_mean_Z[iTest]+3);
             h->Fit(f_fit, "RQWLN");

            // Phase 3: release sigma
            if (iTest == 0) {
                f_fit->ReleaseParameter(5); f_fit->ReleaseParameter(6);
                f_fit->ReleaseParameter(11); f_fit->ReleaseParameter(12);
                f_fit->SetParLimits(5,0.5,6.0); f_fit->SetParLimits(6,0.5,6.0);
                f_fit->SetParLimits(11,0.5,6.0); f_fit->SetParLimits(12,0.5,6.0);
            } else {
                f_fit->ReleaseParameter(5); f_fit->ReleaseParameter(6);
                f_fit->ReleaseParameter(11); f_fit->ReleaseParameter(12);
                f_fit->SetParLimits(5,0.5,6.0); f_fit->SetParLimits(6,0.5,6.0);
                f_fit->SetParLimits(11,0.5,6); f_fit->SetParLimits(12,0.5,6);
            }
             h->Fit(f_fit, "LWQN");

            // final fit for this toy
            TFitResultPtr fitResult = h->Fit(f_fit, "SLWQN"); 
            if (fitResult->Status() != 0) {
                h->Fit(f_fit, "SLEWQN");
            }


            //save results 

            if (fitResult->Status() == 0) {
                Double_t mW = f_fit->GetParameter(4);
    Double_t mZ = f_fit->GetParameter(10);
                h_results_mean_W->Fill(mW);
               h_results_mean_Z->Fill(mZ);
    
    Double_t mMean = (mW + mZ) / 2.0;
    h_results_M_mean->Fill(mMean);
            }

            delete h;
            delete f_fit;
            
            if(iToy % 10 == 0) cout << "Completed " << iToy << " / " << nToys << " toys.." << endl;
        }

          //print results
        cout << "\n========================================" << endl;
        cout << "MEAN RESULTS AFTER " << nToys << " TOYS [" << title << "]" << endl;
        cout << "========================================" << endl;
        cout << "W Mean Mean = " << h_results_mean_W->GetMean() << " +/- " << h_results_mean_W->GetMeanError() << " GeV" << endl;
        cout << "W Mean RMS = " << h_results_mean_W->GetRMS() << " GeV" << endl;
        cout << "Z Mean Mean = " << h_results_mean_Z->GetMean() << " +/- " << h_results_mean_Z->GetMeanError() << " GeV" << endl;
        cout << "Z Mean RMS = " << h_results_mean_Z->GetRMS() << " GeV" << endl;
        cout << "M_mean Finale (Mean of Toys) = " << h_results_M_mean->GetMean() << " +/- " << h_results_M_mean->GetMeanError() << " GeV" << endl;
        cout << "M_mean Resolution (RMS del Toy) = " << h_results_M_mean->GetRMS() << " GeV" << endl;
        cout << "========================================\n" << endl;

        // Distribution of results
        TCanvas *c_res = new TCanvas(Form("c_res_%d", iTest), "Distribution fitted parameters", 1200, 500);
        c_res->Divide(2,1);
        
        c_res->cd(1);
        h_results_mean_W->SetLineColor(kBlue);
        h_results_mean_W->Draw();
        
        c_res->cd(2);
        h_results_mean_Z->SetLineColor(kGreen+2);
        h_results_mean_Z->Draw();
        
        c_res->SaveAs(Form("/eos/home-f/fzavatti/eventi6_R0.2_unq6/figPDFB/toy_results_bin%d_idealideal.pdf", iTest+1));

       // distribution of the mean result
TCanvas *c_media_single = new TCanvas(Form("c_media_single_%d", iTest), "Distribution of fitted Average ", 800, 600);
c_media_single->SetGrid();
c_media_single->SetLeftMargin(0.14); 
c_media_single->SetBottomMargin(0.12);

h_results_M_mean->SetStats(0);

h_results_M_mean->GetXaxis()->SetTitle("Fitted #bar{M} [GeV]");
h_results_M_mean->GetYaxis()->SetTitle("Toy MC / Bin");
h_results_M_mean->GetXaxis()->SetTitleSize(0.045);
h_results_M_mean->GetYaxis()->SetTitleSize(0.045);

h_results_M_mean->SetLineColor(kRed);
h_results_M_mean->SetLineWidth(3);
h_results_M_mean->SetFillColorAlpha(kRed, 0.1); 
h_results_M_mean->SetTitle("");

h_results_M_mean->Draw("HIST");

TLegend *leg_media = new TLegend(0.55, 0.68, 0.85, 0.88);
leg_media->SetTextSize(0.032);     
leg_media->SetTextFont(42);         
leg_media->SetFillStyle(0);        

Double_t mean_val = h_results_M_mean->GetMean();
Double_t mean_err = h_results_M_mean->GetMeanError();
Double_t rms_val  = h_results_M_mean->GetRMS();
Double_t rms_err  = h_results_M_mean->GetRMSError();

leg_media->AddEntry(h_results_M_mean, Form("Dati Toy MC (%d fit)", (Int_t)h_results_M_mean->GetEntries()), "f");
leg_media->AddEntry((TObject*)0, Form("Media = %.3f #pm %.3f GeV", mean_val, mean_err), "");
leg_media->AddEntry((TObject*)0, Form("#sigma_{Hist} = %.3f #pm %.3f GeV", rms_val, rms_err), "");
leg_media->SetHeader(Form("Bin %d - %s", iTest+1, plotTitles[iTest].Data()), "C"); 

leg_media->Draw();

c_media_single->SaveAs(Form("/eos/home-f/fzavatti/eventi6_R0.2_unq6/figPDFB/toy_distribution_media_bin%dideal.pdf", iTest+1));

delete leg_media;
delete c_media_single;

 delete f_gen; delete h_results_mean_W; delete h_results_mean_Z; 
    }
}