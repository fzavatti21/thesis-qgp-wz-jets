// creating4binsgroomed.C
// Reads WZgen.cc output ROOT files from a text list and fills
// 2D histograms of SoftDrop groomed mass vs N2 energy correlator,
// split by pT bin and by event type (W-matched, Z-matched, QCD background).
// The 2D histograms are then projected onto the mass axis with a N2 < 0.18
// cut and passed as input to the four-phase fit (fit4binsfase1-4B.C).
//
// Usage:
//   Edit the default arguments (fileList and output_name) to point
//   to your directory, then run:
//   root creating.C
//
// Input:  text file listing one ROOT file path per line (from pygen.cc)
// Output: ROOT file with all 2D and 1D histograms
//
// NOTE: also update the mkdir path below to match your directory.
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TString.h>

#include <fstream>
#include <iostream>
#include <memory>

using namespace std;

void creating(
    // NOTE: change this path to match your output directory
    TString fileList = "/eos/home-f/fzavatti/eventi6_R0.2_bkgCONsub400pclesSOFTz02_unq/cucu.txt",
    TString output_name = "/eos/home-f/fzavatti/eventi6_R0.2_bkgCONsub400pclesSOFTz02_unq/WZ_histograms4bins.root"
) {
    
    cout << "CREATING W and Z GROOMED MASS HISTOGRAMS" << endl;
    cout << "Reading file list: " << fileList << endl;
    cout << "Output file: " << output_name << endl;
    
    gROOT->SetBatch(kTRUE);
    gSystem->Exec("mkdir -p /eos/home-f/fzavatti/eventi6_R0.2_bkgCONsub400pclesSOFTz02_unq");

    
    // CREATE Histograms
    
    
    TH1D* h_mFat_all = new TH1D("mFat_all", 
                                 " Groomed Mass (all);mass [GeV];Events", 
                                 150, 0, 150);
    h_mFat_all->Sumw2();
    
    TH1D* h_mFat_n2cut = new TH1D("mFat_n2cut", 
                                   " Groomed Mass (N_{2} < 0.18);mass [GeV];Events", 
                                   150, 0, 150);
    h_mFat_n2cut->Sumw2();
    
    TH2F* h2D_mFat_n2 = new TH2F("mFat_n2_2D", 
                                  " Groomed Mass vs N_{2};mass [GeV];N_{2}",
                                  150, 0, 150,
                                  150, 0, 0.6);
    h2D_mFat_n2->Sumw2();
    
    TH2F* h2D_mFat_n2bkg = new TH2F("mFat_n2_2Dbkg", 
                                  " Groomed Mass bkg vs N_{2};mass [GeV];N_{2}",
                                  150, 0, 150,
                                  150, 0, 0.6);
    h2D_mFat_n2bkg->Sumw2();
    
    TH2F* h2D_mFat_n2W = new TH2F("mFat_n2_2DW", 
                                  " Groomed Mass W vs N_{2};mass [GeV];N_{2}",
                                  150, 0, 150,
                                  150, 0, 0.6);
    h2D_mFat_n2W->Sumw2();
    
    TH2F* h2D_mFat_n2Z = new TH2F("mFat_n2_2DZ", 
                                  " Groomed Mass Z vs N_{2};mass [GeV];N_{2}",
                                  150, 0, 150,
                                  150, 0, 0.6);
    h2D_mFat_n2Z->Sumw2();
    
    TH1D* h_mFat_sig = new TH1D("mFat_sig", 
                                 " Groomed Mass (sig);mass [GeV];Events", 
                                 150, 0, 150);
    h_mFat_sig->Sumw2();

    TH1D* h_mFat_bkg = new TH1D("mFat_bkg", 
                                 " Groomed Mass (bkg);mass [GeV];Events", 
                                 150, 0, 150);
    h_mFat_bkg->Sumw2();
    
    // 2D bkg in bin di pT  [150-350, 350-500, 500-700, 700-1000]
    TH2F* h2D_mFat_n2bkg_pt150_350  = new TH2F("mFat_n2_2Dbkg_pt150_350", "",150,0,150,150,0,0.6);
    TH2F* h2D_mFat_n2bkg_pt350_500  = new TH2F("mFat_n2_2Dbkg_pt350_500", "",150,0,150,150,0,0.6);
    TH2F* h2D_mFat_n2bkg_pt500_700  = new TH2F("mFat_n2_2Dbkg_pt500_700", "",150,0,150,150,0,0.6);
    TH2F* h2D_mFat_n2bkg_pt700_1000 = new TH2F("mFat_n2_2Dbkg_pt700_1000","",150,0,150,150,0,0.6);

    // 2D W in bin di pT
    TH2F* h2D_mFat_n2W_pt150_350    = new TH2F("mFat_n2_2DW_pt150_350",   "",150,0,150,150,0,0.6);
    TH2F* h2D_mFat_n2W_pt350_500    = new TH2F("mFat_n2_2DW_pt350_500",   "",150,0,150,150,0,0.6);
    TH2F* h2D_mFat_n2W_pt500_700    = new TH2F("mFat_n2_2DW_pt500_700",   "",150,0,150,150,0,0.6);
    TH2F* h2D_mFat_n2W_pt700_1000   = new TH2F("mFat_n2_2DW_pt700_1000",  "",150,0,150,150,0,0.6);

    // 2D Z in bin di pT
    TH2F* h2D_mFat_n2Z_pt150_350    = new TH2F("mFat_n2_2DZ_pt150_350",   "",150,0,150,150,0,0.6);
    TH2F* h2D_mFat_n2Z_pt350_500    = new TH2F("mFat_n2_2DZ_pt350_500",   "",150,0,150,150,0,0.6);
    TH2F* h2D_mFat_n2Z_pt500_700    = new TH2F("mFat_n2_2DZ_pt500_700",   "",150,0,150,150,0,0.6);
    TH2F* h2D_mFat_n2Z_pt700_1000   = new TH2F("mFat_n2_2DZ_pt700_1000",  "",150,0,150,150,0,0.6);

    // 2D inclusivi in bin di pT
    TH2F* h2D_mFat_n2_pt150_350     = new TH2F("mFat_n2_2D_pt150_350",    "",150,0,150,150,0,0.6);
    TH2F* h2D_mFat_n2_pt350_500     = new TH2F("mFat_n2_2D_pt350_500",    "",150,0,150,150,0,0.6);
    TH2F* h2D_mFat_n2_pt500_700     = new TH2F("mFat_n2_2D_pt500_700",    "",150,0,150,150,0,0.6);
    TH2F* h2D_mFat_n2_pt700_1000    = new TH2F("mFat_n2_2D_pt700_1000",   "",150,0,150,150,0,0.6);

    h2D_mFat_n2bkg_pt150_350->Sumw2();  h2D_mFat_n2bkg_pt350_500->Sumw2();
    h2D_mFat_n2bkg_pt500_700->Sumw2();  h2D_mFat_n2bkg_pt700_1000->Sumw2();
    h2D_mFat_n2W_pt150_350->Sumw2();    h2D_mFat_n2W_pt350_500->Sumw2();
    h2D_mFat_n2W_pt500_700->Sumw2();    h2D_mFat_n2W_pt700_1000->Sumw2();
    h2D_mFat_n2Z_pt150_350->Sumw2();    h2D_mFat_n2Z_pt350_500->Sumw2();
    h2D_mFat_n2Z_pt500_700->Sumw2();    h2D_mFat_n2Z_pt700_1000->Sumw2();
    h2D_mFat_n2_pt150_350->Sumw2();     h2D_mFat_n2_pt350_500->Sumw2();
    h2D_mFat_n2_pt500_700->Sumw2();     h2D_mFat_n2_pt700_1000->Sumw2();

    
    // Counters
    
    Long64_t nEventsTot = 0;
    Long64_t nEventsWithMass = 0;
    Long64_t nEventsN2cut = 0;
    int nFilesProcessed = 0;
    int nFilesTotal = 0;
    
    
    // Loop over files
    
    ifstream list(fileList.Data());
    if (!list.is_open()) {
        cout << "ERROR: Cannot open file list: " << fileList << endl;
        return;
    }
    
    TString fname;
    
    while (list >> fname) {
        nFilesTotal++;
        
        if (fname.IsWhitespace() || fname.BeginsWith("#")) continue;
        
        cout << "\rProcessing file " << nFilesTotal << ": " << fname << flush;
        
        unique_ptr<TFile> f(TFile::Open(fname));
        if (!f || f->IsZombie()) {
            cout << "\n  WARNING: Cannot open " << fname << endl;
            continue;
        }
        
        TTree* tEvent = (TTree*)f->Get("EventTree");
        if (!tEvent) {
            cout << "\n  WARNING: No EventTree in " << fname << endl;
            continue;
        }
        
        
        // SETUP BRANCHES
        
        Float_t mFat_groomed;
        Float_t pFat; 
        Float_t n2Fat;
        Float_t evweight;
        Float_t ev_WMatchedFatJet;
        Float_t ev_ZMatchedFatJet;
        Float_t hofatjet;

        if (!tEvent->GetBranch("mFat_groomed")) {
            cout << "\n  WARNING: No branch 'mFat_groomed' in EventTree" << endl;
            continue;
        }
        if (!tEvent->GetBranch("n2Fat")) {
            cout << "\n  WARNING: No branch 'n2Fat' in EventTree" << endl;
            continue;
        }
        if (!tEvent->GetBranch("evweight")) {
            cout << "\n  WARNING: No branch 'evweight' in EventTree" << endl;
            continue;
        }
        
        tEvent->SetBranchAddress("mFat_groomed", &mFat_groomed);
        tEvent->SetBranchAddress("pFat", &pFat);
        tEvent->SetBranchAddress("n2Fat", &n2Fat);
        tEvent->SetBranchAddress("evweight", &evweight);
        tEvent->SetBranchAddress("ZMatchedFatJet", &ev_ZMatchedFatJet);
        tEvent->SetBranchAddress("WMatchedFatJet", &ev_WMatchedFatJet);
        tEvent->SetBranchAddress("hofatjet", &hofatjet);

        
        // Event loop
        
        Long64_t nEntries = tEvent->GetEntries();
        
        for (Long64_t ie = 0; ie < nEntries; ++ie) {
            tEvent->GetEntry(ie);
            nEventsTot++;
            
            if (mFat_groomed > 0 && hofatjet == 1) {
                nEventsWithMass++;
                
                h_mFat_all->Fill(mFat_groomed, evweight);
                
                if ((ev_ZMatchedFatJet == 1 || ev_WMatchedFatJet == 1) && hofatjet == 1)
                    h_mFat_sig->Fill(mFat_groomed, evweight);
                else
                    h_mFat_bkg->Fill(mFat_groomed, evweight);

                if (n2Fat >= 0 && hofatjet == 1) {
                    h2D_mFat_n2->Fill(mFat_groomed, n2Fat, evweight);

                    // Inclusive pT bins (all event types: W, Z, QCD)
                    if      (pFat >= 150  && pFat < 350)  h2D_mFat_n2_pt150_350->Fill(mFat_groomed, n2Fat, evweight);
                    else if (pFat >= 350  && pFat < 500)  h2D_mFat_n2_pt350_500->Fill(mFat_groomed, n2Fat, evweight);
                    else if (pFat >= 500  && pFat < 700)  h2D_mFat_n2_pt500_700->Fill(mFat_groomed, n2Fat, evweight);
                    else if (pFat >= 700  && pFat < 1000) h2D_mFat_n2_pt700_1000->Fill(mFat_groomed, n2Fat, evweight);

// Split by event type (W-matched, Z-matched, QCD background), integrated over pT
                    if      (ev_WMatchedFatJet == 1)  h2D_mFat_n2W->Fill(mFat_groomed, n2Fat, evweight);
                    else if (ev_ZMatchedFatJet == 1)  h2D_mFat_n2Z->Fill(mFat_groomed, n2Fat, evweight);
                    else                              h2D_mFat_n2bkg->Fill(mFat_groomed, n2Fat, evweight);

// Split by event type AND pT bin
                    if (pFat >= 150 && pFat < 350) {
                        if      (ev_WMatchedFatJet == 1)  h2D_mFat_n2W_pt150_350->Fill(mFat_groomed, n2Fat, evweight);
                        else if (ev_ZMatchedFatJet == 1)  h2D_mFat_n2Z_pt150_350->Fill(mFat_groomed, n2Fat, evweight);
                        else                              h2D_mFat_n2bkg_pt150_350->Fill(mFat_groomed, n2Fat, evweight);
                    }
                    else if (pFat >= 350 && pFat < 500) {
                        if      (ev_WMatchedFatJet == 1)  h2D_mFat_n2W_pt350_500->Fill(mFat_groomed, n2Fat, evweight);
                        else if (ev_ZMatchedFatJet == 1)  h2D_mFat_n2Z_pt350_500->Fill(mFat_groomed, n2Fat, evweight);
                        else                              h2D_mFat_n2bkg_pt350_500->Fill(mFat_groomed, n2Fat, evweight);
                    }
                    else if (pFat >= 500 && pFat < 700) {
                        if      (ev_WMatchedFatJet == 1)  h2D_mFat_n2W_pt500_700->Fill(mFat_groomed, n2Fat, evweight);
                        else if (ev_ZMatchedFatJet == 1)  h2D_mFat_n2Z_pt500_700->Fill(mFat_groomed, n2Fat, evweight);
                        else                              h2D_mFat_n2bkg_pt500_700->Fill(mFat_groomed, n2Fat, evweight);
                    }
                    else if (pFat >= 700 && pFat < 1000) {
                        if      (ev_WMatchedFatJet == 1)  h2D_mFat_n2W_pt700_1000->Fill(mFat_groomed, n2Fat, evweight);
                        else if (ev_ZMatchedFatJet == 1)  h2D_mFat_n2Z_pt700_1000->Fill(mFat_groomed, n2Fat, evweight);
                        else                              h2D_mFat_n2bkg_pt700_1000->Fill(mFat_groomed, n2Fat, evweight);
                    }
                }

                if (n2Fat >= 0 && n2Fat < 0.18 && hofatjet == 1) {
                    h_mFat_n2cut->Fill(mFat_groomed, evweight);
                    nEventsN2cut++;
                }
            }
        }
        
        nFilesProcessed++;
    }
    
    cout << endl;
    list.close();
    
    
    // Print Summary
    
    cout << "\n========================================" << endl;
    cout << "SUMMARY" << endl;
    cout << "========================================" << endl;
    cout << "Files in list: " << nFilesTotal << endl;
    cout << "Files processed: " << nFilesProcessed << endl;
    cout << "Total events: " << nEventsTot << endl;
    cout << "Events with mFat_groomed > 0: " << nEventsWithMass << endl;
    cout << "Events with N2 < 0.18: " << nEventsN2cut << endl;
    cout << "----------------------------------------" << endl;
    cout << "Histogram entries:" << endl;
    cout << "  h_mFat_all: " << h_mFat_all->GetEntries() << endl;
    cout << "  h_mFat_n2cut: " << h_mFat_n2cut->GetEntries() << endl;
    cout << "  h2D_mFat_n2: " << h2D_mFat_n2->GetEntries() << endl;
    cout << "----------------------------------------" << endl;
    cout << "h_mFat_all statistics:" << endl;
    cout << "  Mean: " << h_mFat_all->GetMean() << " GeV" << endl;
    cout << "  RMS: " << h_mFat_all->GetRMS() << " GeV" << endl;
    cout << "========================================\n" << endl;
    
    // Write histograms to output file 

    TFile* output = new TFile(output_name, "RECREATE");
    
    h_mFat_all->Write();
    h_mFat_n2cut->Write();
    h2D_mFat_n2->Write();
    h2D_mFat_n2bkg->Write();
    h2D_mFat_n2W->Write();
    h2D_mFat_n2Z->Write();
    h_mFat_bkg->Write();
    h_mFat_sig->Write();

    h2D_mFat_n2bkg_pt150_350->Write();  h2D_mFat_n2bkg_pt350_500->Write();
    h2D_mFat_n2bkg_pt500_700->Write();  h2D_mFat_n2bkg_pt700_1000->Write();
    h2D_mFat_n2W_pt150_350->Write();    h2D_mFat_n2W_pt350_500->Write();
    h2D_mFat_n2W_pt500_700->Write();    h2D_mFat_n2W_pt700_1000->Write();
    h2D_mFat_n2Z_pt150_350->Write();    h2D_mFat_n2Z_pt350_500->Write();
    h2D_mFat_n2Z_pt500_700->Write();    h2D_mFat_n2Z_pt700_1000->Write();
    h2D_mFat_n2_pt150_350->Write();     h2D_mFat_n2_pt350_500->Write();
    h2D_mFat_n2_pt500_700->Write();     h2D_mFat_n2_pt700_1000->Write();

    output->Close();
    
    // CLEANUP
    delete h_mFat_all;
    delete h_mFat_n2cut;
    delete h2D_mFat_n2;
    delete h2D_mFat_n2bkg;
    delete h2D_mFat_n2W;
    delete h2D_mFat_n2Z;
    delete h_mFat_bkg;
    delete h_mFat_sig;

    delete h2D_mFat_n2bkg_pt150_350;  delete h2D_mFat_n2bkg_pt350_500;
    delete h2D_mFat_n2bkg_pt500_700;  delete h2D_mFat_n2bkg_pt700_1000;
    delete h2D_mFat_n2W_pt150_350;    delete h2D_mFat_n2W_pt350_500;
    delete h2D_mFat_n2W_pt500_700;    delete h2D_mFat_n2W_pt700_1000;
    delete h2D_mFat_n2Z_pt150_350;    delete h2D_mFat_n2Z_pt350_500;
    delete h2D_mFat_n2Z_pt500_700;    delete h2D_mFat_n2Z_pt700_1000;
    delete h2D_mFat_n2_pt150_350;     delete h2D_mFat_n2_pt350_500;
    delete h2D_mFat_n2_pt500_700;     delete h2D_mFat_n2_pt700_1000;

    delete output;
}
