#include "myLib.h"

using namespace std;

void mkDir(TString dir){
    TString makedir = "mkdir -p " + dir + "/";
    const char *mkDIR = makedir.Data();
    gSystem->Exec(mkDIR);
}

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

    TString extraText = " Sim.";
    TString lumi = "Pythia Tune 4C, #sqrt{s} = 5.02 TeV, #hat{p}_{T}^{min}=150 GeV";

    //cmsColl->DrawLatex(0.20-0.005, 0.942, "");
    //cmsPrel->DrawLatex(0.26-0.005, 0.942, extraText);
    lumiText->DrawLatex(0.97-0.005, 0.942, lumi);
}

void DrawLatex(Float_t x, Float_t y, Int_t color, const char* text, Float_t textSize){
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC();
    latex->SetTextSize(textSize);
    latex->SetTextColor(color);
    latex->SetTextFont(42);
    latex->Draw();
}

void canvasSetup(TCanvas* c){
    c->cd();
    c->SetCanvasSize(600, 800);
    //c->SetCanvasSize(500, 500);

    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.05);
    c->SetTopMargin(0.07);
    c->SetBottomMargin(0.1);
    c->SetFrameLineWidth(2);
    
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelSize(0.04, "XYZ");
    gStyle->SetTitleSize(0.05, "XYZ");
    gStyle->SetTitleOffset(1.2, "Y");
    gStyle->SetOptStat(0);
    gStyle->SetLineWidth(2);
}

void histSetup(TH1* histo, const TString& x_label, const TString& y_label, Color_t color, bool fill){
    //AXIS LABEL
    // ... x axis
    histo->SetTitle("");

    histo->GetXaxis()->SetTitle(x_label);
    histo->GetXaxis()->SetTitleSize(0.04);
    histo->GetXaxis()->SetLabelSize(0.04);
    // histo->GetXaxis()->SetRangeUser(80,180);
    histo->GetXaxis()->CenterTitle();
    // ... y axis
    histo->GetYaxis()->SetTitle(y_label);
    histo->GetYaxis()->SetTitleSize(0.04);
    histo->GetYaxis()->SetLabelSize(0.04);

    //WIDTH & COLOR                                                                                                                                                               
    histo->SetLineWidth(2);
    histo->SetLineColor(color);
    if(fill) {
        // histo->SetFillColor(color);
        histo->SetFillColorAlpha(color, 0.5);
        histo->SetFillStyle(1001);
    }

    // cout<< "Hist: "<< histo->GetName()<<", fill = "<<histo->GetFillColor()<< ", colore = "<< color<<endl;
    histo->SetMarkerStyle(kFullCircle);
    //histo->SetMarkerStyle(24);
    histo->SetMarkerColor(color);
    histo->SetMarkerSize(0.5);
    //histo->SetErrorSize(3);

    
    // AUTO Y RANGE (FIXED STYLE)..QUESTA PARTE CREA PROBLEMI DEVI STARE ATTENTA
   // double ymax = histo->GetMaximum();
   // histo->GetYaxis()->SetRangeUser(0., ymax * 1.2);


}

#include <vector>

void RatioPlot_severalHists(
    TCanvas* c, 
    vector<TH1*> hists, 
    vector<TString> legends,
    vector <double> leg_coord, 
    Double_t ratio_y_inf, 
    Double_t ratio_y_sup, 
    TString ratio_leg,
    Double_t xMin,
    Double_t xMax,
    bool logy_index
) {
    //c->SetCanvasSize(1000, 1200);
    c->SetCanvasSize(1000, 1200);
   //double upperHeight = 0.595;
    //double lowerHeight = 1.0 - upperHeight; // = 0.405
    //c->SetCanvasSize(600, 800);
    c->SetMargin(0, 0, 0, 0); // Rimuove i margini dalla canvas
  
    double upperHeight = 0.7;
    double lowerHeight = 0.3;
    double leftMargin = 0.20;
    double y_offset = 0.92;
    
    TGaxis::SetMaxDigits(3);
    TGaxis::SetExponentOffset(-0.07, 0.0, "y");

    // **Upper pad**
    TPad* up_pad = new TPad("up_pad", "", 0., lowerHeight, 1., 1.);


    //up_pad->SetLogy();
    

    up_pad->SetBottomMargin(0);
    up_pad->SetTopMargin(0.07);
    up_pad->SetLeftMargin(leftMargin);
    up_pad->SetRightMargin(0.05);
    up_pad->Draw();
    up_pad->cd();
    // gPad->SetLogy(1);

    // **Plot degli istogrammi**
    double max_value = -1.;
    for (size_t i = 0; i < hists.size(); i++) {
        hists[i]->SetTitle("");
        hists[i]->SetStats(0);
        max_value = TMath::Max(max_value, hists[i]->GetMaximum());
      /*  if (hists[i]->GetFillColor() < 0.5) {
            if (i == 0) hists[i]->Draw("PE1"); 
            else        hists[i]->Draw("PE1 same");
            // if (i == 0) hists[i]->Draw("hist"); 
            // else        hists[i]->Draw("hist same");
        }
        else{
            if (i == 0) hists[i]->Draw("hist"); 
            else        hists[i]->Draw("hist same");
        }*/
        
          if (i == 0) {
        // === curva di riferimento: RIEMPITA ===
        hists[i]->SetFillStyle(1001);
        hists[i]->SetFillColorAlpha(hists[i]->GetLineColor(), 0.30);
        hists[i]->SetMarkerStyle(0); // niente marker
        hists[i]->Draw("HIST");
    }
    else {
        // === curve di confronto: PALLINI ===
       /*hists[i]->SetFillStyle(0);
        hists[i]->SetMarkerStyle(20 + (i - 1)); // marker diversi
        hists[i]->SetMarkerSize(0.6);
        hists[i]->SetLineWidth(2);*/
        //hists[i]->Draw("E1 P SAME");
        hists[i]->SetFillStyle(1001);
        hists[i]->SetFillColorAlpha(hists[i]->GetLineColor(), 0.30);
        hists[i]->SetMarkerStyle(0); // niente marker
        hists[i]->Draw("HIST SAME");

    }
    }


    hists[0]->GetYaxis()->ChangeLabel(1, -1, -1, -1, -1, -1, " ");

    //hists[0]->GetYaxis()->SetRangeUser(1e-10, 1e-6); //°°°°°°°

    //hists[0]->GetYaxis()->SetRangeUser(1e-7, max_value * 1.4); //°°°°°°°
    
    if (logy_index)
    hists[0]->GetYaxis()->SetRangeUser(1e-7, max_value * 50.0);
else
    hists[0]->GetYaxis()->SetRangeUser(0.0, max_value * 1.4);



    if (logy_index) up_pad->SetLogy(1);

    // hists[0]->GetYaxis()->SetRangeUser(0, 0.25);
    // hists[0]->GetYaxis()->SetRangeUser(0, 0.03);
    hists[0]->GetYaxis()->SetLabelFont(43);
    hists[0]->GetYaxis()->SetLabelSize(38);
    hists[0]->GetYaxis()->SetTitleSize(0.06);
    hists[0]->GetYaxis()->SetTitleOffset(y_offset);
    // hists[0]->GetXaxis()->SetRangeUser(80, 180);
    // hists[0]->GetXaxis()->SetRangeUser(-1.8, 1.8);
    hists[0]->GetXaxis()->SetRangeUser(xMin,xMax);
    // hists[1]->GetXaxis()->SetRangeUser(xMin,xMax);
    // hists[2]->GetXaxis()->SetRangeUser(xMin,xMax);

    hists[0]->GetXaxis()->SetTickLength(0.04);

    // if (strstr(hists[0]->GetXaxis()->GetTitle(), "theta") != nullptr) hists[0]->GetXaxis()->SetRangeUser(1.6, 3.);

    // **Legenda**

    TLegend* legend = new TLegend(leg_coord[0], leg_coord[1], leg_coord[2], leg_coord[3]);

    //TLegend* legend = new TLegend(0.6, 0.7, 0.96, 0.92);
    for (size_t i = 0; i < hists.size(); i++) legend->AddEntry(hists[i], legends[i], "f");

//legend->AddEntry((TObject*)0, "p_{T}^{jet}=40-800 GeV", "");
//legend->AddEntry((TObject*)0, " |#eta| < 2.2", "");

    
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextFont(43);
    legend->SetTextSize(30);
    legend->Draw("same");

    // **Lower pad**
    c->cd();
    TPad* ratio_pad = new TPad("ratio_pad", "", 0., 0.0, 1., lowerHeight);
    ratio_pad->SetTopMargin(0);
    ratio_pad->SetBottomMargin(0.5);
    ratio_pad->SetLeftMargin(leftMargin);
    ratio_pad->SetRightMargin(0.05);
    ratio_pad->SetGridy();
    ratio_pad->Draw();
    ratio_pad->cd();

    vector<TH1F*> h_ratios;
    for (size_t i = 1; i < hists.size(); i++) {

        TH1F* h_ratio = (TH1F*)hists[i]->Clone(Form("h_ratio_%zu", i));
        h_ratio->Divide(hists[0]);
        h_ratio->SetMinimum(ratio_y_inf);
        h_ratio->SetMaximum(ratio_y_sup);
        h_ratios.push_back(h_ratio);
        h_ratio->SetLineColor(hists[i]->GetLineColor());
h_ratio->SetMarkerColor(hists[i]->GetLineColor());

    }

    h_ratios[0]->GetYaxis()->SetTitle(ratio_leg);
    // h_ratios[1]->GetYaxis()->SetTitle(ratio_leg);
    h_ratios[0]->GetYaxis()->SetLabelFont(43);
    h_ratios[0]->GetYaxis()->SetLabelSize(30);
    h_ratios[0]->GetYaxis()->SetTitleSize(30);
    h_ratios[0]->GetYaxis()->SetTitleOffset(y_offset + 1.2);
    h_ratios[0]->GetYaxis()->SetTickLength(0.0375);
    
 
    h_ratios[0]->GetXaxis()->SetTitle(hists[0]->GetXaxis()->GetTitle());



    h_ratios[0]->GetXaxis()->SetLabelSize(0.1);
    h_ratios[0]->GetXaxis()->SetTitleFont(43);
    h_ratios[0]->GetXaxis()->SetTitleSize(33);
    h_ratios[0]->GetXaxis()->SetTitleOffset(1.2);
    h_ratios[0]->GetXaxis()->SetTickLength(0.09);
    h_ratios[0]->GetXaxis()->CenterTitle();
    // if (strstr(h_ratios[0]->GetXaxis()->GetTitle(), "theta") != nullptr) h_ratios[0]->GetXaxis()->SetRangeUser(1.6, 3.);
    // h_ratios[0]->GetXaxis()->SetRangeUser(-1.8, 1.8);
    // h_ratios[0]->GetXaxis()->SetRangeUser(80, 180);
    // h_ratios[0]->GetXaxis()->SetRangeUser(0, 80);
    
    double max_value_ratio = 0;
    for(size_t i = 1; i < h_ratios.size(); i++) max_value_ratio = TMath::Max(max_value_ratio, h_ratios[i]->GetMaximum());
    // h_ratios[0]->GetYaxis()->SetRangeUser(0, max_value_ratio * 2);
    h_ratios[0]->GetXaxis()->SetRangeUser(xMin,xMax);


    h_ratios[0]->Draw("PE");
    for (size_t i = 1; i < h_ratios.size(); i++) h_ratios[i]->Draw("PE same");

    // **Legenda nel pad inferiore**
    TLegend* legend_ratio = new TLegend(0.25, 0.82, 0.40, 0.97);

    for (size_t i = 1; i < hists.size(); i++) {
        TString ratio_legend = legends[i] + " / " + legends[0];
        legend_ratio->AddEntry(h_ratios[i - 1], ratio_legend, "lep");
    }
    legend_ratio->SetFillStyle(0);
    legend_ratio->SetBorderSize(0);
    legend_ratio->SetTextFont(43);
    legend_ratio->SetTextSize(22); // Un po' più piccola rispetto alla legenda superiore
    if(hists.size() > 1) legend_ratio->Draw("same");
    ratio_pad->Modified();
    ratio_pad->Update();
    up_pad->cd();
} //RatioPlot_severalHists()

