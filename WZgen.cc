// W/Z+jet production simulation with phenomenological heavy-ion quenching model.
// Generates pp events with Pythia8 at sqrt(s)=5.02 TeV, embeds thermal background
// to simulate OO or PbPb collisions, applies constituent subtraction, and computes
// jet substructure observables (N2, tau21, SoftDrop mass) on large-radius jets (R=0.8).
// Output: ROOT TTree with event-level and jet-level variables.
//
// Dependencies: Pythia8, FastJet + contrib, ROOT
// Usage: ./pygen <tune> <seed> <nEvts> <jetR> <massflag> <ns>

#include "fastjet/Selector.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/ClusterSequenceAreaBase.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/contrib/ConstituentSubtractor.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/ModifiedMassDropTagger.hh"
#include "fastjet/contrib/Recluster.hh"
#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>
#include <string>
#include <cstring>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include "Pythia8/Pythia.h"
#include "TTree.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TList.h"
#include "TVector3.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TString.h"
#include "TRandom3.h"
#include "TH1D.h"
#include <ctime>
#include <cstdio>
#include <valarray>
#include "TSystem.h"
#include "fastjet/contrib/EnergyCorrelator.hh"



const double gev_to_fm       = 0.1973;
// Medium length: L = A^{1/3} fm
// For OO:   A=16 : L = 16^{1/3} ≈ 2.52 fm  : current setting
// For PbPb: A=208 : L = 208^{1/3} ≈ 5.93 fm  : replace 16 with 208
const double L_medium        = TMath::Power(16, 1.0/3.0);
// quenching percentage
// For OO:   Q = 0.03  : current setting
// For PbPb: Q = 0.10 : replace 0.03 with 0.10
const double quenching_normal= 0.03;
// Transport coefficient q-hat [GeV^2/fm]
// For OO:   q-hat = 1.0  : current setting
// For PbPb: q-hat = 4.0 : replace 1.0 with 4.0
const double qhat            = 1.0;

using namespace Pythia8;

double Calculate_pX(double pT, double eta, double phi){ return(pT*TMath::Cos(phi)); }
double Calculate_pY(double pT, double eta, double phi){ return(pT*TMath::Sin(phi)); }
double Calculate_pZ(double pT, double eta, double phi){ return(pT*TMath::SinH(eta)); }
double Calculate_E (double pT, double eta, double phi){
    double pZ = Calculate_pZ(pT,eta,phi);
    return(TMath::Sqrt(pT*pT + pZ*pZ));
}

bool EtaCut(fastjet::PseudoJet fjJet, double etaMin, double etaMax){
    if(fjJet.eta()>etaMax || fjJet.eta()<etaMin) return false;
    else return true;
}

Double_t RelativePhi(Double_t mphi, Double_t vphi){
    if      (vphi < -TMath::Pi()) vphi += TMath::TwoPi();
    else if (vphi >  TMath::Pi()) vphi -= TMath::TwoPi();
    if      (mphi < -TMath::Pi()) mphi += TMath::TwoPi();
    else if (mphi >  TMath::Pi()) mphi -= TMath::TwoPi();
    Double_t dphi = mphi-vphi;
    if      (dphi < -TMath::Pi()) dphi += TMath::TwoPi();
    else if (dphi >  TMath::Pi()) dphi -= TMath::TwoPi();
    return dphi;
}

Double_t CalcAngleWFrame(const Event& event, int id1, int id2, int iW){
    double pz_W    = event[iW].pz();
    double m_W     = event[iW].m();
    double E_long_W= TMath::Sqrt(pz_W*pz_W + m_W*m_W);
    double beta    = pz_W/E_long_W;
    double gamma   = E_long_W/m_W;
    double px1=event[id1].px(), py1=event[id1].py(), pz1=event[id1].pz(), m1=event[id1].m();
    double px2=event[id2].px(), py2=event[id2].py(), pz2=event[id2].pz(), m2=event[id2].m();
    double E_long_1= TMath::Sqrt(pz1*pz1+m1*m1);
    double E_long_2= TMath::Sqrt(pz2*pz2+m2*m2);
    double pz1_new = gamma*(pz1-beta*E_long_1);
    double pz2_new = gamma*(pz2-beta*E_long_2);
    double dotProd = px1*px2+py1*py2+pz1_new*pz2_new;
    double norm1   = TMath::Sqrt(px1*px1+py1*py1+pz1_new*pz1_new);
    double norm2   = TMath::Sqrt(px2*px2+py2*py2+pz2_new*pz2_new);
    if(norm1==0 || norm2==0) return 0;
    return TMath::ACos(dotProd/(norm1*norm2));
}

void IterativeDeclustering(fastjet::PseudoJet jet, double jetR, double weight){
    double erad2=0,r2=0,kt2=0,kt1=0,erad1=0,r1=0;
    int flagtag=0,flagsubjetSD=0,flag1=0,flag2=0,tagpart=-1;
    double jet_radius_ca=100.0;
    fastjet::JetDefinition jet_def(fastjet::genkt_algorithm,jet_radius_ca,0,
        static_cast<fastjet::RecombinationScheme>(0),fastjet::Best);
    try{
        std::vector<fastjet::PseudoJet> particles=jet.constituents();
        fastjet::ClusterSequence cs_gen(particles,jet_def);
        std::vector<fastjet::PseudoJet> output_jets=cs_gen.inclusive_jets(0);
        output_jets=sorted_by_pt(output_jets);
        fastjet::PseudoJet jj=output_jets[0];
        fastjet::PseudoJet j1,j2;
        while(jj.has_parents(j1,j2)){
            flagtag=0;
            if(j1.perp()<j2.perp()) std::swap(j1,j2);
            double delta_R=j1.delta_R(j2);
            double cut=1;
            double kt=j2.perp()*TMath::Sin(delta_R);
            double erad=j1.e();
            double zg=j2.perp()/(j1.perp()+j2.perp());
            if(kt>cut){ erad1=erad; r1=delta_R; flag1=flagtag; kt1=j2.perp()/j1.perp(); }
            if(zg>0.1 && flagsubjetSD==0){ erad2=erad; r2=delta_R; flag2=flagtag; kt2=kt; flagsubjetSD=1; }
            jj=j1;
        }
    } catch(fastjet::Error){}
}


int main(int argc, char* argv[])

{
    Int_t cislo=-1, tune=-1, massflag=0, ns=0;
    Float_t pthatmin=5, pthatmax=5020, jetR=0.4;
    double efficienzaHF=0.8;

    if(argc!=7){
        cout<<"Usage:"<<endl<<"./pygen <PythiaTune> <Number> <nEvts> <jetR> <massflag> <ns>"<<endl;
        return 0;
    }
    tune    =atoi(argv[1]);
    cislo   =atoi(argv[2]);
    Int_t nEvents=atoi(argv[3]);
    // Note: jetR is accepted as input argument but not used in the analysis.
    // All jet finding is performed with jetRfat=0.8 (hardcoded below).
    // jetR is kept in the argument list to allow future extensions without changing the calling scripts.
    jetR    =atof(argv[4]);
    massflag=atof(argv[5]);
    ns      =atoi(argv[6]);

    TTree *tEvent = new TTree("EventTree","Event-level variables");

    float ev_eventID, ev_pTHat, Hjet_pt;
    float ev_ptW, ev_etaW, ev_phiW;
    float id1, id2;
    float ev_thetaqq, ev_dRqq, ev_mqq;
    float evweight;
    float evtau21Fat, mFat, pFat;
    float ev_boostW, ev_WdecayTime, ev_WcoherenceTime, ev_WtotalTime;
    float c2Fat, d2Fat, n2Fat, n3Fat;
    float n2Fat_pure;  
    float ev_isW, ev_hasHFFatJet, ev_WMatchedFatJet, hofatjet;
    float ev_ptZ, ev_etaZ, ev_phiZ;
    float ev_boostZ, ev_ZdecayTime, ev_ZcoherenceTime, ev_ZtotalTime;
    float ev_isZ, ev_thetaqq_Z, ev_dRqq_Z, ev_mqq_Z, ev_ZMatchedFatJet;
    float mFat_pure,ev_WMatchedprobe, ev_ZMatchedprobe, hoprobe;
    float mFat_groomed,mFatpure_groomed;

    tEvent->Branch("eventID",        &ev_eventID,        "eventID/F");
    tEvent->Branch("pTHat",          &ev_pTHat,          "pTHat/F");
    tEvent->Branch("ptHJet",         &Hjet_pt,           "ptHJet/F");
    tEvent->Branch("ptW",            &ev_ptW,            "ptW/F");
    tEvent->Branch("etaW",           &ev_etaW,           "etaW/F");
    tEvent->Branch("phiW",           &ev_phiW,           "phiW/F");
    tEvent->Branch("id1",            &id1,               "id1/F");
    tEvent->Branch("id2",            &id2,               "id2/F");
    tEvent->Branch("thetaqq",        &ev_thetaqq,        "thetaqq/F");
    tEvent->Branch("dRqq",           &ev_dRqq,           "dRqq/F");
    tEvent->Branch("mqq",            &ev_mqq,            "mqq/F");
    tEvent->Branch("evweight",       &evweight,          "evweight/F");
    tEvent->Branch("mFat",           &mFat,              "mFat/F");
    tEvent->Branch("evtau21Fat",     &evtau21Fat,        "evtau21Fat/F");
    tEvent->Branch("pFat",           &pFat,              "pFat/F");
    tEvent->Branch("isWEvent",       &ev_isW,            "isWEvent/F");
    tEvent->Branch("hasHFFatJet",    &ev_hasHFFatJet,    "hasHFFatJet/F");
    tEvent->Branch("WMatchedFatJet", &ev_WMatchedFatJet, "WMatchedFatJet/F");
    tEvent->Branch("boostW",         &ev_boostW,         "boostW/F");
    tEvent->Branch("WdecayTime",     &ev_WdecayTime,     "WdecayTime/F");
    tEvent->Branch("WcoherenceTime", &ev_WcoherenceTime, "WcoherenceTime/F");
    tEvent->Branch("WtotalTime",     &ev_WtotalTime,     "WtotalTime/F");
    tEvent->Branch("c2Fat",          &c2Fat,             "c2Fat/F");
    tEvent->Branch("d2Fat",          &d2Fat,             "d2Fat/F");
    tEvent->Branch("n2Fat",          &n2Fat,             "n2Fat/F");
    tEvent->Branch("n3Fat",          &n3Fat,             "n3Fat/F");
    tEvent->Branch("n2Fat_pure",     &n2Fat_pure,        "n2Fat_pure/F");
    tEvent->Branch("ptZ",            &ev_ptZ,            "ptZ/F");
    tEvent->Branch("etaZ",           &ev_etaZ,           "etaZ/F");
    tEvent->Branch("phiZ",           &ev_phiZ,           "phiZ/F");
    tEvent->Branch("boostZ",         &ev_boostZ,         "boostZ/F");
    tEvent->Branch("ZdecayTime",     &ev_ZdecayTime,     "ZdecayTime/F");
    tEvent->Branch("ZcoherenceTime", &ev_ZcoherenceTime, "ZcoherenceTime/F");
    tEvent->Branch("ZtotalTime",     &ev_ZtotalTime,     "ZtotalTime/F");
    tEvent->Branch("isZEvent",       &ev_isZ,            "isZEvent/F");
    tEvent->Branch("thetaqq_Z",      &ev_thetaqq_Z,      "thetaqq_Z/F");
    tEvent->Branch("dRqq_Z",         &ev_dRqq_Z,         "dRqq_Z/F");
    tEvent->Branch("mqq_Z",          &ev_mqq_Z,          "mqq_Z/F");
    tEvent->Branch("ZMatchedFatJet", &ev_ZMatchedFatJet, "ZMatchedFatJet/F");
    tEvent->Branch("hofatjet",       &hofatjet,          "hofatjet/F");
    tEvent->Branch("mFat_pure", &mFat_pure, "mFat_pure/F");
     tEvent->Branch("ZMatchedprobe", &ev_ZMatchedprobe, "ZMatchedprobe/F");
    tEvent->Branch("WMatchedprobe", &ev_WMatchedprobe, "WMatchedprobe/F");
    tEvent->Branch("hoprobe",       &hoprobe,          "hoprobe/F");
    tEvent->Branch("mFat_groomed", &mFat_groomed, "mFat_groomed/F");
        tEvent->Branch("mFatpure_groomed", &mFatpure_groomed, "mFatpure_groomed/F");



    double trackEtaCut  =2.4;
    double trackLowPtCut=0.;
    

    Int_t mypart=-1;
    TString name;
    Pythia pythia;
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eCM = 5020.");
    pythia.readString("Tune:pp = 5");
    pythia.readString("Random:setSeed = on");
    pythia.readString(Form("Random:seed = %d",cislo));
    pythia.readString("PartonLevel:MPI = off");
    pythia.readString("PartonLevel:ISR = off");
    // Process selection via massflag argument:
// 1: gg->qqbar + qq->qq (light quark QCD)
// 2: gg->qqbar only
// 3: qq->qq only
// 4: gg->gg only
// 5: gg->ccbar (charm)
// 6: HardQCD + W/Z+parton (SIGNAL: used in the thesis analysis)
// 7: qg->qg + qq->qq with hadron decays off
    if(massflag==1){ pythia.readString("HardQCD:gg2qqbar = on"); pythia.readString("HardQCD:qq2qq = on"); }
    if(massflag==2)  pythia.readString("HardQCD:gg2qqbar = on");
    if(massflag==3)  pythia.readString("HardQCD:qq2qq = on");
    if(massflag==4)  pythia.readString("HardQCD:gg2gg = on");
    if(massflag==5)  pythia.readString("HardQCD:gg2ccbar = on");
    if(massflag==6){
        pythia.readString("HardQCD:all = on");
        pythia.readString("WeakBosonAndParton:qqbar2Wg = on");
        pythia.readString("WeakBosonAndParton:qg2Wq = on");
        pythia.readString("WeakBosonAndParton:qqbar2gmZg = on");
        pythia.readString("WeakBosonAndParton:qg2gmZq = on");
    }
    if(massflag==7){
        pythia.readString("HardQCD:qg2qg = on");
        pythia.readString("HardQCD:qq2qq = on");
        pythia.readString("HadronLevel:decay = off");
    }
    if(ns==0) pythia.readString("TimeShower:nGluonToQuark=3");

    pythia.readString("PhaseSpace:pTHatMin = 150"); 
    // Phase space bias to enhance high-pT events:
// events are weighted by (pThat/ref)^pow to sample the hard tail efficiently.
// The weight is stored in evweight and must be used in all histograms.
    pythia.readString("PhaseSpace:bias2Selection = on");
    pythia.readString("PhaseSpace:bias2SelectionPow = 5.");
    pythia.readString("PhaseSpace:bias2SelectionRef = 8.");
    pythia.init();

    fastjet::contrib::NsubjettinessRatio nSub21_beta1(2,1,fastjet::contrib::OnePass_KT_Axes(),fastjet::contrib::UnnormalizedMeasure(1));
    fastjet::contrib::Nsubjettiness     nSub1_beta1  (1,  fastjet::contrib::OnePass_KT_Axes(),fastjet::contrib::UnnormalizedMeasure(1));
    fastjet::contrib::Nsubjettiness     nSub2_beta1  (2,  fastjet::contrib::OnePass_KT_Axes(),fastjet::contrib::UnnormalizedMeasure(1));
    fastjet::contrib::EnergyCorrelatorC2 C2_beta1(1.0);
    fastjet::contrib::EnergyCorrelatorD2 D2_beta1(1.0);
    fastjet::contrib::EnergyCorrelatorN2 N2_beta1(1.0);
    fastjet::contrib::EnergyCorrelatorN3 N3_beta1(1.0);

    double etamin_Sig=-trackEtaCut+jetR;
    double etamax_Sig=-etamin_Sig;
    fastjet::Strategy            strategy    =fastjet::Best;
    fastjet::RecombinationScheme recombScheme=fastjet::E_scheme;
    fastjet::JetDefinition *jetDefAKT_Sig=new fastjet::JetDefinition(fastjet::antikt_algorithm,jetR,    recombScheme,strategy);
    double jetRfat=0.8;
    double etamin_Fat=-trackEtaCut+jetRfat;
    double etamax_Fat=-etamin_Fat;
    fastjet::JetDefinition *jetDefAKT_Fat=new fastjet::JetDefinition(fastjet::antikt_algorithm,jetRfat, recombScheme,strategy);
    fastjet::GhostedAreaSpec ghostareaspec(trackEtaCut,1,0.01);
    fastjet::AreaType        areaType=fastjet::active_area_explicit_ghosts;
    fastjet::AreaDefinition *areaDef =new fastjet::AreaDefinition(areaType,ghostareaspec);

    std::vector<fastjet::PseudoJet> fjInputs1;
    std::vector<fastjet::PseudoJet> fjInputs2;
    std::vector<fastjet::PseudoJet> fjInputsTot;

    Int_t count=0;
    int nNoPair=0,nLowPtW=0;
    int nEvents_with_W=0,nEvents_with_Z=0;
    int nParticles_from_W_reduced=0,nParticles_from_W_no=0;
    int nParticles_from_Z_reduced=0,nParticles_from_Z_no=0;
    int nParticles_normal=0;

    // Thermal background parameters (tuned to reproduce OO/PbPb multiplicity):
    // naverage: mean number of thermal particles per event
    // sigman:   Gaussian spread around the mean
    // beta_thermal: slope parameter of the Gamma pT distribution [GeV]

    double naverage = 400.;
    double sigman   = 12.;
    TRandom3 *r3 = new TRandom3(0);
    double beta_thermal = 0.34;
    bool doTuning = false;

    TF1 *f_pT = new TF1("f_pT", "TMath::GammaDist(x, 2, 0, [0])", 0, 20);
    f_pT->SetParameter(0, beta_thermal);
    TF1 *f_eta = new TF1("f_eta", "1", -trackEtaCut, trackEtaCut);
    TF1 *f_phi = new TF1("f_phi", "1", -TMath::Pi(), TMath::Pi());

    TH1F *hDeltaPt = new TH1F("hDeltaPt","deltaPt;#DeltapT (GeV);counts",200,-100.,100.);
    TH1F *hT_pT  = new TH1F("hT_pT",  "pT termali;pT (GeV);counts", 100, 0., 10.);
    TH1F *hT_eta = new TH1F("hT_eta", "eta termali;#eta;counts",     100, -3., 3.);
    TH1F *hT_phi = new TH1F("hT_phi", "phi termali;#phi;counts",     100, -4., 4.);

    for(int i=0; i<nEvents; i++)
    {
        double fourvec[4];
        double fourvecp[4];
        if(!pythia.next()) continue;
        count++;

        double weight=pythia.info.weight();
        double pthat =pythia.info.pTHat();
        evweight=weight;

        fjInputs1.resize(0);
        fjInputs2.resize(0);
        fjInputsTot.resize(0);

        if(pythia.event.size()==0) continue;

        double ptW=-1.,phiW=-999.,etaW=-999.,boostW=-1;
        double WdecayTime=-1,WcoherenceTime=-1,WtotalTime=-1;
        int Windex=-1;
        double ptZ=-1.,phiZ=-999.,etaZ=-999.,boostZ=-1;
        double ZdecayTime=-1,ZcoherenceTime=-1,ZtotalTime=-1;
        int Zindex=-1;
        double dPhiqq=-999.,thetaqq=-999.,Mqq=-999.,dRqq=-999.;
        double dPhiqq_Z=-999.,thetaqq_Z=-999.,Mqq_Z=-999.,dRqq_Z=-999.;

        for(int ii=0; ii<pythia.event.size(); ++ii){
            int absID=TMath::Abs(pythia.event[ii].id());
            if(absID==24 && Windex==-1){ Windex=ii; ptW=pythia.event[ii].pT(); etaW=pythia.event[ii].eta(); phiW=pythia.event[ii].phi(); }
            if(absID==23 && Zindex==-1){ Zindex=ii; ptZ=pythia.event[ii].pT(); etaZ=pythia.event[ii].eta(); phiZ=pythia.event[ii].phi(); }
        }
        ev_isW=(Windex>=0)?1:0;
        ev_isZ=(Zindex>=0)?1:0;

        if(ev_isW){
            // W decay time: tau = gamma * hbar/Gamma_W, sampled from exponential distribution
            // Gamma_W = 2.085 GeV (W decay width), boost = mT/m
            // hbar_c = 0.1973 GeV*fm converts to fm
            const auto &W=pythia.event[Windex];
            double mT_W=W.mT(),m_W=W.m();
            boostW=mT_W/m_W;
            double W_tau=1.0/2.085,hbar_c=0.1973;
            double W_tau_mean=boostW*W_tau*hbar_c;
            TRandom3 *W_rndm=new TRandom3(cislo+i);
            WdecayTime=W_rndm->Exp(W_tau_mean);
            delete W_rndm;
            int index=Windex; int d1=0,d2=0;
            do{ d1=pythia.event[index].daughter1(); d2=pythia.event[index].daughter2(); index=d1; }
            while(TMath::Abs(pythia.event[index].id())==24);
            dPhiqq=RelativePhi(pythia.event[d1].phi(),pythia.event[d2].phi());
            double E1=pythia.event[d1].e(),px1=pythia.event[d1].px(),py1=pythia.event[d1].py(),pz1=pythia.event[d1].pz();
            double E2=pythia.event[d2].e(),px2=pythia.event[d2].px(),py2=pythia.event[d2].py(),pz2=pythia.event[d2].pz();
            double Etot=E1+E2,pxtot=px1+px2,pytot=py1+py2,pztot=pz1+pz2;
            double M2=Etot*Etot-(pxtot*pxtot+pytot*pytot+pztot*pztot);
            Mqq=(M2>0?TMath::Sqrt(M2):0);
            double deta=pythia.event[d1].eta()-pythia.event[d2].eta();
            dRqq=TMath::Sqrt(deta*deta+dPhiqq*dPhiqq);
            thetaqq=CalcAngleWFrame(pythia.event,d1,d2,Windex);
            // Decoherence time: tau_d = (12 / (qhat * theta_qq^2))^{1/3}
            // qhat divided by hbar_c^2 converts from GeV^2/fm to fm^3 
            // result is in fm
            double qhat_loc=1.0,hbar_c_sq=0.1973*0.1973;
            if(thetaqq>0) WcoherenceTime=TMath::Power(12.0/(4*qhat_loc/hbar_c_sq*thetaqq*thetaqq),1.0/3.0);
            else WcoherenceTime=-1;
            WtotalTime=WdecayTime+WcoherenceTime;
        }

        // Identical procedure for the Z boson (Gamma_Z = 2.495 GeV)
        if(ev_isZ){
            const auto &Z=pythia.event[Zindex];
            double mT_Z=Z.mT(),m_Z=Z.m();
            boostZ=mT_Z/m_Z;
            double Z_tau=1.0/2.495,hbar_c=0.1973;
            double Z_tau_mean=boostZ*Z_tau*hbar_c;
            TRandom3 *Z_rndm=new TRandom3(cislo+i+10000);
            ZdecayTime=Z_rndm->Exp(Z_tau_mean);
            delete Z_rndm;
            int index=Zindex; int d1=0,d2=0;
            do{ d1=pythia.event[index].daughter1(); d2=pythia.event[index].daughter2(); index=d1; }
            while(pythia.event[index].id()==23);
            dPhiqq_Z=RelativePhi(pythia.event[d1].phi(),pythia.event[d2].phi());
            double E1=pythia.event[d1].e(),px1=pythia.event[d1].px(),py1=pythia.event[d1].py(),pz1=pythia.event[d1].pz();
            double E2=pythia.event[d2].e(),px2=pythia.event[d2].px(),py2=pythia.event[d2].py(),pz2=pythia.event[d2].pz();
            double Etot=E1+E2,pxtot=px1+px2,pytot=py1+py2,pztot=pz1+pz2;
            double M2=Etot*Etot-(pxtot*pxtot+pytot*pytot+pztot*pztot);
            Mqq_Z=(M2>0?TMath::Sqrt(M2):0);
            double deta=pythia.event[d1].eta()-pythia.event[d2].eta();
            dRqq_Z=TMath::Sqrt(deta*deta+dPhiqq_Z*dPhiqq_Z);
            thetaqq_Z=CalcAngleWFrame(pythia.event,d1,d2,Zindex);
            double qhat_loc=1.0,hbar_c_sq=0.1973*0.1973;
            if(thetaqq_Z>0) ZcoherenceTime=TMath::Power(12.0/(4*qhat_loc/hbar_c_sq*thetaqq_Z*thetaqq_Z),1.0/3.0);
            else ZcoherenceTime=-1;
            ZtotalTime=ZdecayTime+ZcoherenceTime;
        }

        if(Windex>=0) nEvents_with_W++;
        if(Zindex>=0) nEvents_with_Z++;

        const double rescale_mu   =quenching_normal;
        const double rescale_sigma=0.15*rescale_mu;
        TRandom3 *rescale_rndm=new TRandom3(cislo+i+20000);

        for(int j=0; j<pythia.event.size(); j++){
            if(!pythia.event[j].isFinal()) continue;
            if(pythia.event[j].pT()<trackLowPtCut) continue;
            if(TMath::Abs(pythia.event[j].eta())>trackEtaCut) continue;

            bool fromW=false,fromZ=false;
            // Walk up the decay chain to check if this final-state particle originates from a W (id=24) or Z (id=23).
            // The vector is used as a growing stack: mothers are added as we go up.
            std::vector<int> mothers;
            mothers.push_back(j);
            for(unsigned int iMother=0; iMother<mothers.size(); ++iMother){
                int currentID=TMath::Abs(pythia.event[mothers[iMother]].id());
                if(currentID==24){ fromW=true; break; }
                if(currentID==23){ fromZ=true; break; }
                int m1=pythia.event[mothers[iMother]].mother1();
                int m2=pythia.event[mothers[iMother]].mother2();
                if(m1>2) mothers.push_back(m1);
                if(m2!=m1 && m2>2) mothers.push_back(m2);
            }

            double quenching_factor=1.0,rescale=0.0;

            // Quenching scenario: uncomment the desired block 

            //  SCENARIO 1: UNQUENCHED (pp baseline) 
            // Uncomment this block and comment Scenarios 2 and 3
            /*
            fourvecp[0]=pythia.event[j].px();
            fourvecp[1]=pythia.event[j].py();
            fourvecp[2]=pythia.event[j].pz();
            fourvecp[3]=pythia.event[j].e();
            */

// SCENARIO 2: PARTIAL QUENCHING (default) 
// W/Z daughters quenched only after decay + decoherence time
// Uncomment this block and comment Scenarios 1 and 3
if(fromW){
    if(WtotalTime>L_medium){ rescale=0.0; nParticles_from_W_no++; }
    else{
        double frac=(L_medium-WtotalTime)/L_medium;
        rescale=rescale_rndm->Gaus(quenching_normal*frac,rescale_sigma);
        rescale=(rescale>0)?rescale:0.0; nParticles_from_W_reduced++;
    }
    quenching_factor=1.0-rescale;
}
else if(fromZ){
    if(ZtotalTime>L_medium){ rescale=0.0; nParticles_from_Z_no++; }
    else{
        double frac=(L_medium-ZtotalTime)/L_medium;
        rescale=rescale_rndm->Gaus(quenching_normal*frac,rescale_sigma);
        rescale=(rescale>0)?rescale:0.0; nParticles_from_Z_reduced++;
    }
    quenching_factor=1.0-rescale;
}
else{
    rescale=rescale_rndm->Gaus(rescale_mu,rescale_sigma);
    rescale=(rescale>0)?rescale:0.0;
    quenching_factor=1.0-rescale;
    nParticles_normal++;
}
fourvecp[0]=pythia.event[j].px()*quenching_factor;
fourvecp[1]=pythia.event[j].py()*quenching_factor;
fourvecp[2]=pythia.event[j].pz()*quenching_factor;
fourvecp[3]=pythia.event[j].e() *quenching_factor;

// SCENARIO 3: TOTAL QUENCHING 
// All particles quenched equally, no time structure
// Uncomment this block and comment Scenarios 1 and 2
/*
rescale=rescale_rndm->Gaus(rescale_mu,rescale_sigma);
rescale=(rescale>0)?rescale:0.0;
quenching_factor=1.0-rescale;
fourvecp[0]=pythia.event[j].px()*quenching_factor;
fourvecp[1]=pythia.event[j].py()*quenching_factor;
fourvecp[2]=pythia.event[j].pz()*quenching_factor;
fourvecp[3]=pythia.event[j].e() *quenching_factor;
*/


            fastjet::PseudoJet PythiaParticlep(fourvecp);
            PythiaParticlep.set_user_index(0);
            int absID=TMath::Abs(pythia.event[j].id());
            // Tag heavy-flavour hadrons (user_index=1) for HF contamination studies:
// B mesons: 511-544, D mesons: 411-434
// B baryons: 5112-5553, D baryons: 4112-4443
            if(absID>=511  && absID<545)  PythiaParticlep.set_user_index(1);
            if(absID>=411  && absID<435)  PythiaParticlep.set_user_index(1);
            if(absID>=5112 && absID<5554) PythiaParticlep.set_user_index(1);
            if(absID>=4112 && absID<4444) PythiaParticlep.set_user_index(1);

            fjInputs1.push_back(PythiaParticlep);
        }
        delete rescale_rndm;

        // ── PROBE: hardest jet da fjInputs1 (Pythia puro) ──────────
        fastjet::PseudoJet probe;
        bool hasProbe=false;
        
          // 1. CREATE THE PROBE (Pure Jet, without background)
fastjet::ClusterSequence cs_probe(fjInputs1, *jetDefAKT_Fat);
std::vector<fastjet::PseudoJet> probeJets = sorted_by_pt(cs_probe.inclusive_jets(0.0)); 
            


        //  Generate TERMAL Background
        fjInputsTot=fjInputs1;
        Int_t nThermalParticles=(Int_t)r3->Gaus(naverage,sigman);
        if(nThermalParticles<0) nThermalParticles=0;

        for(int j=0; j<nThermalParticles; j++){
            double pT  = f_pT->GetRandom();
            double eta = f_eta->GetRandom();
            double phi = f_phi->GetRandom();
            if(pT<trackLowPtCut) continue;
            hT_pT ->Fill(pT);
            hT_eta->Fill(eta);
            hT_phi->Fill(phi);
            fourvec[0]=Calculate_pX(pT,eta,phi);
            fourvec[1]=Calculate_pY(pT,eta,phi);
            fourvec[2]=Calculate_pZ(pT,eta,phi);
            fourvec[3]=Calculate_E (pT,eta,phi);
            fastjet::PseudoJet ThermalParticle(fourvec);
            ThermalParticle.set_user_index(0);
            fjInputs2.push_back(ThermalParticle);
            fjInputsTot.push_back(ThermalParticle);
        }

        //  Background Subtraction 
        fastjet::JetMedianBackgroundEstimator bge;
        fastjet::Selector BGSelector=fastjet::SelectorAbsEtaMax(2.4);
        fastjet::JetDefinition jetDefBG(fastjet::kt_algorithm,jetRfat,recombScheme,strategy);
        fastjet::AreaDefinition fAreaDefBG(fastjet::active_area_explicit_ghosts,ghostareaspec);
        // Background estimation from thermal particles only (fjInputs2),
// not from the full event, to avoid biasing rho with hard Pythia particles.
        fastjet::ClusterSequenceArea clustSeqBG(fjInputs2,jetDefBG,fAreaDefBG);
        std::vector<fastjet::PseudoJet> BGJets=clustSeqBG.inclusive_jets();
        bge.set_selector(BGSelector);
        bge.set_jets(BGJets);

        fastjet::contrib::ConstituentSubtractor subtractor(&bge);
        subtractor.set_common_bge_for_rho_and_rhom(true);
        subtractor.set_max_standardDeltaR(0.45);

        // Subtraction on fjInputsTot : fjInputs3
        std::vector<fastjet::PseudoJet> fjInputs3=subtractor.subtract_event(fjInputsTot,2.4);
    

        //  JET FINDING on fjInputs3 
        fastjet::GhostedAreaSpec New_ghost_spec(2.4,1,0.05);
        fastjet::AreaDefinition  New_fAreaDef(fastjet::passive_area,New_ghost_spec);
        fastjet::ClusterSequenceArea New_clustSeq(fjInputs3,*jetDefAKT_Fat,New_fAreaDef);
        std::vector<fastjet::PseudoJet> PythiaJetsFat=New_clustSeq.inclusive_jets(0);
        PythiaJetsFat=sorted_by_pt(PythiaJetsFat);

        //  ANALISIS HARDEST FAT JET 
        try{
            if(PythiaJetsFat.size()==0) continue;

            evtau21Fat=-1; mFat=-1; pFat=-1;mFat_groomed = -1; mFatpure_groomed = -1;
            ev_hasHFFatJet=0; ev_WMatchedFatJet=0; ev_ZMatchedFatJet=0;
             ev_WMatchedprobe=0; ev_ZMatchedprobe=0;

            c2Fat=-1; d2Fat=-1; n2Fat=-1; n3Fat=-1; hofatjet=0,hoprobe=0;  // inizializza a -1
     probe = probeJets[0];

            if(probeJets.size() > 0 && TMath::Abs(probeJets[0].eta()) <= etamax_Fat && probeJets[0].perp() >= 150.){
                 mFat_pure = probeJets[0].m();
        n2Fat_pure = N2_beta1(probeJets[0]);
                        hasProbe=true;
                        hoprobe=1;

    } else {
        n2Fat_pure = -1;
        mFat_pure  = -1;
    }
    
            fastjet::PseudoJet hardestFat=PythiaJetsFat[0];

            if(TMath::Abs(hardestFat.eta())<=etamax_Fat && hardestFat.perp()>=150.){

                hofatjet=1;
                evtau21Fat=nSub21_beta1(hardestFat);
                mFat=hardestFat.m();
                pFat=hardestFat.perp();

                c2Fat=C2_beta1(hardestFat);
                d2Fat=D2_beta1(hardestFat);
                n2Fat=N2_beta1(hardestFat);  
                n3Fat=N3_beta1(hardestFat);
 

                //SOFT DROP MASS BLOCK
                fastjet::contrib::SoftDrop sd(0, 0.2);  // beta=0, zcut=0.2
                fastjet::PseudoJet groomedJet = sd(hardestFat);
                if (groomedJet != 0) mFat_groomed = groomedJet.m();
                else                 mFat_groomed = -1;

                fastjet::PseudoJet groomedProbe = sd(probe);
                if (groomedProbe != 0) mFatpure_groomed = groomedProbe.m();
                else                 mFatpure_groomed = -1;


                int nHF_in_hardest=0;
                const auto &constituents=hardestFat.constituents();
                for(const auto &p:constituents){
                    if(p.user_index()==1) nHF_in_hardest++;
                }
                ev_hasHFFatJet=(nHF_in_hardest==1)?1:0;

                 if(ev_isW==1){
                    double dEta=etaW-hardestFat.eta();
                    double dPhi=RelativePhi(phiW,hardestFat.phi());
                    double dR=TMath::Sqrt(dEta*dEta+dPhi*dPhi);
                    // Match W to hardest fat jet within dR < R_fat = 0.8
                    if(dR<0.8) ev_WMatchedFatJet=1;
                    
                }

                if(ev_isZ==1){
                    double dEta=etaZ-hardestFat.eta();
                    double dPhi=RelativePhi(phiZ,hardestFat.phi());
                    double dR=TMath::Sqrt(dEta*dEta+dPhi*dPhi);
                    if(dR<0.8) ev_ZMatchedFatJet=1;
                }


                
            }
            if(hoprobe==1 && ev_isW==1){
            double dEtap=etaW-probe.eta();
                    double dPhip=RelativePhi(phiW,probe.phi());
                    double dRp=TMath::Sqrt(dEtap*dEtap+dPhip*dPhip);
                    if(dRp<0.8) ev_WMatchedprobe=1;
            }
            if(hoprobe==1 && ev_isZ==1){
                    double dEtap=etaZ-probe.eta();
                    double dPhip=RelativePhi(phiZ,probe.phi());
                    double dRp=TMath::Sqrt(dEtap*dEtap+dPhip*dPhip);
                    if(dRp<0.8) ev_ZMatchedprobe=1;
                    }


                  
            ev_eventID=i; ev_pTHat=pthat;
            ev_ptW=ptW; ev_etaW=etaW; ev_phiW=phiW;
            ev_thetaqq=thetaqq; ev_mqq=Mqq; ev_dRqq=dRqq;
            ev_boostW=boostW; ev_WdecayTime=WdecayTime;
            ev_WcoherenceTime=WcoherenceTime; ev_WtotalTime=WtotalTime;
            ev_ptZ=ptZ; ev_etaZ=etaZ; ev_phiZ=phiZ;
            ev_boostZ=boostZ; ev_ZdecayTime=ZdecayTime;
            ev_ZcoherenceTime=ZcoherenceTime; ev_ZtotalTime=ZtotalTime;
            ev_thetaqq_Z=thetaqq_Z; ev_mqq_Z=Mqq_Z; ev_dRqq_Z=dRqq_Z;
            id1=TMath::Abs(pythia.event[5].id());
            id2=TMath::Abs(pythia.event[6].id());
            Hjet_pt=0;
            
            tEvent->Fill();

        } catch(fastjet::Error){}

    } // end event loop 

    pythia.stat();

    // Output path on CERN EOS. Adapt to local directory if running outside CERN.
    TString outDir =Form("/eos/home-f/fzavatti/eventi%d_R%.1f_OOquenchdelayed",massflag,jetR);
    gSystem->mkdir(outDir,kTRUE);
    TString outPath=Form("%s/Output_single_num%d_jetR%.2f_proc%d_gs%d_withUE.root",outDir.Data(),cislo,jetR,massflag,ns);

    TFile *outFile=TFile::Open(outPath,"RECREATE");
    tEvent->Write();
    hDeltaPt->Write();
    hT_pT->Write();
    hT_eta->Write();
    hT_phi->Write();
    outFile->Write();
    outFile->Close();
    delete outFile;

    delete r3;
    delete f_pT;
    delete f_eta;
    delete f_phi;

    return 0;
}
