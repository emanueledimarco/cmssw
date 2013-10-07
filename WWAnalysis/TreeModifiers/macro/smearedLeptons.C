#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1F.h"
#include "scales2.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"

#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooFitResult.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace RooFit;

int ch = 1;
int eventmax = 1000;

void makeRecoPlots() {

    init(false);

    TH1F* hist_re_z1mass = new TH1F("hist_re_z1mass", "", 40, 40., 120.);
    TH1F* hist_re_z2mass = new TH1F("hist_re_z2mass", "", 50, 10., 110.);
    TH1F* hist_re_mass   = new TH1F("hist_re_mass"  , "", 50, 100., 150.);
    TH1F* hist_re_resp   = new TH1F("hist_re_resp"  , "", 40, 0.6, 1.4);

    TFile* file  = new TFile("/home/avartak/CMS/Higgs/Moriond/CMSSW_5_3_3_patch3/src/WWAnalysis/AnalysisStep/trees/hzzTree.root");
    TTree* tree  = (TTree*)file->Get("zz4lTree/probe_tree");

    Float_t channel = 0.;
    Float_t mass    = 0.;
    Float_t z1mass  = 0.;
    Float_t z2mass  = 0.;
    Float_t numsim  = 0.;
    Float_t mcl1pt  = 0.;
    Float_t mcl2pt  = 0.;
    Float_t mcl3pt  = 0.;
    Float_t mcl4pt  = 0.;
    Float_t mcl1eta = 0.;
    Float_t mcl2eta = 0.;
    Float_t mcl3eta = 0.;
    Float_t mcl4eta = 0.;
    Float_t l1pt    = 0.;
    Float_t l2pt    = 0.;
    Float_t l3pt    = 0.;
    Float_t l4pt    = 0.;
    Float_t l1eta   = 0.;
    Float_t l2eta   = 0.;
    Float_t l3eta   = 0.;
    Float_t l4eta   = 0.;

    TBranch *bchannel = tree->GetBranch("channel");
    TBranch *bmass    = tree->GetBranch("mass");
    TBranch *bz1mass  = tree->GetBranch("z1mass");
    TBranch *bz2mass  = tree->GetBranch("z2mass");
    TBranch *bnumsim  = tree->GetBranch("numTrueInteractions");
    TBranch *bmcl1pt  = tree->GetBranch("mc_l1pt");
    TBranch *bmcl2pt  = tree->GetBranch("mc_l2pt");
    TBranch *bmcl3pt  = tree->GetBranch("mc_l3pt");
    TBranch *bmcl4pt  = tree->GetBranch("mc_l4pt");
    TBranch *bmcl1eta = tree->GetBranch("mc_l1eta");
    TBranch *bmcl2eta = tree->GetBranch("mc_l2eta");
    TBranch *bmcl3eta = tree->GetBranch("mc_l3eta");
    TBranch *bmcl4eta = tree->GetBranch("mc_l4eta");
    TBranch *bl1pt    = tree->GetBranch("l1pt");
    TBranch *bl2pt    = tree->GetBranch("l2pt");
    TBranch *bl3pt    = tree->GetBranch("l3pt");
    TBranch *bl4pt    = tree->GetBranch("l4pt");
    TBranch *bl1eta   = tree->GetBranch("l1eta");
    TBranch *bl2eta   = tree->GetBranch("l2eta");
    TBranch *bl3eta   = tree->GetBranch("l3eta");
    TBranch *bl4eta   = tree->GetBranch("l4eta");

    bchannel->SetAddress(&channel);
    bmass->SetAddress(&mass);
    bz1mass->SetAddress(&z1mass);
    bz2mass->SetAddress(&z2mass);
    bmass->SetAddress(&mass);
    bnumsim->SetAddress (&numsim);
    bmcl1pt->SetAddress (&mcl1pt);
    bmcl2pt->SetAddress (&mcl2pt);
    bmcl3pt->SetAddress (&mcl3pt);
    bmcl4pt->SetAddress (&mcl4pt);
    bmcl1eta->SetAddress(&mcl1eta);
    bmcl2eta->SetAddress(&mcl2eta);
    bmcl3eta->SetAddress(&mcl3eta);
    bmcl4eta->SetAddress(&mcl4eta);
    bl1pt->SetAddress (&l1pt);
    bl2pt->SetAddress (&l2pt);
    bl3pt->SetAddress (&l3pt);
    bl4pt->SetAddress (&l4pt);
    bl1eta->SetAddress(&l1eta);
    bl2eta->SetAddress(&l2eta);
    bl3eta->SetAddress(&l3eta);
    bl4eta->SetAddress(&l4eta);

    int eventcount = 0;

    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; i++){
        bchannel->GetEvent(i);
        bmass->GetEvent(i);
        bz1mass->GetEvent(i);
        bz2mass->GetEvent(i);
        bnumsim->GetEvent(i);
        bmcl1pt->GetEvent(i);
        bmcl2pt->GetEvent(i);
        bmcl3pt->GetEvent(i);
        bmcl4pt->GetEvent(i);
        bmcl1eta->GetEvent(i);
        bmcl2eta->GetEvent(i);
        bmcl3eta->GetEvent(i);
        bmcl4eta->GetEvent(i);
        bl1pt->GetEvent(i);
        bl2pt->GetEvent(i);
        bl3pt->GetEvent(i);
        bl4pt->GetEvent(i);
        bl1eta->GetEvent(i);
        bl2eta->GetEvent(i);
        bl3eta->GetEvent(i);
        bl4eta->GetEvent(i);

        if (channel != ch) continue;

        float wgt = getPUWeight(numsim);

        if (mcl1pt == 0 || mcl2pt == 0 || mcl3pt == 0 || mcl4pt == 0) continue;

        eventcount++;
        if (eventcount>eventmax) break;

        if (z1mass>40 && z1mass<120 && z2mass>12 && z2mass<120) {
            hist_re_z1mass->Fill(z1mass, wgt);
            hist_re_z2mass->Fill(z2mass, wgt);
            hist_re_mass  ->Fill(mass  , wgt);
            hist_re_resp  ->Fill((l1pt*cosh(l1eta))/(mcl1pt*cosh(mcl1eta))  , wgt);
            hist_re_resp  ->Fill((l2pt*cosh(l2eta))/(mcl2pt*cosh(mcl2eta))  , wgt);
            hist_re_resp  ->Fill((l3pt*cosh(l3eta))/(mcl3pt*cosh(mcl3eta))  , wgt);
            hist_re_resp  ->Fill((l4pt*cosh(l4eta))/(mcl2pt*cosh(mcl2eta))  , wgt);
        }

    }

    TFile* outfile = new TFile("validationReco.root", "RECREATE");
    hist_re_z1mass->Write();
    hist_re_z2mass->Write();
    hist_re_mass->Write();
    hist_re_resp->Write();
    outfile->Close();


}




void makeGenPlots() {

    init(false);

    TH1F* hist_mc_z1mass = new TH1F("hist_mc_z1mass", "", 40, 40., 120.);
    TH1F* hist_mc_z2mass = new TH1F("hist_mc_z2mass", "", 50, 10., 110.);
    TH1F* hist_mc_mass   = new TH1F("hist_mc_mass"  , "", 50, 100., 150.);
    TH1F* hist_mc_resp   = new TH1F("hist_mc_resp"  , "", 40, 0.6, 1.4);

    TFile* file  = new TFile("/home/avartak/CMS/Higgs/Moriond/CMSSW_5_3_3_patch3/src/WWAnalysis/AnalysisStep/trees/gentree.root");
    //TFile* file  = new TFile("/home/avartak/CMS/Higgs/Moriond/CMSSW_5_3_3_patch3/src/WWAnalysis/AnalysisStep/trees/hzzTree.root");
    TTree* tree  = (TTree*)file->Get("zz4lTree/probe_tree");

    Float_t channel = 0.;
    Float_t mass    = 0.;
    Float_t z1mass  = 0.;
    Float_t z2mass  = 0.;
    Float_t mcl1pt  = 0.;
    Float_t mcl2pt  = 0.;
    Float_t mcl3pt  = 0.;
    Float_t mcl4pt  = 0.;
    Float_t mcl1eta = 0.;
    Float_t mcl2eta = 0.;
    Float_t mcl3eta = 0.;
    Float_t mcl4eta = 0.;
    Float_t mcl1phi = 0.;
    Float_t mcl2phi = 0.;
    Float_t mcl3phi = 0.;
    Float_t mcl4phi = 0.;
    Float_t numsim  = 0.;

    TBranch *bchannel = tree->GetBranch("channel");
    TBranch *bmass    = tree->GetBranch("mass");
    TBranch *bz1mass  = tree->GetBranch("z1mass");
    TBranch *bz2mass  = tree->GetBranch("z2mass");
    TBranch *bmcl1pt  = tree->GetBranch("mc_l1pt");
    TBranch *bmcl2pt  = tree->GetBranch("mc_l2pt");
    TBranch *bmcl3pt  = tree->GetBranch("mc_l3pt");
    TBranch *bmcl4pt  = tree->GetBranch("mc_l4pt");
    TBranch *bmcl1eta = tree->GetBranch("mc_l1eta");
    TBranch *bmcl2eta = tree->GetBranch("mc_l2eta");
    TBranch *bmcl3eta = tree->GetBranch("mc_l3eta");
    TBranch *bmcl4eta = tree->GetBranch("mc_l4eta");
    TBranch *bmcl1phi = tree->GetBranch("mc_l1phi");
    TBranch *bmcl2phi = tree->GetBranch("mc_l2phi");
    TBranch *bmcl3phi = tree->GetBranch("mc_l3phi");
    TBranch *bmcl4phi = tree->GetBranch("mc_l4phi");
    TBranch *bnumsim = tree->GetBranch("numTrueInteractions");

    bchannel->SetAddress(&channel);
    bmass->SetAddress(&mass);
    bz1mass->SetAddress(&z1mass);
    bz2mass->SetAddress(&z2mass);
    bmass->SetAddress(&mass);
    bmcl1pt->SetAddress(&mcl1pt);
    bmcl2pt->SetAddress(&mcl2pt);
    bmcl3pt->SetAddress(&mcl3pt);
    bmcl4pt->SetAddress(&mcl4pt);
    bmcl1eta->SetAddress(&mcl1eta);
    bmcl2eta->SetAddress(&mcl2eta);
    bmcl3eta->SetAddress(&mcl3eta);
    bmcl4eta->SetAddress(&mcl4eta);
    bmcl1phi->SetAddress(&mcl1phi);
    bmcl2phi->SetAddress(&mcl2phi);
    bmcl3phi->SetAddress(&mcl3phi);
    bmcl4phi->SetAddress(&mcl4phi);
    bnumsim->SetAddress(&numsim);


    TFile *elshapefile = new TFile("electronResponse.root");
    TFile *mushapefile = new TFile("muonResponse.root");

    TH2F* hist_mu_mean  = (TH2F*)mushapefile->Get("hist_mean"); 
    TH2F* hist_mu_sigma = (TH2F*)mushapefile->Get("hist_sigma"); 
    TH2F* hist_mu_aL    = (TH2F*)mushapefile->Get("hist_aL"); 
    TH2F* hist_mu_aR    = (TH2F*)mushapefile->Get("hist_aR"); 
    TH2F* hist_mu_nL    = (TH2F*)mushapefile->Get("hist_nL"); 
    TH2F* hist_mu_nR    = (TH2F*)mushapefile->Get("hist_nR"); 

    hist_mu_mean ->SetName("hist_mu_mean");
    hist_mu_sigma->SetName("hist_mu_sigma");
    hist_mu_aL   ->SetName("hist_mu_aL");
    hist_mu_aR   ->SetName("hist_mu_aR");
    hist_mu_nL   ->SetName("hist_mu_nL");
    hist_mu_nR   ->SetName("hist_mu_nR");

    TH2F* hist_el_mean  = (TH2F*)elshapefile->Get("hist_mean"); 
    TH2F* hist_el_sigma = (TH2F*)elshapefile->Get("hist_sigma"); 
    TH2F* hist_el_aL    = (TH2F*)elshapefile->Get("hist_aL"); 
    TH2F* hist_el_aR    = (TH2F*)elshapefile->Get("hist_aR"); 
    TH2F* hist_el_nL    = (TH2F*)elshapefile->Get("hist_nL"); 
    TH2F* hist_el_nR    = (TH2F*)elshapefile->Get("hist_nR"); 

    hist_el_mean ->SetName("hist_el_mean");
    hist_el_sigma->SetName("hist_el_sigma");
    hist_el_aL   ->SetName("hist_el_aL");
    hist_el_aR   ->SetName("hist_el_aR");
    hist_el_nL   ->SetName("hist_el_nL");
    hist_el_nR   ->SetName("hist_el_nR");

    int eventcount = 0;

    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; i++){
        bchannel->GetEvent(i);
        bmass->GetEvent(i);
        bz1mass->GetEvent(i);
        bz2mass->GetEvent(i);
        bmcl1pt->GetEvent(i);
        bmcl2pt->GetEvent(i);
        bmcl3pt->GetEvent(i);
        bmcl4pt->GetEvent(i);
        bmcl1eta->GetEvent(i);
        bmcl2eta->GetEvent(i);
        bmcl3eta->GetEvent(i);
        bmcl4eta->GetEvent(i);
        bmcl1phi->GetEvent(i);
        bmcl2phi->GetEvent(i);
        bmcl3phi->GetEvent(i);
        bmcl4phi->GetEvent(i);
        bnumsim->GetEvent(i);

        if (mcl1pt<5 || mcl2pt<5 || mcl3pt<5 || mcl4pt<5) continue;
        if (mcl1eta>2.5 || mcl2eta>2.5 || mcl3eta>2.5 || mcl4eta>2.5) continue;
        if (mcl1eta<-2.5 || mcl2eta<-2.5 || mcl3eta<-2.5 || mcl4eta<-2.5) continue;

        if (channel != ch) continue;        

        float wgt = getPUWeight(numsim);
        wgt = 1.0;

        float mean1  = 0.0;
        float sigma1 = 0.0;
        float aL1    = 0.0;
        float aR1    = 0.0;
        float nL1    = 0.0;
        float nR1    = 0.0;
        
        float mean2  = 0.0;
        float sigma2 = 0.0;
        float aL2    = 0.0;
        float aR2    = 0.0;
        float nL2    = 0.0;
        float nR2    = 0.0;
        
        float mean3  = 0.0;
        float sigma3 = 0.0;
        float aL3    = 0.0;
        float aR3    = 0.0;
        float nL3    = 0.0;
        float nR3    = 0.0;
        
        float mean4  = 0.0;
        float sigma4 = 0.0;
        float aL4    = 0.0;
        float aR4    = 0.0;
        float nL4    = 0.0;
        float nR4    = 0.0;

        if (channel == 0) {
            mean1  = hist_mu_mean ->GetBinContent(hist_mu_mean ->FindBin((float)min(mcl1pt, float(49.9)), fabs(mcl1eta)));
            sigma1 = hist_mu_sigma->GetBinContent(hist_mu_sigma->FindBin((float)min(mcl1pt, float(49.9)), fabs(mcl1eta)));
            aL1    = hist_mu_aL   ->GetBinContent(hist_mu_aL   ->FindBin((float)min(mcl1pt, float(49.9)), fabs(mcl1eta)));
            aR1    = hist_mu_aR   ->GetBinContent(hist_mu_aR   ->FindBin((float)min(mcl1pt, float(49.9)), fabs(mcl1eta)));
            nL1    = hist_mu_nL   ->GetBinContent(hist_mu_nL   ->FindBin((float)min(mcl1pt, float(49.9)), fabs(mcl1eta)));
            nR1    = hist_mu_nR   ->GetBinContent(hist_mu_nR   ->FindBin((float)min(mcl1pt, float(49.9)), fabs(mcl1eta)));

            mean2  = hist_mu_mean ->GetBinContent(hist_mu_mean ->FindBin((float)min(mcl2pt, float(49.9)), fabs(mcl2eta)));
            sigma2 = hist_mu_sigma->GetBinContent(hist_mu_sigma->FindBin((float)min(mcl2pt, float(49.9)), fabs(mcl2eta)));
            aL2    = hist_mu_aL   ->GetBinContent(hist_mu_aL   ->FindBin((float)min(mcl2pt, float(49.9)), fabs(mcl2eta)));
            aR2    = hist_mu_aR   ->GetBinContent(hist_mu_aR   ->FindBin((float)min(mcl2pt, float(49.9)), fabs(mcl2eta)));
            nL2    = hist_mu_nL   ->GetBinContent(hist_mu_nL   ->FindBin((float)min(mcl2pt, float(49.9)), fabs(mcl2eta)));
            nR2    = hist_mu_nR   ->GetBinContent(hist_mu_nR   ->FindBin((float)min(mcl2pt, float(49.9)), fabs(mcl2eta)));

            mean3  = hist_mu_mean ->GetBinContent(hist_mu_mean ->FindBin((float)min(mcl3pt, float(49.9)), fabs(mcl3eta)));
            sigma3 = hist_mu_sigma->GetBinContent(hist_mu_sigma->FindBin((float)min(mcl3pt, float(49.9)), fabs(mcl3eta)));
            aL3    = hist_mu_aL   ->GetBinContent(hist_mu_aL   ->FindBin((float)min(mcl3pt, float(49.9)), fabs(mcl3eta)));
            aR3    = hist_mu_aR   ->GetBinContent(hist_mu_aR   ->FindBin((float)min(mcl3pt, float(49.9)), fabs(mcl3eta)));
            nL3    = hist_mu_nL   ->GetBinContent(hist_mu_nL   ->FindBin((float)min(mcl3pt, float(49.9)), fabs(mcl3eta)));
            nR3    = hist_mu_nR   ->GetBinContent(hist_mu_nR   ->FindBin((float)min(mcl3pt, float(49.9)), fabs(mcl3eta)));

            mean4  = hist_mu_mean ->GetBinContent(hist_mu_mean ->FindBin((float)min(mcl4pt, float(49.9)), fabs(mcl4eta)));
            sigma4 = hist_mu_sigma->GetBinContent(hist_mu_sigma->FindBin((float)min(mcl4pt, float(49.9)), fabs(mcl4eta)));
            aL4    = hist_mu_aL   ->GetBinContent(hist_mu_aL   ->FindBin((float)min(mcl4pt, float(49.9)), fabs(mcl4eta)));
            aR4    = hist_mu_aR   ->GetBinContent(hist_mu_aR   ->FindBin((float)min(mcl4pt, float(49.9)), fabs(mcl4eta)));
            nL4    = hist_mu_nL   ->GetBinContent(hist_mu_nL   ->FindBin((float)min(mcl4pt, float(49.9)), fabs(mcl4eta)));
            nR4    = hist_mu_nR   ->GetBinContent(hist_mu_nR   ->FindBin((float)min(mcl4pt, float(49.9)), fabs(mcl4eta)));
        }
       
        else if (channel == 1) {
            mean1  = hist_el_mean ->GetBinContent(hist_el_mean ->FindBin((float)min(mcl1pt, float(99.9)), fabs(mcl1eta)));
            sigma1 = hist_el_sigma->GetBinContent(hist_el_sigma->FindBin((float)min(mcl1pt, float(99.9)), fabs(mcl1eta)));
            aL1    = hist_el_aL   ->GetBinContent(hist_el_aL   ->FindBin((float)min(mcl1pt, float(99.9)), fabs(mcl1eta)));
            aR1    = hist_el_aR   ->GetBinContent(hist_el_aR   ->FindBin((float)min(mcl1pt, float(99.9)), fabs(mcl1eta)));
            nL1    = hist_el_nL   ->GetBinContent(hist_el_nL   ->FindBin((float)min(mcl1pt, float(99.9)), fabs(mcl1eta)));
            nR1    = hist_el_nR   ->GetBinContent(hist_el_nR   ->FindBin((float)min(mcl1pt, float(99.9)), fabs(mcl1eta)));

            mean2  = hist_el_mean ->GetBinContent(hist_el_mean ->FindBin((float)min(mcl2pt, float(99.9)), fabs(mcl2eta)));
            sigma2 = hist_el_sigma->GetBinContent(hist_el_sigma->FindBin((float)min(mcl2pt, float(99.9)), fabs(mcl2eta)));
            aL2    = hist_el_aL   ->GetBinContent(hist_el_aL   ->FindBin((float)min(mcl2pt, float(99.9)), fabs(mcl2eta)));
            aR2    = hist_el_aR   ->GetBinContent(hist_el_aR   ->FindBin((float)min(mcl2pt, float(99.9)), fabs(mcl2eta)));
            nL2    = hist_el_nL   ->GetBinContent(hist_el_nL   ->FindBin((float)min(mcl2pt, float(99.9)), fabs(mcl2eta)));
            nR2    = hist_el_nR   ->GetBinContent(hist_el_nR   ->FindBin((float)min(mcl2pt, float(99.9)), fabs(mcl2eta)));

            mean3  = hist_el_mean ->GetBinContent(hist_el_mean ->FindBin((float)min(mcl3pt, float(99.9)), fabs(mcl3eta)));
            sigma3 = hist_el_sigma->GetBinContent(hist_el_sigma->FindBin((float)min(mcl3pt, float(99.9)), fabs(mcl3eta)));
            aL3    = hist_el_aL   ->GetBinContent(hist_el_aL   ->FindBin((float)min(mcl3pt, float(99.9)), fabs(mcl3eta)));
            aR3    = hist_el_aR   ->GetBinContent(hist_el_aR   ->FindBin((float)min(mcl3pt, float(99.9)), fabs(mcl3eta)));
            nL3    = hist_el_nL   ->GetBinContent(hist_el_nL   ->FindBin((float)min(mcl3pt, float(99.9)), fabs(mcl3eta)));
            nR3    = hist_el_nR   ->GetBinContent(hist_el_nR   ->FindBin((float)min(mcl3pt, float(99.9)), fabs(mcl3eta)));

            mean4  = hist_el_mean ->GetBinContent(hist_el_mean ->FindBin((float)min(mcl4pt, float(99.9)), fabs(mcl4eta)));
            sigma4 = hist_el_sigma->GetBinContent(hist_el_sigma->FindBin((float)min(mcl4pt, float(99.9)), fabs(mcl4eta)));
            aL4    = hist_el_aL   ->GetBinContent(hist_el_aL   ->FindBin((float)min(mcl4pt, float(99.9)), fabs(mcl4eta)));
            aR4    = hist_el_aR   ->GetBinContent(hist_el_aR   ->FindBin((float)min(mcl4pt, float(99.9)), fabs(mcl4eta)));
            nL4    = hist_el_nL   ->GetBinContent(hist_el_nL   ->FindBin((float)min(mcl4pt, float(99.9)), fabs(mcl4eta)));
            nR4    = hist_el_nR   ->GetBinContent(hist_el_nR   ->FindBin((float)min(mcl4pt, float(99.9)), fabs(mcl4eta)));
        }

        else if (channel == 2) {
            mean1  = hist_el_mean ->GetBinContent(hist_el_mean ->FindBin((float)min(mcl1pt, float(99.9)), fabs(mcl1eta)));
            sigma1 = hist_el_sigma->GetBinContent(hist_el_sigma->FindBin((float)min(mcl1pt, float(99.9)), fabs(mcl1eta)));
            aL1    = hist_el_aL   ->GetBinContent(hist_el_aL   ->FindBin((float)min(mcl1pt, float(99.9)), fabs(mcl1eta)));
            aR1    = hist_el_aR   ->GetBinContent(hist_el_aR   ->FindBin((float)min(mcl1pt, float(99.9)), fabs(mcl1eta)));
            nL1    = hist_el_nL   ->GetBinContent(hist_el_nL   ->FindBin((float)min(mcl1pt, float(99.9)), fabs(mcl1eta)));
            nR1    = hist_el_nR   ->GetBinContent(hist_el_nR   ->FindBin((float)min(mcl1pt, float(99.9)), fabs(mcl1eta)));

            mean2  = hist_el_mean ->GetBinContent(hist_el_mean ->FindBin((float)min(mcl2pt, float(99.9)), fabs(mcl2eta)));
            sigma2 = hist_el_sigma->GetBinContent(hist_el_sigma->FindBin((float)min(mcl2pt, float(99.9)), fabs(mcl2eta)));
            aL2    = hist_el_aL   ->GetBinContent(hist_el_aL   ->FindBin((float)min(mcl2pt, float(99.9)), fabs(mcl2eta)));
            aR2    = hist_el_aR   ->GetBinContent(hist_el_aR   ->FindBin((float)min(mcl2pt, float(99.9)), fabs(mcl2eta)));
            nL2    = hist_el_nL   ->GetBinContent(hist_el_nL   ->FindBin((float)min(mcl2pt, float(99.9)), fabs(mcl2eta)));
            nR2    = hist_el_nR   ->GetBinContent(hist_el_nR   ->FindBin((float)min(mcl2pt, float(99.9)), fabs(mcl2eta)));

            mean3  = hist_mu_mean ->GetBinContent(hist_mu_mean ->FindBin((float)min(mcl3pt, float(49.9)), fabs(mcl3eta)));
            sigma3 = hist_mu_sigma->GetBinContent(hist_mu_sigma->FindBin((float)min(mcl3pt, float(49.9)), fabs(mcl3eta)));
            aL3    = hist_mu_aL   ->GetBinContent(hist_mu_aL   ->FindBin((float)min(mcl3pt, float(49.9)), fabs(mcl3eta)));
            aR3    = hist_mu_aR   ->GetBinContent(hist_mu_aR   ->FindBin((float)min(mcl3pt, float(49.9)), fabs(mcl3eta)));
            nL3    = hist_mu_nL   ->GetBinContent(hist_mu_nL   ->FindBin((float)min(mcl3pt, float(49.9)), fabs(mcl3eta)));
            nR3    = hist_mu_nR   ->GetBinContent(hist_mu_nR   ->FindBin((float)min(mcl3pt, float(49.9)), fabs(mcl3eta)));

            mean4  = hist_mu_mean ->GetBinContent(hist_mu_mean ->FindBin((float)min(mcl4pt, float(49.9)), fabs(mcl4eta)));
            sigma4 = hist_mu_sigma->GetBinContent(hist_mu_sigma->FindBin((float)min(mcl4pt, float(49.9)), fabs(mcl4eta)));
            aL4    = hist_mu_aL   ->GetBinContent(hist_mu_aL   ->FindBin((float)min(mcl4pt, float(49.9)), fabs(mcl4eta)));
            aR4    = hist_mu_aR   ->GetBinContent(hist_mu_aR   ->FindBin((float)min(mcl4pt, float(49.9)), fabs(mcl4eta)));
            nL4    = hist_mu_nL   ->GetBinContent(hist_mu_nL   ->FindBin((float)min(mcl4pt, float(49.9)), fabs(mcl4eta)));
            nR4    = hist_mu_nR   ->GetBinContent(hist_mu_nR   ->FindBin((float)min(mcl4pt, float(49.9)), fabs(mcl4eta)));
        }

        else  {
            mean1  = hist_mu_mean ->GetBinContent(hist_mu_mean ->FindBin((float)min(mcl1pt, float(49.9)), fabs(mcl1eta)));
            sigma1 = hist_mu_sigma->GetBinContent(hist_mu_sigma->FindBin((float)min(mcl1pt, float(49.9)), fabs(mcl1eta)));
            aL1    = hist_mu_aL   ->GetBinContent(hist_mu_aL   ->FindBin((float)min(mcl1pt, float(49.9)), fabs(mcl1eta)));
            aR1    = hist_mu_aR   ->GetBinContent(hist_mu_aR   ->FindBin((float)min(mcl1pt, float(49.9)), fabs(mcl1eta)));
            nL1    = hist_mu_nL   ->GetBinContent(hist_mu_nL   ->FindBin((float)min(mcl1pt, float(49.9)), fabs(mcl1eta)));
            nR1    = hist_mu_nR   ->GetBinContent(hist_mu_nR   ->FindBin((float)min(mcl1pt, float(49.9)), fabs(mcl1eta)));

            mean2  = hist_mu_mean ->GetBinContent(hist_mu_mean ->FindBin((float)min(mcl2pt, float(49.9)), fabs(mcl2eta)));
            sigma2 = hist_mu_sigma->GetBinContent(hist_mu_sigma->FindBin((float)min(mcl2pt, float(49.9)), fabs(mcl2eta)));
            aL2    = hist_mu_aL   ->GetBinContent(hist_mu_aL   ->FindBin((float)min(mcl2pt, float(49.9)), fabs(mcl2eta)));
            aR2    = hist_mu_aR   ->GetBinContent(hist_mu_aR   ->FindBin((float)min(mcl2pt, float(49.9)), fabs(mcl2eta)));
            nL2    = hist_mu_nL   ->GetBinContent(hist_mu_nL   ->FindBin((float)min(mcl2pt, float(49.9)), fabs(mcl2eta)));
            nR2    = hist_mu_nR   ->GetBinContent(hist_mu_nR   ->FindBin((float)min(mcl2pt, float(49.9)), fabs(mcl2eta)));

            mean3  = hist_el_mean ->GetBinContent(hist_el_mean ->FindBin((float)min(mcl3pt, float(99.9)), fabs(mcl3eta)));
            sigma3 = hist_el_sigma->GetBinContent(hist_el_sigma->FindBin((float)min(mcl3pt, float(99.9)), fabs(mcl3eta)));
            aL3    = hist_el_aL   ->GetBinContent(hist_el_aL   ->FindBin((float)min(mcl3pt, float(99.9)), fabs(mcl3eta)));
            aR3    = hist_el_aR   ->GetBinContent(hist_el_aR   ->FindBin((float)min(mcl3pt, float(99.9)), fabs(mcl3eta)));
            nL3    = hist_el_nL   ->GetBinContent(hist_el_nL   ->FindBin((float)min(mcl3pt, float(99.9)), fabs(mcl3eta)));
            nR3    = hist_el_nR   ->GetBinContent(hist_el_nR   ->FindBin((float)min(mcl3pt, float(99.9)), fabs(mcl3eta)));

            mean4  = hist_el_mean ->GetBinContent(hist_el_mean ->FindBin((float)min(mcl4pt, float(99.9)), fabs(mcl4eta)));
            sigma4 = hist_el_sigma->GetBinContent(hist_el_sigma->FindBin((float)min(mcl4pt, float(99.9)), fabs(mcl4eta)));
            aL4    = hist_el_aL   ->GetBinContent(hist_el_aL   ->FindBin((float)min(mcl4pt, float(99.9)), fabs(mcl4eta)));
            aR4    = hist_el_aR   ->GetBinContent(hist_el_aR   ->FindBin((float)min(mcl4pt, float(99.9)), fabs(mcl4eta)));
            nL4    = hist_el_nL   ->GetBinContent(hist_el_nL   ->FindBin((float)min(mcl4pt, float(99.9)), fabs(mcl4eta)));
            nR4    = hist_el_nR   ->GetBinContent(hist_el_nR   ->FindBin((float)min(mcl4pt, float(99.9)), fabs(mcl4eta)));
        }

        /*
        std::cout << mcl1pt << " " << mcl1eta <<  std::endl;
        std::cout << mcl2pt << " " << mcl2eta <<  std::endl;
        std::cout << mcl3pt << " " << mcl3eta <<  std::endl;
        std::cout << mcl4pt << " " << mcl4eta << std::endl;
        */

        RooRealVar res1("res1", "res1", 1, 0, 10);
        RooRealVar res2("res2", "res2", 1, 0, 10);
        RooRealVar res3("res3", "res3", 1, 0, 10);
        RooRealVar res4("res4", "res4", 1, 0, 10);

        RooRealVar meanvar1("meanvar1", "", mean1);
        RooRealVar meanvar2("meanvar2", "", mean2);
        RooRealVar meanvar3("meanvar3", "", mean3);
        RooRealVar meanvar4("meanvar4", "", mean4);

        RooRealVar sigmavar1("sigmavar1", "", sigma1);
        RooRealVar sigmavar2("sigmavar2", "", sigma2);
        RooRealVar sigmavar3("sigmavar3", "", sigma3);
        RooRealVar sigmavar4("sigmavar4", "", sigma4);

        RooRealVar aLvar1("aLvar1", "", aL1);
        RooRealVar aLvar2("aLvar2", "", aL2);
        RooRealVar aLvar3("aLvar3", "", aL3);
        RooRealVar aLvar4("aLvar4", "", aL4);

        RooRealVar aRvar1("aRvar1", "", aR1);
        RooRealVar aRvar2("aRvar2", "", aR2);
        RooRealVar aRvar3("aRvar3", "", aR3);
        RooRealVar aRvar4("aRvar4", "", aR4);

        RooRealVar nLvar1("nLvar1", "", nL1);
        RooRealVar nLvar2("nLvar2", "", nL2);
        RooRealVar nLvar3("nLvar3", "", nL3);
        RooRealVar nLvar4("nLvar4", "", nL4);

        RooRealVar nRvar1("nRvar1", "", nR1);
        RooRealVar nRvar2("nRvar2", "", nR2);
        RooRealVar nRvar3("nRvar3", "", nR3);
        RooRealVar nRvar4("nRvar4", "", nR4);

        RooDoubleCB pdf1("pdf1", "pdf1", res1, meanvar1, sigmavar1, aLvar1, nLvar1, aRvar1, nRvar1);
        RooDoubleCB pdf2("pdf2", "pdf2", res2, meanvar2, sigmavar2, aLvar2, nLvar2, aRvar2, nRvar2);
        RooDoubleCB pdf3("pdf3", "pdf3", res3, meanvar3, sigmavar3, aLvar3, nLvar3, aRvar3, nRvar3);
        RooDoubleCB pdf4("pdf4", "pdf4", res4, meanvar4, sigmavar4, aLvar4, nLvar4, aRvar4, nRvar4);

        RooDataSet* dset1 = pdf1.generate(res1, 1);
        RooDataSet* dset2 = pdf2.generate(res2, 1);
        RooDataSet* dset3 = pdf3.generate(res3, 1);
        RooDataSet* dset4 = pdf4.generate(res4, 1);

        float ran1 = 0.0;
        float ran2 = 0.0;
        float ran3 = 0.0;
        float ran4 = 0.0;

        ran1 = dset1->get(0)->getRealValue("res1");
        ran2 = dset2->get(0)->getRealValue("res2");
        ran3 = dset3->get(0)->getRealValue("res3");
        ran4 = dset4->get(0)->getRealValue("res4");

        /*
        std::cout << mean1 << " " <<  sigma1 << " " <<  aL1 << " " <<  aR1 << " " <<  nL1 << " " <<  nR1 << " " << ran1 << " " <<  std::endl;
        std::cout << mean2 << " " <<  sigma2 << " " <<  aL2 << " " <<  aR2 << " " <<  nL2 << " " <<  nR2 << " " << ran2 << " " <<  std::endl;
        std::cout << mean3 << " " <<  sigma3 << " " <<  aL3 << " " <<  aR3 << " " <<  nL3 << " " <<  nR3 << " " << ran3 << " " <<  std::endl;
        std::cout << mean4 << " " <<  sigma4 << " " <<  aL4 << " " <<  aR4 << " " <<  nL4 << " " <<  nR4 << " " << ran4 << " " <<  std::endl << std::endl;
        */

        float l1px = mcl1pt*cos(mcl1phi);
        float l2px = mcl2pt*cos(mcl2phi);
        float l3px = mcl3pt*cos(mcl3phi);
        float l4px = mcl4pt*cos(mcl4phi);

        float l1py = mcl1pt*sin(mcl1phi);
        float l2py = mcl2pt*sin(mcl2phi);
        float l3py = mcl3pt*sin(mcl3phi);
        float l4py = mcl4pt*sin(mcl4phi);

        float l1pz = mcl1pt*sinh(mcl1eta);
        float l2pz = mcl2pt*sinh(mcl2eta);
        float l3pz = mcl3pt*sinh(mcl3eta);
        float l4pz = mcl4pt*sinh(mcl4eta);

        float l1p  = mcl1pt*cosh(mcl1eta);
        float l2p  = mcl2pt*cosh(mcl2eta);
        float l3p  = mcl3pt*cosh(mcl3eta);
        float l4p  = mcl4pt*cosh(mcl4eta);

        l1px *= ran1;
        l1py *= ran1;
        l1pz *= ran1;
        l1p  *= ran1;

        l2px *= ran2;
        l2py *= ran2;
        l2pz *= ran2;
        l2p  *= ran2;

        l3px *= ran3;
        l3py *= ran3;
        l3pz *= ran3;
        l3p  *= ran3;

        l4px *= ran4;
        l4py *= ran4;
        l4pz *= ran4;
        l4p  *= ran4;

        float sl1pt = sqrt(l1px*l1px + l1py*l1py);
        float sl2pt = sqrt(l2px*l2px + l2py*l2py);
        float sl3pt = sqrt(l3px*l3px + l3py*l3py);
        float sl4pt = sqrt(l4px*l4px + l4py*l4py);

        float z1px = l1px + l2px;
        float z1py = l1py + l2py;
        float z1pz = l1pz + l2pz;
        float z1p  = l1p  + l2p;

        float z2px = l3px + l4px;
        float z2py = l3py + l4py;
        float z2pz = l3pz + l4pz;
        float z2p  = l3p  + l4p;

        float hpx = z1px + z2px;
        float hpy = z1py + z2py;
        float hpz = z1pz + z2pz;
        float hp  = z1p  + z2p;

        float z1m  = sqrt(z1p*z1p - (z1px*z1px + z1py*z1py + z1pz*z1pz));
        float z2m  = sqrt(z2p*z2p - (z2px*z2px + z2py*z2py + z2pz*z2pz));
        float hm   = sqrt(hp*hp   - (hpx*hpx   + hpy*hpy   + hpz*hpz));

        if (channel == 0 && (sl1pt<5 || sl2pt<5 || sl3pt<5 || sl4pt<5 || fabs(mcl1eta)>2.4 || fabs(mcl2eta)>2.4 || fabs(mcl3eta)>2.4 || fabs(mcl4eta)>2.4 || z1mass<40 || z1mass>120 || z2mass<12 || z2mass>120)) continue;
        if (channel == 1 && (sl1pt<7 || sl2pt<7 || sl3pt<7 || sl4pt<7 || fabs(mcl1eta)>2.5 || fabs(mcl2eta)>2.5 || fabs(mcl3eta)>2.5 || fabs(mcl4eta)>2.5 || z1mass<40 || z1mass>120 || z2mass<12 || z2mass>120)) continue;
        if (channel == 2 && (sl1pt<7 || sl2pt<7 || sl3pt<5 || sl4pt<5 || fabs(mcl1eta)>2.5 || fabs(mcl2eta)>2.5 || fabs(mcl3eta)>2.4 || fabs(mcl4eta)>2.4 || z1mass<40 || z1mass>120 || z2mass<12 || z2mass>120)) continue;
        if (channel == 3 && (sl1pt<5 || sl2pt<5 || sl3pt<7 || sl4pt<7 || fabs(mcl1eta)>2.4 || fabs(mcl2eta)>2.4 || fabs(mcl3eta)>2.5 || fabs(mcl4eta)>2.5 || z1mass<40 || z1mass>120 || z2mass<12 || z2mass>120)) continue;

        float wgt1 = 1.0;
        float wgt2 = 1.0;
        float wgt3 = 1.0;
        float wgt4 = 1.0;

        if (channel == 0) {
            wgt1 *= getPR(sl1pt, mcl1eta, 13.);
            wgt2 *= getPR(sl2pt, mcl2eta, 13.);
            wgt3 *= getPR(sl3pt, mcl3eta, 13.);
            wgt4 *= getPR(sl4pt, mcl4eta, 13.);
        }

        else if (channel == 1) {
            wgt1 *= getPR(sl1pt, mcl1eta, 11.);
            wgt2 *= getPR(sl2pt, mcl2eta, 11.);
            wgt3 *= getPR(sl3pt, mcl3eta, 11.);
            wgt4 *= getPR(sl4pt, mcl4eta, 11.);
        }
   
        else if (channel == 2) {
            wgt1 *= getPR(sl1pt, mcl1eta, 11.);
            wgt2 *= getPR(sl2pt, mcl2eta, 11.);
            wgt3 *= getPR(sl3pt, mcl3eta, 13.);
            wgt4 *= getPR(sl4pt, mcl4eta, 13.);
        }

        else {
            wgt1 *= getPR(sl1pt, mcl1eta, 13.);
            wgt2 *= getPR(sl2pt, mcl2eta, 13.);
            wgt3 *= getPR(sl3pt, mcl3eta, 11.);
            wgt4 *= getPR(sl4pt, mcl4eta, 11.);
        }

        wgt *= wgt1;
        wgt *= wgt2;
        wgt *= wgt3;
        wgt *= wgt4;

        hist_mc_z1mass->Fill(z1m, wgt);
        hist_mc_z2mass->Fill(z2m, wgt);
        hist_mc_mass  ->Fill(hm , wgt);
        hist_mc_resp->Fill(ran1, wgt1);
        hist_mc_resp->Fill(ran2, wgt2);
        hist_mc_resp->Fill(ran3, wgt3);
        hist_mc_resp->Fill(ran4, wgt4);


        eventcount++;
        if (eventcount>eventmax) break;

    }

    TFile* outfile = new TFile("validationGen.root", "RECREATE");
    hist_mc_z1mass->Write();
    hist_mc_z2mass->Write();
    hist_mc_mass->Write();
    hist_mc_resp->Write();
    outfile->Close();

}

void smearedLeptons() {

    makeRecoPlots();
    //makeGenPlots();

}
