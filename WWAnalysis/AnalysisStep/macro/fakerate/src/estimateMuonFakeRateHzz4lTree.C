#define estimateMuonFakeRateHzz4lTree_cxx
#include "include/estimateMuonFakeRateHzz4lTree.h"
#include "include/EfficiencyEvaluator.hh"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <iostream>

using namespace std;

Bool_t estimateMuonFakeRateHzz4lTree::passRefMuSel() {
  
  float abseta=fabs(eta);

  MuonEffectiveArea::MuonEffectiveAreaTarget effAreaTarget_ = MuonEffectiveArea::kMuEAData2011;
  MuonEffectiveArea::MuonEffectiveAreaType effAreaGamma_   = MuonEffectiveArea::kMuGammaIso04;
  MuonEffectiveArea::MuonEffectiveAreaType effAreaNeutralHad_ = MuonEffectiveArea::kMuNeutralHadronIso04;
  MuonEffectiveArea::MuonEffectiveAreaType effAreaGammaAndNeutralHad_ =  MuonEffectiveArea::kMuGammaAndNeutralHadronIso04;

  float eff_area_ga  = MuonEffectiveArea::GetMuonEffectiveArea(effAreaGamma_, abseta, effAreaTarget_);
  float eff_area_nh  = MuonEffectiveArea::GetMuonEffectiveArea(effAreaNeutralHad_, abseta, effAreaTarget_);
  float eff_area_ganh = MuonEffectiveArea::GetMuonEffectiveArea(effAreaGammaAndNeutralHad_, abseta, effAreaTarget_);

  float iso = pfIsoChHad04;
  iso += max<float>(0.,pfIsoNHad04_NoEA+pfIsoPhoton04_NoEA - eff_area_ganh*rhoAA);

  return (pfid && iso/pt<0.40 && fabs(sip3d)<4.0);
} 

Bool_t estimateMuonFakeRateHzz4lTree::passMvaMuSel(bool pfnopuiso) {
  
  // id part
  if (!pfid || fabs(sip3d)>=4.0) return false;
  // iso part
  float abseta=fabs(eta);
  float bdtiso = (pfnopuiso) ? bdtIsoPnp : bdtIsoDz;
  
  if(globalmu && trackermu) {
    if(pt<10 && abseta<1.5) return bdtiso>-0.593;
    if(pt>=10 && abseta<1.5) return bdtiso>0.337;
    if(pt<10 && abseta>=1.5) return bdtiso>-0.767;
    if(pt>=10 && abseta>=1.5) return bdtiso>0.410;
  }
  if(!globalmu && trackermu) return bdtiso>-0.989;
  if(globalmu && !trackermu) return bdtiso>-0.995;
  std::cout << "This should never happen! passMvaMuSel() returning false." << std::endl;
  return false;
} 

void estimateMuonFakeRateHzz4lTree::Loop(const char *outname)
{
  if (fChain == 0) return;

  std::vector<TString> hmuID;
  hmuID.push_back("hzzPfIso");  // reference pf iso
  hmuID.push_back("hzzMvaPfIso");  // duncan's optimization for HZZ

  // -----------------------------------------------------------------------
  // study vs eta

  Float_t LowerEta[5];
  LowerEta[0]=0.0;
  LowerEta[1]=1.0;
  LowerEta[2]=1.479;  
  LowerEta[3]=2.0;
  LowerEta[4]=2.5;
  TH1F *RecoEtaHighPt   = new TH1F( "RecoEtaHighPt",    "reconstructed #eta", 4, LowerEta);
  TH1F *RecoEtaLowPt    = new TH1F( "RecoEtaLowPt",     "reconstructed #eta", 4, LowerEta);
  
  // eta, high pT
  std::vector<TH1F*> hmuEtaHighPt;
  for (int i=0;i<(int)hmuID.size();++i) {
    TH1F* aHisto = new TH1F( "hmu"+TString(hmuID[i])+"EtaHighPt",   "HZZ BDT ID #eta", 4, LowerEta);     hmuEtaHighPt.push_back(aHisto);
  }

  // eta, low pT
  std::vector<TH1F*> hmuEtaLowPt;
  for (int i=0;i<(int)hmuID.size();++i) {
    TH1F* aHisto = new TH1F( "hmu"+TString(hmuID[i])+"EtaLowPt",   "HZZ BDT ID #eta", 4, LowerEta);
    hmuEtaLowPt.push_back(aHisto);
  }

  // -----------------------------------------------------------------------
  // study vs pT
  Float_t LowerPt[11] = {0.0,5.0,7.0,10.0,15.0,20.0,25.0,30.0,40.0,50.0,80.0};
  TH1F *RecoPtBarrel   = new TH1F( "RecoPtBarrel",    "reconstructed p_{T} (GeV)", 10, LowerPt);
  TH1F *RecoPtEndcap   = new TH1F( "RecoPtEndcap",    "reconstructed p_{T} (GeV)", 10, LowerPt);

  // to have the full picture in the barrel
  std::vector<TH1F*> hmuPtBarrel;
  for (int i=0;i<(int)hmuID.size();++i) {
    TH1F* aHisto = new TH1F( "hmu"+TString(hmuID[i])+"PtBarrel", "HZZ BDT ID #eta",   10, LowerPt );
    hmuPtBarrel.push_back(aHisto);
  }

  // to have the full picture in the endcap
  std::vector<TH1F*> hmuPtEndcap;
  for (int i=0;i<(int)hmuID.size();++i) {
    TH1F* aHisto = new TH1F( "hmu"+TString(hmuID[i])+"PtEndcap", "HZZ BDT ID #eta",   10, LowerPt);
    hmuPtEndcap.push_back(aHisto);
  }

  // -----------------------------------------------------------------------
  // study vs PU
  Float_t LowerPU[11];
  LowerPU[0] = 1;
  LowerPU[1] = 4;
  for(int i=2;i<6;++i) LowerPU[i]=i+3;
  LowerPU[6] = 10;
  LowerPU[7] = 15;
  LowerPU[8] = 25;
  LowerPU[9] = 30;
  LowerPU[10] = 35;

  TH1F *RecoPUBarrel   = new TH1F( "RecoPUBarrel",   "reconstructed nPU", 10, LowerPU);
  TH1F *RecoPUEndcap   = new TH1F( "RecoPUEndcap",   "reconstructed nPU", 10, LowerPU);

  //  barrel
  std::vector<TH1F*> hmuPUBarrel;
  for (int i=0;i<(int)hmuID.size();++i) {
    TH1F* aHisto = new TH1F( "hmu"+TString(hmuID[i])+"PUBarrel", "HZZ BDT ID #eta",   10, LowerPU );
    hmuPUBarrel.push_back(aHisto);
  }

  // endcap
  std::vector<TH1F*> hmuPUEndcap;
  for (int i=0;i<(int)hmuID.size();++i) {
    TH1F* aHisto = new TH1F( "hmu"+TString(hmuID[i])+"PUEndcap", "HZZ BDT ID #eta",   10, LowerPU );
    hmuPUEndcap.push_back(aHisto);
  }

  // loop on events
  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%100000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    // fill the denominator: take only the highest pT denominator candidate
    float etaFake = fabs(eta);
    float etFake  = pt;
    bool isInEB   = fabs(eta)<1.479;
    bool isInEE   = !isInEB;
    bool highPt   = (pt>10.);
    bool lowPt    = (pt<=10.);

    // pass the event selection: strict Z mass and fakeable lepton is required to have deltaR > 0.8 w.r.t. the Z1 leptons
    TVector3 l1p,l2p,lp;
    l1p.SetPtEtaPhi(l1pt,l1eta,l1phi);
    l2p.SetPtEtaPhi(l2pt,l2eta,l2phi);
    lp.SetPtEtaPhi(pt,eta,phi);
    if(fabs(zmass-91.1876)>10 || lp.DeltaR(l1p)<0.8 || lp.DeltaR(l2p)<0.8) continue;

    // pass the denominator object
    if(etFake<5 || etaFake>2.4 || fabs(sip3d)>4) continue;

    // filling
    if (highPt) RecoEtaHighPt -> Fill(etaFake);  //, theWeight); 
    if (lowPt)  RecoEtaLowPt  -> Fill(etaFake);  // , theWeight);

    if (isInEB) {
      RecoPtBarrel -> Fill(etFake);  //, theWeight);
      RecoPUBarrel -> Fill(numvertices);  //, theWeight);
    }
    if (isInEE) {
      RecoPtEndcap -> Fill(etFake);  //, theWeight);
      RecoPUEndcap -> Fill(numvertices);  //, theWeight);
    }


    // fill the numerator(s)
    // === PF-muon ID, combined Iso ===
    if(passRefMuSel()) {
      if (highPt) hmuEtaHighPt[khzzPfIso]->Fill(etaFake);
      if (lowPt)  hmuEtaLowPt[khzzPfIso] ->Fill(etaFake);
      if (isInEB) { 
	hmuPtBarrel[khzzPfIso] ->Fill(etFake);
	hmuPUBarrel[khzzPfIso] ->Fill(numvertices);
      }
      if (isInEE) {
	hmuPtEndcap[khzzPfIso] ->Fill(etFake);
	hmuPUEndcap[khzzPfIso] ->Fill(numvertices);
      }
    }

    // === PF-muon ID, MVA Iso ===
    if(passMvaMuSel(true)) {
      if (highPt) hmuEtaHighPt[khzzMvaPfIso]->Fill(etaFake);
      if (lowPt)  hmuEtaLowPt[khzzMvaPfIso] ->Fill(etaFake);
      if (isInEB) { 
	hmuPtBarrel[khzzMvaPfIso] ->Fill(etFake);
	hmuPUBarrel[khzzMvaPfIso] ->Fill(numvertices);
      }
      if (isInEE) {
	hmuPtEndcap[khzzMvaPfIso] ->Fill(etFake);
	hmuPUEndcap[khzzMvaPfIso] ->Fill(numvertices);
      }
    }

  } // loop over events

    // saving efficiency histos
  // === as a function of eta ===
  char filename[500];
  sprintf(filename,"%s-MuonMisidEtaHighPt.root",outname);
  EfficiencyEvaluator MuonEffEtaHighPt(filename);
  MuonEffEtaHighPt.AddNumerator(RecoEtaHighPt);
  for (int icut=0;icut<(int)hmuID.size();++icut){
    MuonEffEtaHighPt.AddNumerator(hmuEtaHighPt[icut]);
  }

  MuonEffEtaHighPt.SetDenominator(RecoEtaHighPt);
  MuonEffEtaHighPt.ComputeEfficiencies();
  MuonEffEtaHighPt.SetTitle("fake rate vs pT");
  MuonEffEtaHighPt.SetXaxisTitle("electron p_{T}");
  MuonEffEtaHighPt.SetYaxisTitle("Fake rate");
  MuonEffEtaHighPt.SetYaxisMin(0.0);
  MuonEffEtaHighPt.Write();

  sprintf(filename,"%s-MuonMisidEtaLowPt.root",outname);
  EfficiencyEvaluator MuonEffEtaLowPt(filename);
  MuonEffEtaLowPt.AddNumerator(RecoEtaLowPt);
  for (int icut=0;icut<(int)hmuID.size();++icut){
    MuonEffEtaLowPt.AddNumerator(hmuEtaLowPt[icut]);
  }

  MuonEffEtaLowPt.SetDenominator(RecoEtaLowPt);
  MuonEffEtaLowPt.ComputeEfficiencies();
  MuonEffEtaLowPt.SetTitle("fake rate vs eta");
  MuonEffEtaLowPt.SetXaxisTitle("electron #eta");
  MuonEffEtaLowPt.SetYaxisTitle("Fake rate");
  MuonEffEtaLowPt.SetYaxisMin(0.0);
  MuonEffEtaLowPt.Write();

  // === as a function of pt ===
  sprintf(filename,"%s-MuonMisidPtBarrel.root",outname);
  EfficiencyEvaluator MuonEffPtBarrel(filename);
  MuonEffPtBarrel.AddNumerator(RecoPtBarrel);
  for (int icut=0;icut<(int)hmuID.size();++icut){
    MuonEffPtBarrel.AddNumerator(hmuPtBarrel[icut]);
  }

  MuonEffPtBarrel.SetDenominator(RecoPtBarrel);
  MuonEffPtBarrel.ComputeEfficiencies();
  MuonEffPtBarrel.SetTitle("fake rate vs pT");
  MuonEffPtBarrel.SetXaxisTitle("electron pT");
  MuonEffPtBarrel.SetYaxisTitle("Fake rate");
  MuonEffPtBarrel.SetYaxisMin(0.0);
  MuonEffPtBarrel.Write();

  sprintf(filename,"%s-MuonMisidPtEndcap.root",outname);
  EfficiencyEvaluator MuonEffPtEndcap(filename);
  MuonEffPtEndcap.AddNumerator(RecoPtEndcap);
  for (int icut=0;icut<(int)hmuID.size();++icut){
    MuonEffPtEndcap.AddNumerator(hmuPtEndcap[icut]);
  }

  MuonEffPtEndcap.SetDenominator(RecoPtEndcap);
  MuonEffPtEndcap.ComputeEfficiencies();
  MuonEffPtEndcap.SetTitle("fake rate vs pT");
  MuonEffPtEndcap.SetXaxisTitle("electron pT");
  MuonEffPtEndcap.SetYaxisTitle("Fake rate");
  MuonEffPtEndcap.SetYaxisMin(0.0);
  MuonEffPtEndcap.Write();

  // === as a function of PU ===
  sprintf(filename,"%s-MuonMisidPUBarrel.root",outname);
  EfficiencyEvaluator MuonEffPUBarrel(filename);
  MuonEffPUBarrel.AddNumerator(RecoPUBarrel);
  for (int icut=0;icut<(int)hmuID.size();++icut){
    MuonEffPUBarrel.AddNumerator(hmuPUBarrel[icut]);
  }

  MuonEffPUBarrel.SetDenominator(RecoPUBarrel);
  MuonEffPUBarrel.ComputeEfficiencies();
  MuonEffPUBarrel.SetTitle("fake rate vs vertices");
  MuonEffPUBarrel.SetXaxisTitle("# vertices");
  MuonEffPUBarrel.SetYaxisTitle("Fake rate");
  MuonEffPUBarrel.SetYaxisMin(0.0);
  MuonEffPUBarrel.Write();

  sprintf(filename,"%s-MuonMisidPUEndcap.root",outname);
  EfficiencyEvaluator MuonEffPUEndcap(filename);
  MuonEffPUEndcap.AddNumerator(RecoPUEndcap);
  for (int icut=0;icut<(int)hmuID.size();++icut){
    MuonEffPUEndcap.AddNumerator(hmuPUEndcap[icut]);
  }

  MuonEffPUEndcap.SetDenominator(RecoPUEndcap);
  MuonEffPUEndcap.ComputeEfficiencies();
  MuonEffPUEndcap.SetTitle("fake rate vs vertices");
  MuonEffPUEndcap.SetXaxisTitle("# vertices");
  MuonEffPUEndcap.SetYaxisTitle("Fake rate");
  MuonEffPUEndcap.SetYaxisMin(0.0);
  MuonEffPUEndcap.Write();

}
