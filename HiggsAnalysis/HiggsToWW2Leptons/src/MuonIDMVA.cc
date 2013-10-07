#include <TFile.h>
#include "../interface/MuonIDMVA.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

using namespace reco;

//--------------------------------------------------------------------------------------------------
MuonIDMVA::MuonIDMVA() :
fMethodname("BDTG method"),
fIsInitialized(kFALSE)
{
  // Constructor.
  for(UInt_t i=0; i<6; ++i) {
    fTMVAReader[i] = 0;
  }
}


//--------------------------------------------------------------------------------------------------
MuonIDMVA::~MuonIDMVA()
{
  for(UInt_t i=0; i<6; ++i) {
    if (fTMVAReader[i]) delete fTMVAReader[i];
  }
}

//--------------------------------------------------------------------------------------------------
void MuonIDMVA::Initialize( std::string methodName,
                            std::string Subdet0Pt10To14p5Weights , 
                            std::string Subdet1Pt10To14p5Weights , 
                            std::string Subdet0Pt14p5To20Weights,
                            std::string Subdet1Pt14p5To20Weights, 
                            std::string Subdet0Pt20ToInfWeights, 
                            std::string Subdet1Pt20ToInfWeights,
                            MuonIDMVA::MVAType type) {
  
  fIsInitialized = kTRUE;
  
  fMethodname = methodName;
    
  for(UInt_t i=0; i<6; ++i) {
    if (fTMVAReader[i]) delete fTMVAReader[i];

    fTMVAReader[i] = new TMVA::Reader( "!Color:!Silent:Error" );  
    fTMVAReader[i]->SetVerbose(kTRUE);

    if (type == kIDIsoCombinedDetIso) {
      fTMVAReader[i]->AddVariable( "TkNchi2",              &fMVAVar_MuTkNchi2               );
      fTMVAReader[i]->AddVariable( "GlobalNchi2",          &fMVAVar_MuGlobalNchi2           );
      fTMVAReader[i]->AddVariable( "NValidHits",           &fMVAVar_MuNValidHits            );
      fTMVAReader[i]->AddVariable( "NTrackerHits",         &fMVAVar_MuNTrackerHits          );
      fTMVAReader[i]->AddVariable( "NPixelHits",           &fMVAVar_MuNPixelHits            );
      fTMVAReader[i]->AddVariable( "NMatches",             &fMVAVar_MuNMatches              );
      fTMVAReader[i]->AddVariable( "D0",                   &fMVAVar_MuD0                    );      
      fTMVAReader[i]->AddVariable( "IP3d",                 &fMVAVar_MuIP3d                  );      
      fTMVAReader[i]->AddVariable( "IP3dSig",              &fMVAVar_MuIP3dSig               );      
      fTMVAReader[i]->AddVariable( "TrkKink",              &fMVAVar_MuTrkKink               );      
      fTMVAReader[i]->AddVariable( "SegmentCompatibility", &fMVAVar_MuSegmentCompatibility  );      
      fTMVAReader[i]->AddVariable( "CaloCompatibility",    &fMVAVar_MuCaloCompatibility     );      
      fTMVAReader[i]->AddVariable( "HadEnergyOverPt",      &fMVAVar_MuHadEnergyOverPt       );      
      fTMVAReader[i]->AddVariable( "EmEnergyOverPt",       &fMVAVar_MuEmEnergyOverPt        );      
      fTMVAReader[i]->AddVariable( "HadS9EnergyOverPt",    &fMVAVar_MuHadS9EnergyOverPt     );      
      fTMVAReader[i]->AddVariable( "EmS9EnergyOverPt",     &fMVAVar_MuEmS9EnergyOverPt      );      
      fTMVAReader[i]->AddVariable( "TrkIso03OverPt",       &fMVAVar_MuTrkIso03OverPt        );
      fTMVAReader[i]->AddVariable( "EMIso03OverPt",        &fMVAVar_MuEMIso03OverPt         );
      fTMVAReader[i]->AddVariable( "HadIso03OverPt",       &fMVAVar_MuHadIso03OverPt        );
      fTMVAReader[i]->AddVariable( "TrkIso05OverPt",       &fMVAVar_MuTrkIso05OverPt        );
      fTMVAReader[i]->AddVariable( "EMIso05OverPt",        &fMVAVar_MuEMIso05OverPt         );
      fTMVAReader[i]->AddVariable( "HadIso05OverPt",       &fMVAVar_MuHadIso05OverPt        );
    
    }
    
    if (i==0) fTMVAReader[i]->BookMVA(fMethodname , Subdet0Pt10To14p5Weights );
    if (i==1) fTMVAReader[i]->BookMVA(fMethodname , Subdet1Pt10To14p5Weights );
    if (i==2) fTMVAReader[i]->BookMVA(fMethodname , Subdet0Pt14p5To20Weights );
    if (i==3) fTMVAReader[i]->BookMVA(fMethodname , Subdet1Pt14p5To20Weights );
    if (i==4) fTMVAReader[i]->BookMVA(fMethodname , Subdet0Pt20ToInfWeights  );
    if (i==5) fTMVAReader[i]->BookMVA(fMethodname , Subdet1Pt20ToInfWeights  );

  }

  std::cout << "Muon ID MVA Initialization\n";
  std::cout << "MethodName : " << fMethodname << " , type == " << type << std::endl;
  std::cout << "Load weights file : " << Subdet0Pt10To14p5Weights << std::endl;
  std::cout << "Load weights file : " << Subdet1Pt10To14p5Weights << std::endl;
  std::cout << "Load weights file : " << Subdet0Pt14p5To20Weights << std::endl;
  std::cout << "Load weights file : " << Subdet1Pt14p5To20Weights << std::endl;
  std::cout << "Load weights file : " << Subdet0Pt20ToInfWeights << std::endl;
  std::cout << "Load weights file : " << Subdet1Pt20ToInfWeights << std::endl;

}


//--------------------------------------------------------------------------------------------------
Double_t MuonIDMVA::MVAValue(Double_t MuPt , Double_t MuEta,
                             Double_t                   MuTkNchi2, 
                             Double_t                   MuGlobalNchi2, 
                             Double_t                   MuNValidHits, 
                             Double_t                   MuNTrackerHits, 
                             Double_t                   MuNPixelHits, 
                             Double_t                   MuNMatches, 
                             Double_t                   MuD0, 
                             Double_t                   MuIP3d, 
                             Double_t                   MuIP3dSig, 
                             Double_t                   MuTrkKink, 
                             Double_t                   MuSegmentCompatibility, 
                             Double_t                   MuCaloCompatibility, 
                             Double_t                   MuHadEnergyOverPt, 
                             Double_t                   MuHoEnergyOverPt, 
                             Double_t                   MuEmEnergyOverPt, 
                             Double_t                   MuHadS9EnergyOverPt, 
                             Double_t                   MuHoS9EnergyOverPt, 
                             Double_t                   MuEmS9EnergyOverPt,
                             Double_t                   MuTrkIso03OverPt,
                             Double_t                   MuEMIso03OverPt,
                             Double_t                   MuHadIso03OverPt,
                             Double_t                   MuTrkIso05OverPt,
                             Double_t                   MuEMIso05OverPt,
                             Double_t                   MuHadIso05OverPt,
                             Bool_t                     printDebug                            
  ) {
  
  if (!fIsInitialized) { 
    std::cout << "Error: MuonIDMVA not properly initialized.\n"; 
    return -9999;
  }

  Int_t subdet = 0;
  if (fabs(MuEta) < 1.479) subdet = 0;
  else subdet = 1;
  Int_t ptBin = 0;
  if (MuPt > 14.5) ptBin = 1;
  if (MuPt > 20.0) ptBin = 2;

  
  //set all input variables
  fMVAVar_MuTkNchi2              = MuTkNchi2; 
  fMVAVar_MuGlobalNchi2          = MuGlobalNchi2; 
  fMVAVar_MuNValidHits           = MuNValidHits; 
  fMVAVar_MuNTrackerHits         = MuNTrackerHits; 
  fMVAVar_MuNPixelHits           = MuNPixelHits;  
  fMVAVar_MuNMatches             = MuNMatches; 
  fMVAVar_MuD0                   = MuD0; 
  fMVAVar_MuIP3d                 = MuIP3d; 
  fMVAVar_MuIP3dSig              = MuIP3dSig; 
  fMVAVar_MuTrkKink              = MuTrkKink; 
  fMVAVar_MuSegmentCompatibility = MuSegmentCompatibility; 
  fMVAVar_MuCaloCompatibility    = MuCaloCompatibility; 
  fMVAVar_MuHadEnergyOverPt      = MuHadEnergyOverPt; 
  fMVAVar_MuHoEnergyOverPt       = MuHoEnergyOverPt; 
  fMVAVar_MuEmEnergyOverPt       = MuEmEnergyOverPt; 
  fMVAVar_MuHadS9EnergyOverPt    = MuHadS9EnergyOverPt; 
  fMVAVar_MuHoS9EnergyOverPt     = MuHoS9EnergyOverPt; 
  fMVAVar_MuEmS9EnergyOverPt     = MuEmS9EnergyOverPt; 
  fMVAVar_MuTrkIso03OverPt       = MuTrkIso03OverPt; 
  fMVAVar_MuEMIso03OverPt        = MuEMIso03OverPt; 
  fMVAVar_MuHadIso03OverPt       = MuHadIso03OverPt; 
  fMVAVar_MuTrkIso05OverPt       = MuTrkIso05OverPt; 
  fMVAVar_MuEMIso05OverPt        = MuEMIso05OverPt; 
  fMVAVar_MuHadIso05OverPt       = MuHadIso05OverPt; 

  Double_t mva = -9999;  
  TMVA::Reader *reader = 0;

  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 0 && ptBin == 1) MVABin = 2;
  if (subdet == 1 && ptBin == 1) MVABin = 3;
  if (subdet == 0 && ptBin == 2) MVABin = 4;
  if (subdet == 1 && ptBin == 2) MVABin = 5;
  assert(MVABin >= 0 && MVABin <= 5);
  reader = fTMVAReader[MVABin];
                                                
  mva = reader->EvaluateMVA( fMethodname );

  if (printDebug) {
    std::cout << "Debug Muon MVA: "
	 << MuPt << " " << MuEta << " --> MVABin " << MVABin << " : "     
	 << fMVAVar_MuTkNchi2              << " " 
	 << fMVAVar_MuGlobalNchi2          << " " 
	 << fMVAVar_MuNValidHits           << " " 
	 << fMVAVar_MuNTrackerHits         << " " 
	 << fMVAVar_MuNPixelHits           << " "  
	 << fMVAVar_MuNMatches             << " " 
	 << fMVAVar_MuD0                   << " " 
	 << fMVAVar_MuIP3d                 << " " 
	 << fMVAVar_MuIP3dSig              << " " 
	 << fMVAVar_MuTrkKink              << " " 
	 << fMVAVar_MuSegmentCompatibility << " " 
	 << fMVAVar_MuCaloCompatibility    << " " 
	 << fMVAVar_MuHadEnergyOverPt      << " " 
	 << fMVAVar_MuHoEnergyOverPt       << " " 
	 << fMVAVar_MuEmEnergyOverPt       << " " 
	 << fMVAVar_MuHadS9EnergyOverPt    << " " 
	 << fMVAVar_MuHoS9EnergyOverPt     << " " 
	 << fMVAVar_MuEmS9EnergyOverPt     << " " 
	 << fMVAVar_MuTrkIso03OverPt   << " " 
	 << fMVAVar_MuEMIso03OverPt   << " " 
	 << fMVAVar_MuHadIso03OverPt   << " " 
	 << fMVAVar_MuTrkIso05OverPt   << " " 
	 << fMVAVar_MuEMIso05OverPt   << " " 
	 << fMVAVar_MuHadIso05OverPt   << " " 
	 << " === : === "
	 << mva 
	 << std::endl;
  }

  return mva;
}



//--------------------------------------------------------------------------------------------------
double MuonIDMVA::MVAValue(const reco::Muon *mu, const reco::Vertex vertex, 
                             const TransientTrackBuilder *transientTrackBuilder,
                             double Rho,
                             bool printDebug) {
  
  if (!fIsInitialized) { 
    std::cout << "Error: MuonIDMVA not properly initialized.\n"; 
    return -9999;
  }

  TrackRef muTrk = mu->track();
  if (muTrk.isNull()) {
    muTrk = mu->standAloneMuon();
  }

  double muNchi2 = 0.0; 
  if (mu->combinedMuon().isNonnull()) { 
    muNchi2 = mu->combinedMuon()->chi2() / (Double_t)mu->combinedMuon()->ndof(); 
  } else if (mu->standAloneMuon().isNonnull()) { 
    muNchi2 = mu->standAloneMuon()->chi2() / (Double_t)mu->standAloneMuon()->ndof(); 
  } else if (mu->track().isNonnull()) { 
    muNchi2 = mu->track()->chi2() / (Double_t)mu->track()->ndof(); 
  }

  double rho = 0;
  if (!(isnan(Rho) || isinf(Rho))) rho = Rho;

  Int_t subdet = 0;
  if (fabs(muTrk->eta()) < 1.479) subdet = 0;
  else subdet = 1;
  Int_t ptBin = 0;
  if (muTrk->pt() > 14.5) ptBin = 1;
  if (muTrk->pt() > 20.0) ptBin = 2;

  //set all input variables
  fMVAVar_MuTkNchi2              = muTrk->chi2() / (Double_t)muTrk->ndof();
  fMVAVar_MuGlobalNchi2          = muNchi2;
  fMVAVar_MuNValidHits           = mu->globalTrack().isNonnull() ? mu->globalTrack()->hitPattern().numberOfValidMuonHits() : 0;
  fMVAVar_MuNTrackerHits         = muTrk->numberOfValidHits();
  fMVAVar_MuNPixelHits           = muTrk->hitPattern().numberOfValidPixelHits();
  fMVAVar_MuNMatches             = mu->numberOfMatches();

  fMVAVar_MuD0                   = (-1.0)*muTrk->dxy(vertex.position());

  const reco::TransientTrack &tt = transientTrackBuilder->build(muTrk); 
  const double thesign   = ( (muTrk->dxy(vertex.position()))   >=0 ) ? 1. : -1.;
  const std::pair<bool,Measurement1D> &ip3dpv =  IPTools::absoluteImpactParameter3D(tt,vertex);
  if (ip3dpv.first) {
    fMVAVar_MuIP3d                 = thesign*ip3dpv.second.value();
    fMVAVar_MuIP3dSig              = thesign*ip3dpv.second.value() / ip3dpv.second.error();
  }
  else {
    fMVAVar_MuIP3d                 = -99;
    fMVAVar_MuIP3dSig              = -99999;
  }

  fMVAVar_MuTrkKink              = mu->combinedQuality().trkKink;


  fMVAVar_MuSegmentCompatibility = muon::segmentCompatibility(*mu, reco::Muon::SegmentAndTrackArbitration);
  fMVAVar_MuCaloCompatibility    = mu->caloCompatibility();

  fMVAVar_MuHadEnergyOverPt      = (mu->calEnergy().had - rho*MuonEffectiveArea(kMuHadEnergy,muTrk->eta()))/muTrk->pt();
  fMVAVar_MuEmEnergyOverPt       = (mu->calEnergy().em - rho*MuonEffectiveArea(kMuEmEnergy,muTrk->eta()))/muTrk->pt();
  fMVAVar_MuHadS9EnergyOverPt    = (mu->calEnergy().hadS9 - rho*MuonEffectiveArea(kMuHadS9Energy,muTrk->eta()))/muTrk->pt();
  fMVAVar_MuEmS9EnergyOverPt     = (mu->calEnergy().emS9 - rho*MuonEffectiveArea(kMuEmS9Energy,muTrk->eta()))/muTrk->pt();
  fMVAVar_MuTrkIso03OverPt       = (mu->isolationR03().sumPt - rho*MuonEffectiveArea(kMuTrkIso03,muTrk->eta()))/muTrk->pt();
  fMVAVar_MuEMIso03OverPt        = (mu->isolationR03().emEt - rho*MuonEffectiveArea(kMuEMIso03,muTrk->eta()))/muTrk->pt();
  fMVAVar_MuHadIso03OverPt       = (mu->isolationR03().hadEt - rho*MuonEffectiveArea(kMuHadIso03,muTrk->eta()))/muTrk->pt();
  fMVAVar_MuTrkIso05OverPt       = (mu->isolationR05().sumPt - rho*MuonEffectiveArea(kMuTrkIso05,muTrk->eta()))/muTrk->pt();
  fMVAVar_MuEMIso05OverPt        = (mu->isolationR05().emEt - rho*MuonEffectiveArea(kMuEMIso05,muTrk->eta()))/muTrk->pt();
  fMVAVar_MuHadIso05OverPt       = (mu->isolationR05().hadEt - rho*MuonEffectiveArea(kMuHadIso05,muTrk->eta()))/muTrk->pt();

  double mva = -9999;  
  TMVA::Reader *reader = 0;

  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 0 && ptBin == 1) MVABin = 2;
  if (subdet == 1 && ptBin == 1) MVABin = 3;
  if (subdet == 0 && ptBin == 2) MVABin = 4;
  if (subdet == 1 && ptBin == 2) MVABin = 5;
  assert(MVABin >= 0 && MVABin <= 5);
  reader = fTMVAReader[MVABin];
                                                
  mva = reader->EvaluateMVA( fMethodname );

  if (printDebug) {
    std::cout << "Debug Muon MVA: "
              << mu->pt() << " " << mu->eta() << " " << mu->phi() << " : "
              << muTrk->pt() << " " << muTrk->eta() << " --> MVABin " << MVABin << " : "     
              << fMVAVar_MuTkNchi2              << " " 
              << fMVAVar_MuGlobalNchi2          << " " 
              << fMVAVar_MuNValidHits           << " " 
              << fMVAVar_MuNTrackerHits         << " " 
              << fMVAVar_MuNPixelHits           << " "  
              << fMVAVar_MuNMatches             << " " 
              << fMVAVar_MuD0                   << " " 
              << fMVAVar_MuIP3d                 << " " 
              << fMVAVar_MuIP3dSig              << " " 
              << fMVAVar_MuTrkKink              << " " 
              << fMVAVar_MuSegmentCompatibility << " " 
              << fMVAVar_MuCaloCompatibility    << " " 
              << fMVAVar_MuHadEnergyOverPt      << " " 
              << fMVAVar_MuHoEnergyOverPt       << " " 
              << fMVAVar_MuEmEnergyOverPt       << " " 
              << fMVAVar_MuHadS9EnergyOverPt    << " " 
              << fMVAVar_MuHoS9EnergyOverPt     << " " 
              << fMVAVar_MuEmS9EnergyOverPt     << " " 
              << fMVAVar_MuChargedIso03OverPt   << " " 
              << fMVAVar_MuNeutralIso03OverPt   << " " 
              << fMVAVar_MuChargedIso04OverPt   << " " 
              << fMVAVar_MuNeutralIso04OverPt   << " " 
              << fMVAVar_MuTrkIso03OverPt   << " " 
              << fMVAVar_MuEMIso03OverPt   << " " 
              << fMVAVar_MuHadIso03OverPt   << " " 
              << fMVAVar_MuTrkIso05OverPt   << " " 
              << fMVAVar_MuEMIso05OverPt   << " " 
              << fMVAVar_MuHadIso05OverPt   << " " 
              << " === : === "
              << mva 
              << std::endl;
  }

  return mva;
}


double MuonIDMVA::MuonEffectiveArea(MuonEffectiveAreaType type, double Eta) {

  double EffectiveArea = 0;
  if (fabs(Eta) < 1.0) {
    if (type == kMuChargedIso03) EffectiveArea = 0.000;
    if (type == kMuNeutralIso03) EffectiveArea = 0.080;
    if (type == kMuChargedIso04) EffectiveArea = 0.000;
    if (type == kMuNeutralIso04) EffectiveArea = 0.163;
    if (type == kMuHadEnergy)    EffectiveArea = 0.000;
    if (type == kMuHoEnergy)     EffectiveArea = 0.000;
    if (type == kMuEmEnergy)     EffectiveArea = 0.000;
    if (type == kMuHadS9Energy)  EffectiveArea = 0.016;
    if (type == kMuHoS9Energy)   EffectiveArea = 0.000;
    if (type == kMuEmS9Energy)   EffectiveArea = 0.000;
    if (type == kMuTrkIso03)     EffectiveArea = 0.000;
    if (type == kMuEMIso03)      EffectiveArea = 0.080;
    if (type == kMuHadIso03)     EffectiveArea = 0.025;
    if (type == kMuTrkIso05)     EffectiveArea = 0.000;
    if (type == kMuEMIso05)      EffectiveArea = 0.290;
    if (type == kMuHadIso05)     EffectiveArea = 0.091;
  } else if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) {
    if (type == kMuChargedIso03) EffectiveArea = 0.000;
    if (type == kMuNeutralIso03) EffectiveArea = 0.083;
    if (type == kMuChargedIso04) EffectiveArea = 0.000;
    if (type == kMuNeutralIso04) EffectiveArea = 0.168;
    if (type == kMuHadEnergy)    EffectiveArea = 0.005;
    if (type == kMuHoEnergy)     EffectiveArea = 0.000;
    if (type == kMuEmEnergy)     EffectiveArea = 0.000;
    if (type == kMuHadS9Energy)  EffectiveArea = 0.041;
    if (type == kMuHoS9Energy)   EffectiveArea = 0.000;
    if (type == kMuEmS9Energy)   EffectiveArea = 0.000;
    if (type == kMuTrkIso03)     EffectiveArea = 0.000;
    if (type == kMuEMIso03)      EffectiveArea = 0.043;
    if (type == kMuHadIso03)     EffectiveArea = 0.028;
    if (type == kMuTrkIso05)     EffectiveArea = 0.000;
    if (type == kMuEMIso05)      EffectiveArea = 0.184;
    if (type == kMuHadIso05)     EffectiveArea = 0.106;
  } else if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) {
    if (type == kMuChargedIso03) EffectiveArea = 0.000;
    if (type == kMuNeutralIso03) EffectiveArea = 0.060;
    if (type == kMuChargedIso04) EffectiveArea = 0.000;
    if (type == kMuNeutralIso04) EffectiveArea = 0.131;
    if (type == kMuHadEnergy)    EffectiveArea = 0.020;
    if (type == kMuHoEnergy)     EffectiveArea = 0.000;
    if (type == kMuEmEnergy)     EffectiveArea = 0.000;
    if (type == kMuHadS9Energy)  EffectiveArea = 0.072;
    if (type == kMuHoS9Energy)   EffectiveArea = 0.000;
    if (type == kMuEmS9Energy)   EffectiveArea = 0.000;
    if (type == kMuTrkIso03)     EffectiveArea = 0.000;
    if (type == kMuEMIso03)      EffectiveArea = 0.025;
    if (type == kMuHadIso03)     EffectiveArea = 0.036;
    if (type == kMuTrkIso05)     EffectiveArea = 0.000;
    if (type == kMuEMIso05)      EffectiveArea = 0.124;
    if (type == kMuHadIso05)     EffectiveArea = 0.140;
  } else if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.25 ) {
    if (type == kMuChargedIso03) EffectiveArea = 0.000;
    if (type == kMuNeutralIso03) EffectiveArea = 0.066;
    if (type == kMuChargedIso04) EffectiveArea = 0.000;
    if (type == kMuNeutralIso04) EffectiveArea = 0.149;
    if (type == kMuHadEnergy)    EffectiveArea = 0.056;
    if (type == kMuHoEnergy)     EffectiveArea = 0.000;
    if (type == kMuEmEnergy)     EffectiveArea = 0.000;
    if (type == kMuHadS9Energy)  EffectiveArea = 0.148;
    if (type == kMuHoS9Energy)   EffectiveArea = 0.000;
    if (type == kMuEmS9Energy)   EffectiveArea = 0.000;
    if (type == kMuTrkIso03)     EffectiveArea = 0.000;
    if (type == kMuEMIso03)      EffectiveArea = 0.025;
    if (type == kMuHadIso03)     EffectiveArea = 0.050;
    if (type == kMuTrkIso05)     EffectiveArea = 0.000;
    if (type == kMuEMIso05)      EffectiveArea = 0.120;
    if (type == kMuHadIso05)     EffectiveArea = 0.186;
  } else if (fabs(Eta) >= 2.25 && fabs(Eta) < 2.4 ) {
    if (type == kMuChargedIso03) EffectiveArea = 0.000;
    if (type == kMuNeutralIso03) EffectiveArea = 0.098;
    if (type == kMuChargedIso04) EffectiveArea = 0.000;
    if (type == kMuNeutralIso04) EffectiveArea = 0.200;
    if (type == kMuHadEnergy)    EffectiveArea = 0.093;
    if (type == kMuHoEnergy)     EffectiveArea = 0.000;
    if (type == kMuEmEnergy)     EffectiveArea = 0.000;
    if (type == kMuHadS9Energy)  EffectiveArea = 0.260;
    if (type == kMuHoS9Energy)   EffectiveArea = 0.000;
    if (type == kMuEmS9Energy)   EffectiveArea = 0.000;
    if (type == kMuTrkIso03)     EffectiveArea = 0.000;
    if (type == kMuEMIso03)      EffectiveArea = 0.027;
    if (type == kMuHadIso03)     EffectiveArea = 0.060;
    if (type == kMuTrkIso05)     EffectiveArea = 0.000;
    if (type == kMuEMIso05)      EffectiveArea = 0.139;
    if (type == kMuHadIso05)     EffectiveArea = 0.228;
  }
  return EffectiveArea;
}

