#include <TFile.h>
#include "../interface/ElectronIDMVA.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

using namespace reco;

//--------------------------------------------------------------------------------------------------
ElectronIDMVA::ElectronIDMVA() :
fMethodname("BDTG method"),
fIsInitialized(kFALSE)
{
  // Constructor.
  for(UInt_t i=0; i<6; ++i) {
    fTMVAReader[i] = 0;
  }
}



//--------------------------------------------------------------------------------------------------
ElectronIDMVA::~ElectronIDMVA()
{
  for(UInt_t i=0; i<6; ++i) {
    if (fTMVAReader[i]) delete fTMVAReader[i];
  }
}

//--------------------------------------------------------------------------------------------------
void ElectronIDMVA::Initialize( std::string methodName,
                                std::string Subdet0Pt10To20Weights , 
                                std::string Subdet1Pt10To20Weights , 
                                std::string Subdet2Pt10To20Weights,
                                std::string Subdet0Pt20ToInfWeights,
                                std::string Subdet1Pt20ToInfWeights, 
                                std::string Subdet2Pt20ToInfWeights,
                                ElectronIDMVA::MVAType type) {

  fIsInitialized = kTRUE;
  fMVAType = type;

  fMethodname = methodName;
    
  for(UInt_t i=0; i<6; ++i) {
    if (fTMVAReader[i]) delete fTMVAReader[i];

    fTMVAReader[i] = new TMVA::Reader( "!Color:!Silent:Error" );  
    fTMVAReader[i]->SetVerbose(kTRUE);

    if (type == kBaseline) {
      fTMVAReader[i]->AddVariable( "SigmaIEtaIEta",         &fMVAVar_EleSigmaIEtaIEta         );
      fTMVAReader[i]->AddVariable( "DEtaIn",                &fMVAVar_EleDEtaIn                );
      fTMVAReader[i]->AddVariable( "DPhiIn",                &fMVAVar_EleDPhiIn                );
      fTMVAReader[i]->AddVariable( "FBrem",                 &fMVAVar_EleFBrem                 );
      fTMVAReader[i]->AddVariable( "SigmaIPhiIPhi",         &fMVAVar_EleSigmaIPhiIPhi         );
      fTMVAReader[i]->AddVariable( "NBrem",                 &fMVAVar_EleNBrem                 );
      fTMVAReader[i]->AddVariable( "OneOverEMinusOneOverP", &fMVAVar_EleOneOverEMinusOneOverP );      
    }
    
    if (type == kNoIPInfo) {
      fTMVAReader[i]->AddVariable( "SigmaIEtaIEta",         &fMVAVar_EleSigmaIEtaIEta         );
      fTMVAReader[i]->AddVariable( "DEtaIn",                &fMVAVar_EleDEtaIn                );
      fTMVAReader[i]->AddVariable( "DPhiIn",                &fMVAVar_EleDPhiIn                );
      fTMVAReader[i]->AddVariable( "FBrem",                 &fMVAVar_EleFBrem                 );
      fTMVAReader[i]->AddVariable( "EOverP",                &fMVAVar_EleEOverP                );
      fTMVAReader[i]->AddVariable( "ESeedClusterOverPout",  &fMVAVar_EleESeedClusterOverPout  );
      fTMVAReader[i]->AddVariable( "SigmaIPhiIPhi",         &fMVAVar_EleSigmaIPhiIPhi         );
      fTMVAReader[i]->AddVariable( "NBrem",                 &fMVAVar_EleNBrem                 );
      fTMVAReader[i]->AddVariable( "OneOverEMinusOneOverP", &fMVAVar_EleOneOverEMinusOneOverP );      
      fTMVAReader[i]->AddVariable( "ESeedClusterOverPIn",   &fMVAVar_EleESeedClusterOverPIn   );
    }
    if (type == kWithIPInfo) {
      fTMVAReader[i]->AddVariable( "SigmaIEtaIEta",         &fMVAVar_EleSigmaIEtaIEta         );
      fTMVAReader[i]->AddVariable( "DEtaIn",                &fMVAVar_EleDEtaIn                );
      fTMVAReader[i]->AddVariable( "DPhiIn",                &fMVAVar_EleDPhiIn                );
      fTMVAReader[i]->AddVariable( "D0",                    &fMVAVar_EleD0                    );
      fTMVAReader[i]->AddVariable( "FBrem",                 &fMVAVar_EleFBrem                 );
      fTMVAReader[i]->AddVariable( "EOverP",                &fMVAVar_EleEOverP                );
      fTMVAReader[i]->AddVariable( "ESeedClusterOverPout",  &fMVAVar_EleESeedClusterOverPout  );
      fTMVAReader[i]->AddVariable( "SigmaIPhiIPhi",         &fMVAVar_EleSigmaIPhiIPhi         );
      fTMVAReader[i]->AddVariable( "NBrem",                 &fMVAVar_EleNBrem                 );
      fTMVAReader[i]->AddVariable( "OneOverEMinusOneOverP", &fMVAVar_EleOneOverEMinusOneOverP );      
      fTMVAReader[i]->AddVariable( "ESeedClusterOverPIn",   &fMVAVar_EleESeedClusterOverPIn   );
      fTMVAReader[i]->AddVariable( "IP3d",                  &fMVAVar_EleIP3d                  );
      fTMVAReader[i]->AddVariable( "IP3dSig",               &fMVAVar_EleIP3dSig               );
    }
    if (type == kIDIsoCombined) {
      fTMVAReader[i]->AddVariable( "SigmaIEtaIEta",         &fMVAVar_EleSigmaIEtaIEta            );
      fTMVAReader[i]->AddVariable( "DEtaIn",                &fMVAVar_EleDEtaIn                   );
      fTMVAReader[i]->AddVariable( "DPhiIn",                &fMVAVar_EleDPhiIn                   );
      fTMVAReader[i]->AddVariable( "D0",                    &fMVAVar_EleD0                       );
      fTMVAReader[i]->AddVariable( "FBrem",                 &fMVAVar_EleFBrem                    );
      fTMVAReader[i]->AddVariable( "EOverP",                &fMVAVar_EleEOverP                   );
      fTMVAReader[i]->AddVariable( "ESeedClusterOverPout",  &fMVAVar_EleESeedClusterOverPout     );
      fTMVAReader[i]->AddVariable( "SigmaIPhiIPhi",         &fMVAVar_EleSigmaIPhiIPhi            );
      fTMVAReader[i]->AddVariable( "OneOverEMinusOneOverP", &fMVAVar_EleOneOverEMinusOneOverP    );      
      fTMVAReader[i]->AddVariable( "ESeedClusterOverPIn",   &fMVAVar_EleESeedClusterOverPIn      );
      fTMVAReader[i]->AddVariable( "IP3d",                  &fMVAVar_EleIP3d                     );
      fTMVAReader[i]->AddVariable( "IP3dSig",               &fMVAVar_EleIP3dSig                  );

      fTMVAReader[i]->AddVariable( "GsfTrackChi2OverNdof",  &fMVAVar_EleGsfTrackChi2OverNdof     );
      fTMVAReader[i]->AddVariable( "dEtaCalo",              &fMVAVar_EledEtaCalo                 );
      fTMVAReader[i]->AddVariable( "dPhiCalo",              &fMVAVar_EledPhiCalo                 );
      fTMVAReader[i]->AddVariable( "R9",                    &fMVAVar_EleR9                       );
      fTMVAReader[i]->AddVariable( "SCEtaWidth",            &fMVAVar_EleSCEtaWidth               );
      fTMVAReader[i]->AddVariable( "SCPhiWidth",            &fMVAVar_EleSCPhiWidth               );
      fTMVAReader[i]->AddVariable( "CovIEtaIPhi",           &fMVAVar_EleCovIEtaIPhi              );
      if (i == 2 || i == 5) {
        fTMVAReader[i]->AddVariable( "PreShowerOverRaw",      &fMVAVar_ElePreShowerOverRaw       );
      }
      fTMVAReader[i]->AddVariable( "ChargedIso03",          &fMVAVar_EleChargedIso03OverPt       );
      fTMVAReader[i]->AddVariable( "NeutralHadronIso03",    &fMVAVar_EleNeutralHadronIso03OverPt );
      fTMVAReader[i]->AddVariable( "GammaIso03",            &fMVAVar_EleGammaIso03OverPt         );
      fTMVAReader[i]->AddVariable( "ChargedIso04",          &fMVAVar_EleChargedIso04OverPt       );
      fTMVAReader[i]->AddVariable( "NeutralHadronIso04",    &fMVAVar_EleNeutralHadronIso04OverPt );
      fTMVAReader[i]->AddVariable( "GammaIso04",            &fMVAVar_EleGammaIso04OverPt         );
    }

    if (i==0) fTMVAReader[i]->BookMVA(fMethodname , Subdet0Pt10To20Weights );
    if (i==1) fTMVAReader[i]->BookMVA(fMethodname , Subdet1Pt10To20Weights );
    if (i==2) fTMVAReader[i]->BookMVA(fMethodname , Subdet2Pt10To20Weights );
    if (i==3) fTMVAReader[i]->BookMVA(fMethodname , Subdet0Pt20ToInfWeights );
    if (i==4) fTMVAReader[i]->BookMVA(fMethodname , Subdet1Pt20ToInfWeights );
    if (i==5) fTMVAReader[i]->BookMVA(fMethodname , Subdet2Pt20ToInfWeights );

  }

  std::cout << "Electron ID MVA Initialization\n";
  std::cout << "MethodName : " << fMethodname << " , type == " << type << std::endl;
  std::cout << "Load weights file : " << Subdet0Pt10To20Weights << std::endl;
  std::cout << "Load weights file : " << Subdet1Pt10To20Weights << std::endl;
  std::cout << "Load weights file : " << Subdet2Pt10To20Weights << std::endl;
  std::cout << "Load weights file : " << Subdet0Pt20ToInfWeights << std::endl;
  std::cout << "Load weights file : " << Subdet1Pt20ToInfWeights << std::endl;
  std::cout << "Load weights file : " << Subdet2Pt20ToInfWeights << std::endl;

}


//--------------------------------------------------------------------------------------------------
Double_t ElectronIDMVA::MVAValue(Double_t ElePt , Double_t EleSCEta,
                                 Double_t EleSigmaIEtaIEta,
                                 Double_t EleDEtaIn,
                                 Double_t EleDPhiIn,
                                 Double_t EleHoverE,
                                 Double_t EleD0,
                                 Double_t EleFBrem,
                                 Double_t EleEOverP,
                                 Double_t EleESeedClusterOverPout,
                                 Double_t EleSigmaIPhiIPhi,
                                 Double_t EleNBrem,
                                 Double_t EleOneOverEMinusOneOverP,
                                 Double_t EleESeedClusterOverPIn,
                                 Double_t EleIP3d,
                                 Double_t EleIP3dSig
  ) {
  
  if (!fIsInitialized) { 
    std::cout << "Error: ElectronIDMVA not properly initialized.\n"; 
    return -9999;
  }

  Int_t subdet = 0;
  if (fabs(EleSCEta) < 1.0) subdet = 0;
  else if (fabs(EleSCEta) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (ElePt > 20.0) ptBin = 1;
  
  //set all input variables
  fMVAVar_EleSigmaIEtaIEta = EleSigmaIEtaIEta;
  fMVAVar_EleDEtaIn = EleDEtaIn;
  fMVAVar_EleDPhiIn = EleDPhiIn;
  fMVAVar_EleHoverE = EleHoverE;
  fMVAVar_EleD0 = EleD0;
  fMVAVar_EleFBrem = EleFBrem;
  fMVAVar_EleEOverP = EleEOverP;
  fMVAVar_EleESeedClusterOverPout = EleESeedClusterOverPout;
  fMVAVar_EleSigmaIPhiIPhi = EleSigmaIPhiIPhi;
  fMVAVar_EleNBrem = EleNBrem;
  fMVAVar_EleOneOverEMinusOneOverP = EleOneOverEMinusOneOverP;
  fMVAVar_EleESeedClusterOverPIn = EleESeedClusterOverPIn;
  fMVAVar_EleIP3d = EleIP3d;
  fMVAVar_EleIP3dSig = EleIP3dSig;

  Double_t mva = -9999;  
  TMVA::Reader *reader = 0;
  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 2 && ptBin == 0) MVABin = 2;
  if (subdet == 0 && ptBin == 1) MVABin = 3;
  if (subdet == 1 && ptBin == 1) MVABin = 4;
  if (subdet == 2 && ptBin == 1) MVABin = 5;
  assert(MVABin >= 0 && MVABin <= 5);
  reader = fTMVAReader[MVABin];
                                                
  mva = reader->EvaluateMVA( fMethodname );

  return mva;
}



//--------------------------------------------------------------------------------------------------
Double_t ElectronIDMVA::MVAValue(const reco::GsfElectron *ele, const reco::Vertex vertex, 
                                 EcalClusterLazyTools myEcalCluster,
                                 const TransientTrackBuilder *transientTrackBuilder,
                                 const reco::PFCandidateCollection &PFCandidates,
                                 double Rho
) {
  
  if (!fIsInitialized) { 
    std::cout << "Error: ElectronIDMVA not properly initialized.\n"; 
    return -9999;
  }

  //initialize
  fMVAVar_EleSigmaIEtaIEta = -9999.0;
  fMVAVar_EleDEtaIn = -9999.0;
  fMVAVar_EleDPhiIn = -9999.0;
  fMVAVar_EleHoverE = -9999.0;
  fMVAVar_EleD0 = -9999.0;
  fMVAVar_EleFBrem = -9999.0;
  fMVAVar_EleEOverP = -9999.0;
  fMVAVar_EleESeedClusterOverPout = -9999.0;
  fMVAVar_EleSigmaIPhiIPhi = -9999.0;
  fMVAVar_EleNBrem = -9999.0;
  fMVAVar_EleOneOverEMinusOneOverP = -9999.0;
  fMVAVar_EleESeedClusterOverPIn = -9999.0;
  fMVAVar_EleIP3d = -9999.0;
  fMVAVar_EleIP3dSig = -9999.0;
  fMVAVar_EleGsfTrackChi2OverNdof = -9999.0;
  fMVAVar_EledEtaCalo = -9999.0;
  fMVAVar_EledPhiCalo = -9999.0;
  fMVAVar_EleR9 = -9999.0;
  fMVAVar_EleSCEtaWidth = -9999.0;
  fMVAVar_EleSCPhiWidth = -9999.0;
  fMVAVar_EleCovIEtaIPhi = -9999.0;
  fMVAVar_ElePreShowerOverRaw = -9999.0;
  fMVAVar_EleChargedIso03OverPt = -9999.0;
  fMVAVar_EleNeutralHadronIso03OverPt = -9999.0;
  fMVAVar_EleGammaIso03OverPt = -9999.0;
  fMVAVar_EleChargedIso04OverPt = -9999.0;
  fMVAVar_EleNeutralHadronIso04OverPt = -9999.0;
  fMVAVar_EleGammaIso04OverPt = -9999.0;
  
  Int_t subdet = 0;
  if (fabs(ele->superCluster()->eta()) < 1.0) subdet = 0;
  else if (fabs(ele->superCluster()->eta()) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (ele->pt() > 20.0) ptBin = 1;
  
  double rho = 0;
  if (!(isnan(Rho) || isinf(Rho))) rho = Rho;

  double electronTrackZ = 0;
  if (ele->gsfTrack().isNonnull()) {
    electronTrackZ = ele->gsfTrack()->dz(vertex.position());
  } else if (ele->closestCtfTrackRef().isNonnull()) {
    electronTrackZ = ele->closestCtfTrackRef()->dz(vertex.position());
  }

  //set all input variables
  fMVAVar_EleSigmaIEtaIEta = ele->sigmaIetaIeta() ; 
  fMVAVar_EleDEtaIn = ele->deltaEtaSuperClusterTrackAtVtx();
  fMVAVar_EleDPhiIn = ele->deltaPhiSuperClusterTrackAtVtx();
  fMVAVar_EleHoverE = ele->hcalOverEcal(); 

  fMVAVar_EleFBrem = ele->fbrem(); 
  fMVAVar_EleEOverP = ele->eSuperClusterOverP(); 
  fMVAVar_EleESeedClusterOverPout = ele->eSeedClusterOverPout(); 

  //temporary fix for weird electrons with Sigma iPhi iPhi == Nan
  //these occur at the sub-percent level
  std::vector<float> vCov = myEcalCluster.localCovariances(*(ele->superCluster()->seed())) ;
  if (!isnan(vCov[2])) fMVAVar_EleSigmaIPhiIPhi = sqrt (vCov[2]);
  else fMVAVar_EleSigmaIPhiIPhi = ele->sigmaIetaIeta();

  fMVAVar_EleNBrem = ele->basicClustersSize() - 1; 
  fMVAVar_EleOneOverEMinusOneOverP = (1.0/(ele->superCluster()->energy())) - 1.0 / ele->gsfTrack()->p(); 
  fMVAVar_EleESeedClusterOverPIn = ele->superCluster()->seed()->energy() / ele->trackMomentumAtVtx().R(); 

  //d0
  if (ele->gsfTrack().isNonnull()) {
    fMVAVar_EleD0 = (-1.0)*ele->gsfTrack()->dxy(vertex.position()); 
  } else if (ele->closestCtfTrackRef().isNonnull()) {
    fMVAVar_EleD0 = (-1.0)*ele->closestCtfTrackRef()->dxy(vertex.position()); 
  } else {
    fMVAVar_EleD0 = -9999.0;
  }
    
  //default values for IP3D
  fMVAVar_EleIP3d = -999.0; 
  fMVAVar_EleIP3dSig = 0.0;
  if (ele->gsfTrack().isNonnull()) {
    const double gsfsign   = ( (-ele->gsfTrack()->dxy(vertex.position()))   >=0 ) ? 1. : -1.;
      
    const reco::TransientTrack &tt = transientTrackBuilder->build(ele->gsfTrack()); 
    const std::pair<bool,Measurement1D> &ip3dpv =  IPTools::absoluteImpactParameter3D(tt,vertex);
    if (ip3dpv.first) {
      double ip3d = gsfsign*ip3dpv.second.value();
      double ip3derr = ip3dpv.second.error();  
      fMVAVar_EleIP3d = ip3d; 
      fMVAVar_EleIP3dSig = ip3d/ip3derr;
    }
  }


  fMVAVar_EleGsfTrackChi2OverNdof = ele->gsfTrack()->chi2() / (double)ele->gsfTrack()->ndof();
  fMVAVar_EledEtaCalo = ele->deltaEtaSeedClusterTrackAtCalo();
  fMVAVar_EledPhiCalo = ele->deltaPhiSeedClusterTrackAtCalo();
  fMVAVar_EleR9 = myEcalCluster.e3x3(*(ele->superCluster()->seed())) / ele->superCluster()->rawEnergy();
  fMVAVar_EleSCEtaWidth = ele->superCluster()->etaWidth();
  fMVAVar_EleSCPhiWidth = ele->superCluster()->phiWidth();
  fMVAVar_EleCovIEtaIPhi = vCov[1];
  fMVAVar_ElePreShowerOverRaw = ele->superCluster()->preshowerEnergy() / ele->superCluster()->rawEnergy();

  

  double chargedIso03 = 0;
  double gammaIso03 = 0;
  double neutralHadronIso03 = 0;
  double chargedIso04 = 0;
  double gammaIso04 = 0;
  double neutralHadronIso04 = 0;
  for (reco::PFCandidateCollection::const_iterator iP = PFCandidates.begin(); 
       iP != PFCandidates.end(); ++iP) {

    //exclude the electron itself
    if(iP->gsfTrackRef().isNonnull() && ele->gsfTrack().isNonnull() &&
       refToPtr(iP->gsfTrackRef()) == refToPtr(ele->gsfTrack())) continue;
    if(iP->trackRef().isNonnull() && ele->closestCtfTrackRef().isNonnull() &&
       refToPtr(iP->trackRef()) == refToPtr(ele->closestCtfTrackRef())) continue;      

    double dr = sqrt(pow(iP->eta() - ele->eta(),2) + pow(acos(cos(iP->phi() - ele->phi())),2));
    
    //charged PFCand
    if (iP->trackRef().isNonnull()) {
      Double_t deltaZ = fabs(iP->trackRef()->dz(vertex.position()) - electronTrackZ);
      if (deltaZ <= 0.1) {
        if (dr < 0.3) chargedIso03 += iP->pt();
        if (dr < 0.4) chargedIso04 += iP->pt();
      }
    } 
    //PF gamma
    else if (iP->particleId() == reco::PFCandidate::gamma) {
      if (iP->pt() > 0.5 
          && 
          fabs(ele->eta() - iP->eta()) >= 0.025                    
        ) {
        if (dr < 0.3) gammaIso03 += iP->pt();
        if (dr < 0.4) gammaIso04 += iP->pt();
      }
    } else {
      if (dr >= 0.07) {
        if (dr < 0.3) neutralHadronIso03 += iP->pt();
        if (dr < 0.4) neutralHadronIso04 += iP->pt();
      }
    }
  }

  fMVAVar_EleChargedIso03OverPt 
    = (chargedIso03 
       - rho * ElectronEffectiveArea(kEleChargedIso03, ele->superCluster()->eta())) / ele->pt();
  fMVAVar_EleGammaIso03OverPt 
    = (gammaIso03 
       - rho * ElectronEffectiveArea(kEleGammaIso03, ele->superCluster()->eta()) 
       + rho * ElectronEffectiveArea(kEleGammaIsoVetoEtaStrip03,ele->superCluster()->eta()))/ele->pt();
  fMVAVar_EleNeutralHadronIso03OverPt 
    = (neutralHadronIso03
       - rho * ElectronEffectiveArea(kEleNeutralHadronIso03, ele->superCluster()->eta()) 
       + rho * ElectronEffectiveArea(kEleNeutralHadronIso007,ele->superCluster()->eta())) / ele->pt();
  fMVAVar_EleChargedIso04OverPt 
    = (chargedIso04 
       - rho * ElectronEffectiveArea(kEleChargedIso04, ele->superCluster()->eta()))/ele->pt();
  fMVAVar_EleGammaIso04OverPt 
    = (gammaIso04 
       - rho * ElectronEffectiveArea(kEleGammaIso04, ele->superCluster()->eta()) 
       + rho * ElectronEffectiveArea(kEleGammaIsoVetoEtaStrip04,ele->superCluster()->eta()))/ele->pt();
  fMVAVar_EleNeutralHadronIso04OverPt
    = (neutralHadronIso04 
       - rho * ElectronEffectiveArea(kEleNeutralHadronIso04, ele->superCluster()->eta()) 
       + rho * ElectronEffectiveArea(kEleNeutralHadronIso007,ele->superCluster()->eta()))/ele->pt();
  


  Double_t mva = -9999;  
  TMVA::Reader *reader = 0;
  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 2 && ptBin == 0) MVABin = 2;
  if (subdet == 0 && ptBin == 1) MVABin = 3;
  if (subdet == 1 && ptBin == 1) MVABin = 4;
  if (subdet == 2 && ptBin == 1) MVABin = 5;
  assert(MVABin >= 0 && MVABin <= 5);
  reader = fTMVAReader[MVABin];
                                
  mva = reader->EvaluateMVA( fMethodname );

//   ***************************************************
//   For DEBUGGING
//   std::cout << "********************************\n";
//   std::cout << "Electron MVA\n";
//   std::cout << ele->pt() << " " << ele->eta() << " " << ele->phi() << " : " << ele->superCluster()->eta() << " : " << MVABin << std::endl;
//   std::cout << fMVAVar_EleSigmaIEtaIEta << "\n"
//        << fMVAVar_EleDEtaIn << "\n"
//        << fMVAVar_EleDPhiIn << "\n"
//        << fMVAVar_EleHoverE << "\n"
//        << fMVAVar_EleD0 << "\n"
//        << fMVAVar_EleFBrem << "\n"
//        << fMVAVar_EleEOverP << "\n"
//        << fMVAVar_EleESeedClusterOverPout << "\n"
//        << fMVAVar_EleSigmaIPhiIPhi << "\n"
//        << fMVAVar_EleNBrem << "\n"
//        << fMVAVar_EleOneOverEMinusOneOverP << "\n"
//        << fMVAVar_EleESeedClusterOverPIn << "\n"
//        << fMVAVar_EleIP3d << "\n"
//        << fMVAVar_EleIP3dSig << "\n"
//        << std::endl;
//   std::cout << "MVA: " << mva << std::endl;
//   std::cout << "********************************\n";

  return mva;
}


//--------------------------------------------------------------------------------------------------
Double_t ElectronIDMVA::MVAValue(const reco::GsfElectron *ele, 
                                 EcalClusterLazyTools myEcalCluster) {
  
  if (!fIsInitialized) { 
    std::cout << "Error: ElectronIDMVA not properly initialized.\n"; 
    return -9999;
  }
  if (fMVAType == kWithIPInfo) {
    std::cout << "Error: ElectronIDMVA was initialized with type that requires vertex information. "
              << "It is not compatible with this accessor. \n";
    return -9999;
  }

  //initialize
  fMVAVar_EleSigmaIEtaIEta = -9999.0;
  fMVAVar_EleDEtaIn = -9999.0;
  fMVAVar_EleDPhiIn = -9999.0;
  fMVAVar_EleHoverE = -9999.0;
  fMVAVar_EleFBrem = -9999.0;
  fMVAVar_EleEOverP = -9999.0;
  fMVAVar_EleESeedClusterOverPout = -9999.0;
  fMVAVar_EleSigmaIPhiIPhi = -9999.0;
  fMVAVar_EleNBrem = -9999.0;
  fMVAVar_EleOneOverEMinusOneOverP = -9999.0;
  fMVAVar_EleESeedClusterOverPIn = -9999.0;

  Int_t subdet = 0;
  if (fabs(ele->superCluster()->eta()) < 1.0) subdet = 0;
  else if (fabs(ele->superCluster()->eta()) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (ele->pt() > 20.0) ptBin = 1;
  
  //set all input variables
  fMVAVar_EleSigmaIEtaIEta = ele->sigmaIetaIeta() ; 
  fMVAVar_EleDEtaIn = ele->deltaEtaSuperClusterTrackAtVtx();
  fMVAVar_EleDPhiIn = ele->deltaPhiSuperClusterTrackAtVtx();
  fMVAVar_EleHoverE = ele->hcalOverEcal(); 

  fMVAVar_EleFBrem = ele->fbrem(); 
  fMVAVar_EleEOverP = ele->eSuperClusterOverP(); 
  fMVAVar_EleESeedClusterOverPout = ele->eSeedClusterOverPout(); 

  //temporary fix for weird electrons with Sigma iPhi iPhi == Nan
  //these occur at the sub-percent level
  std::vector<float> vCov = myEcalCluster.localCovariances(*(ele->superCluster()->seed())) ;
  if (!isnan(vCov[2])) fMVAVar_EleSigmaIPhiIPhi = sqrt (vCov[2]);
  else fMVAVar_EleSigmaIPhiIPhi = ele->sigmaIetaIeta();

  fMVAVar_EleNBrem = ele->basicClustersSize() - 1; 
  fMVAVar_EleOneOverEMinusOneOverP = (1.0/(ele->superCluster()->energy())) - 1.0 / ele->gsfTrack()->p(); 
  fMVAVar_EleESeedClusterOverPIn = ele->superCluster()->seed()->energy() / ele->trackMomentumAtVtx().R(); 


  Double_t mva = -9999;  
  TMVA::Reader *reader = 0;
  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 2 && ptBin == 0) MVABin = 2;
  if (subdet == 0 && ptBin == 1) MVABin = 3;
  if (subdet == 1 && ptBin == 1) MVABin = 4;
  if (subdet == 2 && ptBin == 1) MVABin = 5;
  assert(MVABin >= 0 && MVABin <= 5);
  reader = fTMVAReader[MVABin];
                                
  mva = reader->EvaluateMVA( fMethodname );

//   ***************************************************
//   For DEBUGGING
//   std::cout << "********************************\n";
//   std::cout << "Electron MVA\n";
//   std::cout << ele->pt() << " " << ele->eta() << " " << ele->phi() << " : " << ele->superCluster()->eta() << " : " << MVABin << std::endl;
//   std::cout << fMVAVar_EleSigmaIEtaIEta << "\n"
//        << fMVAVar_EleDEtaIn << "\n"
//        << fMVAVar_EleDPhiIn << "\n"
//        << fMVAVar_EleHoverE << "\n"
//        << fMVAVar_EleFBrem << "\n"
//        << fMVAVar_EleEOverP << "\n"
//        << fMVAVar_EleESeedClusterOverPout << "\n"
//        << fMVAVar_EleSigmaIPhiIPhi << "\n"
//        << fMVAVar_EleNBrem << "\n"
//        << fMVAVar_EleOneOverEMinusOneOverP << "\n"
//        << fMVAVar_EleESeedClusterOverPIn << "\n"
//        << std::endl;
//   std::cout << "MVA: " << mva << std::endl;
//   std::cout << "********************************\n";

  return mva;
}


double ElectronIDMVA::ElectronEffectiveArea(ElectronEffectiveAreaType type, double Eta) {

  double EffectiveArea = 0;

  if (fabs(Eta) < 1.0) {
    if (type == kEleChargedIso03) EffectiveArea = 0.000;
    if (type == kEleNeutralHadronIso03) EffectiveArea = 0.017;
    if (type == kEleGammaIso03) EffectiveArea = 0.045;
    if (type == kEleGammaIsoVetoEtaStrip03) EffectiveArea = 0.014;
    if (type == kEleChargedIso04) EffectiveArea = 0.000;
    if (type == kEleNeutralHadronIso04) EffectiveArea = 0.034;
    if (type == kEleGammaIso04) EffectiveArea = 0.079;
    if (type == kEleGammaIsoVetoEtaStrip04) EffectiveArea = 0.014;
    if (type == kEleNeutralHadronIso007) EffectiveArea = 0.000;
    if (type == kEleHoverE) EffectiveArea = 0.00016;
    if (type == kEleHcalDepth1OverEcal) EffectiveArea = 0.00016;
    if (type == kEleHcalDepth2OverEcal) EffectiveArea = 0.00000;    
  } else if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) {
    if (type == kEleChargedIso03) EffectiveArea = 0.000;
    if (type == kEleNeutralHadronIso03) EffectiveArea = 0.025;
    if (type == kEleGammaIso03) EffectiveArea = 0.052;
    if (type == kEleGammaIsoVetoEtaStrip03) EffectiveArea = 0.030;
    if (type == kEleChargedIso04) EffectiveArea = 0.000;
    if (type == kEleNeutralHadronIso04) EffectiveArea = 0.050;
    if (type == kEleGammaIso04) EffectiveArea = 0.073;
    if (type == kEleGammaIsoVetoEtaStrip04) EffectiveArea = 0.030;
    if (type == kEleNeutralHadronIso007) EffectiveArea = 0.000;
    if (type == kEleHoverE) EffectiveArea = 0.00022;
    if (type == kEleHcalDepth1OverEcal) EffectiveArea = 0.00022;
    if (type == kEleHcalDepth2OverEcal) EffectiveArea = 0.00000;    
  } else if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) {
    if (type == kEleChargedIso03) EffectiveArea = 0.000;
    if (type == kEleNeutralHadronIso03) EffectiveArea = 0.030;
    if (type == kEleGammaIso03) EffectiveArea = 0.170;
    if (type == kEleGammaIsoVetoEtaStrip03) EffectiveArea = 0.134;
    if (type == kEleChargedIso04) EffectiveArea = 0.000;
    if (type == kEleNeutralHadronIso04) EffectiveArea = 0.060;
    if (type == kEleGammaIso04) EffectiveArea = 0.187;
    if (type == kEleGammaIsoVetoEtaStrip04) EffectiveArea = 0.134;
    if (type == kEleNeutralHadronIso007) EffectiveArea = 0.000;
    if (type == kEleHoverE) EffectiveArea = 0.00030;
    if (type == kEleHcalDepth1OverEcal) EffectiveArea = 0.00026;
    if (type == kEleHcalDepth2OverEcal) EffectiveArea = 0.00002;        
  } else if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.25 ) {
    if (type == kEleChargedIso03) EffectiveArea = 0.000;
    if (type == kEleNeutralHadronIso03) EffectiveArea = 0.022;
    if (type == kEleGammaIso03) EffectiveArea = 0.623;
    if (type == kEleGammaIsoVetoEtaStrip03) EffectiveArea = 0.516;
    if (type == kEleChargedIso04) EffectiveArea = 0.000;
    if (type == kEleNeutralHadronIso04) EffectiveArea = 0.055;
    if (type == kEleGammaIso04) EffectiveArea = 0.659;
    if (type == kEleGammaIsoVetoEtaStrip04) EffectiveArea = 0.517;
    if (type == kEleNeutralHadronIso007) EffectiveArea = 0.000;
    if (type == kEleHoverE) EffectiveArea = 0.00054;
    if (type == kEleHcalDepth1OverEcal) EffectiveArea = 0.00045;
    if (type == kEleHcalDepth2OverEcal) EffectiveArea = 0.00003;
  } else if (fabs(Eta) >= 2.25 && fabs(Eta) < 2.5 ) {
    if (type == kEleChargedIso03) EffectiveArea = 0.000;
    if (type == kEleNeutralHadronIso03) EffectiveArea = 0.018;
    if (type == kEleGammaIso03) EffectiveArea = 1.198;
    if (type == kEleGammaIsoVetoEtaStrip03) EffectiveArea = 1.049;
    if (type == kEleChargedIso04) EffectiveArea = 0.000;
    if (type == kEleNeutralHadronIso04) EffectiveArea = 0.073;
    if (type == kEleGammaIso04) EffectiveArea = 1.258;
    if (type == kEleGammaIsoVetoEtaStrip04) EffectiveArea = 1.051;
    if (type == kEleNeutralHadronIso007) EffectiveArea = 0.000;
    if (type == kEleHoverE) EffectiveArea = 0.00082;
    if (type == kEleHcalDepth1OverEcal) EffectiveArea = 0.00066;
    if (type == kEleHcalDepth2OverEcal) EffectiveArea = 0.00004;
  }
    
  return EffectiveArea;  
}


