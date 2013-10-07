//--------------------------------------------------------------------------------------------------
// $Id $
//
// ElectronIDMVA
//
// Helper Class for applying MVA electron ID selection
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef HIGGSANALYSIS_ElectronIDMVA_H
#define HIGGSANALYSIS_ElectronIDMVA_H

#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

class ElectronIDMVA {
  public:
    ElectronIDMVA();
    ~ElectronIDMVA(); 

    enum MVAType {
      kBaseline = 0,      // SigmaIEtaIEta, DEtaIn, DPhiIn, FBrem, SigmaIPhiIPhi, NBrem, 
                          // OneOverEMinusOneOverP
      kNoIPInfo,          // kBaseline + EOverP, ESeedClusterOverPout, ESeedClusterOverPIn
      kWithIPInfo,        // kV2 + d0 , IP3d, IP3dSig
      kIDIsoCombined      // ID variables , Iso03 , Iso04
    };

    enum ElectronEffectiveAreaType {
      kEleChargedIso03, 
      kEleNeutralHadronIso03, 
      kEleGammaIso03, 
      kEleGammaIsoVetoEtaStrip03, 
      kEleChargedIso04, 
      kEleNeutralHadronIso04, 
      kEleGammaIso04, 
      kEleGammaIsoVetoEtaStrip04, 
      kEleNeutralHadronIso007, 
      kEleHoverE, 
      kEleHcalDepth1OverEcal, 
      kEleHcalDepth2OverEcal    
    };
    
    void     Initialize(std::string methodName,
                        std::string Subdet0Pt10To20Weights , 
                        std::string Subdet1Pt10To20Weights , 
                        std::string Subdet2Pt10To20Weights,
                        std::string Subdet0Pt20ToInfWeights, 
                        std::string Subdet1Pt20ToInfWeights, 
                        std::string Subdet2Pt20ToInfWeights,
                        ElectronIDMVA::MVAType type );
    Bool_t   IsInitialized() const { return fIsInitialized; }

    void            SetPrintMVADebug(bool b) { fPrintMVADebug = b; }
    static Double_t ElectronEffectiveArea(ElectronEffectiveAreaType type, Double_t Eta);
    
    double MVAValue(const reco::GsfElectron *ele, const reco::Vertex vertex, 
                    EcalClusterLazyTools myEcalCluster,
                    const TransientTrackBuilder *transientTrackBuilder,
                    const reco::PFCandidateCollection &PFCandidates,
                    double Rho);
    double MVAValue(const reco::GsfElectron *ele, 
                    EcalClusterLazyTools myEcalCluster);
    double MVAValue(double ElePt , double EleSCEta,
                    double EleSigmaIEtaIEta,
                    double EleDEtaIn,
                    double EleDPhiIn,
                    double EleHoverE,
                    double EleD0,
                    double EleFBrem,
                    double EleEOverP,
                    double EleESeedClusterOverPout,
                    double EleSigmaIPhiIPhi,
                    double EleNBrem,
                    double EleOneOverEMinusOneOverP,
                    double EleESeedClusterOverPIn,
                    double EleIP3d,
                    double EleIP3dSig );
      Double_t MVAValue(Double_t ElePt , Double_t EleEta,
                        Double_t PileupEnergyDensity,
                        Double_t EleSigmaIEtaIEta,
                        Double_t EleDEtaIn,
                        Double_t EleDPhiIn,
                        Double_t EleHoverE,
                        Double_t EleD0,
                        Double_t EleDZ,
                        Double_t EleFBrem,
                        Double_t EleEOverP,
                        Double_t EleESeedClusterOverPout,
                        Double_t EleSigmaIPhiIPhi,
                        Double_t EleNBrem,
                        Double_t EleOneOverEMinusOneOverP,
                        Double_t EleESeedClusterOverPIn,
                        Double_t EleIP3d,
                        Double_t EleIP3dSig,
                        Double_t EleGsfTrackChi2OverNdof,
                        Double_t EledEtaCalo,
                        Double_t EledPhiCalo,
                        Double_t EleR9,
                        Double_t EleSCEtaWidth,
                        Double_t EleSCPhiWidth,
                        Double_t EleCovIEtaIPhi,
                        Double_t ElePreShowerOverRaw,
                        Double_t EleChargedIso03,
                        Double_t EleNeutralHadronIso03,
                        Double_t EleGammaIso03,
                        Double_t EleChargedIso04,
                        Double_t EleNeutralHadronIso04,
                        Double_t EleGammaIso04,
                        Bool_t printDebug = kFALSE );



  protected:
    TMVA::Reader             *fTMVAReader[6];
    std::string               fMethodname;
    MVAType                   fMVAType;
    
    Bool_t                    fIsInitialized;
    Bool_t                    fPrintMVADebug;
    
    Float_t                   fMVAVar_EleSigmaIEtaIEta; 
    Float_t                   fMVAVar_EleDEtaIn; 
    Float_t                   fMVAVar_EleDPhiIn; 
    Float_t                   fMVAVar_EleHoverE; 
    Float_t                   fMVAVar_EleD0; 
    Float_t                   fMVAVar_EleFBrem; 
    Float_t                   fMVAVar_EleEOverP; 
    Float_t                   fMVAVar_EleESeedClusterOverPout; 
    Float_t                   fMVAVar_EleSigmaIPhiIPhi; 
    Float_t                   fMVAVar_EleNBrem; 
    Float_t                   fMVAVar_EleOneOverEMinusOneOverP; 
    Float_t                   fMVAVar_EleESeedClusterOverPIn; 
    Float_t                   fMVAVar_EleIP3d; 
    Float_t                   fMVAVar_EleIP3dSig; 
    Float_t                   fMVAVar_EleGsfTrackChi2OverNdof;
    Float_t                   fMVAVar_EledEtaCalo;
    Float_t                   fMVAVar_EledPhiCalo;
    Float_t                   fMVAVar_EleR9;
    Float_t                   fMVAVar_EleSCEtaWidth;
    Float_t                   fMVAVar_EleSCPhiWidth;
    Float_t                   fMVAVar_EleCovIEtaIPhi;
    Float_t                   fMVAVar_ElePreShowerOverRaw;
    Float_t                   fMVAVar_EleChargedIso03OverPt;
    Float_t                   fMVAVar_EleNeutralHadronIso03OverPt;
    Float_t                   fMVAVar_EleGammaIso03OverPt;
    Float_t                   fMVAVar_EleChargedIso04OverPt;
    Float_t                   fMVAVar_EleNeutralHadronIso04OverPt;
    Float_t                   fMVAVar_EleGammaIso04OverPt;

};

#endif
