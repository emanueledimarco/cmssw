//--------------------------------------------------------------------------------------------------
// $Id $
//
// MuonIDMVA
//
// Helper Class for applying MVA muon ID selection
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef HIGGSANALYSIS_MuonIDMVA_H
#define HIGGSANALYSIS_MuonIDMVA_H

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

class MuonIDMVA {
  public:
    MuonIDMVA();
    ~MuonIDMVA(); 

    enum MVAType {
      kIDIsoCombinedDetIso
    };

    enum MuonEffectiveAreaType {
      kMuChargedIso03, 
      kMuNeutralIso03, 
      kMuChargedIso04, 
      kMuNeutralIso04, 
      kMuHadEnergy, 
      kMuHoEnergy, 
      kMuEmEnergy, 
      kMuHadS9Energy, 
      kMuHoS9Energy, 
      kMuEmS9Energy,
      kMuTrkIso03, 
      kMuEMIso03, 
      kMuHadIso03, 
      kMuTrkIso05, 
      kMuEMIso05, 
      kMuHadIso05 
    };

    void  Initialize(std::string methodName,
                     std::string Subdet0Pt10To14p5Weights , 
                     std::string Subdet1Pt10To14p5Weights , 
                     std::string Subdet0Pt14p5To20Weights,
                     std::string Subdet1Pt14p5To20Weights, 
                     std::string Subdet0Pt20ToInfWeights, 
                     std::string Subdet1Pt20ToInfWeights,
                     MuonIDMVA::MVAType type);
    Bool_t   IsInitialized() const { return fIsInitialized; }
    
    void            SetPrintMVADebug(bool b) { fPrintMVADebug = b; }
    static double MuonEffectiveArea(MuonEffectiveAreaType type, double Eta);

    double MVAValue(const reco::Muon *mu, const reco::Vertex vertex, 
                    const TransientTrackBuilder *transientTrackBuilder,
                    double Rho,
                    bool printDebug);
    
    double MVAValue( double MuPt , double MuEta,
                     double                   MuTkNchi2, 
                     double                   MuGlobalNchi2, 
                     double                   MuNValidHits, 
                     double                   MuNTrackerHits, 
                     double                   MuNPixelHits, 
                     double                   MuNMatches, 
                     double                   MuD0, 
                     double                   MuIP3d, 
                     double                   MuIP3dSig, 
                     double                   MuTrkKink, 
                     double                   MuSegmentCompatibility, 
                     double                   MuCaloCompatibility, 
                     double                   MuHadEnergyOverPt, 
                     double                   MuHoEnergyOverPt, 
                     double                   MuEmEnergyOverPt, 
                     double                   MuHadS9EnergyOverPt, 
                     double                   MuHoS9EnergyOverPt, 
                     double                   MuEmS9EnergyOverPt, 
                     double                   MuTrkIso03OverPt,
                     double                   MuEMIso03OverPt,
                     double                   MuHadIso03OverPt,
                     double                   MuTrkIso05OverPt,
                     double                   MuEMIso05OverPt,
                     double                   MuHadIso05OverPt,
                     bool                     printDebug = kFALSE
      );


  protected:      
    TMVA::Reader             *fTMVAReader[6];
    std::string               fMethodname;
    MVAType                   fMVAType;

    Bool_t                    fIsInitialized;
    Bool_t                    fPrintMVADebug;

    Float_t                   fMVAVar_MuTkNchi2; 
    Float_t                   fMVAVar_MuGlobalNchi2; 
    Float_t                   fMVAVar_MuNValidHits; 
    Float_t                   fMVAVar_MuNTrackerHits; 
    Float_t                   fMVAVar_MuNPixelHits; 
    Float_t                   fMVAVar_MuNMatches; 
    Float_t                   fMVAVar_MuD0; 
    Float_t                   fMVAVar_MuIP3d; 
    Float_t                   fMVAVar_MuIP3dSig; 
    Float_t                   fMVAVar_MuTrkKink; 
    Float_t                   fMVAVar_MuSegmentCompatibility; 
    Float_t                   fMVAVar_MuCaloCompatibility; 
    Float_t                   fMVAVar_MuHadEnergyOverPt; 
    Float_t                   fMVAVar_MuHoEnergyOverPt; 
    Float_t                   fMVAVar_MuEmEnergyOverPt; 
    Float_t                   fMVAVar_MuHadS9EnergyOverPt; 
    Float_t                   fMVAVar_MuHoS9EnergyOverPt; 
    Float_t                   fMVAVar_MuEmS9EnergyOverPt; 
    Float_t                   fMVAVar_MuChargedIso03OverPt;
    Float_t                   fMVAVar_MuNeutralIso03OverPt;
    Float_t                   fMVAVar_MuChargedIso04OverPt;
    Float_t                   fMVAVar_MuNeutralIso04OverPt;
    Float_t                   fMVAVar_MuTrkIso03OverPt;
    Float_t                   fMVAVar_MuEMIso03OverPt;
    Float_t                   fMVAVar_MuHadIso03OverPt;
    Float_t                   fMVAVar_MuTrkIso05OverPt;
    Float_t                   fMVAVar_MuEMIso05OverPt;
    Float_t                   fMVAVar_MuHadIso05OverPt;
      
};

#endif
