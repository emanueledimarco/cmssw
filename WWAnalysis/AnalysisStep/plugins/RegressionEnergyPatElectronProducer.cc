#include "WWAnalysis/AnalysisStep/interface/RegressionEnergyPatElectronProducer.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include <iostream>

// get ROOT_VERSION_CODE
#include <TROOT.h> 


using namespace edm ;
using namespace std ;
using namespace reco ;
using namespace pat ;

RegressionEnergyPatElectronProducerCalifornia::RegressionEnergyPatElectronProducerCalifornia( const edm::ParameterSet & cfg )
{

  produces<ElectronCollection>();

  inputPatElectrons = cfg.getParameter<edm::InputTag>("inputPatElectronsTag");
  rhoInputTag_ = cfg.getParameter<edm::InputTag>("rhoCollection");
  energyRegressionType = cfg.getParameter<uint32_t>("energyRegressionType");
  regressionInputFile = cfg.getParameter<std::string>("regressionInputFile");
  debug = cfg.getUntrackedParameter<bool>("debug");


  //****************************************************************************************
  //set up regression calculator
  //****************************************************************************************
  regressionEvaluator = new ElectronEnergyRegressionEvaluate();

  //set regression type
  ElectronEnergyRegressionEvaluate::ElectronEnergyRegressionType type = ElectronEnergyRegressionEvaluate::kNoTrkVar;
  if (energyRegressionType == 1) type = ElectronEnergyRegressionEvaluate::kNoTrkVar;
  else if (energyRegressionType == 2) type = ElectronEnergyRegressionEvaluate::kWithTrkVarV1;
  else if (energyRegressionType == 3) type = ElectronEnergyRegressionEvaluate::kWithTrkVarV2;

  //load weights and initialize
  regressionEvaluator->initialize(regressionInputFile.c_str(),type);

}
 
RegressionEnergyPatElectronProducerCalifornia::~RegressionEnergyPatElectronProducerCalifornia()
{
  delete regressionEvaluator;
}

void RegressionEnergyPatElectronProducerCalifornia::produce( edm::Event & event, const edm::EventSetup & setup )
{

  assert(regressionEvaluator->isInitialized());


  //**************************************************************************
  //Get Number of Vertices
  //**************************************************************************
  Handle<reco::VertexCollection> hVertexProduct;
  event.getByLabel(edm::InputTag("goodPrimaryVertices"),hVertexProduct);
  const reco::VertexCollection inVertices = *(hVertexProduct.product());  

  // loop through all vertices
  Int_t nvertices = 0;
  for (reco::VertexCollection::const_iterator inV = inVertices.begin(); 
       inV != inVertices.end(); ++inV) {

    //pass these vertex cuts
    if (inV->ndof() >= 4
        && inV->position().Rho() <= 2.0
        && fabs(inV->z()) <= 24.0
      ) {
      nvertices++;
    }
  }

  //**************************************************************************
  //Get Rho
  //**************************************************************************
  double rho = 0;
  Handle<double> hRhoKt6PFJets;
  event.getByLabel(rhoInputTag_, hRhoKt6PFJets);
  rho = (*hRhoKt6PFJets);

  edm::Handle<edm::View<reco::Candidate> > oldElectrons ;
  event.getByLabel(inputPatElectrons,oldElectrons) ;

  std::auto_ptr<ElectronCollection> electrons( new ElectronCollection ) ;

  ElectronCollection::const_iterator electron ;
  ElectronCollection::iterator ele ;
  // first clone the initial collection
  for(edm::View<reco::Candidate>::const_iterator ele=oldElectrons->begin(); ele!=oldElectrons->end(); ++ele) {    
    const pat::ElectronRef elecsRef = edm::RefToBase<reco::Candidate>(oldElectrons,ele-oldElectrons->begin()).castTo<pat::ElectronRef>();
    pat::Electron clone = *edm::RefToBase<reco::Candidate>(oldElectrons,ele-oldElectrons->begin()).castTo<pat::ElectronRef>();
    electrons->push_back(clone);
  }

  for
    ( ele = electrons->begin() ;
      ele != electrons->end() ;
      ++ele ) 
  {     


    bool printDebug = debug;
    if (printDebug) {
      std::cout << "***********************************************************************\n";
      std::cout << "Run Lumi Event: " << event.id().run() << " " << event.luminosityBlock() << " " << event.id().event() << "\n";
      std::cout << "Pat Electron : " << ele->pt() << " " << ele->eta() << " " << ele->phi() << "\n";
    }

    // apply regression energy
    Double_t FinalMomentum = 0;
    Double_t FinalMomentumError = 0;
    Double_t RegressionMomentum = 0;
    Double_t RegressionMomentumError = 0;

    //compute spp and sep
    double spp = 0;
    if (!isnan(ele->userFloat("covIPhiIPhi"))) spp = sqrt (ele->userFloat("covIPhiIPhi"));
    double sep;
    if (ele->sigmaIetaIeta()*spp > 0) {
      sep = ele->userFloat("covIEtaIPhi") / (ele->sigmaIetaIeta()*spp);
    } else if (ele->userFloat("covIEtaIPhi") > 0) {
      sep = 1.0;
    } else {
      sep = -1.0;
    }

    if (energyRegressionType == 1) {
      RegressionMomentum = regressionEvaluator->regressionValueNoTrkVar( ele->userFloat("rawEnergy"),
                                                                         ele->userFloat("eta"),
                                                                         ele->userFloat("phi"),
                                                                         ele->userFloat("e3x3") / ele->userFloat("rawEnergy"),
                                                                         ele->userFloat("etaWidth"),
                                                                         ele->userFloat("phiWidth"),
                                                                         ele->userInt("nBC"),
                                                                         ele->hadronicOverEm(),
                                                                         rho, 
                                                                         nvertices, 
                                                                         ele->userFloat("seedClusterEta"),
                                                                         ele->userFloat("seedClusterPhi"),
                                                                         ele->userFloat("seedClusterEnergy"),
                                                                         ele->userFloat("e3x3"),
                                                                         ele->userFloat("e5x5"),
                                                                         ele->sigmaIetaIeta(),
                                                                         spp,
                                                                         sep,
                                                                         ele->userFloat("eMax"),
                                                                         ele->userFloat("e2nd"),
                                                                         ele->userFloat("eTop"),
                                                                         ele->userFloat("eBottom"),
                                                                         ele->userFloat("eLeft"),
                                                                         ele->userFloat("eRight"),
                                                                         ele->userFloat("e2x5Max"),
                                                                         ele->userFloat("e2x5Top"),
                                                                         ele->userFloat("e2x5Bottom"),
                                                                         ele->userFloat("e2x5Left"),
                                                                         ele->userFloat("e2x5Right"),
                                                                         ele->userFloat("seedieta"),
                                                                         ele->userFloat("seediphi"),
                                                                         ele->userFloat("etacryseed"),
                                                                         ele->userFloat("phicryseed"),
                                                                         ele->userFloat("esEnergy")/ele->userFloat("rawEnergy"),
                                                                         printDebug);
      RegressionMomentumError = regressionEvaluator->regressionUncertaintyNoTrkVar( ele->userFloat("rawEnergy"),
                                                                                    ele->userFloat("eta"),
                                                                                    ele->userFloat("phi"),
                                                                                    ele->userFloat("e3x3") / ele->userFloat("rawEnergy"),
                                                                                    ele->userFloat("etaWidth"),
                                                                                    ele->userFloat("phiWidth"),
                                                                                    ele->userInt("nBC"),
                                                                                    ele->hadronicOverEm(),
                                                                                    rho, 
                                                                                    nvertices, 
                                                                                    ele->userFloat("seedClusterEta"),
                                                                                    ele->userFloat("seedClusterPhi"),
                                                                                    ele->userFloat("seedClusterEnergy"),
                                                                                    ele->userFloat("e3x3"),
                                                                                    ele->userFloat("e5x5"),
                                                                                    ele->sigmaIetaIeta(),
                                                                                    spp,
                                                                                    sep,
                                                                                    ele->userFloat("eMax"),
                                                                                    ele->userFloat("e2nd"),
                                                                                    ele->userFloat("eTop"),
                                                                                    ele->userFloat("eBottom"),
                                                                                    ele->userFloat("eLeft"),
                                                                                    ele->userFloat("eRight"),
                                                                                    ele->userFloat("e2x5Max"),
                                                                                    ele->userFloat("e2x5Top"),
                                                                                    ele->userFloat("e2x5Bottom"),
                                                                                    ele->userFloat("e2x5Left"),
                                                                                    ele->userFloat("e2x5Right"),
                                                                                    ele->userFloat("seedieta"),
                                                                                    ele->userFloat("seediphi"),
                                                                                    ele->userFloat("etacryseed"),
                                                                                    ele->userFloat("phicryseed"),
                                                                                    ele->userFloat("esEnergy")/ele->userFloat("rawEnergy"),
                                                                                    printDebug);
      
      //Combine with track momentum measurement      
      int elClass = ele->classification() ;
      float TrackMomentum  = ele->trackMomentumAtVtx().R() ;
      float TrackMomentumError = ele->trackMomentumError();

      //initial      
      if (TrackMomentumError/TrackMomentum > 0.5 && RegressionMomentumError/RegressionMomentum <= 0.5) {
        FinalMomentum = RegressionMomentum;    FinalMomentumError = RegressionMomentumError;
      }
      else if (TrackMomentumError/TrackMomentum <= 0.5 && RegressionMomentumError/RegressionMomentum > 0.5){
        FinalMomentum = TrackMomentum;  FinalMomentumError = TrackMomentumError;
      }
      else if (TrackMomentumError/TrackMomentum > 0.5 && RegressionMomentumError/RegressionMomentum > 0.5) {
        if (TrackMomentumError/TrackMomentum < RegressionMomentumError/RegressionMomentum) {
          FinalMomentum = TrackMomentum; FinalMomentumError = TrackMomentumError;
        }
        else{
          FinalMomentum = RegressionMomentum; FinalMomentumError = RegressionMomentumError;
        }
      }
      // then apply the combination algorithm
      else {

        // calculate E/p and corresponding error
        float eOverP = RegressionMomentum / TrackMomentum;
        float errorEOverP = sqrt(
          (RegressionMomentumError/TrackMomentum)*(RegressionMomentumError/TrackMomentum) +
          (RegressionMomentum*TrackMomentumError/TrackMomentum/TrackMomentum)*
          (RegressionMomentum*TrackMomentumError/TrackMomentum/TrackMomentum));

#if ROOT_VERSION_CODE >=  ROOT_VERSION(5,30,00)
	// ======== HERE IS THE 5XY (NEW) E-P COMBINATION	
        bool eleIsNotInCombination = false ;
        if ( (eOverP  > 1 + 2.5*errorEOverP) || (eOverP  < 1 - 2.5*errorEOverP) || (eOverP < 0.8) || (eOverP > 1.3) )
        { eleIsNotInCombination = true ; }
        if (eleIsNotInCombination)
        {
          if (eOverP > 1)
          { FinalMomentum = RegressionMomentum ; FinalMomentumError = RegressionMomentumError ; }
          else
          {
            if (elClass == reco::GsfElectron::GOLDEN)
            { FinalMomentum = RegressionMomentum; FinalMomentumError = RegressionMomentumError; }
            if (elClass == reco::GsfElectron::BIGBREM)
            {
              if (RegressionMomentum<36)
              { FinalMomentum = TrackMomentum ; FinalMomentumError = TrackMomentumError ; }
              else
              { FinalMomentum = RegressionMomentum ; FinalMomentumError = RegressionMomentumError ; }
            }
            if (elClass == reco::GsfElectron::BADTRACK) 
            { FinalMomentum = RegressionMomentum; FinalMomentumError = RegressionMomentumError ; }
            if (elClass == reco::GsfElectron::SHOWERING)
            {
              if (RegressionMomentum<30)
              { FinalMomentum = TrackMomentum ; FinalMomentumError = TrackMomentumError; }
              else
              { FinalMomentum = RegressionMomentum; FinalMomentumError = RegressionMomentumError;}
            }
            if (elClass == reco::GsfElectron::GAP)
            {
              if (RegressionMomentum<60)
              { FinalMomentum = TrackMomentum ; FinalMomentumError = TrackMomentumError ; }
              else
              { FinalMomentum = RegressionMomentum; FinalMomentumError = RegressionMomentumError ; }
            }
          }
        }
 
        else
        {
          // combination
          FinalMomentum = (RegressionMomentum/RegressionMomentumError/RegressionMomentumError + TrackMomentum/TrackMomentumError/TrackMomentumError) /
            (1/RegressionMomentumError/RegressionMomentumError + 1/TrackMomentumError/TrackMomentumError);
          float FinalMomentumVariance = 1 / (1/RegressionMomentumError/RegressionMomentumError + 1/TrackMomentumError/TrackMomentumError);
          FinalMomentumError = sqrt(FinalMomentumVariance);
        }
	// ======== FINISH 5XY E-P COMBINATION
#else
	// ======== HERE IS THE 4XY (NEW) E-P COMBINATION
	if ( eOverP  > 1 + 2.5*errorEOverP )
	  {
	    FinalMomentum = RegressionMomentum; FinalMomentumError = RegressionMomentumError;
	    if ((elClass==reco::GsfElectron::GOLDEN) && ele->isEB() && (eOverP<1.15))
	      {
		if (RegressionMomentum<15) {FinalMomentum = TrackMomentum ; FinalMomentumError = TrackMomentumError;}
	      }
	  }
	else if ( eOverP < 1 - 2.5*errorEOverP )
	  {
	    FinalMomentum = RegressionMomentum; FinalMomentumError = RegressionMomentumError;
	    if (elClass==reco::GsfElectron::SHOWERING)
	      {
		if (ele->isEB())
		  {
		    if(RegressionMomentum<18) {FinalMomentum = TrackMomentum; FinalMomentumError = TrackMomentumError;}
		  }
		else if (ele->isEE())
		  {
		    if(RegressionMomentum<13) {FinalMomentum = TrackMomentum; FinalMomentumError = TrackMomentumError;}
		  }
		else
		  { edm::LogWarning("ElectronMomentumCorrector::correct")<<"nor barrel neither endcap electron ?!" ; }
	      }
	    else if (ele->isGap())
	      {
		if(RegressionMomentum<60) {FinalMomentum = TrackMomentum; FinalMomentumError = TrackMomentumError;}
	      }
	  }
	else 
	  {
	    // combination
	    FinalMomentum = (RegressionMomentum/RegressionMomentumError/RegressionMomentumError + TrackMomentum/TrackMomentumError/TrackMomentumError) /
	      (1/RegressionMomentumError/RegressionMomentumError + 1/TrackMomentumError/TrackMomentumError);
	    float FinalMomentumVariance = 1 / (1/RegressionMomentumError/RegressionMomentumError + 1/TrackMomentumError/TrackMomentumError);
	    FinalMomentumError = sqrt(FinalMomentumVariance);
	  }
	// ======== FINISH 4XY E-P COMBINATION
#endif
      }

      if(printDebug) {
        cout << "E-P combination: " << RegressionMomentum << " , " << TrackMomentum << " --> " << FinalMomentum << "\n";
      }

    } else if (energyRegressionType == 2) {
      RegressionMomentum = regressionEvaluator->regressionValueWithTrkVarV1( 
                                                                     ele->userFloat("rawEnergy"),
                                                                     ele->userFloat("eta"),
                                                                     ele->userFloat("phi"),
                                                                     ele->userFloat("e3x3") / ele->userFloat("rawEnergy"),
                                                                     ele->userFloat("etaWidth"),
                                                                     ele->userFloat("phiWidth"),
                                                                     ele->userInt("nBC"),
                                                                     ele->hadronicOverEm(),
                                                                     rho, 
                                                                     nvertices, 
                                                                     ele->userFloat("seedClusterEta"),
                                                                     ele->userFloat("seedClusterPhi"),
                                                                     ele->userFloat("seedClusterEnergy"),
                                                                     ele->userFloat("e3x3"),
                                                                     ele->userFloat("e5x5"),
                                                                     ele->sigmaIetaIeta(),
                                                                     spp,
                                                                     sep,
                                                                     ele->userFloat("eMax"),
                                                                     ele->userFloat("e2nd"),
                                                                     ele->userFloat("eTop"),
                                                                     ele->userFloat("eBottom"),
                                                                     ele->userFloat("eLeft"),
                                                                     ele->userFloat("eRight"),
                                                                     ele->userFloat("e2x5Max"),
                                                                     ele->userFloat("e2x5Top"),
                                                                     ele->userFloat("e2x5Bottom"),
                                                                     ele->userFloat("e2x5Left"),
                                                                     ele->userFloat("e2x5Right"),
                                                                     ele->userFloat("seedieta"),
                                                                     ele->userFloat("seediphi"),
                                                                     ele->userFloat("etacryseed"),
                                                                     ele->userFloat("phicryseed"),
                                                                     ele->userFloat("esEnergy")/ele->userFloat("rawEnergy"),
                                                                     ele->ecalDrivenSeed(),
                                                                     ele->trackMomentumAtVtx().R(),
                                                                     fmax(ele->fbrem(),-1.0),
                                                                     ele->charge(),
                                                                     fmin(ele->eSuperClusterOverP(), 20.0),
                                                                     ele->trackMomentumError(),
#if ROOT_VERSION_CODE >=  ROOT_VERSION(5,30,00)
                                                                     ele->correctedEcalEnergyError(),
#else
                                                                     ele->ecalEnergyError(),
#endif

                                                                     ele->classification(),
                                                                     printDebug);
      RegressionMomentumError = regressionEvaluator->regressionUncertaintyWithTrkVarV1(
                                                                                ele->userFloat("rawEnergy"),
                                                                                ele->userFloat("eta"),
                                                                                ele->userFloat("phi"),
                                                                                ele->userFloat("e3x3") / ele->userFloat("rawEnergy"),
                                                                                ele->userFloat("etaWidth"),
                                                                                ele->userFloat("phiWidth"),
                                                                                ele->userInt("nBC"),
                                                                                ele->hadronicOverEm(),
                                                                                rho, 
                                                                                nvertices, 
                                                                                ele->userFloat("seedClusterEta"),
                                                                                ele->userFloat("seedClusterPhi"),
                                                                                ele->userFloat("seedClusterEnergy"),
                                                                                ele->userFloat("e3x3"),
                                                                                ele->userFloat("e5x5"),
                                                                                ele->sigmaIetaIeta(),
                                                                                spp,
                                                                                sep,
                                                                                ele->userFloat("eMax"),
                                                                                ele->userFloat("e2nd"),
                                                                                ele->userFloat("eTop"),
                                                                                ele->userFloat("eBottom"),
                                                                                ele->userFloat("eLeft"),
                                                                                ele->userFloat("eRight"),
                                                                                ele->userFloat("e2x5Max"),
                                                                                ele->userFloat("e2x5Top"),
                                                                                ele->userFloat("e2x5Bottom"),
                                                                                ele->userFloat("e2x5Left"),
                                                                                ele->userFloat("e2x5Right"),
                                                                                ele->userFloat("seedieta"),
                                                                                ele->userFloat("seediphi"),
                                                                                ele->userFloat("etacryseed"),
                                                                                ele->userFloat("phicryseed"),
                                                                                ele->userFloat("esEnergy")/ele->userFloat("rawEnergy"),
                                                                                ele->ecalDrivenSeed(),
                                                                                ele->trackMomentumAtVtx().R(),
                                                                                fmax(ele->fbrem(),-1.0),
                                                                                ele->charge(),
                                                                                fmin(ele->eSuperClusterOverP(), 20.0),
                                                                                ele->trackMomentumError(),
#if ROOT_VERSION_CODE >=  ROOT_VERSION(5,30,00)
                                                                     ele->correctedEcalEnergyError(),
#else
                                                                     ele->ecalEnergyError(),
#endif
                                                                                ele->classification(),
                                                                                printDebug);
      FinalMomentum = RegressionMomentum;
      FinalMomentumError = RegressionMomentumError;
    } else if (energyRegressionType == 3) {

      RegressionMomentum = regressionEvaluator->regressionValueWithTrkVarV2( 
                                                                     ele->userFloat("rawEnergy"),
                                                                     ele->userFloat("eta"),
                                                                     ele->userFloat("phi"),
                                                                     ele->userFloat("e3x3") / ele->userFloat("rawEnergy"),
                                                                     ele->userFloat("etaWidth"),
                                                                     ele->userFloat("phiWidth"),
                                                                     ele->userInt("nBC"),
                                                                     ele->hadronicOverEm(),
                                                                     rho, 
                                                                     nvertices, 
                                                                     ele->userFloat("seedClusterEta"),
                                                                     ele->userFloat("seedClusterPhi"),
                                                                     ele->userFloat("seedClusterEnergy"),
                                                                     ele->userFloat("e3x3"),
                                                                     ele->userFloat("e5x5"),
                                                                     ele->sigmaIetaIeta(),
                                                                     spp,
                                                                     sep,
                                                                     ele->userFloat("eMax"),
                                                                     ele->userFloat("e2nd"),
                                                                     ele->userFloat("eTop"),
                                                                     ele->userFloat("eBottom"),
                                                                     ele->userFloat("eLeft"),
                                                                     ele->userFloat("eRight"),
                                                                     ele->userFloat("e2x5Max"),
                                                                     ele->userFloat("e2x5Top"),
                                                                     ele->userFloat("e2x5Bottom"),
                                                                     ele->userFloat("e2x5Left"),
                                                                     ele->userFloat("e2x5Right"),
                                                                     ele->userFloat("seedieta"),
                                                                     ele->userFloat("seediphi"),
                                                                     ele->userFloat("etacryseed"),
                                                                     ele->userFloat("phicryseed"),
                                                                     ele->userFloat("esEnergy")/ele->userFloat("rawEnergy"),
                                                                     ele->ecalDrivenSeed(),
                                                                     ele->trackMomentumAtVtx().R(),
                                                                     fmax(ele->fbrem(),-1.0),
                                                                     ele->charge(),
                                                                     fmin(ele->eSuperClusterOverP(), 20.0),
                                                                     ele->trackMomentumError(),
#if ROOT_VERSION_CODE >=  ROOT_VERSION(5,30,00)
                                                                     ele->correctedEcalEnergyError(),
#else
                                                                     ele->ecalEnergyError(),
#endif
                                                                     ele->classification(),
                                                                     fmin(fabs(ele->deltaEtaSuperClusterTrackAtVtx()), 0.6),
                                                                     ele->deltaPhiSuperClusterTrackAtVtx(),
                                                                     ele->deltaEtaSeedClusterTrackAtCalo(),
                                                                     ele->deltaPhiSeedClusterTrackAtCalo(),
                                                                     ele->gsfTrack()->chi2() / ele->gsfTrack()->ndof(),
                                                                     5, //(ele->closestCtfTrackRef().isNonnull() ? ele->closestCtfTrackRef()->hitPattern().trackerLayersWithMeasurement() : -1), 
                                                                     fmin(ele->eEleClusterOverPout(),20.0),
                                                                     printDebug);
      RegressionMomentumError = regressionEvaluator->regressionUncertaintyWithTrkVarV2(
                                                                                ele->userFloat("rawEnergy"),
                                                                                ele->userFloat("eta"),
                                                                                ele->userFloat("phi"),
                                                                                ele->userFloat("e3x3") / ele->userFloat("rawEnergy"),
                                                                                ele->userFloat("etaWidth"),
                                                                                ele->userFloat("phiWidth"),
                                                                                ele->userInt("nBC"),
                                                                                ele->hadronicOverEm(),
                                                                                rho, 
                                                                                nvertices, 
                                                                                ele->userFloat("seedClusterEta"),
                                                                                ele->userFloat("seedClusterPhi"),
                                                                                ele->userFloat("seedClusterEnergy"),
                                                                                ele->userFloat("e3x3"),
                                                                                ele->userFloat("e5x5"),
                                                                                ele->sigmaIetaIeta(),
                                                                                spp,
                                                                                sep,
                                                                                ele->userFloat("eMax"),
                                                                                ele->userFloat("e2nd"),
                                                                                ele->userFloat("eTop"),
                                                                                ele->userFloat("eBottom"),
                                                                                ele->userFloat("eLeft"),
                                                                                ele->userFloat("eRight"),
                                                                                ele->userFloat("e2x5Max"),
                                                                                ele->userFloat("e2x5Top"),
                                                                                ele->userFloat("e2x5Bottom"),
                                                                                ele->userFloat("e2x5Left"),
                                                                                ele->userFloat("e2x5Right"),
                                                                                ele->userFloat("seedieta"),
                                                                                ele->userFloat("seediphi"),
                                                                                ele->userFloat("etacryseed"),
                                                                                ele->userFloat("phicryseed"),
                                                                                ele->userFloat("esEnergy")/ele->userFloat("rawEnergy"),
                                                                                ele->ecalDrivenSeed(),
                                                                                ele->trackMomentumAtVtx().R(),
                                                                                fmax(ele->fbrem(),-1.0),
                                                                                ele->charge(),
                                                                                fmin(ele->eSuperClusterOverP(), 20.0),
                                                                                ele->trackMomentumError(),
#if ROOT_VERSION_CODE >=  ROOT_VERSION(5,30,00)
                                                                     ele->correctedEcalEnergyError(),
#else
                                                                     ele->ecalEnergyError(),
#endif
                                                                                ele->classification(),
                                                                                fmin(fabs(ele->deltaEtaSuperClusterTrackAtVtx()), 0.6),
                                                                                ele->deltaPhiSuperClusterTrackAtVtx(),
                                                                                ele->deltaEtaSeedClusterTrackAtCalo(),
                                                                                ele->deltaPhiSeedClusterTrackAtCalo(),
                                                                                ele->gsfTrack()->chi2() / ele->gsfTrack()->ndof(),
                                                                                5, // (ele->closestCtfTrackRef().isNonnull() ? ele->closestCtfTrackRef()->hitPattern().trackerLayersWithMeasurement() : -1), 
                                                                                fmin(ele->eEleClusterOverPout(),20.0),
                                                                                printDebug);
      FinalMomentum = RegressionMomentum;
      FinalMomentumError = RegressionMomentumError;
    } else {
      cout << "Error: RegressionType = " << energyRegressionType << " is not supported.\n";
    }
    
    math::XYZTLorentzVector oldMomentum = ele->p4();
    math::XYZTLorentzVector newMomentum = math::XYZTLorentzVector
      ( oldMomentum.x()*FinalMomentum/oldMomentum.t(),
        oldMomentum.y()*FinalMomentum/oldMomentum.t(),
        oldMomentum.z()*FinalMomentum/oldMomentum.t(),
        FinalMomentum ) ;


    if (ele->core()->ecalDrivenSeed()) {              
      ele->correctEcalEnergy(RegressionMomentum, RegressionMomentumError);
      ele->correctMomentum(newMomentum,ele->trackMomentumError(),FinalMomentumError);
    }

  }
  event.put(electrons) ;
     
}
  
  
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
DEFINE_FWK_MODULE(RegressionEnergyPatElectronProducerCalifornia);
