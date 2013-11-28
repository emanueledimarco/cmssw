#include "HiggsAnalysis/HiggsToWW2Leptons/plugins/HWWPreselectionMarker.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include <iostream>

HWWPreselectionMarker::HWWPreselectionMarker(const edm::ParameterSet& iConfig)
{
	using namespace edm;
	using namespace reco;
	
	muonslabel_=iConfig.getParameter<InputTag>( "MuonLabel");
	electronslabel_=iConfig.getParameter<InputTag>( "ElectronLabel");
	metlabel_=iConfig.getParameter<InputTag>( "MetLabel");
	jetslabel_=iConfig.getParameter<InputTag>( "JetLabel");
	
	leptonPtMinMin_=iConfig.getParameter<double>("LeptonPtMinMin");
	leptonPtMaxMin_=iConfig.getParameter<double>("LeptonPtMaxMin");        
	leptonEtaMax_=iConfig.getParameter<double>("LeptonEtaMax");      
	leptonChargeCombination_=iConfig.getParameter<double>("LeptonChargeCombination");  
	metMin_=iConfig.getParameter<double>("MetMin");  
 	invMassMin_=iConfig.getParameter<double>("InvMassMin");      
	produces<bool>();
}

HWWPreselectionMarker::~HWWPreselectionMarker()
{
}

void HWWPreselectionMarker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{       
        using namespace edm;
        using namespace reco;
	using namespace std;

        Handle<MuonCollection> muons;
        Handle<GsfElectronCollection> electrons;
        Handle<METCollection> met;
        Handle<CaloJetCollection> jets;		
        std::auto_ptr<bool> pOut(new bool);
        *pOut=true;
        
        iEvent.getByLabel(muonslabel_, muons);
        iEvent.getByLabel(electronslabel_, electrons);
        iEvent.getByLabel(metlabel_, met);
        iEvent.getByLabel(jetslabel_, jets);
        std::vector<const RecoCandidate*> leptons;

	//  ----- Lepton selection 
        GsfElectronCollection::const_iterator electron;
        for(electron=electrons->begin(); electron != electrons->end(); electron++)
        {
	  if(electron->pt()>=leptonPtMinMin_ && fabs(electron->eta())<=leptonEtaMax_)
	     {
                        const RecoCandidate *lepton=&(*electron);
                        leptons.push_back(lepton);
	     }
        }
        MuonCollection::const_iterator muon;
        for(muon=muons->begin(); muon != muons->end(); muon++)
        {
	  if(muon->pt()>= leptonPtMinMin_&& fabs(muon->eta())<=leptonEtaMax_)
	    {
                        const RecoCandidate *lepton=&(*muon);
                        leptons.push_back(lepton);
	    }
        }
	// event counter
	selectedEvents[0]++;
  	// find hardest lepton
        if(leptons.size()<2)
                *pOut=false;
        else
	// 2 good leptons found
        { 

	  selectedEvents[1]++;

	  const RecoCandidate *lepton1 = leptons[0];
	  const RecoCandidate *lepton2 = leptons[1];

          //choose hardest leptons as lepton1 and lepton
	  if(leptons.size() > 2) {
	    std::vector<const RecoCandidate*>::const_iterator lepton;
	    for(lepton=leptons.begin(); lepton != leptons.end(); lepton++){
	      if((*lepton)->pt() > lepton1->pt()) {
		if((*lepton)->pt() > lepton2->pt()) {
		  if(lepton1->pt() > lepton2->pt() && (*lepton)->pt()!=lepton1->pt()) lepton2 =(*lepton);
		  else	if((*lepton)->pt() != lepton2->pt()) lepton1 = (*lepton);
		}
		else if ((*lepton)->pt() != lepton2->pt()) lepton1 = (*lepton);
	      }
	      else if((*lepton)->pt() > lepton2->pt() && (*lepton)->pt()!= lepton1->pt())
		lepton2 = *lepton;
	    }
	  }

	  // one lepton with pt>20 GeV
	  if(lepton1->pt() <leptonPtMaxMin_  && lepton2->pt()<leptonPtMaxMin_) 
	    *pOut=false;
	  else selectedEvents[2]++;

	  // unlike signed leptons
          if(lepton1->charge()*lepton2->charge() !=leptonChargeCombination_)
	    *pOut=false;
	  else selectedEvents[3]++;                
	  
	  //met cut  
	  reco::METCollection::const_iterator met_iter;
	  for(met_iter=met->begin(); met_iter != met->end(); met_iter++)
	    {
	      if(met_iter->pt() < metMin_)
		*pOut=false;
	      else selectedEvents[4]++;
	    }

	  // invariant mass cut (lepton mass is neglected)               
	  double entot = lepton1->energy()+lepton2->energy();
	  double pxtot = lepton1->px()+lepton2->px();
	  double pytot = lepton1->py()+lepton2->py();
	  double pztot = lepton1->pz()+lepton2->pz();
	  double inmass = entot*entot - pxtot*pxtot - pytot*pytot - pztot*pztot;
	  double inmass_sqrt=sqrt(inmass);
	  if(inmass_sqrt<invMassMin_)*pOut=false;
	  else selectedEvents[5]++;
	}
        if(*pOut==true) selectedEvents[6]++;
        iEvent.put(pOut);
}

void HWWPreselectionMarker::beginJob(const edm::EventSetup&)
{
  for(int i=0;i<7;i++) selectedEvents[i]=0;
}

void HWWPreselectionMarker::endJob() {
  std::string cutnames[7];
  cutnames[0]="# # tot. events          : ";
  cutnames[1]="# of Leptons > 1         : ";
  cutnames[2]="one lepton with pt>20GeV : ";
  cutnames[3]="Lepton charge product    : ";
  cutnames[4]="MET                      : ";
  cutnames[5]="Inv Mass                 : ";
  cutnames[6]="All cuts                 : ";
  std::cout<<" =================================== "<<std::endl;
  for(int i=0; i<7; ++i) 
	std::cout<<" #"<<i <<" "<< cutnames[i] << selectedEvents[i]<<std::endl;
  std::cout<<" =================================== "<<std::endl;}

//define this as a plug-in
DEFINE_FWK_MODULE(HWWPreselectionMarker);
