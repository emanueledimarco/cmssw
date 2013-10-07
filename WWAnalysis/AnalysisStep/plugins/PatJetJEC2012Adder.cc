#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <DataFormats/PatCandidates/interface/Jet.h>
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CommonTools/Utils/interface/PtComparator.h"
#include<vector>
#include<TVector3.h>


class PatJetJEC2012Adder : public edm::EDProducer {
    public:
        explicit PatJetJEC2012Adder(const edm::ParameterSet&);
        ~PatJetJEC2012Adder();

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

	bool 		isData_;
        edm::InputTag jetTag_;
        edm::InputTag rhoTag_;
};

PatJetJEC2012Adder::PatJetJEC2012Adder(const edm::ParameterSet& iConfig) :
    jetTag_(iConfig.getParameter<edm::InputTag>("jetTag")),
    rhoTag_(iConfig.getParameter<edm::InputTag>("rhoTag"))
{
    isData_= iConfig.getUntrackedParameter<bool>("isData",false);
    produces<pat::JetCollection>();  
}




void PatJetJEC2012Adder::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    edm::Handle<edm::View<reco::Candidate> > jetH;
    iEvent.getByLabel(jetTag_,jetH);

    edm::Handle<double> hRho;
    iEvent.getByLabel(rhoTag_,hRho);
    float rho = *hRho;

    std::auto_ptr<pat::JetCollection> pOut(new pat::JetCollection);

  const char *base=getenv("CMSSW_BASE");
  std::string baseFolder(base);
  baseFolder += "/src/WWAnalysis/AnalysisStep/data/";
   // New 2012 JEC
   // Create the JetCorrectorParameter objects, the order does not matter.
        JetCorrectorParameters *ResJetPar,*L3JetPar,*L2JetPar,*L1JetPar;
        std::vector<JetCorrectorParameters> vPar;
        if (isData_){
   ResJetPar = new JetCorrectorParameters(baseFolder+"GR_P_V42_AN3_L2L3Residual_AK5PF.txt");
   L3JetPar  = new JetCorrectorParameters(baseFolder+"GR_P_V42_AN3_L3Absolute_AK5PF.txt");
   L2JetPar  = new JetCorrectorParameters(baseFolder+"GR_P_V42_AN3_L2Relative_AK5PF.txt");
   L1JetPar  = new JetCorrectorParameters(baseFolder+"GR_P_V42_AN3_L1FastJet_AK5PF.txt");

   vPar.push_back(*L1JetPar);
   vPar.push_back(*L2JetPar);
   vPar.push_back(*L3JetPar);
   vPar.push_back(*ResJetPar);
        }
        else{
   L3JetPar  = new JetCorrectorParameters(baseFolder+"START53_V15_L3Absolute_AK5PF.txt");
   L2JetPar  = new JetCorrectorParameters(baseFolder+"START53_V15_L2Relative_AK5PF.txt");
   L1JetPar  = new JetCorrectorParameters(baseFolder+"START53_V15_L1FastJet_AK5PF.txt");
   vPar.push_back(*L1JetPar);
   vPar.push_back(*L2JetPar);
   vPar.push_back(*L3JetPar);
        }

   FactorizedJetCorrector *JetCorrector = new FactorizedJetCorrector(vPar);
	//int count=0;
    
//    for(edm::View<pat::Jet>::const_iterator itJet=jetH->begin(); itJet!=jetH->end(); ++itJet) {    
    for(edm::View<reco::Candidate>::const_iterator itJet=jetH->begin(); itJet!=jetH->end(); ++itJet) {

//        pat::Jet clone = *edm::RefToBase<pat::Jet>(itJet).castTo<pat::JetRef>();

        pat::Jet clone = *edm::RefToBase<reco::Candidate>(jetH,itJet-jetH->begin()).castTo<pat::JetRef>();

	//std::cout<<"orig pt: "<<(*itJet).pt()<<" clone: "<<clone.pt()<<std::endl;
	//std::cout<<"count: "<<++count<<std::endl;
//        pat::JetRef jetRef = edm::RefToBase<reco::Candidate>(jetH,itJet-jetH->begin()).castTo<pat::JetRef>();

   double correction(1.0);
   reco::Jet mijetRAW = clone.correctedJet("Uncorrected");
   JetCorrector->setJetEta(clone.eta());
   JetCorrector->setJetPt(mijetRAW.pt());
   JetCorrector->setJetA(clone.jetArea());
   JetCorrector->setRho(rho);
//   cout <<"et: "<<mijetRAW.et()<<" pt: "<<mijetRAW.pt()<<endl;

    correction = JetCorrector->getCorrection();

	clone.setP4(mijetRAW.p4());
        clone.scaleEnergy (correction); 
        pOut->push_back(clone);
//	std::cout<<"correction: "<<correction<<std::endl;
    }

 NumericSafeGreaterByPt<pat::Jet> compJets;
 // reorder corrected jets
  std::sort (pOut->begin (), pOut->end (), compJets);
 // put corrected jet collection into event
    iEvent.put(pOut);

 delete JetCorrector;
 delete L3JetPar;
 delete L2JetPar;
 delete L1JetPar;
        if(isData_) delete ResJetPar;

}

PatJetJEC2012Adder::~PatJetJEC2012Adder() { }
void PatJetJEC2012Adder::beginJob() { }
void PatJetJEC2012Adder::endJob() { }
DEFINE_FWK_MODULE(PatJetJEC2012Adder);
