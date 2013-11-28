#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include<vector>

class PatElectronLHAdder : public edm::EDProducer {
    public:
        explicit PatElectronLHAdder(const edm::ParameterSet&);
        ~PatElectronLHAdder();

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        edm::InputTag electronTag_;
        std::string lhFile_;
};

PatElectronLHAdder::PatElectronLHAdder(const edm::ParameterSet& iConfig) :
        electronTag_(iConfig.getUntrackedParameter<edm::InputTag>("electronTag")),
        lhFile_(iConfig.getUntrackedParameter<std::string>("lhFile")) {
  produces<pat::ElectronCollection>();  
}



void PatElectronLHAdder::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    edm::Handle<edm::View<reco::Candidate> > electrons;
    iEvent.getByLabel(electronTag_,electrons);
    std::auto_ptr<pat::ElectronCollection> pOut(new pat::ElectronCollection);

    std::cerr << "Warning, this shit ain't implemented yet yo" << std::endl;
    

    for(edm::View<reco::Candidate>::const_iterator ele=electrons->begin(); ele!=electrons->end(); ++ele){    

        pat::Electron clone = *edm::RefToBase<reco::Candidate>(electrons,ele-electrons->begin()).castTo<pat::ElectronRef>();
// 
//         LikelihoodMeasurements measurements;
//         Utils anaUtils;
//         bool inEB=anaUtils.fiducialFlagECAL(fiducialFlagsEle[eleIndex],isEB);
//         measurements.pt = clone->pt();
//         if(clone->isEB() && fabs(clone->eta())<1.0) measurements.subdet = 0;
//         else if (clone->isEB() && fabs(clone->eta())>=1.0) measurements.subdet = 1;
//         else measurements.subdet = 2;
//         measurements.deltaPhi = clone->deltaPhiSuperClusterTrackAtVtx();
//         measurements.deltaEta = clone->deltaEtaSuperClusterTrackAtVtx();
//         measurements.eSeedClusterOverPout = clone->eSeedClusterOverP();
//         measurements.eSuperClusterOverP = clone->eSuperClusterOverP();
//         measurements.hadronicOverEm = clone->hadronicOverEm();
//         measurements.sigmaIEtaIEta = clone->sigmaIetaIeta();
//         measurements.sigmaIPhiIPhi = SigmaiPiP(eleIndex);
//         measurements.fBrem = fbremEle[eleIndex];
//         measurements.nBremClusters = nbremsEle[eleIndex];
//         int gsftrack = gsfTrackIndexEle[eleIndex];
//         TVector3 pIn(pxGsfTrack[gsftrack],pyGsfTrack[gsftrack],pzGsfTrack[gsftrack]);
//         measurements.OneOverEMinusOneOverP = 1./(eSuperClusterOverPEle[eleIndex]*pIn.Mag()) - 1./pIn.Mag();
//         return lh.resultLog(measurements);
// 

//         clone.AddUserFloat(lh);
        pOut->push_back(clone);

    }
    iEvent.put(pOut);
}

PatElectronLHAdder::~PatElectronLHAdder() { } 
void PatElectronLHAdder::beginJob() { }
void PatElectronLHAdder::endJob() { }
DEFINE_FWK_MODULE(PatElectronLHAdder);
