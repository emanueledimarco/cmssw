#include <vector>
#include <algorithm>
#include <cmath>

#include "WWAnalysis/DataFormats/interface/SkimEvent2L.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/ShallowClonePtrCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "WWAnalysis/AnalysisStep/interface/CompositeCandMassResolution.h"


class SkimEvent2LProducer : public edm::EDProducer {
    public:
        SkimEvent2LProducer(const edm::ParameterSet &iConfig) ;
        virtual void produce(edm::Event &iEvent, const edm::EventSetup &iSetup) ;
    private:
        edm::InputTag src_;
        edm::InputTag pfMet_;
        edm::InputTag vertices_;
        bool isMC, doGenMatch_;
        edm::InputTag genMatch_;
        bool doMassRes_, doExtendedMassRes_;
        CompositeCandMassResolution massRes_;
};

SkimEvent2LProducer::SkimEvent2LProducer(const edm::ParameterSet &iConfig) :
    src_(iConfig.getParameter<edm::InputTag>("src")),
    pfMet_(iConfig.getParameter<edm::InputTag>("pfMet")),
    vertices_(iConfig.getParameter<edm::InputTag>("vertices")),
    isMC(iConfig.existsAs<bool>("isMC")?iConfig.getParameter<bool>("isMC"):false),
    doGenMatch_(isMC && iConfig.existsAs<edm::InputTag>("genMatch")),
    genMatch_(doGenMatch_ ? iConfig.getParameter<edm::InputTag>("genMatch") : edm::InputTag("FAKE")),
    doMassRes_(iConfig.existsAs<bool>("doMassRes")?iConfig.getParameter<bool>("doMassRes"):false),
    doExtendedMassRes_(iConfig.existsAs<bool>("doExtendedMassRes")?iConfig.getParameter<bool>("doExtendedMassRes"):false)
{
    produces<std::vector<reco::SkimEvent2L> >();
}

void
SkimEvent2LProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
    if (doMassRes_) massRes_.init(iSetup);

    edm::Handle<reco::CandidateView> src; 
    iEvent.getByLabel(src_, src);

    edm::Handle<reco::VertexCollection> vertices; 
    iEvent.getByLabel(vertices_, vertices);
    edm::Handle<reco::PFMETCollection> pfMet; 
    iEvent.getByLabel(pfMet_, pfMet);
    edm::Handle<std::vector<PileupSummaryInfo> > puH;
    if (isMC) iEvent.getByLabel("addPileupInfo", puH);

    edm::Handle<edm::Association<reco::GenParticleCollection> > genMatch; 
    if (doGenMatch_) iEvent.getByLabel(genMatch_, genMatch);

    std::auto_ptr<std::vector<reco::SkimEvent2L> > out(new std::vector<reco::SkimEvent2L>());
    for (reco::CandidateView::const_iterator it = src->begin(), ed= src->end(); it != ed; ++it) {
        const reco::CompositeCandidate *srcz = dynamic_cast<const reco::CompositeCandidate *>(&*it);
        if (srcz == 0) throw cms::Exception("CorruptData") << "Found something that is not a composite candidate" << std::endl;
        const pat::CompositeCandidate *patz = dynamic_cast<const pat::CompositeCandidate *>(srcz);
        if (patz == 0) {
            const reco::Candidate *d0 = srcz->daughter(0);
            if (d0 == 0) throw cms::Exception("CorruptData") << "Null daughter(0)\n";
            out->push_back( reco::SkimEvent2L(*srcz) );
        } else {
            const reco::Candidate *d0 = patz->daughter(0);
            if (d0 == 0) throw cms::Exception("CorruptData") << "Null daughter(0)\n";
            out->push_back( reco::SkimEvent2L(*patz) );
        }
        reco::SkimEvent2L &z = out->back();
        z.setVertices(vertices);
        z.setPFMet(pfMet);
        if (isMC) z.setPileupInfo(puH);
        if (doGenMatch_) z.setGenMatches(*genMatch);
        if (doExtendedMassRes_) {
            std::vector<double> errs;
            z.addUserFloat("massErr", massRes_.getMassResolutionWithComponents(z, errs));
            for (unsigned int i = 0; i < errs.size(); ++i) {
                char buff[20]; sprintf(buff, "massErr[%u]", i);
                z.addUserFloat(buff, errs[i]);
            }
        } else if (doMassRes_) {
            z.addUserFloat("massErr", massRes_.getMassResolution(z));
        } 
        
    }

    iEvent.put(out);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SkimEvent2LProducer);
