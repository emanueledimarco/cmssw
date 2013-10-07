#include <memory>
#include <vector>
#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "WWAnalysis/AnalysisStep/interface/TriggerBitChecker.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include <map>
#include <string>
#include <iostream>

class GioTriggerEmulator : public edm::EDFilter {
    public:
        explicit GioTriggerEmulator(const edm::ParameterSet&);
        ~GioTriggerEmulator();

    private:
        virtual void beginJob() ;
        virtual bool filter(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        edm::InputTag trigger_;
        edm::InputTag muons_, electrons_;
        std::string doubleMu_, doubleEl_, tripleEl_, mueg_;
        int runForMC_;

        std::map<std::string, int> indexCache_; unsigned int lastRun_;
        bool checkPath(const std::string &path, const edm::TriggerResults &result, const edm::Event&);
        bool HLT_MuX_MuY_emulated(const edm::View<pat::Muon> &muons, const std::string &l3filter, double ptMax, double ptMin);
        
};

GioTriggerEmulator::GioTriggerEmulator(const edm::ParameterSet& iConfig) :
    trigger_(iConfig.getParameter<edm::InputTag>("trigger")),
    muons_(iConfig.getParameter<edm::InputTag>("muons")),
    doubleMu_(iConfig.getParameter<std::string>("doubleMu")),
    doubleEl_(iConfig.getParameter<std::string>("doubleEl")),
    tripleEl_(iConfig.getParameter<std::string>("tripleEl")),
    mueg_(iConfig.getParameter<std::string>("mueg")),
    runForMC_(iConfig.getParameter<uint32_t>("runForMC")),
    lastRun_(0)
{ 
} 


bool GioTriggerEmulator::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    edm::Handle<edm::TriggerResults> triggerResults;
    iEvent.getByLabel(trigger_,triggerResults);
    int run = iEvent.isRealData() ? iEvent.id().run() : runForMC_;
    if (doubleMu_ != "none") {
        if (run <= 165208) { 
            if (checkPath("HLT_DoubleMu7_v", *triggerResults, iEvent)) {
                if (doubleMu_ != "any") {
                    edm::Handle<edm::View<pat::Muon> > muons;
                    iEvent.getByLabel(muons_, muons);
                    if (doubleMu_ == "Mu17_Mu8") {
                        if (HLT_MuX_MuY_emulated(*muons, "hltDiMuonL3PreFiltered7", 17, 8)) return true;
                    } else if (doubleMu_ == "Mu13_Mu8") {
                        if (HLT_MuX_MuY_emulated(*muons, "hltDiMuonL3PreFiltered7", 13, 8)) return true;
                    } else return true;
                } else return true;
            } 
        } else if (run <= 178419) {
            if (checkPath("HLT_Mu13_Mu8_v", *triggerResults, iEvent)) {
                if (doubleMu_ != "any") {
                    edm::Handle<edm::View<pat::Muon> > muons;
                    iEvent.getByLabel(muons_, muons);
                    if (doubleMu_ == "Mu17_Mu8") {
                        if (HLT_MuX_MuY_emulated(*muons, "hltDiMuonL3p5PreFiltered8", 17, 8) || 
                            HLT_MuX_MuY_emulated(*muons, "hltDiMuonL3PreFiltered8",   17, 8)) return true;
                    } else return true;
                } else return true;
            } 
        } else {
            if (checkPath("HLT_Mu17_Mu8_v", *triggerResults, iEvent)) return true;
        }
    } else {
        if (checkPath("HLT_Mu17_Mu8_v", *triggerResults, iEvent)) return true;
    }
    if (doubleEl_ == "Ele17_Ele8_Tight") {
        if (checkPath("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v", *triggerResults, iEvent) ||
            checkPath("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v", *triggerResults, iEvent)) {
            return true;
        }
    } 
    else if (doubleEl_ == "any") {
        if (run <= 170901) {
            if (checkPath("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v", *triggerResults, iEvent)) return true;
        } else {
            if (checkPath("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v", *triggerResults, iEvent) ||
                checkPath("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v", *triggerResults, iEvent)) {
                return true;
            }
        }
    }
    if (tripleEl_ != "none") {
        if (checkPath(tripleEl_, *triggerResults, iEvent)) return true;
    }
    if (mueg_ == "HLT_Mu17_Ele8") {
        if (run <= 175972) { 
            if (checkPath("HLT_Mu17_Ele8_CaloIdL_v", *triggerResults, iEvent)) return true;
        }
        else {
            if (checkPath("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v", *triggerResults, iEvent)) return true;
        } 
    }
    else if (mueg_ == "HLT_Mu8_Ele17") {
        if (run <= 167913) { 
            if (checkPath("HLT_Mu8_Ele17_CaloIdL_v", *triggerResults, iEvent)) return true;
        }
        else {
            if (checkPath("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v", *triggerResults, iEvent)) return true;
        } 
    }
    else if (mueg_ == "HLT_Mu17_8_Ele8_17" || mueg_ == "HLT_Mu8_17_Ele17_8") {
        if (run <= 175972 && checkPath("HLT_Mu17_Ele8_CaloIdL_v", *triggerResults, iEvent)) return true;
        else if (run > 175972 && checkPath("HLT_Mu17_Ele8_CaloIdL_v", *triggerResults, iEvent)) return true; 
        else if (run <= 167913 && checkPath("HLT_Mu8_Ele17_CaloIdL_v", *triggerResults, iEvent)) return true;
        else if (run > 167913 && checkPath("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v", *triggerResults, iEvent)) return true;        
    }    

    return false;
}

bool GioTriggerEmulator::checkPath(const std::string &path, const edm::TriggerResults &result, const edm::Event&event) {
    if (event.id().run() != lastRun_) indexCache_.clear();
    int &idx = indexCache_[path];
    if (idx == 0) {
        const edm::TriggerNames &names = event.triggerNames(result);
        for (size_t i = 0, n = names.size(); i < n; ++i) {
            const std::string &thispath = names.triggerName(i);
            if (thispath.find(path) == 0) { idx = i+1; break;}
        }
        if (idx == 0) idx = -1; // if missing once, assume it stays missing for this run
        lastRun_ = event.id().run();
        //std::cout << " path " << path << " is at index " << idx << " for run " << lastRun_ << std::endl;
    }
    return idx > 0 && result.accept(idx-1);
}

bool GioTriggerEmulator::HLT_MuX_MuY_emulated(const edm::View<pat::Muon> &muons, const std::string &l3filter, double ptMax, double ptMin) {
    int nHigh = 0, nLow = 0;
    for (edm::View<pat::Muon>::const_iterator it = muons.begin(), ed = muons.end(); it != ed; ++it) {
        const pat::TriggerObjectStandAlone *match = it->triggerObjectMatchByFilter(l3filter);
        if (match) {
            if (match->pt() > ptMax) nHigh++;
            if (match->pt() > ptMin) nLow++;
        }
        if (nHigh >= 1 && nLow >= 2) return true;
    }
    return false;
}

GioTriggerEmulator::~GioTriggerEmulator() { }
void GioTriggerEmulator::beginJob() { }
void GioTriggerEmulator::endJob() { }
DEFINE_FWK_MODULE(GioTriggerEmulator);


