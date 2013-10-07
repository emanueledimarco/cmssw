#ifndef MYPRODUCERS_MCDumper_MCDumper_H
#define MYPRODUCERS_MCDumper_MCDumper_H

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "WWAnalysis/AnalysisStep/interface/LatinoReWeighting.h"


#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/METReco/interface/GenMET.h"

#include <TTree.h>
#include <vector>
#include <string>


class MCDumper : public edm::EDAnalyzer {
    public:
        explicit MCDumper(const edm::ParameterSet&);
        ~MCDumper();

    private:
        virtual void beginJob() ;
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        edm::InputTag GenJetTag_;
        edm::InputTag GenMetTag_;

        void FindDaughterParticles(const reco::Candidate** pCurrent, std::vector<const reco::Candidate*>* pFinal = 0) const;

        TFileService *fs_;
        TTree *tree_;

        //---- Higgs and decay ----
        float MHiggs_;
        float PtHiggs_;
        float EtaHiggs_;
        float PhiHiggs_;
        float PdgIdHiggs_;

        float V1pt_;
        float V1eta_;
        float V1phi_;
        float V1pdgId_;
        float V2pt_;
        float V2eta_;
        float V2phi_;
        float V2pdgId_;

        float f1V1pt_;
        float f1V1eta_;
        float f1V1phi_;
        float f1V1pdgId_;
        float f2V1pt_;
        float f2V1eta_;
        float f2V1phi_;
        float f2V1pdgId_;

        float f1V2pt_;
        float f1V2eta_;
        float f1V2phi_;
        float f1V2pdgId_;
        float f2V2pt_;
        float f2V2eta_;
        float f2V2phi_;
        float f2V2pdgId_;

        //---- GenJet ----
        float jet1pt_;
        float jet1eta_;
        float jet1phi_;

        float jet2pt_;
        float jet2eta_;
        float jet2phi_;

        float jet3pt_;
        float jet3eta_;
        float jet3phi_;

        float jet4pt_;
        float jet4eta_;
        float jet4phi_;

        float jet5pt_;
        float jet5eta_;
        float jet5phi_;


        float metpt_;
        float metphi_;

        reco::GenParticleRefVector genParticles_;
};

#endif
