#ifndef WWANALYSIS_ANALYSISSTEP_MYSINGLEDEPOSIT_H
#define WWANALYSIS_ANALYSISSTEP_MYSINGLEDEPOSIT_H

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "PhysicsTools/IsolationAlgos/interface/EventDependentAbsVeto.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "DataFormats/Common/interface/OwnVector.h"

#include<string>

//
// class declaration
//

class MySingleDeposit {
    public:
        enum Mode { Sum, SumRelative, Sum2, Sum2Relative, Max, MaxRelative, Count, NearestDR };
        MySingleDeposit(const edm::ParameterSet &p);
        void cleanup() ;
        void open(const edm::Event &iEvent, const edm::EventSetup &iSetup) ;
        double compute(const reco::CandidateBaseRef &cand) ;
        const reco::IsoDepositMap & map() { return *hDeps_; }
        const std::string & getLabel() const {return labelName_;}
        bool isNumber(const std::string&);
        double toNumber(const std::string&);
    private:
        Mode mode_;
        edm::InputTag src_;
        double deltaR_;
        bool   usesFunction_;
        double weight_;
        StringObjectFunction<reco::Candidate> weightExpr_;
        reco::isodeposit::AbsVetos vetos_;
        reco::isodeposit::EventDependentAbsVetos evdepVetos_; // note: these are a subset of the above. Don't delete twice!
        bool   skipDefaultVeto_; 
        edm::Handle<reco::IsoDepositMap> hDeps_; // transient
        std::string labelName_;
};

#endif
