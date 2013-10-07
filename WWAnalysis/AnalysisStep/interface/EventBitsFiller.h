#ifndef WWANALYSIS_ANALYSISSTEP_EVENTFILLER_H
#define WWANALYSIS_ANALYSISSTEP_EVENTFILLER_H

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "WWAnalysis/DataFormats/interface/SkimEvent.h"

#include <stdio.h>
#include <sstream>
#include <vector>
#include <map>
#include <bitset>
#include <algorithm>
#include <limits>
#include <functional>

const size_t MAX = 30;

class EventBitsFiller {
    public:
    

        typedef std::bitset<MAX> bits;

        class AdvancedEventBits {
            public: 
                AdvancedEventBits(const int &b) :
                    passedCuts_(b), run_(0), lumi_(0), evt_(0), instance_(0), weight_(1) { }
                AdvancedEventBits(const bits &b) :
                    passedCuts_(b), run_(0), lumi_(0), evt_(0), instance_(0), weight_(1) { }
                AdvancedEventBits(const int &b, const unsigned int &r, const unsigned int &l, const unsigned int &e, const size_t &i) :
                    passedCuts_(b), run_(r), lumi_(l), evt_(e), instance_(i), weight_(1) { }
                AdvancedEventBits(const unsigned int &r, const unsigned int &l, const unsigned int &e, const size_t &i) :
                    passedCuts_(0), run_(r), lumi_(l), evt_(e), instance_(i), weight_(1) { }

                const unsigned int run() const {return run_;}
                const unsigned int event() const {return evt_;}
                const unsigned int lumi() const {return lumi_;}
                const unsigned int instance() const {return instance_;}

                void addEventVariable(const float &f) { eventVariables_.push_back( f ); }
                void clearEventVariables() { eventVariables_.clear(); }
                float weight() { return weight_; }
                void clearWeight() { weight_ = 1; }
                void addWeight(const float &w) { weight_ *= w; }
                void turnOn(const size_t &bit) { passedCuts_[bit] = true; }

                float getEventVariable(size_t i) const { return ( (i<eventVariables_.size()) ? eventVariables_[i]:-9999. ); }
                const bits & getBits() const {return passedCuts_;}
                bool hasSameBits(const AdvancedEventBits& rhs) const;

                friend std::ostream& operator<<(std::ostream& out, const AdvancedEventBits &evt) {
                    out << std::setw(MAX+1) << evt.getBits() << std::setw(10) << evt.run() << std::setw(10) << evt.lumi() << std::setw(10) << evt.event() << std::setw(5) << evt.instance();
                    return out;
                }

                AdvancedEventBits operator&(const AdvancedEventBits &rhs) const;
                AdvancedEventBits operator|(const AdvancedEventBits &rhs) const;
                const size_t lowestZero() const ;

                bool operator==(const AdvancedEventBits& rhs) const;
                bool operator>(const AdvancedEventBits& rhs) const;
                bool operator<(const AdvancedEventBits& rhs) const;
                bool operator>=(const AdvancedEventBits& rhs) const;
                bool operator<=(const AdvancedEventBits& rhs) const;


            private:
                bits passedCuts_;
                unsigned int  run_;
                unsigned int  lumi_;
                unsigned int  evt_;
                size_t        instance_;
                float         weight_;
                std::vector<float> eventVariables_;
        };

        struct findHypo : std::unary_function<AdvancedEventBits, bool> {
            unsigned int  run_;
            unsigned int  lumi_;
            unsigned int  evt_;
            size_t        instance_;
            findHypo(AdvancedEventBits const & r) : run_(r.run()), lumi_(r.lumi()), evt_(r.event()), instance_(r.instance())  {}
            bool operator()(AdvancedEventBits const& rhs) const {
                return (run_==rhs.run() && lumi_==rhs.lumi() && evt_==rhs.event() && instance_==rhs.instance());
            }
        };

        struct findEvent : std::unary_function<AdvancedEventBits, bool> {
            unsigned int  run_;
            unsigned int  lumi_;
            unsigned int  evt_;
            findEvent(AdvancedEventBits const & r) : run_(r.run()), lumi_(r.lumi()), evt_(r.event()) {}
            bool operator()(AdvancedEventBits const& rhs) const {
                return (run_==rhs.run() && lumi_==rhs.lumi() && evt_==rhs.event() );
            }
        };


        typedef std::vector<AdvancedEventBits> HypoList;
        typedef std::map<std::string,HypoList> HypoSummary;

        EventBitsFiller(const edm::ParameterSet &p, const TFileDirectory &, const std::string &);
        ~EventBitsFiller() { }

        const bool isEventSummaryPopulated() const {return eventSummaryPopulated_;}
        const std::vector<int> getEventYieldVector(const std::string& str);
        void setTotalNumberOfCuts(const std::string &s, const size_t &n) { numCuts_[s] = n<MAX?n:MAX; }

        void operator()(edm::EventBase const * evt,const std::string &str, const size_t &inst, const size_t &cut, const reco::SkimEvent &, const std::vector<float>&);
        void printHypoSummary();
        void printEventSummary();
        void printHypoList(const std::string &str,const size_t &cut);
        void printEventList(const std::string &str,const size_t &cut);
        void printFuckingEverything();
        void setCutLabels(const std::string&, const std::vector<std::string> &l);
        void fillHypoVariables(AdvancedEventBits &hypo, const reco::SkimEvent &se);
        void setWeightVariables(AdvancedEventBits &hypo, const std::vector<float> &se);
        void writeAllYieldHists();
        void writeAllNMinus1Plots();
        void writeAllByCutPlots();
        void writeYieldHist(const std::string &str="");
        void writeNMinus1Plots(const std::string &str);
        void writeByCutPlots(const std::string &str);

    private:

        void setEventSummaryPopulated(const bool &b=true) { eventSummaryPopulated_ = b; }
        void populateEventSummary();

        HypoSummary hypoSummary_;
        HypoSummary eventSummary_;
        bool eventSummaryPopulated_;
        std::vector<AdvancedEventBits> cutMasks_;
        std::map<std::string,size_t> numCuts_;
        std::map<std::string,std::vector<std::string> > cutLabels_;
        edm::ParameterSet histParams_;
        std::vector<StringObjectFunction<reco::SkimEvent> > stringFunctions_;
        TFileDirectory fd_;
        std::string sample_;

};


#endif
