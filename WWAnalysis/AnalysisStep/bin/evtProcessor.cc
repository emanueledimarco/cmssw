#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "PhysicsTools/SelectorUtils/interface/RunLumiSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"

//#if !defined(__CINT__) || !defined(__MAKECINT__)
//Headers for the data items
//#include <WWAnalysis/DataFormats/interface/SkimEvent.h>
//#endif

#include "WWAnalysis/DataFormats/interface/SkimEvent.h"

#define DEBUG 0

#include <TSystem.h>
#include "TDirectory.h"
#include "TPad.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2.h"
#include "TF1.h"
#include "TLine.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "Math/VectorUtil.h"

#include <stdio.h>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <limits>
#include <functional>
#include <bitset>
using namespace std;
const size_t MAX = 30;
typedef bitset<MAX> bits;

class EventFiller {
    public:

        class MyHypoStruct {
            public: 
                MyHypoStruct(const int &b) :
                    passedCuts_(b), run_(0), lumi_(0), evt_(0), instance_(0) { if(DEBUG) cerr << " MyHypoStruct::MyHypoStruct(int)" << endl; }
                MyHypoStruct(const bits &b) :
                    passedCuts_(b), run_(0), lumi_(0), evt_(0), instance_(0) { if(DEBUG) cerr << " MyHypoStruct::MyHypoStruct(bits)" << endl; }
                MyHypoStruct(const int &b, const unsigned int &r, const unsigned int &l, const unsigned int &e, const size_t &i) :
                    passedCuts_(b), run_(r), lumi_(l), evt_(e), instance_(i) { if(DEBUG) cerr << " MyHypoStruct::MyHypoStruct(b,r,l,e,i)" << endl; }
                MyHypoStruct(const unsigned int &r, const unsigned int &l, const unsigned int &e, const size_t &i) :
                    passedCuts_(0), run_(r), lumi_(l), evt_(e), instance_(i) { if(DEBUG) cerr << " MyHypoStruct::MyHypoStruct(r,l,e,i)" << endl; }

                const unsigned int run() const {return run_;}
                const unsigned int event() const {return evt_;}
                const unsigned int lumi() const {return lumi_;}
                const unsigned int instance() const {return instance_;}

                void addEventVariable(const float &f) { eventVariables_.push_back( f ); }
                void clearEventVariables() { eventVariables_.clear(); }
                float getEventVariable(size_t i) { return ( (i<eventVariables_.size()) ? eventVariables_[i]:-9999. ); }
                void turnOn(const size_t &bit) { passedCuts_[bit] = true; }

                friend ostream& operator<<(ostream& out, const MyHypoStruct &evt) {
                    out << setw(MAX+1) << evt.getBits() << setw(10) << evt.run() << setw(10) << evt.lumi() << setw(10) << evt.event() << setw(5) << evt.instance();
                    return out;
                }
                MyHypoStruct operator&(const MyHypoStruct &rhs) const {
                    return (getBits()&rhs.getBits());
                }
                MyHypoStruct operator|(const MyHypoStruct &rhs) const {
                    return (getBits()|rhs.getBits());
                }
                const size_t lowestZero() const { 
                    for(size_t i = 0;i<MAX;++i) if(passedCuts_[i]==false) return i;
                    return MAX;
                }
                bool hasSameBits(const MyHypoStruct& rhs) const {
                    return ( getBits() == rhs.getBits() );
                }
                bool operator==(const MyHypoStruct& rhs) const {
                    return ( run()==rhs.run() && lumi()==rhs.lumi() && event()==rhs.event());
                }
                bool operator>(const MyHypoStruct& rhs) const {
                    return ( lowestZero() > rhs.lowestZero() );
                }
                bool operator<(const MyHypoStruct& rhs) const {
                    return ( lowestZero() < rhs.lowestZero() );
                }
                bool operator>=(const MyHypoStruct& rhs) const {
                    return ( lowestZero() >= rhs.lowestZero() );
                }
                bool operator<=(const MyHypoStruct& rhs) const {
                    return ( lowestZero() <= rhs.lowestZero() );
                }
                const bits & getBits() const {return passedCuts_;}
            private:
                bits passedCuts_;
                unsigned int  run_;
                unsigned int  lumi_;
                unsigned int  evt_;
                size_t        instance_;
                vector<float> eventVariables_;
        };

        struct findHypo : std::unary_function<MyHypoStruct, bool> {
            unsigned int  run_;
            unsigned int  lumi_;
            unsigned int  evt_;
            size_t        instance_;
            findHypo(MyHypoStruct const & r) : run_(r.run()), lumi_(r.lumi()), evt_(r.event()), instance_(r.instance())  {}
            bool operator()(MyHypoStruct const& rhs) const {
                return (run_==rhs.run() && lumi_==rhs.lumi() && evt_==rhs.event() && instance_==rhs.instance());
            }
        };

        struct findEvent : std::unary_function<MyHypoStruct, bool> {
            unsigned int  run_;
            unsigned int  lumi_;
            unsigned int  evt_;
            findEvent(MyHypoStruct const & r) : run_(r.run()), lumi_(r.lumi()), evt_(r.event()) {}
            bool operator()(MyHypoStruct const& rhs) const {
                return (run_==rhs.run() && lumi_==rhs.lumi() && evt_==rhs.event() );
            }
        };


        typedef vector<MyHypoStruct> HypoList;
        typedef map<string,HypoList> HypoSummary;

        EventFiller() : 
                eventSummaryPopulated_(false), 
                numCuts_(MAX), 
                cutLabels_(MAX,"") { 
            cutMasks_.push_back(1);
            for(size_t i=1;i<MAX;++i) {
                cutMasks_.push_back( cutMasks_[i-1]|(1<<i) );
                if(DEBUG) cerr << "cutMasks_[" << i << "] == " << cutMasks_[i] << endl;
            }
        }
        EventFiller(const edm::ParameterSet &p) : 
                eventSummaryPopulated_(false), 
                numCuts_(MAX), 
                cutLabels_(MAX,"") ,
                histParams_(p) { 

            cutMasks_.push_back(1);
            for(size_t i=1;i<MAX;++i) {
                cutMasks_.push_back( cutMasks_[i-1]|(1<<i) );
                if(DEBUG) cerr << "cutMasks_[" << i << "] == " << cutMasks_[i] << endl;
            }

            vector<string> plots = histParams_.getParameterNamesForType<edm::ParameterSet>();
            string var;
            for(size_t i=0;i<plots.size();i++) {
                var = histParams_.getParameter<edm::ParameterSet>(plots[i]).getParameter<string>("variable");
                stringFunctions_.push_back( StringObjectFunction<reco::SkimEvent>( var ) );
            }
        }
        ~EventFiller() { }

        void operator()(edm::EventBase const * evt,const string &str, const size_t &inst, const size_t &cut, const reco::SkimEvent &);
        void printHypoSummary();
        void printEventSummary();
        const bool isEventSummaryPopulated() const {return eventSummaryPopulated_;}
        void printHypoList(const string &str,const size_t &cut);
        void printEventList(const string &str,const size_t &cut);
        void printFuckingEverything();
        void writeYieldHistToFile(const string &fn, const string &str="");
        void writeYieldHistsToFile(const string &fn);
        void setTotalNumberOfCuts(const size_t &n) { numCuts_ = n<MAX?n:MAX; }
        void setCutLabels(const vector<string> &l);
        const vector<int> getEventYieldVector(const string& str);
        void fillHypoVariables(MyHypoStruct &hypo, const reco::SkimEvent &se);
        void writeNMinus1PlotsToFile(const string &fn, const string &str);
        void writeAllNMinus1PlotsToFile(const string &fn);

    private:

        void setEventSummaryPopulated(const bool &b=true) { eventSummaryPopulated_ = b; }
        void populateEventSummary();

        HypoSummary hypoSummary_;
        HypoSummary eventSummary_;
        bool eventSummaryPopulated_;
        vector<MyHypoStruct> cutMasks_;
        size_t numCuts_;
        vector<string> cutLabels_;
        edm::ParameterSet histParams_;
        vector<StringObjectFunction<reco::SkimEvent> > stringFunctions_;

};

void EventFiller::setCutLabels(const vector<string> &l) {
    
    for(size_t i=0;i<l.size()&&i<cutLabels_.size();++i) {
        cutLabels_[i] = l[i];
    }

}

void EventFiller::fillHypoVariables(MyHypoStruct &hypo, const reco::SkimEvent &se) {
    
    hypo.clearEventVariables();
    for(size_t i=0;i<stringFunctions_.size();i++) {
        hypo.addEventVariable( stringFunctions_[i](se) );
    }

}

void EventFiller::operator()(edm::EventBase const * evt,const string &str, const size_t &inst, const size_t &cut, const reco::SkimEvent &se ) {

    setEventSummaryPopulated(false);

    MyHypoStruct temp(
        evt->eventAuxiliary().run(),
        evt->eventAuxiliary().luminosityBlock(),
        evt->eventAuxiliary().event(),
        inst
    );

    HypoList::iterator thisIt = find_if( hypoSummary_[str].begin(), hypoSummary_[str].end(), findHypo(temp) );
    if(thisIt == hypoSummary_[str].end()) {
        hypoSummary_[str].push_back(temp);
        thisIt = hypoSummary_[str].end()-1;
        fillHypoVariables(*thisIt,se);
    }
    thisIt->turnOn(cut);
}

void EventFiller::printHypoSummary() {

    for(HypoSummary::iterator sumIt=hypoSummary_.begin();sumIt!=hypoSummary_.end();++sumIt) {
        cout << "Hypotheses Breakdown for " << sumIt->first << endl;
        vector<int> perCut(numCuts_,0);
        for(HypoList::iterator evtIt=sumIt->second.begin(); evtIt!=sumIt->second.end();++evtIt) {
            for(size_t i=0;i<numCuts_;++i) {
                if( ((*evtIt)&cutMasks_[i]).hasSameBits(cutMasks_[i]) ) perCut[i]++;
            }
        }
        for(size_t i=0;i<numCuts_;++i) {
            cout << setw(3) << i << setw(7) << perCut[i] << endl;
        }
    }
    
}

void EventFiller::populateEventSummary() {

    eventSummary_.clear();
    for(HypoSummary::iterator sumIt=hypoSummary_.begin();sumIt!=hypoSummary_.end();++sumIt) {
        for(HypoList::iterator evtIt=sumIt->second.begin(); evtIt!=sumIt->second.end();++evtIt) {
            HypoList::iterator curPos = find_if(eventSummary_[sumIt->first].begin(),eventSummary_[sumIt->first].end(),findEvent(*evtIt));
            if( curPos == eventSummary_[sumIt->first].end() ) {
                eventSummary_[sumIt->first].push_back(*evtIt);
            } else {
                if( (*evtIt) > (*curPos) ) (*curPos) = (*evtIt);
            }
        }
    }
    setEventSummaryPopulated();
}

const vector<int> EventFiller::getEventYieldVector(const string& str) {
    if( !isEventSummaryPopulated() ) populateEventSummary();
    vector<int> perCut(numCuts_,0);
    for(HypoList::iterator evtIt=eventSummary_[str].begin(); evtIt!=eventSummary_[str].end();++evtIt) {
        for(size_t i=0;i<numCuts_;++i) {
            if( ((*evtIt)&cutMasks_[i]).hasSameBits(cutMasks_[i]) ) perCut[i]++;
        }
    }
    return perCut;
}

void EventFiller::printEventSummary() {

    if( !isEventSummaryPopulated() ) populateEventSummary();

    for(HypoSummary::iterator sumIt=eventSummary_.begin();sumIt!=eventSummary_.end();++sumIt) {
        cout << "Event Breakdown for " << sumIt->first << endl;
        const vector<int> &perCut = getEventYieldVector(sumIt->first);
        for(size_t i=0;i<perCut.size();++i) {
            cout << setw(3) << i << setw(7) << perCut[i] << endl;
        }
    }
    
}

void EventFiller::printHypoList(const string &str, const size_t &cut) {

    cout << "Hypo List for cut " << cut << " and hypothesis " << str << endl;
    for(HypoList::iterator evtIt=hypoSummary_[str].begin(); evtIt!=hypoSummary_[str].end();++evtIt) {
        if( ((*evtIt)&cutMasks_[cut]).hasSameBits(cutMasks_[cut]) ) 
            cout << "Hypo: " << setw(10) << str << ", Cut: " << setw(4) << cut << ", Event: " << *evtIt << endl;
    }

}

void EventFiller::printEventList(const string &str, const size_t &cut) {

    if( !isEventSummaryPopulated() ) populateEventSummary();

    cout << "Event List for cut " << cut << " and hypothesis " << str << endl;
    for(HypoList::iterator evtIt=eventSummary_[str].begin(); evtIt!=eventSummary_[str].end();++evtIt) {
        if( ((*evtIt)&cutMasks_[cut]).hasSameBits(cutMasks_[cut]) ) 
            cout << "Hypo: " << setw(10) << str << ", Cut: " << setw(4) << cut << ", Event: " << *evtIt << endl;
    }

}

void EventFiller::printFuckingEverything() {

    if( !isEventSummaryPopulated() ) populateEventSummary();
    printEventSummary();
    for(HypoSummary::iterator sumIt=eventSummary_.begin();sumIt!=eventSummary_.end();++sumIt) {
        for(size_t i=0;i<numCuts_;++i) {
            printEventList(sumIt->first,i);
        }
    }

}

void EventFiller::writeYieldHistToFile(const string &fn, const string &str) {

    if( !isEventSummaryPopulated() ) populateEventSummary();
    TFile *f = TFile::Open(fn.c_str(),"UPDATE");
    TH1F *h;

    if( !f->Get("yields") ) f->mkdir( "yields" );
    f->cd("yields");

    h = new TH1F(str.c_str(),str.c_str(),numCuts_,0,numCuts_);
    h->Sumw2();

    const vector<int> &perCut = getEventYieldVector(str);
    double total = 0;
    for(size_t i=0;i<perCut.size();++i) {
        h->SetBinContent(i+1,perCut[i]);
        h->SetBinError(i+1,TMath::Sqrt(perCut[i]));
        h->GetXaxis()->SetBinLabel(i+1,cutLabels_[i].c_str());
        h->LabelsOption("v");
        total+=perCut[i];
    }
    h->SetEntries(total);
    
    h->Write(0,TObject::kOverwrite);
    f->Close();

}

void EventFiller::writeNMinus1PlotsToFile(const string &fn, const string &str) {

    if( !isEventSummaryPopulated() ) populateEventSummary();
    TFile *f = TFile::Open(fn.c_str(),"UPDATE");

    if( !f->Get("nminus1") ) f->mkdir( "nminus1" );
    if( !f->Get(("nminus1/"+str).c_str()) ) {
        TDirectory *g = (TDirectory*)f->Get("nminus1");
        g->mkdir( str.c_str() );
    }
    f->cd( ("nminus1/"+str).c_str() );

    vector<string> plots = histParams_.getParameterNamesForType<edm::ParameterSet>();

    vector<vector<TH1F*> > allHists(plots.size());
    edm::ParameterSet plot;

    ostringstream os;
    for(size_t i=0;i<plots.size();++i) {
        for(size_t cut=0;cut<numCuts_;++cut) {

            plot = histParams_.getParameter<edm::ParameterSet>(plots[i]);
            os << plot.getParameter<string>("variable") << setw(2) << setfill('0') << cut;
            allHists[i].push_back( new TH1F(
                os.str().c_str(),
                plot.getParameter<string>("variable").c_str(),
                plot.getParameter<uint>("nbins"),
                plot.getParameter<double>("low"),
                plot.getParameter<double>("high")
            ));
            allHists[i].back()->GetXaxis()->SetTitle( plot.getParameter<string>("xtitle").c_str() );
            allHists[i].back()->Sumw2();
            os.clear(); os.str("");
        }
    }

    for(HypoList::iterator evtIt=eventSummary_[str].begin(); evtIt!=eventSummary_[str].end();++evtIt) {
        for(size_t i=0;i<plots.size();++i) {
            for(size_t cut=0;cut<numCuts_;++cut) {
                if( ((*evtIt)|(MyHypoStruct( bits((1<<cut)) )) ).lowestZero() == numCuts_) {
//                     cout << "I suppose that I never get here huh?" << ((*evtIt)|(MyHypoStruct( bits((1<<cut)) )) ).lowestZero() << endl;
                    allHists[i][cut]->Fill( evtIt->getEventVariable(i) );
                }
            }
        }
    }

    for(size_t i=0;i<plots.size();++i) {
        for(size_t cut=0;cut<numCuts_;++cut) {
            allHists[i][cut]->Write(0,TObject::kOverwrite);
        }
    }

    f->Close();

}

void EventFiller::writeAllNMinus1PlotsToFile(const string &fn) {

    // once for each hypothesis
    if( !isEventSummaryPopulated() ) populateEventSummary();
    for(HypoSummary::iterator sumIt=eventSummary_.begin();sumIt!=eventSummary_.end();++sumIt)  {
        writeNMinus1PlotsToFile(fn,sumIt->first);
    }

}

void EventFiller::writeYieldHistsToFile(const string &fn) {

    // once for each hypothesis
    if( !isEventSummaryPopulated() ) populateEventSummary();
    for(HypoSummary::iterator sumIt=eventSummary_.begin();sumIt!=eventSummary_.end();++sumIt) 
        writeYieldHistToFile(fn,sumIt->first);

}


//cut bits
enum cutBits {
    SKIMMED = 0,
    FIDUCIAL,
    ID,
    ISO,
    CONVERSION,
    IP,
    LEPTONVETO,
    MET,
    MLL,
    MZ,
    PROJMET,
    JETVETO,
    SOFTMU,
    BJETS,
    MLL2,
    PTMAX,
    PTMIN,
    DPHI
};


int main(int argc,char* argv[]) {

    gSystem->Load( "libFWCoreFWLite" );
    AutoLibraryLoader::enable();

    if ( argc < 2 ) {
        cout << "Usage : " << argv[0] << " [parameters.py]" << endl;
        return 0;
    }

    // Get the python configuration
    PythonProcessDesc builder(argv[1]);
    edm::ParameterSet const& allPars = 
        builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("yieldParams");

    // output name, to be handled here or by the Event Describer thingy
    string outputFileName(  allPars.getParameter<string>("fileOutName") );

    // if allPars has the right parameter, it will filter events by lumi
    RunLumiSelector myLumiSel(allPars);

    vector<string> cutLabels;
    cutLabels.push_back("10/10");
    cutLabels.push_back("20/10");
    cutLabels.push_back("Lepton ID");
    cutLabels.push_back("Lepton ISO");
    cutLabels.push_back("Conversion");
    cutLabels.push_back("d0/dZ");
    cutLabels.push_back("Lepton Veto");
    cutLabels.push_back("#slash{E}_{T}");
    cutLabels.push_back("m_{ll} > 12");
    cutLabels.push_back("Z-veto");
    cutLabels.push_back("p#slash{E}_{T}");
    cutLabels.push_back("Jet Veto");
    cutLabels.push_back("Soft #mu Veto");
    cutLabels.push_back("b-tag Veto");
    cutLabels.push_back("m_{ll} < 50");
    cutLabels.push_back("p_{T}^{MAX} > 30");
    cutLabels.push_back("p_{T}^{MIN} > 25");
    cutLabels.push_back("#Delta#phi < 60");

    // selection variables
    double etaMu = 0;
    double etaEl = 0;
    double ptMin = 0;
    double ptMax = 0;
    double d0 = 0;
    double dZ = 0;
    double isoEl = 0;
    double isoMu = 0;
    double met = 0;
    double mll = 0;
    double mZ = 0;
    double pMet = 0;
    double jetPt = 0;
    double jetEta = 0;
    bool useJEC = false;
    int nCentralJet = 0;
    int nSoftMu = 0;
    int nExtraLep = 0;
    double bValue = 0;
    int nBtagJets = 0;
    double mllMaxFinal = 9999999;
    double ptMinFinal  = 0;
    double ptMaxFinal  = 0;
    double deltaPhiLL  = 0; 

    EventFiller eventFiller(allPars.getParameter<edm::ParameterSet>("histParams"));

    // vector of strings, each name is a SkimEvent stored in the file
    edm::ParameterSet selectionParams = allPars.getParameter<edm::ParameterSet>("selectionParams");
    vector<string> hypoTypes = selectionParams.getParameterNamesForType<edm::ParameterSet>();
    vector<fwlite::Handle<vector<reco::SkimEvent> >* > hypoHandles;
    for(vector<string>::const_iterator itHypo=hypoTypes.begin();itHypo!=hypoTypes.end();++itHypo) {
        hypoHandles.push_back( new fwlite::Handle<vector<reco::SkimEvent> >() );
    }

    // vector of strings, each describes an input file
    edm::ParameterSet sampleInputParams = allPars.getParameter<edm::ParameterSet>("inputParams");
    vector<string> sampleInputs = sampleInputParams.getParameterNamesForType<edm::ParameterSet>();


    //loop on datasets
    for(vector<string>::const_iterator itSample=sampleInputs.begin();itSample!=sampleInputs.end();++itSample) {

        // parameter set describing the input file
        edm::ParameterSet input = sampleInputParams.getParameter<edm::ParameterSet>(*itSample);
        fwlite::ChainEvent ev(input.getParameter<vector<string> >("files"));

        // variable just for letting user know wtf is going on
        int evtCount = 0;
        int maxCount = input.getParameter<int>("maxEvents");
        if( maxCount == -1 ) maxCount = numeric_limits<int>::max();

        //loop on events                                     // and fail if not part of lumi list
        if(ev.size()) for( ev.toBegin(); ! ev.atEnd(); ++ev) if(myLumiSel(ev)) {

            //stop if the user is tired o this stuff
            if(evtCount > maxCount) break;
            // let the user know, ya know?
            if( evtCount++%1000 == 0 ) cerr << "Processing " << *itSample << ": " << setw(10) << evtCount << endl;

            //loop on hypothesis
            for(size_t hypoI=0;hypoI<hypoHandles.size();++hypoI) {

                // get the skim event from the file
                hypoHandles[hypoI]->getByLabel(ev,hypoTypes[hypoI].c_str());

                // update the cut variables
                edm::ParameterSet thisSelection = selectionParams.getParameter<edm::ParameterSet>(hypoTypes[hypoI]);
                etaMu = thisSelection.getParameter<double>("etaMu");
                etaEl = thisSelection.getParameter<double>("etaEl");
                ptMax = thisSelection.getParameter<double>("ptMax");
                ptMin = thisSelection.getParameter<double>("ptMin");
                d0 = thisSelection.getParameter<double>("d0");
                dZ = thisSelection.getParameter<double>("dZ");
                isoEl = thisSelection.getParameter<double>("isoEl");
                isoMu = thisSelection.getParameter<double>("isoMu");
                met = thisSelection.getParameter<double>("met");
                mll = thisSelection.getParameter<double>("mll");
                mZ = thisSelection.getParameter<double>("mZ");
                pMet = thisSelection.getParameter<double>("pMet");
                jetPt = thisSelection.getParameter<double>("jetPt");
                jetEta = thisSelection.getParameter<double>("jetEta");
                useJEC = thisSelection.getParameter<bool>("useJEC");
                nCentralJet = thisSelection.getParameter<int>("nCentralJet");
                nSoftMu = thisSelection.getParameter<int>("nSoftMu");
                nExtraLep = thisSelection.getParameter<int>("nExtraLep");
                bValue = thisSelection.getParameter<double>("bValue");
                nBtagJets = thisSelection.getParameter<int>("nBtagJets");
                mllMaxFinal = thisSelection.getParameter<double>("mllMaxFinal");
                ptMinFinal = thisSelection.getParameter<double>("ptMinFinal");
                ptMaxFinal = thisSelection.getParameter<double>("ptMaxFinal");
                deltaPhiLL = thisSelection.getParameter<double>("deltaPhiLL");

                //loop on no of hypos in event
                for(vector<reco::SkimEvent>::const_iterator mySkimEvent = hypoHandles[hypoI]->ptr()->begin();
                        mySkimEvent != hypoHandles[hypoI]->ptr()->end(); mySkimEvent++){

                       size_t instance = mySkimEvent - hypoHandles[hypoI]->ptr()->begin();
                       eventFiller(ev.event(),hypoTypes[hypoI],instance,SKIMMED,*mySkimEvent);

                    if( mySkimEvent->q(0)*mySkimEvent->q(1) < 0 && 
                       !mySkimEvent->isSTA(0) && !mySkimEvent->isSTA(1) && 
                        mySkimEvent->leptEtaCut(etaMu, etaEl) && 
                        mySkimEvent->ptMin() > ptMin && 
                        mySkimEvent->ptMax() > ptMax ) {

                        eventFiller(ev.event(),hypoTypes[hypoI],instance,FIDUCIAL,*mySkimEvent);
                        eventFiller(ev.event(),hypoTypes[hypoI],instance,ID,*mySkimEvent);
                        eventFiller(ev.event(),hypoTypes[hypoI],instance,ISO,*mySkimEvent);
                    }

//                     if( mySkimEvent->passesIDV1(0) && mySkimEvent->passesIDV1(1) ) {
//                         eventFiller(ev.event(),hypoTypes[hypoI],instance,ID,*mySkimEvent);
//                     }

//                     if( !(abs(mySkimEvent->pdgId(0)) == 11 && mySkimEvent->allIso(0)/mySkimEvent->pt(0) >= isoEl) && 
//                         !(abs(mySkimEvent->pdgId(1)) == 11 && mySkimEvent->allIso(1)/mySkimEvent->pt(1) >= isoEl) && 
//                         !(abs(mySkimEvent->pdgId(0)) == 13 && mySkimEvent->allIso(0)/mySkimEvent->pt(0) >= isoMu) && 
//                         !(abs(mySkimEvent->pdgId(1)) == 13 && mySkimEvent->allIso(1)/mySkimEvent->pt(1) >= isoMu) ) {
// 
//                         eventFiller(ev.event(),hypoTypes[hypoI],instance,ISO,*mySkimEvent);
// 
//                     }

                    if( mySkimEvent->passesConversion(0) && mySkimEvent->passesConversion(1) ) {
                        eventFiller(ev.event(),hypoTypes[hypoI],instance,CONVERSION,*mySkimEvent);
                    }

                    if( fabs(mySkimEvent->d0Reco(0)) < d0 && fabs(mySkimEvent->d0Reco(1)) < d0 &&
                        fabs(mySkimEvent->dZReco(0)) < dZ && fabs(mySkimEvent->dZReco(1)) < dZ ) {
                        eventFiller(ev.event(),hypoTypes[hypoI],instance,IP,*mySkimEvent);
                    }

                    if( mySkimEvent->nExtraLep(10) <= nExtraLep ) {
                        eventFiller(ev.event(),hypoTypes[hypoI],instance,LEPTONVETO,*mySkimEvent);
                    }

                    if( mySkimEvent->tcMet() > met) {
                        eventFiller(ev.event(),hypoTypes[hypoI],instance,MET,*mySkimEvent);
                    }

                    if( mySkimEvent->mll()>mll) {
                        eventFiller(ev.event(),hypoTypes[hypoI],instance,MLL,*mySkimEvent);
                    }

                    if( fabs(mySkimEvent->mll()-91.1876)>mZ ) {
                        eventFiller(ev.event(),hypoTypes[hypoI],instance,MZ,*mySkimEvent);
                    }

                    if( mySkimEvent->projTcMet() > pMet ) {
                        eventFiller(ev.event(),hypoTypes[hypoI],instance,PROJMET,*mySkimEvent);
                    }

                    if( mySkimEvent->nCentralJets( jetPt, jetEta, useJEC) <= nCentralJet) {
                        eventFiller(ev.event(),hypoTypes[hypoI],instance,JETVETO,*mySkimEvent);
                    }

                    if( mySkimEvent->nSoftMu(3.) <= nSoftMu ) {
                        eventFiller(ev.event(),hypoTypes[hypoI],instance,SOFTMU,*mySkimEvent);
                    }

                    if( mySkimEvent->bTaggedJetsUnder( jetPt, bValue) <= nBtagJets ) {
                        eventFiller(ev.event(),hypoTypes[hypoI],instance,BJETS,*mySkimEvent);
                    }

                    if(mySkimEvent->mll() < mllMaxFinal ) {
                        eventFiller(ev.event(),hypoTypes[hypoI],instance,MLL2,*mySkimEvent);
                    }

                    if(mySkimEvent->ptMax() > ptMaxFinal ) {
                        eventFiller(ev.event(),hypoTypes[hypoI],instance,PTMAX,*mySkimEvent);
                    }

                    if(mySkimEvent->ptMin() > ptMinFinal ) {
                        eventFiller(ev.event(),hypoTypes[hypoI],instance,PTMIN,*mySkimEvent);
                    }

                    if(mySkimEvent->dPhill() < deltaPhiLL ) {
                        eventFiller(ev.event(),hypoTypes[hypoI],instance,DPHI,*mySkimEvent);
                    }


                } //end of loop over SkimEvent

            } // end loop over hypothesis

        } // end loop over edm::events 

        if( input.getParameter<bool>("printEvents") ) eventFiller.printFuckingEverything();
        eventFiller.setTotalNumberOfCuts(18);
        eventFiller.setCutLabels(cutLabels);
        eventFiller.writeYieldHistsToFile(outputFileName);
        eventFiller.writeAllNMinus1PlotsToFile(outputFileName);



    } //end loop over input datasets



    return 0;
}


