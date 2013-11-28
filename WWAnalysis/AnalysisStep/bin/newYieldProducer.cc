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

#include <WWAnalysis/DataFormats/interface/SkimEvent.h>


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
using namespace std;

#include <bitset>

//Attention, now that we be doin N-1 shizzy, this number has to 
//be changed to reflect the number of cuts we actually be makin
const int maxCuts = 26;
typedef std::bitset<maxCuts> mybit;
const mybit allOn = 0xFFFF;

void fillHists(const reco::SkimEvent&,std::vector<TH1F*> &);
void setupHists(const edm::ParameterSet&,std::vector<TH1F*> &, const string&);
void writeHistsToFile(TFile *f, const std::string &folder,const std::string &postfix,std::vector<TH1F*> &);
void setAxisLabels(const TH1*);
TH1F *getNewHistFromPSet( const edm::ParameterSet &pset, const string& );

struct EvtSummary {
    unsigned int run;
    unsigned int lumi;
    unsigned int evt;
    int cut;
    string hypo;

    friend ostream& operator<<(ostream& out, const EvtSummary &evt) {
        out << setw(9) << evt.run << setw(15) << evt.lumi << setw(20) << evt.evt << setw(10) << evt.hypo << setw(4) << evt.cut;
        return out;
    }

    bool operator==(const EvtSummary& rhs) const {
        return ( run==rhs.run && lumi==rhs.lumi && evt==rhs.evt );
    }

    EvtSummary(unsigned int r, unsigned int l, unsigned int e, int c, string h) :
        run(r), lumi(l), evt(e), cut(c), hypo(h) { }
};


int main(int argc,char* argv[]) {

    vector<mybit> cutMasks;
    cutMasks.push_back(mybit(1));
    for(int i=1;i<maxCuts;++i) cutMasks.push_back( cutMasks[i-1]|(mybit(1<<i)) );

    gROOT->Reset();  
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);

    // load framework libraries
    gSystem->Load( "libFWCoreFWLite" );
    AutoLibraryLoader::enable();

    if ( argc < 2 ) {
        std::cout << "Usage : " << argv[0] << " [parameters.py]" << std::endl;
        return 0;
    }

    // Get the python configuration
    PythonProcessDesc builder(argv[1]);
    edm::ParameterSet const& allPars = 
        builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("yieldParams");

    string outputFileName(  allPars.getParameter<string>("fileOutName") );
//     TFile *outputFile = TFile::Open(outputFileName.c_str(),"UPDATE");
//     if(!outputFile || outputFile->IsZombie() ) {
//         cerr << "Dagnabit, your output file sucks ..." << endl;
//         return 1;
//     }
//     outputFile->mkdir("raw");

//     double scaleToData = allPars.getParameter<double>("integratedLumi");
//     reco::SkimEvent::setupJEC(
//         allPars.getParameter<string>("l2File"),
//         allPars.getParameter<string>("l3File"),
//         allPars.getParameter<string>("resFile")
//     );

    RunLumiSelector myLumiSel(allPars);

    map<string,vector<TH1F*> > histsRaw;
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

    // final cuts
    double mllMaxFinal = 9999999;
    double ptMinFinal  = 0;
    double ptMaxFinal  = 0;
    double pMetFinal   = 0; 
    double metOverPtLL = 0;
    double deltaPhiLL  = 9999999;

    edm::ParameterSet selectionParams = allPars.getParameter<edm::ParameterSet>("selectionParams");
    vector<string> hypoTypes = selectionParams.getParameterNamesForType<edm::ParameterSet>();

    edm::ParameterSet sampleInputParams = allPars.getParameter<edm::ParameterSet>("inputParams");
    vector<string> sampleInputs = sampleInputParams.getParameterNamesForType<edm::ParameterSet>();


//     ostringstream os;
//       //samp     //hypo //cut  //plot
//     map<string,map<string,vector<vector<TH1F*> > > > cutByCutHists;
//     map<string,map<string,vector<vector<TH1F*> > > > nMinus1Hists;
//     TDirectory *curDir;
//     for(vector<string>::const_iterator itSample=sampleInputs.begin();itSample!=sampleInputs.end();++itSample) {
//         curDir = outputFile->mkdir(itSample->c_str());
//         for(vector<string>::const_iterator itHypo=hypoTypes.begin();itHypo!=hypoTypes.end();++itHypo) {
//             curDir = curDir->mkdir(itHypo->c_str());
//             curDir->mkdir("byCut");
//             curDir->mkdir("nMinus1");
//             for(int i=0;i<maxCuts;++i) {
//                 cutByCutHists[*itSample][*itHypo].push_back( vector<TH1F*>() );
//                 nMinus1Hists[*itSample][*itHypo].push_back( vector<TH1F*>() );
//                 os.clear(); os.str("");
//                 os << "_cutbycut_" << setw(2) << setfill('0') << i << "_" << *itSample << "_" <<  *itHypo;
//                 setupHists(allPars.getParameter<edm::ParameterSet>("histParams"),cutByCutHists[*itSample][*itHypo][i],os.str());
//                 os.clear(); os.str("");
//                 os << "_nminus1_" << setw(2) << setfill('0') << i << "_" << *itSample << "_" <<  *itHypo;
//                 setupHists(allPars.getParameter<edm::ParameterSet>("histParams"),nMinus1Hists[*itSample][*itHypo][i],os.str());
//             }
//             curDir = curDir->GetMotherDir();
//         }
//         curDir = curDir->mkdir("all");
//         curDir->mkdir("byCut");
//         curDir->mkdir("nMinus1");
//         for(int i=0;i<maxCuts;++i) {
//             cutByCutHists[*itSample]["all"].push_back( vector<TH1F*>() );
//             nMinus1Hists[*itSample]["all"].push_back( vector<TH1F*>() );
//             os.clear(); os.str("");
//             os << "_cutbycut_" << setw(2) << setfill('0') << i << "_" << *itSample << "_" <<  "all";
//             setupHists(allPars.getParameter<edm::ParameterSet>("histParams"),cutByCutHists[*itSample]["all"][i],os.str());
//             os.clear(); os.str("");
//             os << "_nminus1_" << setw(2) << setfill('0') << i << "_" << *itSample << "_" <<  "all";
//             setupHists(allPars.getParameter<edm::ParameterSet>("histParams"),nMinus1Hists[*itSample]["all"][i],os.str());
//         }
//         curDir = curDir->GetMotherDir();
//     }
// 
    //loop on datasets
    for(vector<string>::const_iterator itSample=sampleInputs.begin();itSample!=sampleInputs.end();++itSample) {

        edm::ParameterSet input = sampleInputParams.getParameter<edm::ParameterSet>(*itSample);
        fwlite::ChainEvent ev(input.getParameter<vector<string> >("files"));

        map<string,map<size_t,int> > counterEvents;
        map<string,mybit > passer;
        vector<EvtSummary> eventList;

        int evtCount = 0;
        //loop on events
        if(ev.size()) for( ev.toBegin(); ! ev.atEnd(); ++ev) if(myLumiSel(ev)) {

            if( evtCount++%1000 == 0 ) cerr << "Processing " << *itSample << ": " << setw(10) << evtCount << endl;


            //loop on hypothesis
            passer["all"] = 0;
            for(vector<string>::const_iterator itHypo=hypoTypes.begin();itHypo!=hypoTypes.end();++itHypo) {

                fwlite::Handle<std::vector<reco::SkimEvent> > skimEventH;
                skimEventH.getByLabel(ev,itHypo->c_str());

                edm::ParameterSet thisSelection = selectionParams.getParameter<edm::ParameterSet>(*itHypo);
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
                // final selection cuts
                mllMaxFinal = thisSelection.getParameter<double>("mllMaxFinal");
                ptMinFinal = thisSelection.getParameter<double>("ptMinFinal");
                ptMaxFinal = thisSelection.getParameter<double>("ptMaxFinal");
                pMetFinal = thisSelection.getParameter<double>("pMetFinal");
                metOverPtLL = thisSelection.getParameter<double>("metOverPtLL");
                deltaPhiLL = thisSelection.getParameter<double>("deltaPhiLL");

                passer[*itHypo] = 0;
                //loop on no of hypos in event
                for(vector<reco::SkimEvent>::const_iterator mySkimEvent = skimEventH.ptr()->begin();
                        mySkimEvent != skimEventH.ptr()->end(); mySkimEvent++){


                    int plotNum=0;
                    passer[*itHypo][plotNum]=true;
                    passer["all"][plotNum]=true;

//                     plotNum++;
//                     if( mySkimEvent->hasGoodVertex() ) {
//                         passer[*itHypo][plotNum]=true;
//                         passer["all"][plotNum]=true;
//                     }
// 
//                     plotNum++;
//                     if( mySkimEvent->q(0)*mySkimEvent->q(1) < 0 ) {
//                         passer[*itHypo][plotNum]=true;
//                         passer["all"][plotNum]=true;
//                     }
// 
//                     plotNum++;
//                     if( !mySkimEvent->isSTA(0) && !mySkimEvent->isSTA(1) ) {
//                         passer[*itHypo][plotNum]=true;
//                         passer["all"][plotNum]=true;
//                     }
// 
//                     plotNum++;
//                     if( mySkimEvent->leptEtaCut(etaMu, etaEl) ) {
//                         passer[*itHypo][plotNum]=true;
//                         passer["all"][plotNum]=true;
//                     }
// 
//                     plotNum++;
//                     if( mySkimEvent->ptMin() > ptMin ) {
//                         passer[*itHypo][plotNum]=true;
//                         passer["all"][plotNum]=true;
//                     }
// 
//                     plotNum++;
//                     if( mySkimEvent->ptMax() > ptMax ) {
//                         passer[*itHypo][plotNum]=true;
//                         passer["all"][plotNum]=true;
//                     }
// 
//                     plotNum++;
//                     if( mySkimEvent->passesIDV1(0) && mySkimEvent->passesIDV1(1) ) {
//                         passer[*itHypo][plotNum]=true;
//                         passer["all"][plotNum]=true;
//                     }
// 
//                     plotNum++;
//                     if( mySkimEvent->q(0)*mySkimEvent->q(1) < 0 ) {
//                         passer[*itHypo][plotNum]=true;
//                         passer["all"][plotNum]=true;
//                     }

                    plotNum++;
                    if( mySkimEvent->q(0)*mySkimEvent->q(1) < 0 && 
                       !mySkimEvent->isSTA(0) && !mySkimEvent->isSTA(1) && 
                        mySkimEvent->leptEtaCut(etaMu, etaEl) && 
                        mySkimEvent->ptMin() > ptMin && 
                        mySkimEvent->ptMax() > ptMax ) {
                        passer[*itHypo][plotNum]=true;
                        passer["all"][plotNum]=true;
                    }

//                     plotNum++;
//                     if( mySkimEvent->passesIDV1(0) && mySkimEvent->passesIDV1(1) ) {
//                         passer[*itHypo][plotNum]=true;
//                         passer["all"][plotNum]=true;
//                     }

                    plotNum++;
                    if( !(abs(mySkimEvent->pdgId(0)) == 11 && mySkimEvent->allIso(0)/mySkimEvent->pt(0) >= isoEl) && 
                        !(abs(mySkimEvent->pdgId(1)) == 11 && mySkimEvent->allIso(1)/mySkimEvent->pt(1) >= isoEl) && 
                        !(abs(mySkimEvent->pdgId(0)) == 13 && mySkimEvent->allIso(0)/mySkimEvent->pt(0) >= isoMu) && 
                        !(abs(mySkimEvent->pdgId(1)) == 13 && mySkimEvent->allIso(1)/mySkimEvent->pt(1) >= isoMu) ) {
                        passer[*itHypo][plotNum]=true;
                        passer["all"][plotNum]=true;
                    }

                    plotNum++;
                    if( mySkimEvent->passesConversion(0) && mySkimEvent->passesConversion(1) ) {
                        passer[*itHypo][plotNum]=true;
                        passer["all"][plotNum]=true;
                    }

                    plotNum++;
                    if( fabs(mySkimEvent->d0Reco(0)) < d0 && fabs(mySkimEvent->d0Reco(1)) < d0 &&
                        fabs(mySkimEvent->dZReco(0)) < dZ && fabs(mySkimEvent->dZReco(1)) < dZ ) {
                        passer[*itHypo][plotNum]=true;
                        passer["all"][plotNum]=true;
                    }

                    plotNum++;
                    if( mySkimEvent->nExtraLep(10) <= nExtraLep ) {
                        passer[*itHypo][plotNum]=true;
                        passer["all"][plotNum]=true;
                    }

                    plotNum++;
                    if( mySkimEvent->tcMet() > met) {
                        passer[*itHypo][plotNum]=true;
                        passer["all"][plotNum]=true;
                    }

                    plotNum++;
                    if( mySkimEvent->mll()>mll) {
                        passer[*itHypo][plotNum]=true;
                        passer["all"][plotNum]=true;
                    }

                    plotNum++;
                    if( fabs(mySkimEvent->mll()-91.1876)>mZ ) {
                        passer[*itHypo][plotNum]=true;
                        passer["all"][plotNum]=true;
                    }

                    plotNum++;
                    if( mySkimEvent->projTcMet() > pMet ) {
                        passer[*itHypo][plotNum]=true;
                        passer["all"][plotNum]=true;
                    }

                    plotNum++;
                    if( mySkimEvent->nCentralJets( jetPt, jetEta, useJEC) <= nCentralJet) {
//                         fillHists(*mySkimEvent,cutByCutHists[*itSample][*itHypo][plotNum]);
//                         fillHists(*mySkimEvent,cutByCutHists[*itSample]["all"][plotNum]);
                        passer[*itHypo][plotNum]=true;
                        passer["all"][plotNum]=true;
                    }

                    plotNum++;
                    if( mySkimEvent->nSoftMu(3.) <= nSoftMu ) {
//                         fillHists(*mySkimEvent,cutByCutHists[*itSample][*itHypo][plotNum]);
//                         fillHists(*mySkimEvent,cutByCutHists[*itSample]["all"][plotNum]);
                        passer[*itHypo][plotNum]=true;
                        passer["all"][plotNum]=true;
                    }

                    plotNum++;
                    if( mySkimEvent->bTaggedJetsUnder( jetPt, bValue) <= nBtagJets ) {
//                         fillHists(*mySkimEvent,cutByCutHists[*itSample][*itHypo][plotNum]);
//                         fillHists(*mySkimEvent,cutByCutHists[*itSample]["all"][plotNum]);
                        passer[*itHypo][plotNum]=true;
                        passer["all"][plotNum]=true;
                    }

                    //final selection cuts
                    plotNum++;
                    if(mySkimEvent->mll() < mllMaxFinal ) {
//                         fillHists(*mySkimEvent,cutByCutHists[*itSample][*itHypo][plotNum]);
//                         fillHists(*mySkimEvent,cutByCutHists[*itSample]["all"][plotNum]);
                        passer[*itHypo][plotNum]=true;
                        passer["all"][plotNum]=true;
                    }

                    plotNum++;
                    if(mySkimEvent->ptMax() > ptMaxFinal ) {
//                         fillHists(*mySkimEvent,cutByCutHists[*itSample][*itHypo][plotNum]);
//                         fillHists(*mySkimEvent,cutByCutHists[*itSample]["all"][plotNum]);
                        passer[*itHypo][plotNum]=true;
                        passer["all"][plotNum]=true;
                    }

                    plotNum++;
                    if(mySkimEvent->ptMin() > ptMinFinal ) {
//                         fillHists(*mySkimEvent,cutByCutHists[*itSample][*itHypo][plotNum]);
//                         fillHists(*mySkimEvent,cutByCutHists[*itSample]["all"][plotNum]);
                        passer[*itHypo][plotNum]=true;
                        passer["all"][plotNum]=true;
                    }

                    plotNum++;
                    if(mySkimEvent->dPhill() < deltaPhiLL ) {
//                         fillHists(*mySkimEvent,cutByCutHists[*itSample][*itHypo][plotNum]);
//                         fillHists(*mySkimEvent,cutByCutHists[*itSample]["all"][plotNum]);
                        passer[*itHypo][plotNum]=true;
                        passer["all"][plotNum]=true;
                    }


                } //end of loop over SkimEvent

                EvtSummary tempEvt(
                    ev.eventAuxiliary().run(),
                    ev.eventAuxiliary().luminosityBlock(),
                    ev.eventAuxiliary().event(), 0, *itHypo
                );
                vector<EvtSummary>::iterator myEvt = find(eventList.begin(),eventList.end(),tempEvt);
                if( myEvt == eventList.end() ) {
                    eventList.push_back(tempEvt);
                    myEvt = eventList.end()-1;
                }

                for(int j=0;j<maxCuts;++j) {
                    if( (passer[*itHypo]&cutMasks[j]) == cutMasks[j] && skimEventH->size()==1) {
                        counterEvents[*itHypo][j]++;
//                         fillHists(skimEventH->at(0),cutByCutHists[*itSample][*itHypo][j]);
//                         fillHists(skimEventH->at(0),cutByCutHists[*itSample]["all"][j]);
                    }
//                     if( (passer[*itHypo]|mybit(1<<j)) == allOn && skimEventH->size()==1) {
//                         fillHists(skimEventH->at(0),nMinus1Hists[*itSample][*itHypo][j]);
//                         fillHists(skimEventH->at(0),nMinus1Hists[*itSample]["all"][j]);
//                     }
                    if(j > myEvt->cut) {
                        myEvt->cut = j;
                        myEvt->hypo = *itHypo;
                    }
                }

//                 for(int j=0;j<maxCuts;++j) {
//                     if( (passer[*itHypo]&cutMasks[j]) == cutMasks[j] ) {
//                         counterEvents[*itHypo][j]++;
//                     }
//                 }
            } // end loop over hypothesis

            for(int j=0;j<maxCuts;++j) {
                if( (passer["all"]&cutMasks[j]) == cutMasks[j]) {
                    counterEvents["all"][j]++;
                } 
            }

        } // end loop over edm::events 

        double scaleTo1pb1 = input.getParameter<double>("scale");

//         outputFile->cd("raw");
//         histsRaw["all"].push_back(new TH1F(("all_"+*itSample).c_str(),("all_"+*itSample).c_str(),maxCuts,0,maxCuts));
//         histsRaw["all"].back()->Sumw2();

        cout << "Counts: " << *itSample << " == " << "all" << endl;
        std::map<size_t,int>::const_iterator ite = counterEvents["all"].begin();
        for(;ite!=counterEvents["all"].end();++ite) {
//             for(int j=0;j<ite->second;++j) {
//                 outputFile->cd("raw");
//                 histsRaw["all"].back()->Fill(ite->first);
//             }
            cout << setw(3) << ite->first << setw(7) << ite->second << setw(12) << setprecision(6) << ite->second * scaleTo1pb1 << endl;
        }

        for(vector<string>::const_iterator itHypo=hypoTypes.begin();itHypo!=hypoTypes.end();++itHypo) {

//             outputFile->cd("raw");
//             histsRaw[*itHypo].push_back(new TH1F((*itHypo+"_"+*itSample).c_str(),(*itHypo+"_"+*itSample).c_str(),maxCuts,0,maxCuts));
//             histsRaw[*itHypo].back()->Sumw2();

            cout << "Counts: " << *itSample << " == " << *itHypo << endl;
            for(ite = counterEvents[*itHypo].begin();ite!=counterEvents[*itHypo].end();++ite) {
//                 for(int j=0;j<ite->second;++j) {
//                     outputFile->cd("raw");
//                     histsRaw[*itHypo].back()->Fill(ite->first);
//                 }
                cout << setw(3) << ite->first << setw(7) << ite->second << setw(12) << setprecision(6) << ite->second * scaleTo1pb1 << endl;
            }
        }

        vector<EvtSummary>::const_iterator itEvt;
//         if(input.getParameter<bool>("printEvents")) {
//             cout << "Printing event list for: " << *itSample << endl;
//             cout << "===================================" << endl;
//             for(itEvt=eventList.begin();itEvt!=eventList.end();++itEvt) 
//                 cout << *itEvt;
//                 for(int k=0;k<=itEvt->cut;k++) 
//                     printf("%03d ",k);
//                 cout << endl;
//         }
        if(input.getParameter<bool>("printEvents")) {
            cout << "Printing event list for: " << *itSample << endl;
            cout << "===================================" << endl;
            for(itEvt=eventList.begin();itEvt!=eventList.end();++itEvt) 
                if(itEvt->cut == 14) cout << *itEvt << endl;
        }

    } //end loop over input datasets

//     vector<TH1F*>::const_iterator itHist;
//     outputFile->cd("raw");
//     for(itHist=histsRaw["all"].begin();itHist!=histsRaw["all"].end();++itHist) {
//         setAxisLabels(*itHist);
//         (*itHist)->LabelsOption("v");
//         (*itHist)->Write((*itHist)->GetName(),TObject::kOverwrite);
//     }
//     for(vector<string>::const_iterator itHypo=hypoTypes.begin();itHypo!=hypoTypes.end();++itHypo) {
//         for(itHist=histsRaw[*itHypo].begin();itHist!=histsRaw[*itHypo].end();++itHist) {
//             setAxisLabels(*itHist);
//             (*itHist)->LabelsOption("v");
//             (*itHist)->Write((*itHist)->GetName(),TObject::kOverwrite);
//         }
//     }

//     string folder(""), postfix("");
//     for(vector<string>::const_iterator itSample=sampleInputs.begin();itSample!=sampleInputs.end();++itSample) {
//         for(vector<string>::const_iterator itHypo=hypoTypes.begin();itHypo!=hypoTypes.end();++itHypo) {
//             for(int i=0;i<maxCuts;++i) {
//                 os.clear(); os.str("");
//                 os << "_" << setw(2) << setfill('0') << i;
//                 postfix = os.str();
//                 os.clear(); os.str("");
//                 os << *itSample << "/" << *itHypo << "/byCut";
//                 folder =  os.str();
//                 writeHistsToFile(outputFile,folder,postfix,cutByCutHists[*itSample][*itHypo][i]);
//                 os.clear(); os.str("");
//                 os << *itSample << "/" << *itHypo << "/nMinus1";
//                 folder =  os.str();
//                 writeHistsToFile(outputFile,folder,postfix,nMinus1Hists[*itSample][*itHypo][i]);
//             }
//         }
//         for(int i=0;i<maxCuts;++i) {
//             os.clear(); os.str("");
//             os << "_" << setw(2) << setfill('0') << i;
//             postfix = os.str();
//             os.clear(); os.str("");
//             os << *itSample << "/" << "all" << "/byCut";
//             folder =  os.str();
//             writeHistsToFile(outputFile,folder,postfix,cutByCutHists[*itSample]["all"][i]);
//             os.clear(); os.str("");
//             os << *itSample << "/" << "all" << "/nMinus1";
//             folder =  os.str();
//             writeHistsToFile(outputFile,folder,postfix,nMinus1Hists[*itSample]["all"][i]);
//         }
//     }


//     outputFile->Close();

    return 0;
}


void fillHists(const reco::SkimEvent& skimEvent,std::vector<TH1F*> &hists) {

    for(std::vector<TH1F*>::iterator itHist=hists.begin();itHist!=hists.end();++itHist) {
        (*itHist)->Fill(StringObjectFunction<reco::SkimEvent>( (*itHist)->GetTitle() )(skimEvent));
    }


}


void setupHists(const edm::ParameterSet& pset,std::vector<TH1F*> &hists, const std::string& postfix) {
    
    vector<string> plots = pset.getParameterNamesForType<edm::ParameterSet>();
    for(vector<string>::const_iterator itPlot=plots.begin();itPlot!=plots.end();++itPlot) {
        hists.push_back( getNewHistFromPSet(pset.getParameter<edm::ParameterSet>(*itPlot), postfix) );
    }
        

}

void writeHistsToFile(TFile *f, const std::string &folder,const std::string &postfix,std::vector<TH1F*> &hists) {

//     f->mkdir(folder.c_str());
    f->cd(folder.c_str());
    for(std::vector<TH1F*>::const_iterator itHist=hists.begin();itHist!=hists.end();++itHist) {
        (*itHist)->Write( ((*itHist)->GetTitle()+postfix).c_str() );
    }

}

TH1F *getNewHistFromPSet( const edm::ParameterSet &pset, const std::string &postfix ) {
    TH1F *temp = new TH1F(
        (pset.getParameter<string>("variable")+postfix).c_str(),
        (pset.getParameter<string>("variable")).c_str(),
        pset.getParameter<uint>("nbins"),
        pset.getParameter<double>("low"),
        pset.getParameter<double>("high")
    );
    temp->GetXaxis()->SetTitle( pset.getParameter<string>("xtitle").c_str() );
    return temp;
}


void setAxisLabels(const TH1 *h) {
    h->GetXaxis()->SetBinLabel( 1,"Initial Skim");
    h->GetXaxis()->SetBinLabel( 2,"OS,20/10,Fid");
    h->GetXaxis()->SetBinLabel( 3,"ID");
    h->GetXaxis()->SetBinLabel( 4,"Iso");
    h->GetXaxis()->SetBinLabel( 5,"Conversion");
    h->GetXaxis()->SetBinLabel( 6,"IP");
    h->GetXaxis()->SetBinLabel( 7,"Extra Lepton");
    h->GetXaxis()->SetBinLabel( 8,"#{E}_{T}");
    h->GetXaxis()->SetBinLabel( 9,"M_{ll} 1");
    h->GetXaxis()->SetBinLabel(10,"Z Veto");
    h->GetXaxis()->SetBinLabel(11,"p#slash{E}_{T}");
    h->GetXaxis()->SetBinLabel(12,"Jet Veto");
    h->GetXaxis()->SetBinLabel(13,"Soft Muon");
    h->GetXaxis()->SetBinLabel(14,"Top Veto");
    h->GetXaxis()->SetBinLabel(15,"M_{ll} 2");
    h->GetXaxis()->SetBinLabel(16,"p_{T} Max");
    h->GetXaxis()->SetBinLabel(17,"p_{T} Min");
    h->GetXaxis()->SetBinLabel(18,"d#Phi_{ll}");
}

