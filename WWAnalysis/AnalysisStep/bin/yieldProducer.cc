#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "PhysicsTools/SelectorUtils/interface/RunLumiSelector.h"

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

#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
using namespace std;

const int maxCuts = 30;

void setAxisLabels(const TH1 *h) {
    h->GetXaxis()->SetBinLabel( 1,"Initial Skim");
    h->GetXaxis()->SetBinLabel( 2,"20/10");
    h->GetXaxis()->SetBinLabel( 3,"20/20");
    h->GetXaxis()->SetBinLabel( 4,"D0");
    h->GetXaxis()->SetBinLabel( 5,"DZ");
    h->GetXaxis()->SetBinLabel( 6,"Iso");
    h->GetXaxis()->SetBinLabel( 7,"ID");
    h->GetXaxis()->SetBinLabel( 8,"Conversion");
    h->GetXaxis()->SetBinLabel( 9,"#slash{E}_{T}");
    h->GetXaxis()->SetBinLabel(10,"M_{ll}");
    h->GetXaxis()->SetBinLabel(11,"Z Veto");
    h->GetXaxis()->SetBinLabel(12,"p#slash{E}_{T}");
    h->GetXaxis()->SetBinLabel(13,"Jet Veto");
    h->GetXaxis()->SetBinLabel(14,"Soft Muon");
    h->GetXaxis()->SetBinLabel(15,"Extra Lepton");
    h->GetXaxis()->SetBinLabel(16,"Top Veto");
    h->GetXaxis()->SetBinLabel(17,"M_{ll} < ");
    h->GetXaxis()->SetBinLabel(18,"p_{T} Max");
    h->GetXaxis()->SetBinLabel(19,"p_{T} Min");
    h->GetXaxis()->SetBinLabel(20,"pMet 2");
    h->GetXaxis()->SetBinLabel(21,"d#Phi_{ll}");
}

struct EvtSummary {
    unsigned int run;
    unsigned int lumi;
    unsigned int evt;
    size_t cut;
    string hypo;

    friend ostream& operator<<(ostream& out, const EvtSummary &evt) {
        out << setw(9) << evt.run << setw(15) << evt.lumi << setw(20) << evt.evt << setw(10) << evt.hypo << setw(4) << evt.cut;
        return out;
    }

    bool operator==(const EvtSummary& rhs) const {
        return ( run==rhs.run && lumi==rhs.lumi && evt==rhs.evt );
    }

    EvtSummary(unsigned int r, unsigned int l, unsigned int e, size_t c, string h) :
        run(r), lumi(l), evt(e), cut(c), hypo(h) { }
};

double scaleFactor(double pt1, double eta1, double pt2, double eta2, int type, int nsel);

int main(int argc,char* argv[]) {

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
    TFile *outputFile = TFile::Open(outputFileName.c_str(),"UPDATE");
    if(!outputFile || outputFile->IsZombie() ) {
        cerr << "Dagnabit, your output file sucks ..." << endl;
        return 1;
    }
    outputFile->mkdir("raw");
    outputFile->mkdir("scaled");

//     double scaleToData = allPars.getParameter<double>("integratedLumi");
    reco::SkimEvent::setupJEC(
        allPars.getParameter<string>("l2File"),
        allPars.getParameter<string>("l3File"),
        allPars.getParameter<string>("resFile")
    );

    RunLumiSelector myLumiSel(allPars);

    map<string,vector<TH1F*> > histsScaled;
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

    //loop on datasets
    edm::ParameterSet sampleInputParams = allPars.getParameter<edm::ParameterSet>("inputParams");
    vector<string> sampleInputs = sampleInputParams.getParameterNamesForType<edm::ParameterSet>();
    for(vector<string>::const_iterator itSample=sampleInputs.begin();itSample!=sampleInputs.end();++itSample) {

        edm::ParameterSet input = sampleInputParams.getParameter<edm::ParameterSet>(*itSample);
        fwlite::ChainEvent ev(input.getParameter<vector<string> >("files"));

        map<string,map<size_t,int> > counterEvents;
        map<string,map<size_t,bool> > passer;
        vector<EvtSummary> eventList;

        int evtCount = 0;
        //loop on events
        if(ev.size()) for( ev.toBegin(); ! ev.atEnd(); ++ev) if(myLumiSel(ev)) {

            if( evtCount++%1000 == 0 ) cerr << "Processing " << *itSample << ": " << setw(10) << evtCount << endl;


            //loop on hypothesis
            passer["all"].clear();
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

                passer[*itHypo].clear();
                //loop on no of hypos in event
                for(vector<reco::SkimEvent>::const_iterator mySkimEvent = skimEventH.ptr()->begin();
                        mySkimEvent != skimEventH.ptr()->end(); mySkimEvent++){


                    int i=0;
                    passer[*itHypo][i++]=true;
//                     if( vetoMuTriggered && ( mySkimEvent->isMuTriggered(0) || mySkimEvent->isMuTriggered(1) ) ) continue;
//                     if( !vetoMuTriggered && !mySkimEvent->isMuTriggered(0) && !mySkimEvent->isMuTriggered(1) ) continue;
//                     if( !mySkimEvent->isEcalSeeded(0) ) continue;
//                     if( !mySkimEvent->isEcalSeeded(1) ) continue;
                    if( mySkimEvent->q(0)*mySkimEvent->q(1) >= 0 ) continue;
                    if( mySkimEvent->isSTA(0) ) continue;
                    if( mySkimEvent->isSTA(1) ) continue;
                    if( mySkimEvent->ptMin() <= ptMin ) continue;
                    if( mySkimEvent->ptMax() <= ptMax ) continue;
                    if( !(mySkimEvent->leptEtaCut(etaMu, etaEl)) ) continue;
                    passer[*itHypo][i++]=true;

                    if( ! mySkimEvent->hasGoodVertex() ) continue;
                    passer[*itHypo][i++]=true;

                    if( fabs(mySkimEvent->d0Reco(0)) >= d0 ) continue;
                    if( fabs(mySkimEvent->d0Reco(1)) >= d0 ) continue;
                    passer[*itHypo][i++]=true;

                    if( fabs(mySkimEvent->dZReco(0)) >= dZ ) continue;
                    if( fabs(mySkimEvent->dZReco(1)) >= dZ ) continue;
                    passer[*itHypo][i++]=true;

                    if( abs(mySkimEvent->pdgId(0)) == 11 && mySkimEvent->allIso(0)/mySkimEvent->pt(0) >= isoEl ) continue;
                    if( abs(mySkimEvent->pdgId(1)) == 11 && mySkimEvent->allIso(1)/mySkimEvent->pt(1) >= isoEl ) continue;
                    if( abs(mySkimEvent->pdgId(0)) == 13 && mySkimEvent->allIso(0)/mySkimEvent->pt(0) >= isoMu ) continue;
                    if( abs(mySkimEvent->pdgId(1)) == 13 && mySkimEvent->allIso(1)/mySkimEvent->pt(1) >= isoMu ) continue;
                    passer[*itHypo][i++]=true;

//                     if( !(mySkimEvent->passesIDV1(0) && mySkimEvent->passesIDV1(1)) ) continue;
//                     passer[*itHypo][i++]=true;

                    if( !(mySkimEvent->passesConversion(0) && mySkimEvent->passesConversion(1)) ) continue;
                    passer[*itHypo][i++]=true;

                    if( mySkimEvent->tcMet() <= met) continue;
                    passer[*itHypo][i++]=true;

//                     if( !(mySkimEvent->triggerMatchingCut()) ) continue;
//                     passer[*itHypo][i++]=true;

                    if(!(mySkimEvent->mll()>mll) ) continue;
                    passer[*itHypo][i++]=true;

                    if(!( fabs(mySkimEvent->mll()-91.1876)>mZ) ) continue;
                    passer[*itHypo][i++]=true;

                    if( mySkimEvent->projTcMet() <= pMet ) continue;
                    passer[*itHypo][i++]=true;

                    if( !(mySkimEvent->nCentralJets( jetPt, jetEta, useJEC) <= nCentralJet) ) continue;
                    passer[*itHypo][i++]=true;

                    if(mySkimEvent->nSoftMu(3.)>nSoftMu ) continue;
                    passer[*itHypo][i++]=true;

                    if(mySkimEvent->nExtraLep(ptMin)>nExtraLep ) continue;
                    passer[*itHypo][i++]=true;

                    if(mySkimEvent->bTaggedJetsUnder( jetPt, bValue) >nBtagJets ) continue;
                    passer[*itHypo][i++]=true;

                    //final selection cuts
                    if(mySkimEvent->mll()>=mllMaxFinal ) continue;
                    passer[*itHypo][i++]=true;

                    if(mySkimEvent->ptMax() < ptMaxFinal ) continue;
                    passer[*itHypo][i++]=true;

                    if(mySkimEvent->ptMin() < ptMinFinal ) continue;
                    passer[*itHypo][i++]=true;

                    if(mySkimEvent->projTcMet() < pMetFinal ) continue;
                    passer[*itHypo][i++]=true;

                    if(mySkimEvent->tcMet()/mySkimEvent->pTll() < metOverPtLL ) continue;
                    passer[*itHypo][i++]=true;

                    if(mySkimEvent->dPhill() >= deltaPhiLL ) continue;
                    passer[*itHypo][i++]=true;


                } //end of loop over SkimEvent

                EvtSummary tempEvt(
                        ev.eventAuxiliary().run(),
                        ev.eventAuxiliary().luminosityBlock(),
                        ev.eventAuxiliary().event(), 0, *itHypo
                        );
//                 vector<EvtSummary>::iterator myEvt = find(eventList.begin(),eventList.end(),tempEvt);
                vector<EvtSummary>::iterator myEvt;
//                 if( myEvt == eventList.end() ) {
                eventList.push_back(tempEvt);
                myEvt = eventList.end()-1;
//                 }

                std::map<size_t,bool>::const_iterator itp = passer[*itHypo].begin();
                for(;itp!=passer[*itHypo].end();++itp) {
                    if(itp->second) {
                        counterEvents[*itHypo][itp->first]++;
// 	                fillHistos(*itHypo,met[itp->first],projMet[itp->first],ptmax[itp->first],ptmin[itp->first],mll[itp->first],
// 		                dMllMz[itp->first],nSoftMu[itp->first],nJet[itp->first],metOverPTll[itp->first],dPhill[itp->first],
// 		                d0,d1,d2,d3,d4,d5,i6,i7,d8,d9);
                        passer["all"][itp->first] = true;
                    } 
                    if(itp->first > myEvt->cut) {
                        myEvt->cut = itp->first;
                        myEvt->hypo = *itHypo;
                    }
                }
            } // end loop over hypothesis

            std::map<size_t,bool>::const_iterator itp = passer["all"].begin();
            for(;itp!=passer["all"].end();++itp) {
                if(itp->second) {
                    //If I get here, itp->first represents which cut it passed
// 	                fillHistos("all",met[itp->first],projMet[itp->first],ptmax[itp->first],ptmin[itp->first],mll[itp->first],
// 		                dMllMz[itp->first],nSoftMu[itp->first],nJet[itp->first],metOverPTll[itp->first],dPhill[itp->first],
// 		                d0,d1,d2,d3,d4,d5,i6,i7,d8,d9);
                    counterEvents["all"][itp->first]++;
                } 
            }

        } // end loop over edm::events 

        double scaleTo1pb1 = input.getParameter<double>("scale");

        outputFile->cd("scaled");
        histsScaled["all"].push_back(new TH1F(("all_"+*itSample).c_str(),("all_"+*itSample).c_str(),maxCuts,0,maxCuts));
        histsScaled["all"].back()->Sumw2();
        outputFile->cd("raw");
        histsRaw["all"].push_back(new TH1F(("all_"+*itSample).c_str(),("all_"+*itSample).c_str(),maxCuts,0,maxCuts));
        histsRaw["all"].back()->Sumw2();

        std::map<size_t,int>::const_iterator ite = counterEvents["all"].begin();
        for(;ite!=counterEvents["all"].end();++ite) {
            for(int j=0;j<ite->second;++j) {
                outputFile->cd("scaled");
                histsScaled["all"].back()->Fill(ite->first,scaleTo1pb1);
                outputFile->cd("raw");
                histsRaw["all"].back()->Fill(ite->first);
            }
        }

        for(vector<string>::const_iterator itHypo=hypoTypes.begin();itHypo!=hypoTypes.end();++itHypo) {

            outputFile->cd("scaled");
            histsScaled[*itHypo].push_back(new TH1F((*itHypo+"_"+*itSample).c_str(),(*itHypo+"_"+*itSample).c_str(),maxCuts,0,maxCuts));
            histsScaled[*itHypo].back()->Sumw2();
            outputFile->cd("raw");
            histsRaw[*itHypo].push_back(new TH1F((*itHypo+"_"+*itSample).c_str(),(*itHypo+"_"+*itSample).c_str(),maxCuts,0,maxCuts));
            histsRaw[*itHypo].back()->Sumw2();

            for(ite = counterEvents[*itHypo].begin();ite!=counterEvents[*itHypo].end();++ite) {
                for(int j=0;j<ite->second;++j) {
                    outputFile->cd("scaled");
                    histsScaled[*itHypo].back()->Fill(ite->first,scaleTo1pb1);
                    outputFile->cd("raw");
                    histsRaw[*itHypo].back()->Fill(ite->first);
                }
            }
        }

        vector<EvtSummary>::const_iterator itEvt;
        if(input.getParameter<bool>("printEvents")) {
            cout << "Printing event list for: " << *itSample << endl;
            cout << "===================================" << endl;
            for(itEvt=eventList.begin();itEvt!=eventList.end();++itEvt) 
                if(itEvt->cut == 21) cout << *itEvt << endl;
        }

    } //end loop over input datasets

    vector<TH1F*>::const_iterator itHist;
    outputFile->cd("scaled");
    for(itHist=histsScaled["all"].begin();itHist!=histsScaled["all"].end();++itHist) {
        setAxisLabels(*itHist);
        (*itHist)->LabelsOption("v");
        (*itHist)->Write((*itHist)->GetName(),TObject::kOverwrite);
    }
    outputFile->cd("raw");
    for(itHist=histsRaw["all"].begin();itHist!=histsRaw["all"].end();++itHist) {
        setAxisLabels(*itHist);
        (*itHist)->LabelsOption("v");
        (*itHist)->Write((*itHist)->GetName(),TObject::kOverwrite);
    }
    outputFile->cd("scaled");
    for(vector<string>::const_iterator itHypo=hypoTypes.begin();itHypo!=hypoTypes.end();++itHypo) {
        for(itHist=histsScaled[*itHypo].begin();itHist!=histsScaled[*itHypo].end();++itHist) {
            setAxisLabels(*itHist);
            (*itHist)->LabelsOption("v");
            (*itHist)->Write((*itHist)->GetName(),TObject::kOverwrite);
        }
    }
    outputFile->cd("raw");
    for(vector<string>::const_iterator itHypo=hypoTypes.begin();itHypo!=hypoTypes.end();++itHypo) {
        for(itHist=histsRaw[*itHypo].begin();itHist!=histsRaw[*itHypo].end();++itHist) {
            setAxisLabels(*itHist);
            (*itHist)->LabelsOption("v");
            (*itHist)->Write((*itHist)->GetName(),TObject::kOverwrite);
        }
    }

    outputFile->Close();

    return 0;
}


double scaleFactor(double pt1, double eta1, double pt2, double eta2, int type, int nsel){
  // type == 10/11/12/13 ->mm/ee/em/me
  // hardcoded, it's not used much
  int syst = 0;
  double scaleE[2] = {0.963, 0.977};
  if(syst == 1){
    scaleE[0] = scaleE[0] - 0.025;
    scaleE[1] = scaleE[1] - 0.029;
  };
  double scaleM[2] = {1.000, 1.000};
  if(syst == 1){
    scaleM[0] = scaleM[0] + 0.002;
    scaleM[1] = scaleM[1] + 0.002;
  };
  double weight = 1.0;
  
  if(nsel == 0){ // electron scale factor
    if     (type == 11){ 
      if(fabs(eta1) < 1.479) weight = weight * scaleE[0];
      else                   weight = weight * scaleE[1];
      if(fabs(eta2) < 1.479) weight = weight * scaleE[0];
      else                   weight = weight * scaleE[1];
    }
    else if(type == 12){ 
      if(fabs(eta1) < 1.479) weight = weight * scaleE[0];
      else                   weight = weight * scaleE[1];
    }
    else if(type == 13){ 
      if(fabs(eta2) < 1.479) weight = weight * scaleE[0];
      else                   weight = weight * scaleE[1];
    }
  }
  else if(nsel == 1){ // muon scale factor
    if     (type == 10){ 
      if(fabs(eta1) < 1.479) weight = weight * scaleM[0];
      else                   weight = weight * scaleM[1];
      if(fabs(eta2) < 1.479) weight = weight * scaleM[0];
      else                   weight = weight * scaleM[1];
    }
    else if(type == 12){ 
      if(fabs(eta2) < 1.479) weight = weight * scaleM[0];
      else                   weight = weight * scaleM[1];
    }
    else if(type == 13){ 
      if(fabs(eta1) < 1.479) weight = weight * scaleM[0];
      else                   weight = weight * scaleM[1];
    }
  }
  else if(nsel == 2){ // luminosity scale factor
    if     (type == 1) weight = weight * 1.4286;
    else if(type == 2) weight = weight * 0.9963;
    else if(type >= 3) weight = weight * 0.8106;
  }
  else if(nsel == 3){ // momentum scale factor
    if     (type == 10){ 
      if(pt1*0.99 < 20) weight = 0.0;
      if(pt2*0.99 < 20) weight = 0.0;
    }
    else if(type == 11){ 
      if     (fabs(eta1) <  1.479 && pt1*0.98 < 20) weight = 0.0;
      else if(fabs(eta1) >= 1.479 && pt1*0.96 < 20) weight = 0.0;
      if     (fabs(eta2) <  1.479 && pt2*0.98 < 20) weight = 0.0;
      else if(fabs(eta2) >= 1.479 && pt2*0.96 < 20) weight = 0.0;
    }
    else if(type == 12){ 
      if     (fabs(eta1) <  1.479 && pt1*0.98 < 20) weight = 0.0;
      else if(fabs(eta1) >= 1.479 && pt1*0.96 < 20) weight = 0.0;
      if(pt2*0.99 < 20) weight = 0.0;
    }
    else if(type == 13){ 
      if(pt1*0.99 < 20) weight = 0.0;
      if     (fabs(eta2) <  1.479 && pt2*0.98 < 20) weight = 0.0;
      else if(fabs(eta2) >= 1.479 && pt2*0.96 < 20) weight = 0.0;
    }
  }
  
  return weight;
}
