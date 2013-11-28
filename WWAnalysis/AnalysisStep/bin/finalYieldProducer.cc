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

const int maxCuts = 1;
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

        int evtCount = 0;

        outputFile->cd("scaled");
        histsRaw["wwelel"].push_back(new TH1F(("wwelel_"+*itSample).c_str(),("wwelel_"+*itSample).c_str(),maxCuts,0.5,maxCuts+0.5));
        histsRaw["wwelmu"].push_back(new TH1F(("wwelmu_"+*itSample).c_str(),("wwelmu_"+*itSample).c_str(),maxCuts,0.5,maxCuts+0.5));
        histsRaw["wwmumu"].push_back(new TH1F(("wwmumu_"+*itSample).c_str(),("wwmumu_"+*itSample).c_str(),maxCuts,0.5,maxCuts+0.5));
        histsRaw["all"].push_back(new TH1F(("all_"+*itSample).c_str(),("all_"+*itSample).c_str(),maxCuts,0.5,maxCuts+0.5));

        histsRaw["wwelel"].back()->Sumw2();
        histsRaw["wwelmu"].back()->Sumw2();
        histsRaw["wwmumu"].back()->Sumw2();
        histsRaw["all"].back()->Sumw2();

        outputFile->cd("raw");
        histsScaled["wwelel"].push_back(new TH1F(("wwelel_"+*itSample).c_str(),("wwelel_"+*itSample).c_str(),maxCuts,0.5,maxCuts+0.5));
        histsScaled["wwelmu"].push_back(new TH1F(("wwelmu_"+*itSample).c_str(),("wwelmu_"+*itSample).c_str(),maxCuts,0.5,maxCuts+0.5));
        histsScaled["wwmumu"].push_back(new TH1F(("wwmumu_"+*itSample).c_str(),("wwmumu_"+*itSample).c_str(),maxCuts,0.5,maxCuts+0.5));
        histsScaled["all"].push_back(new TH1F(("all_"+*itSample).c_str(),("all_"+*itSample).c_str(),maxCuts,0.5,maxCuts+0.5));

        histsScaled["wwelel"].back()->Sumw2();
        histsScaled["wwelmu"].back()->Sumw2();
        histsScaled["wwmumu"].back()->Sumw2();
        histsScaled["all"].back()->Sumw2();

        //loop on events
        if(ev.size()) for( ev.toBegin(); ! ev.atEnd(); ++ev) if(myLumiSel(ev)) {

            if( evtCount++%1000 == 0 ) cerr << "Processing " << *itSample << ": " << setw(10) << evtCount << endl;


            //loop on hypothesis
            double thisScale1 = 0 , thisScale2 = 0;
            map<string,bool> pass;
            pass["wwelel"] = false;
            pass["wwelmu"] = false;
            pass["wwmumu"] = false;
            pass["all"] = false;
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

                //loop on no of hypos in event
                for(vector<reco::SkimEvent>::const_iterator mySkimEvent = skimEventH.ptr()->begin();
                        mySkimEvent != skimEventH.ptr()->end(); mySkimEvent++){


                    if( mySkimEvent->q(0)*mySkimEvent->q(1) >= 0 ) continue;
                    if( mySkimEvent->isSTA(0) ) continue;
                    if( mySkimEvent->isSTA(1) ) continue;
                    if( mySkimEvent->ptMin() <= ptMin ) continue;
                    if( mySkimEvent->ptMax() <= ptMax ) continue;
                    if( !(mySkimEvent->leptEtaCut(etaMu, etaEl)) ) continue;
                    if( ! mySkimEvent->hasGoodVertex() ) continue;
                    if( fabs(mySkimEvent->d0Reco(0)) >= d0 ) continue;
                    if( fabs(mySkimEvent->d0Reco(1)) >= d0 ) continue;
                    if( fabs(mySkimEvent->dZReco(0)) >= dZ ) continue;
                    if( fabs(mySkimEvent->dZReco(1)) >= dZ ) continue;
                    if( abs(mySkimEvent->pdgId(0)) == 11 && mySkimEvent->allIso(0)/mySkimEvent->pt(0) >= isoEl ) continue;
                    if( abs(mySkimEvent->pdgId(1)) == 11 && mySkimEvent->allIso(1)/mySkimEvent->pt(1) >= isoEl ) continue;
                    if( abs(mySkimEvent->pdgId(0)) == 13 && mySkimEvent->allIso(0)/mySkimEvent->pt(0) >= isoMu ) continue;
                    if( abs(mySkimEvent->pdgId(1)) == 13 && mySkimEvent->allIso(1)/mySkimEvent->pt(1) >= isoMu ) continue;
//                     if( !(mySkimEvent->passesIDV1(0) && mySkimEvent->passesIDV1(1)) ) continue;
                    if( !(mySkimEvent->passesConversion(0) && mySkimEvent->passesConversion(1)) ) continue;
                    if( mySkimEvent->tcMet() <= met) continue;
                    if(!(mySkimEvent->mll()>mll) ) continue;
                    if(!( fabs(mySkimEvent->mll()-91.1876)>mZ) ) continue;
                    if( mySkimEvent->projTcMet() <= pMet ) continue;
                    if( !(mySkimEvent->nCentralJets( jetPt, jetEta, useJEC) <= nCentralJet) ) continue;
                    if(mySkimEvent->nSoftMu(3.)>nSoftMu ) continue;
                    if(mySkimEvent->nExtraLep(ptMin)>nExtraLep ) continue;
                    if(mySkimEvent->bTaggedJetsUnder( jetPt, bValue) >nBtagJets ) continue;
                    if(mySkimEvent->mll()>=mllMaxFinal ) continue;
                    if(mySkimEvent->ptMax() < ptMaxFinal ) continue;
                    if(mySkimEvent->ptMin() < ptMinFinal ) continue;
                    if(mySkimEvent->projTcMet() < pMetFinal ) continue;
                    if(mySkimEvent->tcMet()/mySkimEvent->pTll() < metOverPtLL ) continue;
                    if(mySkimEvent->dPhill() >= deltaPhiLL ) continue;

                    pass[*itHypo]=true;
                    pass["all"]=true;
                    thisScale1 = scaleFactor(mySkimEvent->ptByPt(0), mySkimEvent->etaByPt(0), mySkimEvent->ptByPt(1), mySkimEvent->etaByPt(1), mySkimEvent->mitType(), 0);
                    thisScale2 = scaleFactor(mySkimEvent->ptByPt(0), mySkimEvent->etaByPt(0), mySkimEvent->ptByPt(1), mySkimEvent->etaByPt(1), mySkimEvent->nGoodVertices(), 2);


                } //end of loop over SkimEvent

            } // end loop over hypothesis

            double thisScale = 1;
            if(input.getParameter<bool>("isMC")) {
                thisScale *= thisScale1;
                if(input.getParameter<bool>("isSignal")) thisScale *= 0.945;
                else                                     thisScale *= thisScale2;
            }

            if(thisScale) cout << "This is the scale: " << thisScale << endl;
            if(thisScale && !pass["all"]) cerr << "oh shit man" << endl;
            if(pass["wwelel"]) histsRaw["wwelel"].back()->Fill(1,1);
            if(pass["wwelmu"]) histsRaw["wwelmu"].back()->Fill(1,1);
            if(pass["wwmumu"]) histsRaw["wwmumu"].back()->Fill(1,1);
            if(pass["all"   ]) histsRaw["all"   ].back()->Fill(1,1);
            if(pass["wwelel"]) histsScaled["wwelel"].back()->Fill(1,thisScale);
            if(pass["wwelmu"]) histsScaled["wwelmu"].back()->Fill(1,thisScale);
            if(pass["wwmumu"]) histsScaled["wwmumu"].back()->Fill(1,thisScale);
            if(pass["all"   ]) histsScaled["all"   ].back()->Fill(1,thisScale);

        } // end loop over edm::events 

    } //end loop over input datasets

    vector<TH1F*>::const_iterator itHist;

    outputFile->cd("scaled");
    for(itHist=histsScaled["all"].begin();itHist!=histsScaled["all"].end();++itHist) {
        (*itHist)->Write((*itHist)->GetName(),TObject::kOverwrite);
    }
    outputFile->cd("raw");
    for(itHist=histsRaw["all"].begin();itHist!=histsRaw["all"].end();++itHist) {
        (*itHist)->Write((*itHist)->GetName(),TObject::kOverwrite);
    }
    outputFile->cd("scaled");
    for(vector<string>::const_iterator itHypo=hypoTypes.begin();itHypo!=hypoTypes.end();++itHypo) {
        for(itHist=histsScaled[*itHypo].begin();itHist!=histsScaled[*itHypo].end();++itHist) {
            (*itHist)->Write((*itHist)->GetName(),TObject::kOverwrite);
        }
    }
    outputFile->cd("raw");
    for(vector<string>::const_iterator itHypo=hypoTypes.begin();itHypo!=hypoTypes.end();++itHypo) {
        for(itHist=histsRaw[*itHypo].begin();itHist!=histsRaw[*itHypo].end();++itHist) {
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
