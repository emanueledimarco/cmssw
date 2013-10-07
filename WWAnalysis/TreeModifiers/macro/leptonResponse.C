#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1F.h"
#include "scales2.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"

#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooFitResult.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace RooFit;

struct DoubleCBParams {

    float mean;
    float sigma;
    float aL;
    float nL;
    float aR;
    float nR;

};

DoubleCBParams fitResponse(TTree* tree, float ptmin, float ptmax, float pzpmin, float pzpmax, std::string filename) {

    float boundlow = 0.8;
    float boundhigh = 1.2;

    if (ptmin<10 && pzpmin > 0.8945) {
        boundlow = 0.4;
        boundhigh = 2.0;
    }

    RooRealVar res ("res"    , "res"    , 1.0, boundlow, boundhigh);
    RooRealVar pt  ("genpt"  , "genpt"  , 0, 10000);
    RooRealVar p   ("genp"   , "genp"   , 0, 10000);
    RooRealVar pz  ("genpz"  , "genpz"  , 0, 10000);
    RooRealVar wgt ("wgt"    , "wgt"    , 0, 1000);

    RooRealVar mean  ("mean" , "mean" , 1.0 , 0.6, 1.4);
    RooRealVar sigma ("sigma", "sigma", 0.05, 0.0, 0.5);
    RooRealVar aL    ("aL"   , "aL"   , 1   , 0.0, 10);
    RooRealVar aR    ("aR"   , "aR"   , 1   , 0.0, 10);
    RooRealVar nL    ("nL"   , "nL"   , 1   , 0.0, 10);
    RooRealVar nR    ("nR"   , "nR"   , 1   , 0.0, 10);
    //RooRealVar nL    ("nL"   , "nL"   , 5);
    //RooRealVar nR    ("nR"   , "nR"   , 5);

    //nL.setConstant(kTRUE);
    //nR.setConstant(kTRUE);

    std::stringstream cutss;
    cutss << "genpt>" << ptmin << " && genpt<" << ptmax << " && " << "abs(genpz/genp)>" << pzpmin << " && abs(genpz/genp)<" << pzpmax; 

    RooDataSet dataset("dataset", "dataset", tree, RooArgSet(res, pt, p, pz, wgt), cutss.str().c_str(), "");

    RooDoubleCB pdf("pdf", "pdf", res, mean, sigma, aL, nL, aR, nR);
    pdf.fitTo(dataset, SumW2Error(kTRUE));

    RooPlot *frame = res.frame(50);
    dataset.plotOn(frame,RooFit::DataError(RooAbsData::SumW2));
    pdf.plotOn(frame);
    TCanvas c1;
    c1.cd();
    frame->Draw();
    c1.SaveAs((filename+".pdf").c_str());
    c1.SaveAs((filename+".png").c_str());


    DoubleCBParams dcbp;
    dcbp.mean  = mean.getVal();
    dcbp.sigma = sigma.getVal();
    dcbp.aL    = aL.getVal();
    dcbp.aR    = aR.getVal();
    dcbp.nL    = nL.getVal();
    dcbp.nR    = nR.getVal();
         
    return dcbp;

}

void leptonResponse() {

    init(false);    

    TFile* file  = new TFile("/home/avartak/CMS/Higgs/Moriond/CMSSW_5_3_3_patch3/src/WWAnalysis/AnalysisStep/trees/MuonFits/hzzTree.root");
    TTree* tree  = (TTree*)file->Get("electronTree/probe_tree");

    Float_t fres  = 0.;
    Float_t fwgt  = 0.;
    Float_t fp    = 0.;
    Float_t fpgen = 0.;
    Int_t   fnus  = 0 ;
    
    TBranch *bres  = tree->Branch("res", &fres, "res/F");
    TBranch *bwgt  = tree->Branch("wgt", &fwgt, "wgt/F");
    TBranch *bp    = tree->GetBranch("p");
    TBranch *bpgen = tree->GetBranch("genp");
    TBranch *bnus  = tree->GetBranch("event_nTruePV");
    
    bp->SetAddress(&fp);
    bpgen->SetAddress(&fpgen);
    bnus->SetAddress(&fnus);
    
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; i++){
        bp->GetEvent(i);
        bpgen->GetEvent(i);
        bnus->GetEvent(i);
        fres = fp/fpgen;
        fwgt = getPUWeight((float)fnus);
        //fwgt = 1.0;
        bres->Fill();
        bwgt->Fill();
    }


    float ptBins[12] = {5,7,10,15,20,25,30,35,40,50,70,100};
    float pzpBins[10] = {0,0.3799489623,0.6640367703,0.833654607,0.894540039,0.9163871674,0.9468060128,0.9704519366,0.9800963963,0.9866142982}; 
    float etaBins[10] = {0,0.4,0.8,1.2,1.4442,1.566,1.8,2.1,2.3,2.5};

    /*
    TH2F* hist_mean  = new TH2F("hist_mean" , "mean" , 11, ptBins, 9, etaBins);
    TH2F* hist_sigma = new TH2F("hist_sigma", "sigma", 11, ptBins, 9, etaBins);
    TH2F* hist_aL    = new TH2F("hist_aL"   , "aL"   , 11, ptBins, 9, etaBins);
    TH2F* hist_nL    = new TH2F("hist_nL"   , "nL"   , 11, ptBins, 9, etaBins);
    TH2F* hist_aR    = new TH2F("hist_aR"   , "aR"   , 11, ptBins, 9, etaBins);
    TH2F* hist_nR    = new TH2F("hist_nR"   , "nR"   , 11, ptBins, 9, etaBins);

    for (std::size_t i = 9; i < 11; i++) {
        for (std::size_t j = 0; j < 9; j++) {
            std::stringstream fnss;
            fnss << "fit_pt_" << i << "_" << i+1 << "_pzp_" << j << "_" << j+1; 
            DoubleCBParams param = fitResponse(tree, ptBins[i], ptBins[i+1], pzpBins[j], pzpBins[j+1], fnss.str());
            std::cout << "##############################################################" << std::endl;
            std::cout << "pT > " << ptBins[i] << " and pT < " << ptBins[i+1] << std::endl;    
            std::cout << "pz/p > " << pzpBins[j] << " and pz/p < " << pzpBins[j+1] << std::endl;    
            std::cout << "mean = " << param.mean << std::endl;
            std::cout << "sigma = " << param.sigma << std::endl;
            std::cout << "aL = " << param.aL << std::endl;
            std::cout << "nL = " << param.nL << std::endl;
            std::cout << "aR = " << param.aR << std::endl;
            std::cout << "nR = " << param.nR << std::endl;
            std::cout << std::endl;

            hist_mean ->SetBinContent(i+1, j+1, param.mean);
            hist_sigma->SetBinContent(i+1, j+1, param.sigma);
            hist_aL   ->SetBinContent(i+1, j+1, param.aL);
            hist_aR   ->SetBinContent(i+1, j+1, param.aR);
            hist_nL   ->SetBinContent(i+1, j+1, param.nL);
            hist_nR   ->SetBinContent(i+1, j+1, param.nR);
        }
    }

    std::stringstream electronoutss;
    electronoutss << "electronResponse2.root";

    TFile* electronoutfile = new TFile(electronoutss.str().c_str(), "RECREATE");
    hist_mean->Write();
    hist_sigma->Write();
    hist_aL->Write();
    hist_aR->Write();
    hist_nL->Write();
    hist_nR->Write();
    electronoutfile->Close();
    */

    fitResponse(tree, ptBins[10], ptBins[11], pzpBins[3], pzpBins[4], "fit");
}
