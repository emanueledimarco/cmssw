#ifndef FAKERATECALCULATOR_H
#define FAKERATECALCULATOR_H

#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TBranch.h>
#include <TLorentzVector.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <utility>
#include <sstream>
#include <cmath>
#include <algorithm>
#include "EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h"
#include "Muon/MuonAnalysisTools/interface/MuonEffectiveArea.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "scales2.h"

class FakeRateCalculator {

    private:
        std::vector<std::pair<int, int> > runeventinfo;    

        std::vector<float> mmbarrel;
        std::vector<float> mmendcap;
        std::vector<float> eebarrel;
        std::vector<float> eeendcap;
        std::vector<float> embarrel;
        std::vector<float> emendcap;
        std::vector<float> mebarrel;
        std::vector<float> meendcap;

        std::vector<float> mmbarrelerr;
        std::vector<float> mmendcaperr;
        std::vector<float> eebarrelerr;
        std::vector<float> eeendcaperr;
        std::vector<float> embarrelerr;
        std::vector<float> emendcaperr;
        std::vector<float> mebarrelerr;
        std::vector<float> meendcaperr;

        bool dofsr;

        int mode;

        float muboundary;
        float elboundary;

        TH2D* hfake0, *hfake1, *hfake2, *hfake3, *hall0, *hall1, *hall2, *hall3, *hpass0, *hpass1, *hpass2, *hpass3;

        float getMass(float l1pt, float l1eta, float l1phi, float l2pt, float l2eta, float l2phi) {
            float px1 = l1pt*cos(l1phi);
            float py1 = l1pt*sin(l1phi);
            float pz1 = l1pt*sinh(l1eta);
            float p1  = sqrt(px1*px1 + py1*py1 + pz1*pz1);
            
            float px2 = l2pt*cos(l2phi);
            float py2 = l2pt*sin(l2phi);
            float pz2 = l2pt*sinh(l2eta);
            float p2  = sqrt(px2*px2 + py2*py2 + pz2*pz2);
            
            TLorentzVector v1(px1, py1, pz1, p1);
            TLorentzVector v2(px2, py2, pz2, p2);
            
            v1 += v2;
            
            return v1.M();
        }

        float getMass2(float l1pt, float l1eta, float l1phi, float l2pt, float l2eta, float l2phi, float l3pt, float l3eta, float l3phi) {
            float px1 = l1pt*cos(l1phi);
            float py1 = l1pt*sin(l1phi);
            float pz1 = l1pt*sinh(l1eta);
            float p1  = sqrt(px1*px1 + py1*py1 + pz1*pz1);

            float px2 = l2pt*cos(l2phi);
            float py2 = l2pt*sin(l2phi);
            float pz2 = l2pt*sinh(l2eta);
            float p2  = sqrt(px2*px2 + py2*py2 + pz2*pz2);

            float px3 = l3pt*cos(l3phi);
            float py3 = l3pt*sin(l3phi);
            float pz3 = l3pt*sinh(l3eta);
            float p3  = sqrt(px3*px3 + py3*py3 + pz3*pz3);

            TLorentzVector v1(px1, py1, pz1, p1);
            TLorentzVector v2(px2, py2, pz2, p2);
            TLorentzVector v3(px3, py3, pz3, p3);

            v1 += v2;
            v1 += v3;

            return v1.M();
        }


        void dofakes(std::string path, bool do7TeV, float zmin, float zmax, float drcut, float mcut, int mde=0) {
            Float_t mubins[] = {0.0,10.0,15.0,20.0,30.0,40.0,50.0,80.0};
            //Float_t elbins[] = {0.0,10.0,15.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0};
            Float_t elbins[] = {0.0,10.0,20.0,30.0,40.0,50.0,60.0,80.0};
            UInt_t mubinnum = 7;
            UInt_t elbinnum = 7;
            Float_t ybins[] = {0.,1.,2.};
            UInt_t ybinnum = 2;

            hfake0 = new TH2D("hfake0", "", mubinnum, mubins, ybinnum, ybins);
            hall0  = new TH2D("hall0",  "", mubinnum, mubins, ybinnum, ybins);
            hpass0 = new TH2D("hpass0", "", mubinnum, mubins, ybinnum, ybins);

            hfake1 = new TH2D("hfake1", "", elbinnum, elbins, ybinnum, ybins);
            hall1  = new TH2D("hall1",  "", elbinnum, elbins, ybinnum, ybins);
            hpass1 = new TH2D("hpass1", "", elbinnum, elbins, ybinnum, ybins);
            
            hfake2 = new TH2D("hfake2", "", mubinnum, mubins, ybinnum, ybins);
            hall2  = new TH2D("hall2",  "", mubinnum, mubins, ybinnum, ybins);
            hpass2 = new TH2D("hpass2", "", mubinnum, mubins, ybinnum, ybins);

            hfake3 = new TH2D("hfake3", "", elbinnum, elbins, ybinnum, ybins);
            hall3  = new TH2D("hall3",  "", elbinnum, elbins, ybinnum, ybins);
            hpass3 = new TH2D("hpass3", "", elbinnum, elbins, ybinnum, ybins);
            
            hall0  ->Sumw2();
            hpass0 ->Sumw2();
            hfake0 ->Sumw2();
            hall1  ->Sumw2();
            hpass1 ->Sumw2();
            hfake1 ->Sumw2();
            hall2  ->Sumw2();
            hpass2 ->Sumw2();
            hfake2 ->Sumw2();
            hall3  ->Sumw2();
            hpass3 ->Sumw2();
            hfake3 ->Sumw2();

            domufakes(path, do7TeV, zmin, zmax, drcut, mcut);
            doelfakes(path, do7TeV, zmin, zmax, drcut, mcut);

            if (mde == 0) {
                hall0->Add(hall2);
                hpass0->Add(hpass2);

                hall1->Add(hall3);
                hpass1->Add(hpass3);

                hfake0 ->Divide(hpass0,  hall0);
                hfake1 ->Divide(hpass1,  hall1);

                for (std::size_t k = 0; k < mubinnum; k++) {
                    mmbarrel[k]    = hfake0->GetBinContent(k+1, 1);
                    mmendcap[k]    = hfake0->GetBinContent(k+1, 2);
                    mmbarrelerr[k] = hfake0->GetBinError(k+1, 1);
                    mmendcaperr[k] = hfake0->GetBinError(k+1, 2);
                    embarrel[k]    = hfake0->GetBinContent(k+1, 1);
                    emendcap[k]    = hfake0->GetBinContent(k+1, 2);
                    embarrelerr[k] = hfake0->GetBinError(k+1, 1);
                    emendcaperr[k] = hfake0->GetBinError(k+1, 2);
                }

                for (std::size_t k = 0; k < elbinnum; k++) {
                    eebarrel[k]    = hfake1->GetBinContent(k+1, 1);
                    eeendcap[k]    = hfake1->GetBinContent(k+1, 2);
                    eebarrelerr[k] = hfake1->GetBinError(k+1, 1);
                    eeendcaperr[k] = hfake1->GetBinError(k+1, 2);
                    mebarrel[k]    = hfake1->GetBinContent(k+1, 1);
                    meendcap[k]    = hfake1->GetBinContent(k+1, 2);
                    mebarrelerr[k] = hfake1->GetBinError(k+1, 1);
                    meendcaperr[k] = hfake1->GetBinError(k+1, 2);
                }
            }

            else {
                hfake0 ->Divide(hpass0,  hall0);
                hfake1 ->Divide(hpass1,  hall1);
                hfake2 ->Divide(hpass2,  hall2);
                hfake3 ->Divide(hpass3,  hall3);

                for (std::size_t k = 0; k < mubinnum; k++) {
                    mmbarrel[k]    = hfake0->GetBinContent(k+1, 1);
                    mmendcap[k]    = hfake0->GetBinContent(k+1, 2);
                    mmbarrelerr[k] = hfake0->GetBinError(k+1, 1);
                    mmendcaperr[k] = hfake0->GetBinError(k+1, 2);
                    embarrel[k]    = hfake2->GetBinContent(k+1, 1);
                    emendcap[k]    = hfake2->GetBinContent(k+1, 2);
                    embarrelerr[k] = hfake2->GetBinError(k+1, 1);
                    emendcaperr[k] = hfake2->GetBinError(k+1, 2);
                }

                for (std::size_t k = 0; k < elbinnum; k++) {
                    eebarrel[k]    = hfake1->GetBinContent(k+1, 1);
                    eeendcap[k]    = hfake1->GetBinContent(k+1, 2);
                    eebarrelerr[k] = hfake1->GetBinError(k+1, 1);
                    eeendcaperr[k] = hfake1->GetBinError(k+1, 2);
                    mebarrel[k]    = hfake3->GetBinContent(k+1, 1);
                    meendcap[k]    = hfake3->GetBinContent(k+1, 2);
                    mebarrelerr[k] = hfake3->GetBinError(k+1, 1);
                    meendcaperr[k] = hfake3->GetBinError(k+1, 2);
                }
            }
        }

        void domufakes(std::string path, bool do7TeV, float zmin, float zmax, float drcut, float mcut) {
            TFile* file = TFile::Open(path.c_str());
            TTree* tree = (TTree*)file->Get("zllmtree/probe_tree");
            
            TBranch *bmet     = tree->GetBranch("met");
            TBranch *bzmass   = tree->GetBranch("zmass");
            TBranch *bl1pt    = tree->GetBranch("l1pt");
            TBranch *bl1eta   = tree->GetBranch("l1eta");
            TBranch *bl1phi   = tree->GetBranch("l1phi");
            TBranch *bl1pdgId = tree->GetBranch("l1pdgId");
            TBranch *bl2pt    = tree->GetBranch("l2pt");
            TBranch *bl2eta   = tree->GetBranch("l2eta");
            TBranch *bl2phi   = tree->GetBranch("l2phi");
            TBranch *bl2pdgId = tree->GetBranch("l2pdgId");
            TBranch *bpt      = tree->GetBranch("pt");
            TBranch *beta     = tree->GetBranch("eta");
            TBranch *bphi     = tree->GetBranch("phi");
            TBranch *bid      = tree->GetBranch("newID");
            TBranch *biso     = dofsr ? tree->GetBranch("pfCombRelIso04EACorr_WithFSR") : tree->GetBranch("pfCombRelIso04EACorr");
            TBranch *brun     = tree->GetBranch("run");
            TBranch *blum     = tree->GetBranch("lumi");
            TBranch *bevt     = tree->GetBranch("event");
            
            float met     = 0.0;
            float l1pt    = 0.0;
            float l1eta   = 0.0;
            float l1phi   = 0.0;
            float l1pdgId = 0.0;
            float l2pt    = 0.0;
            float l2eta   = 0.0;
            float l2phi   = 0.0;
            float l2pdgId = 0.0;
            float zmass   = 0.0;
            float pt      = 0.0;
            float eta     = 0.0;
            float phi     = 0.0;
            float iso     = 0.0;
            int   id      = 0;
            int   run     = 0;
            int   lum     = 0;
            int   evt     = 0;
            
            bmet      ->SetAddress(&met);
            bzmass    ->SetAddress(&zmass);
            bl1pt     ->SetAddress(&l1pt);
            bl1eta    ->SetAddress(&l1eta);
            bl1phi    ->SetAddress(&l1phi);
            bl1pdgId  ->SetAddress(&l1pdgId);
            bl2pt     ->SetAddress(&l2pt);
            bl2eta    ->SetAddress(&l2eta);
            bl2phi    ->SetAddress(&l2phi);
            bl2pdgId  ->SetAddress(&l2pdgId);
            bpt       ->SetAddress(&pt);
            beta      ->SetAddress(&eta);
            bphi      ->SetAddress(&phi);
            bid       ->SetAddress(&id);
            biso      ->SetAddress(&iso);
            brun      ->SetAddress(&run);
            blum      ->SetAddress(&lum);
            bevt      ->SetAddress(&evt);
            
            for (int i = 0; i < tree->GetEntries(); i++) {
                bmet     ->GetEvent(i);
                bzmass   ->GetEvent(i);
                bl1pt    ->GetEvent(i);
                bl1eta   ->GetEvent(i);
                bl1phi   ->GetEvent(i);
                bl1pdgId ->GetEvent(i);
                bl2pt    ->GetEvent(i);
                bl2eta   ->GetEvent(i);
                bl2phi   ->GetEvent(i);
                bl2pdgId ->GetEvent(i);
                bpt      ->GetEvent(i);
                beta     ->GetEvent(i);
                bphi     ->GetEvent(i);
                bid      ->GetEvent(i);
                biso     ->GetEvent(i);
                brun     ->GetEvent(i);
                blum     ->GetEvent(i);
                bevt     ->GetEvent(i);
            
                bool existsAlready = false;
                for (std::size_t k = 0; k < runeventinfo.size(); k++) {
                    if (run == runeventinfo[k].first && evt == runeventinfo[k].second) existsAlready = true;
                }
       
                if (existsAlready) continue;            
                runeventinfo.push_back(std::pair<int, int>(run, evt));
 
                float dR1 = reco::deltaR(l1eta, l1phi, eta, phi); 
                float dR2 = reco::deltaR(l2eta, l2phi, eta, phi); 
            
                float m1 = getMass(l1pt, l1eta, l1phi, pt, eta, phi); 
                float m2 = getMass(l2pt, l2eta, l2phi, pt, eta, phi); 
            
                if (zmass>zmin && zmass<zmax && dR1 > drcut && dR2 > drcut && m1>mcut && m2>mcut && fabs(l1pdgId) == 13) {
                    float etafill = -1.;    
                    if (fabs(eta) < muboundary) etafill = 0.5;
                    else etafill = 1.5;
    
                    hall0->Fill(pt,etafill);
                    if (id>0 && iso<0.4) {
                        hpass0->Fill(pt,etafill);
                    }
                }
                if (zmass>zmin && zmass<zmax && dR1 > drcut && dR2 > drcut && m1>mcut && m2>mcut && fabs(l1pdgId) == 11) {
                    float etafill = -1.;    
                    if (fabs(eta) < muboundary) etafill = 0.5;
                    else etafill = 1.5;
    
                    hall2->Fill(pt,etafill);
                    if (id>0 && iso<0.4) {
                        hpass2->Fill(pt,etafill);
                    }
                }

            
            }
           
            file->Close();
            delete file;
        }


        void doelfakes(std::string path, bool do7TeV, float zmin, float zmax, float drcut, float mcut) {
            TFile* file = TFile::Open(path.c_str());
            TTree* tree = (TTree*)file->Get("zlletree/probe_tree");
            
            TBranch *bmet     = tree->GetBranch("met");
            TBranch *bzmass   = tree->GetBranch("zmass");
            TBranch *bl1pt    = tree->GetBranch("l1pt");
            TBranch *bl1eta   = tree->GetBranch("l1eta");
            TBranch *bl1phi   = tree->GetBranch("l1phi");
            TBranch *bl1pdgId = tree->GetBranch("l1pdgId");
            TBranch *bl2pt    = tree->GetBranch("l2pt");
            TBranch *bl2eta   = tree->GetBranch("l2eta");
            TBranch *bl2phi   = tree->GetBranch("l2phi");
            TBranch *bl2pdgId = tree->GetBranch("l2pdgId");
            TBranch *bpt      = tree->GetBranch("pt");
            TBranch *beta     = tree->GetBranch("eta");
            TBranch *bphi     = tree->GetBranch("phi");
            TBranch *bid      = tree->GetBranch("newID");
            TBranch *biso     = dofsr ? tree->GetBranch("pfCombRelIso04EACorr_WithFSR") : tree->GetBranch("pfCombRelIso04EACorr");
            TBranch *brun     = tree->GetBranch("run");
            TBranch *blum     = tree->GetBranch("lumi");
            TBranch *bevt     = tree->GetBranch("event");
           
            float met     = 0.0;
            float l1pt    = 0.0;
            float l1eta   = 0.0;
            float l1phi   = 0.0;
            float l1pdgId = 0.0;
            float l2pt    = 0.0;
            float l2eta   = 0.0;
            float l2phi   = 0.0;
            float l2pdgId = 0.0;
            float zmass   = 0.0;
            float pt      = 0.0;
            float eta     = 0.0;
            float phi     = 0.0;
            float iso     = 0.0;
            int   id      = 0;
            int   run     = 0;
            int   lum     = 0;
            int   evt     = 0;

            bmet      ->SetAddress(&met);
            bzmass    ->SetAddress(&zmass);
            bl1pt     ->SetAddress(&l1pt);
            bl1eta    ->SetAddress(&l1eta);
            bl1phi    ->SetAddress(&l1phi);
            bl1pdgId  ->SetAddress(&l1pdgId);
            bl2pt     ->SetAddress(&l2pt);
            bl2eta    ->SetAddress(&l2eta);
            bl2phi    ->SetAddress(&l2phi);
            bl2pdgId  ->SetAddress(&l2pdgId);
            bpt       ->SetAddress(&pt);
            beta      ->SetAddress(&eta);
            bphi      ->SetAddress(&phi);
            bid       ->SetAddress(&id);
            biso      ->SetAddress(&iso);
            brun      ->SetAddress(&run);
            blum      ->SetAddress(&lum);
            bevt      ->SetAddress(&evt);
 
            for (int i = 0; i < tree->GetEntries(); i++) {
                bmet     ->GetEvent(i);
                bzmass   ->GetEvent(i);
                bl1pt    ->GetEvent(i);
                bl1eta   ->GetEvent(i);
                bl1phi   ->GetEvent(i);
                bl1pdgId ->GetEvent(i);
                bl2pt    ->GetEvent(i);
                bl2eta   ->GetEvent(i);
                bl2phi   ->GetEvent(i);
                bl2pdgId ->GetEvent(i);
                bpt      ->GetEvent(i);
                beta     ->GetEvent(i);
                bphi     ->GetEvent(i);
                bid      ->GetEvent(i);
                biso     ->GetEvent(i);
                brun     ->GetEvent(i);
                blum     ->GetEvent(i);
                bevt     ->GetEvent(i);
            
                bool existsAlready = false;
                for (std::size_t k = 0; k < runeventinfo.size(); k++) {
                    if (run == runeventinfo[k].first && evt == runeventinfo[k].second) existsAlready = true;
                }
       
                if (existsAlready) continue;            
                runeventinfo.push_back(std::pair<int, int>(run, evt));
 
                float dR1 = reco::deltaR(l1eta, l1phi, eta, phi); 
                float dR2 = reco::deltaR(l2eta, l2phi, eta, phi); 
            
                float m1 = getMass(l1pt, l1eta, l1phi, pt, eta, phi); 
                float m2 = getMass(l2pt, l2eta, l2phi, pt, eta, phi); 
            
                if (zmass>zmin && zmass<zmax && dR1 > drcut && dR2 > drcut && m1>mcut && m2>mcut && fabs(l1pdgId) == 11) {
                    float etafill = -1.;    
                    if (fabs(eta) < elboundary) etafill = 0.5;
                    else etafill = 1.5;
    
                    hall1->Fill(pt,etafill);
                    if (id>0 && iso<0.4) {
                        hpass1->Fill(pt,etafill);
                    }
                }
                if (zmass>zmin && zmass<zmax && dR1 > drcut && dR2 > drcut && m1>mcut && m2>mcut && fabs(l1pdgId) == 13) {
                    float etafill = -1.;    
                    if (fabs(eta) < elboundary) etafill = 0.5;
                    else etafill = 1.5;
    
                    hall3->Fill(pt,etafill);
                    if (id>0 && iso<0.4) {
                        hpass3->Fill(pt,etafill);
                    }
                }
            
            }
           
            file->Close();
            delete file;
        }

    public:
        FakeRateCalculator():
            mmbarrel(std::vector<float>(11, 0.0)),
            mmendcap(std::vector<float>(11, 0.0)),
            eebarrel(std::vector<float>(11, 0.0)),
            eeendcap(std::vector<float>(11, 0.0)),
            embarrel(std::vector<float>(11, 0.0)),
            emendcap(std::vector<float>(11, 0.0)),
            mebarrel(std::vector<float>(11, 0.0)),
            meendcap(std::vector<float>(11, 0.0)),

            mmbarrelerr(std::vector<float>(11, 0.0)),
            mmendcaperr(std::vector<float>(11, 0.0)),
            eebarrelerr(std::vector<float>(11, 0.0)),
            eeendcaperr(std::vector<float>(11, 0.0)),
            embarrelerr(std::vector<float>(11, 0.0)),
            emendcaperr(std::vector<float>(11, 0.0)),
            mebarrelerr(std::vector<float>(11, 0.0)),
            meendcaperr(std::vector<float>(11, 0.0)),

            dofsr(false),
 
            mode(0),

            muboundary(1.479),
            elboundary(1.479) 
        {}

        FakeRateCalculator(std::string path, bool do7TeV, float zmin, float zmax, float drcut, float mcut, bool df):
            mmbarrel(std::vector<float>(11, 0.0)),
            mmendcap(std::vector<float>(11, 0.0)),
            eebarrel(std::vector<float>(11, 0.0)),
            eeendcap(std::vector<float>(11, 0.0)),
            embarrel(std::vector<float>(11, 0.0)),
            emendcap(std::vector<float>(11, 0.0)),
            mebarrel(std::vector<float>(11, 0.0)),
            meendcap(std::vector<float>(11, 0.0)),

            mmbarrelerr(std::vector<float>(11, 0.0)),
            mmendcaperr(std::vector<float>(11, 0.0)),
            eebarrelerr(std::vector<float>(11, 0.0)),
            eeendcaperr(std::vector<float>(11, 0.0)),
            embarrelerr(std::vector<float>(11, 0.0)),
            emendcaperr(std::vector<float>(11, 0.0)),
            mebarrelerr(std::vector<float>(11, 0.0)),
            meendcaperr(std::vector<float>(11, 0.0)),

            dofsr(df),

            mode(0),

            muboundary(1.479),
            elboundary(1.479) 
        {
            init(do7TeV);
            dofakes(path, do7TeV, zmin, zmax, drcut, mcut);
        }

        FakeRateCalculator(std::string path, bool do7TeV, float zmin, float zmax, float drcut, float mcut, bool df, int md):
            mmbarrel(std::vector<float>(11, 0.0)),
            mmendcap(std::vector<float>(11, 0.0)),
            eebarrel(std::vector<float>(11, 0.0)),
            eeendcap(std::vector<float>(11, 0.0)),
            embarrel(std::vector<float>(11, 0.0)),
            emendcap(std::vector<float>(11, 0.0)),
            mebarrel(std::vector<float>(11, 0.0)),
            meendcap(std::vector<float>(11, 0.0)),

            mmbarrelerr(std::vector<float>(11, 0.0)),
            mmendcaperr(std::vector<float>(11, 0.0)),
            eebarrelerr(std::vector<float>(11, 0.0)),
            eeendcaperr(std::vector<float>(11, 0.0)),
            embarrelerr(std::vector<float>(11, 0.0)),
            emendcaperr(std::vector<float>(11, 0.0)),
            mebarrelerr(std::vector<float>(11, 0.0)),
            meendcaperr(std::vector<float>(11, 0.0)),

            dofsr(df),

            mode(md),

            muboundary(1.479),
            elboundary(1.479) 
        {
            init(do7TeV);
            dofakes(path, do7TeV, zmin, zmax, drcut, mcut, mode);
        }

        float getPromptRate(float pt, float eta, float id) {
            return getPR(pt, eta, id);
        }

        float getFakeRate(float pt, float eta, float ch) {
            int bin = 0;
          
            if (pt>0  && pt<= 5) bin = 0;
            if (pt>5  && pt<= 7) bin = 1;
            if (pt>7  && pt<=10) bin = 2;
            if (pt>10 && pt<=15) bin = 3;
            if (pt>15 && pt<=20) bin = 4;
            if (pt>20 && pt<=25) bin = 5;
            if (pt>25 && pt<=30) bin = 6;
            if (pt>30 && pt<=40) bin = 7;
            if (pt>40 && pt<=50) bin = 8;
            if (pt>50 && pt<=80) bin = 9;
            if (pt>80)           bin = 9;

            
            if      (fabs(ch) == 0 && fabs(eta)<muboundary) return mmbarrel[bin];
            else if (fabs(ch) == 0 && fabs(eta)>muboundary) return mmendcap[bin];
            else if (fabs(ch) == 1 && fabs(eta)<elboundary) return eebarrel[bin];
            else if (fabs(ch) == 1 && fabs(eta)>elboundary) return eeendcap[bin];
            else if (fabs(ch) == 2 && fabs(eta)<muboundary) return embarrel[bin];
            else if (fabs(ch) == 2 && fabs(eta)>muboundary) return emendcap[bin];
            else if (fabs(ch) == 3 && fabs(eta)<elboundary) return mebarrel[bin];
            else if (fabs(ch) == 3 && fabs(eta)>elboundary) return meendcap[bin];
            else return 0.0;
            
        }

        float getFakeRateErr(float pt, float eta, float ch) {
            int bin = 0;

            if (pt>0  && pt<= 5) bin = 0;
            if (pt>5  && pt<= 7) bin = 1;
            if (pt>7  && pt<=10) bin = 2;
            if (pt>10 && pt<=15) bin = 3;
            if (pt>15 && pt<=20) bin = 4;
            if (pt>20 && pt<=25) bin = 5;
            if (pt>25 && pt<=30) bin = 6;
            if (pt>30 && pt<=40) bin = 7;
            if (pt>40 && pt<=50) bin = 8;
            if (pt>50 && pt<=80) bin = 9;
            if (pt>80)           bin = 9;

            
            if      (fabs(ch) == 0 && fabs(eta)<muboundary) return mmbarrelerr[bin];
            else if (fabs(ch) == 0 && fabs(eta)>muboundary) return mmendcaperr[bin];
            else if (fabs(ch) == 1 && fabs(eta)<elboundary) return eebarrelerr[bin];
            else if (fabs(ch) == 2 && fabs(eta)<elboundary) return embarrelerr[bin];
            else if (fabs(ch) == 2 && fabs(eta)>muboundary) return emendcaperr[bin];
            else if (fabs(ch) == 3 && fabs(eta)<elboundary) return mebarrelerr[bin];
            else if (fabs(ch) == 3 && fabs(eta)>elboundary) return meendcaperr[bin];
            else return 0.0;
            
        }

        void printMuonFakeRate(std::ostream& out) {

            out << "---- Muon Muon Barrel ----" << std::endl;
            for (int i = 0; i < hfake0->GetNbinsX(); i++) {
                out << "pT [ " << int(hfake0->GetBinLowEdge(i+1)) << " - " << int(hfake0->GetBinLowEdge(i+1))+int(hfake0->GetBinWidth(i+1)) << "] : " << mmbarrel[i]  << " +/- " << mmbarrelerr[i]  << endl;
            }
            out << std::endl;
            out << "---- Muon Muon Endcap ----" << std::endl;
            for (int i = 0; i < hfake0->GetNbinsX(); i++) {
                out << "pT [ " << int(hfake0->GetBinLowEdge(i+1)) << " - " << int(hfake0->GetBinLowEdge(i+1))+int(hfake0->GetBinWidth(i+1)) << "] : " << mmendcap[i]  << " +/- " << mmendcaperr[i]  << endl;
            }
            out << std::endl;

            out << "---- Electron Muon Barrel ----" << std::endl;
            for (int i = 0; i < hfake2->GetNbinsX(); i++) {
                out << "pT [ " << int(hfake2->GetBinLowEdge(i+1)) << " - " << int(hfake2->GetBinLowEdge(i+1))+int(hfake2->GetBinWidth(i+1)) << "] : " << embarrel[i]  << " +/- " << embarrelerr[i]  << endl;
            }
            out << std::endl;
            out << "---- Electron Muon Endcap ----" << std::endl;
            for (int i = 0; i < hfake2->GetNbinsX(); i++) {
                out << "pT [ " << int(hfake2->GetBinLowEdge(i+1)) << " - " << int(hfake2->GetBinLowEdge(i+1))+int(hfake2->GetBinWidth(i+1)) << "] : " << emendcap[i]  << " +/- " << emendcaperr[i]  << endl;
            }
            out << std::endl;

        }

        void printElectronFakeRate(std::ostream& out) {
            out << "---- Electron Electron Barrel ----" << std::endl;
            for (int i = 0; i < hfake1->GetNbinsX(); i++) {
                out << "pT [ " << int(hfake1->GetBinLowEdge(i+1)) << " - " << int(hfake1->GetBinLowEdge(i+1))+int(hfake1->GetBinWidth(i+1)) << "] : " << eebarrel[i]  << " +/- " << eebarrelerr[i]  << endl;
            }
            out << std::endl;
            out << "---- Electron Electron Endcap ----" << std::endl;
            for (int i = 0; i < hfake1->GetNbinsX(); i++) {
                out << "pT [ " << int(hfake1->GetBinLowEdge(i+1)) << " - " << int(hfake1->GetBinLowEdge(i+1))+int(hfake1->GetBinWidth(i+1)) << "] : " << eeendcap[i]  << " +/- " << eeendcaperr[i]  << endl;
            }
            out << std::endl;

            out << "---- Muon ElectronBarrel ----" << std::endl;
            for (int i = 0; i < hfake3->GetNbinsX(); i++) {
                out << "pT [ " << int(hfake3->GetBinLowEdge(i+1)) << " - " << int(hfake3->GetBinLowEdge(i+1))+int(hfake3->GetBinWidth(i+1)) << "] : " << mebarrel[i]  << " +/- " << mebarrelerr[i]  << endl;
            }
            out << std::endl;
            out << "---- Muon Electron Endcap ----" << std::endl;
            for (int i = 0; i < hfake3->GetNbinsX(); i++) {
                out << "pT [ " << int(hfake3->GetBinLowEdge(i+1)) << " - " << int(hfake3->GetBinLowEdge(i+1))+int(hfake3->GetBinWidth(i+1)) << "] : " << meendcap[i]  << " +/- " << meendcaperr[i]  << endl;
            }
            out << std::endl;

        }


        TH2D* getFakeRateHist(int ch, int mde) {
            if (mde == 0) {
                if (ch == 0 || ch == 2) return hfake0;
                else return hfake1;
            }

            else {
                if      (ch == 0) return hfake0;
                else if (ch == 1) return hfake1;
                else if (ch == 2) return hfake2;
                else              return hfake3;
            }
        } 
};

#endif 
