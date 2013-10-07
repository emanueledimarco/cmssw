#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPad.h>
#include <TPaveText.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <algorithm>
#include <TF1.h>
#include <../FitMaker.h>

#define YEAR 2011

enum {
  kMuRef = 0,
  kMuMvaIso
};

enum {
  kEleRef = 0,
  kEleMvaLoose,
  kEleMvaTight
};

enum {
  kFF = 0,
  kMatrix
};

#define METHOD kMatrix

float getFRbin(float pT, float eta, TH2* myH, int nsigma) {
  float fr=-1.;
  int   xBins = myH->GetXaxis()->GetNbins();
  float xMin  = myH->GetXaxis()->GetBinLowEdge(1);
  float xMax  = myH->GetXaxis()->GetBinUpEdge(xBins);
  int   yBins = myH->GetYaxis()->GetNbins();
  float yMin  = myH->GetYaxis()->GetBinLowEdge(1);
  float yMax  = myH->GetYaxis()->GetBinUpEdge(yBins);
  int theBin = myH->FindBin(pT, fabs(eta));
  if (pT>xMin && pT<xMax && fabs(eta)>yMin && fabs(eta)<yMax) {
    fr = myH->GetBinContent(theBin);
    fr += nsigma * myH->GetBinError(theBin);
  } else if(pT>xMax && fabs(eta)<yMax) { 
    int theOvf = myH->FindBin(xMax-1.0,fabs(eta));
    fr = myH->GetBinContent(theOvf);
    fr += nsigma * myH->GetBinError(theOvf);
  } else if(pT<xMax && fabs(eta)>yMax) { 
    int theOvf = myH->FindBin(pT,fabs(eta)-0.1);
    fr = myH->GetBinContent(theOvf);
    fr += nsigma * myH->GetBinError(theOvf);
  } else if(pT>xMax && fabs(eta)>yMax) {
    int theOvf = myH->FindBin(xMax-0.1,fabs(eta)-0.1);
    fr = myH->GetBinContent(theOvf);
    fr += nsigma * myH->GetBinError(theOvf);
  } else cout << "fuck you" << endl;

  return fr;
}

float computeSF(float pt, float abseta, bool isMu, int elwp, int muwp, 
		TH2 *elfake, TH2 *mufake, int nsigma, bool usefunction=false) {

  if(!usefunction) {
    if(isMu) return getFRbin(pt, abseta, mufake, nsigma);
    else return getFRbin(pt, abseta, elfake, nsigma);
  } else {
  // 2012 fake rates parameterizations (PromptReco)
  float mu2012_bar_par0[2] = { 0.088146, 0.172263};
  float mu2012_bar_par1[2] = { 0.813482, 0.267368};
  float mu2012_bar_par2[2] = { 0.100546, 0.099457};
  float mu2012_bar_par3[2] = { 0.078343, 0.026675};
    
  float mu2012_end_par0[2] = { 0.097366, 0.233072 };
  float mu2012_end_par1[2] = { 0.829411, 0.332927 };
  float mu2012_end_par2[2] = { 0.10532,  0.126052 };
  float mu2012_end_par3[2] = { 0.081078, 0.020154 };
    
  float el2012_bar_par0[3] = { 0.0545295, 0.0736869, 0.0560685 };
  float el2012_bar_par1[3] = { 1.04825,   1.17273,   0.77083   };
  float el2012_bar_par2[3] = { 0.214458,  0.18919,   0.149952  };
  float el2012_bar_par3[3] = { 0.194721,  0.165607,  0.131008  };
  
  float el2012_end_par0[3] = { 0.172428,  0.192236,  0.18395  };
  float el2012_end_par1[3] = { 2.59677,   1.6312,    0.592073 };
  float el2012_end_par2[3] = { 0.26902,   0.252914,  0.236726 };
  float el2012_end_par3[3] = { 0.217469,  0.175229,  0.100146 };

  // 2011 fake rates parameterizations (Jan16 rereco)
  float mu2011_bar_par0[2] = { 0.0442711, 0.0924824 };
  float mu2011_bar_par1[2] = { 1.6023,    0.0 };
  float mu2011_bar_par2[2] = { 0.398586,  0.101534 };
  float mu2011_bar_par3[2] = { 0.386322,  0.0943444 };
  
  float mu2011_end_par0[2] = { 0.0364204, 0.13092 };
  float mu2011_end_par1[2] = { 0.562334,  1.40327 };
  float mu2011_end_par2[2] = { 0.171554,  0.0345104 };
  float mu2011_end_par3[2] = { 0.167615,  0.0291992 };
  
  float el2011_bar_par0[3] = { 0.0786844, 0.0868894, 0.072289 };
  float el2011_bar_par1[3] = { 1.84035,   1.3331,    1.74931  };
  float el2011_bar_par2[3] = { 0.251348,  0.25139,   0.208941 };
  float el2011_bar_par3[3] = { 0.224511,  0.211936,  0.189535 };
 
  float el2011_end_par0[3] = { 0.207228,  0.221917, 0.173104 };
  float el2011_end_par1[3] = { 3.79389,   0.75421,  1.95527  };
  float el2011_end_par2[3] = { 0.326206,  1.76081,  0.229775 };
  float el2011_end_par3[3] = { 0.261076,  0.155266, 0.176875 };

  if (YEAR == 2012) {
    if (!isMu && abseta<1.479)   return el2012_bar_par0[elwp] + el2012_bar_par1[elwp]*(exp(-el2012_bar_par2[elwp]*pt) - exp(-el2012_bar_par3[elwp]*pt));
    if (!isMu && abseta>1.479)   return el2012_end_par0[elwp] + el2012_end_par1[elwp]*(exp(-el2012_end_par2[elwp]*pt) - exp(-el2012_end_par3[elwp]*pt));
    
    if (isMu && abseta<1.479)    return mu2012_bar_par0[muwp] + mu2012_bar_par1[muwp]*(exp(-mu2012_bar_par2[muwp]*pt) - exp(-mu2012_bar_par3[muwp]*pt));
    if (isMu && abseta>1.479)    return mu2012_end_par0[muwp] + mu2012_end_par1[muwp]*(exp(-mu2012_end_par2[muwp]*pt) - exp(-mu2012_end_par3[muwp]*pt));
  } else { 
    if (!isMu && abseta<1.479)   return el2011_bar_par0[elwp] + el2011_bar_par1[elwp]*(exp(-el2011_bar_par2[elwp]*pt) - exp(-el2011_bar_par3[elwp]*pt));
    if (!isMu && abseta>1.479)   return el2011_end_par0[elwp] + el2011_end_par1[elwp]*(exp(-el2011_end_par2[elwp]*pt) - exp(-el2011_end_par3[elwp]*pt));
    
    if (isMu && abseta<1.479)    return mu2011_bar_par0[muwp] + mu2011_bar_par1[muwp]*(exp(-mu2011_bar_par2[muwp]*pt) - exp(-mu2011_bar_par3[muwp]*pt));
    if (isMu && abseta>1.479)    return mu2011_end_par0[muwp] + mu2011_end_par1[muwp]*(exp(-mu2011_end_par2[muwp]*pt) - exp(-mu2011_end_par3[muwp]*pt));
  }
  }
  return 0;
}

float dozx(int ch, int elwp, int muwp, bool barecount) {

    TFile *elfakefile, *mufakefile;
    elfakefile = mufakefile = 0;
    if(elwp==kEleRef) {
      elfakefile = TFile::Open("src/elfrmapPfIso.root");
    } else if(elwp==kEleMvaLoose) {
      elfakefile = TFile::Open("src/elfrmapMvaLoose.root");
    } else if(elwp==kEleMvaTight) {
      elfakefile = TFile::Open("src/elfrmapMvaTight.root");
    } else {
      cout << "WP not allowed for ele." << endl;
    }
    if(muwp==kMuRef) {
      mufakefile = TFile::Open("src/mufrmapPfIso.root");
    } else if(muwp==kMuMvaIso) {
      mufakefile = TFile::Open("src/mufrmapMvaPfIso.root");
    }  else {
      cout << "WP not allowed for mu." << endl;
    }
    TH2D *elfake = (TH2D*)elfakefile->Get("hfake");
    TH2D *mufake = (TH2D*)mufakefile->Get("hfake");

    // prompt rate (using efficiency in 1st approximation)
    TFile *elprfile = TFile::Open("src/hel_eff_2011.root");
    TFile *muprfile = TFile::Open("src/scale_factors_withTH2.30-05-12.root");
    TH2F *elpr = (TH2F*)elprfile->Get("heff");
    TH2F *mupr = (TH2F*)muprfile->Get("TH2D_ALL_2012");

    TChain* tree = new TChain("zxTree/probe_tree");
    tree->Add("results_data/hzzTree_data2011.root");

    TBranch *bl3pt      = tree->GetBranch("l3pt");
    TBranch *bl3eta     = tree->GetBranch("l3eta");
    TBranch *bl4pt      = tree->GetBranch("l4pt");
    TBranch *bl4eta     = tree->GetBranch("l4eta");

    TBranch *bmass      = tree->GetBranch("mass");
    TBranch *bz1mass    = tree->GetBranch("z1mass");
    TBranch *bz2mass    = tree->GetBranch("z2mass");
    TBranch *bchannel   = tree->GetBranch("channel");
    TBranch *bl3pdgId   = tree->GetBranch("l3pdgId");
    TBranch *bl4pdgId   = tree->GetBranch("l4pdgId");
    TBranch *bl3idNew   = tree->GetBranch("l3idNew");
    TBranch *bl4idNew   = tree->GetBranch("l4idNew");

    TBranch *bl3pfIsoComb04EACorr = tree->GetBranch("l3pfIsoComb04EACorr");
    TBranch *bl4pfIsoComb04EACorr = tree->GetBranch("l4pfIsoComb04EACorr");

    TBranch *bevent     = tree->GetBranch("event");
    TBranch *brun       = tree->GetBranch("run");

    float l3pt      = 0.0;
    float l3eta     = 0.0;
    float l4pt      = 0.0;
    float l4eta     = 0.0;

    float mass      = 0.0;
    float z1mass    = 0.0;
    float z2mass    = 0.0;
    float channel   = 0.0;
    float l3pdgId   = 0.0;
    float l4pdgId   = 0.0;
    float l3idNew   = 0.0;
    float l4idNew   = 0.0;
    float l3pfIsoComb04EACorr = 0.0;
    float l4pfIsoComb04EACorr = 0.0;

    long   event     = 0.0;
    long   run       = 0.0;

    bl3pt   ->SetAddress(&l3pt);
    bl3eta  ->SetAddress(&l3eta);
    bl4pt   ->SetAddress(&l4pt);
    bl4eta  ->SetAddress(&l4eta);

    bmass   ->SetAddress(&mass);
    bz1mass ->SetAddress(&z1mass);
    bz2mass ->SetAddress(&z2mass);
    bchannel->SetAddress(&channel);
    bl3pdgId ->SetAddress(&l3pdgId);
    bl4pdgId ->SetAddress(&l4pdgId);
    bl3idNew ->SetAddress(&l3idNew);
    bl4idNew ->SetAddress(&l4idNew);

    bl3pfIsoComb04EACorr ->SetAddress(&l3pfIsoComb04EACorr);
    bl4pfIsoComb04EACorr ->SetAddress(&l4pfIsoComb04EACorr);

    bevent  ->SetAddress(&event);
    brun    ->SetAddress(&run);

    ZXFitMaker zxfit("xzfit",120,100,600);

    float yield = 0.0;
    std::vector<int> usedevents;
    for (int i = 0; i < tree->GetEntries(); i++) {

      bl3pt   ->GetEvent(i);
      bl3eta  ->GetEvent(i);
      bl4pt   ->GetEvent(i);
      bl4eta  ->GetEvent(i);

      bmass   ->GetEvent(i);
      bz1mass ->GetEvent(i);
      bz2mass ->GetEvent(i);
      bchannel->GetEvent(i);
      bl3pdgId ->GetEvent(i);
      bl4pdgId ->GetEvent(i);
      bl3idNew ->GetEvent(i);
      bl4idNew ->GetEvent(i);

      bl3pfIsoComb04EACorr ->GetEvent(i);
      bl4pfIsoComb04EACorr ->GetEvent(i);
      
      bevent  ->GetEvent(i);
      brun    ->GetEvent(i);

      vector<int>::iterator it;
      it = find (usedevents.begin(), usedevents.end(), event);
      if(it!=usedevents.end()) continue;

      usedevents.push_back(event);

      float z1min = 40;
      float z2min = 12;
      
      // float z1min = 50;
      // float z2min = 12;
	
        if (z1mass<z1min || z2mass<z2min) continue;
        if (channel != ch) continue;
        if (mass<100 || mass>600) continue;
      	//if (l3pdgId!=l4pdgId) continue; // using SS leptons
	if(METHOD==kFF) {
	  bool osfail = l3pdgId==(-1*l4pdgId) && (l3idNew==0 || l3pfIsoComb04EACorr/l3pt>0.40) && (l4idNew==0 || l4pfIsoComb04EACorr/l4pt>0.40);
	  if (!osfail) continue;
	  //if (z1mass>120 || z2mass>120) continue;
	  
	  float sf1 = 0.0; 
	  float sf2 = 0.0; 
	  
	  bool isMu = false;
	  if (channel == 0 || channel == 2) isMu = true;
	  
	  float osss = 0.0;
	  // if (channel == 0) osss = 1.28;    
	  // if (channel == 1) osss = 0.93;    
	  // if (channel == 2 || channel == 3) osss = 0.94;    
	  if (channel == 0) osss = 1.0;    
	  if (channel == 1) osss = 1.0;    
	  if (channel == 2 || channel == 3) osss = 1.0;    
	  
	  sf1 = computeSF(l3pt, fabs(l3eta), isMu, elwp, muwp, elfake, mufake, 0);
	  sf2 = computeSF(l4pt, fabs(l4eta), isMu, elwp, muwp, elfake, mufake, 0);
	  
	  //yield += (sf11*sf21 + 0.5*sf11*sf22 + 0.5*sf21*sf12)*osss;  // PRL triangle
	  //yield += sf1*sf2*osss;  // ss
	  if(barecount) yield += 1.0;
	  else yield += sf1/(1-sf1)*sf2/(1-sf2); // osfail
	}
	if(METHOD==kMatrix) {

	  float sf1 = 0.0; 
	  float sf2 = 0.0; 
	  float p1 = 0.0;
	  float p2 = 0.0;

	  bool isMu = false;
	  if (channel == 0 || channel == 2) isMu = true;

	  sf1 = computeSF(l3pt, fabs(l3eta), isMu, elwp, muwp, elfake, mufake, 0);
	  sf2 = computeSF(l4pt, fabs(l4eta), isMu, elwp, muwp, elfake, mufake, 0);
	  p1 = computeSF(l3pt, fabs(l3eta), isMu, elwp, muwp, elpr, mupr, 0);
	  p2 = computeSF(l4pt, fabs(l4eta), isMu, elwp, muwp, elpr, mupr, 0);

	  float eps1=sf1/(1-sf1);
	  float eps2=sf2/(1-sf2);
	  // float eta1=(1-p1)/p1;
	  // float eta2=(1-p2)/p2;
	  float eta1, eta2;
	  eta1=eta2=0.;
	  
	  float denom=1/(1-eps1*eta1)/(1-eps2*eta2);
	  
	  bool l3t=l3idNew==1 && l3pfIsoComb04EACorr/l3pt<0.40;
	  bool l4t=l4idNew==1 && l4pfIsoComb04EACorr/l4pt<0.40;

	  float wfp,wff,wpp;
	  wfp=wff=wpp=0.0;

	  // t00
	  if( (l3pdgId==-l4pdgId) && (!l3t) && (!l4t) ) {
	    wfp=-2*eps1*eps2 * denom;
	    wff=eps1*eps2 * denom;
	    wpp=eps1*eps2 * denom;
	    if(barecount) yield += 1.0;	  
	  }
	  // t10
	  if( (l3pdgId==-l4pdgId) && l3t && (!l4t) ) {
	    wfp=(eps2+eps2*eps1*eta2) * denom;
	    wff=-1*eps1*eps2*eta2 * denom;
	    wpp=(-eps1) * denom;
	    if(barecount) yield += 1.0;	  
	  }
	  // t01
	  if( (l3pdgId==-l4pdgId) && (!l3t) && l4t ) {
	    wfp = eps2*denom;
	    wfp=(eps1*eps2*eta1 + eps1) * denom;
	    wff=-1*eps1*eps2*eta1 * denom;
	    wpp=(-eps2) * denom;
	    if(barecount) yield += 1.0;	  
	  }
	  // t11
	  if( (l3pdgId==-l4pdgId) && l3t && l4t ) {
	    // wfp =(-eps1*eta1-eps2*eta2) * denom;
	    // wff = eps1*eps2*eta1*eta2 * denom;
	    // wpp = denom;
	    //	    if(barecount) yield += 1.0;	  
	  }
	  float weight = wfp+wff;
	  zxfit.add(mass,weight);
	  if(!barecount) yield += weight;
	}
    }

    return yield;
    delete tree;
}

void zx(int elwp, int muwp) {

    float yield_ee = 0.0;
    float yield_mm = 0.0;
    float yield_em = 0.0;

    float crcount_ee = 0.0;
    float crcount_mm = 0.0;
    float crcount_em = 0.0;

    yield_ee += dozx(1, elwp, muwp,false);
    yield_mm += dozx(0, elwp, muwp,false);
    yield_em += dozx(2, elwp, muwp,false);
    yield_em += dozx(3, elwp, muwp,false);

    crcount_ee += dozx(1, elwp, muwp,true);
    crcount_mm += dozx(0, elwp, muwp,true);
    crcount_em += dozx(2, elwp, muwp,true);
    crcount_em += dozx(3, elwp, muwp,true);

    cout << "yield(4e)   = " << yield_ee << " (CR = " << crcount_ee << " events)" << std::endl;
    cout << "yield(4mu)  = " << yield_mm << " (CR = " << crcount_mm << " events)" << std::endl;
    cout << "yeld(2e2mu) = " << yield_em << " (CR = " << crcount_em << " events)" << std::endl;
}


void all() {
  std::cout << "--- elpfiso / mupfiso ---" << std::endl;
  zx(kEleRef,kMuRef);
  std::cout << "-------------------------" << std::endl;

  std::cout << "--- elpfiso / muMvapfisoLoose ---" << std::endl;
  zx(kEleRef,kMuMvaIso);
  std::cout << "-------------------------" << std::endl;

  std::cout << "--- elpfisoLoose / mupfiso ---" << std::endl;
  zx(kEleMvaLoose,kMuRef);
  std::cout << "-------------------------" << std::endl;

  std::cout << "--- elpfisoLoose / muMvapfisoLoose ---" << std::endl;
  zx(kEleMvaLoose,kMuMvaIso);
  std::cout << "-------------------------" << std::endl;

  std::cout << "--- elpfisoTight / mupfiso ---" << std::endl;
  zx(kEleMvaTight,kMuRef);
  std::cout << "-------------------------" << std::endl;

  std::cout << "--- elpfisoTight / muMvapfisoLoose ---" << std::endl;
  zx(kEleMvaTight,kMuMvaIso);
  std::cout << "-------------------------" << std::endl;


}
