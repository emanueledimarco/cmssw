#include <string>
#include <sstream>
#include <vector>
#include <iostream>

#include "TFile.h"
#include "TProfile2D.h"

using namespace std;

void templatesFromLaserAllCrystals(const char *dqmfile) {

  TFile *tfile = TFile::Open(dqmfile);
  string dirEB = "/DQMData/Run 202299/EcalBarrel/Run summary/EBLaserTask/Laser1";
  string dirEE = "/DQMData/Run 202299/EcalEndcap/Run summary/EELaserTask/Laser1";

  tfile->cd(dirEB.c_str());
  vector<TProfile2D> templates;
  for(int ism=1; ism<=36; ++ism) {
    stringstream ebname;
    if(ism<10) ebname << "EBLT shape EB+0" << ism << " L1";
    else if(ism<19) ebname << "EBLT shape EB+" << ism << " L1";
    else if(ism<28) ebname << "EBLT shape EB-0" << ism-18 << " L1";
    else ebname << "EBLT shape EB-" << ism-18 << " L1";
    cout << "ism = " << ism << "  SM name = " << ebname.str() << endl;
    TProfile2D *shape = (TProfile2D*)gDirectory->Get(ebname.str().c_str());
    templates.push_back(*shape);
  }

  tfile->cd(dirEE.c_str());
  for(int ism=-9; ism<=9; ++ism) {
    if(ism==0) continue;
    stringstream eename;
    if(ism<0) eename << "EELT shape EE-0" << abs(ism) << " L1";
    else eename << "EELT shape EE+0" << ism << " L1";
    cout << "ism = " << ism << "  SM name = " << eename.str() << endl;
    TProfile2D *shape = (TProfile2D*)gDirectory->Get(eename.str().c_str());
    templates.push_back(*shape);
  }

  TFile *strippedfile = TFile::Open("templates.root","recreate");
  strippedfile->cd();
  
  for(int t=0; t<(int)templates.size(); ++t) templates[t].Write();
  
  strippedfile->Close();

}

void templatesFromLaserSingleCrystal(const char *dqmfile) {

  TFile *tfile = TFile::Open(dqmfile);
  string dirEB = "/DQMData/Run 202299/EcalBarrel/Run summary/EBLaserClient";
  string dirEE = "/DQMData/Run 202299/EcalEndcap/Run summary/EELaserClient";

  tfile->cd(dirEB.c_str());
  vector<TH1F> templates;
  for(int ism=1; ism<=36; ++ism) {
    stringstream ebname;
    if(ism<10) ebname << "EBLT laser shape L1 EB+0" << ism;
    else if(ism<19) ebname << "EBLT laser shape L1 EB+" << ism;
    else if(ism<28) ebname << "EBLT laser shape L1 EB-0" << ism-18;
    else ebname << "EBLT laser shape L1 EB-" << ism-18;
    cout << "ism = " << ism << "  SM name = " << ebname.str() << endl;
    TH1F *shape = (TH1F*)gDirectory->Get(ebname.str().c_str());
    string newname = ebname.str().replace(5,5,"single");
    shape->SetName(newname.c_str());
    templates.push_back(*shape);
  }

  tfile->cd(dirEE.c_str());
  for(int ism=-9; ism<=9; ++ism) {
    if(ism==0) continue;
    stringstream eename;
    if(ism<0) eename << "EELT laser shape L1 EE-0" << abs(ism);
    else eename << "EELT laser shape L1 EE+0" << ism;
    cout << "ism = " << ism << "  SM name = " << eename.str() << endl;
    TH1F *shape = (TH1F*)gDirectory->Get(eename.str().c_str());
    string newname = eename.str().replace(5,5,"single");
    shape->SetName(newname.c_str());
    templates.push_back(*shape);
  }

  TFile *strippedfile = TFile::Open("templates.root","update");
  strippedfile->cd();
  
  for(int t=0; t<(int)templates.size(); ++t) templates[t].Write();
  
  strippedfile->Close();

}


void makeTemplates(const char *dqmfile="~/Work/data/ecalreco/DQM_V0013_EcalBarrel_R000202299.root") {
  templatesFromLaserAllCrystals(dqmfile);
  templatesFromLaserSingleCrystal(dqmfile);
}
