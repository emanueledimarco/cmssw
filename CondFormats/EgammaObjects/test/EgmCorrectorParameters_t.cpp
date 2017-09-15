#include "Utilities/Testing/interface/CppUnit_testdriver.icpp"
#include "cppunit/extensions/HelperMacros.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "CondFormats/EgammaObjects/interface/EgmCorrectorParameters.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cstdlib>

#include "TRandom3.h"
#include "TBenchmark.h"

using std::cout;
using std::endl;
using std::flush;
using std::vector;
using std::string;
using std::setw;
using std::ofstream;

class testEgmCorrectorParameters: public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(testEgmCorrectorParameters);
  CPPUNIT_TEST(setUp);
  CPPUNIT_TEST(setupCorrector);
  CPPUNIT_TEST_SUITE_END();

public:
  testEgmCorrectorParameters() {}
  ~testEgmCorrectorParameters() {}
  void setUp();
  void tearDown();
  void setupCorrector();
  void destroyCorrector();
  void generateFiles();

  inline void loadbar3(unsigned int x, unsigned int n, unsigned int w = 50, unsigned int freq = 100, string prefix = "") {
    if ( (x != n) && (x % (n/freq) != 0) ) return;
    float ratio  =  x/(float)n;
    int   c      =  ratio * w;

    cout << prefix << std::fixed << setw(8) << std::setprecision(0) << (ratio*100) << "% [";
    for (int x=0; x<c; x++) cout << "=";
    for (unsigned int x=c; x<w; x++) cout << " ";
    cout << "] (" << x << "/" << n << ")\r" << flush;
  }

private:
  string filename;
  TBenchmark* m_benchmark;
  EgmCorrectorParameters* electronScalePar;
  vector<EgmCorrectorParameters> vPar;
  vector<float> fX;
  vector<float> veta;
  vector<float> vr9;
  vector<float> vpt;
};

  ///registration of the test so that the runner can find it
  CPPUNIT_TEST_SUITE_REGISTRATION(testEgmCorrectorParameters);

void testEgmCorrectorParameters::setUp() {
  m_benchmark = new TBenchmark();
  veta = {-2.5, -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, 
          -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 
          0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 
          1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5};
  for(unsigned int ir9=0; ir9<2; ir9++) {
    vr9.push_back(ir9);
  }
  vpt = {5,6,7,8,9,10,11,12,13,14,15,17,20,23,27,30,35,40,45,57,72,90,120,
         150,200,300,400,550,750,1000,6500};

  generateFiles();
}

void testEgmCorrectorParameters::tearDown() {
  if(m_benchmark!=nullptr)
    delete m_benchmark;
}

void testEgmCorrectorParameters::setupCorrector() {
  string path = "CondFormats/EgammaObjects/data/";
  try {
    string to_try = path + ("testEgmCorrectorParameters_EtaR9Et.txt");
    edm::FileInPath strFIP(to_try);
    filename = strFIP.fullPath();
    std::cout << "Setting up corrector with the file " << filename << std::endl;
  }
  catch (edm::Exception ex) {
    throw ex;
  }

  std::cout << "Loading... " << std::endl;
  electronScalePar = new EgmCorrectorParameters(filename);
  vPar.push_back(*electronScalePar);
}

void testEgmCorrectorParameters::destroyCorrector() {
  delete electronScalePar;
}

void testEgmCorrectorParameters::generateFiles() {
  string path = std::getenv("CMSSW_BASE");
  path+="/src/CondFormats/EgammaObjects/data/";
  string name3D = "testEgmCorrectorParameters_EtaR9Et.txt";

  std::cout << "Writing dummy corrections file in " << name3D << std::endl;
  ofstream txtFile3D;
  txtFile3D.open((path+name3D).c_str());
  CPPUNIT_ASSERT(txtFile3D.is_open());
  txtFile3D << "{3 eta r9 ET 3 Correction StatUnc SystUnc}" << endl;
  float statUnc = 1e-04;
  float systUnc = 1e-03;
  for(unsigned int ieta=0; ieta<veta.size()-1; ieta++) {
    for(unsigned int ir9=0; ir9<vr9.size()-1; ir9++) {
      for(unsigned int ipt=0; ipt<vpt.size()-1; ipt++) {
        txtFile3D << std::setw(15) << veta[ieta] << std::setw(15) << veta[ieta+1]
                  << std::setw(15) << vr9[ir9] << std::setw(15) << vr9[ir9+1];
        txtFile3D << std::setw(15) << vpt[ipt] << std::setw(15);
        if(ipt+1==vpt.size()-1) txtFile3D << "6500";
        else                       txtFile3D << vpt[ipt+1];
        txtFile3D << std::setw(15) << gRandom->Gaus(1.0,1e-03) 
                  << std::setw(15) << statUnc << std::setw(15) << systUnc << std::endl;
      }
    }
  }
  txtFile3D.close();
  std::cout << "Done." << std::endl;
}

