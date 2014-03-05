#include "SimCalorimetry/EcalSimAlgos/interface/EcalSimParameterMap.h"
#include "SimCalorimetry/EcalSimAlgos/interface/APDShape.h"
#include "SimCalorimetry/EcalSimAlgos/interface/EBShape.h"
#include "SimCalorimetry/EcalSimAlgos/interface/EEShape.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <iostream>
#include <iomanip>
#include <vector>

#include "TROOT.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TMath.h"
#include "TFile.h"

int main ()
{  
  edm::MessageDrop::instance()->debugEnabled = false;

  const EcalSimParameterMap parameterMap ;

  const APDShape theAPDShape( 74.5, 40.5 ) ;
  const EBShape theEBShape ;
  const EEShape theEEShape ;

  const int nsamp = EcalShapeBase::k1NSecBinsTotal;
  const int tconv = EcalShapeBase::kNBinsPerNSec;

  const unsigned int histsiz = nsamp*tconv;

  TFile *fileShape = TFile::Open("EcalShapes.root","recreate");

  for( unsigned int i ( 0 ) ; i != 3 ; ++i )
    {  

      const DetId id ( 0 == i || 2 == i ?
                       (DetId) EBDetId(1,1) :
                       (DetId) EEDetId(1,50,1) ) ;

      const EcalShapeBase* theShape ( 0 == i ?
                                      (EcalShapeBase*) &theEBShape :
                                      ( 1 == i ?
                                        (EcalShapeBase*) &theEEShape :
                                        (EcalShapeBase*) &theAPDShape  ) ) ;

      const double ToM ( theShape->timeOfMax()  ) ;
      const double T0  ( theShape->timeOfThr()  ) ;
      const double rT  ( theShape->timeToRise() ) ;


      // standard display of the implemented shape function
      const int csize = EcalShapeBase::k1NSecBinsTotal;
      TCanvas * showShape = new TCanvas("showShape","showShape",csize,csize);

      const std::string name ( 0 == i ? "Barrel" : 
                               ( 1 == i ? "Endcap" :
                                 "APD" ) ) ;

      std::cout << "\n ********************* "
                << name 
                << "************************" ; 

      std::cout << "\n Maximum time from tabulated values = " 
                << std::fixed    << std::setw(6)   
                << std::setprecision(2) << ToM << std::endl ;

      std::cout << "\n Tzero from tabulated values        = " 
                << std::fixed    << std::setw(6)   
                << std::setprecision(2) << T0 << std::endl ;

      std::cout << "\n Rising time from tabulated values  = " 
                << std::fixed    << std::setw(6)   
                << std::setprecision(2) << rT << std::endl;

      // signal used with the nominal parameters and no jitter

      std::cout << "\n computed ECAL " << name 
                << " pulse shape (LHC timePhaseShift = 1) \n" << std::endl;
      const double tzero = rT - ( parameterMap.simParameters( id ).binOfMaximum() - 1. )*25. ;
      double x = tzero ;

      const std::string title ( "Computed Ecal " + name + " MGPA shape" ) ;
      const std::string nameh  ( "Ecal" + name + "Shape" ) ;

      double maxT = ((double)nsamp) / 25.;
      TH1F* shape = new TH1F( nameh.c_str(), title.c_str(), nsamp, 0., maxT ) ;
      double y = 0.;

      for( unsigned int i ( 0 ) ; i != histsiz ; ++i ) 
        {
          y  = (*theShape)(x);
          double scaledX = ((double)(x-tzero))/25.;
          shape->Fill(scaledX,(double)y);
          std::cout << " time (ns) = "  << std::fixed    << std::setw(6)         << std::setprecision(2) << x-tzero 
                    << " shape = "      << std::setw(11) << std::setprecision(5) << y << std::endl;
          x = x+1./(double)tconv;
        }
      
      // remove the rounding spike
      for (int i = 2 ; i< shape->GetNbinsX(); ++i) {
        if(shape->GetBinContent(i)>shape->GetBinContent(i-1)+0.5) shape->SetBinContent(i,(shape->GetBinContent(i-1)+shape->GetBinContent(i+1))/2.);
      }

      std::vector<TH1F*> oot_shapes;
      // replicate it for -10bx +3bx
      for (int i = -10; i<4; ++i) {
        if(i==0) continue;
        char buf[100],titleo[100];
        char sign[1];
        if(i>0) sprintf(sign,"+"); else sprintf(sign,"");
        sprintf(buf,"%s_oot%s%dbx",nameh.c_str(),sign,i);
        sprintf(titleo,"%s OOT %s%d BX",title.c_str(),sign,i);
        TH1F *ootshape = new TH1F(buf, titleo, nsamp, i, maxT+i);
        for(int b=0; b<nsamp; ++b) ootshape->SetBinContent(b,shape->GetBinContent(b));
        oot_shapes.push_back(ootshape);
      }

      for( unsigned int iSample ( 0 ) ; iSample != 10 ; ++iSample ) 
        {
          std::cout << (*theShape)(tzero + iSample*25.0) << std::endl; 
        }

      showShape->cd();
      gPad->SetGrid();
      shape->GetXaxis()->SetNdivisions(10,kFALSE);
      shape->Draw();

      const std::string fname ( name + "EcalShapeUsed.pdf" ) ;
      showShape->SaveAs( fname.c_str() );

      fileShape->cd();
      shape->Write();
      for(int i=0; i<(int)oot_shapes.size(); ++i) { 
        oot_shapes[i]->Write();
        delete oot_shapes[i];
      }
        
      delete shape;
      delete showShape;

    }

  fileShape->Close();

  return 0;
} 
