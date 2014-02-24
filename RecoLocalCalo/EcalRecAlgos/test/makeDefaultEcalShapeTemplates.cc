#include "SimCalorimetry/EcalSimAlgos/interface/EcalSimParameterMap.h"
#include "SimCalorimetry/EcalSimAlgos/interface/APDShape.h"
#include "SimCalorimetry/EcalSimAlgos/interface/EBShape.h"
#include "SimCalorimetry/EcalSimAlgos/interface/EEShape.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include<iostream>
#include<iomanip>

#include "TROOT.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TF1.h"

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

      TH1F* shape = new TH1F( "shape", title.c_str(), nsamp, 0., (float) nsamp / 25. ) ;
      double y = 0.;
      double dy = 0.;

      for( unsigned int i ( 0 ) ; i != histsiz ; ++i ) 
        {
          y  = (*theShape)(x);
          shape->Fill((float)(x-tzero)/25.,(float)y);
          std::cout << " time (ns) = "  << std::fixed    << std::setw(6)         << std::setprecision(2) << x-tzero 
                    << " shape = "      << std::setw(11) << std::setprecision(5) << y << std::endl;
          x = x+1./(double)tconv;
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

      delete shape;
      delete showShape;

    }

  return 0;
} 
