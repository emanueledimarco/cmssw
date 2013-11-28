#include "HiggsAnalysis/HiggsToWW2Leptons/plugins/HWWUtils.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include <iostream>

void setMomentum (TLorentzVector & myvector, const reco::Candidate & gen)
{
    myvector.SetPx (gen.px ()) ;
    myvector.SetPy (gen.py ()) ;
    myvector.SetPz (gen.pz ()) ;
    myvector.SetE (gen.energy ()) ;
}


// --------------------------------------------------------------------


std::pair<HWWjetIt,HWWjetIt>	
findTagJets (HWWjetIt begin, HWWjetIt end,
             double jetPtMin, double jetEtaMax) 
{
	
  std::pair<HWWjetIt,HWWjetIt> tagJets (begin,begin) ;
  double maxInvMass = 0. ;

  //PG find the tagging jets

  //PG first loop over jets
  for (HWWjetIt firstJet = begin ; 
		 firstJet != end ; 
		 ++firstJet ) 
  {
    if (firstJet->pt () < jetPtMin || 
            fabs (firstJet->eta ()) > jetEtaMax) continue ;
    math::XYZTLorentzVector firstLV = firstJet->p4 () ;
    //PG second loop over jets
    for (HWWjetIt secondJet = firstJet + 1 ; 
        	       secondJet != end ; 
        	       ++secondJet ) {
      if (secondJet->pt () < jetPtMin || 
              fabs (secondJet->eta ()) > jetEtaMax) continue ;
      math::XYZTLorentzVector sumLV = secondJet->p4 () + firstLV ;
      if (sumLV.M () > maxInvMass) {
        		    maxInvMass = sumLV.M () ;
        		    tagJets.first = firstJet ;
        		    tagJets.second = secondJet ;
      }
    } //PG second loop over jets
  } //PG first loop over jets

  return tagJets ;
}

double deltaPhi (double phi1, double phi2)
{
  double deltaphi=fabs(phi1-phi2);
  if (deltaphi > 2*M_PI) deltaphi -= 2*M_PI;
  if (deltaphi > M_PI)   deltaphi  = 2*M_PI-deltaphi;
  return deltaphi;
}
