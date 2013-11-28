#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "HiggsAnalysis/HiggsToWW2Leptons/plugins/EleTrackerIsolationAlgo.h"
#include "CLHEP/Vector/ThreeVector.h"

using namespace std;
EleTrackerIsolationAlgo::EleTrackerIsolationAlgo (){}

EleTrackerIsolationAlgo::EleTrackerIsolationAlgo (const GsfElectron *gsfEle, const TrackCollection trackColl) : 
  _myGsfEle(gsfEle),
  _tracks(trackColl)   
{
  _extRadius = 0.200;
  _intRadius = 0.015;
}

EleTrackerIsolationAlgo::~EleTrackerIsolationAlgo (){}


void EleTrackerIsolationAlgo::setExtRadius (float extRadius){_extRadius = extRadius; }

void EleTrackerIsolationAlgo::setIntRadius (float intRadius){_intRadius = intRadius; }

// sum of pt of tracks (within a given cone) / electron pt
float EleTrackerIsolationAlgo::getPtTracks () const
{
  float dummyPt = 0 ;

  CLHEP::Hep3Vector elePAtVtx(_myGsfEle->px(), _myGsfEle->py(), _myGsfEle->pz()); 
  float ele_pt  = _myGsfEle->pt();
  float ele_lip = _myGsfEle->vz();   

  for(TrackCollection::const_iterator this_track = _tracks.begin(); this_track != _tracks.end(); this_track++ ){ 
    
    CLHEP::Hep3Vector trackPAtVtx(this_track->px(),this_track->py(),this_track->pz());
    float this_pt  = trackPAtVtx.perp();

    // only tracks from the same vertex as the electron
    float this_lip = this_track->vz();
    if ( fabs(this_lip - ele_lip) > 0.2 ){ continue; }

    double dr = elePAtVtx.deltaR(trackPAtVtx);
    if ( fabs(dr) < _extRadius && fabs(dr) > _intRadius ){ dummyPt += this_pt; } 
    
  } //end loop over tracks		       
  
  // sum tracks pt / ele pt
  dummyPt = dummyPt/ele_pt; 

  return dummyPt;
}

// minimum distance from tracks upper a given pt cut - without veto
float EleTrackerIsolationAlgo::minDeltaR(float minPt) const
{
  float minDR = 100000. ;

  CLHEP::Hep3Vector elePAtVtx(_myGsfEle->px(), _myGsfEle->py(), _myGsfEle->pz()); 
  float ele_lip = _myGsfEle->vz();   

  for(TrackCollection::const_iterator this_track = _tracks.begin(); this_track != _tracks.end(); this_track++ ){ 
    CLHEP::Hep3Vector trackPAtVtx(this_track->px(),this_track->py(),this_track->pz());
    float this_pt  = trackPAtVtx.perp();

    if (this_pt < minPt){ continue;} 

    // only tracks from the same vertex as the electron
    float this_lip = this_track->vz();
    if ( fabs(this_lip - ele_lip) > 0.2 ){ continue; }
    
    double dr = elePAtVtx.deltaR(trackPAtVtx);
    if ( fabs(dr) < minDR ){ minDR = dr; }
 
  } //end loop over tracks		       
  
  return minDR;
}

// minimum distance from tracks upper a given pt cut - with veto
float EleTrackerIsolationAlgo::minDeltaR_withVeto(float minPt) const
{
  float minDR = 100000. ;

  CLHEP::Hep3Vector elePAtVtx(_myGsfEle->px(), _myGsfEle->py(), _myGsfEle->pz()); 
  float ele_lip = _myGsfEle->vz();   

  for(TrackCollection::const_iterator this_track = _tracks.begin(); this_track != _tracks.end(); this_track++ ){ 
    
    CLHEP::Hep3Vector trackPAtVtx(this_track->px(),this_track->py(),this_track->pz());
    float this_pt  = trackPAtVtx.perp();
    if (this_pt < minPt){ continue;} 

    // only tracks from the same vertex as the electron
    float this_lip = this_track->vz();
    if ( fabs(this_lip - ele_lip) > 0.2 ){ continue; }
    
    double dr = elePAtVtx.deltaR(trackPAtVtx);
    if ( (fabs(dr) < minDR) && (fabs(dr)>_intRadius) ){ minDR = dr; }
 
  } //end loop over tracks		       
  
  return minDR;
}

bool EleTrackerIsolationAlgo::isIsolated (float ptCut) const
{
  bool dummyIsolation = true ;
  
  if (EleTrackerIsolationAlgo::getPtTracks() > ptCut ) // default = 0.05 
    dummyIsolation = false ;
  
  return dummyIsolation ;
}
