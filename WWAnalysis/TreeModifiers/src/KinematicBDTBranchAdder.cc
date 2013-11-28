#include "WWAnalysis/TreeModifiers/interface/KinematicBDTBranchAdder.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRotation.h"

KinematicBDTBranchAdder::KinematicBDTBranchAdder(const edm::ParameterSet& pset) : BranchAdder(pset){
  std::string datasetConf_     = pset.getParameter<std::string>("datasetConf"); 
  std::string weightfile     = pset.getParameter<std::string>("weightfile");   
  
      BDTReader = new TMVA::Reader( "V" );
      BDTReader->AddVariable( "costheta1",         &MVAInputVar_costheta1);
      BDTReader->AddVariable( "costheta2",         &MVAInputVar_costheta2);
      BDTReader->AddVariable( "costhetastar",      &MVAInputVar_costhetastar);
      BDTReader->AddVariable( "Phi",		   &MVAInputVar_Phi);
      BDTReader->AddVariable( "Phi1",		   &MVAInputVar_Phi1);
      BDTReader->AddVariable( "mZ1",		   &MVAInputVar_mZ1);
      BDTReader->AddVariable( "mZ2",		   &MVAInputVar_mZ2);
      BDTReader->AddVariable( "ZZpt/m4l",	   &MVAInputVar_ZZpt);
      BDTReader->AddVariable( "ZZdotZ1/(m4l*mZ1)", &MVAInputVar_zzdotz1);
      BDTReader->AddVariable( "ZZdotZ2/(m4l*mZ2)", &MVAInputVar_zzdotz2);
      BDTReader->AddVariable( "ZZptCosDphiZ1pt",   &MVAInputVar_ZZptZ1ptCosDphi);
      BDTReader->AddVariable( "ZZptCosDphiZ2pt",   &MVAInputVar_ZZptZ2ptCosDphi);
      BDTReader->AddVariable( "Z1pt/m4l",	   &MVAInputVar_reduced_Z1pt);
      BDTReader->AddVariable( "Z2pt/m4l",	   &MVAInputVar_reduced_Z2pt);
      BDTReader->AddVariable( "ZZy",		   &MVAInputVar_ZZy);
      // add spectators
      BDTReader->AddSpectator("m4l",               &MVAInputVar_m4l);
      // initialize
      BDTReader->BookMVA("BDTG", weightfile.c_str());

}

void KinematicBDTBranchAdder::initialize(TTree& tree) {
  tree.SetBranchAddress("l1pt",&l1pt);
  tree.SetBranchAddress("l1eta",&l1eta);
  tree.SetBranchAddress("l1phi",&l1phi);
  tree.SetBranchAddress("l1pdgId",&l1pdgId);
  tree.SetBranchAddress("l2pt",&l2pt);
  tree.SetBranchAddress("l2eta",&l2eta);
  tree.SetBranchAddress("l2phi",&l2phi);
  tree.SetBranchAddress("l2pdgId",&l2pdgId);
  tree.SetBranchAddress("pho1pt",&pho1pt);
  tree.SetBranchAddress("pho1eta",&pho1eta);
  tree.SetBranchAddress("pho1phi",&pho1phi);

  tree.SetBranchAddress("l3pt",&l3pt);
  tree.SetBranchAddress("l3eta",&l3eta);
  tree.SetBranchAddress("l3phi",&l3phi);
  tree.SetBranchAddress("l3pdgId",&l3pdgId);
  tree.SetBranchAddress("l4pt",&l4pt);
  tree.SetBranchAddress("l4eta",&l4eta);
  tree.SetBranchAddress("l4phi",&l4phi);
  tree.SetBranchAddress("l4pdgId",&l4pdgId);
  tree.SetBranchAddress("pho2pt",&pho2pt);
  tree.SetBranchAddress("pho2eta",&pho2eta);
  tree.SetBranchAddress("pho2phi",&pho2phi);

}


float KinematicBDTBranchAdder::calculateValue(TTree& tree,int entry, int id, double& value)  {
  tree.GetEntry(entry);

  TLorentzVector lep1vec;
  TLorentzVector lep2vec;
  TLorentzVector lep3vec;
  TLorentzVector lep4vec;
  double mass1 = 0; 
  double mass2 = 0;
  double mass3 = 0;
  double mass4 = 0;
  if (abs(l1pdgId) == 11) mass1 = 0.51099892e-3;
  else if (abs(l1pdgId) == 13) mass1 = 105.658369e-3;
  if (abs(l2pdgId) == 11) mass2 = 0.51099892e-3;
  else if (abs(l2pdgId) == 13) mass2 = 105.658369e-3;
  if (abs(l3pdgId) == 11) mass3 = 0.51099892e-3;
  else if (abs(l3pdgId) == 13) mass3 = 105.658369e-3;
  if (abs(l4pdgId) == 11) mass4 = 0.51099892e-3;
  else if (abs(l4pdgId) == 13) mass4 = 105.658369e-3;


  //Order Leptons by charge
  //lep1 and lep3 are negatively charged
  //lep2 and lep4 are positively charged

  if (l1pdgId > 0) {
    lep1vec.SetPtEtaPhiM(l1pt, l1eta, l1phi, mass1);
    lep2vec.SetPtEtaPhiM(l2pt, l2eta, l2phi, mass2);
  } else {
    lep1vec.SetPtEtaPhiM(l2pt, l2eta, l2phi, mass2);
    lep2vec.SetPtEtaPhiM(l1pt, l1eta, l1phi, mass1);
  }

  if (l3pdgId > 0) {
    lep3vec.SetPtEtaPhiM(l3pt, l3eta, l3phi, mass3);
    lep4vec.SetPtEtaPhiM(l4pt, l4eta, l4phi, mass4);
  } else {
    lep3vec.SetPtEtaPhiM(l4pt, l4eta, l4phi, mass4);
    lep4vec.SetPtEtaPhiM(l3pt, l3eta, l3phi, mass3);
  }

  //do FSR recovery
  TLorentzVector pho1vec;
  TLorentzVector pho2vec;
  if (pho1pt > 0) {
    pho1vec.SetPtEtaPhiM(pho1pt, pho1eta, pho1phi, 0);
    if (reco::deltaR(lep1vec.Eta(),lep1vec.Phi(), pho1eta,pho1phi) < reco::deltaR(lep2vec.Eta(),lep2vec.Phi(), pho1eta,pho1phi) ) {
      lep1vec += pho1vec;
    } else {
      lep2vec += pho1vec;
    }
  }
  if (pho2pt > 0) {
    pho2vec.SetPtEtaPhiM(pho2pt, pho2eta, pho2phi, 0);
    if (reco::deltaR(lep3vec.Eta(),lep3vec.Phi(), pho2eta,pho2phi) < reco::deltaR(lep4vec.Eta(),lep4vec.Phi(), pho2eta,pho2phi) ) {
      lep3vec += pho2vec;
    } else {
      lep4vec += pho2vec;
    }
  }
  

  //Compute Angles
  TLorentzVector z1vec = lep1vec+lep2vec;
  TLorentzVector z2vec = lep3vec+lep4vec;
  TLorentzVector zzvec = z1vec+z2vec;
  
  TVector3 Xframe  = zzvec.BoostVector();
  TVector3 Z1frame = z1vec.BoostVector();
  TVector3 Z2frame = z2vec.BoostVector();

  //  
  // costheta(l1m,X) in Z1 frame
  // --------------------------------------------------------------------------------------------------
  TLorentzVector vecl1m_in_Z1frame( lep1vec ); // make a copies of the vectors 
  TLorentzVector vec4l_in_Z1frame( zzvec );   // 
  vecl1m_in_Z1frame.Boost( -1*Z1frame );        // now go to the Z1 rest frame 
  vec4l_in_Z1frame.Boost( -1*Z1frame );         //
  MVAInputVar_costheta1 
    = vecl1m_in_Z1frame.Vect().Dot(-1*vec4l_in_Z1frame.Vect()) / 
    (vecl1m_in_Z1frame.Vect().Mag()*vec4l_in_Z1frame.Vect().Mag());

  //  
  // costheta(l3,X) in Z2 frame
  // --------------------------------------------------------------------------------------------------
  TLorentzVector vecl2m_in_Z2frame( lep3vec ); // make copies of the vectors
  TLorentzVector vec4l_in_Z2frame( zzvec );   // 
  vecl2m_in_Z2frame.Boost( -1*Z2frame );        // now go to the Z2 rest frame 
  vec4l_in_Z2frame.Boost( -1*Z2frame );  
  MVAInputVar_costheta2 = 
    vecl2m_in_Z2frame.Vect().Dot(-1*vec4l_in_Z2frame.Vect()) / 
    (vecl2m_in_Z2frame.Vect().Mag()*vec4l_in_Z2frame.Vect().Mag());


  //  
  // Phi(Z1,Z2) in X frame 
  // --------------------------------------------------------------------------------------------------

  // this is essentially what they do and it reproduces their numbers, 
  // it's the angle between l1 and the z2 decay plane in the Z1 frame
  TLorentzVector vecz2_in_Z1frame( z2vec );      //
  TLorentzVector vecl2m_in_Z1frame( lep3vec );    // 
  TLorentzVector vecl2p_in_Z1frame( lep4vec );    // 
  vecz2_in_Z1frame.Boost( -1*Z1frame );            //
  vecl2m_in_Z1frame.Boost( -1*Z1frame );           // 
  vecl2p_in_Z1frame.Boost( -1*Z1frame );           // 

  //
  // define the Z2 coordinate system in the Z1 frame
  // x : opposite the Z2 dir, lies in the Z2 decay plane
  // z : orthogonal to the Z2 decay plane
  // y : orthogonal to x&z, lies in the Z2 decay plane
  TVector3 nnewX = ( -1*vecz2_in_Z1frame.Vect()).Unit(); 
  TVector3 nnewZ = (vecl2m_in_Z1frame.Vect().Cross(vecl2p_in_Z1frame.Vect())).Unit();
  TVector3 nnewY = nnewZ.Cross(nnewX).Unit();       
  TRotation rotation2;
  rotation2 = (rotation2.RotateAxes( nnewX, nnewY, nnewZ )).Inverse();
  TVector3 vecl1m_in_Z2coords_Z1frame  = vecl1m_in_Z1frame.Vect().Transform(rotation2);

  //
  // axes are redfined so that Phi is the angle between the Z1 orthogonal 
  // and the Z2 decay plane as given by y
  //
  TRotation rotation3;
  TVector3 nnnewX( 0,1,0 ); 
  TVector3 nnnewY( 0,0,1 ); 
  TVector3 nnnewZ( 1,0,0 );
  rotation3 = (rotation3.RotateAxes( nnnewX, nnnewY, nnnewZ )).Inverse();
  TVector3 vecl1m_in_Z2coords_redefined_Z1frame = vecl1m_in_Z2coords_Z1frame.Transform(rotation3);
  MVAInputVar_Phi = vecl1m_in_Z2coords_redefined_Z1frame.Phi(); 


  //
  // rather, I'd expect something like this ...
  //
  TLorentzVector vecl1p_in_Xframe( lep1vec ); // make a copies of the vectors 
  TLorentzVector vecl1m_in_Xframe( lep2vec ); // 
  TLorentzVector vecl2p_in_Xframe( lep3vec ); // 
  TLorentzVector vecl2m_in_Xframe( lep4vec ); // 

  vecl1p_in_Xframe.Boost( -1*Xframe );         // go to the X rest frame 
  vecl1m_in_Xframe.Boost( -1*Xframe );         // 
  vecl2p_in_Xframe.Boost( -1*Xframe );         // 
  vecl2m_in_Xframe.Boost( -1*Xframe );         // 

  TVector3 Z1ortho_in_Xframe =                                         // define the decay planes via the orthogonals 
    vecl1m_in_Xframe.Vect().Cross(vecl1p_in_Xframe.Vect()); 
  TVector3 Z2ortho_in_Xframe = 
    vecl2m_in_Xframe.Vect().Cross(vecl2p_in_Xframe.Vect());

  //  
  // costheta*(Z1,pp) in X frame 
  // --------------------------------------------------------------------------------------------------
  TLorentzVector vecz1_in_Xframe( z1vec );        // make a copies of the vectors   
  TLorentzVector vecpp_in_Xframe( 0., 0., 1., 0. ); // define pp axis

  vecz1_in_Xframe.Boost( -1*Xframe );         // go to the X rest frame 
  vecpp_in_Xframe.Boost( -1*Xframe );         // 

  MVAInputVar_costhetastar = 
    vecz1_in_Xframe.Vect().Dot(vecpp_in_Xframe.Vect()) /
    (vecz1_in_Xframe.Vect().Mag()*vecpp_in_Xframe.Vect().Mag());


  //  
  // Phi1(Z1,pp) in X frame 
  // --------------------------------------------------------------------------------------------------

  // define the Z1 decay plane in the Xframe
  // y: orthogonal to the Z1 decay plane 
  // x: lies in the Z1 decay plane, but orthogonal to Z1
  // z: z1 in Xframe 
  TVector3 newY = (vecl1p_in_Xframe.Vect().Unit().Cross(vecl1m_in_Xframe.Vect().Unit())).Unit(); 
  TVector3 newX = (newY.Cross(vecz1_in_Xframe.Vect().Unit())).Unit(); 
  TVector3 newZ( (vecz1_in_Xframe.Vect()).Unit() ); 

  TRotation rotation;
  rotation = (rotation.RotateAxes( newX, newY, newZ )).Inverse();
  TVector3 vecpp( 0., 0., 1. );
  TVector3 vecpp_in_Z1coords_Xframe  = vecpp.Transform(rotation);
  MVAInputVar_Phi1 = vecpp_in_Z1coords_Xframe.Phi(); 


  //  
  // More kinematics variables
  // --------------------------------------------------------------------------------------------------
  TVector3 zzvec3 = zzvec.Vect();
  TVector3 z1vec3 = z1vec.Vect();
  TVector3 z2vec3 = z2vec.Vect();
  TVector2 zzvecxy = zzvec3.XYvector();
  TVector2 z1vecxy = z1vec3.XYvector();
  TVector2 z2vecxy = z2vec3.XYvector();
  MVAInputVar_ZZptZ1ptCosDphi = TMath::Cos(zzvecxy.DeltaPhi( z1vecxy));
  MVAInputVar_ZZptZ2ptCosDphi = TMath::Cos(zzvecxy.DeltaPhi( z2vecxy));
  MVAInputVar_reduced_Z1pt = z1vec.Pt() / zzvec.M();
  MVAInputVar_reduced_Z2pt = z2vec.Pt() / zzvec.M();
  MVAInputVar_ZZy = zzvec.Rapidity();
  MVAInputVar_m4l = zzvec.M();

  value = BDTReader->EvaluateMVA("BDTG");

  return 0;
}


