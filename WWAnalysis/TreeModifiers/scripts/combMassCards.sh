#!/bin/bash

dirs='cards/c1DNoMassErr cards/c1DMassErr cards/c2DMassErr'

echo "Combining m(4l) cards..."
for dir in $dirs 
  do
  cd $dir
  pwd=$PWD
  echo $pwd
  twod='2D'
  let ndim=1
  if [[ $dir =~ $twod ]]
      then
      let ndim=2
  fi
# sed -i 's|card|'$pwd'/card|g' card*txt
  combineCards.py -S hzz4l_7TeV_4mu=card_$ndim\D_m126_7TeV_4mu.txt     hzz4l_8TeV_4mu=card_$ndim\D_m126_8TeV_4mu.txt     > card_$ndim\D_m126_4mu.txt
  combineCards.py -S hzz4l_7TeV_4e=card_$ndim\D_m126_7TeV_4e.txt       hzz4l_8TeV_4mu=card_$ndim\D_m126_8TeV_4e.txt      > card_$ndim\D_m126_4e.txt
  combineCards.py -S hzz4l_7TeV_2e2mu=card_$ndim\D_m126_7TeV_2e2mu.txt hzz4l_8TeV_2e2mu=card_$ndim\D_m126_8TeV_2e2mu.txt > card_$ndim\D_m126_2e2mu.txt
  
  combineCards.py -S hzz4l_7TeV_4mu=card_$ndim\D_m126_7TeV_4mu.txt hzz4l_7TeV_4e=card_$ndim\D_m126_7TeV_4e.txt hzz4l_7TeV_2e2mu=card_$ndim\D_m126_7TeV_2e2mu.txt \
      hzz4l_8TeV_4mu=card_$ndim\D_m126_8TeV_4mu.txt hzz4l_8TeV_4e=card_$ndim\D_m126_8TeV_4e.txt hzz4l_8TeV_2e2mu=card_$ndim\D_m126_8TeV_2e2mu.txt \
      > card_$ndim\D_m126_comb.txt
  
  text2workspace.py card_$ndim\D_m126_4mu.txt    -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingHiggsMass --PO higgsMassRange=120,130 -o FloatMass_4mu_hzz.root
  text2workspace.py card_$ndim\D_m126_4e.txt     -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingHiggsMass --PO higgsMassRange=120,130 -o FloatMass_4e_hzz.root
  text2workspace.py card_$ndim\D_m126_2e2mu.txt  -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingHiggsMass --PO higgsMassRange=120,130 -o FloatMass_2e2mu_hzz.root
  
  text2workspace.py card_$ndim\D_m126_comb.txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingHiggsMass --PO higgsMassRange=120,130 -o FloatMass_comb_hzz.root
  cd -
done;

echo "done"
