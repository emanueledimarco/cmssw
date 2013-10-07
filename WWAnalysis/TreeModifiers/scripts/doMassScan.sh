#!/bin/bash

dirs='cards/c1DNoMassErr cards/c1DMassErr cards/c2DMassErr'
channels='4mu 4e 2e2mu comb'

for dir in $dirs
  do
  echo $dir
  twod='2D'
  nomasserr='NoMassErr'
  let ndim=1
  if [[ $dir =~ $twod ]]
      then
      let ndim=3
  elif [[ $dir =~ $nomasserr ]]
      then
      let ndim=1
  else
      let ndim=2
  fi
  for channel in $channels 
    do
    echo "   scanning channel HZZ: " $channel
    echo "#!/bin/bash" > SCAN$ndim\DFast$channel\.src
    echo "cd " $PWD >> SCAN$ndim\DFast$channel\.src
    echo "eval \`scram ru -sh\`" >> SCAN$ndim\DFast$channel\.src
    echo "combine -M MultiDimFit $PWD/$dir/FloatMass_"$channel"_hzz.root -m 125.8 -P MH --floatOtherPOI=1 --algo=grid --points=200 -n SCAN"$ndim"D"$channel"Fast" >> SCAN$ndim\DFast$channel\.src
    echo "mv higgsCombine*.MultiDimFit*root " $PWD >> SCAN$ndim\DFast$channel\.src
    bsub -q 8nh -J $ndim\D$channel -o SCAN$ndim\DFast$channel\.log source SCAN$ndim\DFast$channel\.src
  done
done

echo "submitted. Running in batch"

