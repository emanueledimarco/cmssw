#!/bin/env bash
whoami=$(whoami)
wait=0
nsub=0
for mass in $(ls -d m*)
do
  cd $mass
  for mode in `ls -d *`
  do
    cd $mode
    echo $mode
    sampleName=$(echo $mode|awk -F"." '{print $2}')
    echo "cmsRun $sampleName.py >& $sampleName.log &"
    cmsRun $sampleName.py >& $sampleName.log &
    nsub=$(($nsub+1))
    echo "#job submitted: " $nsub
    cd ..
    wait=1
    while [ $wait -eq 1 ]
    do
      njobs=$(ps aux |grep $whoami |grep cmsRun|grep -v grep|wc -l)
      if [[ $njobs -gt 8 ]];then #have to wait
	  echo "njobs running: "$njobs
	  echo "will sleep for 10s.."
	  sleep 10s
      else
	  wait=0
      fi
    done
  done
  cd ..
done
#while [ $wait -eq 1 ]
#do
#  date
#  sleep 3s
#done

