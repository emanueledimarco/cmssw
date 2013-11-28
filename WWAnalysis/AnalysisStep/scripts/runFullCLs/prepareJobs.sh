#!/bin/bash


prefix="../datacards/hww-4.63fb.mH"
#type="comb_0j1j2j"
type="TYPE"
postfix="${type}_cb.txt"

ncores=8

### TO CONVERT FROM OUR OLD FORMAT TO THE NEW ONE
#for x in limitSummary_cutBased_*; do mv $x ${x/y_cutBased/y}; done;
#for x in limitSummary_*; do mv $x ${x/.txt/_cb.txt}; done;
#mv limitSummary_comb_cb.txt limitSummary_comb_0j1j2j_cb.txt

if [ "$1" == "10" -o "$1" == "init" ]
then

while read line
do
    mass=`echo $line|awk '{print $1}'`
    max=`echo $line|awk '{if($2>$3){print $2*4.}else{print $3*4.};}'`
    min=`echo $line|awk '{if($2<$3){print $2*0.2}else{print $3*0.2};}'`
    card=$prefix$mass.$postfix
    rm -rf m$mass
    mkdir m$mass
    cd m$mass
    #FAST (about 2000 crab jobs)
    #python ../makeGridUsingCrab.py ../$card $min $max -n 20 -j 100 -T 8  -t 100 -m $mass -r --log --glidein --server -o ${type}_m$mass
    #SLOW (about 200 crab jobs)
    python ../makeGridUsingCrab.py ../$card $min $max -n 25 -j 15 -T 80  -t 100 -m $mass -r --log --glidein --server -o ${type}_m$mass
    cd ..
done < ../limitSummary_$postfix

fi



if [ "$1" == "2" -o "$1" == "submit" ]
then
    for x in m???; do
	cd $x
	ls
	rm -rf crab_0_$type_*;
	crab -create -submit -cfg ${type}_m???.cfg >& submit.log &
	cd ..
	while [ `ps | grep python | wc -l` -ge $ncores ] ; do sleep 10; done
    done
fi

if [ "$1" == "2" -o "$1" == "onlySubmit" ]
then
    for x in m???; do
	cd $x
	ls
	crab -submit -cfg ${type}_m???.cfg >& submit.log &
	cd ..
	while [ `ps | grep python | wc -l` -ge $ncores ] ; do sleep 10; done
    done
fi


if [ "$1" == "3" -o "$1" == "statusGet" ]
then
    for x in m???; do
	cd $x
	crab -status -get >& get.log &
	cd ..
	while [ `ps | grep python | wc -l` -ge $ncores ] ; do sleep 10; done
    done
fi

if [ "$1" == "4" -o "$1" == "status" ]
then
    for x in m???; do
	cd $x
	crab -status  >& status.log &
	cd ..
	while [ `ps | grep python | wc -l` -ge $ncores ] ; do sleep 10; done
    done
fi


if [ "$1" == "5" -o "$1" == "print" ]
then
    for x in m???; do
	cd $x
	ls *.cfg
	cat status.log |grep " Y "
	cat status.log |grep " N "
	cd ..
    done
fi



if [ "$1" == "6" -o "$1" == "merge"  ]
then
    for x in m???; do hadd tree_$x.root $x/crab_*/res/*.root; done;
fi

if [ "$1" == "7" -o "$1" == "limits" ]
then
    for x in m???; do
	mass=`echo $x|awk -Fm '{print $2}'` 
	card=$prefix$mass.$postfix
	combine $card -M HybridNew -m $mass --freq --grid=tree_$x.root >& obs_$x.log &
	combine $card -M HybridNew -m $mass --freq --grid=tree_$x.root --expectedFromGrid 0.5    >& exp_$x.log &
	combine $card -M HybridNew -m $mass --freq --grid=tree_$x.root --expectedFromGrid 0.16   >& 68d_$x.log &
	combine $card -M HybridNew -m $mass --freq --grid=tree_$x.root --expectedFromGrid 0.84   >& 68u_$x.log &
	combine $card -M HybridNew -m $mass --freq --grid=tree_$x.root --expectedFromGrid 0.0275 >& 95d_$x.log &
	combine $card -M HybridNew -m $mass --freq --grid=tree_$x.root --expectedFromGrid 0.975  >& 95u_$x.log &
	while [ `ps | grep combine | wc -l` -ge $ncores ] ; do sleep 10; done
    done
fi


if [ "$1" == "8" -o "$1" == "tables" ]
then
    echo -e "mass \t obs \t exp \t 68d \t 68u \t 95d \t 95u" > limitSummary_final_$type.txt
    echo "-------------------------------------------------------" >> limitSummary_final_$type.txt
    for x in m???; do
	mass=`echo $x|awk -Fm '{print $2}'` 
	obs=`cat obs_$x.log|grep Limit|tail -n1 | awk '{print $4}'`
	exp=`cat exp_$x.log|grep Limit|tail -n1 | awk '{print $4}'`
	b68d=`cat 68d_$x.log|grep Limit|tail -n1 | awk '{print $4}'`
	b68u=`cat 68u_$x.log|grep Limit|tail -n1 | awk '{print $4}'`
	b95d=`cat 95d_$x.log|grep Limit|tail -n1 | awk '{print $4}'`
	b95u=`cat 95u_$x.log|grep Limit|tail -n1 | awk '{print $4}'`
	printf "%3.0f \t %1.2f \t %1.2f \t %1.2f \t %1.2f \t %1.2f \t %1.2f \n" $mass $obs $exp ${b68d} ${b68u} ${b95d} ${b95u} >> limitSummary_final_$type.txt
    done
fi

