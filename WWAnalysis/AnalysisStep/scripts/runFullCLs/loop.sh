#!/bin/bash

#configurations="comb_0j1j2j"
#configurations="of_0j"
#configurations="sf_0j"
#configurations="of_1j"
#configurations="sf_1j"
#configurations="comb_2j"

#configurations="comb_0j"
configurations="comb_1j"

#configurations="comb_0j1j2j of_0j sf_0j of_1j sf_1j comb_2j"

if [ "$1" == "init" ]
then
    for x in $configurations;
    do
	mkdir $x
	cat prepareJobs.sh |sed "s#TYPE#$x#" > $x/prepareJobs.sh
	cp makeGridUsingCrab.py $x/
	cd $x	
	source prepareJobs.sh init
	cd ..
    done
fi


if [ "$1" == "submit" ]
then
    for x in $configurations;
    do
	cat prepareJobs.sh |sed "s#TYPE#$x#" > $x/prepareJobs.sh
	cd $x	
	source prepareJobs.sh submit
	cd ..
    done
fi

if [ "$1" == "onlySubmit" ]
then
    for x in $configurations;
    do
	cat prepareJobs.sh |sed "s#TYPE#$x#" > $x/prepareJobs.sh
	cd $x	
	source prepareJobs.sh onlySubmit
	cd ..
    done
fi


if [ "$1" == "statusGet" ]
then
    for x in $configurations;
    do
	cat prepareJobs.sh |sed "s#TYPE#$x#" > $x/prepareJobs.sh
	cd $x	
	source prepareJobs.sh statusGet
	cd ..
    done
fi

if [ "$1" == "status" ]
then
    for x in $configurations;
    do
	cat prepareJobs.sh |sed "s#TYPE#$x#" > $x/prepareJobs.sh
	cd $x	
	source prepareJobs.sh status
	cd ..
    done
fi

if [ "$1" == "print" ]
then
    for x in $configurations;
    do
	cat prepareJobs.sh |sed "s#TYPE#$x#" > $x/prepareJobs.sh
	cd $x	
	source prepareJobs.sh print
	cd ..
    done
fi


if [ "$1" == "merge" ]
then
    for x in $configurations;
    do
	cat prepareJobs.sh |sed "s#TYPE#$x#" > $x/prepareJobs.sh
	cd $x	
	source prepareJobs.sh merge
	cd ..
    done
fi


if [ "$1" == "limits" ]
then
    for x in $configurations;
    do
	cat prepareJobs.sh |sed "s#TYPE#$x#" > $x/prepareJobs.sh
	cd $x	
	source prepareJobs.sh limits
	cd ..
    done
fi

if [ "$1" == "tables" ]
then
    for x in $configurations;
    do
	cat prepareJobs.sh |sed "s#TYPE#$x#" > $x/prepareJobs.sh
	cd $x	
	source prepareJobs.sh tables
	cd ..
    done
fi

if [ "$1" == "plots" ]
then
    mkdir finalPlots
    for x in $configurations;
    do
	cp $x/limitSummary_final_${x}.txt finalPlots
    done
    cp makePlots.sh finalPlots
    cp PlotLimit.C finalPlots
    cp tdrstyle.C finalPlots
    cd finalPlots
    source makePlots.sh
    cd ..
fi






