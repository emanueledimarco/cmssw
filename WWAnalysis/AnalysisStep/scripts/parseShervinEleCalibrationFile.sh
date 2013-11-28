#!/bin/bash


lastNewCategory="null";
while read line; do
    category=$(echo $line |awk '{print  $1}');
    runSX=$(echo $line |awk '{print  $3}');
    runDX=$(echo $line |awk '{print  $4}');
    scale=$(echo $line |awk '{print  $5}');

    if [ "$category" != "$lastNewCategory" ] 
    then
	if [ "$lastNewCategory" == "null" ]; 
	then echo -n "if "; 
	else echo -n "} else if "; 
	fi;
	    
	case $category in
	    "EB-absEta_0_1-bad" )
		echo "(electron.isEB() && fabs(electron.superCluster()->eta())<1 && r9<0.94) {";;
	    "EB-absEta_0_1-gold" )
		echo "(electron.isEB() && fabs(electron.superCluster()->eta())<1 && r9>=0.94) {";;
	    "EB-absEta_1_1.4442-bad" )
		echo "(electron.isEB() && fabs(electron.superCluster()->eta())>=1 && r9<0.94) {";;
	    "EB-absEta_1_1.4442-gold" )
		echo "(electron.isEB() && fabs(electron.superCluster()->eta())>=1 && r9>=0.94) {";;
	    "EE-absEta_1.566_2-bad" )
		echo "(electron.isEE() && fabs(electron.superCluster()->eta())<2 && r9<0.94) {";;
	    "EE-absEta_1.566_2-gold" )
		echo "(electron.isEE() && fabs(electron.superCluster()->eta())<2 && r9>=0.94) {";;	    
	    "EE-absEta_2_2.5-bad" )
		echo "(electron.isEE() && fabs(electron.superCluster()->eta())>=2 && r9<0.94) {";;
	    "EE-absEta_2_2.5-gold" )
		echo "(electron.isEE() && fabs(electron.superCluster()->eta())>=2 && r9>=0.94) {";;	    
	esac			    
	lastNewCategory="$category";
    fi
    echo -e "\t if (run>="$runSX" && run<="$runDX") scale = "$scale";" ; 
done < shervin.txt
echo "}"