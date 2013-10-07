#!/bin/bash

listOfModuleToDump=""
listOfModuleToDump+=" hltFilter"
listOfModuleToDump+=" primaryVertexFilter"
listOfModuleToDump+=" noscraping"
listOfModuleToDump+=" incompleteECALReadout"
listOfModuleToDump+=" elmuPtFilter"
listOfModuleToDump+=" elelPtFilter"
listOfModuleToDump+=" mumuPtFilter"
#listOfModuleToDump=$listOfModuleToDump" \"TrigReport Event\""

echo "listOfModuleToDump: " $listOfModuleToDump

total=0;
passed=0;

listOfFiles=$(ls $1/*.stdout)

for module in $listOfModuleToDump; do
    Run=0;
    Pass=0;
    
    for file in $listOfFiles; do   
	Count=$(grep TrigReport $file| grep $module |tail -n1 |awk '{print $3" "$4}')
	Run=$(( $Run+$(echo $Count|awk '{print $1}') ))
	Pass=$(( $Pass+$(echo $Count|awk '{print $2}') ))
    done

    echo $module" Run: "  $Run
    echo $module" Pass: " $Pass
    echo $module" efficiency: " $(echo "scale=2; $Pass*100 / $Run" |bc)"%"
    echo;
done



Run=0;
Pass=0;
    
for file in $listOfFiles; do   
    Count=$(grep "TrigReport Event" $file|tail -n1 |awk '{print $5" "$8}')
    Run=$(( $Run+$(echo $Count|awk '{print $1}') ))
    Pass=$(( $Pass+$(echo $Count|awk '{print $2}') ))
done

echo "TOTAL Run: "  $Run
echo "TOTAL Pass: " $Pass
echo "TOTAL efficiency: " $(echo "scale=2; $Pass*100 / $Run" |bc)"%"
