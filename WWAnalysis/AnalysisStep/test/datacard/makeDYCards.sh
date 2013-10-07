#!/usr/bin/env bash

for C in ee mm; do for J in 0 1; do
    echo " ---- DYCard_${C}_${J}j.txt ----"
    echo "#Yield  Err(Data)   Err(MC) - both to be gmM" | tee DYCard_${C}_${J}j.txt
    # grep "N.Z->$C, all cuts, ${J}j" ../step3/dyEstimate.cuts*txt | \
    grep "N.Z->$C, all cuts, ${J}j" $* | \
        # perl -npe 's/.*cuts(\d+).*=\s+(\S+) .*-\s+(\S+).*-\s+(\S+).*/$1  $2 $3 $4/' | \
        # fix this, now it needs to know the prefix in the file name to the mass norm or cuts in mine and gio's case
        perl -npe 's/.*norm(\d+).*=\s+(\S+) .*-\s+(\S+).*-\s+(\S+).*/$1  $2 $3 $4/' | \
        tee -a DYCard_${C}_${J}j.txt
    done
done
