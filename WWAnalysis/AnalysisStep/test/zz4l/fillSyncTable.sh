if [[ "$1" == "" ]]; then echo "Usage: $1 logfile"; exit 1; fi
test -f $1 || echo "Usage: $1 logfile";

LOG=/dev/shm/$(basename $(mktemp) )

function getline { 
    MATCH=$1; shift
    ROWS=$1; shift;
    PASS=$1; shift
    grep "PATH\|$MATCH" $LOG | \
        perl -npe 'm/PATH/ and chomp' | \
        head -n $ROWS | \
        awk "{ORS=\"/\"; print \$$(($PASS+4))}" | sed 's#/$#\n#'
}

~gpetrucc/sh/skimreport --path=zzP -a $1 > $LOG
echo -n "Begin:               "; getline trigger      1 0
echo -n "Trigger:             "; getline trigger      1 1
echo -n "3a1 oneZ:            "; getline selectedZ1   3 0
echo -n "3b1 good Z1:         "; getline selectedZ1   3 1
echo -n "4a. 4 lepts:         "; getline oneZZ        4 0
echo -n "4b. one ZZ:          "; getline oneZZ        4 1
echo -n "4d. 4 < mz2 < 120    "; getline selectedZZs1 4 1
echo -n "5. pt 20/10          "; getline selectedZZs2 4 1
echo -n "6. 4/4 mll > 4       "; getline selectedZZs3 4 1
echo -n "7. m > 70            "; getline selectedZZs4 4 1
echo -n "8. m > 100, mz2 > 12 "; getline selectedZZs5 4 1
echo -n "9. MELA > 0.1        "; getline selectedZZs6 4 1

~gpetrucc/sh/skimreport --path=zzPath_1FSR -a $1 > $LOG
echo -n "9. + 1 FSR           "; getline countZZ1FSR 4 1
~gpetrucc/sh/skimreport --path=zzPath_2FSR -a $1 > $LOG
echo -n "9. + 2 FSR           "; getline countZZ2FSR 4 1

~gpetrucc/sh/skimreport --path=vbfPath -a $1 > $LOG
echo -n "9 VBF a. >= 1 jets   "; getline selectedZZsVBF1 4 1
echo -n "9 VBF b. == 2 jets   "; getline selectedZZsVBF3 4 1
echo -n "9 VBF c. deta, mjj   "; getline selectedZZsVBF4 4 1
rm $LOG
