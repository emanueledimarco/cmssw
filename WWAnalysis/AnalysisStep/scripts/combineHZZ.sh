#!/bin/bash    

createJob() {
if [ -d $4 ]; then 
    echo "$4 exists"
else
    mkdir $4
fi
cd $4
echo "#!/bin/bash" > `echo "job"$5".sh"`
echo "cd $PWD" >> `echo "job"$5".sh"`
echo "export SCRAM_ARCH=slc5_amd64_gcc462" >> `echo "job"$5".sh"`
echo "eval \`scramv1 runtime -sh\`" >> `echo "job"$5".sh"`

cp ../`echo "card*m"$5"*"` ./
if [ $1 == '7TeV' -o $1 == '8TeV' ]; then
    combineCards.py `echo "card_"$2"_m"$5"_"$1"_4mu.txt"` `echo "card_"$2"_m"$5"_"$1"_4e.txt"` `echo "card_"$2"_m"$5"_"$1"_2e2mu.txt"` > `echo "card_"$2"_m"$5"_"$1".txt"`  
elif [ $1 == '78TeV' ]; then
    combineCards.py `echo "card_"$2"_m"$5"_7TeV_4mu.txt"` `echo "card_"$2"_m"$5"_7TeV_4e.txt"` `echo "card_"$2"_m"$5"_7TeV_2e2mu.txt"` `echo "card_"$2"_m"$5"_8TeV_4mu.txt"` `echo "card_"$2"_m"$5"_8TeV_4e.txt"` `echo "card_"$2"_m"$5"_8TeV_2e2mu.txt"`> `echo "card_"$2"_m"$5"_"$1".txt"`
fi

if [ $3 == 'limits' ]; then
    echo "combine -M Asymptotic -m "$5" card_"$2"_m"$5"_"$1".txt > limits_"$2"_m"$5"_"$1".txt" >> `echo "job"$5".sh"`
elif [ $3 == 'signif'  ]; then echo "combine -M ProfileLikelihood -m "$5" --signif card_"$2"_m"$5"_"$1".txt --rMax=$6 > signif_"$2"_m"$5"_"$1".txt" >> `echo "job"$5".sh"`
elif [ $3 == 'fit'  ]; then echo "combine -M MaxLikelihoodFit -m "$5" --justFit card_"$2"_m"$5"_"$1".txt --rMin -2 --rMax 2 > fit_"$2"_m"$5"_"$1".txt" >> `echo "job"$5".sh"`
else echo "combine -M ProfileLikelihood --significance -m "$5" card_"$2"_m"$5"_"$1".txt -t -1 --expectSignal=1 --toysFreq > expsignif_"$2"_m"$5"_"$1".txt" >> `echo "job"$5".sh"`
fi
chmod 755 `echo "job"$5".sh"`

bsub -q 8nh -W 60 -M 200000 `echo "job"$5".sh"`
cd ../
}


for X in $(seq 110 1 160)
    do
        createJob $1 $2 $3 job $X 1
    done

for X in $(seq 162 2 180)
    do
        createJob $1 $2 $3 job $X 1
    done

for X in $(seq 182 2 290)
    do
        createJob $1 $2 $3 job $X 1
    done

for X in $(seq 295 5 350)
    do
        createJob $1 $2 $3 job $X 1
    done

for X in $(seq 360 10 400)
    do
        createJob $1 $2 $3 job $X 5
    done

for X in $(seq 420 20 600)
    do
        createJob $1 $2 $3 job $X 10
    done

for X in $(seq 620 20 1000)
    do
        createJob $1 $2 $3 job $X 50
    done

