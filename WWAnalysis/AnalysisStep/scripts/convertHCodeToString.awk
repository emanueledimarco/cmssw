#!/bin/awk -f 

NF==1 {
    a=int($1%1000)
    b=int((($1-100000)/1000%10));
    c = ""
    if(b == 1) {c="GluGluTo";d="2L2Nu"; }
    if(b == 2) {c="GluGluTo";d="LNuTauNu"; }
    if(b == 3) {c="GluGluTo";d="2Tau2Nu"; }
    if(b == 4) {c="VBF_";d="2L2Nu"; }
    if(b == 5) {c="VBF_";d="LNuTauNu"; }
    if(b == 6) {c="VBF_";d="2Tau2Nu"; }

    e=int(($1%100000)/10000)
    if(e == 0) {f="NoPU";}
    if(e == 1) {f="PU";}
    printf "%sHToWWTo%s_M-%i\n",c,d,a
} 

NF==0 {
    print
}
