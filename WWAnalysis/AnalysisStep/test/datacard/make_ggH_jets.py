import os, re
from math import *

def file2map(file):
    map = {}
    for line in open(file,"r"):
        fields = [x.strip() for x in line.split()]
        if len(fields) > 1 and re.match("\d+", fields[0]):
            map[float(fields[0])] = [float(y) for y in fields[1:] if re.match(r"-?\d+.*",y) ]
    return map

SYST_PATH = os.getenv('CMSSW_BASE')+'/src/WWAnalysis/AnalysisStep/test/datacard/'
YR_ggH  = file2map(SYST_PATH+'YR-ggH.txt')
YR_sigma  = dict([(m, xs) for m,(xs,xs_hi,xs_lo,sca_hi,sca_lo,pdf_hi,pdf_lo) in YR_ggH.items()] )
YR_kappa  = dict([(m, sqrt((1+0.01*sca_hi)/(1+0.01*sca_lo))) for m,(xs,xs_hi,xs_lo,sca_hi,sca_lo,pdf_hi,pdf_lo) in YR_ggH.items()] )

for X in 450, 550:
    YR_sigma[X] = sqrt(YR_sigma[X-10]*YR_sigma[X+10])
    YR_kappa[X] = 0.5*(YR_kappa[X-10]+YR_kappa[X+10])

MCFM_1jA = file2map(SYST_PATH+'MCFM_HWW_1jets_noBR.Nominal.txt')
MCFM_1jB = file2map(SYST_PATH+'MCFM_HWW_1jets.band.Nominal.txt')
MCFM_1jC = file2map(SYST_PATH+'MCFM_HWW_1jets.band.Fdn_Rdn.txt')

MCFM_2j = file2map(SYST_PATH+'MCFM_HWW_2jets.band.Fdn_Rdn.txt')

MCFM_1j_kappa = dict([(m, kappa)                for m,(cen,min,max,kappa) in MCFM_1jC.items()])
MCFM_1j_sigma = {}
for m in MCFM_1jA.keys(): MCFM_1j_sigma[m] = 1e-3*MCFM_1jA[m][0] * (MCFM_1jC[m][0]/MCFM_1jB[m][0])

MCFM_2j_kappa = dict([(m, kappa)    for m,(cen,min,max,kappa) in MCFM_2j.items()])
MCFM_2j_sigma = dict([(m, 1e-3*cen) for m,(cen,min,max,kappa) in MCFM_2j.items()])

masses = MCFM_2j_sigma.keys(); masses.sort()
print "%3s   %8s  %8s  %8s    %4s  %4s  %4s" % ('M_H', 'xs(tot)', 'xs(>=1)', 'xs(>=2)', 'ktot', 'k>=1', 'k>=2')
for m in masses:
    print "%3.0f   %8.2f  %8.2f  %8.2f    %4.2f  %4.2f  %4.2f" % (m, YR_sigma[m], MCFM_1j_sigma[m], MCFM_2j_sigma[m], YR_kappa[m], MCFM_1j_kappa[m], MCFM_2j_kappa[m])

print ""
out = open("ggH_jetBins.txt","w")
out.write("mH   f_0   f_1   f_2   k_1   k_2\n")
print     "mH   f_0   f_1   f_2   k_1   k_2  " 
for m in masses:
    f0 = (YR_sigma[m]      - MCFM_1j_sigma[m])/YR_sigma[m]
    f1 = (MCFM_1j_sigma[m] - MCFM_2j_sigma[m])/YR_sigma[m]
    f2 = 1-f0-f1
    print     "%3d %5.2f %5.2f %5.2f  %4.2f  %4.2f  " % (m, f0, f1, f2, MCFM_1j_kappa[m], MCFM_2j_kappa[m])
    out.write("%3d %5.2f %5.2f %5.2f  %4.2f  %4.2f\n" % (m, f0, f1, f2, MCFM_1j_kappa[m], MCFM_2j_kappa[m]))
