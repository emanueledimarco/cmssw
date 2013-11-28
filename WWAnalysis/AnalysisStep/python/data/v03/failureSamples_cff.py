import FWCore.ParameterSet.Config as cms
from glob import glob

# Prefix where we expect idXX/*.root files to be
prefix = '/nfs/bluearc/group/trees/ww/V03_WWEventSkimProd_Newest/'

failureSamples = {
    # id          name                       scale to 1/pb   isMC       Enhancements
    'id101160.T00': ['ggH1602L2Nu_T00',   0.000003777578691 ,  False, [], 8.405E+00 , False],
    'id101160.T01': ['ggH1602L2Nu_T01',   0.000003748807287 ,  False, [], 8.405E+00 , False],
    'id101160.T02': ['ggH1602L2Nu_T02',   0.000003748807287 ,  False, [], 8.405E+00 , False],
    'id101160.T03': ['ggH1602L2Nu_T03',   0.000003748807287 ,  False, [], 8.405E+00 , False],
    'id101160.T04': ['ggH1602L2Nu_T04',   0.000003748807287 ,  False, [], 8.405E+00 , False],
    'id101160.T05': ['ggH1602L2Nu_T05',   0.000003748807287 ,  False, [], 8.405E+00 , False],
    'id101160.T06': ['ggH1602L2Nu_T06',   0.000003748807287 ,  False, [], 8.405E+00 , False],
    'id101160.T07': ['ggH1602L2Nu_T07',   0.000003748807287 ,  False, [], 8.405E+00 , False],
    'id101160.T08': ['ggH1602L2Nu_T08',   0.000003748807287 ,  False, [], 8.405E+00 , False],
    'id101160.T09': ['ggH1602L2Nu_T09',   0.000003748807287 ,  False, [], 8.405E+00 , False],
    'id101160.T10': ['ggH1602L2Nu_T10',   0.000003748807287 ,  False, [], 8.405E+00 , False],
    'id101160.T11': ['ggH1602L2Nu_T11',   0.000003748807287 ,  False, [], 8.405E+00 , False],
    'id101160.T12': ['ggH1602L2Nu_T12',   0.000003748807287 ,  False, [], 8.405E+00 , False],
    'id101160.T13': ['ggH1602L2Nu_T13',   0.000003748807287 ,  False, [], 8.405E+00 , False],
}


#loop over the looper
for key in failureSamples:
    failureSamples[key][3] += [ '%s'%x for x in glob(prefix+key+'/*.root') ]







































































































