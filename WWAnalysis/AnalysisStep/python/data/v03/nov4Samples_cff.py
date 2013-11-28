import FWCore.ParameterSet.Config as cms
from glob import glob

# Prefix where we expect idXX/*.root files to be
prefix = '/nfs/bluearc/group/trees/ww/V03_WWEventSkimProd_Newest/'

nov4Samples = {
      #id    #name      #scale to 1pb^-1 #is data  #skim event files
    'id65': ['Mu2010ANov4'        , 1, True, [] , 1.0, False],
    'id66': ['EG2010ANov4'        , 1, True, [] , 1.0, False],
    'id67': ['Mu2010BNov4'        , 1, True, [] , 1.0, False],
    'id68': ['Electron2010BNov4'  , 1, True, [] , 1.0, False],
}


for key in nov4Samples:
    nov4Samples[key][3] += [ '%s'%x for x in glob(prefix+key+'/*.root') ]























