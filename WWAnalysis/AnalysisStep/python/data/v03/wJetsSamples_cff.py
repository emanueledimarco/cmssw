import FWCore.ParameterSet.Config as cms
from glob import glob

# Prefix where we expect idXX/*.root files to be
prefix = '/nfs/bluearc/group/trees/ww/V03_WWEventSkimProd_Newest/'

wJetsSamples = {
      #id    #name      #scale to 1pb^-1  #skim event # is data files
    'id30': ['W0JetsD6T        ',   5.51E-03,    False,    [], 1.0, False],
    #'id31': ['W0JetsZ2         ',   DIV/0! ,    False,    [], 1.0, False],  Sample is still not ready!
    'id32': ['W1Jets0to100Z2   ',   4.18E-03,    False,    [], 1.0, False],
    'id33': ['W1Jets100to300Z2 ',   3.59E-04,    False,    [], 1.0, False],
    'id34': ['W1Jets300to800Z2 ',   5.74E-06,    False,    [], 1.0, False],
    'id35': ['W1Jets800to1600Z2',   1.12E-08,    False,    [], 1.0, False],
    'id36': ['W2Jets0to100Z2   ',   1.17E-03,    False,    [], 1.0, False],
    'id37': ['W2Jets100to300Z2 ',   8.13E-04,    False,    [], 1.0, False],
    'id38': ['W2Jets300to800Z2 ',   1.34E-05,    False,    [], 1.0, False],
    'id39': ['W2Jets800to1600Z2',   9.90E-08,    False,    [], 1.0, False],
    'id40': ['W3Jets0to100Z2   ',   5.15E-04,    False,    [], 1.0, False],
    'id41': ['W3Jets100to300Z2 ',   6.69E-04,    False,    [], 1.0, False],
    'id42': ['W3Jets300to800Z2 ',   1.97E-05,    False,    [], 1.0, False],
    'id43': ['W3Jets800to1600Z2',   1.81E-07,    False,    [], 1.0, False],
    'id44': ['W4Jets0to100Z2   ',   6.73E-04,    False,    [], 1.0, False],
    'id45': ['W4Jets100to300Z2 ',   1.13E-03,    False,    [], 1.0, False],
    'id46': ['W4Jets300to800Z2 ',   2.63E-05,    False,    [], 1.0, False],
    'id47': ['W4Jets800to1600Z2',   9.15E-08,    False,    [], 1.0, False],
    'id48': ['W5Jets0to100Z2   ',   1.03E-03,    False,    [], 1.0, False],
    'id49': ['W5Jets100to300Z2 ',   3.37E-04,    False,    [], 1.0, False],
    'id50': ['W5Jets300to800Z2 ',   1.44E-05,    False,    [], 1.0, False],
    'id51': ['W5Jets800to1600Z2',   6.24E-08,    False,    [], 1.0, False],
}


#samples haven't been skimevented
for key in wJetsSamples:
    wJetsSamples[key][3] += [ '%s'%x for x in glob(prefix+key+'/*.root') ]























