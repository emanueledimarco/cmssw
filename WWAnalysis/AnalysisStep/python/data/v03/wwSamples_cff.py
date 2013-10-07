import FWCore.ParameterSet.Config as cms
from glob import glob

# Prefix where we expect idXX/*.root files to be
prefix = '/nfs/bluearc/group/trees/ww/V03_WWEventSkimProd_Newest/'

wwSamples = {
      #id    #name      #scale to 1pb^-1 #is data  #skim event files #enhancement #is signal
    'id01' : ['wz'      , 8.29E-06 ,   False,     [], 1.0, False],
    'id02' : ['zz'      , 2.81E-06 ,   False,     [], 1.0, False],
    'id03' : ['dyee'    , 8.33E-04 ,   False,     [], 1.0, False],
    'id04' : ['dymm'    , 8.33E-04 ,   False,     [], 1.0, False],
    'id05' : ['dytt'    , 8.35E-04 ,   False,     [], 1.0, False],
    'id06' : ['wm'      , 2.27E-03 ,   False,     [], 1.0, False],
    'id07' : ['we'      , 2.72E-03 ,   False,     [], 1.0, False],
    'id08' : ['wt'      , 2.62E-03 ,   False,     [], 1.0, False],
#     'id09' : ['ttbar'   , 1.66E-05 ,   False,     [], 1.0, False],
    'id10' : ['ts'      , 3.01E-06 ,   False,     [], 1.0, False],
    'id11' : ['tt'      , 4.32E-05 ,   False,     [], 1.0, False],
    'id12' : ['tw'      , 2.14E-05 ,   False,     [], 1.0, False],
#     'id13' : ['ww'      , 4.10E-05 ,   False,     [], 1.0, False],
#     'id14' : ['vv'      , 4.98E-06 ,   False,     [], 1.0, False],
    'id14b': ['vvFilt'  , 6.08E-06 ,   False,     [], 1.0, False],
#     'id15' : ['ttall'   , #DIV/0!  ,   False,     [], 1.0, False],
#     'id16' : ['wwNoPU'  , 4.10E-05 ,   False,     [], 1.0, False],
    'id17' : ['dyee10'  , 1.95E-03 ,   False,     [], 1.0, False],
    'id18' : ['dymm10'  , 1.96E-03 ,   False,     [], 1.0, False],
    'id19' : ['dytt10'  , 1.95E-03 ,   False,     [], 1.0, False],
                          
#     'id21' : ['phov'    , 1.50E-04 ,   False,     [], 1.0, False],
    'id22' : ['ggWW'    , 1.39E-06 ,   False,     [], 1.0, False],
    'id23' : ['ttjets'  , 1.39E-04 ,   False,     [], 1.0, False],
#     'id24' : ['dyToLL'  , 2.80E-04 ,   False,     [], 1.0, False],
}

wwGroupings = { 'ww':   ['vvFilt','ggWW'],
                'diBoson':['wz','zz'],
                'top':    ['ttjets','ts','tw','tt'],
                'zJets':  ['dyee','dymm','dytt','dyee10','dymm10','dytt10'],
                'wJets':  ['we','wm','wt'] }

# wwGroupings = { 'qqWW':   ['vvFilt'],
#                 'ggWW':   ['ggWW'],
#                 'diBoson':['wz','zz'],
#                 'top':    ['ttjets','ts','tw','tt'],
#                 'zJets':  ['dyee','dymm','dytt','dyee10','dymm10','dytt10'],
#                 'wJets':  ['we','wm','wt'] }

#loop over the looper
for key in wwSamples:
    wwSamples[key][3] += [ '%s'%x for x in glob(prefix+key+'/*.root') ]


