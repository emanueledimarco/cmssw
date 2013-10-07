#!/usr/bin/python
import sys
import os
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-r", "--recompile",
                  action="store_true", dest="recompile", default=False,
                  help="recompile the fake rate estimator")
parser.add_option("-d", "--drawonly",
                  action="store_true", dest="drawonly", default=False,
                  help="do not run. Only draw the results based on existing ROOT files")

(options, args) = parser.parse_args()

if(options.drawonly == False):
    if(options.recompile):
        print 'recompiling...'
        os.system('make clean')
        os.system('make fakerate')
        print 'done complining.'

os.system('./fakerate fr')
os.system('hadd -f fakerates_zll1e.root fr_zlle-ElectronMisid*')
os.system('hadd -f fakerates_zll1m.root fr_zllm-MuonMisid*')
os.system('rm -f fr_zlle-ElectronMisid*')
os.system('rm -f fr_zllm-MuonMisid*')        

from ROOT import gROOT
gROOT.LoadMacro('src/drawFR.cc+')

# compute electron FR
from ROOT import drawElectronIdsUnbiased
drawElectronIdsUnbiased()
from ROOT import fitElectronIdsUnbiased
fitElectronIdsUnbiased(0)
fitElectronIdsUnbiased(1)
fitElectronIdsUnbiased(2)

# compute muon FR
from ROOT import drawMuonIdsUnbiased
drawMuonIdsUnbiased()
print "And now fit with thr paramteric function..."
from ROOT import fitMuonIdsUnbiased
fitMuonIdsUnbiased(0)
fitMuonIdsUnbiased(1)

# make the 2D map and save in a ROOT file
gROOT.LoadMacro('src/doFRMaps.cc+')
from ROOT import doEleFRMapsHzz4l
doEleFRMapsHzz4l(0)
doEleFRMapsHzz4l(1)
from ROOT import doMuFRMapsHzz4l
doMuFRMapsHzz4l(0)
doMuFRMapsHzz4l(1)

os.system('mv elfr*root src')
os.system('mv mufr*root src')
