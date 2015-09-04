#!/usr/bin/python
import sys
import os
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-v", "--vtx", dest="vtxbinning", default="false",
                  help="do vtx,eta bins instead of default pt,eta bins")
(options, args) = parser.parse_args()

inputfdata="zeetree.root"
outputpref="Z2012_"

from ROOT import gROOT
gROOT.LoadMacro('FitZMassScaleAndResolution.C+')
from ROOT import FitZMassScaleAndResolution

if not options.vtxbinning:
    # do the fits in eta,pt bins of the electrons
    for ptb in range(0,4):
        print "pt bin = ",ptb
        for etab in range(0,2):
            print "   eta bin = ",etab
            for p4kind in range(0,4):
                print " local reco type = ",p4kind
                outputfdata="data"+outputpref+"PtBin"+str(ptb)+"_EtaBin"+str(etab)+"_Reco"+str(p4kind)
                print "        fitting DATA..."
                FitZMassScaleAndResolution(p4kind,inputfdata,outputfdata,ptb,etab,-1)
                print "        output is in file "+outputfdata
else:
    # do the fits in eta, nvtx bins 
    for vtxb in range(0,4):
        print "nvtx bin = ",vtxb
        for etab in range(0,2):
            print "   eta bin = ",etab
            for p4kind in range(0,4):
                print " local reco type = ",p4kind
                outputfdata="data"+outputpref+"VtxBin"+str(vtxb)+"_EtaBin"+str(etab)+"_Reco"+str(p4kind)
                print "        fitting DATA..."
                FitZMassScaleAndResolution(p4kind,inputfdata,outputfdata,-1,etab,vtxb)
                print "        output is in file "+outputfdata
