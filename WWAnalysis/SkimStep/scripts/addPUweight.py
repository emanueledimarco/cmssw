#!/usr/bin/env python

## Example:
#  ./addPUweight.py    --data puDATA.root --mc puMC.root --inputTreeFileName in.root -outputTreeFileName bubu.root --histo pileup --name puW2012 --treeName latino --kindPU  trpu
#  ./addPUweight.py    -d     puDATA.root -m   puMC.root -i                  in.root -o                  bubu.root -h      pileup -n     puW2012 -t         latino --k       trpu


import optparse
import os
import ROOT
import numpy
from ROOT import *
from ROOT import std
import math

import re
from math import *
from pprint import pprint


#
#     
#        \  |   ___|         \         |      |       _ \   |   |     ___ \    _ \ _ | ___ \  
#       |\/ |  |            _ \     _` |   _` |      |   |  |   |        ) |  |   |  |    ) | 
#       |   |  |           ___ \   (   |  (   |      ___/   |   |       __/   |   |  |   __/  
#      _|  _| \____|     _/    _\ \__,_| \__,_|     _|     \___/      _____| \___/  _| _____| 
#
#                                                                                             



def openTFile(path, option=''):
    f =  ROOT.TFile.Open(path,option)
    if not f.__nonzero__() or not f.IsOpen():
        raise NameError('File '+path+' not open')
    return f

def getHist(file, hist):
    h = file.Get(hist)
    if not h.__nonzero__():
        raise NameError('Histogram '+hist+' doesn\'t exist in '+str(file))
    return h




###############################################################################################
###############################################################################################
###############################################################################################

class addPUweight:
    def __init__(self):
        self.inputDATAFileName = ''
        self.inputMCFileName = ''
        self.inputTreeFileName = ''
        self.outputTreeFileName = ''
        self.inFile = None
        self.outFile = None
        self.inDATAFile = None
        self.inMCFile   = None
        self.puScaleDATAhisto = None
        self.puScaleMChisto   = None
        self.HistoName = ''
        self.ttree = None
        self.oldttree = None
        self.nentries = 0
        self.treeDir = ''
        self.weightName = ''
        self.kindPU = ''


        
    def __del__(self):
        if self.outFile:
            self.outFile.Write()
            self.outFile.Close()



    def openOriginalTFile(self):
        filename = self.inputTreeFileName
        print 'opening file: '+filename
        self.inFile = ROOT.TFile.Open(filename)
        filename = self.inputDATAFileName
        print 'opening file: '+filename
        self.inDATAFile = ROOT.TFile.Open(filename)
        filename = self.inputMCFileName
        print 'opening file: '+filename
        self.inMCFile = ROOT.TFile.Open(filename)

    def openOutputTFile(self):
        filename = self.outputTreeFileName
        print 'output file name is : '+filename
        dir = os.path.dirname(filename)
        if not os.path.exists(dir):
            print 'directory does not exist yet...'
            print 'creating output directory: '+dir
            os.system('mkdir -p '+dir)
        print 'creating output file: '+filename
        self.outFile = ROOT.TFile.Open(filename,'recreate')

###############################################################################################
##  _____ _                    _____             
## /  __ \ |                  |_   _|            
## | /  \/ | ___  _ __   ___    | |_ __ ___  ___ 
## | |   | |/ _ \| '_ \ / _ \   | | '__/ _ \/ _ \
## | \__/\ | (_) | | | |  __/   | | | |  __/  __/
##  \____/_|\___/|_| |_|\___|   \_/_|  \___|\___|
##                                              
    def cloneTree(self):
        ## clone the tree
        self.inFile.ls()
        oldTree = self.inFile.Get(self.treeDir)
        ## do not clone the branches which should be scaled
        ## i.e. set status to 0            
        if oldTree.GetBranch(self.weightName) :
            oldTree.SetBranchStatus(self.weightName,0)

        newTree = oldTree.CloneTree(0)
        nentries = oldTree.GetEntriesFast()
        print 'Tree with '+str(nentries)+' entries cloned...'
        self.nentries = nentries
        
        self.ttree = newTree
        ## BUT keep all branches "active" in the old tree
        oldTree.SetBranchStatus('*'  ,1)

        self.oldttree = oldTree
        for branch in self.ttree.GetListOfBranches():
            print branch



###############################################################################################
##
##       _ \   |   |                     _)         |      |   
##      |   |  |   |     \ \  \   /  _ \  |   _` |  __ \   __| 
##      ___/   |   |      \ \  \ /   __/  |  (   |  | | |  |   
##     _|     \___/        \_/\_/  \___| _| \__, | _| |_| \__| 
##                                          |___/              
##
##

    def puWeight(self):
        print 'PU scaling'
        puW = numpy.ones(1, dtype=numpy.float32)
        self.ttree.Branch(self.weightName,puW,self.weightName+"/F")

        self.inDATAFile = openTFile(self.inputDATAFileName)
        self.inMCFile   = openTFile(self.inputMCFileName)

        self.puScaleDATAhisto = getHist(self.inDATAFile,self.HistoName)
        self.puScaleMChisto   = getHist(self.inMCFile,self.HistoName)


        data_nBin = self.puScaleDATAhisto.GetNbinsX()
        data_minValue = self.puScaleDATAhisto.GetXaxis().GetXmin()
        data_maxValue = self.puScaleDATAhisto.GetXaxis().GetXmax()
        data_dValue = (data_maxValue - data_minValue) / data_nBin

        mc_nBin = self.puScaleMChisto.GetNbinsX()
        mc_minValue = self.puScaleMChisto.GetXaxis().GetXmin()
        mc_maxValue = self.puScaleMChisto.GetXaxis().GetXmax()
        mc_dValue = (mc_maxValue - mc_minValue) / mc_nBin
  
        ratio = mc_dValue/data_dValue
        nBin = data_nBin
        minValue = data_minValue
        maxValue = data_maxValue
        dValue = data_dValue
 
        print " ratio = "
        print ratio
 
        if (mc_dValue/data_dValue - (int) (mc_dValue/data_dValue)) != 0 :
          print " ERROR:: incompatible intervals!"
          exit(0);
        
 
        puScaleDATA   = std.vector(float)()
        puScaleMCtemp = std.vector(float)()
        puScaleMC     = std.vector(float)()

        ################
        # remove last bin -> peak in DATA distribution
        # nBin = nBin-1
        ################   

        for iBin in range(0, nBin):
            puScaleDATA.push_back(self.puScaleDATAhisto.GetBinContent(iBin+1))
            mcbin = int(floor(iBin / ratio))
            puScaleMCtemp.push_back(self.puScaleMChisto.GetBinContent(mcbin+1))

 
        integralDATA = 0.
        integralMC   = 0.
 
        for iBin in range(0, nBin):
            integralDATA += puScaleDATA.at(iBin)
            integralMC   += puScaleMCtemp.at(iBin)

        print " integralDATA = " + "%.3f" %integralDATA
        print " integralMC   = " + "%.3f" %integralMC
 
        for iBin in range(0, nBin):
            puScaleMC.push_back( puScaleMCtemp.at(iBin) * integralDATA / integralMC) 


 
        nentries = self.nentries
        print 'total number of entries: '+str(nentries)
        i = 0
        for ientry in range(0,nentries):
            i+=1
            self.oldttree.GetEntry(ientry)

            ## print event count
            step = 50000
            if i > 0 and i%step == 0.:
                print str(i)+' events processed.'

##
## ----------------------------------------------------------------
##
            puW[0] = 1.
            
            ibin = -1

## true pu reweighting 
            if self.kindPU == "trpu" : 
               ibin = int(self.oldttree.trpu / dValue)
## in time pu reweighting (observed)
            if self.kindPU == "itpu" : 
               ibin = int(self.oldttree.itpu / dValue)

            if ibin < puScaleDATA.size() :
               if puScaleMC.at(ibin) != 0 :
                  puW[0] = 1. * puScaleDATA.at(ibin) / puScaleMC.at(ibin)
               else :
                  puW[0] = 1.
            else :
               ibin = puScaleDATA.size()-1
               if puScaleMC.at(ibin) != 0 :
                  puW[0] = 1. * puScaleDATA.at(ibin) / puScaleMC.at(ibin)
               else :
                  puW[0] = 1.

#            print puW

            self.ttree.Fill()


###############################################################################################
##                  _       
##                 (_)      
##  _ __ ___   __ _ _ _ __  
## | '_ ` _ \ / _` | | '_ \ 
## | | | | | | (_| | | | | |
## |_| |_| |_|\__,_|_|_| |_|
##                         

def main():
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    
    parser.set_defaults(overwrite=False)
    
    parser.add_option('-d', '--data'     ,   dest='inputDATAFileName',   help='Name of the input *.root file with pu histo from data',)
    parser.add_option('-m', '--mc'       ,   dest='inputMCFileName'  ,   help='Name of the input *.root file with pu histo from mc',)
    parser.add_option('-H', '--histoName',   dest='HistoName'        ,   help='Histogram name',)

    parser.add_option('-i', '--inputTreeFileName' ,  dest='inputTreeFileName' ,  help='Name of the input  *.root file',)
    parser.add_option('-o', '--outputTreeFileName',  dest='outputTreeFileName',  help='Name of the output *.root file',)
    parser.add_option('-t', '--treeName'          ,  dest='treeDir'           ,  help='Name of the tree ',)
    parser.add_option('-k', '--kindPU'            ,  dest='kindPU'            ,  help='kind of PU reweighting: trpu (= true pu), itpu (= in time pu, that is observed!)',)




    parser.add_option('-n', '--name' ,  dest='weightName' ,  help='Name of the branch of PU weight',)

#    
#    --data puDATA.root --mc puMC.root --inputTreeFileName in.root -outputTreeFileName bubu.root --histo pileup --name puW2012 --treeName latino --kindPU  trpu
#    -d     puDATA.root -m   puMC.root -i                  in.root -o                  bubu.root -h      pileup -n     puW2012 -t         latino --k       trpu
#    

    (opt, args) = parser.parse_args()

    if opt.inputDATAFileName is None:
        parser.error('No input DATA pu file defined')
    if opt.inputMCFileName is None:
        parser.error('No input MC pu file defined')
    if opt.inputTreeFileName is None:
        parser.error('No input file defined')
    if opt.outputTreeFileName is None:
        parser.error('No output file defined')
    if opt.HistoName is None:
        parser.error('No Histogram name defined')
    if opt.weightName is None:
        parser.error('No Weigth branch name defined')

##     sys.argv.append('-b')
##     ROOT.gROOT.SetBatch()

    w = addPUweight()
    w.inputDATAFileName = opt.inputDATAFileName
    w.inputMCFileName = opt.inputMCFileName
    w.HistoName = opt.HistoName

    w.inputTreeFileName  = opt.inputTreeFileName
    w.outputTreeFileName = opt.outputTreeFileName
    w.treeDir            = opt.treeDir
    w.weightName = opt.weightName
    w.kindPU = opt.kindPU



#    print s.systArgument

    w.openOriginalTFile()
    w.openOutputTFile()
    w.cloneTree()
    
    w.puWeight()
    


    print 'Job finished...'


if __name__ == '__main__':
    main()



