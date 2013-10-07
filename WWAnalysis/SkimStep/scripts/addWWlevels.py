#!/bin/env python



## Example:
#  ./addPUweight.py    --data puDATA.root --mc puMC.root --inputTreeFileName in.root -outputTreeFileName bubu.root --histo pileup --name puW2012 --treeName latino
#  ./addPUweight.py    -d     puDATA.root -m   puMC.root -i                  in.root -o                  bubu.root -h      pileup -n     puW2012 -t         latino


import optparse
import os
import sys
import ROOT
import numpy
import math

import re
from math import *
from pprint import pprint

wwBase=[
    'trigger==1.',
    'pfmet>20.',
    'mll>12',
    'zveto==1',
    'mpmet>20.',
    '(njet==0 || njet==1 || (dphilljetjet<pi/180.*165. || !sameflav )  )',
    'bveto_mu==1',
    'nextra==0',
    '(bveto_ip==1 &&  (njet != 1  || nbjet==0) && ((njet<2 || njet>3) || (jetbjpb1<=1.05 && jetbjpb2<=1.05)) )',
    'ptll>45.',
    '( !sameflav || ( (njet!=0 || dymva1>0.60) && (njet!=1 || dymva1>0.30) && ( njet==0 || njet==1 || (pfmet > 45.0)) ) )'
]


zerojet = 'njet == 0'
onejet  = 'njet == 1'
vbf     = '(njet >= 2 && njet <= 3 && (jetpt3 <= 30 || !(jetpt3 > 30 && (  (jeteta1-jeteta3 > 0 && jeteta2-jeteta3 < 0) || (jeteta2-jeteta3 > 0 && jeteta1-jeteta3 < 0))))) '

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

class WWSelAdder:
    def __init__(self):
        self.inputTreeFileName = ''
        self.outputTreeFileName = ''
        self.inFile = None
        self.outFile = None
        self.tree = None
        self.oldtree = None
        self.nentries = 0
        self.treeDir = ''

        self.branches = dict(
            wwsel=' && '.join(wwBase),
            wwsel0j=' && '.join(wwBase+[zerojet]),
            wwsel1j=' && '.join(wwBase+[onejet]),
            wwsel2j=' && '.join(wwBase+[vbf]),
        )

        
    def connect(self):
        filename = self.inputTreeFileName
        print 'Opening file: '+filename
        self.inFile = openTFile(filename)
        filename = self.outputTreeFileName
        print 'Output file name is : '+filename
        dir = os.path.dirname(filename)
        if not os.path.exists(dir) and dir:
            print 'directory does not exist yet...'
            print 'creating output directory: '+dir
            os.system('mkdir -p '+dir)
        print 'Creating output file: '+filename
        self.outFile = openTFile(filename,'recreate')

    def disconnect(self):
        self.inFile.Close()
        self.outFile.Write()
        self.outFile.Close()

###############################################################################################
##  _____ _                    _____             
## /  __ \ |                  |_   _|            
## | /  \/ | ___  _ __   ___    | |_ __ ___  ___ 
## | |   | |/ _ \| '_ \ / _ \   | | '__/ _ \/ _ \
## | \__/\ | (_) | | | |  __/   | | | |  __/  __/
##  \____/_|\___/|_| |_|\___|   \_/_|  \___|\___|
##                                              
    def makeClone(self):
        ## clone the tree
#         self.inFile.ls()
        oldTree = self.inFile.Get(self.treeDir)
#         print oldTree
        ## do not clone the branches which will be replaced
        ## i.e. set status to 0            
        oldbranches = [ b.GetName() for b in oldTree.GetListOfBranches() ]
        for branch in self.branches.iterkeys():
            if branch not in oldbranches: continue
            oldTree.SetBranchStatus(branch,0)

        newTree = oldTree.CloneTree(0)
        nentries = oldTree.GetEntries()
        print 'Tree with '+str(nentries)+' entries cloned...'
        self.nentries = nentries
        
        self.tree = newTree
        ## BUT keep all branches "active" in the old tree
        oldTree.SetBranchStatus('*'  ,1)

        self.oldtree = oldTree
#         for branch in self.tree.GetListOfBranches():
#             print branch



###############################################################################################
#  _       ___       _______      __           __  _              __     
# | |     / / |     / / ___/___  / /___  _____/ /_(_)____  ____  / /     
# | | /| / /| | /| / /\__ \/ _ \/ // _ \/ ___/ __/ // __ \/ __ \/ /      
# | |/ |/ / | |/ |/ /___/ /  __/ //  __/ /__/ /_/ // /_/ / / / /_/       
# |__/|__/  |__/|__//____/\___/_/ \___/\___/\__/_/ \____/_/ /_(_)        
#                                                                

    def add(self):

        formulas = [ (numpy.zeros(1, dtype = numpy.int32),ROOT.TTreeFormula(k,v,self.oldtree)) for k,v in self.branches.iteritems() ]
        print 'Adding WW selection flags'

        for flag,formula in formulas:
            name = formula.GetName()
            self.tree.Branch(name,flag,name+'/I')

        nentries = self.nentries
        print 'total number of entries: '+str(nentries)
        i = 0
        for ientry in xrange(0,nentries):
            i+=1
            self.oldtree.GetEntry(ientry)

            ## print event count
            step = 50000
            if i > 0 and i%step == 0.:
                print str(i)+' events processed.'

            import warnings

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                for flag,formula in formulas:
                    flag[0] = formula.EvalInstance()

            self.tree.Fill()



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
    
    parser.add_option('-i', '--inputTreeFileName' ,  dest='inputTreeFileName' ,  help='Name of the input  *.root file',)
    parser.add_option('-o', '--outputTreeFileName',  dest='outputTreeFileName',  help='Name of the output *.root file',)
    parser.add_option('-t', '--treeName'          ,  dest='treeDir'           ,  help='Name of the tree ',)


    (opt, args) = parser.parse_args()

    sys.argv.append('-b')
    ROOT.gROOT.SetBatch()

    if opt.inputTreeFileName is None:
        parser.error('No input file defined')
    if opt.outputTreeFileName is None:
        parser.error('No output file defined')


    adder = WWSelAdder()

    adder.inputTreeFileName  = opt.inputTreeFileName
    adder.outputTreeFileName = opt.outputTreeFileName
    adder.treeDir            = opt.treeDir


#    print s.systArgument

    adder.connect()
    adder.makeClone()
    adder.add()
    adder.disconnect()
    


    print 'Job finished...'


if __name__ == '__main__':
    main()



