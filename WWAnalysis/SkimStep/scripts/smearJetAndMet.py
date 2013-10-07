#!/usr/bin/env python

## Example:
#  ./smearJetAndMet.py    --inputTreeFileName in.root -outputTreeFileName bubu.root  --kindMCDATA  MC2012
#  ./smearJetAndMet.py     -i                  in.root -o                  bubu.root --k           MC2012


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
#                                                        |        |                          |       \  |  ____| __ __| 
#             __|  __ `__ \    _ \   _` |   __|          |   _ \  __|       _` |  __ \    _` |      |\/ |  __|      |   
#           \__ \  |   |   |   __/  (   |  |         \   |   __/  |        (   |  |   |  (   |      |   |  |        |   
#           ____/ _|  _|  _| \___| \__,_| _|        \___/  \___| \__|     \__,_| _|  _| \__,_|     _|  _| _____|   _|   
#                                               
#           
#
#                                                                                             

                                                               







def deltaPhi(l1,l2):
    dphi = fabs(l1.DeltaPhi(l2))
#    dphi = fabs(ROOT.Math.VectorUtil.DeltaPhi(l1.p4(),l2.p4()))
    return dphi


# pfjet resolutions. taken from AN-2010-371
def ErrEt(Et,Eta) :

    InvPerr2 = 0.0

    N = 0.0
    S = 0.0
    C = 0.0
    m = 0.0

    if fabs(Eta) < 0.5 :
        N = 3.96859
        S = 0.18348
        C = 0.
        m = 0.62627
    elif fabs(Eta) < 1.0 :
        N = 3.55226
        S = 0.24026
        C = 0.
        m = 0.52571
    elif fabs(Eta) < 1.5 :
        N = 4.54826
        S = 0.22652
        C = 0.
        m = 0.58963
    elif fabs(Eta) < 2.0 :
        N = 4.62622
        S = 0.23664
        C = 0.
        m = 0.48738
    elif fabs(Eta) < 3.0 :
        N = 2.53324
        S = 0.34306
        C = 0.
        m = 0.28662
    elif fabs(Eta) < 5.0 :
        N = 2.95397
        S = 0.11619
        C = 0.
        m = 0.96086

    # this is the absolute resolution (squared), not sigma(pt)/pt
    # so have to multiply by pt^2, thats why m+1 instead of m-1
    InvPerr2 =  (N * fabs(N) ) + (S * S) * pow(Et, m+1) + (C * C) * Et * Et

    return sqrt(InvPerr2)/Et;




def openTFile(path, option=''):
    f =  ROOT.TFile.Open(path,option)
    if not f.__nonzero__() or not f.IsOpen():
        raise NameError('File '+path+' not open')
    return f

def deltaPhi(l1,l2):
    dphi = fabs(l1.DeltaPhi(l2))
#    dphi = fabs(ROOT.Math.VectorUtil.DeltaPhi(l1.p4(),l2.p4()))
    return dphi



###############################################################################################
###############################################################################################
###############################################################################################

class JetAndMetSmear:
    def __init__(self):
        self.inputTreeFileName = ''
        self.outputTreeFileName = ''
        self.inFile = None
        self.outFile = None
        self.HistoName = ''
        self.ttree = None
        self.oldttree = None
        self.nentries = 0
        self.treeDir = ''
        self.kindMCDATA = ''  # MC2011, MC2012, DATA2011, DATA2012
        self.getDYMVAV0j0 = None
        self.getDYMVAV0j1 = None
        self.getDYMVAV0j0 = None
        self.getDYMVAV1j1 = None
        self.var1 = numpy.ones(1, dtype=numpy.float32)
        self.var2 = numpy.ones(1, dtype=numpy.float32)
        self.var3 = numpy.ones(1, dtype=numpy.float32)
        self.var4 = numpy.ones(1, dtype=numpy.float32)
        self.var5 = numpy.ones(1, dtype=numpy.float32)
        self.var6 = numpy.ones(1, dtype=numpy.float32)
        self.var7 = numpy.ones(1, dtype=numpy.float32)
        self.var8 = numpy.ones(1, dtype=numpy.float32)
        self.var9 = numpy.ones(1, dtype=numpy.float32)
        self.var10 = numpy.ones(1, dtype=numpy.float32)
        self.var11 = numpy.ones(1, dtype=numpy.float32)
        self.var12 = numpy.ones(1, dtype=numpy.float32)
        self.var13 = numpy.ones(1, dtype=numpy.float32)
        self.var14 = numpy.ones(1, dtype=numpy.float32)
        self.var15 = numpy.ones(1, dtype=numpy.float32)
        self.var16 = numpy.ones(1, dtype=numpy.float32)
        self.var17 = numpy.ones(1, dtype=numpy.float32)
        self.var18 = numpy.ones(1, dtype=numpy.float32)
        self.var19 = numpy.ones(1, dtype=numpy.float32)
        self.var20 = numpy.ones(1, dtype=numpy.float32)
        self.var21 = numpy.ones(1, dtype=numpy.float32)
        self.var22 = numpy.ones(1, dtype=numpy.float32)

        
    def __del__(self):
        if self.outFile:
            self.outFile.Write()
            self.outFile.Close()



    def openOriginalTFile(self):
        filename = self.inputTreeFileName
        print 'opening file: '+filename
        self.inFile = ROOT.TFile.Open(filename)

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


    def createDYMVA(self):
        self.getDYMVAV0j0 = TMVA.Reader();
        self.getDYMVAV0j1 = TMVA.Reader();
        self.getDYMVAV1j0 = TMVA.Reader();
        self.getDYMVAV1j1 = TMVA.Reader();

        self.getDYMVAV0j0.AddVariable("met",           (self.var1))
        self.getDYMVAV0j0.AddVariable("trackMet",      (self.var2))
        self.getDYMVAV0j0.AddVariable("jet1pt",        (self.var3))
        self.getDYMVAV0j0.AddVariable("metSig",        (self.var4))
        self.getDYMVAV0j0.AddVariable("dPhiDiLepJet1", (self.var5))
        self.getDYMVAV0j0.AddVariable("dPhiJet1MET",   (self.var6))
        self.getDYMVAV0j0.AddVariable("mt",            (self.var7))

        self.getDYMVAV0j1.AddVariable("met",           (self.var1))
        self.getDYMVAV0j1.AddVariable("trackMet",      (self.var2))
        self.getDYMVAV0j1.AddVariable("jet1pt",        (self.var3))
        self.getDYMVAV0j1.AddVariable("metSig",        (self.var4))
        self.getDYMVAV0j1.AddVariable("dPhiDiLepJet1", (self.var5))
        self.getDYMVAV0j1.AddVariable("dPhiJet1MET",   (self.var6))
        self.getDYMVAV0j1.AddVariable("mt",            (self.var7))

        self.getDYMVAV1j0.AddVariable("pmet",           (self.var1))
        self.getDYMVAV1j0.AddVariable("pTrackMet",      (self.var2))
        self.getDYMVAV1j0.AddVariable("nvtx",           (self.var3))
        self.getDYMVAV1j0.AddVariable("dilpt",          (self.var4))
        self.getDYMVAV1j0.AddVariable("jet1pt",         (self.var5))
        self.getDYMVAV1j0.AddVariable("metSig",         (self.var6))
        self.getDYMVAV1j0.AddVariable("dPhiDiLepJet1",  (self.var7))
        self.getDYMVAV1j0.AddVariable("dPhiDiLepMET",   (self.var8))
        self.getDYMVAV1j0.AddVariable("dPhiJet1MET",    (self.var9))
        self.getDYMVAV1j0.AddVariable("recoil",         (self.var10))
        self.getDYMVAV1j0.AddVariable("mt",             (self.var11))

        self.getDYMVAV1j1.AddVariable("pmet",           (self.var1))
        self.getDYMVAV1j1.AddVariable("pTrackMet",      (self.var2))
        self.getDYMVAV1j1.AddVariable("nvtx",           (self.var3))
        self.getDYMVAV1j1.AddVariable("dilpt",          (self.var4))
        self.getDYMVAV1j1.AddVariable("jet1pt",         (self.var5))
        self.getDYMVAV1j1.AddVariable("metSig",         (self.var6))
        self.getDYMVAV1j1.AddVariable("dPhiDiLepJet1",  (self.var7))
        self.getDYMVAV1j1.AddVariable("dPhiDiLepMET",   (self.var8))
        self.getDYMVAV1j1.AddVariable("dPhiJet1MET",    (self.var9))
        self.getDYMVAV1j1.AddVariable("recoil",         (self.var10))
        self.getDYMVAV1j1.AddVariable("mt",             (self.var11))


        baseCMSSW = string(getenv("CMSSW_BASE"))
        self.getDYMVAV0j0.BookMVA("BDTB",baseCMSSW+string("/src/DYMvaInCMSSW/GetDYMVA/data/TMVA_0j_BDTB.weights.xml"))
        self.getDYMVAV0j1.BookMVA("BDTB",baseCMSSW+string("/src/DYMvaInCMSSW/GetDYMVA/data/TMVA_1j_BDTB.weights.xml"))
        self.getDYMVAV1j0.BookMVA("BDTG",baseCMSSW+string("/src/DYMvaInCMSSW/GetDYMVA/data/TMVA_BDTG_0j_MCtrain.weights.xml"))
        self.getDYMVAV1j1.BookMVA("BDTG",baseCMSSW+string("/src/DYMvaInCMSSW/GetDYMVA/data/TMVA_BDTG_1j_MCtrain.weights.xml"))





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
        oldTree.SetBranchStatus("dphillmet",0)
        oldTree.SetBranchStatus("dphilmet",0)
        oldTree.SetBranchStatus("dphilmet1",0)
        oldTree.SetBranchStatus("dphilmet2",0)
        oldTree.SetBranchStatus("mth",0)
        oldTree.SetBranchStatus("mtw1",0)
        oldTree.SetBranchStatus("mtw2",0)
        oldTree.SetBranchStatus("pfmet",0)
        oldTree.SetBranchStatus("pfmetphi",0)
        oldTree.SetBranchStatus("ppfmet",0)
        oldTree.SetBranchStatus("mpmet",0)
#        oldTree.SetBranchStatus("dymva0",0)
#        oldTree.SetBranchStatus("dymva1",0)
                     
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
##          
##             \  |  ____| __ __|                   
##            |\/ |  __|      |      
##            |   |  |        |      
##           _|  _| _____|   _|      
##                                   
##          
##
##
##

    def applyMetSmearing(self):
        print 'MET Jet smearing'

        dphillmet = numpy.ones(1, dtype=numpy.float32)
        dphilmet  = numpy.ones(1, dtype=numpy.float32)
        dphilmet1 = numpy.ones(1, dtype=numpy.float32)
        dphilmet2 = numpy.ones(1, dtype=numpy.float32)
        mth       = numpy.ones(1, dtype=numpy.float32)
        mtw1      = numpy.ones(1, dtype=numpy.float32)
        mtw2      = numpy.ones(1, dtype=numpy.float32)
        pfmet     = numpy.ones(1, dtype=numpy.float32)
        pfmetphi  = numpy.ones(1, dtype=numpy.float32)
        ppfmet    = numpy.ones(1, dtype=numpy.float32)
        mpmet     = numpy.ones(1, dtype=numpy.float32)
        dymva0    = numpy.ones(1, dtype=numpy.float32)
        dymva1    = numpy.ones(1, dtype=numpy.float32)

        self.ttree.Branch("dphillmet",dphillmet,"dphillmet/F")
        self.ttree.Branch("dphilmet" ,dphilmet ,"dphilmet/F")
        self.ttree.Branch("dphilmet1",dphilmet1,"dphilmet1/F")
        self.ttree.Branch("dphilmet2",dphilmet2,"dphilmet2/F")
        self.ttree.Branch("mth"      ,mth      ,"mth/F")
        self.ttree.Branch("mtw1"     ,mtw1     ,"mtw1/F")
        self.ttree.Branch("mtw2"     ,mtw2     ,"mtw2/F")
        self.ttree.Branch("pfmet"    ,pfmet    ,"pfmet/F")
        self.ttree.Branch("pfmetphi" ,pfmetphi ,"pfmetphi/F")
        self.ttree.Branch("ppfmet"   ,ppfmet   ,"ppfmet/F")
        self.ttree.Branch("mpmet"    ,mpmet    ,"mpmet/F")
#        self.ttree.Branch("dymva0"   ,dymva0   ,"dymva0/F")
#        self.ttree.Branch("dymva1"   ,dymva1   ,"dymva1/F")


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





            metx = self.oldttree.pfmet * cos (self.oldttree.pfmetphi)
            mety = self.oldttree.pfmet * sin (self.oldttree.pfmetphi)

            dmetx = 0.0
            dmety = 0.0

#            if self.kindMCDATA == "DATA2011" :
#            if self.kindMCDATA == "MC2011" :
#            if self.kindMCDATA == "DATA2012" :
#            if self.kindMCDATA == "MC2012" :
            if self.kindMCDATA == "NOCHANGE" :
                dmetx = 0
                dmety = 0


            #
            ##
            ###
            ####
            ####
            
            dpxsum = 0
            dpysum = 0
            
            if self.oldttree.jetpt1>0 :
                oldpt  = self.oldttree.jetpt1
                oldphi = self.oldttree.jetphi1
                oldeta = self.oldttree.jeteta1
#                oldp   = oldpt / sin (theta)
#                theta = 2.0 * atan ( exp (- oldeta))

                smear = 1.00000

                if  (fabs(oldeta)<=0.5) :
                    smear = 1.052  
                if  ((fabs(oldeta)<=1.1) and (fabs(oldeta)>0.5)) :
                    smear = 1.057  
                if  ((fabs(oldeta)<=1.7) and (fabs(oldeta)>1.1)) :
                    smear = 1.096
                if  ((fabs(oldeta)<=2.3) and (fabs(oldeta)>1.7)) :
                    smear = 1.134
                if  (fabs(oldeta)>2.3) :
                    smear = 1.288

                sigma = ErrEt(oldpt,oldeta)

#                print sigma
#                addres = (sigma * (smear - 1.))
#                print addres
#                print " ~~~~ "
#
#                deltap  = gRandom.Gaus(oldp,  (smear - 1.) * oldp )
#                deltapt = gRandom.Gaus(oldpt, sigma * (smear - 1.) * oldpt)
                deltapt = gRandom.Gaus(oldpt, sigma * sqrt(smear*smear - 1.) * oldpt)
                deltapt = deltapt - oldpt

                dpx = deltapt * cos (oldphi)
                dpy = deltapt * sin (oldphi)

                dpxsum = dpxsum + dpx
                dpysum = dpysum + dpy

#           see if second jet
            if self.oldttree.jetpt2>0 :
                oldpt  = self.oldttree.jetpt2
                oldphi = self.oldttree.jetphi2
                oldeta = self.oldttree.jeteta2

                smear = 1.00000

                if  (fabs(oldeta)<=0.5) :
                    smear = 1.052  
                if  ((fabs(oldeta)<=1.1) and (fabs(oldeta)>0.5)) :
                    smear = 1.057  
                if  ((fabs(oldeta)<=1.7) and (fabs(oldeta)>1.1)) :
                    smear = 1.096
                if  ((fabs(oldeta)<=2.3) and (fabs(oldeta)>1.7)) :
                    smear = 1.134
                if  (fabs(oldeta)>2.3) :
                    smear = 1.288

                sigma = ErrEt(oldpt,oldeta)

#                deltap  = gRandom.Gaus(oldp,  (smear - 1.) * oldp )
                #deltapt = gRandom.Gaus(oldpt, sigma * (smear - 1.) * oldpt)
                deltapt = gRandom.Gaus(oldpt, sigma * sqrt(smear*smear - 1.) * oldpt)
                deltapt = deltapt - oldpt

                dpx = deltapt * cos (oldphi)
                dpy = deltapt * sin (oldphi)

                dpxsum = dpxsum + dpx
                dpysum = dpysum + dpy


#          see if third jet
            if self.oldttree.jetpt3>0 :
                oldpt  = self.oldttree.jetpt3
                oldphi = self.oldttree.jetphi3
                oldeta = self.oldttree.jeteta3

                smear = 1.00000

                if  (fabs(oldeta)<=0.5) :
                    smear = 1.052  
                if  ((fabs(oldeta)<=1.1) and (fabs(oldeta)>0.5)) :
                    smear = 1.057  
                if  ((fabs(oldeta)<=1.7) and (fabs(oldeta)>1.1)) :
                    smear = 1.096
                if  ((fabs(oldeta)<=2.3) and (fabs(oldeta)>1.7)) :
                    smear = 1.134
                if  (fabs(oldeta)>2.3) :
                    smear = 1.288

                sigma = ErrEt(oldpt,oldeta)

#                deltap  = gRandom.Gaus(oldp,  (smear - 1.) * oldp )
                #deltapt = gRandom.Gaus(oldpt, sigma * (smear - 1.) * oldpt)
                deltapt = gRandom.Gaus(oldpt, sigma * sqrt(smear*smear - 1.) * oldpt)
                deltapt = deltapt - oldpt

                dpx = deltapt * cos (oldphi)
                dpy = deltapt * sin (oldphi)

                dpxsum = dpxsum + dpx
                dpysum = dpysum + dpy

#          see if fourth jet
            if self.oldttree.jetpt4>0 :
                oldpt  = self.oldttree.jetpt4
                oldphi = self.oldttree.jetphi4
                oldeta = self.oldttree.jeteta4

                smear = 1.00000

                if  (fabs(oldeta)<=0.5) :
                    smear = 1.052  
                if  ((fabs(oldeta)<=1.1) and (fabs(oldeta)>0.5)) :
                    smear = 1.057  
                if  ((fabs(oldeta)<=1.7) and (fabs(oldeta)>1.1)) :
                    smear = 1.096
                if  ((fabs(oldeta)<=2.3) and (fabs(oldeta)>1.7)) :
                    smear = 1.134
                if  (fabs(oldeta)>2.3) :
                    smear = 1.288

                sigma = ErrEt(oldpt,oldeta)

#                deltap  = gRandom.Gaus(oldp,  (smear - 1.) * oldp )
                #deltapt = gRandom.Gaus(oldpt, sigma * (smear - 1.) * oldpt)
                deltapt = gRandom.Gaus(oldpt, sigma * sqrt(smear*smear - 1.) * oldpt)
                deltapt = deltapt - oldpt

                dpx = deltapt * cos (oldphi)
                dpy = deltapt * sin (oldphi)

                dpxsum = dpxsum + dpx
                dpysum = dpysum + dpy


            newmetx = metx - dpxsum
            newmety = mety - dpysum

            ####
            ####
            ###
            ##
            #

        
            newmet = sqrt(newmetx*newmetx + newmety*newmety)

            met = ROOT.TLorentzVector()
            met.SetPxPyPzE(newmetx, newmety, 0, newmet)           
        
            pfmet    [0] = met.Pt()
            pfmetphi[0] = met.Phi()


            l1 = ROOT.TLorentzVector()
            l2 = ROOT.TLorentzVector()
            l1.SetPtEtaPhiM(self.oldttree.pt1, self.oldttree.eta1, self.oldttree.phi1, 0)
            l2.SetPtEtaPhiM(self.oldttree.pt2, self.oldttree.eta2, self.oldttree.phi2, 0)

            dphillmet[0] = deltaPhi(l1+l2,met)
            dphilmet1[0] = deltaPhi(l1   ,met)
            dphilmet2[0] = deltaPhi(   l2,met)
            dphilmet [0] = min (dphilmet1[0], dphilmet2[0])

            mth      [0] = sqrt( 2 * self.oldttree.ptll * pfmet[0] * ( 1 - cos(dphillmet[0]) ) );
            mtw1     [0] = sqrt( 2 * self.oldttree.pt1  * pfmet[0] * ( 1 - cos(dphilmet1[0]) ) );
            mtw2     [0] = sqrt( 2 * self.oldttree.pt2  * pfmet[0] * ( 1 - cos(dphilmet2[0]) ) );


            if dphilmet [0] < math.pi/2.:
               ppfmet   [0] = pfmet[0] * sin(dphilmet [0])
            else :
               ppfmet   [0] = pfmet[0]
 

            mpmet[0] = min(ppfmet[0],self.oldttree.pchmet)


            #### correct dymva ####


#            dphilljet1  = numpy.ones(1, dtype=numpy.float32)
#            dphimetjet1 = numpy.ones(1, dtype=numpy.float32)
#
#            jetpt1 = self.oldttree.jetpt1
#            if self.oldttree.jetpt1 < 15 :
#                jetpt1 = 15
#                dphilljet1[0]  = -0.1;
#                dphimetjet1[0] = -0.1;
#            else :
#                jet1 = ROOT.TLorentzVector()
#                jet1.SetPtEtaPhiM(self.oldttree.jetpt1, self.oldttree.jeteta1, self.oldttree.jetphi1, 0)
#                dphilljet1[0]  = deltaPhi(l1+l2,jet1)
#                dphimetjet1[0] = deltaPhi(met,jet1)
#
#            self.var1[0] =  pfmet[0] 
#            self.var2[0] =  self.oldttree.chmet
#            self.var3[0] =  jetpt1
#            self.var4[0] =  self.oldttree.pfmetSignificance #pfMetSignificance
#            self.var5[0] =  dphilljet1[0]  
#            self.var6[0] =  dphimetjet1[0]
#            self.var7[0] =  mth[0]
# 
#
#            if self.oldttree.njet == 0:
#                  dymva0[0] = self.getDYMVAV0j0.EvaluateMVA("BDTB")
#            elif self.oldttree.njet == 1:
#                  dymva0[0] = self.getDYMVAV0j1.EvaluateMVA("BDTB")
#            else :
#                  dymva0[0] = -999
#
#            recoil = numpy.ones(1, dtype=numpy.float32)
#            px_rec = pfmet[0] * cos(pfmetphi[0]) + (l1+l2).Px(); 
#            py_rec = pfmet[0] * sin(pfmetphi[0]) + (l1+l2).Py(); 
#
#            recoil[0] = sqrt(px_rec*px_rec + py_rec*py_rec)
#
#            self.var1[0]  =  ppfmet[0]  
#            self.var2[0]  =  self.oldttree.pchmet
#            self.var3[0]  =  self.oldttree.nvtx
#            self.var4[0]  =  self.oldttree.ptll
#            self.var5[0]  =  jetpt1
#            self.var6[0]  =  self.oldttree.pfmetMEtSig #pfMetMEtSig
#            self.var7[0]  =  dphilljet1[0]  
#            self.var8[0]  =  dphillmet[0]
#            self.var9[0]  =  dphimetjet1[0]  
#            self.var10[0] =  recoil[0]  
#            self.var11[0] =  mth[0]
#
#
#            if self.oldttree.njet == 0:
#                  dymva1[0] = self.getDYMVAV1j0.EvaluateMVA("BDTG")
#            elif self.oldttree.njet == 1:
#                  dymva1[0] = self.getDYMVAV1j1.EvaluateMVA("BDTG")
#            else :
#                  dymva1[0] = -999
#


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
    
    parser.add_option('-i', '--inputTreeFileName' ,  dest='inputTreeFileName' ,  help='Name of the input  *.root file',)
    parser.add_option('-o', '--outputTreeFileName',  dest='outputTreeFileName',  help='Name of the output *.root file',)
    parser.add_option('-t', '--treeName'          ,  dest='treeDir'           ,  help='Name of the tree ',)
    parser.add_option('-k', '--kindMCDATA'        ,  dest='kindMCDATA'        ,  help='kind of sample: MC2011, MC2012, DATA2011, DATA2012',)


    

    (opt, args) = parser.parse_args()

    if opt.inputTreeFileName is None:
        parser.error('No input file defined')
    if opt.outputTreeFileName is None:
        parser.error('No output file defined')
    if opt.kindMCDATA is None:
        parser.error('No kindMCDATA branch name defined')

##     sys.argv.append('-b')
##     ROOT.gROOT.SetBatch()

    w = JetAndMetSmear()

    w.inputTreeFileName  = opt.inputTreeFileName
    w.outputTreeFileName = opt.outputTreeFileName
    w.treeDir            = opt.treeDir
    w.kindMCDATA         = opt.kindMCDATA




#    print s.systArgument

    w.createDYMVA()

    w.openOriginalTFile()
    w.openOutputTFile()
    w.cloneTree()
    
    w.applyMetSmearing()
    


    print 'Job finished...'


if __name__ == '__main__':
    main()



