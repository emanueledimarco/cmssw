# print sys.argv
# args = sys.argv[1:]
import getopt
import sys

from WWAnalysis.AnalysisStep.data.v03.nov4Samples_cff       import nov4Samples
from WWAnalysis.AnalysisStep.data.v03.sep172010BSamples_cff import sep172010BSamples
from WWAnalysis.AnalysisStep.data.v03.wJetsSamples_cff      import wJetsSamples
from WWAnalysis.AnalysisStep.data.v03.wwSamples_cff         import wwSamples
from WWAnalysis.AnalysisStep.data.v03.higgsSamples_cff      import higgsSamples
from WWAnalysis.AnalysisStep.data.v03.failureSamples_cff    import failureSamples


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:],
                                    "an:vswjgtc:f:o:h",
                                   ["all",
                                    "files-per-job=",
                                    "nov4Samples",
                                    "sep172010BSamples",
                                    "wwSamples",
                                    "wJetsSamples",
                                    "higgsSamples",
                                    "failureSamples",
                                    "cut-sel=",
                                    "filterMass",
                                    "output-dir=",
                                    "help"])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage(sys.argv[0])
        sys.exit(1)

    
    doAll = False
    doNov4 = False
    doSep17 = False
    doWW = True
    doWJet = False
    doHiggs = False
    doFailure = False
    filterMass = ""
    filesPer = 1
    outputDir = "data/"
    cutSel="ww"

    for o, a in opts:
        if o in ('-a','--all'):
            doAll = True
        elif o in ('-n', '--files-per-job'):
            filesPer = int(a)
        elif o in ('-v', '--nov4Samples'):
            doNov4 = True;
            doFailure = doHiggs = doSep17 = doWW = doWJet = False;
        elif o in ('-s', '--sep172010BSamples'):
            doSep17 = True;
            doFailure = doHiggs = doNov4 = doWW = doWJet = False;
        elif o in ('-w', '--wwSamples'):
            doWW = True;
            doFailure = doHiggs = doNov4 = doSep17 = doWJet = False;
        elif o in ('-j', '--wJetsSamples'):
            doWJet = True;
            doFailure = doHiggs = doNov4 = doSep17 = doWW = False;
        elif o in ('-g', '--higgsSamples'):
            doHiggs = True;
            doFailure = doWJet = doNov4 = doSep17 = doWW = False;
        elif o in ('-t', '--failureSamples'):
            doFailure = True;
            doHiggs = doWJet = doNov4 = doSep17 = doWW = False;
        elif o in ('-c', '--cut-sel'):
            cutSel = a
        elif o in ('-f', '--filterMass'):
            filterMass = a
        elif o in ('-o', '--output-dir'):
            outputDir = a
        elif o in ('-h', '--help'):
            usage(sys.argv[0])
            sys.exit(0)  
        else:
            assert False, "unhandled option"

    if doAll:
        doHiggs = doNov4 = doSep17 = doWW = doWJet = True;

    
    if doNov4:
        makeFiles(nov4Samples,filesPer,cutSel,filterMass,outputDir) 
    if doSep17:
        makeFiles(sep172010BSamples,filesPer,cutSel,filterMass,outputDir) 
    if doWJet:
        makeFiles(wJetsSamples,filesPer,cutSel,filterMass,outputDir) 
    if doWW:
        makeFiles(wwSamples,filesPer,cutSel,filterMass,outputDir) 
    if doHiggs:
        makeFiles(higgsSamples,filesPer,cutSel,filterMass,outputDir) 
    if doFailure:
        makeFiles(failureSamples,filesPer,cutSel,filterMass,outputDir) 




def makeFiles(myDicts,filesPer,cutSel,filterMass,outputDir):
    for key in myDicts:
        print key,cutSel,filterMass
#         if not filterMass or key[5:] in cutSel:
        if  filterMass == "" or filterMass in key:
            for i in range(0,len(myDicts[key][3]),filesPer):
    
                f = open( "job."+key+"."+myDicts[key][0]+"."+str(i)+".py" , 'w')
                f.write(prefix.format(cutSel))
                if( myDicts[key][2] ) : f.write(dataLumiString)
    
                files = myDicts[key][3][i:min(i+filesPer,len(myDicts[key][3]))]
                fileString = ""
                for file in files:
                    fileString += "'" + str(file) + "'" + ", "
                thisString =  'process.yieldParams.inputParams = cms.PSet(\n'
                thisString += '    %s = cms.PSet(\n' % myDicts[key][0]
                thisString += '        files = cms.vstring(%s),\n'  % fileString
                thisString += '        scale = cms.double(%e),\n' % myDicts[key][1]
                thisString += '        printEvents = cms.bool(%s),\n' % str(myDicts[key][2])
                thisString += '        isMC = cms.bool(%s),\n' % str(not myDicts[key][2])
                thisString += '        isSignal = cms.bool(%s),\n' % str(myDicts[key][5])
                thisString += '    )\n'
                thisString += ')\n'
    
                f.write(thisString)
                f.write("process.yieldParams.fileOutName = cms.string('%s/%s')\n" % ((outputDir),("job."+key+"."+myDicts[key][0]+"."+str(i)+".root")))
                f.close()


prefix = """
import FWCore.ParameterSet.Config as cms

from FWCore.PythonUtilities.LumiList import LumiList
process = cms.Process("HistoProducer")
process.load("WWAnalysis.AnalysisStep.yieldProducer_cfi")
process.yieldParams = process.FWLiteParams.clone()

from WWAnalysis.AnalysisStep.cutPSets_cfi import *

swapOutPSetValues(process.yieldParams.selectionParams.wwelel,{0})
process.yieldParams.selectionParams.wwmumu = process.yieldParams.selectionParams.wwelel.clone()
process.yieldParams.selectionParams.wwelmu = process.yieldParams.selectionParams.wwelel.clone()
swapOutPSetValues(process.yieldParams.selectionParams.wwelmu,oppositeFlavor)

"""

dataLumiString = """
lumis = LumiList(filename = '/nfs/bluearc/group/jsons/Skim_Mu_And_El_byGiovanni/Final.json')
process.yieldParams.lumisToProcess = cms.untracked.VLuminosityBlockRange()
process.yieldParams.lumisToProcess = lumis.getCMSSWString().split(',')
process.yieldParams.resFile = 'WWAnalysis/AnalysisStep/data/START38_V13_AK5PF_L2L3Residual.txt'

"""





def usage(a):
    print """
Usage:

%s [OPTIONS] 

Options:

 -h, --help              This message.
 -a, --all               Use all available files (trumps all) (default False)
 -n, --files-per-job=    Number of files per job  (default 1)
 -v, --nov4Samples       Use Nov 4th (default False)
 -s, --sep172010BSamples Use Sep 17th (default False)
 -w, --wwSamples         Use WW and Backgrounds (default True)
 -g, --higgsSamples      Use Higgs Samples (default False)
 -t, --failureSamples    Use Failure modes samples (default False)
 -j, --wJetsSamples      Use WJets Samples (default False)
 -c, --cut-sel=          which sel of cuts (ww by default, could be h160 e.g.)
 -f, --filterMass        only make the samples if they mactch the mass of the cut sel
 """ % a


if __name__ == "__main__":
    main()

