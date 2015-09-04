import FWCore.ParameterSet.Config as cms
import re
process = cms.Process("PEDS")

######## configure here #######
calcPedestalMean = False
run = 211831
barrel = True
FEDused = "EB+13"
Gain = 12
###############################

templateOutputFile = str("templates_run"+str(run)+"_"+FEDused+"_Gain"+str(Gain)+".root")
templateOutputFile = re.sub("\+","p", templateOutputFile)
templateOutputFile = re.sub("\-","m", templateOutputFile)

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.Geometry.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.XMLFromDBSource.label = cms.string("Extended2015ZeroMaterial")
process.GlobalTag.globaltag = 'POSTLS172_V7::All'

process.load('DPGAnalysis.EcalTools.pulseDump_cfi')

# Run2013A 
# sitestring = 'root://cms-xrd-global.cern.ch/'
# process.source = cms.Source("PoolSource",
#                             duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
#                             # run 211831, low PU run at 2.76 TeV in 2013 
#                             fileNames = cms.untracked.vstring( sitestring+'/store/data/Run2013A/PPPhoton/RAW/v1/000/211/831/00000/0C09133A-A276-E211-ACA3-E0CB4E55365C.root',
#                                                                sitestring+'/store/data/Run2013A/PPPhoton/RAW/v1/000/211/831/00000/1423C319-AA76-E211-B8B4-BCAEC5329716.root',
#                                                                sitestring+'/store/data/Run2013A/PPPhoton/RAW/v1/000/211/831/00000/1A9DBCEB-C276-E211-AB2E-E0CB4E4408E3.root',
#                                                                sitestring+'/store/data/Run2013A/PPPhoton/RAW/v1/000/211/831/00000/26AB1CFE-C276-E211-84CC-003048D2BEA8.root',
#                                                                sitestring+'/store/data/Run2013A/PPPhoton/RAW/v1/000/211/831/00000/289C17F2-C276-E211-A6CF-0025901D624A.root',
#                                                                sitestring+'/store/data/Run2013A/PPPhoton/RAW/v1/000/211/831/00000/308B3A4B-AF76-E211-8D76-001D09F24EE3.root',
#                                                                sitestring+'/store/data/Run2013A/PPPhoton/RAW/v1/000/211/831/00000/3486B9F0-C276-E211-B27D-0025901D6288.root',
#                                                                sitestring+'/store/data/Run2013A/PPPhoton/RAW/v1/000/211/831/00000/489747ED-C276-E211-8890-0025901D5C88.root',
#                                                                sitestring+'/store/data/Run2013A/PPPhoton/RAW/v1/000/211/831/00000/4C1661F1-C276-E211-ADD5-5404A63886B4.root',
#                                                                sitestring+'/store/data/Run2013A/PPPhoton/RAW/v1/000/211/831/00000/4EAC03EE-C276-E211-A527-0025901D5D78.root',
#                                                                sitestring+'/store/data/Run2013A/PPPhoton/RAW/v1/000/211/831/00000/546625F2-A576-E211-AD8B-003048F024E0.root',
#                                                                sitestring+'/store/data/Run2013A/PPPhoton/RAW/v1/000/211/831/00000/5C776549-AF76-E211-812E-0019B9F72CE5.root',
#                                                                sitestring+'/store/data/Run2013A/PPPhoton/RAW/v1/000/211/831/00000/7A0D87F5-A576-E211-922B-00237DDBE49C.root',
#                                                                sitestring+'/store/data/Run2013A/PPPhoton/RAW/v1/000/211/831/00000/7CD30F64-C376-E211-BC54-0025901D5DF4.root',
#                                                                sitestring+'/store/data/Run2013A/PPPhoton/RAW/v1/000/211/831/00000/7EAA71ED-C276-E211-8E79-5404A63886CE.root',
#                                                                sitestring+'/store/data/Run2013A/PPPhoton/RAW/v1/000/211/831/00000/82239BEC-C276-E211-AB79-0025901D62A6.root',
#                                                                sitestring+'/store/data/Run2013A/PPPhoton/RAW/v1/000/211/831/00000/9E08094B-AF76-E211-A45A-001D09F28EA3.root',
#                                                                sitestring+'/store/data/Run2013A/PPPhoton/RAW/v1/000/211/831/00000/B4E67CED-C276-E211-9BF5-003048F118C2.root',
#                                                                sitestring+'/store/data/Run2013A/PPPhoton/RAW/v1/000/211/831/00000/D0CC5EEB-C276-E211-9EDB-5404A63886EB.root',
#                                                                sitestring+'/store/data/Run2013A/PPPhoton/RAW/v1/000/211/831/00000/DE3D95F4-C276-E211-8A52-003048D2BC5C.root',
#                                                                sitestring+'/store/data/Run2013A/PPPhoton/RAW/v1/000/211/831/00000/E6E326F1-C276-E211-8D57-0025901D6268.root',
#                                                                sitestring+'/store/data/Run2013A/PPPhoton/RAW/v1/000/211/831/00000/E865D7EE-C276-E211-9BEE-0025901D623C.root',
#                                                                sitestring+'/store/data/Run2013A/PPPhoton/RAW/v1/000/211/831/00000/F09165EB-C276-E211-B971-BCAEC53296F8.root',
#                                                                sitestring+'/store/data/Run2013A/PPPhoton/RAW/v1/000/211/831/00000/F823A949-AF76-E211-886C-001D09F24D8A.root',
#                                                                sitestring+'/store/data/Run2013A/PPPhoton/RAW/v1/000/211/831/00000/F8BB23F1-C276-E211-A244-0025901D5DF4.root')
#                             )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('/store/data/Run2012A/DoubleElectron/RAW-RECO/ZElectron-22Jan2013-v1/30000/B06652D0-7170-E211-9BD0-00304867BFBC.root')
                            )

process.TFileService = cms.Service("TFileService", fileName = cms.string("templates.root") )

# start from RAW format for more flexibility
process.raw2digi_step = cms.Sequence(process.RawToDigi)

process.p = cms.Path(process.ecalDigis * 
                     process.pulseDump )
