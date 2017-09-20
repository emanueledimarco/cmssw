import FWCore.ParameterSet.Config as cms

from PhysicsTools.NanoAOD.nano_cff import *

nanoAOD = cms.Sequence(nanoSequence)
nanoAOD_MC = cms.Sequence(nanoSequenceMC)
