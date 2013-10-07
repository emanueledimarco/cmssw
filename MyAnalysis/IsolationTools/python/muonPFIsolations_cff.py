import FWCore.ParameterSet.Config as cms
import MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi

### DR=0.1
# Charged Hadron isolation
muonPFIsoChHad01 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoChHad01.pfTypes = cms.untracked.vint32(1)
muonPFIsoChHad01.deltaR = 0.1

# Neutral Hadron isolation
muonPFIsoNHad01 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoNHad01.pfTypes = cms.untracked.vint32(5)
muonPFIsoNHad01.deltaR = 0.1

# Photon isolation
muonPFIsoPhoton01 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoPhoton01.pfTypes = cms.untracked.vint32(4)
muonPFIsoPhoton01.deltaR = 0.1


### DR=0.2
# Charged Hadron isolation
muonPFIsoChHad02 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoChHad02.pfTypes = cms.untracked.vint32(1)
muonPFIsoChHad02.deltaR = 0.2

# Neutral Hadron isolation
muonPFIsoNHad02 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoNHad02.pfTypes = cms.untracked.vint32(5)
muonPFIsoNHad02.deltaR = 0.2

# Photon isolation
muonPFIsoPhoton02 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoPhoton02.pfTypes = cms.untracked.vint32(4)
muonPFIsoPhoton02.deltaR = 0.2


### DR=0.3
# Charged Hadron isolation
muonPFIsoChHad03 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoChHad03.pfTypes = cms.untracked.vint32(1)
muonPFIsoChHad03.deltaR = 0.3

# Neutral Hadron isolation
muonPFIsoNHad03 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoNHad03.pfTypes = cms.untracked.vint32(5)
muonPFIsoNHad03.deltaR = 0.3

# Photon isolation
muonPFIsoPhoton03 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoPhoton03.pfTypes = cms.untracked.vint32(4)
muonPFIsoPhoton03.deltaR = 0.3


### DR=0.4
# Charged Hadron isolation
muonPFIsoChHad04 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoChHad04.pfTypes = cms.untracked.vint32(1)
muonPFIsoChHad04.deltaR = 0.4

# Neutral Hadron isolation
muonPFIsoNHad04 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoNHad04.pfTypes = cms.untracked.vint32(5)
muonPFIsoNHad04.deltaR = 0.4

# Photon isolation
muonPFIsoPhoton04 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoPhoton04.pfTypes = cms.untracked.vint32(4)
muonPFIsoPhoton04.deltaR = 0.4


### DR=0.5
# Charged Hadron isolation
muonPFIsoChHad05 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoChHad05.pfTypes = cms.untracked.vint32(1)
muonPFIsoChHad05.deltaR = 0.5

# Neutral Hadron isolation
muonPFIsoNHad05 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoNHad05.pfTypes = cms.untracked.vint32(5)
muonPFIsoNHad05.deltaR = 0.5

# Photon isolation
muonPFIsoPhoton05 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoPhoton05.pfTypes = cms.untracked.vint32(4)
muonPFIsoPhoton05.deltaR = 0.5


### DR=0.6
# Charged Hadron isolation
muonPFIsoChHad06 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoChHad06.pfTypes = cms.untracked.vint32(1)
muonPFIsoChHad06.deltaR = 0.6

# Neutral Hadron isolation
muonPFIsoNHad06 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoNHad06.pfTypes = cms.untracked.vint32(5)
muonPFIsoNHad06.deltaR = 0.6

# Photon isolation
muonPFIsoPhoton06 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoPhoton06.pfTypes = cms.untracked.vint32(4)
muonPFIsoPhoton06.deltaR = 0.6


### DR=0.7
# Charged Hadron isolation
muonPFIsoChHad07 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoChHad07.pfTypes = cms.untracked.vint32(1)
muonPFIsoChHad07.deltaR = 0.7

# Neutral Hadron isolation
muonPFIsoNHad07 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoNHad07.pfTypes = cms.untracked.vint32(5)
muonPFIsoNHad07.deltaR = 0.7

# Photon isolation
muonPFIsoPhoton07 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoPhoton07.pfTypes = cms.untracked.vint32(4)
muonPFIsoPhoton07.deltaR = 0.7
