import FWCore.ParameterSet.Config as cms
import MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi

### DR=0.1 cone
# Charged Hadron isolation
electronPFIsoChHad01 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoChHad01.pfTypes = cms.untracked.vint32(1)
electronPFIsoChHad01.deltaR = 0.1

# Neutral Hadron isolation
electronPFIsoNHad01 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoNHad01.pfTypes = cms.untracked.vint32(5)
electronPFIsoNHad01.deltaR = 0.1

# Photon isolation
electronPFIsoPhoton01 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoPhoton01.pfTypes = cms.untracked.vint32(4)
electronPFIsoPhoton01.deltaR = 0.1


### DR=0.2 cone
# Charged Hadron isolation
electronPFIsoChHad02 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoChHad02.pfTypes = cms.untracked.vint32(1)
electronPFIsoChHad02.deltaR = 0.2

# Neutral Hadron isolation
electronPFIsoNHad02 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoNHad02.pfTypes = cms.untracked.vint32(5)
electronPFIsoNHad02.deltaR = 0.2

# Photon isolation
electronPFIsoPhoton02 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoPhoton02.pfTypes = cms.untracked.vint32(4)
electronPFIsoPhoton02.deltaR = 0.2


### DR=0.3 cone
# Charged Hadron isolation
electronPFIsoChHad03 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoChHad03.pfTypes = cms.untracked.vint32(1)
electronPFIsoChHad03.deltaR = 0.3

# Neutral Hadron isolation
electronPFIsoNHad03 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoNHad03.pfTypes = cms.untracked.vint32(5)
electronPFIsoNHad03.deltaR = 0.3

# Photon isolation
electronPFIsoPhoton03 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoPhoton03.pfTypes = cms.untracked.vint32(4)
electronPFIsoPhoton03.deltaR = 0.3


### DR=0.4 cone
# Charged Hadron isolation
electronPFIsoChHad04 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoChHad04.pfTypes = cms.untracked.vint32(1)
electronPFIsoChHad04.deltaR = 0.4

# Neutral Hadron isolation
electronPFIsoNHad04 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoNHad04.pfTypes = cms.untracked.vint32(5)
electronPFIsoNHad04.deltaR = 0.4

# Photon isolation
electronPFIsoPhoton04 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoPhoton04.pfTypes = cms.untracked.vint32(4)
electronPFIsoPhoton04.deltaR = 0.4


### DR=0.5 cone
# Charged Hadron isolation
electronPFIsoChHad05 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoChHad05.pfTypes = cms.untracked.vint32(1)
electronPFIsoChHad05.deltaR = 0.5

# Neutral Hadron isolation
electronPFIsoNHad05 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoNHad05.pfTypes = cms.untracked.vint32(5)
electronPFIsoNHad05.deltaR = 0.5

# Photon isolation
electronPFIsoPhoton05 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoPhoton05.pfTypes = cms.untracked.vint32(4)
electronPFIsoPhoton05.deltaR = 0.5


### DR=0.6 cone
# Charged Hadron isolation
electronPFIsoChHad06 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoChHad06.pfTypes = cms.untracked.vint32(1)
electronPFIsoChHad06.deltaR = 0.6

# Neutral Hadron isolation
electronPFIsoNHad06 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoNHad06.pfTypes = cms.untracked.vint32(5)
electronPFIsoNHad06.deltaR = 0.6

# Photon isolation
electronPFIsoPhoton06 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoPhoton06.pfTypes = cms.untracked.vint32(4)
electronPFIsoPhoton06.deltaR = 0.6


### DR=0.7 cone
# Charged Hadron isolation
electronPFIsoChHad07 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoChHad07.pfTypes = cms.untracked.vint32(1)
electronPFIsoChHad07.deltaR = 0.7

# Neutral Hadron isolation
electronPFIsoNHad07 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoNHad07.pfTypes = cms.untracked.vint32(5)
electronPFIsoNHad07.deltaR = 0.7

# Photon isolation
electronPFIsoPhoton07 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoPhoton07.pfTypes = cms.untracked.vint32(4)
electronPFIsoPhoton07.deltaR = 0.7
