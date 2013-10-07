import FWCore.ParameterSet.Config as cms
import MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi

### DR=0.1 cone
# Charged Hadron isolation
electronDirPFIsoChHad01 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronDirPFIsoChHad01.pfTypes = cms.untracked.vint32(1)
electronDirPFIsoChHad01.deltaR = 0.1
electronDirPFIsoChHad01.directional = True

# Neutral Hadron isolation
electronDirPFIsoNHad01 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronDirPFIsoNHad01.pfTypes = cms.untracked.vint32(5)
electronDirPFIsoNHad01.deltaR = 0.1
electronDirPFIsoNHad01.directional = True

# Photon isolation
electronDirPFIsoPhoton01 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronDirPFIsoPhoton01.pfTypes = cms.untracked.vint32(4)
electronDirPFIsoPhoton01.deltaR = 0.1
electronDirPFIsoPhoton01.directional = True


### DR=0.2 cone
# Charged Hadron isolation
electronDirPFIsoChHad02 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronDirPFIsoChHad02.pfTypes = cms.untracked.vint32(1)
electronDirPFIsoChHad02.deltaR = 0.2
electronDirPFIsoChHad02.directional = True

# Neutral Hadron isolation
electronDirPFIsoNHad02 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronDirPFIsoNHad02.pfTypes = cms.untracked.vint32(5)
electronDirPFIsoNHad02.deltaR = 0.2
electronDirPFIsoNHad02.directional = True

# Photon isolation
electronDirPFIsoPhoton02 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronDirPFIsoPhoton02.pfTypes = cms.untracked.vint32(4)
electronDirPFIsoPhoton02.deltaR = 0.2
electronDirPFIsoPhoton02.directional = True


### DR=0.3 cone
# Charged Hadron isolation
electronDirPFIsoChHad03 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronDirPFIsoChHad03.pfTypes = cms.untracked.vint32(1)
electronDirPFIsoChHad03.deltaR = 0.3
electronDirPFIsoChHad03.directional = True

# Neutral Hadron isolation
electronDirPFIsoNHad03 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronDirPFIsoNHad03.pfTypes = cms.untracked.vint32(5)
electronDirPFIsoNHad03.deltaR = 0.3
electronDirPFIsoNHad03.directional = True

# Photon isolation
electronDirPFIsoPhoton03 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronDirPFIsoPhoton03.pfTypes = cms.untracked.vint32(4)
electronDirPFIsoPhoton03.deltaR = 0.3
electronDirPFIsoPhoton03.directional = True


### DR=0.4 cone
# Charged Hadron isolation
electronDirPFIsoChHad04 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronDirPFIsoChHad04.pfTypes = cms.untracked.vint32(1)
electronDirPFIsoChHad04.deltaR = 0.4
electronDirPFIsoChHad04.directional = True

# Neutral Hadron isolation
electronDirPFIsoNHad04 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronDirPFIsoNHad04.pfTypes = cms.untracked.vint32(5)
electronDirPFIsoNHad04.deltaR = 0.4
electronDirPFIsoNHad04.directional = True

# Photon isolation
electronDirPFIsoPhoton04 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronDirPFIsoPhoton04.pfTypes = cms.untracked.vint32(4)
electronDirPFIsoPhoton04.deltaR = 0.4
electronDirPFIsoPhoton04.directional = True


### DR=0.5 cone
# Charged Hadron isolation
electronDirPFIsoChHad05 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronDirPFIsoChHad05.pfTypes = cms.untracked.vint32(1)
electronDirPFIsoChHad05.deltaR = 0.5
electronDirPFIsoChHad05.directional = True

# Neutral Hadron isolation
electronDirPFIsoNHad05 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronDirPFIsoNHad05.pfTypes = cms.untracked.vint32(5)
electronDirPFIsoNHad05.deltaR = 0.5
electronDirPFIsoNHad05.directional = True

# Photon isolation
electronDirPFIsoPhoton05 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronDirPFIsoPhoton05.pfTypes = cms.untracked.vint32(4)
electronDirPFIsoPhoton05.deltaR = 0.5
electronDirPFIsoPhoton05.directional = True


### DR=0.6 cone
# Charged Hadron isolation
electronDirPFIsoChHad06 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronDirPFIsoChHad06.pfTypes = cms.untracked.vint32(1)
electronDirPFIsoChHad06.deltaR = 0.6
electronDirPFIsoChHad06.directional = True

# Neutral Hadron isolation
electronDirPFIsoNHad06 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronDirPFIsoNHad06.pfTypes = cms.untracked.vint32(5)
electronDirPFIsoNHad06.deltaR = 0.6
electronDirPFIsoNHad06.directional = True

# Photon isolation
electronDirPFIsoPhoton06 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronDirPFIsoPhoton06.pfTypes = cms.untracked.vint32(4)
electronDirPFIsoPhoton06.deltaR = 0.6
electronDirPFIsoPhoton06.directional = True


### DR=0.7 cone
# Charged Hadron isolation
electronDirPFIsoChHad07 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronDirPFIsoChHad07.pfTypes = cms.untracked.vint32(1)
electronDirPFIsoChHad07.deltaR = 0.7
electronDirPFIsoChHad07.directional = True

# Neutral Hadron isolation
electronDirPFIsoNHad07 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronDirPFIsoNHad07.pfTypes = cms.untracked.vint32(5)
electronDirPFIsoNHad07.deltaR = 0.7
electronDirPFIsoNHad07.directional = True

# Photon isolation
electronDirPFIsoPhoton07 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronDirPFIsoPhoton07.pfTypes = cms.untracked.vint32(4)
electronDirPFIsoPhoton07.deltaR = 0.7
electronDirPFIsoPhoton07.directional = True



