import FWCore.ParameterSet.Config as cms
import MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi

### DR=0.1 cone
# Charged Hadron isolation
muonDirPFIsoChHad01 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonDirPFIsoChHad01.pfTypes = cms.untracked.vint32(1)
muonDirPFIsoChHad01.deltaR = 0.1
muonDirPFIsoChHad01.directional = True

# Neutral Hadron isolation
muonDirPFIsoNHad01 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonDirPFIsoNHad01.pfTypes = cms.untracked.vint32(5)
muonDirPFIsoNHad01.deltaR = 0.1
muonDirPFIsoNHad01.directional = True

# Photon isolation
muonDirPFIsoPhoton01 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonDirPFIsoPhoton01.pfTypes = cms.untracked.vint32(4)
muonDirPFIsoPhoton01.deltaR = 0.1
muonDirPFIsoPhoton01.directional = True


### DR=0.2 cone
# Charged Hadron isolation
muonDirPFIsoChHad02 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonDirPFIsoChHad02.pfTypes = cms.untracked.vint32(1)
muonDirPFIsoChHad02.deltaR = 0.2
muonDirPFIsoChHad02.directional = True

# Neutral Hadron isolation
muonDirPFIsoNHad02 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonDirPFIsoNHad02.pfTypes = cms.untracked.vint32(5)
muonDirPFIsoNHad02.deltaR = 0.2
muonDirPFIsoNHad02.directional = True

# Photon isolation
muonDirPFIsoPhoton02 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonDirPFIsoPhoton02.pfTypes = cms.untracked.vint32(4)
muonDirPFIsoPhoton02.deltaR = 0.2
muonDirPFIsoPhoton02.directional = True


### DR=0.3 cone
# Charged Hadron isolation
muonDirPFIsoChHad03 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonDirPFIsoChHad03.pfTypes = cms.untracked.vint32(1)
muonDirPFIsoChHad03.deltaR = 0.3
muonDirPFIsoChHad03.directional = True

# Neutral Hadron isolation
muonDirPFIsoNHad03 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonDirPFIsoNHad03.pfTypes = cms.untracked.vint32(5)
muonDirPFIsoNHad03.deltaR = 0.3
muonDirPFIsoNHad03.directional = True

# Photon isolation
muonDirPFIsoPhoton03 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonDirPFIsoPhoton03.pfTypes = cms.untracked.vint32(4)
muonDirPFIsoPhoton03.deltaR = 0.3
muonDirPFIsoPhoton03.directional = True


### DR=0.4 cone
# Charged Hadron isolation
muonDirPFIsoChHad04 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonDirPFIsoChHad04.pfTypes = cms.untracked.vint32(1)
muonDirPFIsoChHad04.deltaR = 0.4
muonDirPFIsoChHad04.directional = True

# Neutral Hadron isolation
muonDirPFIsoNHad04 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonDirPFIsoNHad04.pfTypes = cms.untracked.vint32(5)
muonDirPFIsoNHad04.deltaR = 0.4
muonDirPFIsoNHad04.directional = True

# Photon isolation
muonDirPFIsoPhoton04 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonDirPFIsoPhoton04.pfTypes = cms.untracked.vint32(4)
muonDirPFIsoPhoton04.deltaR = 0.4
muonDirPFIsoPhoton04.directional = True


### DR=0.5 cone
# Charged Hadron isolation
muonDirPFIsoChHad05 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonDirPFIsoChHad05.pfTypes = cms.untracked.vint32(1)
muonDirPFIsoChHad05.deltaR = 0.5
muonDirPFIsoChHad05.directional = True

# Neutral Hadron isolation
muonDirPFIsoNHad05 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonDirPFIsoNHad05.pfTypes = cms.untracked.vint32(5)
muonDirPFIsoNHad05.deltaR = 0.5
muonDirPFIsoNHad05.directional = True

# Photon isolation
muonDirPFIsoPhoton05 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonDirPFIsoPhoton05.pfTypes = cms.untracked.vint32(4)
muonDirPFIsoPhoton05.deltaR = 0.5
muonDirPFIsoPhoton05.directional = True


### DR=0.6 cone
# Charged Hadron isolation
muonDirPFIsoChHad06 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonDirPFIsoChHad06.pfTypes = cms.untracked.vint32(1)
muonDirPFIsoChHad06.deltaR = 0.6
muonDirPFIsoChHad06.directional = True

# Neutral Hadron isolation
muonDirPFIsoNHad06 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonDirPFIsoNHad06.pfTypes = cms.untracked.vint32(5)
muonDirPFIsoNHad06.deltaR = 0.6
muonDirPFIsoNHad06.directional = True

# Photon isolation
muonDirPFIsoPhoton06 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonDirPFIsoPhoton06.pfTypes = cms.untracked.vint32(4)
muonDirPFIsoPhoton06.deltaR = 0.6
muonDirPFIsoPhoton06.directional = True


### DR=0.7 cone
# Charged Hadron isolation
muonDirPFIsoChHad07 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonDirPFIsoChHad07.pfTypes = cms.untracked.vint32(1)
muonDirPFIsoChHad07.deltaR = 0.7
muonDirPFIsoChHad07.directional = True

# Neutral Hadron isolation
muonDirPFIsoNHad07 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonDirPFIsoNHad07.pfTypes = cms.untracked.vint32(5)
muonDirPFIsoNHad07.deltaR = 0.7
muonDirPFIsoNHad07.directional = True

# Photon isolation
muonDirPFIsoPhoton07 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonDirPFIsoPhoton07.pfTypes = cms.untracked.vint32(4)
muonDirPFIsoPhoton07.deltaR = 0.7
muonDirPFIsoPhoton07.directional = True



