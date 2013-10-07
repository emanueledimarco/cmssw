import FWCore.ParameterSet.Config as cms
import MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi

### DR=0.1 cone
# Charged Hadron isolation
muonRadPFIsoChHad01 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonRadPFIsoChHad01.pfTypes = cms.untracked.vint32(1)
muonRadPFIsoChHad01.deltaR = 0.1
muonRadPFIsoChHad01.radial = True

# Neutral Hadron isolation
muonRadPFIsoNHad01 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonRadPFIsoNHad01.pfTypes = cms.untracked.vint32(5)
muonRadPFIsoNHad01.deltaR = 0.1
muonRadPFIsoNHad01.radial = True

# Photon isolation
muonRadPFIsoPhoton01 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonRadPFIsoPhoton01.pfTypes = cms.untracked.vint32(4)
muonRadPFIsoPhoton01.deltaR = 0.1
muonRadPFIsoPhoton01.radial = True


### DR=0.2 cone
# Charged Hadron isolation
muonRadPFIsoChHad02 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonRadPFIsoChHad02.pfTypes = cms.untracked.vint32(1)
muonRadPFIsoChHad02.deltaR = 0.2
muonRadPFIsoChHad02.radial = True

# Neutral Hadron isolation
muonRadPFIsoNHad02 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonRadPFIsoNHad02.pfTypes = cms.untracked.vint32(5)
muonRadPFIsoNHad02.deltaR = 0.2
muonRadPFIsoNHad02.radial = True

# Photon isolation
muonRadPFIsoPhoton02 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonRadPFIsoPhoton02.pfTypes = cms.untracked.vint32(4)
muonRadPFIsoPhoton02.deltaR = 0.2
muonRadPFIsoPhoton02.radial = True


### DR=0.3 cone
# Charged Hadron isolation
muonRadPFIsoChHad03 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonRadPFIsoChHad03.pfTypes = cms.untracked.vint32(1)
muonRadPFIsoChHad03.deltaR = 0.3
muonRadPFIsoChHad03.radial = True

# Neutral Hadron isolation
muonRadPFIsoNHad03 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonRadPFIsoNHad03.pfTypes = cms.untracked.vint32(5)
muonRadPFIsoNHad03.deltaR = 0.3
muonRadPFIsoNHad03.radial = True

# Photon isolation
muonRadPFIsoPhoton03 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonRadPFIsoPhoton03.pfTypes = cms.untracked.vint32(4)
muonRadPFIsoPhoton03.deltaR = 0.3
muonRadPFIsoPhoton03.radial = True


### DR=0.4 cone
# Charged Hadron isolation
muonRadPFIsoChHad04 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonRadPFIsoChHad04.pfTypes = cms.untracked.vint32(1)
muonRadPFIsoChHad04.deltaR = 0.4
muonRadPFIsoChHad04.radial = True

# Neutral Hadron isolation
muonRadPFIsoNHad04 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonRadPFIsoNHad04.pfTypes = cms.untracked.vint32(5)
muonRadPFIsoNHad04.deltaR = 0.4
muonRadPFIsoNHad04.radial = True

# Photon isolation
muonRadPFIsoPhoton04 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonRadPFIsoPhoton04.pfTypes = cms.untracked.vint32(4)
muonRadPFIsoPhoton04.deltaR = 0.4
muonRadPFIsoPhoton04.radial = True


### DR=0.5 cone
# Charged Hadron isolation
muonRadPFIsoChHad05 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonRadPFIsoChHad05.pfTypes = cms.untracked.vint32(1)
muonRadPFIsoChHad05.deltaR = 0.5
muonRadPFIsoChHad05.radial = True

# Neutral Hadron isolation
muonRadPFIsoNHad05 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonRadPFIsoNHad05.pfTypes = cms.untracked.vint32(5)
muonRadPFIsoNHad05.deltaR = 0.5
muonRadPFIsoNHad05.radial = True

# Photon isolation
muonRadPFIsoPhoton05 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonRadPFIsoPhoton05.pfTypes = cms.untracked.vint32(4)
muonRadPFIsoPhoton05.deltaR = 0.5
muonRadPFIsoPhoton05.radial = True


### DR=0.6 cone
# Charged Hadron isolation
muonRadPFIsoChHad06 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonRadPFIsoChHad06.pfTypes = cms.untracked.vint32(1)
muonRadPFIsoChHad06.deltaR = 0.6
muonRadPFIsoChHad06.radial = True

# Neutral Hadron isolation
muonRadPFIsoNHad06 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonRadPFIsoNHad06.pfTypes = cms.untracked.vint32(5)
muonRadPFIsoNHad06.deltaR = 0.6
muonRadPFIsoNHad06.radial = True

# Photon isolation
muonRadPFIsoPhoton06 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonRadPFIsoPhoton06.pfTypes = cms.untracked.vint32(4)
muonRadPFIsoPhoton06.deltaR = 0.6
muonRadPFIsoPhoton06.radial = True


### DR=0.7 cone
# Charged Hadron isolation
muonRadPFIsoChHad07 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonRadPFIsoChHad07.pfTypes = cms.untracked.vint32(1)
muonRadPFIsoChHad07.deltaR = 0.7
muonRadPFIsoChHad07.radial = True

# Neutral Hadron isolation
muonRadPFIsoNHad07 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonRadPFIsoNHad07.pfTypes = cms.untracked.vint32(5)
muonRadPFIsoNHad07.deltaR = 0.7
muonRadPFIsoNHad07.radial = True

# Photon isolation
muonRadPFIsoPhoton07 = MyAnalysis.IsolationTools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonRadPFIsoPhoton07.pfTypes = cms.untracked.vint32(4)
muonRadPFIsoPhoton07.deltaR = 0.7
muonRadPFIsoPhoton07.radial = True



