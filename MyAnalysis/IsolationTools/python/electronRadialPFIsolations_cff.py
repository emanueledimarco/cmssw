import FWCore.ParameterSet.Config as cms
import MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi

### DR=0.1 cone
# Charged Hadron isolation
electronRadPFIsoChHad01 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronRadPFIsoChHad01.pfTypes = cms.untracked.vint32(1)
electronRadPFIsoChHad01.deltaR = 0.1
electronRadPFIsoChHad01.radial = True

# Neutral Hadron isolation
electronRadPFIsoNHad01 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronRadPFIsoNHad01.pfTypes = cms.untracked.vint32(5)
electronRadPFIsoNHad01.deltaR = 0.1
electronRadPFIsoNHad01.radial = True

# Photon isolation
electronRadPFIsoPhoton01 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronRadPFIsoPhoton01.pfTypes = cms.untracked.vint32(4)
electronRadPFIsoPhoton01.deltaR = 0.1
electronRadPFIsoPhoton01.radial = True


### DR=0.2 cone
# Charged Hadron isolation
electronRadPFIsoChHad02 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronRadPFIsoChHad02.pfTypes = cms.untracked.vint32(1)
electronRadPFIsoChHad02.deltaR = 0.2
electronRadPFIsoChHad02.radial = True

# Neutral Hadron isolation
electronRadPFIsoNHad02 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronRadPFIsoNHad02.pfTypes = cms.untracked.vint32(5)
electronRadPFIsoNHad02.deltaR = 0.2
electronRadPFIsoNHad02.radial = True

# Photon isolation
electronRadPFIsoPhoton02 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronRadPFIsoPhoton02.pfTypes = cms.untracked.vint32(4)
electronRadPFIsoPhoton02.deltaR = 0.2
electronRadPFIsoPhoton02.radial = True


### DR=0.3 cone
# Charged Hadron isolation
electronRadPFIsoChHad03 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronRadPFIsoChHad03.pfTypes = cms.untracked.vint32(1)
electronRadPFIsoChHad03.deltaR = 0.3
electronRadPFIsoChHad03.radial = True

# Neutral Hadron isolation
electronRadPFIsoNHad03 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronRadPFIsoNHad03.pfTypes = cms.untracked.vint32(5)
electronRadPFIsoNHad03.deltaR = 0.3
electronRadPFIsoNHad03.radial = True

# Photon isolation
electronRadPFIsoPhoton03 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronRadPFIsoPhoton03.pfTypes = cms.untracked.vint32(4)
electronRadPFIsoPhoton03.deltaR = 0.3
electronRadPFIsoPhoton03.radial = True


### DR=0.4 cone
# Charged Hadron isolation
electronRadPFIsoChHad04 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronRadPFIsoChHad04.pfTypes = cms.untracked.vint32(1)
electronRadPFIsoChHad04.deltaR = 0.4
electronRadPFIsoChHad04.radial = True

# Neutral Hadron isolation
electronRadPFIsoNHad04 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronRadPFIsoNHad04.pfTypes = cms.untracked.vint32(5)
electronRadPFIsoNHad04.deltaR = 0.4
electronRadPFIsoNHad04.radial = True

# Photon isolation
electronRadPFIsoPhoton04 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronRadPFIsoPhoton04.pfTypes = cms.untracked.vint32(4)
electronRadPFIsoPhoton04.deltaR = 0.4
electronRadPFIsoPhoton04.radial = True


### DR=0.5 cone
# Charged Hadron isolation
electronRadPFIsoChHad05 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronRadPFIsoChHad05.pfTypes = cms.untracked.vint32(1)
electronRadPFIsoChHad05.deltaR = 0.5
electronRadPFIsoChHad05.radial = True

# Neutral Hadron isolation
electronRadPFIsoNHad05 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronRadPFIsoNHad05.pfTypes = cms.untracked.vint32(5)
electronRadPFIsoNHad05.deltaR = 0.5
electronRadPFIsoNHad05.radial = True

# Photon isolation
electronRadPFIsoPhoton05 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronRadPFIsoPhoton05.pfTypes = cms.untracked.vint32(4)
electronRadPFIsoPhoton05.deltaR = 0.5
electronRadPFIsoPhoton05.radial = True


### DR=0.6 cone
# Charged Hadron isolation
electronRadPFIsoChHad06 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronRadPFIsoChHad06.pfTypes = cms.untracked.vint32(1)
electronRadPFIsoChHad06.deltaR = 0.6
electronRadPFIsoChHad06.radial = True

# Neutral Hadron isolation
electronRadPFIsoNHad06 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronRadPFIsoNHad06.pfTypes = cms.untracked.vint32(5)
electronRadPFIsoNHad06.deltaR = 0.6
electronRadPFIsoNHad06.radial = True

# Photon isolation
electronRadPFIsoPhoton06 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronRadPFIsoPhoton06.pfTypes = cms.untracked.vint32(4)
electronRadPFIsoPhoton06.deltaR = 0.6
electronRadPFIsoPhoton06.radial = True


### DR=0.7 cone
# Charged Hadron isolation
electronRadPFIsoChHad07 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronRadPFIsoChHad07.pfTypes = cms.untracked.vint32(1)
electronRadPFIsoChHad07.deltaR = 0.7
electronRadPFIsoChHad07.radial = True

# Neutral Hadron isolation
electronRadPFIsoNHad07 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronRadPFIsoNHad07.pfTypes = cms.untracked.vint32(5)
electronRadPFIsoNHad07.deltaR = 0.7
electronRadPFIsoNHad07.radial = True

# Photon isolation
electronRadPFIsoPhoton07 = MyAnalysis.IsolationTools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronRadPFIsoPhoton07.pfTypes = cms.untracked.vint32(4)
electronRadPFIsoPhoton07.deltaR = 0.7
electronRadPFIsoPhoton07.radial = True



