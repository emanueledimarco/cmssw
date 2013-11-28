import FWCore.ParameterSet.Config as cms
import WWAnalysis.Tools.electronPFIsoSingleType_cfi

### DR=0.1 cone
# Charged Hadron isolation
electronPFIsoWithOvlChHad01 = WWAnalysis.Tools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoWithOvlChHad01.pfTypes = cms.untracked.vint32(1,2,3)
electronPFIsoWithOvlChHad01.deltaR = 0.1

# Neutral Hadron isolation
electronPFIsoWithOvlNHad01 = WWAnalysis.Tools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoWithOvlNHad01.pfTypes = cms.untracked.vint32(5,2,3)
electronPFIsoWithOvlNHad01.deltaR = 0.1

# Photon isolation
electronPFIsoWithOvlPhoton01 = WWAnalysis.Tools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoWithOvlPhoton01.pfTypes = cms.untracked.vint32(4,2,3)
electronPFIsoWithOvlPhoton01.deltaR = 0.1


### DR=0.2 cone
# Charged Hadron isolation
electronPFIsoWithOvlChHad02 = WWAnalysis.Tools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoWithOvlChHad02.pfTypes = cms.untracked.vint32(1,2,3)
electronPFIsoWithOvlChHad02.deltaR = 0.2

# Neutral Hadron isolation
electronPFIsoWithOvlNHad02 = WWAnalysis.Tools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoWithOvlNHad02.pfTypes = cms.untracked.vint32(5,2,3)
electronPFIsoWithOvlNHad02.deltaR = 0.2

# Photon isolation
electronPFIsoWithOvlPhoton02 = WWAnalysis.Tools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoWithOvlPhoton02.pfTypes = cms.untracked.vint32(4,2,3)
electronPFIsoWithOvlPhoton02.deltaR = 0.2


### DR=0.3 cone
# Charged Hadron isolation
electronPFIsoWithOvlChHad03 = WWAnalysis.Tools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoWithOvlChHad03.pfTypes = cms.untracked.vint32(1,2,3)
electronPFIsoWithOvlChHad03.deltaR = 0.3

# Neutral Hadron isolation
electronPFIsoWithOvlNHad03 = WWAnalysis.Tools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoWithOvlNHad03.pfTypes = cms.untracked.vint32(5,2,3)
electronPFIsoWithOvlNHad03.deltaR = 0.3

# Photon isolation
electronPFIsoWithOvlPhoton03 = WWAnalysis.Tools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoWithOvlPhoton03.pfTypes = cms.untracked.vint32(4,2,3)
electronPFIsoWithOvlPhoton03.deltaR = 0.3


### DR=0.4 cone
# Charged Hadron isolation
electronPFIsoWithOvlChHad04 = WWAnalysis.Tools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoWithOvlChHad04.pfTypes = cms.untracked.vint32(1,2,3)
electronPFIsoWithOvlChHad04.deltaR = 0.4

# Neutral Hadron isolation
electronPFIsoWithOvlNHad04 = WWAnalysis.Tools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoWithOvlNHad04.pfTypes = cms.untracked.vint32(5,2,3)
electronPFIsoWithOvlNHad04.deltaR = 0.4

# Photon isolation
electronPFIsoWithOvlPhoton04 = WWAnalysis.Tools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoWithOvlPhoton04.pfTypes = cms.untracked.vint32(4,2,3)
electronPFIsoWithOvlPhoton04.deltaR = 0.4


### DR=0.5 cone
# Charged Hadron isolation
electronPFIsoWithOvlChHad05 = WWAnalysis.Tools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoWithOvlChHad05.pfTypes = cms.untracked.vint32(1,2,3)
electronPFIsoWithOvlChHad05.deltaR = 0.5

# Neutral Hadron isolation
electronPFIsoWithOvlNHad05 = WWAnalysis.Tools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoWithOvlNHad05.pfTypes = cms.untracked.vint32(5,2,3)
electronPFIsoWithOvlNHad05.deltaR = 0.5

# Photon isolation
electronPFIsoWithOvlPhoton05 = WWAnalysis.Tools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoWithOvlPhoton05.pfTypes = cms.untracked.vint32(4,2,3)
electronPFIsoWithOvlPhoton05.deltaR = 0.5


### DR=0.6 cone
# Charged Hadron isolation
electronPFIsoWithOvlChHad06 = WWAnalysis.Tools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoWithOvlChHad06.pfTypes = cms.untracked.vint32(1,2,3)
electronPFIsoWithOvlChHad06.deltaR = 0.6

# Neutral Hadron isolation
electronPFIsoWithOvlNHad06 = WWAnalysis.Tools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoWithOvlNHad06.pfTypes = cms.untracked.vint32(5,2,3)
electronPFIsoWithOvlNHad06.deltaR = 0.6

# Photon isolation
electronPFIsoWithOvlPhoton06 = WWAnalysis.Tools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoWithOvlPhoton06.pfTypes = cms.untracked.vint32(4,2,3)
electronPFIsoWithOvlPhoton06.deltaR = 0.6


### DR=0.7 cone
# Charged Hadron isolation
electronPFIsoWithOvlChHad07 = WWAnalysis.Tools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoWithOvlChHad07.pfTypes = cms.untracked.vint32(1,2,3)
electronPFIsoWithOvlChHad07.deltaR = 0.7

# Neutral Hadron isolation
electronPFIsoWithOvlNHad07 = WWAnalysis.Tools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoWithOvlNHad07.pfTypes = cms.untracked.vint32(5,2,3)
electronPFIsoWithOvlNHad07.deltaR = 0.7

# Photon isolation
electronPFIsoWithOvlPhoton07 = WWAnalysis.Tools.electronPFIsoSingleType_cfi.electronPFIsoSingleTypeMapProd.clone()
electronPFIsoWithOvlPhoton07.pfTypes = cms.untracked.vint32(4,2,3)
electronPFIsoWithOvlPhoton07.deltaR = 0.7
