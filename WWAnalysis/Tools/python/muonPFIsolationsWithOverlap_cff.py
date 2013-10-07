import FWCore.ParameterSet.Config as cms
import WWAnalysis.Tools.muonPFIsoSingleType_cfi

### DR=0.1
# Charged Hadron isolation
muonPFIsoWithOvlChHad01 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoWithOvlChHad01.pfTypes = cms.untracked.vint32(1,2,3)
muonPFIsoWithOvlChHad01.deltaR = 0.1

# Neutral Hadron isolation
muonPFIsoWithOvlNHad01 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoWithOvlNHad01.pfTypes = cms.untracked.vint32(5,2,3)
muonPFIsoWithOvlNHad01.deltaR = 0.1

# Photon isolation
muonPFIsoWithOvlPhoton01 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoWithOvlPhoton01.pfTypes = cms.untracked.vint32(4,2,3)
muonPFIsoWithOvlPhoton01.deltaR = 0.1


### DR=0.2
# Charged Hadron isolation
muonPFIsoWithOvlChHad02 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoWithOvlChHad02.pfTypes = cms.untracked.vint32(1,2,3)
muonPFIsoWithOvlChHad02.deltaR = 0.2

# Neutral Hadron isolation
muonPFIsoWithOvlNHad02 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoWithOvlNHad02.pfTypes = cms.untracked.vint32(5,2,3)
muonPFIsoWithOvlNHad02.deltaR = 0.2

# Photon isolation
muonPFIsoWithOvlPhoton02 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoWithOvlPhoton02.pfTypes = cms.untracked.vint32(4,2,3)
muonPFIsoWithOvlPhoton02.deltaR = 0.2


### DR=0.3
# Charged Hadron isolation
muonPFIsoWithOvlChHad03 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoWithOvlChHad03.pfTypes = cms.untracked.vint32(1,2,3)
muonPFIsoWithOvlChHad03.deltaR = 0.3

# Neutral Hadron isolation
muonPFIsoWithOvlNHad03 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoWithOvlNHad03.pfTypes = cms.untracked.vint32(5,2,3)
muonPFIsoWithOvlNHad03.deltaR = 0.3

# Photon isolation
muonPFIsoWithOvlPhoton03 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoWithOvlPhoton03.pfTypes = cms.untracked.vint32(4,2,3)
muonPFIsoWithOvlPhoton03.deltaR = 0.3


### DR=0.4
# Charged Hadron isolation
muonPFIsoWithOvlChHad04 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoWithOvlChHad04.pfTypes = cms.untracked.vint32(1,2,3)
muonPFIsoWithOvlChHad04.deltaR = 0.4

# Neutral Hadron isolation
muonPFIsoWithOvlNHad04 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoWithOvlNHad04.pfTypes = cms.untracked.vint32(5,2,3)
muonPFIsoWithOvlNHad04.deltaR = 0.4

# Photon isolation
muonPFIsoWithOvlPhoton04 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoWithOvlPhoton04.pfTypes = cms.untracked.vint32(4,2,3)
muonPFIsoWithOvlPhoton04.deltaR = 0.4


### DR=0.5
# Charged Hadron isolation
muonPFIsoWithOvlChHad05 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoWithOvlChHad05.pfTypes = cms.untracked.vint32(1,2,3)
muonPFIsoWithOvlChHad05.deltaR = 0.5

# Neutral Hadron isolation
muonPFIsoWithOvlNHad05 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoWithOvlNHad05.pfTypes = cms.untracked.vint32(5,2,3)
muonPFIsoWithOvlNHad05.deltaR = 0.5

# Photon isolation
muonPFIsoWithOvlPhoton05 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoWithOvlPhoton05.pfTypes = cms.untracked.vint32(4,2,3)
muonPFIsoWithOvlPhoton05.deltaR = 0.5


### DR=0.6
# Charged Hadron isolation
muonPFIsoWithOvlChHad06 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoWithOvlChHad06.pfTypes = cms.untracked.vint32(1,2,3)
muonPFIsoWithOvlChHad06.deltaR = 0.6

# Neutral Hadron isolation
muonPFIsoWithOvlNHad06 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoWithOvlNHad06.pfTypes = cms.untracked.vint32(5,2,3)
muonPFIsoWithOvlNHad06.deltaR = 0.6

# Photon isolation
muonPFIsoWithOvlPhoton06 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoWithOvlPhoton06.pfTypes = cms.untracked.vint32(4,2,3)
muonPFIsoWithOvlPhoton06.deltaR = 0.6


### DR=0.7
# Charged Hadron isolation
muonPFIsoWithOvlChHad07 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoWithOvlChHad07.pfTypes = cms.untracked.vint32(1,2,3)
muonPFIsoWithOvlChHad07.deltaR = 0.7

# Neutral Hadron isolation
muonPFIsoWithOvlNHad07 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoWithOvlNHad07.pfTypes = cms.untracked.vint32(5,2,3)
muonPFIsoWithOvlNHad07.deltaR = 0.7

# Photon isolation
muonPFIsoWithOvlPhoton07 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoWithOvlPhoton07.pfTypes = cms.untracked.vint32(4,2,3)
muonPFIsoWithOvlPhoton07.deltaR = 0.7
