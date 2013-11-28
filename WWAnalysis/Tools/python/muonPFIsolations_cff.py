import FWCore.ParameterSet.Config as cms
import WWAnalysis.Tools.muonPFIsoSingleType_cfi

### DR=0.1
# Charged Hadron isolation
muonPFIsoChHad01 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoChHad01.pfTypes = cms.untracked.vint32(1)
muonPFIsoChHad01.deltaR = 0.1
muonPFIsoChHad01.innerDeltaRVeto = 0.0001

# Neutral Hadron isolation
muonPFIsoNHad01 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoNHad01.pfTypes = cms.untracked.vint32(5)
muonPFIsoNHad01.deltaR = 0.1
muonPFIsoNHad01.innerDeltaRVeto = 0.01

# Photon isolation
muonPFIsoPhoton01 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoPhoton01.pfTypes = cms.untracked.vint32(4)
muonPFIsoPhoton01.deltaR = 0.1
muonPFIsoPhoton01.innerDeltaRVeto = 0.01


### DR=0.2
# Charged Hadron isolation
muonPFIsoChHad02 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoChHad02.pfTypes = cms.untracked.vint32(1)
muonPFIsoChHad02.deltaR = 0.2
muonPFIsoChHad02.innerDeltaRVeto = 0.0001

# Neutral Hadron isolation
muonPFIsoNHad02 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoNHad02.pfTypes = cms.untracked.vint32(5)
muonPFIsoNHad02.deltaR = 0.2
muonPFIsoNHad02.innerDeltaRVeto = 0.01

# Photon isolation
muonPFIsoPhoton02 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoPhoton02.pfTypes = cms.untracked.vint32(4)
muonPFIsoPhoton02.deltaR = 0.2
muonPFIsoPhoton02.innerDeltaRVeto = 0.01


### DR=0.3
# Charged Hadron isolation
muonPFIsoChHad03 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoChHad03.pfTypes = cms.untracked.vint32(1)
muonPFIsoChHad03.deltaR = 0.3
muonPFIsoChHad03.innerDeltaRVeto = 0.0001

# Neutral Hadron isolation
muonPFIsoNHad03 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoNHad03.pfTypes = cms.untracked.vint32(5)
muonPFIsoNHad03.deltaR = 0.3
muonPFIsoNHad03.innerDeltaRVeto = 0.01

# Photon isolation
muonPFIsoPhoton03 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoPhoton03.pfTypes = cms.untracked.vint32(4)
muonPFIsoPhoton03.deltaR = 0.3
muonPFIsoPhoton03.innerDeltaRVeto = 0.01


### DR=0.4
# Charged Hadron isolation
muonPFIsoChHad04 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoChHad04.pfTypes = cms.untracked.vint32(1)
muonPFIsoChHad04.deltaR = 0.4
muonPFIsoChHad04.innerDeltaRVeto = 0.0001

# Neutral Hadron isolation
muonPFIsoNHad04 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoNHad04.pfTypes = cms.untracked.vint32(5)
muonPFIsoNHad04.deltaR = 0.4
muonPFIsoNHad04.innerDeltaRVeto = 0.01

# Photon isolation
muonPFIsoPhoton04 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoPhoton04.pfTypes = cms.untracked.vint32(4)
muonPFIsoPhoton04.deltaR = 0.4
muonPFIsoPhoton04.innerDeltaRVeto = 0.01


### DR=0.5
# Charged Hadron isolation
muonPFIsoChHad05 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoChHad05.pfTypes = cms.untracked.vint32(1)
muonPFIsoChHad05.deltaR = 0.5
muonPFIsoChHad05.innerDeltaRVeto = 0.0001

# Neutral Hadron isolation
muonPFIsoNHad05 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoNHad05.pfTypes = cms.untracked.vint32(5)
muonPFIsoNHad05.deltaR = 0.5
muonPFIsoNHad05.innerDeltaRVeto = 0.01

# Photon isolation
muonPFIsoPhoton05 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoPhoton05.pfTypes = cms.untracked.vint32(4)
muonPFIsoPhoton05.deltaR = 0.5
muonPFIsoPhoton05.innerDeltaRVeto = 0.01


### DR=0.6
# Charged Hadron isolation
muonPFIsoChHad06 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoChHad06.pfTypes = cms.untracked.vint32(1)
muonPFIsoChHad06.deltaR = 0.6
muonPFIsoChHad06.innerDeltaRVeto = 0.0001

# Neutral Hadron isolation
muonPFIsoNHad06 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoNHad06.pfTypes = cms.untracked.vint32(5)
muonPFIsoNHad06.deltaR = 0.6
muonPFIsoNHad06.innerDeltaRVeto = 0.01

# Photon isolation
muonPFIsoPhoton06 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoPhoton06.pfTypes = cms.untracked.vint32(4)
muonPFIsoPhoton06.deltaR = 0.6
muonPFIsoPhoton06.innerDeltaRVeto = 0.01


### DR=0.7
# Charged Hadron isolation
muonPFIsoChHad07 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoChHad07.pfTypes = cms.untracked.vint32(1)
muonPFIsoChHad07.deltaR = 0.7
muonPFIsoChHad07.innerDeltaRVeto = 0.0001

# Neutral Hadron isolation
muonPFIsoNHad07 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoNHad07.pfTypes = cms.untracked.vint32(5)
muonPFIsoNHad07.deltaR = 0.7
muonPFIsoNHad07.innerDeltaRVeto = 0.01

# Photon isolation
muonPFIsoPhoton07 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoPhoton07.pfTypes = cms.untracked.vint32(4)
muonPFIsoPhoton07.deltaR = 0.7
muonPFIsoPhoton07.innerDeltaRVeto = 0.01

####=============== DELTA BETA ====================
### DR=0.3
# Charged Hadron isolation
muonPFIsoChHadPU03 = WWAnalysis.Tools.muonPFIsoSingleType_cfi.muonPFIsoSingleTypeMapProd.clone()
muonPFIsoChHadPU03.pfTypes = cms.untracked.vint32(1)
muonPFIsoChHadPU03.deltaR = 0.3
muonPFIsoChHadPU03.innerDeltaRVeto = 0.01
muonPFIsoChHadPU03.minPt = cms.untracked.double(0.5)
muonPFIsoChHadPU03.pfLabel = 'pfPileUp'

### DR=0.4
# Charged Hadron isolation
muonPFIsoChHadPU04 = muonPFIsoChHadPU03.clone(deltaR = 0.4)

