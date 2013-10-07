import FWCore.ParameterSet.Config as cms
from glob import glob

shitSamples = {
    '4150': ['vbfToH150toWWto2L2Nu',0.0002923508],
}

assocProd = {
    '7118': ['wzttH118ToWW',0.0006918764],
    '7120': ['wzttH120ToWW',0.0007238280],
    '7122': ['wzttH122ToWW',0.0008954347],
    '7124': ['wzttH124ToWW',0.0009991516],
    '7126': ['wzttH126ToWW',0.0011041573],
    '7128': ['wzttH128ToWW',0.0012050892],
    '7130': ['wzttH130ToWW',0.0011893260],
    '7135': ['wzttH135ToWW',0.0015132260],
    '7140': ['wzttH140ToWW',0.0015174284],
    '7150': ['wzttH150ToWW',0.0016502763],
    '7160': ['wzttH160ToWW',0.0016590295],
    '7170': ['wzttH170ToWW',0.0014509147],
    '7180': ['wzttH180ToWW',0.0011345486],
    '7190': ['wzttH190ToWW',0.0007897267],
    '7200': ['wzttH200ToWW',0.0006152460],
    '7250': ['wzttH250ToWW',0.0002455993],
    '7300': ['wzttH300ToWW',0.0001154211],
}


sameSignDatasets = {
#     '070': ['WZFull',0.0042670528],
#     '071': ['ZZFull',0.0018300938],
#     '080': ['WJetsToLNuMad',0.3863388729],
    '082': ['WgammaToElNuMad',0.2630577063],
    '083': ['WgammaToMuNuMad',0.2630577063],
    '084': ['WgammaToTauNuMad',0.2630577063],
    '085': ['WGstarToLNu2Mu',0.0160428344],
    '086': ['WGstarToLNu2E' ,0.0208489074],
}

dataSamples = {
    '100': ['SingleElectron2011Av4'],
    '101': ['SingleMuon2011Av4'],
    '102': ['DoubleElectron2011Av4'],
    '103': ['DoubleMuon2011Av4'],
    '104': ['MuEG2011Av4'],
    '120': ['SingleElectron2011Av6'],
    '121': ['SingleMuon2011Av6'],
    '122': ['DoubleElectron2011Av6'],
    '123': ['DoubleMuon2011Av6'],
    '124': ['MuEG2011Av6'],
    '140': ['SingleElectron2011Bv1a'],
    '141': ['SingleMuon2011Bv1a'],
    '142': ['DoubleElectron2011Bv1a'],
    '143': ['DoubleMuon2011Bv1a'],
    '144': ['MuEG2011Bv1a'],
    '150': ['SingleElectron2011AMay10'],
    '151': ['SingleMuon2011AMay10'],
    '152': ['DoubleMuon2011AMay10'],
    '153': ['DoubleElectron2011AMay10'],
    '154': ['MuEG2011AMay10'],
    '160': ['SingleElectron2011AAug05'],
    '161': ['SingleMuon2011AAug05'],
    '162': ['DoubleElectron2011AAug05'],
    '163': ['DoubleMuon2011AAug05'],
    '164': ['MuEG2011AAug05'],
}

allBackgroundSamples = {
    '000': ['WWJets2LMad',0.0040224907],
    '001': ['GluGluToWWTo4L',0.0013949545],
    '011': ['TtWFullDR',0.0096856293],
    '012': ['TbartWFullDR',0.0243926630],
    '013': ['TtFull',0.0107721903],
    '014': ['TbartFull',0.0116785085],
    '015': ['TsFull',0.0122894611],
    '016': ['TbarsFull',0.0104615429],
    '019': ['TTTo2L2Nu2B',0.0016569566],
    '030': ['DYtoElEl',0.0565850053],
    '031': ['DYtoMuMu',0.0561153997],
    '032': ['DYtoTauTau',0.0837185509],
    '033': ['DY10toElEl',0.3389271736],
    '034': ['DY10toMuMu',0.3384486201],
    '035': ['DY10toTauTau',1.5089513601],
    '074': ['WZJetsMad',0.0004893420],
    '071': ['ZZFull',0.0018300938],
    '080': ['WJetsToLNuMad',0.3863388729],
    '082': ['WgammaToElNuMad',0.2630577063],
    '083': ['WgammaToMuNuMad',0.2630577063],
    '084': ['WgammaToTauNuMad',0.2630577063],
    '085': ['WGstarToLNu2Mu',0.0160428344],
    '086': ['WGstarToLNu2E' ,0.0208489074],
}

allOtherBackgrounds = {
    '002': ['WWto2L2NuMCatNLO',0.0088698385],
    '003': ['WWto2L2NuMCatNLOUp',0.0088696742],
    '004': ['WWto2L2NuMCatNLODown',0.0088696413],
    '017': ['TtWFullDS',0.0099167597],
    '018': ['TbartWFullDS',0.0100447226],
    '038': ['ZgammaToElElMad',0.0043042520],
    '039': ['ZgammaToMuMuMad',0.0043354502],
    '040': ['ZgammaToTauTauMad',0.0042718842],
    '041': ['ZgammaToNuNuMad',0.0182223703],
    '037': ['DY50toLLMadD6T',0.0860888229],
    '072': ['GluGluZZ2L2L',0.0000173193],
    '073': ['GluGluZZ4L',0.0000126599],
    '071': ['ZZFull',0.0018300938],
    '075': ['ZZJetsMad',0.0002819586],
    '081': ['VGamma',0.1550354091],
}

allSignalSamples = {
    '1120': ['ggToH120toWWto2L2Nu',0.0010091551],
    '2120': ['ggToH120toWWtoLNuTauNu',0.0010091643],
    '3120': ['ggToH120toWWto2Tau2Nu',0.0005036730],
    '4120': ['vbfToH120toWWto2L2Nu',0.0000766444],
    '5120': ['vbfToH120toWWtoLNuTauNu',0.0000766284],
    '6120': ['vbfToH120toWWto2Tau2Nu',0.0000398429],
    '7120': ['wzttH120ToWW',0.0000759846],
    '1130': ['ggToH130toWWto2L2Nu',0.0018283314],
    '2130': ['ggToH130toWWtoLNuTauNu',0.0018442773],
    '3130': ['ggToH130toWWto2Tau2Nu',0.0009848455],
    '4130': ['vbfToH130toWWto2L2Nu',0.0001486887],
    '5130': ['vbfToH130toWWtoLNuTauNu',0.0001487158],
    '6130': ['vbfToH130toWWto2Tau2Nu',0.0000777812],
    '7130': ['wzttH130ToWW',0.0001248507],
    '1140': ['ggToH140toWWto2L2Nu',0.0026134876],
    '2140': ['ggToH140toWWtoLNuTauNu',0.0025956197],
    '3140': ['ggToH140toWWto2Tau2Nu',0.0014836055],
    '4140': ['vbfToH140toWWto2L2Nu',0.0002236544],
    '5140': ['vbfToH140toWWtoLNuTauNu',0.0002235667],
    '6140': ['vbfToH140toWWto2Tau2Nu',0.0001113159],
    '7140': ['wzttH140ToWW',0.0001592936],
    '1150': ['ggToH150toWWto2L2Nu',0.0031192611],
    '2150': ['ggToH150toWWtoLNuTauNu',0.0031193178],
    '3150': ['ggToH150toWWto2Tau2Nu',0.0015602269],
    '4150': ['vbfToH150toWWto2L2Nu',0.0002923508],
    '5150': ['vbfToH150toWWtoLNuTauNu',0.0002827745],
    '6150': ['vbfToH150toWWto2Tau2Nu',0.0001422942],
    '7150': ['wzttH150ToWW',0.0001732394],
    '1160': ['ggToH160toWWto2L2Nu',0.0036198032],
    '2160': ['ggToH160toWWtoLNuTauNu',0.0035011319],
    '3160': ['ggToH160toWWto2Tau2Nu',0.0017511001],
    '4160': ['vbfToH160toWWto2L2Nu',0.0003365376],
    '5160': ['vbfToH160toWWtoLNuTauNu',0.0003366235],
    '6160': ['vbfToH160toWWto2Tau2Nu',0.0001678863],
    '7160': ['wzttH160ToWW',0.0001741583],
    '1170': ['ggToH170toWWto2L2Nu',0.0032840815],
    '2170': ['ggToH170toWWtoLNuTauNu',0.0031764791],
    '3170': ['ggToH170toWWto2Tau2Nu',0.0017227090],
    '4170': ['vbfToH170toWWto2L2Nu',0.0003319169],
    '5170': ['vbfToH170toWWtoLNuTauNu',0.0003320835],
    '6170': ['vbfToH170toWWto2Tau2Nu',0.0001789436],
    '7170': ['wzttH170ToWW',0.0001523112],
    '1180': ['ggToH180toWWto2L2Nu',0.0026720611],
    '2180': ['ggToH180toWWtoLNuTauNu',0.0026719396],
    '3180': ['ggToH180toWWto2Tau2Nu',0.0013541897],
    '4180': ['vbfToH180toWWto2L2Nu',0.0002945954],
    '5180': ['vbfToH180toWWtoLNuTauNu',0.0002928886],
    '6180': ['vbfToH180toWWto2Tau2Nu',0.0001572636],
    '7180': ['wzttH180ToWW',0.0001191004],
    '1190': ['ggToH190toWWto2L2Nu',0.0020133730],
    '2190': ['ggToH190toWWtoLNuTauNu',0.0019737799],
    '3190': ['ggToH190toWWto2Tau2Nu',0.0009873086],
    '4190': ['vbfToH190toWWto2L2Nu',0.0002303133],
    '5190': ['vbfToH190toWWtoLNuTauNu',0.0002289645],
    '6190': ['vbfToH190toWWto2Tau2Nu',0.0001149111],
    '7190': ['wzttH190ToWW',0.0000829024],
    '1200': ['ggToH200toWWto2L2Nu',0.0016564690],
    '2200': ['ggToH200toWWtoLNuTauNu',0.0016564841],
    '3200': ['ggToH200toWWto2Tau2Nu',0.0008688014],
    '4200': ['vbfToH200toWWto2L2Nu',0.0001988017],
    '5200': ['vbfToH200toWWtoLNuTauNu',0.0002038488],
    '6200': ['vbfToH200toWWto2Tau2Nu',0.0000985437],
    '7200': ['wzttH200ToWW',0.0000645861],
    '1250': ['ggToH250toWWto2L2Nu',0.0010241611],
    '2250': ['ggToH250toWWtoLNuTauNu',0.0010137122],
    '3250': ['ggToH250toWWto2Tau2Nu',0.0005174620],
    '4250': ['vbfToH250toWWto2L2Nu',0.0001267716],
    '5250': ['vbfToH250toWWtoLNuTauNu',0.0001267762],
    '6250': ['vbfToH250toWWto2Tau2Nu',0.0000638278],
    '7250': ['wzttH250ToWW',0.0000257820],
    '1300': ['ggToH300toWWto2L2Nu',0.0007100710],
    '2300': ['ggToH300toWWtoLNuTauNu',0.0007100516],
    '3300': ['ggToH300toWWto2Tau2Nu',0.0003723506],
    '4300': ['vbfToH300toWWto2L2Nu',0.0000875417],
    '5300': ['vbfToH300toWWtoLNuTauNu',0.0000875976],
    '6300': ['vbfToH300toWWto2Tau2Nu',0.0000473345],
    '7300': ['wzttH300ToWW',0.0000121164],
    '1350': ['ggToH350toWWto2L2Nu',0.0006536898],
    '2350': ['ggToH350toWWtoLNuTauNu',0.0006537136],
    '3350': ['ggToH350toWWto2Tau2Nu',0.0003273203],
    '4350': ['vbfToH350toWWto2L2Nu',0.0000611878],
    '5350': ['vbfToH350toWWtoLNuTauNu',0.0000611442],
    '6350': ['vbfToH350toWWto2Tau2Nu',0.0000310497],
    '7350': ['wzttH350ToWW',0],
    '1400': ['ggToH400toWWto2L2Nu',0.0005009410],
    '2400': ['ggToH400toWWtoLNuTauNu',0.0005009319],
    '3400': ['ggToH400toWWto2Tau2Nu',0.0002635247],
    '4400': ['vbfToH400toWWto2L2Nu',0.0000402098],
    '5400': ['vbfToH400toWWtoLNuTauNu',0.0000404553],
    '6400': ['vbfToH400toWWto2Tau2Nu',0.0000200946],
    '7400': ['wzttH400ToWW',0],
    '1450': ['ggToH450toWWto2L2Nu',0.0003200669],
    '2450': ['ggToH450toWWtoLNuTauNu',0.0003200553],
    '3450': ['ggToH450toWWto2Tau2Nu',0.0001715166],
    '4450': ['vbfToH450toWWto2L2Nu',0.0000292606],
    '5450': ['vbfToH450toWWtoLNuTauNu',0.0000292534],
    '6450': ['vbfToH450toWWto2Tau2Nu',0.0000146373],
    '7450': ['wzttH450ToWW',0],
    '1500': ['ggToH500toWWto2L2Nu',0.0002000491],
    '2500': ['ggToH500toWWtoLNuTauNu',0.0002000418],
    '3500': ['ggToH500toWWto2Tau2Nu',0.0001071895],
    '4500': ['vbfToH500toWWto2L2Nu',0.0000220909],
    '5500': ['vbfToH500toWWtoLNuTauNu',0.0000228085],
    '6500': ['vbfToH500toWWto2Tau2Nu',0.0000111392],
    '7500': ['wzttH500ToWW',0],
    '1550': ['ggToH550toWWto2L2Nu',0.0001254853],
    '2550': ['ggToH550toWWtoLNuTauNu',0.0001254762],
    '3550': ['ggToH550toWWto2Tau2Nu',0.0000636584],
    '4550': ['vbfToH550toWWto2L2Nu',0.0000175112],
    '5550': ['vbfToH550toWWtoLNuTauNu',0.0000187858],
    '6550': ['vbfToH550toWWto2Tau2Nu',0.0000091532],
    '7550': ['wzttH550ToWW',0],
    '1600': ['ggToH600toWWto2L2Nu',0.0000791110],
    '2600': ['ggToH600toWWtoLNuTauNu',0.0000791074],
    '3600': ['ggToH600toWWto2Tau2Nu',0.0000428766],
    '4600': ['vbfToH600toWWto2L2Nu',0.0000139997],
    '5600': ['vbfToH600toWWtoLNuTauNu',0.0000138430],
    '6600': ['vbfToH600toWWto2Tau2Nu',0.0000073369],
    '7600': ['wzttH600ToWW',0],
    '7110': ['wzttH110ToWW',0],
    '8110': ['vbfToH110toWWTo2LAndTau2Nu',0.0000717314],
    '9110': ['ggToH110toWWTo2LAndTau2Nu',0.0003412336],
    '7115': ['wzttH115ToWW',0],
    '8115': ['vbfToH115toWWTo2LAndTau2Nu',0.0001229431],
    '9115': ['ggToH115toWWTo2LAndTau2Nu',0.0005497143],
}

wwSamples = {
    '000': ['WWJets2LMad',0.0040224907],
    '001': ['GluGluToWWTo4L',0.0013949545],
}

wwSysSamples = {
    '002': ['WWto2L2NuMCatNLO',0.0088698385],
    '003': ['WWto2L2NuMCatNLOUp',0.0088696742],
    '004': ['WWto2L2NuMCatNLODown',0.0088696413],
}

ggWWSample = {
    '001': ['GluGluToWWTo4L',0.0013949545],
}

ttbarSamples = {
    '019': ['TTTo2L2Nu2B',0.0016569566],
}

topSamples = {
    '011': ['TtWFullDR',0.0096856293],
    '012': ['TbartWFullDR',0.0243926630],
    '013': ['TtFull',0.0107721903],
    '014': ['TbartFull',0.0116785085],
    '015': ['TsFull',0.0122894611],
    '016': ['TbarsFull',0.0104615429],
    '019': ['TTTo2L2Nu2B',0.0016569566],
}

tWSamples = {
    '011': ['TtWFullDR',0.0096856293],
    '012': ['TbartWFullDR',0.0243926630],
    '013': ['TtFull',0.0107721903],
    '014': ['TbartFull',0.0116785085],
    '015': ['TsFull',0.0122894611],
    '016': ['TbarsFull',0.0104615429],
}

tWSysSamples = {
    '017': ['TtWFullDS',0.0099167597],
    '018': ['TbartWFullDS',0.0100447226],
}

dySamples = {
    '030': ['DYtoElEl',0.0565850053],
    '031': ['DYtoMuMu',0.0561153997],
    '032': ['DYtoTauTau',0.0837185509],
    '033': ['DY10toElEl',0.3389271736],
    '034': ['DY10toMuMu',0.3384486201],
    '035': ['DY10toTauTau',1.5089513601],
}

zGamSamples = {
    '038': ['ZgammaToElElMad',0.0043042520],
    '039': ['ZgammaToMuMuMad',0.0043354502],
    '040': ['ZgammaToTauTauMad',0.0042718842],
    '041': ['ZgammaToNuNuMad',0.0182223703],
}

dyMadSamples = {
    '037': ['DY50toLLMadD6T',0.0860888229],
}

wzSamples = {
    '074': ['WZJetsMad',0.0004893420],
}

zzSamples = {
    '071': ['ZZFull',0.0018300938],
}

diBosonSamples = {
    '074': ['WZJetsMad',0.0004893420],
    '071': ['ZZFull',0.0018300938],
}

diBosonSysSamples = {
    '072': ['GluGluZZ2L2L',0.0000173193],
    '073': ['GluGluZZ4L',0.0000126599],
    '071': ['ZZFull',0.0018300938],
    '075': ['ZZJetsMad',0.0002819586],
}

wJetSamples = {
    '080': ['WJetsToLNuMad',0.3863388729],
}

wGamSamples = {
    '082': ['WgammaToElNuMad',0.2630577063],
    '083': ['WgammaToMuNuMad',0.2630577063],
    '084': ['WgammaToTauNuMad',0.2630577063],
}


vGamSamples = {
    '081': ['VGamma',0.1550354091],
}

wGamStar = {
    '085': ['WGstarToLNu2Mu',0.0160428344],
    '086': ['WGstarToLNu2E' ,0.0208489074],
}


h1Samples = {
    '1120': ['ggToH120toWWto2L2Nu',0.0010091551],
    '1130': ['ggToH130toWWto2L2Nu',0.0018283314],
    '1140': ['ggToH140toWWto2L2Nu',0.0026134876],
    '1150': ['ggToH150toWWto2L2Nu',0.0031192611],
    '1160': ['ggToH160toWWto2L2Nu',0.0036198032],
    '1170': ['ggToH170toWWto2L2Nu',0.0032840815],
    '1180': ['ggToH180toWWto2L2Nu',0.0026720611],
    '1190': ['ggToH190toWWto2L2Nu',0.0020133730],
    '1200': ['ggToH200toWWto2L2Nu',0.0016564690],
    '1250': ['ggToH250toWWto2L2Nu',0.0010241611],
    '1300': ['ggToH300toWWto2L2Nu',0.0007100710],
    '1350': ['ggToH350toWWto2L2Nu',0.0006536898],
    '1400': ['ggToH400toWWto2L2Nu',0.0005009410],
    '1450': ['ggToH450toWWto2L2Nu',0.0003200669],
    '1500': ['ggToH500toWWto2L2Nu',0.0002000491],
    '1550': ['ggToH550toWWto2L2Nu',0.0001254853],
    '1600': ['ggToH600toWWto2L2Nu',0.0000791110],
}

h2Samples = {
    '2120': ['ggToH120toWWtoLNuTauNu',0.0010091643],
    '2130': ['ggToH130toWWtoLNuTauNu',0.0018442773],
    '2140': ['ggToH140toWWtoLNuTauNu',0.0025956197],
    '2150': ['ggToH150toWWtoLNuTauNu',0.0031193178],
    '2160': ['ggToH160toWWtoLNuTauNu',0.0035011319],
    '2170': ['ggToH170toWWtoLNuTauNu',0.0031764791],
    '2180': ['ggToH180toWWtoLNuTauNu',0.0026719396],
    '2190': ['ggToH190toWWtoLNuTauNu',0.0019737799],
    '2200': ['ggToH200toWWtoLNuTauNu',0.0016564841],
    '2250': ['ggToH250toWWtoLNuTauNu',0.0010137122],
    '2300': ['ggToH300toWWtoLNuTauNu',0.0007100516],
    '2350': ['ggToH350toWWtoLNuTauNu',0.0006537136],
    '2400': ['ggToH400toWWtoLNuTauNu',0.0005009319],
    '2450': ['ggToH450toWWtoLNuTauNu',0.0003200553],
    '2500': ['ggToH500toWWtoLNuTauNu',0.0002000418],
    '2550': ['ggToH550toWWtoLNuTauNu',0.0001254762],
    '2600': ['ggToH600toWWtoLNuTauNu',0.0000791074],
}

h3Samples = {
    '3120': ['ggToH120toWWto2Tau2Nu',0.0005036730],
    '3130': ['ggToH130toWWto2Tau2Nu',0.0009848455],
    '3140': ['ggToH140toWWto2Tau2Nu',0.0014836055],
    '3150': ['ggToH150toWWto2Tau2Nu',0.0015602269],
    '3160': ['ggToH160toWWto2Tau2Nu',0.0017511001],
    '3170': ['ggToH170toWWto2Tau2Nu',0.0017227090],
    '3180': ['ggToH180toWWto2Tau2Nu',0.0013541897],
    '3190': ['ggToH190toWWto2Tau2Nu',0.0009873086],
    '3200': ['ggToH200toWWto2Tau2Nu',0.0008688014],
    '3250': ['ggToH250toWWto2Tau2Nu',0.0005174620],
    '3300': ['ggToH300toWWto2Tau2Nu',0.0003723506],
    '3350': ['ggToH350toWWto2Tau2Nu',0.0003273203],
    '3400': ['ggToH400toWWto2Tau2Nu',0.0002635247],
    '3450': ['ggToH450toWWto2Tau2Nu',0.0001715166],
    '3500': ['ggToH500toWWto2Tau2Nu',0.0001071895],
    '3550': ['ggToH550toWWto2Tau2Nu',0.0000636584],
    '3600': ['ggToH600toWWto2Tau2Nu',0.0000428766],
}

h4Samples = {
    '4120': ['vbfToH120toWWto2L2Nu',0.0000766444],
    '4130': ['vbfToH130toWWto2L2Nu',0.0001486887],
    '4140': ['vbfToH140toWWto2L2Nu',0.0002236544],
    '4150': ['vbfToH150toWWto2L2Nu',0.0002923508],
    '4160': ['vbfToH160toWWto2L2Nu',0.0003365376],
    '4170': ['vbfToH170toWWto2L2Nu',0.0003319169],
    '4180': ['vbfToH180toWWto2L2Nu',0.0002945954],
    '4190': ['vbfToH190toWWto2L2Nu',0.0002303133],
    '4200': ['vbfToH200toWWto2L2Nu',0.0001988017],
    '4250': ['vbfToH250toWWto2L2Nu',0.0001267716],
    '4300': ['vbfToH300toWWto2L2Nu',0.0000875417],
    '4350': ['vbfToH350toWWto2L2Nu',0.0000611878],
    '4400': ['vbfToH400toWWto2L2Nu',0.0000402098],
    '4450': ['vbfToH450toWWto2L2Nu',0.0000292606],
    '4500': ['vbfToH500toWWto2L2Nu',0.0000220909],
    '4550': ['vbfToH550toWWto2L2Nu',0.0000175112],
    '4600': ['vbfToH600toWWto2L2Nu',0.0000139997],
}

h5Samples = {
    '5120': ['vbfToH120toWWtoLNuTauNu',0.0000766284],
    '5130': ['vbfToH130toWWtoLNuTauNu',0.0001487158],
    '5140': ['vbfToH140toWWtoLNuTauNu',0.0002235667],
    '5150': ['vbfToH150toWWtoLNuTauNu',0.0002827745],
    '5160': ['vbfToH160toWWtoLNuTauNu',0.0003366235],
    '5170': ['vbfToH170toWWtoLNuTauNu',0.0003320835],
    '5180': ['vbfToH180toWWtoLNuTauNu',0.0002928886],
    '5190': ['vbfToH190toWWtoLNuTauNu',0.0002289645],
    '5200': ['vbfToH200toWWtoLNuTauNu',0.0002038488],
    '5250': ['vbfToH250toWWtoLNuTauNu',0.0001267762],
    '5300': ['vbfToH300toWWtoLNuTauNu',0.0000875976],
    '5350': ['vbfToH350toWWtoLNuTauNu',0.0000611442],
    '5400': ['vbfToH400toWWtoLNuTauNu',0.0000404553],
    '5450': ['vbfToH450toWWtoLNuTauNu',0.0000292534],
    '5500': ['vbfToH500toWWtoLNuTauNu',0.0000228085],
    '5550': ['vbfToH550toWWtoLNuTauNu',0.0000187858],
    '5600': ['vbfToH600toWWtoLNuTauNu',0.0000138430],
}

h6Samples = {
    '6120': ['vbfToH120toWWto2Tau2Nu',0.0000398429],
    '6130': ['vbfToH130toWWto2Tau2Nu',0.0000777812],
    '6140': ['vbfToH140toWWto2Tau2Nu',0.0001113159],
    '6150': ['vbfToH150toWWto2Tau2Nu',0.0001422942],
    '6160': ['vbfToH160toWWto2Tau2Nu',0.0001678863],
    '6170': ['vbfToH170toWWto2Tau2Nu',0.0001789436],
    '6180': ['vbfToH180toWWto2Tau2Nu',0.0001572636],
    '6190': ['vbfToH190toWWto2Tau2Nu',0.0001149111],
    '6200': ['vbfToH200toWWto2Tau2Nu',0.0000985437],
    '6250': ['vbfToH250toWWto2Tau2Nu',0.0000638278],
    '6300': ['vbfToH300toWWto2Tau2Nu',0.0000473345],
    '6350': ['vbfToH350toWWto2Tau2Nu',0.0000310497],
    '6400': ['vbfToH400toWWto2Tau2Nu',0.0000200946],
    '6450': ['vbfToH450toWWto2Tau2Nu',0.0000146373],
    '6500': ['vbfToH500toWWto2Tau2Nu',0.0000111392],
    '6550': ['vbfToH550toWWto2Tau2Nu',0.0000091532],
    '6600': ['vbfToH600toWWto2Tau2Nu',0.0000073369],
}

h7Samples = {
    '7110': ['wzttH110ToWW',0],
    '7115': ['wzttH115ToWW',0],
    '7120': ['wzttH120ToWW',0.0000759846],
    '7130': ['wzttH130ToWW',0.0001248507],
    '7140': ['wzttH140ToWW',0.0001592936],
    '7150': ['wzttH150ToWW',0.0001732394],
    '7160': ['wzttH160ToWW',0.0001741583],
    '7170': ['wzttH170ToWW',0.0001523112],
    '7180': ['wzttH180ToWW',0.0001191004],
    '7190': ['wzttH190ToWW',0.0000829024],
    '7200': ['wzttH200ToWW',0.0000645861],
    '7250': ['wzttH250ToWW',0.0000257820],
    '7300': ['wzttH300ToWW',0.0000121164],
    '7350': ['wzttH350ToWW',0],
    '7400': ['wzttH400ToWW',0],
    '7450': ['wzttH450ToWW',0],
    '7500': ['wzttH500ToWW',0],
    '7550': ['wzttH550ToWW',0],
    '7600': ['wzttH600ToWW',0],
}

h8Samples = {
    '8110': ['vbfToH110toWWTo2LAndTau2Nu',0.0000717314],
    '8115': ['vbfToH115toWWTo2LAndTau2Nu',0.0001229431],
}

h9Samples = {
    '9110': ['ggToH110toWWTo2LAndTau2Nu',0.0003412336],
    '9115': ['ggToH115toWWTo2LAndTau2Nu',0.0005497143],
}

h120Samples = {
    '1120': ['ggToH120toWWto2L2Nu',0.0010091551],
    '2120': ['ggToH120toWWtoLNuTauNu',0.0010091643],
    '3120': ['ggToH120toWWto2Tau2Nu',0.0005036730],
    '4120': ['vbfToH120toWWto2L2Nu',0.0000766444],
    '5120': ['vbfToH120toWWtoLNuTauNu',0.0000766284],
    '6120': ['vbfToH120toWWto2Tau2Nu',0.0000398429],
    '7120': ['wzttH120ToWW',0.0000759846],
}

h130Samples = {
    '1130': ['ggToH130toWWto2L2Nu',0.0018283314],
    '2130': ['ggToH130toWWtoLNuTauNu',0.0018442773],
    '3130': ['ggToH130toWWto2Tau2Nu',0.0009848455],
    '4130': ['vbfToH130toWWto2L2Nu',0.0001486887],
    '5130': ['vbfToH130toWWtoLNuTauNu',0.0001487158],
    '6130': ['vbfToH130toWWto2Tau2Nu',0.0000777812],
    '7130': ['wzttH130ToWW',0.0001248507],
}

h140Samples = {
    '1140': ['ggToH140toWWto2L2Nu',0.0026134876],
    '2140': ['ggToH140toWWtoLNuTauNu',0.0025956197],
    '3140': ['ggToH140toWWto2Tau2Nu',0.0014836055],
    '4140': ['vbfToH140toWWto2L2Nu',0.0002236544],
    '5140': ['vbfToH140toWWtoLNuTauNu',0.0002235667],
    '6140': ['vbfToH140toWWto2Tau2Nu',0.0001113159],
    '7140': ['wzttH140ToWW',0.0001592936],
}

h150Samples = {
    '1150': ['ggToH150toWWto2L2Nu',0.0031192611],
    '2150': ['ggToH150toWWtoLNuTauNu',0.0031193178],
    '3150': ['ggToH150toWWto2Tau2Nu',0.0015602269],
    '4150': ['vbfToH150toWWto2L2Nu',0.0002923508],
    '5150': ['vbfToH150toWWtoLNuTauNu',0.0002827745],
    '6150': ['vbfToH150toWWto2Tau2Nu',0.0001422942],
    '7150': ['wzttH150ToWW',0.0001732394],
}

h160Samples = {
    '1160': ['ggToH160toWWto2L2Nu',0.0036198032],
    '2160': ['ggToH160toWWtoLNuTauNu',0.0035011319],
    '3160': ['ggToH160toWWto2Tau2Nu',0.0017511001],
    '4160': ['vbfToH160toWWto2L2Nu',0.0003365376],
    '5160': ['vbfToH160toWWtoLNuTauNu',0.0003366235],
    '6160': ['vbfToH160toWWto2Tau2Nu',0.0001678863],
    '7160': ['wzttH160ToWW',0.0001741583],
}

h170Samples = {
    '1170': ['ggToH170toWWto2L2Nu',0.0032840815],
    '2170': ['ggToH170toWWtoLNuTauNu',0.0031764791],
    '3170': ['ggToH170toWWto2Tau2Nu',0.0017227090],
    '4170': ['vbfToH170toWWto2L2Nu',0.0003319169],
    '5170': ['vbfToH170toWWtoLNuTauNu',0.0003320835],
    '6170': ['vbfToH170toWWto2Tau2Nu',0.0001789436],
    '7170': ['wzttH170ToWW',0.0001523112],
}

h180Samples = {
    '1180': ['ggToH180toWWto2L2Nu',0.0026720611],
    '2180': ['ggToH180toWWtoLNuTauNu',0.0026719396],
    '3180': ['ggToH180toWWto2Tau2Nu',0.0013541897],
    '4180': ['vbfToH180toWWto2L2Nu',0.0002945954],
    '5180': ['vbfToH180toWWtoLNuTauNu',0.0002928886],
    '6180': ['vbfToH180toWWto2Tau2Nu',0.0001572636],
    '7180': ['wzttH180ToWW',0.0001191004],
}

h190Samples = {
    '1190': ['ggToH190toWWto2L2Nu',0.0020133730],
    '2190': ['ggToH190toWWtoLNuTauNu',0.0019737799],
    '3190': ['ggToH190toWWto2Tau2Nu',0.0009873086],
    '4190': ['vbfToH190toWWto2L2Nu',0.0002303133],
    '5190': ['vbfToH190toWWtoLNuTauNu',0.0002289645],
    '6190': ['vbfToH190toWWto2Tau2Nu',0.0001149111],
    '7190': ['wzttH190ToWW',0.0000829024],
}

h200Samples = {
    '1200': ['ggToH200toWWto2L2Nu',0.0016564690],
    '2200': ['ggToH200toWWtoLNuTauNu',0.0016564841],
    '3200': ['ggToH200toWWto2Tau2Nu',0.0008688014],
    '4200': ['vbfToH200toWWto2L2Nu',0.0001988017],
    '5200': ['vbfToH200toWWtoLNuTauNu',0.0002038488],
    '6200': ['vbfToH200toWWto2Tau2Nu',0.0000985437],
    '7200': ['wzttH200ToWW',0.0000645861],
}

h250Samples = {
    '1250': ['ggToH250toWWto2L2Nu',0.0010241611],
    '2250': ['ggToH250toWWtoLNuTauNu',0.0010137122],
    '3250': ['ggToH250toWWto2Tau2Nu',0.0005174620],
    '4250': ['vbfToH250toWWto2L2Nu',0.0001267716],
    '5250': ['vbfToH250toWWtoLNuTauNu',0.0001267762],
    '6250': ['vbfToH250toWWto2Tau2Nu',0.0000638278],
    '7250': ['wzttH250ToWW',0.0000257820],
}

h300Samples = {
    '1300': ['ggToH300toWWto2L2Nu',0.0007100710],
    '2300': ['ggToH300toWWtoLNuTauNu',0.0007100516],
    '3300': ['ggToH300toWWto2Tau2Nu',0.0003723506],
    '4300': ['vbfToH300toWWto2L2Nu',0.0000875417],
    '5300': ['vbfToH300toWWtoLNuTauNu',0.0000875976],
    '6300': ['vbfToH300toWWto2Tau2Nu',0.0000473345],
    '7300': ['wzttH300ToWW',0.0000121164],
}

h350Samples = {
    '1350': ['ggToH350toWWto2L2Nu',0.0006536898],
    '2350': ['ggToH350toWWtoLNuTauNu',0.0006537136],
    '3350': ['ggToH350toWWto2Tau2Nu',0.0003273203],
    '4350': ['vbfToH350toWWto2L2Nu',0.0000611878],
    '5350': ['vbfToH350toWWtoLNuTauNu',0.0000611442],
    '6350': ['vbfToH350toWWto2Tau2Nu',0.0000310497],
    '7350': ['wzttH350ToWW',0],
}

h400Samples = {
    '1400': ['ggToH400toWWto2L2Nu',0.0005009410],
    '2400': ['ggToH400toWWtoLNuTauNu',0.0005009319],
    '3400': ['ggToH400toWWto2Tau2Nu',0.0002635247],
    '4400': ['vbfToH400toWWto2L2Nu',0.0000402098],
    '5400': ['vbfToH400toWWtoLNuTauNu',0.0000404553],
    '6400': ['vbfToH400toWWto2Tau2Nu',0.0000200946],
    '7400': ['wzttH400ToWW',0],
}

h450Samples = {
    '1450': ['ggToH450toWWto2L2Nu',0.0003200669],
    '2450': ['ggToH450toWWtoLNuTauNu',0.0003200553],
    '3450': ['ggToH450toWWto2Tau2Nu',0.0001715166],
    '4450': ['vbfToH450toWWto2L2Nu',0.0000292606],
    '5450': ['vbfToH450toWWtoLNuTauNu',0.0000292534],
    '6450': ['vbfToH450toWWto2Tau2Nu',0.0000146373],
    '7450': ['wzttH450ToWW',0],
}

h500Samples = {
    '1500': ['ggToH500toWWto2L2Nu',0.0002000491],
    '2500': ['ggToH500toWWtoLNuTauNu',0.0002000418],
    '3500': ['ggToH500toWWto2Tau2Nu',0.0001071895],
    '4500': ['vbfToH500toWWto2L2Nu',0.0000220909],
    '5500': ['vbfToH500toWWtoLNuTauNu',0.0000228085],
    '6500': ['vbfToH500toWWto2Tau2Nu',0.0000111392],
    '7500': ['wzttH500ToWW',0],
}

h550Samples = {
    '1550': ['ggToH550toWWto2L2Nu',0.0001254853],
    '2550': ['ggToH550toWWtoLNuTauNu',0.0001254762],
    '3550': ['ggToH550toWWto2Tau2Nu',0.0000636584],
    '4550': ['vbfToH550toWWto2L2Nu',0.0000175112],
    '5550': ['vbfToH550toWWtoLNuTauNu',0.0000187858],
    '6550': ['vbfToH550toWWto2Tau2Nu',0.0000091532],
    '7550': ['wzttH550ToWW',0],
}

h600Samples = {
    '1600': ['ggToH600toWWto2L2Nu',0.0000791110],
    '2600': ['ggToH600toWWtoLNuTauNu',0.0000791074],
    '3600': ['ggToH600toWWto2Tau2Nu',0.0000428766],
    '4600': ['vbfToH600toWWto2L2Nu',0.0000139997],
    '5600': ['vbfToH600toWWtoLNuTauNu',0.0000138430],
    '6600': ['vbfToH600toWWto2Tau2Nu',0.0000073369],
    '7600': ['wzttH600ToWW',0],
}

h115Samples = {
    '7115': ['wzttH115ToWW',0],
    '8115': ['vbfToH115toWWTo2LAndTau2Nu',0.0001229431],
    '9115': ['ggToH115toWWTo2LAndTau2Nu',0.0005497143],
}

h110Samples = {
    '7110': ['wzttH110ToWW',0],
    '8110': ['vbfToH110toWWTo2LAndTau2Nu',0.0000717314],
    '9110': ['ggToH110toWWTo2LAndTau2Nu',0.0003412336],
}

