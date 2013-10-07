import FWCore.ParameterSet.Config as cms
from glob import glob

# Prefix where we expect idXX/*.root files to be
prefix = '/nfs/bluearc/group/trees/ww/V03_WWEventSkimProd_Newest/'

higgsSamples = {
    # id          name                       scale to 1/pb   isMC       Enhancements
    'id101120': ['GluGluToHToWWTo2L2Nu_M120',   0.000001048835009 ,  False, [], 5.753E+00 , True],
    'id101130': ['GluGluToHToWWTo2L2Nu_M130',   0.000001829815319 ,  False, [], 5.863E+00 , True],
    'id101140': ['GluGluToHToWWTo2L2Nu_M140',   0.000002703063408 ,  False, [], 6.430E+00 , True],
    'id101150': ['GluGluToHToWWTo2L2Nu_M150',   0.00000325540878  ,  False, [], 7.245E+00 , True],
    'id101160': ['GluGluToHToWWTo2L2Nu_M160',   0.000003777578691 ,  False, [], 8.405E+00 , True],
    'id101170': ['GluGluToHToWWTo2L2Nu_M170',   0.000003367876121 ,  False, [], 8.662E+00 , True],
    'id101180': ['GluGluToHToWWTo2L2Nu_M180',   0.000002790900049 ,  False, [], 8.610E+00 , True],
    'id101190': ['GluGluToHToWWTo2L2Nu_M190',   0.000002137244856 ,  False, [], 8.373E+00 , True],
    'id101200': ['GluGluToHToWWTo2L2Nu_M200',   0.00000170720294  ,  False, [], 8.289E+00 , True],
    'id101210': ['GluGluToHToWWTo2L2Nu_M210',   0.000002497535769 ,  False, [], 8.1822764 , True],
    'id101220': ['GluGluToHToWWTo2L2Nu_M220',   0.000002234705767 ,  False, [], 8.0458019 , True],
    'id101230': ['GluGluToHToWWTo2L2Nu_M230',   0.000002022070138 ,  False, [], 7.9407256 , True],
    'id101250': ['GluGluToHToWWTo2L2Nu_M250',   0.00000101725358  ,  False, [], 7.655E+00 , True],
    'id101300': ['GluGluToHToWWTo2L2Nu_M300',   0.00000073511357  ,  False, [], 6.795E+00 , True],
    'id101350': ['GluGluToHToWWTo2L2Nu_M350',   0.000000708865814 ,  False, [], 5.547E+00 , True],
    'id101400': ['GluGluToHToWWTo2L2Nu_M400',   0.000000507398744 ,  False, [], 4.429E+00 , True],
    'id101450': ['GluGluToHToWWTo2L2Nu_M450',   0.000000313932996 ,  False, [], 4.195E+00 , True],
    'id101500': ['GluGluToHToWWTo2L2Nu_M500',   0.000000193283112 ,  False, [], 4.167E+00 , True],
    'id101550': ['GluGluToHToWWTo2L2Nu_M550',   0.000000120214038 ,  False, [], 4.228E+00 , True],
    'id101600': ['GluGluToHToWWTo2L2Nu_M600',   0.000000075779319 ,  False, [], 4.350E+00 , True],
    'id102120': ['GluGluToHToWWToLNuTauNu_M120',0.000000998666618 ,  False, [], 5.753E+00 , True],
    'id102130': ['GluGluToHToWWToLNuTauNu_M130',0.000001829881863 ,  False, [], 5.863E+00 , True],
    'id102140': ['GluGluToHToWWToLNuTauNu_M140',0.000002623297739 ,  False, [], 6.430E+00 , True],
    'id102150': ['GluGluToHToWWToLNuTauNu_M150',0.000003183605691 ,  False, [], 7.245E+00 , True],
    'id102160': ['GluGluToHToWWToLNuTauNu_M160',0.000003624777523 ,  False, [], 8.405E+00 , True],
    'id102170': ['GluGluToHToWWToLNuTauNu_M170',0.000003268399498 ,  False, [], 8.662E+00 , True],
    'id102180': ['GluGluToHToWWToLNuTauNu_M180',0.000002770448318 ,  False, [], 8.610E+00 , True],
    'id102190': ['GluGluToHToWWToLNuTauNu_M190',0.000002035089122 ,  False, [], 8.373E+00 , True],
    'id102200': ['GluGluToHToWWToLNuTauNu_M200',0.000001707125335 ,  False, [], 8.289E+00 , True],
    'id102210': ['GluGluToHToWWToLNuTauNu_M210',0.000002497384396 ,  False, [], 8.1822764 , True],
    'id102220': ['GluGluToHToWWToLNuTauNu_M220',0.000002234671905 ,  False, [], 8.0458019 , True],
    'id102230': ['GluGluToHToWWToLNuTauNu_M230',0.000002022223342 ,  False, [], 7.9407256 , True],
    'id102250': ['GluGluToHToWWToLNuTauNu_M250',0.000001017299822 ,  False, [], 7.655E+00 , True],
    'id102300': ['GluGluToHToWWToLNuTauNu_M300',0.000000735126937 ,  False, [], 6.795E+00 , True],
    'id102350': ['GluGluToHToWWToLNuTauNu_M350',0.000000814085828 ,  False, [], 5.547E+00 , True],
    'id102400': ['GluGluToHToWWToLNuTauNu_M400',0.000000505525811 ,  False, [], 4.429E+00 , True],
    'id102450': ['GluGluToHToWWToLNuTauNu_M450',0.000000313935851 ,  False, [], 4.195E+00 , True],
    'id102500': ['GluGluToHToWWToLNuTauNu_M500',0.000000193272568 ,  False, [], 4.167E+00 , True],
    'id102550': ['GluGluToHToWWToLNuTauNu_M550',0.000000120207479 ,  False, [], 4.228E+00 , True],
    'id102600': ['GluGluToHToWWToLNuTauNu_M600',0.00000007577794  ,  False, [], 4.350E+00 , True],
    'id103120': ['GluGluToHToWWTo2Tau2Nu_M120', 0.000000416090915 ,  False, [], 5.753E+00 , True],
    'id103130': ['GluGluToHToWWTo2Tau2Nu_M130', 0.000000762446155 ,  False, [], 5.863E+00 , True],
    'id103140': ['GluGluToHToWWTo2Tau2Nu_M140', 0.000001093130163 ,  False, [], 6.430E+00 , True],
    'id103150': ['GluGluToHToWWTo2Tau2Nu_M150', 0.00000132646218  ,  False, [], 7.245E+00 , True],
    'id103160': ['GluGluToHToWWTo2Tau2Nu_M160', 0.000001510241587 ,  False, [], 8.405E+00 , True],
    'id103170': ['GluGluToHToWWTo2Tau2Nu_M170', 0.000001361911539 ,  False, [], 8.662E+00 , True],
    'id103180': ['GluGluToHToWWTo2Tau2Nu_M180', 0.000001154325479 ,  False, [], 8.610E+00 , True],
    'id103190': ['GluGluToHToWWTo2Tau2Nu_M190', 0.000000847894704 ,  False, [], 8.373E+00 , True],
    'id103200': ['GluGluToHToWWTo2Tau2Nu_M200', 0.00000071132162  ,  False, [], 8.289E+00 , True],
    'id103210': ['GluGluToHToWWTo2Tau2Nu_M210', 0.000000624336647 ,  False, [], 8.1822764 , True],
    'id103220': ['GluGluToHToWWTo2Tau2Nu_M220', 0.000000558651049 ,  False, [], 8.0458019 , True],
    'id103230': ['GluGluToHToWWTo2Tau2Nu_M230', 0.00000050553285  ,  False, [], 7.9407256 , True],
    'id103250': ['GluGluToHToWWTo2Tau2Nu_M250', 0.000000423914759 ,  False, [], 7.655E+00 , True],
    'id103300': ['GluGluToHToWWTo2Tau2Nu_M300', 0.000000306308467 ,  False, [], 6.795E+00 , True],
    'id103350': ['GluGluToHToWWTo2Tau2Nu_M350', 0.000000295375981 ,  False, [], 5.547E+00 , True],
    'id103400': ['GluGluToHToWWTo2Tau2Nu_M400', 0.00000021064533  ,  False, [], 4.429E+00 , True],
    'id103450': ['GluGluToHToWWTo2Tau2Nu_M450', 0.000000130808187 ,  False, [], 4.195E+00 , True],
    'id103500': ['GluGluToHToWWTo2Tau2Nu_M500', 0.000000080531453 ,  False, [], 4.167E+00 , True],
    'id103550': ['GluGluToHToWWTo2Tau2Nu_M550', 0.00000005009116  ,  False, [], 4.228E+00 , True],
    'id103600': ['GluGluToHToWWTo2Tau2Nu_M600', 0.000000031577974 ,  False, [], 4.350E+00 , True],
    'id104120': ['VBF_HToWWTo2L2Nu_M120',       0.000000076065923 ,  False, [], 6.098E-01 , True],
    'id104130': ['VBF_HToWWTo2L2Nu_M130',       0.000000149109518 ,  False, [], 6.357E-01 , True],
    'id104140': ['VBF_HToWWTo2L2Nu_M140',       0.000000226785526 ,  False, [], 7.015E-01 , True],
    'id104150': ['VBF_HToWWTo2L2Nu_M150',       0.000000289563835 ,  False, [], 8.016E-01 , True],
    'id104160': ['VBF_HToWWTo2L2Nu_M160',       0.000000346989969 ,  False, [], 9.395E-01 , True],
    'id104170': ['VBF_HToWWTo2L2Nu_M170',       0.000000341601082 ,  False, [], 9.874E-01 , True],
    'id104180': ['VBF_HToWWTo2L2Nu_M180',       0.000000303386416 ,  False, [], 9.922E-01 , True],
    'id104190': ['VBF_HToWWTo2L2Nu_M190',       0.000000233881502 ,  False, [], 9.792E-01 , True],
    'id104200': ['VBF_HToWWTo2L2Nu_M200',       0.000000225728111 ,  False, [], 9.838E-01 , True],
    'id104250': ['VBF_HToWWTo2L2Nu_M250',       0.000000129831038 ,  False, [], 9.922E-01 , True],
    'id104300': ['VBF_HToWWTo2L2Nu_M300',       0.000000089252544 ,  False, [], 9.937E-01 , True],
    'id104350': ['VBF_HToWWTo2L2Nu_M350',       0.00000006172483  ,  False, [], 1.005E+00 , True],
    'id104400': ['VBF_HToWWTo2L2Nu_M400',       0.000000040649832 ,  False, [], 1.001E+00 , True],
    'id104450': ['VBF_HToWWTo2L2Nu_M450',       0.000000027347591 ,  False, [], 9.960E-01 , True],
    'id104500': ['VBF_HToWWTo2L2Nu_M500',       0.000000020336535 ,  False, [], 9.952E-01 , True],
    'id104550': ['VBF_HToWWTo2L2Nu_M550',       0.000000015644785 ,  False, [], 9.962E-01 , True],
    'id104600': ['VBF_HToWWTo2L2Nu_M600',       0.000000012027343 ,  False, [], 9.980E-01 , True],
    'id105120': ['VBF_HToWWToLNuTauNu_M120',    0.000000076069385 ,  False, [], 6.098E-01 , True],
    'id105130': ['VBF_HToWWToLNuTauNu_M130',    0.00000016674931  ,  False, [], 6.357E-01 , True],
    'id105140': ['VBF_HToWWToLNuTauNu_M140',    0.000000226764879 ,  False, [], 7.015E-01 , True],
    'id105150': ['VBF_HToWWToLNuTauNu_M150',    0.000000289555928 ,  False, [], 8.016E-01 , True],
    'id105160': ['VBF_HToWWToLNuTauNu_M160',    0.000000363106899 ,  False, [], 9.395E-01 , True],
    'id105170': ['VBF_HToWWToLNuTauNu_M170',    0.000000338459909 ,  False, [], 9.874E-01 , True],
    'id105180': ['VBF_HToWWToLNuTauNu_M180',    0.000000300564037 ,  False, [], 9.922E-01 , True],
    'id105190': ['VBF_HToWWToLNuTauNu_M190',    0.000000233149132 ,  False, [], 9.792E-01 , True],
    'id105200': ['VBF_HToWWToLNuTauNu_M200',    0.000000212613162 ,  False, [], 9.838E-01 , True],
    'id105250': ['VBF_HToWWToLNuTauNu_M250',    0.000000129816857 ,  False, [], 9.922E-01 , True],
    'id105300': ['VBF_HToWWToLNuTauNu_M300',    0.000000089265545 ,  False, [], 9.937E-01 , True],
    'id105350': ['VBF_HToWWToLNuTauNu_M350',    0.000000061517872 ,  False, [], 1.005E+00 , True],
    'id105400': ['VBF_HToWWToLNuTauNu_M400',    0.000000040517962 ,  False, [], 1.001E+00 , True],
    'id105450': ['VBF_HToWWToLNuTauNu_M450',    0.00000002725475  ,  False, [], 9.960E-01 , True],
    'id105500': ['VBF_HToWWToLNuTauNu_M500',    0.000000020275962 ,  False, [], 9.952E-01 , True],
    'id105550': ['VBF_HToWWToLNuTauNu_M550',    0.000000015495872 ,  False, [], 9.962E-01 , True],
    'id105600': ['VBF_HToWWToLNuTauNu_M600',    0.000000013115185 ,  False, [], 9.980E-01 , True],
    'id106120': ['VBF_HToWWTo2Tau2Nu_M120',     0.000000031702597 ,  False, [], 6.098E-01 , True],
    'id106130': ['VBF_HToWWTo2Tau2Nu_M130',     0.000000062129724 ,  False, [], 6.357E-01 , True],
    'id106140': ['VBF_HToWWTo2Tau2Nu_M140',     0.000000104990794 ,  False, [], 7.015E-01 , True],
    'id106150': ['VBF_HToWWTo2Tau2Nu_M150',     0.00000012062744  ,  False, [], 8.016E-01 , True],
    'id106160': ['VBF_HToWWTo2Tau2Nu_M160',     0.000000144582225 ,  False, [], 9.395E-01 , True],
    'id106170': ['VBF_HToWWTo2Tau2Nu_M170',     0.000000155540484 ,  False, [], 9.874E-01 , True],
    'id106180': ['VBF_HToWWTo2Tau2Nu_M180',     0.000000125223237 ,  False, [], 9.922E-01 , True],
    'id106190': ['VBF_HToWWTo2Tau2Nu_M190',     0.000000097153729 ,  False, [], 9.792E-01 , True],
    'id106200': ['VBF_HToWWTo2Tau2Nu_M200',     0.000000084408683 ,  False, [], 9.838E-01 , True],
    'id106250': ['VBF_HToWWTo2Tau2Nu_M250',     0.000000054087892 ,  False, [], 9.922E-01 , True],
    'id106300': ['VBF_HToWWTo2Tau2Nu_M300',     0.000000037192954 ,  False, [], 9.937E-01 , True],
    'id106350': ['VBF_HToWWTo2Tau2Nu_M350',     0.000000025637111 ,  False, [], 1.005E+00 , True],
    'id106400': ['VBF_HToWWTo2Tau2Nu_M400',     0.000000016135343 ,  False, [], 1.001E+00 , True],
    'id106450': ['VBF_HToWWTo2Tau2Nu_M450',     0.00000001136354  ,  False, [], 9.960E-01 , True],
    'id106500': ['VBF_HToWWTo2Tau2Nu_M500',     0.000000008453466 ,  False, [], 9.952E-01 , True],
    'id106550': ['VBF_HToWWTo2Tau2Nu_M550',     0.000000006453974 ,  False, [], 9.962E-01 , True],
    'id106600': ['VBF_HToWWTo2Tau2Nu_M600',     0.000000005013626 ,  False, [], 9.980E-01 , True],
}


#loop over the looper
for key in higgsSamples:
    higgsSamples[key][3] += [ '%s'%x for x in glob(prefix+key+'/*.root') ]







































































































