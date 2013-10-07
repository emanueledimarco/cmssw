from WWAnalysis.AnalysisStep.mcAnalysis import *
from optparse import OptionParser
import pickle

def main():
    parser = OptionParser(usage="%prog [options] mcanalysis.root cuts.txt")
    addMCAnalysisOptions(parser)
    parser.add_option("--pickle",      dest="pickle", action="store_true", help="Import for the pickle file instead of re-running");
    parser.add_option("--pickle-file", dest="pickleFile", type="string",  help="name of the file being pickled to/from (dyData.{lumi:d}pb.pkl)");
    (options, args) = parser.parse_args()

    if not options.pickleFile:
        options.pickleFile = "dyData.{lumi:d}pb.pkl".format(lumi=int(options.lumi*1000))

    # if unpickling, we don't give a shit about anything else
    yields = {}
    for c in [EE,MM,OF,SF]: yields[c] = {}
    if options.pickle:
        pkl_file = open(options.pickleFile, 'rb')
        yields = pickle.load(pkl_file)
        pkl_file.close()
    else: 
        # set up the mca
        mca = MCAnalysis(args[0], options)

        # add the first cf
        cf  = CutsFile(args[1],options)

        # and any subsequent cf's
        for cutFile in args[2:]:
            temp = CutsFile(cutFile,options)
            for cut in temp.cuts():
                cf.add(cut[0],cut[1])


        cf.replace('Z veto','Z veto','abs(mll-91.1876)>15').replace('40/20','mpmet','mpmet>37+nvtx/2').replace('Delta.phi.ll','dphillj','dphiveto')
        cf.remove('m_T')

        selections = { }
        selections['yields'   ] =  CutsFile(cf)                                                                         # out of peak
        selections['peak40'   ] =  CutsFile(cf).replace('Z veto','Z peak','abs(mll-91.1876)<=7.5').remove('m_.*ell')    # nominal
        selections['peak30'   ] =  CutsFile(selections['peak40']).replace('mpmet','mpmet','mpmet>30&&mpmet<=37+nvtx/2') # in case stats too low
        selections['peak25'   ] =  CutsFile(selections['peak40']).replace('mpmet','mpmet','mpmet>25&&mpmet<=30')        # for the systematics
        selections['peak20'   ] =  CutsFile(selections['peak40']).replace('mpmet','mpmet','mpmet>20&&mpmet<=25')        # for the systematics

        reports = {}; 
        for S,C in selections.items():
            print "=========== %s ==========" % S
            print C
            reports[S] = mca.getYields(C)
            #reports[S] = mca.getYields(C,subprocesses=True)
            for c in [EE,MM,OF,SF]: 
                yields[c][S] = {
                    'data' : (float(reports[S]['data'][-1][1][c][1]),reports[S]['data'][-1][1][c][2]),
                    'dy'   : (float(sum( reports[S][p][-1][1][c][1]    for p in reports[S].keys() if p[0:2] == 'DY')),
                              sqrt( sum( reports[S][p][-1][1][c][2]**2 for p in reports[S].keys() if p[0:2] == 'DY') ) ),
                    'vv'   : (float(sum( reports[S][p][-1][1][c][1]    for p in reports[S].keys() if p == 'VV' and p[0] != '.')),
                              sqrt( sum( reports[S][p][-1][1][c][2]**2 for p in reports[S].keys() if p == 'VV' and p[0] != '.') ) ),
                    'other': (float(sum( reports[S][p][-1][1][c][1]    for p in reports[S].keys() if p[0:2] != 'DY' and p not in ['data','VV','ggH','vbfH'] and p[0] != '.')),
                              sqrt( sum( reports[S][p][-1][1][c][2]**2 for p in reports[S].keys() if p[0:2] != 'DY' and p not in ['data','VV','ggH','vbfH'] and p[0] != '.') ) ),
                    'all'  : (float(sum( reports[S][p][-1][1][c][1]    for p in reports[S].keys() if p not in ['data','ggH','vbfH'] and p[0] != '.')),
                              sqrt( sum( reports[S][p][-1][1][c][2]**2 for p in reports[S].keys() if p not in ['data','ggH','vbfH'] and p[0] != '.') ) ),
                }

        pkl_file = open(options.pickleFile, 'wb')
        pickle.dump(yields, pkl_file, -1)
        pkl_file.close()

    for n,c in zip(["ee","mm","of","sf"],[EE,MM,OF,SF]):
        divider('-',120)
        print "{0:^30s}".format("All Selections (%s)"%n) 
        divider('-',120)
        for S, y in sorted(yields[c].items()):
            print "Selection %-15s:  n(dy, sim.) = %8.2f +/- %5.2f  n(vv, sim.) = %8.2f +/- %5.2f   n(data) = %8.2f" % ( 
                S, y['dy'][0], y['dy'][1], y['vv'][0], y['vv'][1], y['data'][0] )
        divider('-',120)

    # Let's do some calculating!
    divider()


    nllIn  = yields[SF]['peak40']['data']
    nemIn  = yields[OF]['peak40']['data']
    nOutDY = yields[SF]['yields']['dy']  
    nInDY  = yields[SF]['peak40']['dy']  
    Rll    = div(nOutDY,nInDY)
    nVV    = yields[SF]['peak40']['vv'] 
    kmm2   = div( yields[MM]['peak40']['dy'] , yields[EE]['peak40']['dy'] )
    kmm    = ( sqrt(kmm2[0]), 0.5*sqrt(kmm2[0])*kmm2[1] )
    kll    = (0.5*(kmm[0] + 1./kmm[0]) , kmm[1] * (1.-1./pow(kmm[0],2)) )
    mcExp  = yields[SF][ 'yields']['dy'] 

    fmt = "%8.2f   +/- %8.2f"
    print "nOutDY", fmt % nOutDY
    print "nInDY ", fmt % nInDY 
    print "Rll   ", fmt % Rll   
    print "nllIn ", fmt % nllIn 
    print "nemIn ", fmt % nemIn 
    print "nVV   ", fmt % nVV   
    print "kmm   ", fmt % kmm   
    print "kll   ", fmt % kll   
    print "mcExp ", fmt % mcExp   

    predicted = nDYout(nllIn, tuple([2*x for x in nemIn]), Rll, kll, nVV)
    print "nDYout", fmt % predicted
    print "ratio ", fmt % div( predicted, mcExp )
    divider()

def quadSum(x1,x2,x3=0,x4=0,x5=0,x6=0,x7=0,x8=0):
    return sqrt(x1*x1 + x2*x2 + x3*x3 + x4*x4 + x5*x5 + x6*x6 + x7*x7 + x8*x8)

def nDYout(nDYin, nemu, R, K, nZV):
    val = R[0] * (nDYin[0] - 0.5 * nemu[0] * K[0] - nZV[0])
    err = quadSum(R[1]*(nDYin[0] - 0.5 * nemu[0] * K[0] -nZV[0]),
            R[0] * sqrt(nDYin[0]),
            0.5 * R[0] * K[0] * sqrt(nemu[0]),
            0.5 * R[0] * nemu[0] * K[1],
            R[0] * nZV[1])
    return (val,err)


def divider(str='-',len=45):
    print str*len

def div(a,b):
    return (float(a[0])/b[0],divErr(a,b)) if b!=0 else (0,0)

def divErr(a,b):
    return a[0]/b[0]*sqrt( pow(float(a[1])/a[0],2) + pow(float(b[1])/b[0],2) )

def multErr(a,b):
    return a[0]*b[0]*sqrt( pow(float(a[1])/a[0],2) + pow(float(b[1])/b[0],2) )

def effNDE(num,den,errNum,errDen): 
    """
sigma^2(eff) = [ (sumw2 pass) * (1 - 2*eff) + (sumw2 all) * eff*eff ] / (sumw)^2 = 
             = [ (sumw2 pass) * (1-eff)^2   + (sumw2 fail) * eff^2 ] / (sumw)^2 
And of course sumw2(X) = err(X)^2
              sumw     = den
    """
    if den == 0: return (0,0)
    eff  = num/float(den)
    errPass2 = errNum*errNum;
    errFail2 = errDen*errDen - errPass2
    err = sqrt( errPass2*(1-eff)**2 + errFail2*(eff**2) ) / den
    return ( eff, err )

if __name__ == "__main__":
    main()



### do1 = False if options.zorro else True
### datasets = [
### #     "SingleMuon2011AMay10",     
### #     "SingleMuon2011Av4",     
### #     "DoubleMuon2011AMay10", 
### #     "DoubleMuon2011Av4",
### #     "DoubleElectron2011AMay10", 
### #     "DoubleElectron2011Av4", 
### #     "MuEG2011AMay10",       
### #     "MuEG2011Av4" 
###     "083_SingleElectron2011Av4",
###     "083b_SingleElectron2011Av4",
###     "084_SingleMuon2011Av4",
###     "084b_SingleMuon2011Av4",
###     "085_DoubleElectron2011Av4",
###     "085b_DoubleElectron2011Av4",
###     "086_DoubleMuon2011Av4",
###     "086b_DoubleMuon2011Av4",
###     "087_MuEG2011Av4",
###     "087b_MuEG2011Av4",
###     "089_SingleElectron2011AMay10",
###     "090_SingleMuon2011AMay10",
###     "091_DoubleMuon2011AMay10",
###     "092_DoubleElectron2011AMay10",
###     "093_MuEG2011AMay10",
### ]
### mcsets   = [ 
###     "003.42X_DYtoElEl",
###     "004.42X_DYtoMuMu",
### #     "DYtoMuMu", 
### #     "DYtoElEl", 
###     "DY10toMuMuZ2", 
###     "DY10toElElZ2" 
### ];
### mcmad    = [ "DY50toLLMadD6T" ];
### mcsetsvv = [ "WZtoAny", "ZZtoAny" ]
### runrange = "run > %s && run < %s" % (options.minrun,options.maxrun)
### print runrange
### if options.doVV and not options.subVV: 
###     mcsets += mcsetsvv
###     print "Including WZ->any and ZZ->any in the MC samples."
### elif options.subVV: 
###     print "Including WZ->any and ZZ->any contributions will be subtracted from peak yield"
### if options.mad:
###     print "Including MadGraph samples (for Z peak yields only)"
### ttysData = [ TreeToYield(options.path+"/tree_%s.root" % D, options) for D in datasets ]
### ttysMC   = [ TreeToYield(options.path+"/tree_%s.root" % D, options) for D in mcsets   ] 
### ttysMCVV = [ TreeToYield(options.path+"/tree_%s.root" % D, options) for D in mcsetsvv ] 
### ttysMCM  = [ TreeToYield(options.path+"/tree_%s.root" % D, options) for D in mcmad    if options.mad ] 
### cf0j = CutsFile(args[0],options)
### cf0j.addAfter('trigger','peaking','peaking')
### cf1j       = CutsFile(cf0j).replace('jet veto', 'one jet', 'njet == 1').replace('soft b', 'b-tag veto', 'bveto_ip && nbjet == 0')
### cf0jNoMET  = CutsFile(cf0j).remove('met')
### cf1jNoMET  = CutsFile(cf1j).remove('met')
### cfZ0j      = CutsFile(cf0j).replace('Z veto','Z window','abs(mll-91.1876)<15').remove('m_T|gamma|Delta.phi|m_.*ell|p_T.*(min|max)').replace('40/20','pMET > 40', 'pmmet > 40')
### cfZ1j      = CutsFile(cf1j).replace('Z veto','Z window','abs(mll-91.1876)<15').remove('m_T|gamma|Delta.phi|m_.*ell|p_T.*(min|max)').replace('40/20','pMET > 40', 'pmmet > 40')
### cfZ0jNoMET = CutsFile(cfZ0j).remove('met').remove('pMET')
### cfZ1jNoMET = CutsFile(cfZ1j).remove('met').remove('pMET')
### cf0jData       = CutsFile(cf0j       ).replace('trigger','trigger','trigguil && '+runrange).remove("peaking")
### cf1jData       = CutsFile(cf1j       ).replace('trigger','trigger','trigguil && '+runrange).remove("peaking")
### cf0jNoMETData  = CutsFile(cf0jNoMET  ).replace('trigger','trigger','trigguil && '+runrange).remove("peaking")
### cf1jNoMETData  = CutsFile(cf1jNoMET  ).replace('trigger','trigger','trigguil && '+runrange).remove("peaking")
### cfZ0jData      = CutsFile(cfZ0j      ).replace('trigger','trigger','trigguil && '+runrange).remove("peaking")
### cfZ1jData      = CutsFile(cfZ1j      ).replace('trigger','trigger','trigguil && '+runrange).remove("peaking")
### cfZ0jNoMETData = CutsFile(cfZ0jNoMET ).replace('trigger','trigger','trigguil && '+runrange).remove("peaking")
### cfZ1jNoMETData = CutsFile(cfZ1jNoMET ).replace('trigger','trigger','trigguil && '+runrange).remove("peaking")
### 
### # print 'cf0j          \n', cf0j          , '\n'
### # print 'cf1j          \n', cf1j          , '\n'
### # print 'cf0jNoMET     \n', cf0jNoMET     , '\n'
### # print 'cf1jNoMET     \n', cf1jNoMET     , '\n'
### # print 'cfZ0j         \n', cfZ0j         , '\n'
### # print 'cfZ1j         \n', cfZ1j         , '\n'
### # print 'cfZ0jNoMET    \n', cfZ0jNoMET    , '\n'
### # print 'cfZ1jNoMET    \n', cfZ1jNoMET    , '\n'
### # print 'cf0jData      \n', cf0jData      , '\n'
### # print 'cf1jData      \n', cf1jData      , '\n'
### # print 'cf0jNoMETData \n', cf0jNoMETData , '\n'
### # print 'cf1jNoMETData \n', cf1jNoMETData , '\n'
### # print 'cfZ0jData     \n', cfZ0jData     , '\n'
### # print 'cfZ1jData     \n', cfZ1jData     , '\n'
### # print 'cfZ0jNoMETData\n', cfZ0jNoMETData, '\n'
### # print 'cfZ1jNoMETData\n', cfZ1jNoMETData, '\n'
### # sys.exit()
### 
### ## get yields in peak, with MET cut
### reportDataZ0j = mergeReports([tty.getYields(cfZ0jData) for tty in ttysData])
### reportMCZ0j   = mergeReports([tty.getYields(cfZ0j) for tty in ttysMC])
### if do1: 
###     reportDataZ1j = mergeReports([tty.getYields(cfZ1jData) for tty in ttysData])
###     reportMCZ1j   = mergeReports([tty.getYields(cfZ1j) for tty in ttysMC])
### if options.subVV:
###     reportMCVVZ0j = mergeReports([tty.getYields(cfZ0j) for tty in ttysMCVV])
###     if do1: 
###         reportMCVVZ1j = mergeReports([tty.getYields(cfZ1j) for tty in ttysMCVV])
###     
### if options.mad:
###     reportMCZ0j = mergeReports([reportMCZ0j] + [tty.getYields(cfZ0j) for tty in ttysMCM])
###     for i in range(4): 
###         reportMCZ0j[-1][1][i] = [ reportMCZ0j[-1][1][i][0], 0.5*reportMCZ0j[-1][1][i][1], 0.5*reportMCZ0j[-1][1][i][2] ]
###     if do1:
###         reportMCZ1j = mergeReports([reportMCZ1j] + [tty.getYields(cfZ1j) for tty in ttysMCM])
###         for i in range(4): 
###             reportMCZ1j[-1][1][i] = [ reportMCZ1j[-1][1][i][0], 0.5*reportMCZ1j[-1][1][i][1], 0.5*reportMCZ1j[-1][1][i][2] ]
### print "\nZ peak in Data, WW level no btag, with MET cut, 0 jets"; ttysData[0].prettyPrint(reportDataZ0j)
### print "\nZ peak in Sim., WW level no btag, with MET cut, 0 jets";   ttysMC[0].prettyPrint(reportMCZ0j)
### if options.subVV: 
###     print "\nZ peak in Sim. VV, WW level no btag, with MET cut, 0 jets";
###     ttysMCVV[0].prettyPrint(reportMCVVZ0j)
### if do1:
###     print "\nZ peak in Data, WW level no btag, with MET cut, 1 jets"; ttysData[0].prettyPrint(reportDataZ1j)
###     print "\nZ peak in Sim., WW level no btag, with MET cut, 1 jets";   ttysMC[0].prettyPrint(reportMCZ1j)
###     if options.subVV: 
###         print "\nZ peak in Sim. VV, WW level no btag, with MET cut, 1 jets";
###         ttysMCVV[0].prettyPrint(reportMCVVZ1j)
### Z0jmmData  = reportDataZ0j[-1][1][MM][1]; Z0jeeData  = reportDataZ0j[-1][1][EE][1]; Z0jxxData  = (reportDataZ0j[-1][1][EM][1]+reportDataZ0j[-1][1][ME][1])
### Z0jmmMC    = reportMCZ0j[-1][1][MM][1];   Z0jeeMC    = reportMCZ0j[-1][1][EE][1];   Z0jxxMC    = (reportMCZ0j[-1][1][EM][1]+reportMCZ0j[-1][1][ME][1])
### Z0jmmMCErr = reportMCZ0j[-1][1][MM][2];   Z0jeeMCErr = reportMCZ0j[-1][1][EE][2];   Z0jxxMCErr = hypot(reportMCZ0j[-1][1][EM][2], reportMCZ0j[-1][1][ME][2])
### if options.subVV:
###     Z0jmmMCVV    = reportMCVVZ0j[-1][1][MM][1];   Z0jeeMCVV    = reportMCVVZ0j[-1][1][EE][1];   Z0jxxMCVV    = (reportMCVVZ0j[-1][1][EM][1]+reportMCVVZ0j[-1][1][ME][1])
###     Z0jmmMCVVErr = reportMCVVZ0j[-1][1][MM][2];   Z0jeeMCVVErr = reportMCVVZ0j[-1][1][EE][2];   Z0jxxMCVVErr = hypot(reportMCVVZ0j[-1][1][EM][2], reportMCVVZ0j[-1][1][ME][2])
### if do1:
###     Z1jmmData  = reportDataZ1j[-1][1][MM][1]; Z1jeeData  = reportDataZ1j[-1][1][EE][1]; Z1jxxData  = (reportDataZ1j[-1][1][EM][1]+reportDataZ1j[-1][1][ME][1])
###     Z1jmmMC    = reportMCZ1j[-1][1][MM][1];   Z1jeeMC    = reportMCZ1j[-1][1][EE][1];   Z1jxxMC    = (reportMCZ1j[-1][1][EM][1]+reportMCZ1j[-1][1][ME][1])
###     Z1jmmMCErr = reportMCZ1j[-1][1][MM][2];   Z1jeeMCErr = reportMCZ1j[-1][1][EE][2];   Z1jxxMCErr = hypot(reportMCZ1j[-1][1][EM][2], reportMCZ1j[-1][1][ME][2])
###     if options.subVV:
###         Z1jmmMCVV    = reportMCVVZ1j[-1][1][MM][1];   Z1jeeMCVV    = reportMCVVZ1j[-1][1][EE][1];   Z1jxxMCVV    = (reportMCVVZ1j[-1][1][EM][1]+reportMCVVZ1j[-1][1][ME][1])
###         Z1jmmMCVVErr = reportMCVVZ1j[-1][1][MM][2];   Z1jeeMCVVErr = reportMCVVZ1j[-1][1][EE][2];   Z1jxxMCVVErr = hypot(reportMCVVZ1j[-1][1][EM][2], reportMCVVZ1j[-1][1][ME][2])
### 
### ## Now we compute events in the peak without MET cut, to get efficiencies
### reportDataZNoMET = mergeReports([tty.getYields(cfZ0jNoMETData) for tty in ttysData])
### reportMCZNoMET   = mergeReports([tty.getYields(cfZ0jNoMET) for tty in ttysMC])
### if options.subVV: 
###     reportMCVVZNoMET   = mergeReports([tty.getYields(cfZ0jNoMET) for tty in ttysMCVV])
### print "\nZ peak in Data, WW level no btag, no MET, 0 jets"; ttysData[0].prettyPrint(reportDataZNoMET)
### print "\nZ peak in Sim., WW level no btag, no MET, 0 jets";   ttysMC[0].prettyPrint(reportMCZNoMET)
### if options.subVV: 
###     print "\nZ peak in Sim. VV, WW level no btag, no MET, 0 jets";   
###     ttysMCVV[0].prettyPrint(reportMCVVZNoMET)
### keeData = sqrt(reportDataZNoMET[-1][1][EE][1]/float(reportDataZNoMET[-1][1][MM][1]));
### kmmData = 1.0/keeData;
### keeMC = sqrt(reportMCZNoMET[-1][1][EE][1]/float(reportMCZNoMET[-1][1][MM][1]));
### kmmMC = 1.0/keeMC;
### 
### print ""
### ## Now we compute events in the signal region. we need the weight regularize the yields
### avgweightMM  = ttysMC[0].getAverageWeight(cfZ0jNoMET.allCuts());
### avgweightEE  = ttysMC[1].getAverageWeight(cfZ0jNoMET.allCuts());
### avgweightSMM = avgweightMM #ttysMC[0].getAverageWeight(cf0j.allCuts());
### avgweightSEE = avgweightEE #ttysMC[1].getAverageWeight(cf0j.allCuts());
### print "Average event weight in Z peak w/o MET: ",avgweightMM, "(mm),\t",avgweightEE,"(ee)"
### print "Average event weight in signal region : ",avgweightSMM,"(mm),\t",avgweightSEE,"(ee)"
### 
### ## Do flavour subtraction to get DY-only yield
### Z0jeeSubData = Z0jeeData - 0.5*keeData*Z0jxxData; Z0jeeSubDataErr =  sqrt(Z0jeeData+pow(0.5*keeData,2)*Z0jxxData)
### Z0jeeSubMC   = Z0jeeMC   - 0.5*keeMC  *Z0jxxMC;   Z0jeeSubMCErr   = hypot(Z0jeeMCErr,0.5*keeMC*Z0jxxMCErr)
### Z0jmmSubData = Z0jmmData - 0.5*kmmData*Z0jxxData; Z0jmmSubDataErr =  sqrt(Z0jmmData+pow(0.5*kmmData,2)*Z0jxxData)
### Z0jmmSubMC   = Z0jmmMC   - 0.5*kmmMC  *Z0jxxMC;   Z0jmmSubMCErr   = hypot(Z0jmmMCErr,0.5*kmmMC*Z0jxxMCErr)
### if options.subVV:
###     Z0jeeSubMCVV = Z0jeeMCVV   - 0.5*keeMC*Z0jxxMCVV;   Z0jeeSubMCVVErr   = hypot(Z0jeeMCVVErr,0.5*keeMC*Z0jxxMCVVErr)
###     Z0jmmSubMCVV = Z0jmmMCVV   - 0.5*kmmMC*Z0jxxMCVV;   Z0jmmSubMCVVErr   = hypot(Z0jmmMCVVErr,0.5*kmmMC*Z0jxxMCVVErr)
###     Z0jmmSubData -= Z0jmmSubMCVV; Z0jmmDataErr = sqrt(Z0jmmSubDataErr**2 + Z0jmmSubMCVVErr**2 + (0.1*Z0jmmSubMCVV)**2)
###     Z0jeeSubData -= Z0jeeSubMCVV; Z0jeeDataErr = sqrt(Z0jeeSubDataErr**2 + Z0jeeSubMCVVErr**2 + (0.1*Z0jeeSubMCVV)**2)
### if do1:
###     Z1jeeSubData = Z1jeeData - 0.5*keeData*Z1jxxData; Z1jeeSubDataErr =  sqrt(Z1jeeData+pow(0.5*keeData,2)*Z1jxxData)
###     Z1jeeSubMC   = Z1jeeMC   - 0.5*keeMC  *Z1jxxMC;   Z1jeeSubMCErr   = hypot(Z1jeeMCErr,0.5*keeMC*Z1jxxMCErr)
###     Z1jmmSubData = Z1jmmData - 0.5*kmmData*Z1jxxData; Z1jmmSubDataErr =  sqrt(Z1jmmData+pow(0.5*kmmData,2)*Z1jxxData)
###     Z1jmmSubMC   = Z1jmmMC   - 0.5*kmmMC  *Z1jxxMC;   Z1jmmSubMCErr   = hypot(Z1jmmMCErr,0.5*kmmMC*Z1jxxMCErr)
###     if options.subVV:
###         Z1jeeSubMCVV   = Z1jeeMCVV   - 0.5*keeMC*Z1jxxMCVV;   Z1jeeSubMCVVErr   = hypot(Z1jeeMCVVErr,0.5*keeMC*Z1jxxMCVVErr)
###         Z1jmmSubMCVV   = Z1jmmMCVV   - 0.5*kmmMC*Z1jxxMCVV;   Z1jmmSubMCVVErr   = hypot(Z1jmmMCVVErr,0.5*kmmMC*Z1jxxMCVVErr)
###         Z1jmmSubData -= Z1jmmSubMCVV; Z1jmmDataErr = sqrt(Z1jmmSubDataErr**2 + Z1jmmSubMCVVErr**2 + (0.1*Z1jmmSubMCVV)**2)
###         Z1jeeSubData -= Z1jeeSubMCVV; Z1jeeDataErr = sqrt(Z1jeeSubDataErr**2 + Z1jeeSubMCVVErr**2 + (0.1*Z1jeeSubMCVV)**2)
### 
### report0jMC = mergeReports([tty.getYields(cf0j) for tty in ttysMC])
### print "\nFull selection in Sim., 0 jets";   ttysMC[0].prettyPrint(report0jMC)
### if do1:
###     report1jMC = mergeReports([tty.getYields(cf1j) for tty in ttysMC])
###     print "\nFull selection in Sim., 1 jets";   ttysMC[0].prettyPrint(report1jMC)
### 
### ## Then we try to do the same but as a function of the MET cut
### plot0j  = mergePlots([tty.getPlots("pmmet_0j","pmmet", "50,0.,50.", cf0jNoMET.allCuts()) for tty in ttysMC])
### plotZ0j = mergePlots([tty.getPlots("pmmetZ_0j","pmmet","50,0.,50.",cfZ0jNoMET.allCuts()) for tty in ttysMC])
### totZmm0j = plotZ0j[MM][1].Integral(1,51); totSmm0j = plot0j[MM][1].Integral(1,51);
### totZee0j = plotZ0j[EE][1].Integral(1,51); totSee0j = plot0j[EE][1].Integral(1,51);
### totZme0j = plotZ0j[ME][1].Integral(1,51); totSme0j = plot0j[ME][1].Integral(1,51);
### totZem0j = plotZ0j[EM][1].Integral(1,51); totSem0j = plot0j[ME][1].Integral(1,51);
### plotData0j  = mergePlots([tty.getPlots("pmmetD_0j","pmmet", "50,0.,50.", cf0jNoMETData.allCuts()) for tty in ttysData])
### plotDataZ0j = mergePlots([tty.getPlots("pmmetDZ_0j","pmmet","50,0.,50.",cfZ0jNoMETData.allCuts()) for tty in ttysData])
### totDataZee0j = plotDataZ0j[EE][1].Integral(1,51); totDataSee0j = plotData0j[EE][1].Integral(1,51);
### totDataZem0j = plotDataZ0j[EM][1].Integral(1,51); totDataSem0j = plotData0j[EM][1].Integral(1,51);
### totDataZme0j = plotDataZ0j[ME][1].Integral(1,51); totDataSme0j = plotData0j[ME][1].Integral(1,51); 
### totDataZmm0j = plotDataZ0j[MM][1].Integral(1,51); totDataSmm0j = plotData0j[MM][1].Integral(1,51); 
### if do1:
###     plot1j  = mergePlots([tty.getPlots("pmmet_1j","pmmet", "50,0.,50.", cf1jNoMET.allCuts()) for tty in ttysMC])
###     plotZ1j = mergePlots([tty.getPlots("pmmetZ_1j","pmmet","50,0.,50.",cfZ1jNoMET.allCuts()) for tty in ttysMC])
###     totZmm1j = plotZ1j[MM][1].Integral(1,51); totSmm1j = plot1j[MM][1].Integral(1,51);
###     totZme1j = plotZ1j[ME][1].Integral(1,51); totSme1j = plot1j[ME][1].Integral(1,51);
###     totZem1j = plotZ1j[EM][1].Integral(1,51); totSem1j = plot1j[EM][1].Integral(1,51);
###     totZee1j = plotZ1j[EE][1].Integral(1,51); totSee1j = plot1j[EE][1].Integral(1,51);
###     plotData1j  = mergePlots([tty.getPlots("pmmetD_1j","pmmet", "50,0.,50.", cf1jNoMETData.allCuts()) for tty in ttysData])
###     plotDataZ1j = mergePlots([tty.getPlots("pmmetDZ_1j","pmmet","50,0.,50.",cfZ1jNoMETData.allCuts()) for tty in ttysData])
###     totDataZee1j = plotDataZ1j[EE][1].Integral(1,51); totDataSee1j = plotData1j[EE][1].Integral(1,51);
###     totDataZem1j = plotDataZ1j[EM][1].Integral(1,51); totDataSem1j = plotData1j[EM][1].Integral(1,51);
###     totDataZme1j = plotDataZ1j[ME][1].Integral(1,51); totDataSme1j = plotData1j[ME][1].Integral(1,51);
###     totDataZmm1j = plotDataZ1j[MM][1].Integral(1,51); totDataSmm1j = plotData1j[MM][1].Integral(1,51);
### 
### ## Compute extrapolation factors for nominal pMET cut
### alpha0jee = (report0jMC[-1][1][EE][1]+avgweightSEE)/(reportMCZ0j[-1][1][EE][1]+avgweightEE);
### alpha0jeeErr = min(1.0, sqrt(1.0/(report0jMC[-1][1][EE][1]/avgweightSEE+1) + 1.0/(reportMCZ0j[-1][1][EE][1]/avgweightEE+1)));
### alpha0jmm = (report0jMC[-1][1][MM][1]+avgweightSMM)/(reportMCZ0j[-1][1][MM][1]+avgweightMM);
### alpha0jmmErr = min(1.0, sqrt(1.0/(report0jMC[-1][1][MM][1]/avgweightSMM+1) + 1.0/(reportMCZ0j[-1][1][MM][1]/avgweightMM+1)));
### if do1:
###     alpha1jee = (report1jMC[-1][1][EE][1]+avgweightSEE)/(reportMCZ1j[-1][1][EE][1]+avgweightEE);
###     alpha1jeeErr = min(1.0, sqrt(1.0/(report1jMC[-1][1][EE][1]/avgweightSEE+1) + 1.0/(reportMCZ1j[-1][1][EE][1]/avgweightEE+1)));
###     alpha1jmm = (report1jMC[-1][1][MM][1]+avgweightSMM)/(reportMCZ1j[-1][1][MM][1]+avgweightMM);
###     alpha1jmmErr = min(1.0, sqrt(1.0/(report1jMC[-1][1][MM][1]/avgweightSMM+1) + 1.0/(reportMCZ1j[-1][1][MM][1]/avgweightMM+1)));
### 
### ## Let's do a scan of alpha vs pMET cut
### tfout = ROOT.TFile(options.out, "RECREATE")
### print "\nEstimate of alpha vs pMET, regularized, 0 jets. Data is flavour subtracted."
### print "\tpMET  ------------ MC -------------  |  ------------- DATA ----------  |  ------------ MC -------------  |  ------------- DATA ----------"
### print "\tpMET    Zmm0j    Smm0j   alpha  err  |    Zmm0j    Smm0j   alpha  err  |    Zee0j    See0j   alpha  err  |    Zee0j    See0j   alpha  err"
### gZmm0jMC   = ROOT.TGraphAsymmErrors(50); gZee0jMC   = ROOT.TGraphAsymmErrors(50); gSmm0jMC   = ROOT.TGraphAsymmErrors(50); gSee0jMC   = ROOT.TGraphAsymmErrors(50);
### gZmm0jData = ROOT.TGraphAsymmErrors(50); gZee0jData = ROOT.TGraphAsymmErrors(50); gSmm0jData = ROOT.TGraphAsymmErrors(50); gSee0jData = ROOT.TGraphAsymmErrors(50);
### gAmm0jMC   = ROOT.TGraphAsymmErrors(50); gAee0jMC   = ROOT.TGraphAsymmErrors(50); gAmm0jData = ROOT.TGraphAsymmErrors(50); gAee0jData = ROOT.TGraphAsymmErrors(50)
### for b in range(50):
###     totDataSmm0jSub = max(0, totDataSmm0j - 0.5*kmmData * (totDataSem0j + totDataSme0j))
###     totDataSee0jSub = max(0, totDataSee0j - 0.5*keeData * (totDataSem0j + totDataSme0j))
###     totDataSmm0jErr = sqrt(totDataSmm0j + pow(0.5*kmmData,2) * (totDataSem0j + totDataSme0j))/totDataSmm0jSub if totDataSmm0jSub > 0 else 1.0
###     totDataSee0jErr = sqrt(totDataSee0j + pow(0.5*keeData,2) * (totDataSem0j + totDataSme0j))/totDataSee0jSub if totDataSee0jSub > 0 else 1.0
###     totDataZmm0jSub = max(0, totDataZmm0j - 0.5*kmmData * (totDataZem0j + totDataZme0j))
###     totDataZee0jSub = max(0, totDataZee0j - 0.5*keeData * (totDataZem0j + totDataZme0j))
###     totDataZmm0jErr = sqrt(totDataZmm0j + pow(0.5*kmmData,2) * (totDataZem0j + totDataZme0j))/totDataZmm0jSub if totDataZmm0jSub > 0 else 1.0
###     totDataZee0jErr = sqrt(totDataZee0j + pow(0.5*keeData,2) * (totDataZem0j + totDataZme0j))/totDataZee0jSub if totDataZee0jSub > 0 else 1.0
###     Aee0jData = (totDataSee0jSub+1)/(totDataZee0jSub+1);   Aee0jDataErr = sqrt(1.0/(totDataZee0jSub+1) + 1.0/(totDataSee0jErr+1));
###     Amm0jData = (totDataSmm0jSub+1)/(totDataZmm0jSub+1);   Amm0jDataErr = sqrt(1.0/(totDataZmm0jSub+1) + 1.0/(totDataSmm0jErr+1));
###     totSmm0jSub = max(0, totSmm0j - 0.5*kmmMC * (totSem0j + totSme0j))
###     totSee0jSub = max(0, totSee0j - 0.5*keeMC * (totSem0j + totSme0j))
###     totSmm0jErr = sqrt(totSmm0j + pow(0.5*kmmMC,2) * (totSem0j + totSme0j))/totSmm0jSub if totSmm0jSub > 0 else 1.0
###     totSee0jErr = sqrt(totSee0j + pow(0.5*keeMC,2) * (totSem0j + totSme0j))/totSee0jSub if totSee0jSub > 0 else 1.0
###     totZmm0jSub = max(0, totZmm0j - 0.5*kmmMC * (totZem0j + totZme0j))
###     totZee0jSub = max(0, totZee0j - 0.5*keeMC * (totZem0j + totZme0j))
###     totZmm0jErr = sqrt(totZmm0j + pow(0.5*kmmMC,2) * (totZem0j + totZme0j))/totZmm0jSub if totZmm0jSub > 0 else 1.0
###     totZee0jErr = sqrt(totZee0j + pow(0.5*keeMC,2) * (totZem0j + totZme0j))/totZee0jSub if totZee0jSub > 0 else 1.0
###     Aee0j = (totSee0jSub+avgweightSEE)/(totZee0jSub+avgweightEE);   Aee0jErr = sqrt(1.0/(totZee0jSub/avgweightEE+1) + 1.0/(totSee0jErr/avgweightSEE+1));
###     Amm0j = (totSmm0jSub+avgweightSMM)/(totZmm0jSub+avgweightMM);   Amm0jErr = sqrt(1.0/(totZmm0jSub/avgweightEE+1) + 1.0/(totSmm0jErr/avgweightSMM+1));
###     gZmm0jMC.SetPoint(b, 1.0*b, totZmm0jSub); gZmm0jMC.SetPointError(b, 0, 0, totZmm0jErr*totZmm0jSub, totZmm0jErr*totZmm0jSub);
###     gZee0jMC.SetPoint(b, 1.0*b, totZee0jSub); gZee0jMC.SetPointError(b, 0, 0, totZee0jErr*totZee0jSub, totZee0jErr*totZee0jSub);
###     gSmm0jMC.SetPoint(b, 1.0*b, totSmm0jSub); gSmm0jMC.SetPointError(b, 0, 0, totSmm0jErr*totSmm0jSub, totSmm0jErr*totSmm0jSub);
###     gSee0jMC.SetPoint(b, 1.0*b, totSee0jSub); gSee0jMC.SetPointError(b, 0, 0, totSee0jErr*totSee0jSub, totSee0jErr*totSee0jSub);
###     gAmm0jMC.SetPoint(b, 1.0*b,    Amm0j);    gAmm0jMC.SetPointError(b, 0, 0, Amm0j*Amm0jErr, Amm0j*Amm0jErr);
###     gAee0jMC.SetPoint(b, 1.0*b,    Aee0j);    gAee0jMC.SetPointError(b, 0, 0, Aee0j*Aee0jErr, Aee0j*Aee0jErr);
###     gZmm0jData.SetPoint(b, 1.0*b, totDataZmm0jSub); gZmm0jData.SetPointError(b, 0, 0, totDataZmm0jErr, totDataZmm0jErr);
###     gZee0jData.SetPoint(b, 1.0*b, totDataZee0jSub); gZee0jData.SetPointError(b, 0, 0, totDataZee0jErr, totDataZee0jErr);
###     gSmm0jData.SetPoint(b, 1.0*b, totDataSmm0jSub); gSmm0jData.SetPointError(b, 0, 0, totDataSmm0jErr, totDataSmm0jErr);
###     gSee0jData.SetPoint(b, 1.0*b, totDataSee0jSub); gSee0jData.SetPointError(b, 0, 0, totDataSee0jErr, totDataSee0jErr);
###     gAmm0jData.SetPoint(b, 1.0*b,       Amm0jData); gAmm0jData.SetPointError(b, 0, 0, Amm0jData*Amm0jDataErr, Amm0jData*Amm0jDataErr);
###     gAee0jData.SetPoint(b, 1.0*b,       Aee0jData); gAee0jData.SetPointError(b, 0, 0, Aee0jData*Aee0jDataErr, Aee0jData*Aee0jDataErr);
###     print "\t %2d   %9.1f %8.2f  %5.3f %4.2f |  %9.1f %8.2f  %5.3f %4.2f |  %9.1f %8.2f  %5.3f %4.2f |  %9.1f %8.2f  %5.3f %4.2f" % (
###             b,
###             totZmm0jSub, totSmm0jSub, Amm0j, min(1.0, Amm0jErr),
###             totDataZmm0jSub, totDataSmm0jSub, Amm0jData, min(1.0, Amm0jDataErr),
###             totZee0jSub, totSee0jSub, Aee0j, min(1.0, Aee0jErr),
###             totDataZee0jSub, totDataSee0jSub, Aee0jData, min(1.0, Aee0jDataErr),
###           )
###     totZmm0j -= plotZ0j[MM][1].GetBinContent(b+1); totSmm0j -= plot0j[MM][1].GetBinContent(b+1)
###     totZee0j -= plotZ0j[EE][1].GetBinContent(b+1); totSee0j -= plot0j[EE][1].GetBinContent(b+1)
###     totZem0j -= plotZ0j[EM][1].GetBinContent(b+1); totSem0j -= plot0j[EM][1].GetBinContent(b+1)
###     totZme0j -= plotZ0j[ME][1].GetBinContent(b+1); totSme0j -= plot0j[ME][1].GetBinContent(b+1)
###     totDataZmm0j -= plotDataZ0j[MM][1].GetBinContent(b+1); totDataSmm0j -= plotData0j[MM][1].GetBinContent(b+1)
###     totDataZee0j -= plotDataZ0j[EE][1].GetBinContent(b+1); totDataSee0j -= plotData0j[EE][1].GetBinContent(b+1)
###     totDataZem0j -= plotDataZ0j[EM][1].GetBinContent(b+1); totDataSem0j -= plotData0j[EM][1].GetBinContent(b+1)
###     totDataZme0j -= plotDataZ0j[ME][1].GetBinContent(b+1); totDataSme0j -= plotData0j[ME][1].GetBinContent(b+1)
### gZmm0jMC.SetName("gZmm0jMC"); gSmm0jMC.SetName("gSmm0jMC"); gAmm0jMC.SetName("gAmm0jMC"); 
### gZee0jMC.SetName("gZee0jMC"); gSee0jMC.SetName("gSee0jMC"); gAee0jMC.SetName("gAee0jMC");
### gZmm0jData.SetName("gZmm0jData"); gSmm0jData.SetName("gSmm0jData"); gAmm0jData.SetName("gAmm0jData"); 
### gZee0jData.SetName("gZee0jData"); gSee0jData.SetName("gSee0jData"); gAee0jData.SetName("gAee0jData");
### gZmm0jMC.Write(); gZee0jMC.Write(); gSmm0jMC.Write(); gSee0jMC.Write(); gAmm0jMC.Write(); gAee0jMC.Write();
### gZmm0jData.Write(); gZee0jData.Write(); gSmm0jData.Write(); gSee0jData.Write(); gAmm0jData.Write(); gAee0jData.Write();
### 
### if do1:
###     gZmm1jMC   = ROOT.TGraphAsymmErrors(50); gZee1jMC   = ROOT.TGraphAsymmErrors(50); gSmm1jMC   = ROOT.TGraphAsymmErrors(50); gSee1jMC   = ROOT.TGraphAsymmErrors(50);
###     gZmm1jData = ROOT.TGraphAsymmErrors(50); gZee1jData = ROOT.TGraphAsymmErrors(50); gSmm1jData = ROOT.TGraphAsymmErrors(50); gSee1jData = ROOT.TGraphAsymmErrors(50);
###     gAmm1jMC   = ROOT.TGraphAsymmErrors(50); gAee1jMC   = ROOT.TGraphAsymmErrors(50); gAmm1jData = ROOT.TGraphAsymmErrors(50); gAee1jData = ROOT.TGraphAsymmErrors(50)
###     print "\nEstimate of alpha vs pMET, regularized, 1 jet. Data is flavour subtracted."
###     print "\tpMET  ------------ MC -------------  |  ------------- DATA ----------  |  ------------ MC -------------  |  ------------- DATA ----------"
###     print "\tpMET    Zmm1j    Smm1j   alpha  err  |    Zmm1j    Smm1j   alpha  err  |    Zee1j    See1j   alpha  err  |    Zee1j    See1j   alpha  err"
###     for b in range(50):
###         totDataSmm1jSub = max(0, totDataSmm1j - 0.5*kmmData * (totDataSem1j + totDataSme1j))
###         totDataSee1jSub = max(0, totDataSee1j - 0.5*keeData * (totDataSem1j + totDataSme1j))
###         totDataSmm1jErr = sqrt(totDataSmm1j + pow(0.5*kmmData,2) * (totDataSem1j + totDataSme1j))/totDataSmm1jSub if totDataSmm1jSub > 0 else 1.0
###         totDataSee1jErr = sqrt(totDataSee1j + pow(0.5*keeData,2) * (totDataSem1j + totDataSme1j))/totDataSee1jSub if totDataSee1jSub > 0 else 1.0
###         totDataZmm1jSub = max(0, totDataZmm1j - 0.5*kmmData * (totDataZem1j + totDataZme1j))
###         totDataZee1jSub = max(0, totDataZee1j - 0.5*keeData * (totDataZem1j + totDataZme1j))
###         totDataZmm1jErr = sqrt(totDataZmm1j + pow(0.5*kmmData,2) * (totDataZem1j + totDataZme1j))/totDataZmm1jSub if totDataZmm1jSub > 0 else 1.0
###         totDataZee1jErr = sqrt(totDataZee1j + pow(0.5*keeData,2) * (totDataZem1j + totDataZme1j))/totDataZee1jSub if totDataZee1jSub > 0 else 1.0
###         Aee1jData = (totDataSee1jSub+1)/(totDataZee1jSub+1);   Aee1jDataErr = sqrt(1.0/(totDataZee1jSub+1) + 1.0/(totDataSee1jErr+1));
###         Amm1jData = (totDataSmm1jSub+1)/(totDataZmm1jSub+1);   Amm1jDataErr = sqrt(1.0/(totDataZmm1jSub+1) + 1.0/(totDataSmm1jErr+1));
###         totSmm1jSub = max(0, totSmm1j - 0.5*kmmMC * (totSem1j + totSme1j))
###         totSee1jSub = max(0, totSee1j - 0.5*keeMC * (totSem1j + totSme1j))
###         totSmm1jErr = sqrt(totSmm1j + pow(0.5*kmmMC,2) * (totSem1j + totSme1j))/totSmm1jSub if totSmm1jSub > 0 else 1.0
###         totSee1jErr = sqrt(totSee1j + pow(0.5*keeMC,2) * (totSem1j + totSme1j))/totSee1jSub if totSee1jSub > 0 else 1.0
###         totZmm1jSub = max(0, totZmm1j - 0.5*kmmMC * (totZem1j + totZme1j))
###         totZee1jSub = max(0, totZee1j - 0.5*keeMC * (totZem1j + totZme1j))
###         totZmm1jErr = sqrt(totZmm1j + pow(0.5*kmmMC,2) * (totZem1j + totZme1j))/totZmm1jSub if totZmm1jSub > 0 else 1.0
###         totZee1jErr = sqrt(totZee1j + pow(0.5*keeMC,2) * (totZem1j + totZme1j))/totZee1jSub if totZee1jSub > 0 else 1.0
###         Aee1j = (totSee1jSub+avgweightSEE)/(totZee1jSub+avgweightEE);   Aee1jErr = sqrt(1.0/(totZee1jSub/avgweightEE+1) + 1.0/(totSee1jErr/avgweightSEE+1));
###         Amm1j = (totSmm1jSub+avgweightSMM)/(totZmm1jSub+avgweightMM);   Amm1jErr = sqrt(1.0/(totZmm1jSub/avgweightEE+1) + 1.0/(totSmm1jErr/avgweightSMM+1));
###         gZmm1jMC.SetPoint(b, 1.0*b, totZmm1jSub); gZmm1jMC.SetPointError(b, 0, 0, totZmm1jErr*totZmm1jSub, totZmm1jErr*totZmm1jSub);
###         gZee1jMC.SetPoint(b, 1.0*b, totZee1jSub); gZee1jMC.SetPointError(b, 0, 0, totZee1jErr*totZee1jSub, totZee1jErr*totZee1jSub);
###         gSmm1jMC.SetPoint(b, 1.0*b, totSmm1jSub); gSmm1jMC.SetPointError(b, 0, 0, totSmm1jErr*totSmm1jSub, totSmm1jErr*totSmm1jSub);
###         gSee1jMC.SetPoint(b, 1.0*b, totSee1jSub); gSee1jMC.SetPointError(b, 0, 0, totSee1jErr*totSee1jSub, totSee1jErr*totSee1jSub);
###         gAmm1jMC.SetPoint(b, 1.0*b,    Amm1j);    gAmm1jMC.SetPointError(b, 0, 0, Amm1j*Amm1jErr, Amm1j*Amm1jErr);
###         gAee1jMC.SetPoint(b, 1.0*b,    Aee1j);    gAee1jMC.SetPointError(b, 0, 0, Aee1j*Aee1jErr, Aee1j*Aee1jErr);
###         gZmm1jData.SetPoint(b, 1.0*b, totDataZmm1jSub); gZmm1jData.SetPointError(b, 0, 0, totDataZmm1jErr, totDataZmm1jErr);
###         gZee1jData.SetPoint(b, 1.0*b, totDataZee1jSub); gZee1jData.SetPointError(b, 0, 0, totDataZee1jErr, totDataZee1jErr);
###         gSmm1jData.SetPoint(b, 1.0*b, totDataSmm1jSub); gSmm1jData.SetPointError(b, 0, 0, totDataSmm1jErr, totDataSmm1jErr);
###         gSee1jData.SetPoint(b, 1.0*b, totDataSee1jSub); gSee1jData.SetPointError(b, 0, 0, totDataSee1jErr, totDataSee1jErr);
###         gAmm1jData.SetPoint(b, 1.0*b,       Amm1jData); gAmm1jData.SetPointError(b, 0, 0, Amm1jData*Amm1jDataErr, Amm1jData*Amm1jDataErr);
###         gAee1jData.SetPoint(b, 1.0*b,       Aee1jData); gAee1jData.SetPointError(b, 0, 0, Aee1jData*Aee1jDataErr, Aee1jData*Aee1jDataErr);
###         print "\t %2d   %9.1f %8.2f  %5.3f %4.2f |  %9.1f %8.2f  %5.3f %4.2f |  %9.1f %8.2f  %5.3f %4.2f |  %9.1f %8.2f  %5.3f %4.2f" % (
###                 b,
###                 totZmm1jSub, totSmm1jSub, Amm1j, min(1.0, Amm1jErr),
###                 totDataZmm1jSub, totDataSmm1jSub, Amm1jData, min(1.0, Amm1jDataErr),
###                 totZee1jSub, totSee1jSub, Aee1j, min(1.0, Aee1jErr),
###                 totDataZee1jSub, totDataSee1jSub, Aee1jData, min(1.0, Aee1jDataErr),
###               )
###         totZmm1j -= plotZ1j[MM][1].GetBinContent(b+1); totSmm1j -= plot1j[MM][1].GetBinContent(b+1)
###         totZee1j -= plotZ1j[EE][1].GetBinContent(b+1); totSee1j -= plot1j[EE][1].GetBinContent(b+1)
###         totZem1j -= plotZ1j[EM][1].GetBinContent(b+1); totSem1j -= plot1j[EM][1].GetBinContent(b+1)
###         totZme1j -= plotZ1j[ME][1].GetBinContent(b+1); totSme1j -= plot1j[ME][1].GetBinContent(b+1)
###         totDataZmm1j -= plotDataZ1j[MM][1].GetBinContent(b+1); totDataSmm1j -= plotData1j[MM][1].GetBinContent(b+1)
###         totDataZee1j -= plotDataZ1j[EE][1].GetBinContent(b+1); totDataSee1j -= plotData1j[EE][1].GetBinContent(b+1)
###         totDataZem1j -= plotDataZ1j[EM][1].GetBinContent(b+1); totDataSem1j -= plotData1j[EM][1].GetBinContent(b+1)
###         totDataZme1j -= plotDataZ1j[ME][1].GetBinContent(b+1); totDataSme1j -= plotData1j[ME][1].GetBinContent(b+1)
###     gZmm1jMC.SetName("gZmm1jMC"); gSmm1jMC.SetName("gSmm1jMC"); gAmm1jMC.SetName("gAmm1jMC"); 
###     gZee1jMC.SetName("gZee1jMC"); gSee1jMC.SetName("gSee1jMC"); gAee1jMC.SetName("gAee1jMC");
###     gZmm1jData.SetName("gZmm1jData"); gSmm1jData.SetName("gSmm1jData"); gAmm1jData.SetName("gAmm1jData"); 
###     gZee1jData.SetName("gZee1jData"); gSee1jData.SetName("gSee1jData"); gAee1jData.SetName("gAee1jData");
###     gZmm1jMC.Write(); gZee1jMC.Write(); gSmm1jMC.Write(); gSee1jMC.Write(); gAmm1jMC.Write(); gAee1jMC.Write();
###     gZmm1jData.Write(); gZee1jData.Write(); gSmm1jData.Write(); gSee1jData.Write(); gAmm1jData.Write(); gAee1jData.Write();
### 
### bg0jmm = alpha0jmm*Z0jmmSubData; bg0jee = alpha0jee*Z0jeeSubData;
### if do1:
###     bg1jee = alpha1jee*Z1jeeSubData; bg1jmm = alpha1jmm*Z1jmmSubData;
### 
### print "\n"
### print "\n --- Efficiency factors, from Z peak with no MET cut.  --- "
### print "kee(Data) = sqrt(N(ee)/N(mm)) = %6.4f\tkee(MC) = %6.4f" % (keeData ,keeMC)
### print "kmm(Data) = sqrt(N(mm)/N(ee)) = %6.4f\tkmm(MC) = %6.4f" % (kmmData ,kmmMC)
### 
### if options.subVV:
###     print "\n --- Yields in the peak, with MET cut and VV subtraction --- "
###     print "N(Zmm 0j, Data) = (%5.1f - %5.3f/2 * %5.1f) - %5.2f = %5.2f +/- %4.2f" % (Z0jmmData, kmmData, Z0jxxData, Z0jmmSubMCVV, Z0jmmSubData, Z0jmmSubDataErr)
###     print "N(Zee 0j, Data) = (%5.1f - %5.3f/2 * %5.1f) - %5.2f = %5.2f +/- %4.2f" % (Z0jeeData, keeData, Z0jxxData, Z0jeeSubMCVV, Z0jeeSubData, Z0jeeSubDataErr)
###     print "N(Zmm 0j, VV)   = (%5.1f - %5.3f/2 * %5.1f)         = %5.2f +/- %4.2f" % (Z0jmmMCVV, kmmMC  , Z0jxxMCVV,               Z0jmmSubMCVV, Z0jmmSubMCVVErr  )
###     print "N(Zee 0j, VV)   = (%5.1f - %5.3f/2 * %5.1f)         = %5.2f +/- %4.2f" % (Z0jeeMCVV, keeMC  , Z0jxxMCVV,               Z0jeeSubMCVV, Z0jeeSubMCVVErr  )
###     print "N(Zmm 0j, DY)   = (%5.1f - %5.3f/2 * %5.1f)         = %5.2f +/- %4.2f" % (Z0jmmMC  , kmmMC  , Z0jxxMC  ,               Z0jmmSubMC  , Z0jmmSubMCErr  )
###     print "N(Zee 0j, DY)   = (%5.1f - %5.3f/2 * %5.1f)         = %5.2f +/- %4.2f" % (Z0jeeMC  , keeMC  , Z0jxxMC  ,               Z0jeeSubMC  , Z0jeeSubMCErr  )
###     if do1:
###         print "N(Zmm 1j, Data) = (%5.1f - %5.3f/2 * %5.1f) - %5.2f = %5.2f +/- %4.2f" % (Z1jmmData, kmmData, Z1jxxData, Z1jmmSubMCVV, Z1jmmSubData, Z1jmmSubDataErr)
###         print "N(Zee 1j, Data) = (%5.1f - %5.3f/2 * %5.1f) - %5.2f = %5.2f +/- %4.2f" % (Z1jeeData, keeData, Z1jxxData, Z1jeeSubMCVV, Z1jeeSubData, Z1jeeSubDataErr)
###         print "N(Zmm 1j, VV)   = (%5.1f - %5.3f/2 * %5.1f)         = %5.2f +/- %4.2f" % (Z1jmmMCVV, kmmMC  , Z1jxxMCVV,               Z1jmmSubMCVV, Z1jmmSubMCVVErr  )
###         print "N(Zee 1j, VV)   = (%5.1f - %5.3f/2 * %5.1f)         = %5.2f +/- %4.2f" % (Z1jeeMCVV, keeMC  , Z1jxxMCVV,               Z1jeeSubMCVV, Z1jeeSubMCVVErr  )
###         print "N(Zmm 1j, DY)   = (%5.1f - %5.3f/2 * %5.1f)         = %5.2f +/- %4.2f" % (Z1jmmMC  , kmmMC  , Z1jxxMC  ,               Z1jmmSubMC  , Z1jmmSubMCErr  )
###         print "N(Zee 1j, DY)   = (%5.1f - %5.3f/2 * %5.1f)         = %5.2f +/- %4.2f" % (Z1jeeMC  , keeMC  , Z1jxxMC  ,               Z1jeeSubMC  , Z1jeeSubMCErr  )
### else:
###     print "\n --- Yields in the peak, with MET cut                  --- "
###     print "N(Zmm 0j, Data) = (%5.1f - %5.3f/2 * %5.1f) = %5.2f +/- %4.2f" % (Z0jmmData, kmmData, Z0jxxData, Z0jmmSubData, Z0jmmSubDataErr)
###     print "N(Zmm 0j, Sim.) = (%5.1f - %5.3f/2 * %5.1f) = %5.2f +/- %4.2f" % (Z0jmmMC  , kmmMC  , Z0jxxMC  , Z0jmmSubMC  , Z0jmmSubMCErr  )
###     print "N(Zee 0j, Data) = (%5.1f - %5.3f/2 * %5.1f) = %5.2f +/- %4.2f" % (Z0jeeData, keeData, Z0jxxData, Z0jeeSubData, Z0jeeSubDataErr)
###     print "N(Zee 0j, Sim.) = (%5.1f - %5.3f/2 * %5.1f) = %5.2f +/- %4.2f" % (Z0jeeMC  , keeMC  , Z0jxxMC  , Z0jeeSubMC  , Z0jeeSubMCErr  )
###     if do1:
###         print "N(Zmm 1j, Data) = (%5.1f - %5.3f/2 * %5.1f) = %5.2f +/- %4.2f" % (Z1jmmData, kmmData, Z1jxxData, Z1jmmSubData, Z1jmmSubDataErr)
###         print "N(Zmm 1j, Sim.) = (%5.1f - %5.3f/2 * %5.1f) = %5.2f +/- %4.2f" % (Z1jmmMC  , kmmMC  , Z1jxxMC  , Z1jmmSubMC  , Z1jmmSubMCErr  )
###         print "N(Zee 1j, Data) = (%5.1f - %5.3f/2 * %5.1f) = %5.2f +/- %4.2f" % (Z1jeeData, keeData, Z1jxxData, Z1jeeSubData, Z1jeeSubDataErr)
###         print "N(Zee 1j, Sim.) = (%5.1f - %5.3f/2 * %5.1f) = %5.2f +/- %4.2f" % (Z1jeeMC  , keeMC  , Z1jxxMC  , Z1jeeSubMC  , Z1jeeSubMCErr  )
### 
### 
### print "\n --- Extrapolation factors to final selection, from MC --- "
### print "alpha mm, 0j, bare: %6.2f / %6.2f = %6.2f / %6.2f (effective events) " % (
###          report0jMC[-1][1][MM][1], Z0jmmMC, report0jMC[-1][1][MM][1]/avgweightMM, Z0jmmMC/avgweightMM)
### print "alpha ee, 0j, bare: %6.2f / %6.2f = %6.2f / %6.2f (effective events) " % (
###          report0jMC[-1][1][EE][1], Z0jeeMC, report0jMC[-1][1][EE][1]/avgweightEE, Z0jeeMC/avgweightEE)
### if do1:
###     print "alpha mm, 1j, bare: %6.2f / %6.2f = %6.2f / %6.2f (effective events) " % (
###              report1jMC[-1][1][MM][1], Z1jmmMC, report1jMC[-1][1][MM][1]/avgweightMM, Z1jmmMC/avgweightMM)
###     print "alpha ee, 1j, bare: %6.2f / %6.2f = %6.2f / %6.2f (effective events) " % (
###              report1jMC[-1][1][EE][1], Z1jeeMC, report1jMC[-1][1][EE][1]/avgweightEE, Z1jeeMC/avgweightEE)
### print "alpha mm, 0j, regularized: %7.4f   (rel. uncertainty %6.2f%%)" % ( alpha0jmm, alpha0jmmErr*100 )
### print "alpha ee, 0j, regularized: %7.4f   (rel. uncertainty %6.2f%%)" % ( alpha0jee, alpha0jeeErr*100 )
### if do1:
###     print "alpha ee, 1j, regularized: %7.4f   (rel. uncertainty %6.2f%%)" % ( alpha1jee, alpha1jeeErr*100 )
###     print "alpha mm, 1j, regularized: %7.4f   (rel. uncertainty %6.2f%%)" % ( alpha1jmm, alpha1jmmErr*100 )
### 
### print "\n --- Estimated background in the final selection       --- "
### print "N(Z->mm, all cuts, 0j) = %7.3f +/- %7.3f (stat. data) +/- %7.3f (stat. mc)" % (bg0jmm, alpha0jmm*Z0jmmSubDataErr, bg0jmm*alpha0jmmErr)
### print "N(Z->ee, all cuts, 0j) = %7.3f +/- %7.3f (stat. data) +/- %7.3f (stat. mc)" % (bg0jee, alpha0jee*Z0jeeSubDataErr, bg0jee*alpha0jeeErr)
### if do1:
###     print "N(Z->mm, all cuts, 1j) = %7.3f +/- %7.3f (stat. data) +/- %7.3f (stat. mc)" % (bg1jmm, alpha1jmm*Z1jmmSubDataErr, bg1jmm*alpha1jmmErr)
###     print "N(Z->ee, all cuts, 1j) = %7.3f +/- %7.3f (stat. data) +/- %7.3f (stat. mc)" % (bg1jee, alpha1jee*Z1jeeSubDataErr, bg1jee*alpha1jeeErr)
### 
