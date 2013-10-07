import FWCore.ParameterSet.Config as cms

# from ROOT import ROOT
# from ROOT import TH1F
# from ROOT import TFile
# from ROOT import TTree
# from ROOT import gDirectory
import ROOT

# grab tree from file, pass the file name
# and the tree name/path in the file

colorList = [ 
    ROOT.kBlack,
    ROOT.kRed,
    ROOT.kBlue,
    ROOT.kGreen+2,
    ROOT.kOrange+7,
    ROOT.kMagenta+2,
    ROOT.kGray,
    ROOT.kRed+3,
    ROOT.kAzure+5,
    ROOT.kGreen+4,
    ROOT.kOrange+4,
    ROOT.kMagenta+7
]

legendInfo = [
    [ 0.6,0.7,0.9,0.9,'l' ],
]

myGlobalList = []

def getTree(fn,tn):
    myGlobalList.append(ROOT.TFile(fn,"READ"))
    t = myGlobalList[-1].Get(tn)
    return t

# pass the file name, the tree name/path in the file, 
# number of bins, min max range, the value to plot, and 
# a selection if you want to.  also push over/underflow
# in if you want, on by default; normalize to unit area
# by default also
def getHistFromFile(fn,tn,nbins,min,max,toPlot,sel,ou=True,norm=True):

    #make unique name
    histName = (
        'h_' +
        fn.replace('/','_').replace('.','_').replace('.root','') + 
        tn.replace('/','_').replace('.','_').replace('.root','') + 
        toPlot.translate(None,'`~!@#$%^&*()-+= {}[]\|;:"?/.>,<') +
        sel.GetTitle().translate(None,'`~!@#$%^&*()-+= {}[]\|;:"?/.>,<') )

    #get tree
    t = getTree(fn,tn)

    # if put overflow and underflow in plot
    if ou:
        scale = (max-min)/nbins
        toPlot = 'max(min('+toPlot+','+str(max-0.01*scale)+'),'+str(min+0.01*scale)+')'

    #fill histogram
    binString = '('+str(nbins)+','+str(min)+','+str(max)+')'
    t.Draw(toPlot+'>>'+histName+binString,sel,"GOFF")
    
    h = getattr(ROOT,histName)
    h.Sumw2()
    #if normalize to unit area
    if norm and h.Integral() != 0:
        h.Scale(1.0/h.Integral())

    return h

def getIntegralHist(h1):
    
    h2 = ROOT.TH1F(h1.GetName()+'_eff', h1.GetName()+'_eff',
              h1.GetNbinsX(), h1.GetBinLowEdge(0),
              h1.GetBinLowEdge(h1.GetNbinsX()+1))

    for i in range(0,h1.GetNbinsX()+1): 
        h2.SetBinContent(i,h1.GetIntegral()[i])


    h2.SetTitle(h1.GetTitle())
    h2.GetXaxis().SetTitle(h1.GetXaxis().GetTitle())
    h2.GetYaxis().SetTitle(h1.GetYaxis().GetTitle())
    h2.SetLineWidth(h1.GetLineWidth())
    h2.SetLineStyle(h1.GetLineStyle())
    h2.SetLineColor(h1.GetLineColor())

    return h2

def getClopGraph(h1,h2):
    
    g = ROOT.TGraphAsymmErrors(h1,h2)

#     alpha = (1.0 - .68540158589942957)/2;
#     for i in range(0,g.GetN()):
#         num = h1.GetBinContent(i+1)
#         den = h2.GetBinContent(i+1)
#         eff = 0.0 if den == 0   else float(num)/float(den);
#         lo  = 0.0 if num == 0   else ROOT.Math.beta_quantile(   alpha, num,   den-num+1 )
#         hi  = 1.0 if den == num else ROOT.Math.beta_quantile( 1-alpha, num+1, den-num   )
#         print lo
#         print hi
#         g.SetPoint(i,eff,g.GetY()[i])
#         g.SetPointEYhigh(i,hi)
#         g.SetPointEYlow(i,lo)

    return g


def getEffCurve(hSig,hBkg):
    
    if hSig.GetSize() != hBkg.GetSize():
    	print "Array lengths are not the same!!!" 
    	return ROOT.TGraph()
    
    i=0;
    while hSig.GetArray()[i] == 0:
        i+=1

#     return ROOT.TGraph(hSig.GetSize()-i-1,hBkg.GetArray()[i:],hSig.GetArray()[i:]);
    return ROOT.TGraph(hSig.GetSize()-1,hBkg.GetArray(),hSig.GetArray());

def getEffCurveFromFiles(fn1,tn1,toPlot1,sel1,fn2,tn2,toPlot2,sel2,nbins,min,max,ou=True,norm=True):
    h1 = getHistFromFile(fn1,tn1,nbins,min,max,toPlot1,sel1,ou,norm)
    h2 = getHistFromFile(fn2,tn2,nbins,min,max,toPlot2,sel2,ou,norm)
    h3 = getIntegralHist(h1)
    h4 = getIntegralHist(h2)
    return getEffCurve(h3,h4)

def getLines(fn1,tn1,toPlot1,sel1,sel11,fn2,tn2,toPlot2,sel2,sel21,nbins,min,max,ou=True,norm=True):
    h1 = getHistFromFile(fn1,tn1,nbins,min,max,toPlot1,sel1,True,False)
    h2 = getHistFromFile(fn2,tn2,nbins,min,max,toPlot2,sel2,True,False)
    h3 = getHistFromFile(fn1,tn1,nbins,min,max,toPlot1,sel11,True,False)
    h4 = getHistFromFile(fn2,tn2,nbins,min,max,toPlot2,sel21,True,False)
    x = h4.Integral()/h2.Integral()
    y = h3.Integral()/h1.Integral()
    return [ROOT.TLine(0,y,x,y),ROOT.TLine(x,0,x,y)]


def formatHist(h, i, opt=''):
    h.SetLineColor( colorList[i] )
    h.SetLineWidth( 3 )
    h.SetLineStyle( 2 if opt.find('bkg') != -1 else 1 )

def setupLabels(h,x,y):
    h.GetXaxis().SetTitle(x)
    h.GetYaxis().SetTitle(y)
    h.GetYaxis().SetRangeUser(8.5e-5,1.15)

ls = []
def addLegend(sigHists,i):
    ls.append( ROOT.TLegend(
        legendInfo[i][0],
        legendInfo[i][1],
        legendInfo[i][2],
        legendInfo[i][3]
    ))
    for s in sigHists:
        ls[-1].AddEntry(sigHists[s],s,legendInfo[i][4])
    ls[-1].SetLineWidth(0)
    ls[-1].SetFillColor(0)
    ls[-1].SetLineColor(0)
    ls[-1].Draw()


def createAllHists(c1,plot):
    
    if c1 == None:
        c1 = ROOT.TCanvas("c1","c1",1200,1200)
        c1.SetWindowSize(1200,1200)
        c1.Divide(2,2)

    c1.cd(1).SetLogy()

    sigHists = {}
    bkgHists = {}
    for p in plot[0]:
        sigHists[p[2]] = getHistFromFile(p[0],p[3],p[4],p[5],p[6],p[7],p[8])
        bkgHists[p[2]] = getHistFromFile(p[1],p[3],p[4],p[5],p[6],p[7],p[9])

    temp = sigHists.popitem()
    formatHist(temp[1],0)
    formatHist(bkgHists[temp[0]],0,'bkg')
    setupLabels(temp[1],plot[1],plot[2])
    temp[1].Draw('hist')
    bkgHists[temp[0]].Draw('same hist')

    i = 1
    for sig in sigHists:
        formatHist(sigHists[sig],i)
        formatHist(bkgHists[sig],i,'bkg')
        sigHists[sig].Draw('hist same')
        bkgHists[sig].Draw('hist same')
        i = i+1
    sigHists[temp[0]] = temp[1]
    addLegend(sigHists,0)


#     for bkg in bkgHists:
#         bkgHists[sig].Draw('same')

    return c1


