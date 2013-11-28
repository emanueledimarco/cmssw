/**
  USAGE: 
    DATA-only:    root.exe -b -l -q TnP_ICHEP_MuonID_data_all.root "plotMuonID.cxx(\"data_all\")"
    MC data-like: root.exe -b -l -q TnP_ICHEP_MuonID_datalike_mc.root  "plotMuonID.cxx(\"datalike_mc\")"
    DATA+MC:      root.exe -b -l -q TnP_ICHEP_MuonID_data_all.root TnP_ICHEP_MuonID_datalike_mc.root  "plotMuonID.cxx(\"data_vs_mc\")"

    Tracker Probes vs Calo Probes:
        root.exe -b -l -q TnP_ICHEP_MuonID_FromTK_data_all.root TnP_ICHEP_MuonID_data_all.root  "plotMuonID_ICHEP.cxx(\"data_all\",2)"
        root.exe -b -l -q TnP_ICHEP_MuonID_FromTK_datalike_mc.root TnP_ICHEP_MuonID_datalike_mc.root  "plotMuonID_ICHEP.cxx(\"datalike_mc\",2)"

  REQUIRES:
   1) mkdir -p plots_ichep_dev/muonid/ plots_ichep_dev/muonid_tk/ plots_ichep_dev/muonid_tk_vs_cal/
   2) provide a suitable "tdrStyle.cc" macro or similar
      (by default, it's taken from ~/cpp/tdrstyle.cc;
       if you need one you might want to grab ~gpetrucc/cpp/tdrstyle.cc)
*/
#include <TCanvas.h>
#include <TPad.h>
#include <TFile.h>
TString prefix = "plots/muonid/";
TString basedir  = "tpTreeMuMu";
TString basedir2 = "tpTreeMuMu";
TFile *ref = 0;

#include "plotUtil.cxx"


void plotMuonIDData(bool refFit=1) {
    TString eta = "eta";
    //TString eta = "abseta";
    retitle = "Efficiency";
    const int nplots = 6;
    const char *ids[nplots]    = { "all",        "id",              "iso"       , "hlt"     , "mu15",           "mu9",         };
    const char *titles[nplots] = { "Selection",  "Identification",  "Isolation",  "Trigger",  "Trigger (Mu15)", "Trigger (Mu9)"};
    const int nxvar = 7;
    const char *xnam[nxvar]   = { "abseta", "average", "eta",  "phi",  "vtx",    "jet",          "run" };
    const char *xvar[nxvar]   = { "abseta", "average", "eta",  "phi",  "nVtx",   "pair_nJet",    "run" };
    const char *xtitle[nxvar] = { "|#eta|", "average", "#eta", "#phi", "N(vtx)", "N(jet)",       "run number" };
    for (int i = 0; i < nplots; ++i) {
        for (int j = 0; j < nxvar; ++j) {
            TString idname(TString::Format("%s_%s", ids[i], xnam[j]));
            retitle = TString(titles[i])+" efficiency";
            TDirectory *fit_dir = gFile->GetDirectory(basedir+"/"+idname+"/");
            TString xv = strcmp(xvar[j],"average")!=0 ? xvar[j] : "abseta";
            if (fit_dir == 0) {
                std::cerr << "Measurement " << idname << "(" << titles[i] << ") NOT AVAILABLE" << std::endl;
                continue;
            }
            if (i >= 3) {
                yMin = 0.0;
                doRatioPlot = false;
            } else {
                yMin = 0.8;
                doRatioPlot = true;
            }
            printDataSet(gFile->GetDirectory(basedir+"/"+idname), idname, xv, "pt");
            TDirectory *ref_dir = (ref == 0 ? 0 : ref->GetDirectory(basedir2+"/"+idname+"/"));
            if (ref_dir != 0) {
                printDataSet(ref->GetDirectory(basedir2+"/"+idname),  idname+"_reference", xv, "pt", (refFit ? "fit_eff" : "cnt_eff"));
                if (refFit) {
                    refstack(fit_dir, ref_dir, idname,  xv+"_PLOT");
                } else {
                    mcstack(fit_dir, ref_dir, idname,   xv+"_PLOT");
                }
           } else {
                single(fit_dir, idname, xv+"_PLOT");
            }
            if (1) {
                //doCanvas(fit_dir, 1, 10,  idname+"_bin_%d",  xv+"_bin%d__");
            }
        }
    }

}

void plotMuonID(bool refFit=1,TString dir="muonid") {
    prefix = "plots/"+dir+"/";
    gSystem->mkdir(prefix,true);

    gROOT->ProcessLine(".x /afs/cern.ch/user/g/gpetrucc/cpp/tdrstyle.cc");
    gStyle->SetOptStat(0);
    gStyle->SetTitleYOffset(1.45);
    c1 = new TCanvas("c1","c1");
    c1->SetLeftMargin(0.30);

    if (gROOT->GetListOfFiles()->GetEntries() == 2) {
        ref = (TFile *) gROOT->GetListOfFiles()->At(1);
        ((TFile*) gROOT->GetListOfFiles()->At(0))->cd();
    }

    autoScale = true;
    doRatioPlot = true;
    doDiffPlot = false;
    doPdf    = true;
    doSquare = true; 
    //doSquare = false; 
    yMax = 1.05; yMin = 0.5;
    cropPointsWithErrorsAbove = 0.10;

    datalbl = "Data, 35/pb";
    reflbl  = "Simulation";
    if (dir == "muonid_mc") {
        datalbl = "MC, Fit";
        reflbl  = "MC, Count";
    }
    if (dir == "muonid_test") {
        datalbl = "Data, Fit";
        reflbl  = "Data, Count";
    }

    preliminary = "CMS Preliminary, #sqrt{s} = 7 TeV";
    plotMuonIDData(refFit);
}


