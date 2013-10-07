/**
  USAGE: 
    DATA-only:    root.exe -b -l -q TnP_ICHEP_EleID_data_all.root "plotEleID.cxx(\"data_all\")"
    MC data-like: root.exe -b -l -q TnP_ICHEP_EleID_datalike_mc.root  "plotEleID.cxx(\"datalike_mc\")"
    DATA+MC:      root.exe -b -l -q TnP_ICHEP_EleID_data_all.root TnP_ICHEP_EleID_datalike_mc.root  "plotEleID.cxx(\"data_vs_mc\")"

    Tracker Probes vs Calo Probes:
        root.exe -b -l -q TnP_ICHEP_EleID_FromTK_data_all.root TnP_ICHEP_EleID_data_all.root  "plotEleID_ICHEP.cxx(\"data_all\",2)"
        root.exe -b -l -q TnP_ICHEP_EleID_FromTK_datalike_mc.root TnP_ICHEP_EleID_datalike_mc.root  "plotEleID_ICHEP.cxx(\"datalike_mc\",2)"

  REQUIRES:
   1) mkdir -p plots_ichep_dev/eleid/ plots_ichep_dev/eleid_tk/ plots_ichep_dev/eleid_tk_vs_cal/
   2) provide a suitable "tdrStyle.cc" macro or similar
      (by default, it's taken from ~/cpp/tdrstyle.cc;
       if you need one you might want to grab ~gpetrucc/cpp/tdrstyle.cc)
*/
#include <TCanvas.h>
#include <TPad.h>
#include <TFile.h>
TString prefix = "plots/eleid/";
TString basedir  = "tpTreeElEl";
TString basedir2 = "tpTreeElEl";
TFile *ref = 0;

#include "plotUtil.cxx"


void plotEleIDData(bool refFit=1) {
    TString eta = "eta";
    //TString eta = "abseta";
    retitle = "Efficiency";
    const int nplots = 6;
    const char *ids[nplots]    = { "all",     "id",              "iso"       , "hlt"     , "idNoIso"  , "isoNoId"};
    const char *titles[nplots] = { "Iso*Id",  "Identification",  "Isolation",  "Trigger" , "ID from Gsf"  , "Iso from Gsf"};
    const int nxvar = 6;
    const char *xnam[nxvar]   = { "pt_eta", "average", "eta",   "vtx",    "njet",          "run" };
    const char *xvar[nxvar]   = { "pt",     "average", "eta",   "nVtx",   "pair_nJet",    "run" };
    const char *xtitle[nxvar] = { "pt",     "average", "#eta",  "N(vtx)", "N(jet)",       "run number" };
    for (int i = 0; i < nplots; ++i) {
        for (int j = 0; j < nxvar; ++j) {
            TString idname(TString::Format("%s_%s", ids[i], xnam[j]));
            retitle = TString(titles[i])+" ele efficiency";
            TDirectory *fit_dir = gFile->GetDirectory(basedir+"/"+idname+"/");
            TString xv = strcmp(xvar[j],"average")!=0 ? xvar[j] : "pt";
            if (fit_dir == 0) {
                std::cerr << "Measurement " << idname << "(" << titles[i] << ") NOT AVAILABLE" << std::endl;
                continue;
            }
            yMin = 0.5;
            TDirectory *ref_dir = (ref == 0 ? 0 : ref->GetDirectory(basedir2+"/"+idname+"/"));
            if (strstr("eta|phi|run", xnam[j]) == 0) { // then we have barrel and endcaps
                printDataSet(gFile->GetDirectory(basedir+"/"+idname), idname, xv, "abseta");
                if (ref_dir != 0) {
                    printDataSet(ref->GetDirectory(basedir2+"/"+idname),  idname+"_reference", xv, "abseta", (refFit ? "fit_eff" : "cnt_eff"));
                    if (refFit) {
                        refstack(fit_dir, ref_dir, idname+"_barrel",  xv+"_PLOT_abseta_bin0_");
                        refstack(fit_dir, ref_dir, idname+"_endcap",  xv+"_PLOT_abseta_bin1_");
                    } else {
                        mcstack(fit_dir, ref_dir, idname+"_barrel",   xv+"_PLOT_abseta_bin0_");
                        mcstack(fit_dir, ref_dir, idname+"_endcap",   xv+"_PLOT_abseta_bin1_");
                    }
               } else {
                    single(fit_dir, idname, xv+"_PLOT");
                }
                if (strcmp(xnam[j], "pt_eta") == 0) {
                    doCanvas(fit_dir, 1, 10, idname+"_barrel_pt%d",  "abseta_bin0__pair_nJet_bin0__pt_bin_%d");
                    doCanvas(fit_dir, 1, 10, idname+"_endcap_pt%d",  "abseta_bin1__pair_nJet_bin0__pt_bin_%d");
                } else if (strcmp(xnam[j], "average") == 0) {
                    doCanvas(fit_dir, 1, 1, idname+"_barrel",  "abseta_bin0__");
                    doCanvas(fit_dir, 1, 1, idname+"_endcap",  "abseta_bin1__");
                }
            } else {
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
                    doCanvas(fit_dir, 1, 10,  idname+"_bin_%d",  xv+"_bin%d__");
                }
            }
        }
    }

}

void plotEleID(bool refFit=1,TString dir="eleid") {
    prefix = "plots/"+dir+"/";
    gSystem->mkdir(prefix,true);

    gROOT->ProcessLine(".x /afs/cern.ch/user/g/gpetrucc/cpp/tdrstyle.cc");
    gStyle->SetOptStat(0);
    c1 = new TCanvas("c1","c1");

    if (gROOT->GetListOfFiles()->GetEntries() == 2) {
        ref = (TFile *) gROOT->GetListOfFiles()->At(1);
        ((TFile*) gROOT->GetListOfFiles()->At(0))->cd();
    }

    autoScale = true;
    doRatioPlot = true;
    doDiffPlot = false;
    //doPdf    = true;
    //doSquare = true; 
    doSquare = false; 
    yMax = 1.05; yMin = 0.5;
    datalbl = "Data, 35/pb";
    reflbl  = "MC 3.8.X";
    if (dir == "eleid_mc") {
        datalbl = "MC, Fit";
        reflbl  = "MC, Count";
    }
    if (dir == "eleid_test") {
        datalbl = "Data, Fit";
        reflbl  = "MC, Count";
    }
    //preliminary = "CMS Preliminary,   #sqrt{s} = 7 TeV";
    plotEleIDData(refFit);
}


