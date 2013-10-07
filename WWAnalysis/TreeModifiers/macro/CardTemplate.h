#ifndef CARDTEMPLETE_H
#define CARDTEMPLETE_H

std::string createCardTemplate(bool do7TeV, int channel, bool do1D, std::string workspacefilename, bool isLowMass=true) {

    std::string card = "";
        card += "imax 1\n";
        card += "jmax *\n";
        card += "kmax *\n";
        card += "------------\n";
        card += std::string("shapes * * ") + workspacefilename + " w:$PROCESS\n";
        card += "------------\n";
        card += "bin         BIN\n";
        card += "observation OBS\n";
        card += "------------\n";
        card += "bin     BIN           BIN            BIN           BIN            BIN            BIN             BIN             BIN\n";
        card += "process ggH           qqH            WH            ZH             ttH            bkg_qqzz        bkg_ggzz        bkg_zjets\n";
        card += "process -5            -4             -3            -2             -1             1               2               3\n";
        card += "rate    1             1              1             1              1              BKG_QQZZ_YIELD  BKG_GGZZ_YIELD  BKG_ZJETS_YIELD\n";
        card += "------------\n";
        card += "## mass window [LOW_EDGE, HIGH_EDGE]\n";

    if (do7TeV) {        
        card += "lumi_7TeV                 lnN        1.022   1.022   1.022    1.022   1.022    1.022    1.022    -\n";
    }
    else {        
        card += "lumi_8TeV                 lnN        1.05    1.05    1.05     1.05    1.05     1.05     1.05     -\n";
    }
        card += "pdf_gg                    lnN        GGH_PDF -       -        -       -        -        GGZZ_PDF -\n";
        card += "pdf_qqbar                 lnN         -      VBF_PDF WHI_PDF  ZHI_PDF TTH_PDF  QQZZ_PDF -        -\n";
        card += "pdf_hzz4l_accept          lnN        1.02    1.02    1.02     1.02    1.02     -        -        -\n";
        card += "QCDscale_ggH              lnN        GGH_QCD -       -        -       -        -        -        -\n";
        card += "QCDscale_qqH              lnN        -       VBF_QCD -        -       -        -        -        -\n";
        card += "QCDscale_WH               lnN        -       -       WHI_QCD  -       -        -        -        -\n";
        card += "QCDscale_ZH               lnN        -       -       -        ZHI_QCD -        -        -        -\n";
        card += "QCDscale_ttH              lnN        -       -       -        -       TTH_QCD  -        -        -\n";
        card += "QCDscale_VV               lnN        -       -       -        -       -        QQZZ_QCD -        -\n";
        card += "QCDscale_ggVV             lnN        -       -       -        -       -        -        GGZZ_QCD -\n";
        card += "BRhiggs_ZZ4l              lnN        1.02    1.02    1.02     1.02    1.02     -        -        -\n";
    if (!isLowMass) {        
        card += "sig_gamma_err             param      0        0.05\n";
    }
    if (channel == 0) {
        card += "CMS_eff_m                 lnN        1.015   1.015   1.015    1.015   1.015    1.015    1.015    -\n";
        if (do7TeV) {
        card += "sig_mean_err_4mu_7TeV     param      0        0.001                     \n";
        card += "sig_sigma_err_4mu_7TeV    param      0        0.2                       \n";
        card += "CMS_hzz4mu_Zjets          lnN        -       -       -        -       -        -        -        ZX_SYST\n";
        }
        else {
        card += "sig_mean_err_4mu_8TeV     param      0        0.001                     \n";
        card += "sig_sigma_err_4mu_8TeV    param      0        0.2                       \n";
        card += "CMS_hzz4mu_Zjets          lnN        -       -       -        -       -        -        -        ZX_SYST\n";
        }
    }
    else if (channel == 1) {
        card += "CMS_eff_e                 lnN        1.01    1.01    1.01     1.01     1.01    1.01     1.01     -\n";
        if (do7TeV) {
        card += "sig_mean_err_4e_7TeV      param      0        0.002                     \n";
        card += "sig_sigma_err_4e_7TeV     param      0        0.2                       \n";
        card += "CMS_hzz4e_Zjets           lnN        -       -       -        -        -       -        -        ZX_SYST\n";
        }
        else {
        card += "sig_mean_err_4e_8TeV      param      0        0.002                     \n";
        card += "sig_sigma_err_4e_8TeV     param      0        0.2                       \n";
        card += "CMS_hzz4e_Zjets           lnN        -       -       -        -        -       -        -        ZX_SYST\n";
        }
    }
    else if (channel == 2) {
        card += "CMS_eff_m                 lnN        1.015   1.015   1.015    1.015   1.015    1.015    1.015    -\n";
        card += "CMS_eff_e                 lnN        1.01    1.01    1.01     1.01     1.01    1.01     1.01     -\n";
        if (do7TeV && isLowMass) {
        card += "sig_mean_err_4mu_7TeV     param      0        0.001                      \n";
        card += "sig_sigma_err_4mu_7TeV    param      0        0.2                        \n";
        card += "sig_mean_err_4e_7TeV      param      0        0.002                      \n";
        card += "sig_sigma_err_4e_7TeV     param      0        0.2                        \n";
        card += "CMS_hzz2e2mu_Zjets        lnN        -       -       -        -       -        -        -        ZX_SYST\n";
        }
        else if (!do7TeV && isLowMass) {
        card += "sig_mean_err_4mu_8TeV     param      0        0.001                     \n";
        card += "sig_sigma_err_4mu_8TeV    param      0        0.2                       \n";
        card += "sig_mean_err_4e_8TeV      param      0        0.002                     \n";
        card += "sig_sigma_err_4e_8TeV     param      0        0.2                       \n";
        card += "CMS_hzz2e2mu_Zjets        lnN        -       -       -        -       -        -        -        ZX_SYST\n";
        }
    }
    else if (channel == 3) {
        card += "CMS_eff_m                 lnN        1.015   1.015   1.015    1.015   1.015    1.015    1.015    -\n";
        card += "CMS_eff_e                 lnN        1.01    1.01    1.01     1.01     1.01    1.01     1.01     -\n";
        if (do7TeV && isLowMass) {
        card += "sig_mean_err_4mu_7TeV     param      0        0.001                      \n";
        card += "sig_sigma_err_4mu_7TeV    param      0        0.2                        \n";
        card += "sig_mean_err_4e_7TeV      param      0        0.002                      \n";
        card += "sig_sigma_err_4e_7TeV     param      0        0.2                        \n";
        card += "CMS_hzz2mu2e_Zjets        lnN        -       -       -        -       -        -        -        ZX_SYST\n";
        }
        else if (!do7TeV && isLowMass) {
        card += "sig_mean_err_4mu_8TeV     param      0        0.001                     \n";
        card += "sig_sigma_err_4mu_8TeV    param      0        0.2                       \n";
        card += "sig_mean_err_4e_8TeV      param      0        0.002                     \n";
        card += "sig_sigma_err_4e_8TeV     param      0        0.2                       \n";
        card += "CMS_hzz2mu2e_Zjets        lnN        -       -       -        -       -        -        -        ZX_SYST\n";
        }
    }
    if (!do1D) {
        card += "CMS_zz4l_bkgMELA          param      0       1       [-3,3]             \n";
    }
    return card;
}

std::string createCardTemplateForSignalHypothesis(bool do7TeV, int channel, bool do1D, std::string workspacefilename, bool isLowMass=true) {

    std::string card = "";
        card += "imax 1\n";
        card += "jmax *\n";
        card += "kmax *\n";
        card += "------------\n";
        card += std::string("shapes * * ") + workspacefilename + " w:$PROCESS\n";
        card += "------------\n";
        card += "bin         BIN\n";
        card += "observation OBS\n";
        card += "------------\n";
        card += "bin     BIN   BIN   BIN   BIN   BIN   BIN     BIN     BIN     BIN     BIN      BIN             BIN             BIN\n";
        card += "process ggH   qqH   WH    ZH    ttH   ggH_ALT qqH_ALT WH_ALT  ZH_ALT  ttH_ALT  bkg_qqzz        bkg_ggzz        bkg_zjets\n";
        card += "process -10   -9    -8    -7    -6     -5      -4      -3      -2      -1       1               2               3\n";
        card += "rate    1     1      1     1     1     1       1       1       1       1        BKG_QQZZ_YIELD  BKG_GGZZ_YIELD  BKG_ZJETS_YIELD\n";
        card += "------------\n";
        card += "## mass window [LOW_EDGE, HIGH_EDGE]\n";

    if (do7TeV) {        
        card += "lumi_7TeV                 lnN        1.022   1.022   1.022    1.022   1.022    1.022   1.022   1.022    1.022   1.022     1.022    1.022    -\n";
    }
    else {        
        card += "lumi_8TeV                 lnN        1.05    1.05    1.05     1.05    1.05     1.05    1.05    1.05     1.05    1.05      1.05     1.05     -\n";
    }
        card += "pdf_gg                    lnN        GGH_PDF -       -        -       -        GGH_PDF -       -        -       -         -        GGZZ_PDF -\n";
        card += "pdf_qqbar                 lnN        -       VBF_PDF WHI_PDF  ZHI_PDF TTH_PDF  -       VBF_PDF WHI_PDF  ZHI_PDF TTH_PDF   QQZZ_PDF -        -\n";
        card += "pdf_hzz4l_accept          lnN        1.02    1.02    1.02     1.02    1.02     1.02    1.02    1.02     1.02    1.02      -        -        -\n";
        card += "QCDscale_ggH              lnN        GGH_QCD -       -        -       -        GGH_QCD -       -        -       -         -        -        -\n";
        card += "QCDscale_qqH              lnN        -       VBF_QCD -        -       -        -       VBF_QCD -        -       -         -        -        -\n";
        card += "QCDscale_WH               lnN        -       -       WHI_QCD  -       -        -       -       WHI_QCD  -       -         -        -        -\n";
        card += "QCDscale_ZH               lnN        -       -       -        ZHI_QCD -        -       -       -        ZHI_QCD -         -        -        -\n";
        card += "QCDscale_ttH              lnN        -       -       -        -       TTH_QCD  -       -       -        -       TTH_QCD   -        -        -\n";
        card += "QCDscale_VV               lnN        -       -       -        -       -        -       -       -        -       -         QQZZ_QCD -        -\n";
        card += "QCDscale_ggVV             lnN        -       -       -        -       -        -       -       -        -       -         -        GGZZ_QCD -\n";
        card += "BRhiggs_ZZ4l              lnN        1.02    1.02    1.02     1.02    1.02     1.02    1.02    1.02     1.02    1.02      -        -        -\n";
    if (!isLowMass) {        
        card += "sig_gamma_err             param      0        0.05\n";
    }
    if (channel == 0) {
        card += "CMS_eff_m                 lnN        1.015   1.015   1.015    1.015   1.015    1.015   1.015   1.015    1.015   1.015     1.015    1.015    -\n";
        card += "CMS_hzz4mu_Zjets          lnN        -       -       -        -       -        -       -       -        -       -         -        -        ZX_SYST\n";
        if (do7TeV) {
        card += "sig_mean_err_4mu_7TeV     param      0        0.001                     \n";
        card += "sig_sigma_err_4mu_7TeV    param      0        0.2                       \n";
        }
        else {
        card += "sig_mean_err_4mu_8TeV     param      0        0.001                     \n";
        card += "sig_sigma_err_4mu_8TeV    param      0        0.2                       \n";
        }
    }
    else if (channel == 1) {
        card += "CMS_eff_e                 lnN        1.01    1.01    1.01     1.01     1.01    1.01    1.01    1.01     1.01     1.01     1.01     1.01     -\n";
        card += "CMS_hzz4e_Zjets           lnN        -       -       -        -        -       -       -       -        -        -        -       -     ZX_SYST\n";
        if (do7TeV) {
        card += "sig_mean_err_4e_7TeV      param      0        0.002                     \n";
        card += "sig_sigma_err_4e_7TeV     param      0        0.2                       \n";
        }
        else {
        card += "sig_mean_err_4e_8TeV      param      0        0.002                     \n";
        card += "sig_sigma_err_4e_8TeV     param      0        0.2                       \n";
        }
    }
    else if (channel == 2) {
        card += "CMS_eff_m                 lnN        1.015   1.015   1.015    1.015   1.015    1.015   1.015   1.015    1.015   1.015    1.015    1.015    -\n";
        card += "CMS_eff_e                 lnN        1.01    1.01    1.01     1.01    1.01     1.01    1.01    1.01     1.01    1.01     1.01     1.01     -\n";
        card += "CMS_hzz2e2mu_Zjets        lnN        -       -       -        -       -        -       -       -        -       -        -        -        ZX_SYST\n";
        if (do7TeV && isLowMass) {
        card += "sig_mean_err_4mu_7TeV     param      0        0.001                      \n";
        card += "sig_sigma_err_4mu_7TeV    param      0        0.2                        \n";
        card += "sig_mean_err_4e_7TeV      param      0        0.002                      \n";
        card += "sig_sigma_err_4e_7TeV     param      0        0.2                        \n";
        }
        else if (!do7TeV && isLowMass) {
        card += "sig_mean_err_4mu_8TeV     param      0        0.001                     \n";
        card += "sig_sigma_err_4mu_8TeV    param      0        0.2                       \n";
        card += "sig_mean_err_4e_8TeV      param      0        0.002                     \n";
        card += "sig_sigma_err_4e_8TeV     param      0        0.2                       \n";
        }
    }
    else if (channel == 3) {
        card += "CMS_eff_m                 lnN        1.015   1.015   1.015    1.015   1.015    1.015   1.015   1.015    1.015   1.015    1.015    1.015    -\n";
        card += "CMS_eff_e                 lnN        1.01    1.01    1.01     1.01    1.01     1.01    1.01    1.01     1.01    1.01     1.01     1.01     -\n";
        card += "CMS_hzz2mu2e_Zjets        lnN        -       -       -        -       -        -       -       -        -       -        -        -        ZX_SYST\n";
        if (do7TeV && isLowMass) {
        card += "sig_mean_err_4mu_7TeV     param      0        0.001                      \n";
        card += "sig_sigma_err_4mu_7TeV    param      0        0.2                        \n";
        card += "sig_mean_err_4e_7TeV      param      0        0.002                      \n";
        card += "sig_sigma_err_4e_7TeV     param      0        0.2                        \n";
        }
        else if (!do7TeV && isLowMass) {
        card += "sig_mean_err_4mu_8TeV     param      0        0.001                     \n";
        card += "sig_sigma_err_4mu_8TeV    param      0        0.2                       \n";
        card += "sig_mean_err_4e_8TeV      param      0        0.002                     \n";
        card += "sig_sigma_err_4e_8TeV     param      0        0.2                       \n";
        }
    }
    if (!do1D) {
        card += "CMS_zz4l_bkgMELA          param      0       1       [-3,3]             \n";
    }
    return card;
}

std::string createCardTemplateNoVH(bool do7TeV, int channel, bool do1D, std::string workspacefilename, bool isLowMass=true) {

    std::string card = "";
        card += "imax 1\n";
        card += "jmax *\n";
        card += "kmax *\n";
        card += "------------\n";
        card += std::string("shapes * * ") + workspacefilename + " w:$PROCESS\n";
        card += "------------\n";
        card += "bin         BIN\n";
        card += "observation OBS\n";
        card += "------------\n";
        card += "bin     BIN           BIN            BIN             BIN             BIN\n";
        card += "process ggH           qqH            bkg_qqzz        bkg_ggzz        bkg_zjets\n";
        card += "process -2            -1             1               2               3\n";
        card += "rate    1             1              BKG_QQZZ_YIELD  BKG_GGZZ_YIELD  BKG_ZJETS_YIELD\n";
        card += "------------\n";
        card += "## mass window [LOW_EDGE, HIGH_EDGE]\n";

    if (do7TeV) {        
        card += "lumi_7TeV                 lnN        1.022   1.022    1.022    1.022    -\n";
    }
    else {        
        card += "lumi_8TeV                 lnN        1.05    1.05     1.05     1.05     -\n";
    }
        card += "pdf_gg                    lnN        GGH_PDF -        -        GGZZ_PDF -\n";
        card += "pdf_qqbar                 lnN         -      VBF_PDF  QQZZ_PDF -        -\n";
        card += "pdf_hzz4l_accept          lnN        1.02    1.02     -        -        -\n";
        card += "QCDscale_ggH              lnN        GGH_QCD -        -        -        -\n";
        card += "QCDscale_qqH              lnN        -       VBF_QCD  -        -        -\n";
        card += "QCDscale_WH               lnN        -       -        -        -        -\n";
        card += "QCDscale_ZH               lnN        -       -        -        -        -\n";
        card += "QCDscale_ttH              lnN        -       -        -        -        -\n";
        card += "QCDscale_VV               lnN        -       -        QQZZ_QCD -        -\n";
        card += "QCDscale_ggVV             lnN        -       -        -        GGZZ_QCD -\n";
        card += "BRhiggs_ZZ4l              lnN        1.02    1.02     -        -        -\n";
    if (!isLowMass) {        
        card += "sig_gamma_err             param      0        0.05\n";
    }
    if (channel == 0) {
        card += "CMS_eff_m                 lnN        1.015   1.015    1.015    1.015    -\n";
        card += "CMS_hzz4mu_Zjets          lnN        -       -        -        -        ZX_SYST\n";
        if (do7TeV) {
        card += "sig_mean_err_4mu_7TeV     param      0        0.001                     \n";
        card += "sig_sigma_err_4mu_7TeV    param      0        0.2                       \n";
        }
        else {
        card += "sig_mean_err_4mu_8TeV     param      0        0.001                     \n";
        card += "sig_sigma_err_4mu_8TeV    param      0        0.2                       \n";
        }
    }
    else if (channel == 1) {
        card += "CMS_eff_e                 lnN        1.01    1.01    1.01     1.01     -\n";
        card += "CMS_hzz4e_Zjets           lnN        -       -       -        -        ZX_SYST\n";
        if (do7TeV) {
        card += "sig_mean_err_4e_7TeV      param      0        0.002                     \n";
        card += "sig_sigma_err_4e_7TeV     param      0        0.2                       \n";
        }
        else {
        card += "sig_mean_err_4e_8TeV      param      0        0.002                     \n";
        card += "sig_sigma_err_4e_8TeV     param      0        0.2                       \n";
        }
    }
    else if (channel == 2) {
        card += "CMS_eff_m                 lnN        1.015   1.015    1.015    1.015    -\n";
        card += "CMS_eff_e                 lnN        1.01    1.01    1.01     1.01     -\n";
        card += "CMS_hzz2e2mu_Zjets        lnN        -       -       -        -        ZX_SYST\n";
        if (do7TeV && isLowMass) {
        card += "sig_mean_err_4mu_7TeV     param      0        0.001                      \n";
        card += "sig_sigma_err_4mu_7TeV    param      0        0.2                        \n";
        card += "sig_mean_err_4e_7TeV      param      0        0.002                      \n";
        card += "sig_sigma_err_4e_7TeV     param      0        0.2                        \n";
        }
        else if (!do7TeV && isLowMass) {
        card += "sig_mean_err_4mu_8TeV     param      0        0.001                     \n";
        card += "sig_sigma_err_4mu_8TeV    param      0        0.2                       \n";
        card += "sig_mean_err_4e_8TeV      param      0        0.002                     \n";
        card += "sig_sigma_err_4e_8TeV     param      0        0.2                       \n";
        }
    }
    else if (channel == 3) {
        card += "CMS_eff_m                 lnN        1.015   1.015    1.015    1.015    -\n";
        card += "CMS_eff_e                 lnN        1.01    1.01    1.01     1.01     -\n";
        card += "CMS_hzz2mu2e_Zjets        lnN        -       -       -        -        ZX_SYST\n";
        if (do7TeV && isLowMass) {
        card += "sig_mean_err_4mu_7TeV     param      0        0.001                      \n";
        card += "sig_sigma_err_4mu_7TeV    param      0        0.2                        \n";
        card += "sig_mean_err_4e_7TeV      param      0        0.002                      \n";
        card += "sig_sigma_err_4e_7TeV     param      0        0.2                        \n";
        }
        else if (!do7TeV && isLowMass) {
        card += "sig_mean_err_4mu_8TeV     param      0        0.001                     \n";
        card += "sig_sigma_err_4mu_8TeV    param      0        0.2                       \n";
        card += "sig_mean_err_4e_8TeV      param      0        0.002                     \n";
        card += "sig_sigma_err_4e_8TeV     param      0        0.2                       \n";
        }
    }
    if (!do1D) {
        card += "CMS_zz4l_bkgMELA          param      0       1       [-3,3]             \n";
    }
    return card;
}

std::string createCardTemplateForSignalHypothesisNoVH(bool do7TeV, int channel, bool do1D, std::string workspacefilename, bool isLowMass=true) {

    std::string card = "";
        card += "imax 1\n";
        card += "jmax *\n";
        card += "kmax *\n";
        card += "------------\n";
        card += std::string("shapes * * ") + workspacefilename + " w:$PROCESS\n";
        card += "------------\n";
        card += "bin         BIN\n";
        card += "observation OBS\n";
        card += "------------\n";
        card += "bin     BIN   BIN   BIN     BIN     BIN             BIN             BIN\n";
        card += "process ggH   qqH   ggH_ALT qqH_ALT bkg_qqzz        bkg_ggzz        bkg_zjets\n";
        card += "process -4   -3     -2      -1      1               2               3\n";
        card += "rate    1     1      1       1       BKG_QQZZ_YIELD  BKG_GGZZ_YIELD  BKG_ZJETS_YIELD\n";
        card += "------------\n";
        card += "## mass window [LOW_EDGE, HIGH_EDGE]\n";

    if (do7TeV) {        
        card += "lumi_7TeV                 lnN        1.022   1.022    1.022   1.022   1.022    1.022    -\n";
    }
    else {        
        card += "lumi_8TeV                 lnN        1.05    1.05     1.05    1.05    1.05     1.05     -\n";
    }
        card += "pdf_gg                    lnN        GGH_PDF -        GGH_PDF -       -        GGZZ_PDF -\n";
        card += "pdf_qqbar                 lnN        -       VBF_PDF  -       VBF_PDF QQZZ_PDF -        -\n";
        card += "pdf_hzz4l_accept          lnN        1.02    1.02     1.02    1.02    -        -        -\n";
        card += "QCDscale_ggH              lnN        GGH_QCD -        GGH_QCD -       -        -        -\n";
        card += "QCDscale_qqH              lnN        -       VBF_QCD  -       VBF_QCD -        -        -\n";
        card += "QCDscale_WH               lnN        -       -        -       -       -        -        -\n";
        card += "QCDscale_ZH               lnN        -       -        -       -       -        -        -\n";
        card += "QCDscale_ttH              lnN        -       -        -       -       -        -        -\n";
        card += "QCDscale_VV               lnN        -       -        -       -       QQZZ_QCD -        -\n";
        card += "QCDscale_ggVV             lnN        -       -        -       -       -        GGZZ_QCD -\n";
        card += "BRhiggs_ZZ4l              lnN        1.02    1.02     1.02    1.02    -        -        -\n";
    if (!isLowMass) {        
        card += "sig_gamma_err             param      0        0.05\n";
    }
    if (channel == 0) {
        card += "CMS_eff_m                 lnN        1.015   1.015    1.015   1.015   1.015    1.015    -\n";
        card += "CMS_hzz4mu_Zjets          lnN        -       -        -       -       -        -        ZX_SYST\n";
        if (do7TeV) {
        card += "sig_mean_err_4mu_7TeV     param      0        0.001                     \n";
        card += "sig_sigma_err_4mu_7TeV    param      0        0.2                       \n";
        }
        else {
        card += "sig_mean_err_4mu_8TeV     param      0        0.001                     \n";
        card += "sig_sigma_err_4mu_8TeV    param      0        0.2                       \n";
        }
    }
    else if (channel == 1) {
        card += "CMS_eff_e                 lnN        1.01    1.01    1.01    1.01     1.01     1.01     -\n";
        card += "CMS_hzz4e_Zjets           lnN        -       -       -       -        -       -     ZX_SYST\n";
        if (do7TeV) {
        card += "sig_mean_err_4e_7TeV      param      0        0.002                     \n";
        card += "sig_sigma_err_4e_7TeV     param      0        0.2                       \n";
        }
        else {
        card += "sig_mean_err_4e_8TeV      param      0        0.002                     \n";
        card += "sig_sigma_err_4e_8TeV     param      0        0.2                       \n";
        }
    }
    else if (channel == 2) {
        card += "CMS_eff_m                 lnN        1.015   1.015    1.015   1.015   1.015    1.015    -\n";
        card += "CMS_eff_e                 lnN        1.01    1.01    1.01    1.01     1.01     1.01     -\n";
        card += "CMS_hzz2e2mu_Zjets        lnN        -       -       -       -        -       -     ZX_SYST\n";
        if (do7TeV && isLowMass) {
        card += "sig_mean_err_4mu_7TeV     param      0        0.001                      \n";
        card += "sig_sigma_err_4mu_7TeV    param      0        0.2                        \n";
        card += "sig_mean_err_4e_7TeV      param      0        0.002                      \n";
        card += "sig_sigma_err_4e_7TeV     param      0        0.2                        \n";
        }
        else if (!do7TeV && isLowMass) {
        card += "sig_mean_err_4mu_8TeV     param      0        0.001                     \n";
        card += "sig_sigma_err_4mu_8TeV    param      0        0.2                       \n";
        card += "sig_mean_err_4e_8TeV      param      0        0.002                     \n";
        card += "sig_sigma_err_4e_8TeV     param      0        0.2                       \n";
        }
    }
    else if (channel == 3) {
        card += "CMS_eff_m                 lnN        1.015   1.015    1.015   1.015   1.015    1.015    -\n";
        card += "CMS_eff_e                 lnN        1.01    1.01    1.01    1.01     1.01     1.01     -\n";
        card += "CMS_hzz2mu2e_Zjets        lnN        -       -       -       -        -       -     ZX_SYST\n";
        if (do7TeV && isLowMass) {
        card += "sig_mean_err_4mu_7TeV     param      0        0.001                      \n";
        card += "sig_sigma_err_4mu_7TeV    param      0        0.2                        \n";
        card += "sig_mean_err_4e_7TeV      param      0        0.002                      \n";
        card += "sig_sigma_err_4e_7TeV     param      0        0.2                        \n";
        }
        else if (!do7TeV && isLowMass) {
        card += "sig_mean_err_4mu_8TeV     param      0        0.001                     \n";
        card += "sig_sigma_err_4mu_8TeV    param      0        0.2                       \n";
        card += "sig_mean_err_4e_8TeV      param      0        0.002                     \n";
        card += "sig_sigma_err_4e_8TeV     param      0        0.2                       \n";
        }
    }
    if (!do1D) {
        card += "CMS_zz4l_bkgMELA          param      0       1       [-3,3]             \n";
    }
    return card;
}

std::string createCardTemplateForSignalHypothesisSingleChannel(bool do7TeV, int channel, bool do1D, std::string workspacefilename, bool isLowMass=true) {

    std::string card = "";
        card += "imax 1\n";
        card += "jmax *\n";
        card += "kmax *\n";
        card += "------------\n";
        card += std::string("shapes * * ") + workspacefilename + " w:$PROCESS\n";
        card += "------------\n";
        card += "bin         BIN\n";
        card += "observation OBS\n";
        card += "------------\n";
        card += "bin     BIN             BIN            BIN             BIN             BIN\n";
        card += "process ggH             ggH_ALT        bkg_qqzz        bkg_ggzz        bkg_zjets\n";
        card += "process -2              -1             1               2               3\n";
        card += "rate    SIG_GGH_YIELD   S2G_GGH_YIELD  BKG_QQZZ_YIELD  BKG_GGZZ_YIELD  BKG_ZJETS_YIELD\n";
        card += "------------\n";
        card += "## mass window [LOW_EDGE, HIGH_EDGE]\n";

    if (do7TeV) {
        card += "lumi_7TeV                 lnN        1.022   1.022   1.022    1.022    -\n";
    }
    else {
        card += "lumi_8TeV                 lnN        1.05    1.05    1.05     1.05     -\n";
    }
        card += "QCDscale_WH               lnN        -       -       -        -        -\n";
        card += "QCDscale_ZH               lnN        -       -       -        -        -\n";
        card += "QCDscale_ttH              lnN        -       -       -        -        -\n";
        card += "QCDscale_VV               lnN        -       -       QQZZ_QCD -        -\n";
        card += "QCDscale_ggVV             lnN        -       -       -        GGZZ_QCD -\n";
        card += "BRhiggs_ZZ4l              lnN        1.02    1.02    -        -        -\n";
    if (!isLowMass) {
        card += "sig_gamma_err             param      0        0.05\n";
    }
    if (channel == 0) {
        card += "CMS_eff_m                 lnN        1.015   1.015   1.015    1.015    -\n";
        card += "CMS_hzz4mu_Zjets          lnN        -       -       -        -        ZX_SYST\n";
        if (do7TeV) {
        card += "sig_mean_err_4mu_7TeV     param      0        0.001                     \n";
        card += "sig_sigma_err_4mu_7TeV    param      0        0.2                       \n";
        }
        else {
        card += "sig_mean_err_4mu_8TeV     param      0        0.001                     \n";
        card += "sig_sigma_err_4mu_8TeV    param      0        0.2                       \n";
        }
    }

    else if (channel == 1) {
        card += "CMS_eff_e                 lnN        1.01    1.01     1.01     1.01     -\n";
        card += "CMS_hzz4e_Zjets           lnN        -       -        -       -     ZX_SYST\n";
        if (do7TeV) {
        card += "sig_mean_err_4e_7TeV      param      0        0.002                     \n";
        card += "sig_sigma_err_4e_7TeV     param      0        0.2                       \n";
        }
        else {
        card += "sig_mean_err_4e_8TeV      param      0        0.002                     \n";
        card += "sig_sigma_err_4e_8TeV     param      0        0.2                       \n";
        }
    }
    else if (channel == 2) {
        card += "CMS_eff_m                 lnN        1.015   1.015   1.015    1.015    -\n";
        card += "CMS_eff_e                 lnN        1.01    1.01    1.01     1.01     -\n";
        card += "CMS_hzz2e2mu_Zjets        lnN        -       -       -        -        ZX_SYST\n";
        if (do7TeV && isLowMass) {
        card += "sig_mean_err_4mu_7TeV     param      0        0.001                      \n";
        card += "sig_sigma_err_4mu_7TeV    param      0        0.2                        \n";
        card += "sig_mean_err_4e_7TeV      param      0        0.002                      \n";
        card += "sig_sigma_err_4e_7TeV     param      0        0.2                        \n";
        }
        else if (!do7TeV && isLowMass) {   
        card += "sig_mean_err_4mu_8TeV     param      0        0.001                     \n";
        card += "sig_sigma_err_4mu_8TeV    param      0        0.2                       \n";
        card += "sig_mean_err_4e_8TeV      param      0        0.002                     \n";
        card += "sig_sigma_err_4e_8TeV     param      0        0.2                       \n";
        }
    }
    else if (channel == 3) {
        card += "CMS_eff_m                 lnN        1.015   1.015   1.015    1.015    -\n";
        card += "CMS_eff_e                 lnN        1.01    1.01    1.01     1.01     -\n";
        card += "CMS_hzz2mu2e_Zjets        lnN        -       -       -        -        ZX_SYST\n";
        if (do7TeV && isLowMass) {
        card += "sig_mean_err_4mu_7TeV     param      0        0.001                      \n";
        card += "sig_sigma_err_4mu_7TeV    param      0        0.2                        \n";
        card += "sig_mean_err_4e_7TeV      param      0        0.002                      \n";
        card += "sig_sigma_err_4e_7TeV     param      0        0.2                        \n";
        }
        else if (!do7TeV && isLowMass) {   
        card += "sig_mean_err_4mu_8TeV     param      0        0.001                     \n";
        card += "sig_sigma_err_4mu_8TeV    param      0        0.2                       \n";
        card += "sig_mean_err_4e_8TeV      param      0        0.002                     \n";
        card += "sig_sigma_err_4e_8TeV     param      0        0.2                       \n";
        }
    }
    if (!do1D) {
        card += "CMS_zz4l_bkgAlternativeMELA       param      0       1       [-3,3]             \n";
        card += "CMS_zz4l_bkgMELA       param      0       1       [-3,3]             \n";
    }
    return card;
}

#endif
