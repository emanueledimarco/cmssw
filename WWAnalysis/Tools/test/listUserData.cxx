// vvv----- Needed in 5.2.X to avoid choking on C++11 data
#define __GCCXML__
// ^^^----- Needed in 5.2.X to avoid choking on C++11 data

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"

#if !defined(__CINT__) && !defined(__MAKECINT__)
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include <TString.h>
#include <TFile.h>
#endif

void listUserDataElectron(const char *inputTag, int nEvents=1) { 
    fwlite::Event ev(gFile);
    fwlite::Handle<std::vector<pat::Electron> > handle;
    int nfail = 0, ngood = 0;
    for (ev.toBegin(); !ev.atEnd(); ++ev) {
        handle.getByLabel(ev, inputTag);
        if (handle.failedToGet()) {
            std::cerr  << "Failed to get pat::Electrons " << inputTag << " (at try " << nfail << ")" << std::endl;
            if (nfail++ > 10) break;
        }
        if (!handle->empty()) {
           std::cout << "\n\n  ==============  UserData for Electrons  ==============  " << std::endl;
           const pat::Electron &ele = handle->at(0);
           const std::vector<std::string> & userFloatNames = ele.userFloatNames();
           for (int i = 0, n = userFloatNames.size(); i < n; ++i) {
                std::cout << "userFloat " << userFloatNames[i] << ", value =  " << ele.userFloat(userFloatNames[i]) << std::endl;
           }
           const std::vector<std::string> & userIntNames = ele.userIntNames();
           for (int i = 0, n = userIntNames.size(); i < n; ++i) {
                std::cout << "userInt " << userIntNames[i] << ", value =  " << ele.userInt(userIntNames[i]) << std::endl;
           }
#if !defined(__CINT__) && !defined(__MAKECINT__)
           const std::vector<pat::Electron::IdPair> &  electronIDs = ele.electronIDs();
           for (int i = 0, n = electronIDs.size(); i < n; ++i) {
                std::cout << "electronID " << electronIDs[i].first << ", value =  " << electronIDs[i].second << std::endl;
           }
#endif
           ngood++;
           if (ngood >= nEvents) break;
        }
    }
}

void listUserDataMuon(const char *inputTag, int nEvents=1) { 
    fwlite::Event ev(gFile);
    fwlite::Handle<std::vector<pat::Muon> > handle;
    int nfail = 0, ngood = 0;
    for (ev.toBegin(); !ev.atEnd(); ++ev) {
        handle.getByLabel(ev, inputTag);
        if (handle.failedToGet()) {
            std::cerr  << "Failed to get pat::Muons " << inputTag << " (at try " << nfail << ")" << std::endl;
            if (nfail++ > 10) break;
        }
        if (!handle->empty()) {
           std::cout << "\n\n  ==============  UserData for Muons  ==============  " << std::endl;
           const pat::Muon &ele = handle->at(0);
           const std::vector<std::string> & userFloatNames = ele.userFloatNames();
           for (int i = 0, n = userFloatNames.size(); i < n; ++i) {
                std::cout << "userFloat " << userFloatNames[i] << ", value =  " << ele.userFloat(userFloatNames[i]) << std::endl;
           }
           const std::vector<std::string> & userIntNames = ele.userIntNames();
           for (int i = 0, n = userIntNames.size(); i < n; ++i) {
                std::cout << "userInt " << userIntNames[i] << ", value =  " << ele.userInt(userIntNames[i]) << std::endl;
           }
           ngood++;
           if (ngood >= nEvents) break;
        }
    }

}


void listUserData(const char *particle=0, const char *inputTag=0, int nEvents=1) { 
    if (TString(particle) == "Electrons") {
        listUserDataElectron(inputTag == 0 ? "boostedElectrons" : inputTag, nEvents);
    } else if (TString(particle) == "Muons") {
        listUserDataMuon(inputTag == 0 ? "boostedElectrons" : inputTag, nEvents);
    } else {
        listUserDataMuon("boostedMuons", nEvents);
        listUserDataElectron("boostedElectrons", nEvents);
    }
}
