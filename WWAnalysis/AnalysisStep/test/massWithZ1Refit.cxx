#include <cmath>
#include <cstdio>
#include <TTree.h>
#include <TFile.h>
#include <TMatrixTSym.h>
#include <TVectorT.h>
#include <TLorentzVector.h>
#include <Math/Functor.h>
#include <Minuit2/Minuit2Minimizer.h>

// must be global, temporarily
Float_t mass, massErr, z1mass, z1mll, z1massErr, z1l1massErr, z1l2massErr, l1massErr, l2massErr;;
Float_t l1pdgId, l3pdgId;
Float_t pt[6], eta[6], phi[6]; TLorentzVector p4[6], psum;
Float_t massNew, massErrNew, z1massNew, z1massChi2In, z1massChi2Out;

TLorentzVector p1s, p2s;
double z1Gauss(const double *x) {
    p1s = p4[0]; p2s = p4[1];
    p1s *= (1.0 + x[0]);
    p2s *= (1.0 + x[1]);
    double term1 = ((p1s + p4[1] + p4[2]).M() - z1mass)/z1l1massErr;
    double term2 = ((p2s + p4[0] + p4[2]).M() - z1mass)/z1l2massErr;
    double term3 = ((p1s + p2s + p4[2]).M() - 91.2)/(2.5/2.33);
    return term1*term1 + term2*term2 + term3*term3;
}

double z1NLL(const double *x) {
    p1s = p4[0]; p2s = p4[1];
    p1s *= (1.0 + x[0]);
    p2s *= (1.0 + x[1]);
    double term1 = ((p1s + p4[1] + p4[2]).M() - z1mass)/z1l1massErr;
    double term2 = ((p2s + p4[0] + p4[2]).M() - z1mass)/z1l2massErr;
    double y = ((p1s + p2s + p4[2]).M() - 91.2)/2.5;
    double term3 = 2*log(1 + y*y);
    return term1*term1 + term2*term2 + term3;
}


void massWithZ1Refit(double maxChi2Cut=10.) {
    TTree *in = (TTree*) gFile->Get("zz4lTree/probe_tree");
    in->SetBranchAddress("mass",&mass);
    in->SetBranchAddress("massErr",&massErr);
    in->SetBranchAddress("z1mass",&z1mass);
    in->SetBranchAddress("z1mll",&z1mll);
    in->SetBranchAddress("z1massErr",&z1massErr);
    in->SetBranchAddress("z1l1massErr",&z1l1massErr);
    in->SetBranchAddress("z1l2massErr",&z1l2massErr);
    in->SetBranchAddress("l1massErr",&l1massErr);
    in->SetBranchAddress("l2massErr",&l2massErr);
    in->SetBranchAddress("l1pdgId",&l1pdgId);
    in->SetBranchAddress("l3pdgId",&l3pdgId);
    // load l1,l2,pho1, l3,l4,pho2
    for (int i = 1, j = 0; i <= 4; ++i, j += (j == 1 ? 2 : 1)) {
        in->SetBranchAddress(Form("l%dpt",i),&pt[j]);
        in->SetBranchAddress(Form("l%deta",i),&eta[j]);
        in->SetBranchAddress(Form("l%dphi",i),&phi[j]);
    }
    for (int i = 1, j = 2; i <= 2; ++i, j += 3) {
        in->SetBranchAddress(Form("pho%dpt",i),&pt[j]);
        in->SetBranchAddress(Form("pho%deta",i),&eta[j]);
        in->SetBranchAddress(Form("pho%dphi",i),&phi[j]);
    }

    TString fname = gFile->GetName(); fname.ReplaceAll(".root",".friend.root");
    TFile *frf = new TFile(fname,"RECREATE");
    TTree *out = new TTree("mass_tree","");
    out->Branch("massNew",    &massNew, "F");
    out->Branch("massErrNew", &massErrNew, "F");
    out->Branch("z1massNew", &z1massNew, "F");
    out->Branch("z1massChi2In", &z1massChi2In, "F");
    out->Branch("z1massChi2Out", &z1massChi2Out, "F");

    ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigrad );
    min.SetMaxFunctionCalls(1000000);
    min.SetMaxIterations(100000);
    min.SetTolerance(0.001);

    ROOT::Math::Functor f(&z1Gauss,2); 
    min.SetFunction(f);

    TMatrixTSym<double> cov(2);
    TVectorD jac(2);

    for (unsigned int e = 0, ne = in->GetEntries(); e < ne; ++e) {
        //if (e > 10000) break;
        in->GetEntry(e);
        psum.SetPxPyPzE(0,0,0,0);
        for (int i = 0; i < 6; ++i) {
            double mass = 0;
            if (i % 3 < 2) mass = (abs(i < 3 ? l1pdgId : l3pdgId) == 13) ? 0.105 : 0.0005;
            p4[i].SetPtEtaPhiM(pt[i],eta[i],phi[i],mass);
            if (i >= 2) psum += p4[i];
        }
           
        double zeros[2] =  { 0., 0. };
        z1massChi2In  =  z1Gauss(zeros);
        if (z1massChi2In > maxChi2Cut) {
            massNew = mass; z1massNew = z1mass;
            massErrNew = massErr;
            z1massChi2Out =  z1massChi2In;
            out->Fill();
            continue;
        }
        min.SetLimitedVariable(0,"x",0.,0.01,-8.,8.);
        min.SetLimitedVariable(1,"y",0.,0.01,-8.,8.);
        min.Minimize();

        const double *xs = min.X(); 
        p1s = p4[0]; p2s = p4[1];
        p1s *= (1.0 + xs[0]);
        p2s *= (1.0 + xs[1]);
        //psum = p1s + p2s;

        z1massNew = (p1s+p2s+p4[2]).M();
        massNew = (p1s+p2s+psum).M();

        z1massChi2Out =  z1Gauss(xs);
        // start from the error coming from the other components

        massErrNew = massErr*massErr - l1massErr*l1massErr - l2massErr*l2massErr;

        // then need the rest.
        // now cov(m) = (dm/dx dm/dy) [ cov(x,y) ] * (dm/dx dm/dy)
        cov(0,0) = min.CovMatrix(0,0);
        cov(0,1) = min.CovMatrix(0,1);
        cov(1,1) = min.CovMatrix(1,1);

        p1s = p4[0]; p1s *= (1.0 + xs[0] + 0.001);
        jac(0) = ((p1s+p2s+psum).M() - massNew)/0.001;
        p1s = p4[0]; p1s *= (1.0 + xs[0]);

        p2s = p4[1]; p2s *= (1.0 + xs[1] + 0.001);
        jac(1) = ((p1s+p2s+psum).M() - massNew)/0.001;
        p2s = p4[1]; p2s *= (1.0 + xs[1]);

        massErrNew = sqrt(massErrNew + cov.Similarity(jac));        

        if (e < 50) {
            printf("%s %6.2f +/- %4.2f -> %6.2f (%3.1f s)  zz4l  %6.2f +/- %5.2f -> %6.2f +/- %5.2f (%3.1f s)  chi2 in %5.1f out %5.1f pull ratio %3.1f\n",
                    abs(l1pdgId) == 13 ? "z1mm" : "z1ee",
                    z1mass,z1massErr,z1massNew, fabs(z1massNew-z1mass)/z1massErr,
                    mass,massErr,massNew, massErrNew, fabs(massNew-mass)/massErr,
                    z1massChi2In, z1massChi2Out, (fabs(massNew-126)-fabs(mass-126))/massErr);
        }


        out->Fill();
    }
    
    out->Write();
    frf->Close();
    printf("Friend tree written to %s\n", fname.Data());
}
