void PUcalcMC(std::string nameFile){
  
 TFile* f = new TFile(nameFile.c_str());
 TH1F* h = (TH1F*) f->Get("pileup");
 for (int ibin = 0; ibin < (h->GetNbinsX()-1); ibin++) {
  std::cout << h->GetBinContent(ibin+1)/(h->GetSumOfWeights() - h->GetBinContent(h->GetNbinsX())) << " , ";
 }
 std::cout << std::endl;


 TFile* fNew = new TFile("newPU.root","RECREATE");
 TH1F* h2 = new TH1F("pileup","pileup",h->GetNbinsX(),0,h->GetNbinsX());
 for (int ibin = 0; ibin < (h->GetNbinsX()-1); ibin++) {
   h2->SetBinContent(ibin+1, h->GetBinContent(ibin+1)/(h->GetSumOfWeights() - h->GetBinContent(h->GetNbinsX())) );
 }

 h2->Write();

 std::cout << "end" << std::endl;
 std::cout << std::endl;

}
