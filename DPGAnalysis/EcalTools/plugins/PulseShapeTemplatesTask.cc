#include "DPGAnalysis/EcalTools/plugins//PulseShapeTemplatesTask.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "Geometry/EcalMapping/interface/EcalElectronicsMapping.h"
#include "Geometry/EcalMapping/interface/EcalMappingRcd.h"

#include <iomanip>
#include <algorithm> 

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDigi/interface/EcalDataFrame.h"
#include "DQM/EcalCommon/interface/EcalDQMCommonUtils.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

PulseShapeTemplatesTask::PulseShapeTemplatesTask(const edm::ParameterSet& ps)
{ 
  ebDigiCollectionToken_ = consumes<EBDigiCollection>(ps.getParameter<edm::InputTag>("EBdigiCollection"));
  eeDigiCollectionToken_ = consumes<EEDigiCollection>(ps.getParameter<edm::InputTag>("EEdigiCollection"));

  pedestalMeanFile_ = ps.getParameter<std::string>("pedestalMeanFile");
  templatesFile_ = ps.getParameter<std::string>("templatesFile");
  doBarrel_ = ps.getParameter<bool>("doBarrel");
  calcPedMean_ = ps.getParameter<bool>("calcPedMean");
  feds_ =  ps.getParameter< std::vector<std::string> >("feds");
  gain_ = ps.getParameter<int32_t>("gain"); 
}

void 
PulseShapeTemplatesTask::beginRun(const edm::Run & iRun, const edm::EventSetup & iSetup ) { 
}

void 
PulseShapeTemplatesTask::beginJob() {
  file_ = 0;
  if(pedestalMeanFile_.compare("") && calcPedMean_) file_ = TFile::Open(pedestalMeanFile_.c_str(),"recreate");
  if(templatesFile_.compare("") && !calcPedMean_) file_ = TFile::Open(templatesFile_.c_str(),"recreate");  
}

void
PulseShapeTemplatesTask::analyze(const edm::Event& evt, const edm::EventSetup& es)
{

  using namespace edm;

  ESHandle< EcalElectronicsMapping > ecalmapping;
  es.get< EcalMappingRcd >().get(ecalmapping);
  const EcalElectronicsMapping* TheMapping = ecalmapping.product();

  ecaldqm::setElectronicsMap(TheMapping);

  Handle< EBDigiCollection > pEBDigis;
  Handle< EEDigiCollection > pEEDigis;

  const EBDigiCollection* ebDigis =0;
  const EEDigiCollection* eeDigis =0;


  evt.getByToken( ebDigiCollectionToken_, pEBDigis);
  ebDigis = pEBDigis.product(); // get a ptr to the produc
  edm::LogInfo("EcalUncalibRecHitInfo") << "total # ebDigis: " << ebDigis->size() ;

  evt.getByToken( eeDigiCollectionToken_, pEEDigis);
  eeDigis = pEEDigis.product(); // get a ptr to the product
  edm::LogInfo("EcalUncalibRecHitInfo") << "total # eeDigis: " << eeDigis->size() ;

  const EcalDigiCollection *digis = (doBarrel_) ? (EcalDigiCollection*) ebDigis : (EcalDigiCollection*) eeDigis;

  TFile *pedmeanfile = (calcPedMean_) ? 0 : TFile::Open(pedestalMeanFile_.c_str());
  std::map<DetId, TH1D*> pedmeans_h;

  //  std::cout << "Got the collection of digis with " << digis->size() << " crystals" << std::endl;

  for(EcalDigiCollection::const_iterator digiItr(digis->begin()); digiItr != digis->end(); ++digiItr) {

    DetId id(digiItr->id());

    // EcalDataFrame is not a derived class of edm::DataFrame, but can take edm::DataFrame in the constructor
    EcalDataFrame dataFrame(*digiItr);
    
    int gain(0);
    switch(dataFrame.sample(0).gainId()){
    case 1: gain = 12; break;
    case 2: gain = 6; break;
    case 3: gain = 1; break;
    default: continue;
    }
    if(gain != gain_) continue;         // for the time being, calculate only the gain 12 templates

    if(calcPedMean_) { // calculate the mean of the pedestals and store in a file

      if(pedmeans.count(id)==0) {
        std::vector<double> peds;
        peds.resize(EcalDataFrame::MAXSAMPLES);
        for(int iSample(0); iSample < EcalDataFrame::MAXSAMPLES; iSample++)  peds[iSample] = double(dataFrame.sample(iSample).adc());
        pedmeans[id] = peds;
        nevts_[id] = 0;
        //        std::cout << "inserting new ped for DetId = " << id.rawId() << std::endl;
      } else {
        std::vector<double> &peds = pedmeans[id];
        for(int iSample(0); iSample < EcalDataFrame::MAXSAMPLES; iSample++)  peds[iSample] += double(dataFrame.sample(iSample).adc());        
        nevts_[id] ++;
        //        std::cout << "updating ped for DetId = " << id.rawId() << std::endl;
      }
    } else { // calculate the covariance matrix of the pedestals and store in a file 
      
      int iDCC = ecaldqm::dccId(id) - 1;
      std::string dccname = ecaldqm::smName(iDCC);
      std::vector<std::string>::iterator it= std::find(feds_.begin(), feds_.end(), dccname);
      if(it == feds_.end()) continue;

      TH1D *pedmean = 0;
      if(pedmeans_h.count(id)==0) {
        pedmean = (TH1D*)pedmeanfile->Get(Form("pedmean_%d",id.rawId()));
        pedmeans_h[id] = pedmean;
      } else pedmean = pedmeans_h[id];

      double norm = double(dataFrame.sample(5).adc()) - pedmean->GetBinContent(5+1);
      if(norm==0) continue;

      if(templates.count(id)==0) {
        std::vector<double> templ;
        templ.resize(EcalDataFrame::MAXSAMPLES);
        for(int iSample(0); iSample < EcalDataFrame::MAXSAMPLES; iSample++)  templ[iSample] = (double(dataFrame.sample(iSample).adc()) - pedmean->GetBinContent(iSample+1))/norm;
        templates[id] = templ;
        nevts_[id] = 1;
        //        std::cout << "inserting new ped for DetId = " << id.rawId() << std::endl;
      } else {
        std::vector<double> &templ = templates[id];
        for(int iSample(0); iSample < EcalDataFrame::MAXSAMPLES; iSample++)  templ[iSample] += (double(dataFrame.sample(iSample).adc()) - pedmean->GetBinContent(iSample+1))/norm;
        nevts_[id] ++;
        //        std::cout << "updating ped for DetId = " << id.rawId() << std::endl;
      }
      delete pedmean;
    }
    
  } // loop over digis
  if(!calcPedMean_) pedmeanfile->Close();

}

void 
PulseShapeTemplatesTask::endJob() {
  file_->cd();
    
  if(calcPedMean_) {
    TH1D *hped = new TH1D("hped","",EcalDataFrame::MAXSAMPLES,0,EcalDataFrame::MAXSAMPLES);
    for(std::map<DetId, std::vector<double> >::iterator it=pedmeans.begin(); it!=pedmeans.end(); ++it) {
      //std::cout << "Writing histoigram. DetId = " << it->first.rawId() << " got " << nevts_[it->first] << " events" << std::endl;
      TH1D *ph = (TH1D*)hped->Clone(Form("pedmean_%d",it->first.rawId()));
      for(int iSample(0); iSample < EcalDataFrame::MAXSAMPLES; iSample++) 
        ph->SetBinContent(iSample+1,(it->second)[iSample] / double(nevts_[it->first]));   
      ph->Write();
    }
  } else {
    TH1D *htempl = new TH1D("htempl","",EcalDataFrame::MAXSAMPLES,0,EcalDataFrame::MAXSAMPLES);
    TH1D *htemplAverage =  (TH1D*)htempl->Clone("templates_average");
    double nCry(0);
    for(std::map<DetId, std::vector<double> >::iterator it=templates.begin(); it!=templates.end(); ++it) {
      std::cout << "Writing templates histogram. DetId = " << it->first.rawId() << " got " << nevts_[it->first] << " events" << std::endl;
      TH1D *th = (TH1D*)htempl->Clone(Form("template_%d",it->first.rawId()));
      for(int iSample(0); iSample < EcalDataFrame::MAXSAMPLES; iSample++) 
        th->SetBinContent(iSample+1,(it->second)[iSample] / double(nevts_[it->first]));   
      th->Write();
      if(nevts_[it->first]>0) *htemplAverage = (*htemplAverage) + (*th);
      nCry += 1.;
    }
    htemplAverage->Scale(1./nCry);
    std::cout << "Writing average template histogram (averaged over " << nCry << " crystals" << std::endl;
    htemplAverage->Write();
  }
  file_->Close();
}


DEFINE_FWK_MODULE(PulseShapeTemplatesTask);

