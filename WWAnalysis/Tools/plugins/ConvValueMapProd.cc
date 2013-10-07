#include "WWAnalysis/Tools/interface/ConvValueMapProd.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"



ConvValueMapProd::ConvValueMapProd(const edm::ParameterSet& iConfig) :
    gsfLabel_(iConfig.getUntrackedParameter<edm::InputTag>("gsfLabel")),
    tkLabel_(iConfig.getUntrackedParameter<edm::InputTag>("tkLabel")),
    beamSpotLabel_(iConfig.getUntrackedParameter<edm::InputTag>("beamSpotLabel")),
    conversionLabel_(iConfig.getUntrackedParameter<edm::InputTag>("conversionLabel")) {

    produces<edm::ValueMap<float> >("dist").setBranchAlias("dist");
    produces<edm::ValueMap<float> >("dcot").setBranchAlias("dcot");
    produces<edm::ValueMap<float> >("passVtxConvert").setBranchAlias("passVtxConvert");

}

void ConvValueMapProd::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    edm::Handle<reco::GsfElectronCollection> gsfH;
    iEvent.getByLabel(gsfLabel_,gsfH);

    edm::Handle<reco::TrackCollection> tkH;
    iEvent.getByLabel(tkLabel_,tkH);

    edm::ESHandle<MagneticField> magneticField;
    iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
    const  MagneticField *mField = magneticField.product();
    double bfield = mField->inTesla(GlobalPoint(0.,0.,0.)).z();

    edm::Handle<reco::BeamSpot> bsHandle;
    iEvent.getByLabel(beamSpotLabel_, bsHandle);
    const reco::BeamSpot &thebs = *bsHandle.product();

    edm::Handle<reco::ConversionCollection> hConversions;
    iEvent.getByLabel(conversionLabel_, hConversions);

    std::vector<float> distsV;
    std::auto_ptr<edm::ValueMap<float> > distsM(new edm::ValueMap<float> ());
    edm::ValueMap<float>::Filler distsF(*distsM);

    std::vector<float> dcotsV;
    std::auto_ptr<edm::ValueMap<float> > dcotsM(new edm::ValueMap<float> ());
    edm::ValueMap<float>::Filler dcotsF(*dcotsM);

    std::vector<float> passVtxConvertV;
    std::auto_ptr<edm::ValueMap<float> > passVtxConvertM(new edm::ValueMap<float> ());
    edm::ValueMap<float>::Filler passVtxConvertF(*passVtxConvertM);

    ConversionFinder convFinder;
    ConversionInfo convInfo;

    for(size_t i=0;i<gsfH->size();++i) {

        convInfo = convFinder.getConversionInfo(gsfH->at(i), tkH, bfield);
        distsV.push_back(convInfo.dist());
        dcotsV.push_back(convInfo.dcot());
        passVtxConvertV.push_back( !ConversionTools::hasMatchedConversion(gsfH->at(i),hConversions,thebs.position()) );

    }

    distsF.insert(gsfH,distsV.begin(),distsV.end());
    dcotsF.insert(gsfH,dcotsV.begin(),dcotsV.end());
    passVtxConvertF.insert(gsfH,passVtxConvertV.begin(),passVtxConvertV.end());

    distsF.fill();
    dcotsF.fill();
    passVtxConvertF.fill();
    iEvent.put(distsM,"dist");
    iEvent.put(dcotsM,"dcot");
    iEvent.put(passVtxConvertM,"passVtxConvert");


}

ConvValueMapProd::~ConvValueMapProd() { }
void ConvValueMapProd::beginJob() { }
void ConvValueMapProd::endJob() { }
DEFINE_FWK_MODULE(ConvValueMapProd);
