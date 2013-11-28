#ifndef WWANALYSIS_ANALYSISSTEP_LATINOSREWIEGHT_H
#define WWANALYSIS_ANALYSISSTEP_LATINOSREWIEGHT_H

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "DataFormats/Provenance/interface/ProcessHistory.h"
#include "DataFormats/Provenance/interface/Provenance.h"
#include "FWCore/Common/interface/EventBase.h"

// namespace edm {
//     class LatinoReWeighting : public LumiReWeighting {
//         public:
//             LatinoReWeighting( std::string generatedFile,
//                     std::string dataFile,
//                     std::string GenHistName,
//                     std::string DataHistName)
//                 :LumiReWeighting(generatedFile,dataFile,GenHistName,DataHistName) {}
//             LatinoReWeighting( std::vector< double > MC_distr, std::vector< double > Lumi_distr);
//             LatinoReWeighting( std::vector< float > MC_distr, std::vector< float > Lumi_distr) 
//                 :LumiReWeighting(MC_distr,Lumi_distr) {}
//             LatinoReWeighting ( ) { } 
// 
//             double fullWeight( const edm::EventBase &e ) {
// 
//                 edm::ReleaseVersion TargetRelease1 = edm::ReleaseVersion("\"CMSSW_4_2_2_patch1\"");
//                 edm::ReleaseVersion TargetRelease2 = edm::ReleaseVersion("\"CMSSW_4_2_2_patch2\"");
//                 edm::ProcessHistory PHist = e.processHistory();
//                 edm::ProcessHistory::const_iterator PHist_iter = PHist.begin();
// 
//                 bool isS3 = false;
//                 for(; PHist_iter<PHist.end() ;++PHist_iter) {
//                     edm::ProcessConfiguration PConf = *(PHist_iter);
//                     edm::ReleaseVersion Release =  PConf.releaseVersion() ;
//                     const std::string Process =  PConf.processName();
// 
//                     if((Process=="HLT") && (Release==TargetRelease1 || Release == TargetRelease2)) {
//                         isS3 = true;
//                         break;
//                     }
//                 }
// 
//                 if(isS3) return weight(e);
//                 else     return weight3BX(e);
// 
//             }
// 
//     };
//     LatinoReWeighting::LatinoReWeighting( std::vector< double > MC_distr, std::vector< double > Lumi_distr) {
// 
//         if( MC_distr.size() != Lumi_distr.size() ){   
// 
//             std::cerr <<"ERROR: LumiReWeighting: input vectors have different sizes. Quitting... \n";
//             return;
// 
//         }
// 
//         Int_t NBins = MC_distr.size();
// 
//         weights_ = boost::shared_ptr<TH1> ( new TH1F("luminumer","luminumer",NBins,-0.5, float(NBins)-0.5) );
//         TH1* den = new TH1F("lumidenom","lumidenom",NBins,-0.5, float(NBins)-0.5) ;
// 
//         for(int ibin = 1; ibin<NBins+1; ++ibin ) {
//             weights_->SetBinContent(ibin, Lumi_distr[ibin-1]);
//             den->SetBinContent(ibin,MC_distr[ibin-1]);
//         }
// 
//         // check integrals, make sure things are normalized
// 
//         float deltaH = weights_->Integral();
//         if(fabs(1.0 - deltaH) > 0.02 ) { //*OOPS*...
//             weights_->Scale( 1.0/ weights_->Integral() );
//         }
//         float deltaMC = den->Integral();
//         if(fabs(1.0 - deltaMC) > 0.02 ) {
//             den->Scale(1.0/ den->Integral());
//         }
// 
//         weights_->Divide( den );  // so now the average weight should be 1.0    
// 
//         std::cout << " Lumi/Pileup Reweighting: Computed Weights per In-Time Nint " << std::endl;
// 
//         for(int ibin = 1; ibin<NBins+1; ++ibin){
//             std::cout << "   " << ibin-1 << " " << weights_->GetBinContent(ibin) << std::endl;
//         }
// 
//         weightOOT_init();
// 
//         FirstWarning_ = true;
//         OldLumiSection_ = -1;
// 
//     }
// 
// }
// 

namespace edm {
    class LatinoReWeighting {
        public:
            LatinoReWeighting(std::vector<float> s3Distro,std::vector<float> s4Distro,std::vector<float> dataDistro) : 
                s3(s3Distro,dataDistro),s4(s4Distro,dataDistro) {}
            LatinoReWeighting(std::vector<double> s3Distro,std::vector<double> s4Distro,std::vector<double> dataDistro) : 
                s3(dToF(s3Distro),dToF(dataDistro)),s4(dToF(s4Distro),dToF(dataDistro)) {}
            static std::vector<float> dToF(std::vector<double> d) {
                std::vector<float> f(d.size(),0);
                for(size_t i=0;i<f.size();++i) f[i] = d[i];
                return f;
            }
            double fullWeight( const edm::EventBase &e ) {

                edm::ReleaseVersion TargetRelease1 = edm::ReleaseVersion("\"CMSSW_4_2_2_patch1\"");
                edm::ReleaseVersion TargetRelease2 = edm::ReleaseVersion("\"CMSSW_4_2_2_patch2\"");
                edm::ProcessHistory PHist = e.processHistory();
                edm::ProcessHistory::const_iterator PHist_iter = PHist.begin();

                bool isS3 = false;
                for(; PHist_iter<PHist.end() ;++PHist_iter) {
                    edm::ProcessConfiguration PConf = *(PHist_iter);
                    edm::ReleaseVersion Release =  PConf.releaseVersion() ;
                    const std::string Process =  PConf.processName();

                    if((Process=="HLT") && (Release==TargetRelease1 || Release == TargetRelease2)) {
                        isS3 = true;
                        break;
                    }
                }

                std::cout << isS3 << std::endl;
                if(isS3) return s3.weight(e);
                else     return s4.weight3BX(e);

            }
            double weight( const edm::EventBase &e ) { return s3.weight(e); }
            double weight3BX( const edm::EventBase &e ) { return s4.weight3BX(e); }
        private: 
            LumiReWeighting s3;
            LumiReWeighting s4;

    };
}

#endif
