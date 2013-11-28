#include "WWAnalysis/AnalysisStep/interface/CreateEventHists.h"

#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "PhysicsTools/UtilAlgos/interface/FWLiteAnalyzerWrapper.h"

typedef fwlite::AnalyzerWrapper<CreateEventHists> WrappedFWLiteAnalyzer;

int main(int argc, char **argv) {
    gSystem->Load( "libFWCoreFWLite" );
    AutoLibraryLoader::enable();

    if ( argc < 2 ) {
        std::cout << "Usage : " << argv[0] << " [parameters.py]" << std::endl;
        return 1;
    }
    if( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") ){
        std::cout << " ERROR: ParametersSet 'plot' is missing in your configuration file" << std::endl; 
        return 2;
    }

    edm::ParameterSet process = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");
    WrappedFWLiteAnalyzer analyzer(process,"eventHists","eventHists");

    analyzer.beginJob();
    analyzer.analyze();
    analyzer.endJob();

    return 0;
}
