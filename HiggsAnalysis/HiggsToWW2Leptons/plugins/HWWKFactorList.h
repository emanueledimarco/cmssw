#ifndef DataFormats_ParticleFlowReco_HWWKfactorList_h
#define DataFormats_ParticleFlowReco_HWWKfactorList_h

#include <iostream>
#include <string>
#include <stdexcept>

#include <TH1.h>



/// \brief Resolution Map (resolution as a function of eta and E)
///
/// Basically just a TH1D with text I/O
/// \author Joanna Weng
/// \todo extrapolation
/// \date January 2006
class HWWKfactorList : public TH1D {

 public:

  /// default constructor
  HWWKfactorList() : TH1D() {}
  
  /// create a map from text file mapfile
  HWWKfactorList(const char* name, const char* mapfile);
  
  /// create an empty map and initialize it 
  HWWKfactorList(const char* name, 
		  unsigned nbinspt, double minpt, double maxpt, double value);

  // create a variable binning map and initialize it
  HWWKfactorList(const char* name, unsigned nbinspt, const double* ptbins, double value);
  
  /// create a map from a 1d histogram
  HWWKfactorList(const TH1D& h) : TH1D(h) {}
 

  /// read text file
  bool ReadMapFile(const char* mapfile);

  bool ReadVarBinMapFile( const char* mapfile );

  /// write text file
  /// is not const because mapFile_ will be updated
  bool WriteMapFile(const char* mapfile);

  /// all the information is here, why not just ask?
  int GetBinFromPt(const double &pt);
  double GetWeightFromPt(const double &pt) { return GetBinContent(GetBinFromPt(pt)); }
  
  const char* GetMapFile() const {return mapFile_.c_str();}

  /// print this map
  friend std::ostream& operator<<(std::ostream& out, const HWWKfactorList& rm);


  // get alternative Kfactor
  inline double GetAlterKfactor() { return alternativeK_;};
  inline double GetAlterNNLOKfactor() { return alternativeNNLOK_;};
  
 private:
  static const unsigned lineSize_;
  std::string           mapFile_;

  double                alternativeK_;
  double                alternativeNNLOK_;
};

#endif


