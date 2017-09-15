#ifndef EgmCorrectorParameters_h
#define EgmCorrectorParameters_h

#include "CondFormats/Serialization/interface/Serializable.h"
#include "CondFormats/JetMETObjects/interface/Utilities.h"

#include <string>
#include <vector>
#include <tuple>
#include <algorithm>
#include <functional>
#include <iostream>
#include <memory>
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

class EgmCorrectorParameters 
{
  //---------------- EgmCorrectorParameters class ----------------
  //-- Encapsulates all the information of the parametrization ---
 public:
  //---------------- Definitions class ---------------------------
  //-- Global iformation about the parametrization is kept here --
    class Definitions 
    {
    public:
      //-------- Constructors -------------- 
    Definitions() {}
      Definitions(const std::vector<std::string>& fBinVar, const std::vector<std::string>& fParVar); 
      Definitions(const std::string& fLine); 
      //-------- Member functions ----------
      unsigned nBinVar()                  const {return mBinVar.size(); }
      unsigned nParVar()                  const {return mParVar.size(); }
      std::vector<std::string> parVar()   const {return mParVar;        }
      std::vector<std::string> binVar()   const {return mBinVar;        } 
      std::string parVar(unsigned fIndex) const {return mParVar[fIndex];}
      std::string binVar(unsigned fIndex) const {return mBinVar[fIndex];} 
    private:
      //-------- Member variables ----------
      std::vector<std::string> mParVar;
      std::vector<std::string> mBinVar;
    
      COND_SERIALIZABLE;
    };
    //---------------- Record class --------------------------------
    //-- Each Record holds the properties of a bin ----------------- 
    class Record 
    {
    public:
      //-------- Constructors --------------
    Record() : mNvar(0),mMin(0),mMax(0) {}
    Record(unsigned fNvar, const std::vector<float>& fXMin, const std::vector<float>& fXMax, const std::vector<float>& fParameters) : mNvar(fNvar),mMin(fXMin),mMax(fXMax),mParameters(fParameters) {}
      Record(const std::string& fLine, unsigned fNvar, unsigned fNpar);
      //-------- Member functions ----------
      unsigned nVar()                     const {return mNvar;                      }
      float xMin(unsigned fVar)           const {return mMin[fVar];                 }
      float xMax(unsigned fVar)           const {return mMax[fVar];                 }
      float xMiddle(unsigned fVar)        const {return 0.5*(xMin(fVar)+xMax(fVar));}
      float parameter(unsigned fIndex)    const {return mParameters[fIndex];        }
      std::vector<float> parameters()     const {return mParameters;                }
      unsigned nParameters()              const {return mParameters.size();         }
      bool operator< (const Record& other) const
      {
        if (xMin(0) < other.xMin(0)) return true;
        if (xMin(0) > other.xMin(0)) return false;
        if (xMin(1) < other.xMin(1)) return true;
        if (xMin(1) > other.xMin(1)) return false;
        return (xMin(2) < other.xMin(2));
      }
    private:
      //-------- Member variables ----------
      unsigned           mNvar;
      std::vector<float> mMin;
      std::vector<float> mMax;
      std::vector<float> mParameters;
    
      COND_SERIALIZABLE;
    };
     
    //-------- Constructors --------------
    EgmCorrectorParameters() { valid_ = false;}
    EgmCorrectorParameters(const std::string& fFile);
    EgmCorrectorParameters(const EgmCorrectorParameters::Definitions& fDefinitions,
                           const std::vector<EgmCorrectorParameters::Record>& fRecords) 
      : mDefinitions(fDefinitions),mRecords(fRecords) { valid_ = true;}
    //-------- Member functions ----------
    const Record& record(unsigned fBin)                                                       const {return mRecords[fBin]; }
    const Definitions& definitions()                                                          const {return mDefinitions;   }
    unsigned size()                                                                           const {return mRecords.size();}
    unsigned size(unsigned fVar)                                                              const;
    int binIndex(const std::vector<float>& fX)                                                const;
    int binIndexN(const std::vector<float>& fX)                                               const;
    int neighbourBin(unsigned fIndex, unsigned fVar, bool fNext)                              const;
    std::vector<float> binCenters(unsigned fVar)                                              const;
    void printScreen()                                                                        const;
    void printFile(const std::string& fFileName)                                              const;
    bool isValid() const { return valid_; }

    static const int                                                           MAX_SIZE_DIMENSIONALITY = 3 COND_TRANSIENT;

 private:
    //-------- Member variables ----------
    EgmCorrectorParameters::Definitions                                        mDefinitions;
    std::vector<EgmCorrectorParameters::Record>                                mRecords;
    bool                                                                       valid_; /// is this a valid set?

    COND_SERIALIZABLE;
};
std::ostream& operator<<(std::ostream& out, const EgmCorrectorParameters::Record& fBin);

#endif
