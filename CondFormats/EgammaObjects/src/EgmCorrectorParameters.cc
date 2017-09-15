#include "CondFormats/EgammaObjects/interface/EgmCorrectorParameters.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <algorithm>
#include <cmath>
#include <iterator>

//------------------------------------------------------------------------
//--- EgmCorrectorParameters::Definitions constructor --------------------
//--- takes specific arguments for the member variables ------------------
//------------------------------------------------------------------------
EgmCorrectorParameters::Definitions::Definitions(const std::vector<std::string>& fBinVar, const std::vector<std::string>& fParVar)
{
  for(unsigned i=0;i<fBinVar.size();i++)
    mBinVar.push_back(fBinVar[i]);
  for(unsigned i=0;i<fParVar.size();i++)
    mParVar.push_back(fParVar[i]);
}
//------------------------------------------------------------------------
//--- EgmCorrectorParameters::Definitions constructor --------------------
//--- reads the member variables from a string ---------------------------
//------------------------------------------------------------------------
EgmCorrectorParameters::Definitions::Definitions(const std::string& fLine)
{
  std::vector<std::string> tokens = getTokens(fLine);
  if (!tokens.empty())
    {
      if (tokens.size() < 3)
        {
          std::stringstream sserr;
          sserr<<"(line "<<fLine<<"): less than 4 expected tokens:"<<tokens.size();
          handleError("EgmCorrectorParameters::Definitions",sserr.str());
        }
      unsigned nvar = getUnsigned(tokens[0]);
      unsigned npar = getUnsigned(tokens[nvar+1]);
      for(unsigned i=0;i<nvar;i++)
        mBinVar.push_back(tokens[i+1]);
      for(unsigned i=0;i<npar;i++)
        mParVar.push_back(tokens[nvar+2+i]);
    }
}
//------------------------------------------------------------------------
//--- EgmCorrectorParameters::Record constructor -------------------------
//--- reads the member variables from a string ---------------------------
//------------------------------------------------------------------------
EgmCorrectorParameters::Record::Record(const std::string& fLine,unsigned fNvar,unsigned fNpar) : mMin(0),mMax(0)
{
  mNvar = fNvar;
  // quckly parse the line
  std::vector<std::string> tokens = getTokens(fLine);
  std::cout << "DBG: " << tokens.size() << std::endl;
  if (!tokens.empty())
    {
      if (tokens.size() < 3)
        {
          std::stringstream sserr;
          sserr<<"(line "<<fLine<<"): "<<"three tokens expected, "<<tokens.size()<<" provided.";
          handleError("EgmCorrectorParameters::Record",sserr.str());
        }
      for(unsigned i=0;i<mNvar;i++)
        {
          std::cout << "DBG: token i = " << i << "  " << tokens[i*2] << "  " << tokens[i*2+1] << std::endl;
          mMin.push_back(getFloat(tokens[i*2]));
          mMax.push_back(getFloat(tokens[i*2+1]));
        }
      for (unsigned i = (2*mNvar); i < tokens.size(); ++i)
        mParameters.push_back(getFloat(tokens[i]));
      if (mParameters.size()!=fNpar) {
          std::stringstream sserr;
          sserr<<"the numbers of parameters read in the records field: " << mParameters.size() << " different from the one declared in the header (" 
               << fNpar << "). Will use only the ones in the header.";
          handleError("EgmCorrectorParameters::Record",sserr.str());
      }
    }
}
std::ostream& operator<<(std::ostream& out, const EgmCorrectorParameters::Record& fBin)
{
  for(unsigned j=0;j<fBin.nVar();j++)
    out<<fBin.xMin(j)<<" "<<fBin.xMax(j)<<" ";
  out<<fBin.nParameters()<<" ";
  for(unsigned j=0;j<fBin.nParameters();j++)
    out<<fBin.parameter(j)<<" ";
  out<<std::endl;
  return out;
}
//------------------------------------------------------------------------
//--- EgmCorrectorParameters constructor ---------------------------------
//--- reads the member variables from a string ---------------------------
//------------------------------------------------------------------------
EgmCorrectorParameters::EgmCorrectorParameters(const std::string& fFile)
{
  std::ifstream input(fFile.c_str());
  std::string line;
  std::string currentDefinitions = "";
  while (std::getline(input,line))
    {
      std::string tmp = getDefinitions(line);
      if (!tmp.empty())
        {
          currentDefinitions = tmp;
          continue;
        }
      Definitions definitions(currentDefinitions);
      if (!(definitions.nBinVar()==0))
        mDefinitions = definitions;
      std::cout << "DBG: " << mDefinitions.nBinVar() << " " << mDefinitions.nParVar() << " " << mDefinitions.binVar(0) << std::endl;
      Record record(line,mDefinitions.nBinVar(),mDefinitions.nParVar());
      bool check(true);
      for(unsigned i=0;i<mDefinitions.nBinVar();++i)
        if (record.xMin(i)==0 && record.xMax(i)==0)
          check = false;
      if (record.nParameters() == 0)
        check = false;
      if (check)
        mRecords.push_back(record);
    }
  if (currentDefinitions=="")
    handleError("EgmCorrectorParameters","No definitions found!!!");
  if (mRecords.empty()) mRecords.push_back(Record());
  std::sort(mRecords.begin(), mRecords.end());
  valid_ = true;

}
//------------------------------------------------------------------------
//--- returns the index of the record defined by fX ----------------------
//------------------------------------------------------------------------
int EgmCorrectorParameters::binIndex(const std::vector<float>& fX) const
{
  int result = -1;
  unsigned N = mDefinitions.nBinVar();
  if (N != fX.size())
    {
      std::stringstream sserr;
      sserr<<"# bin variables "<<N<<" doesn't correspont to requested #: "<<fX.size();
      handleError("EgmCorrectorParameters",sserr.str());
    }
  unsigned tmp;
  for (unsigned i = 0; i < size(); ++i)
    {
      tmp = 0;
      for (unsigned j=0;j<N;j++)
        if (fX[j] >= record(i).xMin(j) && fX[j] < record(i).xMax(j))
          tmp+=1;
      if (tmp==N)
        {
          result = i;
          break;
        }
    }
  return result;
}
//------------------------------------------------------------------------
//--- returns the neighbouring bins of fIndex in the direction of fVar ---
//------------------------------------------------------------------------
int EgmCorrectorParameters::neighbourBin(unsigned fIndex, unsigned fVar, bool fNext) const
{
  int result = -1;
  unsigned N = mDefinitions.nBinVar();
  if (fVar >= N)
    {
      std::stringstream sserr;
      sserr<<"# of bin variables "<<N<<" doesn't correspond to requested #: "<<fVar;
      handleError("EgmCorrectorParameters",sserr.str());
    }
  unsigned tmp;
  for (unsigned i = 0; i < size(); ++i)
    {
      tmp = 0;
      for (unsigned j=0;j<fVar;j++)
        if (fabs(record(i).xMin(j)-record(fIndex).xMin(j))<0.0001)
          tmp+=1;
      for (unsigned j=fVar+1;j<N;j++)
        if (fabs(record(i).xMin(j)-record(fIndex).xMin(j))<0.0001)
          tmp+=1;
      if (tmp<N-1)
        continue;
      if (tmp==N-1)
        {
          if (fNext)
            if (fabs(record(i).xMin(fVar)-record(fIndex).xMax(fVar))<0.0001)
              tmp+=1;
          if (!fNext)
            if (fabs(record(i).xMax(fVar)-record(fIndex).xMin(fVar))<0.0001)
              tmp+=1;
        }
      if (tmp==N)
        {
          result = i;
          break;
        }
    }
  return result;
}
//------------------------------------------------------------------------
//--- returns the number of bins in the direction of fVar ----------------
//------------------------------------------------------------------------
unsigned EgmCorrectorParameters::size(unsigned fVar) const
{
  if (fVar >= mDefinitions.nBinVar())
    {
      std::stringstream sserr;
      sserr<<"requested bin variable index "<<fVar<<" is greater than number of variables "<<mDefinitions.nBinVar();
      handleError("EgmCorrectorParameters",sserr.str());
    }
  unsigned result = 0;
  float tmpMin(-9999),tmpMax(-9999);
  for (unsigned i = 0; i < size(); ++i)
    if (record(i).xMin(fVar) > tmpMin && record(i).xMax(fVar) > tmpMax)
      {
        result++;
        tmpMin = record(i).xMin(fVar);
        tmpMax = record(i).xMax(fVar);
      }
  return result;
}
//------------------------------------------------------------------------
//--- returns the vector of bin centers of fVar --------------------------
//------------------------------------------------------------------------
std::vector<float> EgmCorrectorParameters::binCenters(unsigned fVar) const
{
  std::vector<float> result;
  for (unsigned i = 0; i < size(); ++i)
    result.push_back(record(i).xMiddle(fVar));
  return result;
}
//------------------------------------------------------------------------
//--- prints parameters on screen ----------------------------------------
//------------------------------------------------------------------------
void EgmCorrectorParameters::printScreen() const
{
  std::cout<<"--------------------------------------------"<<std::endl;
  std::cout<<"////////  PARAMETERS: //////////////////////"<<std::endl;
  std::cout<<"--------------------------------------------"<<std::endl;
  std::cout<<"Number of binning variables:   "<<definitions().nBinVar()<<std::endl;
  std::cout<<"Names of binning variables:    ";
  for(unsigned i=0;i<definitions().nBinVar();i++)
    std::cout<<definitions().binVar(i)<<" ";
  std::cout<<std::endl;
  std::cout<<"--------------------------------------------"<<std::endl;
  std::cout<<"Number of parameter variables: "<<definitions().nParVar()<<std::endl;
  std::cout<<"Names of parameter variables:  ";
  for(unsigned i=0;i<definitions().nParVar();i++)
    std::cout<<definitions().parVar(i)<<" ";
  std::cout<<std::endl;
  std::cout<<"------- Bin contents -----------------------"<<std::endl;
  for(unsigned i=0;i<size();i++)
    {
      for(unsigned j=0;j<definitions().nBinVar();j++)
        std::cout<<record(i).xMin(j)<<" "<<record(i).xMax(j)<<" ";
      std::cout<<record(i).nParameters()<<" ";
      for(unsigned j=0;j<record(i).nParameters();j++)
        std::cout<<record(i).parameter(j)<<" ";
      std::cout<<std::endl;
    }
}
//------------------------------------------------------------------------
//--- prints parameters on file ----------------------------------------
//------------------------------------------------------------------------
void EgmCorrectorParameters::printFile(const std::string& fFileName) const
{
  std::ofstream txtFile;
  txtFile.open(fFileName.c_str());
  txtFile.setf(std::ios::right);
  txtFile<<"{"<<definitions().nBinVar()<<std::setw(15);
  for(unsigned i=0;i<definitions().nBinVar();i++)
    txtFile<<definitions().binVar(i)<<std::setw(15);
  txtFile<<definitions().nParVar()<<std::setw(15);
  for(unsigned i=0;i<definitions().nParVar();i++)
    txtFile<<definitions().parVar(i)<<std::setw(15);
  for(unsigned i=0;i<size();i++)
    {
      for(unsigned j=0;j<definitions().nBinVar();j++)
        txtFile<<record(i).xMin(j)<<std::setw(15)<<record(i).xMax(j)<<std::setw(15);
      txtFile<<record(i).nParameters()<<std::setw(15);
      for(unsigned j=0;j<record(i).nParameters();j++)
        txtFile<<record(i).parameter(j)<<std::setw(15);
      txtFile<<"\n";
    }
  txtFile.close();
}

#include "FWCore/Utilities/interface/typelookup.h"

TYPELOOKUP_DATA_REG(EgmCorrectorParameters);
