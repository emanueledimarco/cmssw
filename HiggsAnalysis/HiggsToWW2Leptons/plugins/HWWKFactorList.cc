#include "HiggsAnalysis/HiggsToWW2Leptons/plugins/HWWKFactorList.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <functional>
#include <stdexcept>

using namespace std;


const unsigned HWWKfactorList::lineSize_ = 10000;


HWWKfactorList::HWWKfactorList(const char* name, unsigned nbinspt,
                               double minpt, double maxpt, double value)
  : TH1D(name, name, nbinspt, minpt, maxpt), alternativeK_(1.), alternativeNNLOK_(1.)
{

  GetXaxis()->SetTitle("#pt");
  GetYaxis()->SetTitle("Kfactor");
  if(value>0)
  {
    for(int i=1; i<=GetNbinsX(); i++)
    {
      SetBinContent(i, value);
    }
  }

}

HWWKfactorList::HWWKfactorList(const char* name, unsigned nbinspt,
                               const double* ptbins,  double value)
  : TH1D(name, name, nbinspt, ptbins), alternativeK_(1.), alternativeNNLOK_(1.)
{

  GetXaxis()->SetTitle("#pt");
  GetYaxis()->SetTitle("Kfactor");
  if(value>0)
  {
    for(int i=1; i<=GetNbinsX(); i++)
    {
      SetBinContent(i, value);
    }
  }

}

HWWKfactorList::HWWKfactorList(const char* name, const char* mapfile)
  : alternativeK_(1.), alternativeNNLOK_(1.)
{

  SetTitle(mapfile);
  GetXaxis()->SetTitle("#eta");
  GetYaxis()->SetTitle("E");
//   if( ! ReadMapFile(mapfile) )
  if( ! ReadVarBinMapFile(mapfile) )
  {
    string err = "HWWKfactorList::HWWKfactorList : cannot read file ";
    err += mapfile;
    throw invalid_argument(err);
  }

//   std::cout << *this << std::endl;

}


bool HWWKfactorList::WriteMapFile(const char* mapfile)
{
  // open the file
  ofstream outf(mapfile);
  if( !outf.good() )
  {
    cout<<"HWWKfactorList::Write : cannot open file "<<mapfile<<endl;
    return false;
  }

  outf<<(*this)<<endl;
  if(!outf.good() )
  {
    cerr<<"HWWKfactorList::Write : corrupted file "<<mapfile<<endl;
    return false;
  }
  else
  {
    mapFile_ = mapfile;
    return true;
  }
}


bool HWWKfactorList::ReadVarBinMapFile( const char* mapfile ) {

  // open the file
  ifstream inf(mapfile);
  if( !inf.good() )
  {
    return false;
  }
  // first data describes the map: histo bin, max, min
  int nbinspt=0;
  double minpt=0;
  double maxpt=0;
  std::string dummy;
  inf>>nbinspt;inf>>dummy;
  inf>>minpt;inf>>dummy;
  inf>>maxpt;inf>>dummy;
  inf>>this->alternativeK_;inf>>dummy;
  inf>>this->alternativeNNLOK_;inf>>dummy;inf>>dummy;

  std::vector<double> binlow, binhigh, binvalue;

  double low, high, value;
  int k(0);
  for( ;; ++k) { 
    inf >> low >> high >> value;
    if ( !inf.good() )
        break;
//     std::cout << low << "   " << high << " " << value << std::endl;
    binlow.push_back(low);
    binhigh.push_back(high);
    binvalue.push_back(value);

  }

//   std::cout << "k " << k << "   nbins" << nbinspt << std::endl; 

  if ( k != nbinspt ) {
    throw runtime_error( "Wrong number of bins read. Corrupted file?" );
  } 

  double* xbins = new double[nbinspt+1];
  for( int i(0); i<nbinspt; ++i ) {
    if ( i != 0  && binhigh[i-1] != binlow[i] ) 
        throw runtime_error( "Bin edge mismatch. Corrupted file?" );

    xbins[i] = binlow[i];
    xbins[i+1] = binhigh[i];

  }

  this->SetBins(nbinspt, xbins);

  for( int i(0); i<nbinspt; ++i ) {
      this->SetBinContent(i+1,binvalue[i]);
  }

  
  delete [] xbins;
  return true;

}




bool HWWKfactorList::ReadMapFile(const char* mapfile)
{
  // open the file
  ifstream inf(mapfile);
  if( !inf.good() )
  {
    return false;
  }
  // first data describes the map: histo bin, max, min
  int nbinspt=0;
  double minpt=0;
  double maxpt=0;
  std::string dummy;
  inf>>nbinspt;inf>>dummy;
  inf>>minpt;inf>>dummy;
  inf>>maxpt;inf>>dummy;
  inf>>this->alternativeK_;inf>>dummy;
  inf>>this->alternativeNNLOK_;inf>>dummy;inf>>dummy;
  SetBins(nbinspt, minpt, maxpt);

  char s[lineSize_];
  // get position in stream
  int pos=inf.tellg(); int j=0;
  // parse map data
  do
  {
    inf.seekg(pos);
    inf.getline(s,lineSize_);
    pos = inf.tellg();
    //    cout<<"LINE"<< s<<endl;
    if(string(s).empty())
    {
      continue; // remove empty lines
    }
    istringstream lin(s);
    double dataw;

    if (lin.good())
    {
      lin>>dataw;   lin>>dataw;
      //  cout<<"LINE "<<dataw <<endl;
      SetBinContent(j, dataw);
      j++;
    }

  }
  while(inf.good());

  if(inf.eof())
  {
    mapFile_ = mapfile;
    return true;
  }
  else return false;
  mapFile_ = mapfile;
  return true;
}




ostream& operator<<(ostream& outf, const HWWKfactorList& rm)
{

  if(!outf.good() ) return outf;

  // first data describes the map
  outf<<rm.GetNbinsX()<< "\t#nbins" << endl;
  outf<<rm.GetXaxis()->GetXmin()<< "\t#lowedge"<< endl;
  outf<<rm.GetXaxis()->GetXmax()<< "\t#highedge" << endl;

  TAxis* ax = rm.GetXaxis();
  for(int i=0; i<=rm.GetNbinsX(); i++){
    outf<< ax->GetBinLowEdge(i) << '\t' << ax->GetBinUpEdge(i) << '\t' << rm.GetBinContent(i)<<"\n";
  }
  outf<<endl;
  return outf;
}

int HWWKfactorList::GetBinFromPt(const double& pt) {

    if (pt < GetBinLowEdge(1) || pt > GetBinLowEdge(GetNbinsX()+1) )  {
        return GetNbinsX()+1;
    } else {
        return FindBin(pt);
    }

}
