
/*************************************
 * C++ Based class for syst. smearing
 * Original author: A.C. Marini
 * Date: 01/06/2013
 *
 */

#include <cstdio>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>

#include "TMath.h"

#ifndef SYST_H
#define SYST_H

class Bin{
public:
	Bin(){};
	~Bin(){};
	std::pair<float,float> PtRange;	
	std::pair<float,float> EtaRange;	
	std::pair<float,float> RhoRange;	
	std::string tag;
	bool isInside(float pt,float eta,float rho) const
		{
		if( ( pt<=PtRange.second) && (PtRange.first<pt) && 
		    (RhoRange.first<=rho) && (rho<RhoRange.second) && 
		    (EtaRange.first <= eta) && (eta<EtaRange.second) ) return true;
		else { 
			//printf("pt eta rho (%.0f,%.1f,%.1f) outside of (%.0f-%.0f),(%.1f,%.1f),(%.1f,%.1f\n)\n",pt,eta,rho,PtRange.first,PtRange.second,EtaRange.first,EtaRange.second,RhoRange.first,RhoRange.second); //debug
			return false;} 
		}
};

class Parameters
{
public:
	float a;
	float b;
	float lmin;
	float lmax;
}; 

const bool operator<(const Bin&a,const Bin&b)
{
	if( !(a.PtRange.first == b.PtRange.first) ) return (a.PtRange.first < b.PtRange.first);
	if( !(a.PtRange.second == b.PtRange.second) ) return (a.PtRange.second < b.PtRange.second);
	if( !(a.EtaRange.first == b.EtaRange.first) ) return (a.EtaRange.first < b.EtaRange.first);
	if( !(a.EtaRange.second == b.EtaRange.second) ) return (a.EtaRange.second < b.EtaRange.second);
	if( !(a.RhoRange.first == b.RhoRange.first) ) return (a.RhoRange.first < b.RhoRange.first);
	if( !(a.RhoRange.second == b.RhoRange.second) ) return (a.RhoRange.second < b.RhoRange.second);
	return (a.tag < b.tag) ;
	return false; //never here
}

class QGSyst {
public:
	QGSyst();
	~QGSyst();
	float Smear(float pt,float eta, float rho, float x);
	void SetTagger(std::string tag){tagger_=tag;};
	std::string GetTagger(){ return tagger_;};
	void ReadDatabase(std::string fileName);
private:
	std::string tagger_;
	//FILE *database_;
	std::ifstream database_;
	std::map<Bin,Parameters> corrections_;
	float function(float x0, float a ,float b,float min,float max);

};
#endif

#ifndef SYST_C
#define SYST_C
//-------- Constructor -------
QGSyst::QGSyst()
{
	tagger_="QGLHisto";
}
//-------- Destructor -------
QGSyst::~QGSyst()
{
	if(database_.is_open()) database_.close();
}
//----------Read Database ----------
void QGSyst::ReadDatabase(std::string fileName)
{
	database_.open(fileName.c_str(),std::ios::in);
	if( !database_.is_open() ) { std::cerr<<"ERROR: File "<<fileName<<" not open"<<std::endl; return;}
	std::string line;
	while ( database_.good() )
    		{
      		std::getline (database_,line);
      		//cout << line << endl;
		float pt1,pt2,eta1,eta2,rho1,rho2,a,b,lmin,lmax;
		char tag[1023],leadchar;
		sscanf(line.c_str(),"%c",&leadchar);
		if(  (leadchar=='#') || (leadchar=='!')) continue; 
		sscanf(line.c_str(),"%s %f %f %f %f %f %f %f %f %f %f",&tag[0],&pt1,&pt2,&rho1,&rho2,&eta1,&eta2,&a,&b,&lmin,&lmax);
		Bin B;
			B.PtRange=std::pair<float,float>(pt1,pt2);
			B.EtaRange=std::pair<float,float>(eta1,eta2);
			B.RhoRange=std::pair<float,float>(rho1,rho2);
			B.tag=std::string(tag);
		Parameters P;
			P.a=a;P.b=b;P.lmin=lmin;P.lmax=lmax;	
		corrections_[B]=P;
    		}	
	database_.close();
}

float QGSyst::function(float x0, float a ,float b,float min,float max)
{
using namespace TMath;
float x=(x0-min)/(max-min); 
if(x<0)x=0;
if(x>1)x=1;

float x1= (TanH( a* ATanH(2*x-1)+b )/2+.5 ) ;

return x1*(max-min)+min;
 
}

float QGSyst::Smear(float pt,float eta, float rho, float x){
	for(std::map<Bin,Parameters>::iterator it= corrections_.begin();it!=corrections_.end();it++)
		{
		if(it->first.tag.find(tagger_) == std::string::npos)continue;
		if(it->first.isInside(pt,fabs(eta),rho) )
			return function(x,it->second.a,it->second.b,it->second.lmin,it->second.lmax);
		else continue;
		}
	return -99;
}

#endif
