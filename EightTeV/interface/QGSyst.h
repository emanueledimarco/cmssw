
/*************************************
 * C++ Based class for syst. smearing
 * Original author: A.C. Marini
 * Date: 01/06/2013
 *
 */

#ifndef QGSYST_H
#define QGSYST_H

#include <cstdio>
#include <vector>
#include <map>
#include <string>
#include <fstream>



class QGSystBin{
public:
	QGSystBin(){};
	~QGSystBin(){};
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
	bool operator<(const QGSystBin& rhs) const;
};

struct QGSystParameters
{
	float a;
	float b;
	float lmin;
	float lmax;
}; 


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
	std::map<QGSystBin,QGSystParameters> corrections_;
	float function(float x0, float a ,float b,float min,float max);

};
#endif
