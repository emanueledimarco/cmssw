#include "../interface/QGSyst.h"
#include <iostream>
#include "TMath.h"

//
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
		QGSystBin B;
			B.PtRange=std::pair<float,float>(pt1,pt2);
			B.EtaRange=std::pair<float,float>(eta1,eta2);
			B.RhoRange=std::pair<float,float>(rho1,rho2);
			B.tag=std::string(tag);
		QGSystParameters P;
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
	for(std::map<QGSystBin,QGSystParameters>::iterator it= corrections_.begin();it!=corrections_.end();it++)
		{
		if(it->first.tag.find(tagger_) == std::string::npos)continue;
		if(it->first.isInside(pt,fabs(eta),rho) )
			return function(x,it->second.a,it->second.b,it->second.lmin,it->second.lmax);
		else continue;
		}
	return -99;
}

