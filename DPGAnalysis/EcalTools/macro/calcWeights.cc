#include <iostream>


void calcWeights(int firstSample, int lastSample, bool isEB) {

  

  double shape[10];
  double currweights[10];
  if(isEB) {
    shape[0] = 0.00000;
    shape[1] = 0.00000;
    shape[2] = 0.00000;
    shape[3] = 0.01321;
    shape[4] = 0.76730;
    shape[5] = 1.00003;
    shape[6] = 0.88441;
    shape[7] = 0.67007;
    shape[8] = 0.47158;
    shape[9] = 0.31756;
    
    currweights[0] = -0.382846;
    currweights[1] = -0.382846;
    currweights[2] = -0.382846;
    currweights[3] = 0;
    currweights[4] = 0.18463;
    currweights[5] = 0.424762;
    currweights[6] = 0.349301;
    currweights[7] = 0.176883;
    currweights[8] = 0.0129613;
    currweights[9] = 0;

  } else {
    shape[0] = 0.00000;
    shape[1] = 0.00000;
    shape[2] = 0.00000;
    shape[3] = 0.00000;
    shape[4] = 0.70815;
    shape[5] = 1.00008;
    shape[6] = 0.87285;
    shape[7] = 0.65066;
    shape[8] = 0.44798;
    shape[9] = 0.29397;


    currweights[0] = -0.382067;
    currweights[1] = -0.382067;
    currweights[2] = -0.382067;
    currweights[3] = 0;
    currweights[4] = 0.204788;
    currweights[5] = 0.415302;
    currweights[6] = 0.339506;
    currweights[7] = 0.172017;
    currweights[8] = 0.0145885;
    currweights[9] = 0;
  }

  double sampleUsed[10];
  sampleUsed[0]=1;
  sampleUsed[1]=1;
  sampleUsed[2]=1;
  sampleUsed[3]=0;
  sampleUsed[4]=1;
  sampleUsed[5]=1;
  sampleUsed[6]=1;
  sampleUsed[7]=1;
  sampleUsed[8]=1;
  sampleUsed[9]=0;

  double sumwxs=0;

  double norm = 0.;
  for(int iSample=std::max(0,firstSample); iSample<std::min(lastSample,10); ++iSample) {
    norm += pow(shape[iSample],2)*sampleUsed[iSample];
  }

  cout << "sumwxs = " << sumwxs << endl;

  for(int iSample=std::max(0,firstSample); iSample<std::min(lastSample,10); ++iSample) {
    std::cout << "w[" << iSample << "] = " << shape[iSample]/norm * sampleUsed[iSample] << std::endl;
  }


}
