#ifndef PU_H
#define PU_H

#include <vector>
#include <TH1F.h>


TH1F* hist_puweights_2011;
TH1F* hist_puweights_2012;
TH1F* hist_puweights_2012_52;

std::vector<float> puweights_2011;
std::vector<float> puweights_2012;
std::vector<float> puweights_2012_52;

void initpuweights() {
    puweights_2011.push_back(0.212929);
    puweights_2011.push_back(0.0208114);
    puweights_2011.push_back(0.0584048);
    puweights_2011.push_back(0.538898);
    puweights_2011.push_back(1.357);
    puweights_2011.push_back(1.49913);
    puweights_2011.push_back(1.42247);
    puweights_2011.push_back(1.35904);
    puweights_2011.push_back(1.29946);
    puweights_2011.push_back(1.27925);
    puweights_2011.push_back(1.37845);
    puweights_2011.push_back(1.71246);
    puweights_2011.push_back(1.5291);
    puweights_2011.push_back(1.35234);
    puweights_2011.push_back(1.22215);
    puweights_2011.push_back(1.0155);
    puweights_2011.push_back(1.01137);
    puweights_2011.push_back(0.395465);
    puweights_2011.push_back(0.230984);
    puweights_2011.push_back(0.109883);
    puweights_2011.push_back(0.0433739);
    puweights_2011.push_back(0.0111497);
    puweights_2011.push_back(0.00408801);
    puweights_2011.push_back(0.00115678);
    puweights_2011.push_back(0.000365505);
    puweights_2011.push_back(0.000112391);
    puweights_2011.push_back(3.83894e-05);
    puweights_2011.push_back(1.60651e-05);
    puweights_2011.push_back(4.81412e-06);
    puweights_2011.push_back(1.39717e-06);
    puweights_2011.push_back(1.92368e-06);
    puweights_2011.push_back(4.10748e-06);
    puweights_2011.push_back(2.33157e-05);
    puweights_2011.push_back(4.0181e-05);
    puweights_2011.push_back(4.87786e-05);
    puweights_2011.push_back(0.00194128);
    puweights_2011.push_back(8.97414e-05);
    puweights_2011.push_back(1);
    puweights_2011.push_back(1);
    puweights_2011.push_back(0.000162709);
    puweights_2011.push_back(1);
    puweights_2011.push_back(1);
    puweights_2011.push_back(1);
    puweights_2011.push_back(1);
    puweights_2011.push_back(1);
    puweights_2011.push_back(1);
    puweights_2011.push_back(1);
    puweights_2011.push_back(1);
    puweights_2011.push_back(1);
    puweights_2011.push_back(1);

    /*
    puweights_2012.push_back(0.409409);
    puweights_2012.push_back(0.527276);
    puweights_2012.push_back(0.39328);
    puweights_2012.push_back(0.507892);
    puweights_2012.push_back(0.48029);
    puweights_2012.push_back(0.787701);
    puweights_2012.push_back(0.632356);
    puweights_2012.push_back(0.618033);
    puweights_2012.push_back(0.806089);
    puweights_2012.push_back(1.14018);
    puweights_2012.push_back(1.5788);
    puweights_2012.push_back(1.93507);
    puweights_2012.push_back(1.957);
    puweights_2012.push_back(1.73004);
    puweights_2012.push_back(1.46737);
    puweights_2012.push_back(1.28278);
    puweights_2012.push_back(1.18189);
    puweights_2012.push_back(1.13388);
    puweights_2012.push_back(1.12578);
    puweights_2012.push_back(1.14415);
    puweights_2012.push_back(1.16048);
    puweights_2012.push_back(1.1618);
    puweights_2012.push_back(1.15318);
    puweights_2012.push_back(1.13405);
    puweights_2012.push_back(1.09239);
    puweights_2012.push_back(1.01915);
    puweights_2012.push_back(0.914837);
    puweights_2012.push_back(0.786744);
    puweights_2012.push_back(0.644879);
    puweights_2012.push_back(0.502039);
    puweights_2012.push_back(0.371688);
    puweights_2012.push_back(0.263586);
    puweights_2012.push_back(0.18067);
    puweights_2012.push_back(0.120472);
    puweights_2012.push_back(0.0780184);
    puweights_2012.push_back(0.0486113);
    puweights_2012.push_back(0.0289039);
    puweights_2012.push_back(0.0163367);
    puweights_2012.push_back(0.00879674);
    puweights_2012.push_back(0.00456046);
    puweights_2012.push_back(0.0023098);
    puweights_2012.push_back(0.00115977);
    puweights_2012.push_back(0.000583207);
    puweights_2012.push_back(0.000294815);
    puweights_2012.push_back(0.000149865);
    puweights_2012.push_back(7.62892e-05);
    puweights_2012.push_back(3.87537e-05);
    puweights_2012.push_back(1.96105e-05);
    puweights_2012.push_back(9.87744e-06);
    puweights_2012.push_back(4.95418e-06);
    puweights_2012.push_back(2.47913e-06);
    puweights_2012.push_back(1.23919e-06);
    puweights_2012.push_back(6.19751e-07);
    puweights_2012.push_back(3.10125e-07);
    puweights_2012.push_back(1.54934e-07);
    puweights_2012.push_back(7.71425e-08);
    puweights_2012.push_back(3.8182e-08);
    puweights_2012.push_back(1.87455e-08);
    puweights_2012.push_back(9.10765e-09);
    puweights_2012.push_back(9.19802e-09);
    */    

    puweights_2012.push_back(0.246449);
    puweights_2012.push_back(0.319829);
    puweights_2012.push_back(0.332274);
    puweights_2012.push_back(0.378928);
    puweights_2012.push_back(0.324225);
    puweights_2012.push_back(0.571133);
    puweights_2012.push_back(0.445285);
    puweights_2012.push_back(0.431313);
    puweights_2012.push_back(0.594374);
    puweights_2012.push_back(0.911241);
    puweights_2012.push_back(1.3186);
    puweights_2012.push_back(1.68813);
    puweights_2012.push_back(1.75232);
    puweights_2012.push_back(1.56782);
    puweights_2012.push_back(1.33493);
    puweights_2012.push_back(1.164);
    puweights_2012.push_back(1.07775);
    puweights_2012.push_back(1.0522);
    puweights_2012.push_back(1.06937);
    puweights_2012.push_back(1.11167);
    puweights_2012.push_back(1.15048);
    puweights_2012.push_back(1.17265);
    puweights_2012.push_back(1.18182);
    puweights_2012.push_back(1.17864);
    puweights_2012.push_back(1.15641);
    puweights_2012.push_back(1.11001);
    puweights_2012.push_back(1.03986);
    puweights_2012.push_back(0.949639);
    puweights_2012.push_back(0.842786);
    puweights_2012.push_back(0.724516);
    puweights_2012.push_back(0.602474);
    puweights_2012.push_back(0.484649);
    puweights_2012.push_back(0.37676);
    puweights_2012.push_back(0.2826);
    puweights_2012.push_back(0.203926);
    puweights_2012.push_back(0.141043);
    puweights_2012.push_back(0.0935341);
    puweights_2012.push_back(0.059644);
    puweights_2012.push_back(0.0367331);
    puweights_2012.push_back(0.0220279);
    puweights_2012.push_back(0.0129938);
    puweights_2012.push_back(0.00763653);
    puweights_2012.push_back(0.00454065);
    puweights_2012.push_back(0.00277822);
    puweights_2012.push_back(0.00178079);
    puweights_2012.push_back(0.00120881);
    puweights_2012.push_back(0.000871767);
    puweights_2012.push_back(0.000664965);
    puweights_2012.push_back(0.000531206);
    puweights_2012.push_back(0.000439634);
    puweights_2012.push_back(0.000373677);
    puweights_2012.push_back(0.000323738);
    puweights_2012.push_back(0.000284508);
    puweights_2012.push_back(0.000252544);
    puweights_2012.push_back(0.000225388);
    puweights_2012.push_back(0.000201689);
    puweights_2012.push_back(0.000180452);
    puweights_2012.push_back(0.000161084);
    puweights_2012.push_back(0.000143136);
    puweights_2012.push_back(0.000265881);


    puweights_2012_52.push_back(0.0447136);     
    puweights_2012_52.push_back(0.11785);       
    puweights_2012_52.push_back(0.23825);
    puweights_2012_52.push_back(1.08447);
    puweights_2012_52.push_back(0.102575);
    puweights_2012_52.push_back(0.454605);
    puweights_2012_52.push_back(1.79761);
    puweights_2012_52.push_back(4.00271);
    puweights_2012_52.push_back(6.83281);
    puweights_2012_52.push_back(9.83701);
    puweights_2012_52.push_back(10.7966);
    puweights_2012_52.push_back(12.2356);
    puweights_2012_52.push_back(10.0247);
    puweights_2012_52.push_back(8.49395);
    puweights_2012_52.push_back(7.1125);
    puweights_2012_52.push_back(5.69527);
    puweights_2012_52.push_back(4.31256);
    puweights_2012_52.push_back(3.19305);
    puweights_2012_52.push_back(2.42035);
    puweights_2012_52.push_back(1.91666);
    puweights_2012_52.push_back(1.58485);
    puweights_2012_52.push_back(1.36297);
    puweights_2012_52.push_back(1.21166);
    puweights_2012_52.push_back(1.09466);
    puweights_2012_52.push_back(0.978941);
    puweights_2012_52.push_back(0.84653);
    puweights_2012_52.push_back(0.699235);
    puweights_2012_52.push_back(0.548996);
    puweights_2012_52.push_back(0.408673);
    puweights_2012_52.push_back(0.288194);
    puweights_2012_52.push_back(0.193367);
    puweights_2012_52.push_back(0.124653);
    puweights_2012_52.push_back(0.0781124);
    puweights_2012_52.push_back(0.0479268);
    puweights_2012_52.push_back(0.0287763);
    puweights_2012_52.push_back(0.0167744);
    puweights_2012_52.push_back(0.00941834);
    puweights_2012_52.push_back(0.00507877);
    puweights_2012_52.push_back(0.00264364);
    puweights_2012_52.push_back(0.00134612);
    puweights_2012_52.push_back(0.000682678);
    puweights_2012_52.push_back(0.000351412);
    puweights_2012_52.push_back(0.0001864);
    puweights_2012_52.push_back(0.00010259);
    puweights_2012_52.push_back(5.87818e-05);
    puweights_2012_52.push_back(3.5033e-05);
    puweights_2012_52.push_back(2.17116e-05);
    puweights_2012_52.push_back(1.39777e-05);
    puweights_2012_52.push_back(9.36123e-06);
    puweights_2012_52.push_back(6.53328e-06);
    puweights_2012_52.push_back(4.76598e-06);
    puweights_2012_52.push_back(3.64139e-06);
    puweights_2012_52.push_back(2.92018e-06);
    puweights_2012_52.push_back(2.4602e-06);
    puweights_2012_52.push_back(2.17291e-06);
    puweights_2012_52.push_back(2.01107e-06);
    puweights_2012_52.push_back(1.94392e-06);
    puweights_2012_52.push_back(1.9598e-06);
    puweights_2012_52.push_back(2.0583e-06);
    puweights_2012_52.push_back(2.24895e-06);


    hist_puweights_2012 = new TH1F("hist_puweights_2012","",60,0.,60.);

    for(int k=0;k<60;k++){
        hist_puweights_2012->SetBinContent(k+1,puweights_2012[k]);
    }

    hist_puweights_2012_52 = new TH1F("hist_puweights_2012_52","",60,0.,60.);

    for(int k=0;k<60;k++){
        hist_puweights_2012_52->SetBinContent(k+1,puweights_2012_52[k]);
    }

    hist_puweights_2011 = new TH1F("hist_puweights_2011","",50,0.,50.);

    for(int k=0;k<50;k++){
        hist_puweights_2011->SetBinContent(k+1,puweights_2011[k]);
    }
}

float getPUWeight2011(float numsim) {
    return hist_puweights_2011->GetBinContent(hist_puweights_2011->FindBin(numsim));
}

float getPUWeight2012(float numsim, int mode=1) {
    if (mode == 1) return hist_puweights_2012->GetBinContent(hist_puweights_2012->FindBin(numsim));
    else           return hist_puweights_2012_52->GetBinContent(hist_puweights_2012_52->FindBin(numsim));
}

#endif  
