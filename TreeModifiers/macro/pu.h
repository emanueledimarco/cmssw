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
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(0.222451);
    puweights_2012.push_back(0.0658851);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(0.150902);
    puweights_2012.push_back(0.202205);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1.50116);
    puweights_2012.push_back(2.79375);
    puweights_2012.push_back(0.198341);
    puweights_2012.push_back(0.246893);
    puweights_2012.push_back(0.28116);
    puweights_2012.push_back(0.449377);
    puweights_2012.push_back(0.553276);
    puweights_2012.push_back(1.48919);
    puweights_2012.push_back(2.15249);
    puweights_2012.push_back(3.62415);
    puweights_2012.push_back(4.33041);
    puweights_2012.push_back(3.57192);
    puweights_2012.push_back(4.99603);
    puweights_2012.push_back(7.79303);
    puweights_2012.push_back(8.04276);
    puweights_2012.push_back(8.05557);
    puweights_2012.push_back(12.9364);
    puweights_2012.push_back(9.9036);
    puweights_2012.push_back(14.6975);
    puweights_2012.push_back(13.3387);
    puweights_2012.push_back(10.9734);
    puweights_2012.push_back(12.6077);
    puweights_2012.push_back(11.5617);
    puweights_2012.push_back(10.8107);
    puweights_2012.push_back(14.5043);
    puweights_2012.push_back(17.8497);
    puweights_2012.push_back(11.8817);
    puweights_2012.push_back(9.6805);
    puweights_2012.push_back(12.2255);
    puweights_2012.push_back(10.1117);
    puweights_2012.push_back(10.2482);
    puweights_2012.push_back(11.5398);
    puweights_2012.push_back(9.35737);
    puweights_2012.push_back(9.90259);
    puweights_2012.push_back(9.19216);
    puweights_2012.push_back(7.57377);
    puweights_2012.push_back(7.94847);
    puweights_2012.push_back(7.15578);
    puweights_2012.push_back(5.63016);
    puweights_2012.push_back(5.35972);
    puweights_2012.push_back(5.05791);
    puweights_2012.push_back(3.35313);
    puweights_2012.push_back(3.60582);
    puweights_2012.push_back(3.35256);
    puweights_2012.push_back(2.49496);
    puweights_2012.push_back(2.28219);
    puweights_2012.push_back(2.21227);
    puweights_2012.push_back(1.76362);
    puweights_2012.push_back(1.68533);
    puweights_2012.push_back(1.62149);
    puweights_2012.push_back(1.34263);
    puweights_2012.push_back(1.30646);
    puweights_2012.push_back(1.21918);
    puweights_2012.push_back(1.10347);
    puweights_2012.push_back(1.08544);
    puweights_2012.push_back(1.0251);
    puweights_2012.push_back(0.907123);
    puweights_2012.push_back(0.905997);
    puweights_2012.push_back(0.869217);
    puweights_2012.push_back(0.816708);
    puweights_2012.push_back(0.76043);
    puweights_2012.push_back(0.714367);
    puweights_2012.push_back(0.679723);
    puweights_2012.push_back(0.665294);
    puweights_2012.push_back(0.609956);
    puweights_2012.push_back(0.586386);
    puweights_2012.push_back(0.548999);
    puweights_2012.push_back(0.521088);
    puweights_2012.push_back(0.4929);
    puweights_2012.push_back(0.453545);
    puweights_2012.push_back(0.44546);
    puweights_2012.push_back(0.406266);
    puweights_2012.push_back(0.378486);
    puweights_2012.push_back(0.347898);
    puweights_2012.push_back(0.337097);
    puweights_2012.push_back(0.313674);
    puweights_2012.push_back(0.291392);
    puweights_2012.push_back(0.283346);
    puweights_2012.push_back(0.25272);
    puweights_2012.push_back(0.244178);
    puweights_2012.push_back(0.228673);
    puweights_2012.push_back(0.211327);
    puweights_2012.push_back(0.19084);
    puweights_2012.push_back(0.179408);
    puweights_2012.push_back(0.169234);
    puweights_2012.push_back(0.157131);
    puweights_2012.push_back(0.143818);
    puweights_2012.push_back(0.140968);
    puweights_2012.push_back(0.124021);
    puweights_2012.push_back(0.118273);
    puweights_2012.push_back(0.109751);
    puweights_2012.push_back(0.0977754);
    puweights_2012.push_back(0.0967206);
    puweights_2012.push_back(0.0870401);
    puweights_2012.push_back(0.0826372);
    puweights_2012.push_back(0.0746777);
    puweights_2012.push_back(0.0698592);
    puweights_2012.push_back(0.0656062);
    puweights_2012.push_back(0.0601853);
    puweights_2012.push_back(0.057892);
    puweights_2012.push_back(0.0517871);
    puweights_2012.push_back(0.0512109);
    puweights_2012.push_back(0.0465423);
    puweights_2012.push_back(0.0403982);
    puweights_2012.push_back(0.0443631);
    puweights_2012.push_back(0.0399185);
    puweights_2012.push_back(0.0338933);
    puweights_2012.push_back(0.0354274);
    puweights_2012.push_back(0.0310775);
    puweights_2012.push_back(0.0261122);
    puweights_2012.push_back(0.0280878);
    puweights_2012.push_back(0.0264014);
    puweights_2012.push_back(0.021299);
    puweights_2012.push_back(0.0245197);
    puweights_2012.push_back(0.0221076);
    puweights_2012.push_back(0.0189236);
    puweights_2012.push_back(0.0202148);
    puweights_2012.push_back(0.0177248);
    puweights_2012.push_back(0.0163634);
    puweights_2012.push_back(0.0188307);
    puweights_2012.push_back(0.0144512);
    puweights_2012.push_back(0.0134599);
    puweights_2012.push_back(0.0143315);
    puweights_2012.push_back(0.0130668);
    puweights_2012.push_back(0.0108666);
    puweights_2012.push_back(0.0162516);
    puweights_2012.push_back(0.0126035);
    puweights_2012.push_back(0.0102154);
    puweights_2012.push_back(0.0154442);
    puweights_2012.push_back(0.00959973);
    puweights_2012.push_back(0.0106827);
    puweights_2012.push_back(0.0146624);
    puweights_2012.push_back(0.0155156);
    puweights_2012.push_back(0.00761674);
    puweights_2012.push_back(0.0187999);
    puweights_2012.push_back(0.0135013);
    puweights_2012.push_back(0.0160794);
    puweights_2012.push_back(0.0180586);
    puweights_2012.push_back(0.00905508);
    puweights_2012.push_back(0.00869858);
    puweights_2012.push_back(0.0193968);
    puweights_2012.push_back(0.0209201);
    puweights_2012.push_back(0.0084405);
    puweights_2012.push_back(0.0407657);
    puweights_2012.push_back(0.0109116);
    puweights_2012.push_back(0.0262218);
    puweights_2012.push_back(0.0104767);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(0.00658031);
    puweights_2012.push_back(0.0051814);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    */

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
