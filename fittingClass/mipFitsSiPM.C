#include "TROOT.h"
#include "TObject.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TLegend.h"
#include "TMath.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TLatex.h"
#include <cmath>
#include <cstdio>
#include <vector>
#include "SPEfunc.h"

//This function fits the SiPM MIP distribution
void fitSPEMIP(TH1* hfit, int run)
{
    TCanvas c1("c1","c1",800,800);
    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.12);
    gPad->SetRightMargin(0.05);
    gPad->SetLeftMargin(0.14);
    c1.SetLogx();
    c1.SetLogy();
    //hfit->GetXaxis()->SetRangeUser(3, 1000);
    //hfit->SetMinimum(0.0000001);
    //hfit->SetMaximum(200);
    
    char ffname[128];
    
    PPEFunc background(15,false);
    background.h = hfit;
    PPEFunc background_MIP(15,true);
    background_MIP.h = hfit;
                                                      //////////////Parameter Info////////////////////////////////
    double min0  =    100; double max0  =   450;      //p[0]  : pedestal amplitude-------------//               //
    double min1  =      8; double max1  =    20;      //p[1]  : pedestal width-----------------//               //
    double min2  =     50; double max2  =   100;      //p[2]  : Overall Shift------------------//(avg. pedestal)//
    double min3  =      0; double max3  =   0.1;      //p[3]  : background amplidude-----------//(not used)     //
    double min4  =      0; double max4  =    10;      //p[4]  : background width---------------//(not used)     //
    double min5  =  10000; double max5  = 25000;      //p[5]  : overall PE amplitude ----------//               //
    double min6  =      0; double max6  =  0.03;      //p[6]  : poisson mean number of PE------//               //
    double min7  =     35; double max7  =    65;      //p[7]  : PE peak spacing----------------//(Gain)         //
    double min8  =      2; double max8  =    10;      //p[8]  : PE peak width------------------//               //
    double min9  =    0.2; double max9  =   0.8;      //p[9]  : pixel cross-talk probability---//               //
    double min10 =    200; double max10 =  2000;      //p[10] : Total MIP peak area------------//               //
    double min11 =   8000; double max11 =  6000;      //p[11] : Most probable value of Landau--//(MPV)          //
    double min12 =     10; double max12 =   500;      //p[12] : Landau width parameter---------//               //
    double min13 =     20; double max13 =   500;      //p[13] : gaussian width-----------------//               //
    double set0  = 115;                               ////////////////////////////////////////////////////////////
    double set1  = 10;
    double set2  = 55;
    double set3  = 0;
    double set4  = 1;
    double set5  = 20000;
    double set6  = 0.02;
    double set7  =   55;
    double set8  =  3.0;
    double set9  =  0.5;
    double set10 =  500;
    double set11 = 1500;
    double set12 =  150;
    double set13 =  200;
    
    //Fit Pedestal peak
    sprintf(ffname, "sff%s", hfit->GetName());
    TF1* fit = new TF1(ffname, background, 0.0, 100000.0, 10);
    fit->SetParLimits(0, min0, max0); 
    fit->SetParLimits(1, min1, max1); 
    fit->SetParLimits(2, min2, max2); 
    fit->SetParLimits(3, min3, max3); 
    fit->SetParLimits(4, min4, max4); 
    fit->SetParLimits(5, min5, max5); 
    fit->SetParLimits(6, min6, max6); 
    fit->SetParLimits(7, min7, max7); 
    fit->SetParLimits(8, min8, max8); 
    fit->SetParLimits(9, min9, max9); 
    fit->SetParameter(0, set0);
    fit->SetParameter(1, set1);
    fit->SetParameter(2, set2);
    fit->FixParameter(3, set3);
    fit->FixParameter(4, set4);
    fit->FixParameter(5, set5);
    fit->FixParameter(6, set6);
    fit->FixParameter(7, set7);
    fit->FixParameter(8, set8);
    fit->FixParameter(9, set9);
    hfit->Fit(fit, "RNQL", "", 15, 100);
    
    //fit generic background amplitude 
    //sprintf(ffname, "sff2%s", hfit->GetName());
    //TF1* fit2 = new TF1(ffname, background, 0.0, 100000.0, 10);
    //fit2->SetParLimits(0,    0,   6000);
    //fit2->SetParLimits(1,    0,     15);
    //fit2->SetParLimits(2,    0,   1000);
    //fit2->SetParLimits(3,    0, 100000);
    //fit2->SetParLimits(4,    0,   1000);
    //fit2->SetParLimits(5,    0,   1000);
    //fit2->SetParLimits(6,    0,    3.0);
    //fit2->SetParLimits(7,    0,     55);
    //fit2->SetParLimits(8,    0,     20);
    //fit2->SetParLimits(9,    0,    1.0);
    //fit2->FixParameter(0, fit->GetParameter(0));
    //fit2->FixParameter(1, fit->GetParameter(1));
    //fit2->FixParameter(2, fit->GetParameter(2));
    //fit2->SetParameter(3,  100);
    //fit2->SetParameter(4,   50);
    //fit2->FixParameter(5,    3);
    //fit2->FixParameter(6,  0.9);
    //fit2->FixParameter(7,   41);
    //fit2->FixParameter(8, 4.06);
    //fit2->FixParameter(9, 0.15);
    //hfit->Fit(fit2, "RNQL", "", 0, 45);
    //
    ////Try to pin down highside tail
    //int i;
    //for(i =  hfit->GetNbinsX() + 1; i > 0; --i)
    //  {
    //	if(hfit->GetBinContent(i) > 0.00001) break;
    //  }
    //hfit->Fit(fit2, "RNQL", "", hfit->GetBinCenter(i)/5, hfit->GetBinCenter(i));
    //fit2->FixParameter(4, fit2->GetParameter(4));
    //hfit->Fit(fit2, "RNQL", "", 0, 45);

    //fit PE peaks
    sprintf(ffname, "sff3%s", hfit->GetName());
    TF1* fit3 = new TF1(ffname, background, 0.0, 100000.0, 10);
    fit3->SetParLimits(0, min0, max0); 
    fit3->SetParLimits(1, min1, max1); 
    fit3->SetParLimits(2, min2, max2); 
    fit3->SetParLimits(3, min3, max3); 
    fit3->SetParLimits(4, min4, max4); 
    fit3->SetParLimits(5, min5, max5); 
    fit3->SetParLimits(6, min6, max6); 
    fit3->SetParLimits(7, min7, max7); 
    fit3->SetParLimits(8, min8, max8); 
    fit3->SetParLimits(9, min9, max9); 
    fit3->SetParameter(0, fit->GetParameter(0));
    fit3->SetParameter(1, fit->GetParameter(1));
    fit3->SetParameter(2, fit->GetParameter(2));
    fit3->FixParameter(3, fit->GetParameter(3));
    fit3->FixParameter(4, fit->GetParameter(4));
    fit3->SetParameter(5, set5);
    fit3->SetParameter(6, set6);
    fit3->SetParameter(7, set7);
    fit3->SetParameter(8, set8);
    fit3->SetParameter(9, set9);
    hfit->Fit(fit3, "RNQL", "", fit->GetParameter(2)-20, 190);

    //Fit MIP peak
    sprintf(ffname, "sff4%s", hfit->GetName());
    TF1* fit4 = new TF1(ffname, background_MIP, 0.0, 100000.0, 14);
    ///TF1* fit4 = new TF1(ffname, background, 0.0, 100000.0, 14);
    fit4->SetParLimits(0,  min0,  max0); 
    fit4->SetParLimits(1,  min1,  max1); 
    fit4->SetParLimits(2,  min2,  max2); 
    fit4->SetParLimits(3,  min3,  max3); 
    fit4->SetParLimits(4,  min4,  max4); 
    fit4->SetParLimits(5,  min5,  max5); 
    fit4->SetParLimits(6,  min6,  max6); 
    fit4->SetParLimits(7,  min7,  max7); 
    fit4->SetParLimits(8,  min8,  max8); 
    fit4->SetParLimits(9,  min9,  max9);
    fit4->SetParLimits(10, min10, max10); 
    fit4->SetParLimits(11, min11, max11); 
    fit4->SetParLimits(12, min12, max12); 
    fit4->SetParLimits(13, min13, max13); 
    fit4->SetParameter(0,  fit3->GetParameter(0));
    fit4->SetParameter(1,  fit3->GetParameter(1));
    fit4->SetParameter(2,  fit3->GetParameter(2));
    fit4->FixParameter(3,  fit3->GetParameter(3));
    fit4->FixParameter(4,  fit3->GetParameter(4));
    fit4->SetParameter(5,  fit3->GetParameter(5));
    fit4->SetParameter(6,  fit3->GetParameter(6));
    fit4->SetParameter(7,  fit3->GetParameter(7));
    fit4->SetParameter(8,  fit3->GetParameter(8));
    fit4->SetParameter(9,  fit3->GetParameter(9));
    fit4->SetParameter(10, set10);
    fit4->SetParameter(11, set11);
    fit4->SetParameter(12, set12);
    fit4->SetParameter(13, set13);
    hfit->Fit(fit4, "RNQ", "", fit3->GetParameter(2)-20, fit3->GetParameter(2)+4000);

    //Fine Tunning Fit
    sprintf(ffname, "sff4%s", hfit->GetName());
    TF1* fit5 = new TF1(ffname, background_MIP, 0.0, 100000.0, 14);
    ///TF1* fit5 = new TF1(ffname, background, 0.0, 100000.0, 14);
    double down = 0.5;
    double up   = 1.5;
    fit5->SetParLimits(0, down*fit4->GetParameter(0),  up*fit4->GetParameter(0));
    fit5->SetParLimits(1, down*fit4->GetParameter(1),  up*fit4->GetParameter(1));
    fit5->SetParLimits(2, down*fit4->GetParameter(2),  up*fit4->GetParameter(2));
    fit5->SetParLimits(3, down*fit4->GetParameter(3),  up*fit4->GetParameter(3));
    fit5->SetParLimits(4, down*fit4->GetParameter(4),  up*fit4->GetParameter(4));
    fit5->SetParLimits(5, down*fit4->GetParameter(5),  up*fit4->GetParameter(5));
    fit5->SetParLimits(6, down*fit4->GetParameter(6),     fit4->GetParameter(6));
    fit5->SetParLimits(7, down*fit4->GetParameter(7),  up*fit4->GetParameter(7));
    fit5->SetParLimits(8, down*fit4->GetParameter(8),  up*fit4->GetParameter(8));
    fit5->SetParLimits(9, down*fit4->GetParameter(9),  up*fit4->GetParameter(9));
    fit5->SetParLimits(10,down*fit4->GetParameter(10), up*fit4->GetParameter(10));
    fit5->SetParLimits(11,down*fit4->GetParameter(11), up*fit4->GetParameter(11));
    fit5->SetParLimits(12,down*fit4->GetParameter(12), up*fit4->GetParameter(12));
    fit5->SetParLimits(13,down*fit4->GetParameter(13), up*fit4->GetParameter(13));
    fit5->SetParameter(0,  fit4->GetParameter(0));
    fit5->SetParameter(1,  fit4->GetParameter(1));
    fit5->SetParameter(2,  fit4->GetParameter(2));
    fit5->FixParameter(3,  fit4->GetParameter(3));
    fit5->FixParameter(4,  fit4->GetParameter(4));
    fit5->SetParameter(5,  fit4->GetParameter(5));
    fit5->SetParameter(6,  fit4->GetParameter(6));
    fit5->SetParameter(7,  fit4->GetParameter(7));
    fit5->SetParameter(8,  fit4->GetParameter(8));
    fit5->SetParameter(9,  fit4->GetParameter(9));
    fit5->SetParameter(10, fit4->GetParameter(10));
    fit5->SetParameter(11, fit4->GetParameter(11));
    fit5->SetParameter(12, fit4->GetParameter(12));
    fit5->SetParameter(13, fit4->GetParameter(13));
    hfit->Fit(fit5, "RNQ", "", fit3->GetParameter(2)-20, fit4->GetParameter(2)+5000);
    
    ///printf("Chi^2:%10.4f Gain:%10.4f Unc:%10.4f MPV:%10.4f Unc:%10.4f", fit5->GetChisquare(), fit5->GetParameter(6), fit5->GetParError(6), fit5->GetParameter(10), fit5->GetParError(10));
    ///printf("  Run: %i Hname: %s\n",run, hfit->GetName());

    //h.GetXaxis().SetRangeUser(2, 500)
    //std::cout<<"ans: "<<hfit->GetName()<<"  "<<fit5->GetParameter(10)<<"  "<<fit5->GetParError(10)<<"   "<<fit5->GetParameter(6)<<std::endl;
    printf("Chi^2:%10.4f, P0:%10.4f, P1:%10.4f, P2:%10.4f, P5:%10.4f, P6:%10.4f, P7:%10.4f, P8:%10.4f, P9:%10.4f, P10:%10.4f, P11:%10.4f, P12:%10.4f, P13:%10.4f\n", fit5->GetChisquare(), fit5->GetParameter(0), fit5->GetParameter(1), fit5->GetParameter(2), fit5->GetParameter(5), fit5->GetParameter(6), fit5->GetParameter(7), fit5->GetParameter(8), fit5->GetParameter(9), fit5->GetParameter(10), fit5->GetParameter(11), fit5->GetParameter(12), fit5->GetParameter(13));
    hfit->Draw("hist");
    hfit->GetYaxis()->SetTitle("Events /fC");
    hfit->GetXaxis()->SetTitle("Charge [fC]");
    hfit->SetTitleOffset(1,"X");
    hfit->SetTitleOffset(1.2,"Y");
    hfit->SetTitleSize(0.05,"X");
    hfit->SetTitleSize(0.05,"Y");
    fit3->SetLineColor(kRed);
    fit3->SetLineWidth(2);
    fit5->SetLineWidth(2);
    TF1* fit6 = new TF1("BackGround", background, 0.0, 100000.0, 10);
    fit6->FixParameter(0,  fit5->GetParameter(0));
    fit6->FixParameter(1,  fit5->GetParameter(1));
    fit6->FixParameter(2,  fit5->GetParameter(2));
    fit6->FixParameter(3,  fit5->GetParameter(3));
    fit6->FixParameter(4,  fit5->GetParameter(4));
    fit6->FixParameter(5,  fit5->GetParameter(5));
    fit6->FixParameter(6,  fit5->GetParameter(6));
    fit6->FixParameter(7,  fit5->GetParameter(7));
    fit6->FixParameter(8,  fit5->GetParameter(8));
    fit6->FixParameter(9,  fit5->GetParameter(9));

    TF1* fit7 = new TF1("Lang", langaufun, 0.0, 100000.0, 4);
    fit7->FixParameter(0,  fit5->GetParameter(12));
    fit7->FixParameter(1,  fit5->GetParameter(11));
    fit7->FixParameter(2,  fit5->GetParameter(10));
    fit7->FixParameter(3,  fit5->GetParameter(13));
    //fit3->Draw("same");
    fit6->SetLineColor(kBlue);
    fit6->Draw("same");
    fit7->SetLineColor(kGreen+2);
    fit7->Draw("same");
    fit5->Draw("same");
    
    //Make Plots Pretty
    hfit->SetStats(false);
    hfit->SetTitle("");
    
    TLatex* CMSPrelim1 = new TLatex(0.14, 0.91, "CMS #scale[0.9]{#font[52]{Preliminary}}");
    CMSPrelim1->SetNDC();
    CMSPrelim1->SetTextFont(62);

    TLatex* testbeam = new TLatex(0.95, 0.91, "Testbeam 2017");
    testbeam->SetNDC();
    testbeam->SetTextFont(42);
    testbeam->SetTextAlign(31);

    TLatex* SiPMTitle = new TLatex(0.93, 0.86, "SiPM (Silicon Photomultiplier)");
    SiPMTitle->SetNDC();
    SiPMTitle->SetTextFont(42);
    SiPMTitle->SetTextAlign(32);

    TLatex* Muon = new TLatex(0.93, 0.81, "Parasitic Muons");
    Muon->SetNDC();
    Muon->SetTextFont(42);
    Muon->SetTextAlign(32);

    char chan [100];
    int iEta, iPhi, iDepth;
    //sscanf (hfit->GetName(),"beam_adc_%d_%d_%d", &iEta, &iPhi, &iDepth);
    sscanf (hfit->GetName(),"adc_nosub_%d_%d_%d", &iEta, &iPhi, &iDepth);
    //sscanf (hfit->GetName(),"ped_adc_%d_%d_%d", &iEta, &iPhi, &iDepth);
    sprintf (chan, "Channel: %d,%d,%d", iEta, iPhi, iDepth);
    
    TLatex* channel = new TLatex(0.93, 0.76, chan);
    channel->SetNDC();
    channel->SetTextFont(42);
    channel->SetTextAlign(32);

    char mpv [100];
    float intmpv = fit5->GetParameter(11);
    sprintf (mpv,"MPV: %0.3f", intmpv);

    TLatex* MPV = new TLatex(0.93, 0.7, mpv);
    MPV->SetNDC();
    MPV->SetTextFont(42);
    MPV->SetTextAlign(31);

    char gain [100];
    float intgain = fit5->GetParameter(7);
    sprintf (gain,"Gain: %0.3f", intgain);

    TLatex* Gain = new TLatex(0.93, 0.65, gain);
    Gain->SetNDC();
    Gain->SetTextFont(42);
    Gain->SetTextAlign(31);

    CMSPrelim1->Draw();
    testbeam->Draw();
    SiPMTitle->Draw();
    Muon->Draw();
    channel->Draw();
    MPV->Draw();
    Gain->Draw();
    
    char oname[128];
    sprintf(oname, "%s_%i.pdf", hfit->GetName(), run);
    //sprintf(oname, "%s_SiPMRuns_3030to3475.pdf", hfit->GetName());
    //sprintf(oname, "%s_HBRuns_3526to3534.pdf", hfit->GetName());
    c1.Print(oname);
}

int main()
{
    gROOT->SetStyle("Plain");    

    //char fname[] = "Test_run3234_3228.root";
    //char fname[] = "data/analysis_run00%i_EMAP_AC05_01AUG2017_NOMINAL.root";    
    //char fname[] = "All_SiPM_Data.root";
    //char fname[] = "theFile.root";

    //char fname[] = "muonHBRuns.root";
    char fname[] = "../MipCalibration/analysis_run00%i_EMAP_AC07_01SEP2017_Phase1_HB.root";
    //char fname[] = "analysis_run%i_EMAP_AC07_17SEP2017_Phase1_HB_B904.root";
    //char fname[] = "analysis_run00%i_EMAP-6SEP2017_HBHE.root";
    
    //for(int irun = 8849; irun < 8850; ++irun) ///Joe Test
    //for(int irun = 9998; irun < 9999; ++irun) ///All SiPM Data
    for(int irun = 3526; irun < 3571; ++irun)
    //for(int irun = 1000028088; irun < 1000028106; ++irun)
    {
        char fn[128];
        sprintf(fn, fname, irun);
        TFile *fin = TFile::Open(fn);
    
        if(!fin) continue;
    
        TDirectory* din = (TDirectory*)fin->Get("hcalADCHists");
	//TDirectory* din = (TDirectory*)fin;
        
        TIter next(din->GetListOfKeys());
        TObject *obj;
        while((obj = next()))
        {
	    //if(strstr(obj->GetName(), "beam_adc_6") != nullptr || strstr(obj->GetName(), "beam_adc_9") != nullptr || strstr(obj->GetName(), "beam_adc_11") != nullptr || strstr(obj->GetName(), "beam_adc_13") != nullptr)
	    if(strstr(obj->GetName(), "adc_nosub_6") != nullptr || strstr(obj->GetName(), "adc_nosub_9") != nullptr || strstr(obj->GetName(), "adc_nosub_11") != nullptr || strstr(obj->GetName(), "adc_nosub_13") != nullptr)
	    //if(strstr(obj->GetName(), "ped_adc_6") != nullptr || strstr(obj->GetName(), "ped_adc_9") != nullptr || strstr(obj->GetName(), "ped_adc_11") != nullptr || strstr(obj->GetName(), "ped_adc_13") != nullptr)
            {
                TH1* htofit = (TH1*)din->Get(obj->GetName());
                if(htofit->GetMean() > 100+80.0)//50
                {
		    //printf("Hname: %s/n", htofit->GetName());
                    //h = htofit;
                    fitSPEMIP(htofit, irun);
                }
            }
        }
    }
}
