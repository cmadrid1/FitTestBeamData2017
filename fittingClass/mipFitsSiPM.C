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
    //hfit->GetXaxis()->SetRangeUser(10, 400);
    //hfit->SetMinimum(0.0000001);
    //hfit->SetMaximum(200);
    
    char ffname[128];
    double nPeaks = 10;
    
    PPEFunc background(nPeaks,false);
    background.h = hfit;
    PPEFunc background_MIP(nPeaks,true);
    background_MIP.h = hfit;
    hfit->Draw("hist");
                                                       //////////////Parameter Info////////////////////////////////
    double min0  =   100; double max0  =    1000;      //p[0]  : pedestal amplitude-------------//               //
    double min1  =     1; double max1  =      10;      //p[1]  : pedestal width-----------------//               //
    double min2  =    50; double max2  =     100;      //p[2]  : Overall Shift------------------//(avg. pedestal)//
    double min3  =     0; double max3  =    1000;      //p[3]  : background amplidude-----------//               //
    double min4  =     0; double max4  =    1000;      //p[4]  : background width---------------//               //
    double min5  =  1000; double max5  =   25000;      //p[5]  : overall PE amplitude ----------//               //
    double min6  = 0.001; double max6  =    0.03;      //p[6]  : poisson mean number of PE------//               //
    double min7  =    40; double max7  =      55;      //p[7]  : PE peak spacing----------------//(Gain)         //
    double min8  =  0.01; double max8  =      10;      //p[8]  : PE peak width------------------//               //
    double min9  =  0.01; double max9  =    0.95;      //p[9]  : pixel cross-talk probability---//               //
    double min10 =   200; double max10 =    2000;      //p[10] : Total MIP peak area------------//               //
    double min11 =   250; double max11 =    5000;      //p[11] : Most probable value of Landau--//(MPV)          //
    double min12 =    10; double max12 =     500;      //p[12] : Landau width parameter---------//               //
    double min13 =    20; double max13 =     500;      //p[13] : gaussian width-----------------//               //
    double set0  =   115;                              ////////////////////////////////////////////////////////////
    double set1  =     8;
    double set2  =    55;
    double set3  =     0;
    double set4  =     1;
    double set5  = 20000;
    double set6  =  0.02;
    double set7  =    45;
    double set8  =   1.0;
    double set9  =   0.4;
    double set10 =   500;
    double set11 =   800;
    double set12 =   100;
    double set13 =   200;
    
    //Fit Pedestal Peak
    sprintf(ffname, "sff%s", hfit->GetName());
    TF1* fit1 = new TF1(ffname, background, 0.0, 100000.0, 10);
    fit1->SetParLimits(0,  min0, max0); 
    fit1->SetParLimits(1,  min1, max1); 
    fit1->SetParLimits(2,  min2, max2); 
    fit1->SetParLimits(3,  min3, max3); 
    fit1->SetParLimits(4,  min4, max4); 
    fit1->SetParLimits(5,  min5, max5); 
    fit1->SetParLimits(6,  min6, max6); 
    fit1->SetParLimits(7,  min7, max7); 
    fit1->SetParLimits(8,  min8, max8); 
    fit1->SetParLimits(9,  min9, max9); 
    //fit1->SetParLimits(10, min10, max10); 
    //fit1->SetParLimits(11, min11, max11); 
    //fit1->SetParLimits(12, min12, max12); 
    //fit1->SetParLimits(13, min13, max13); 
    fit1->SetParameter(0,  set0);
    fit1->SetParameter(1,  set1);
    fit1->SetParameter(2,  set2);
    fit1->FixParameter(3,  set3);
    fit1->FixParameter(4,  set4);
    fit1->FixParameter(5,  0.01);
    fit1->FixParameter(6,  set6);
    fit1->FixParameter(7,  set7);
    fit1->FixParameter(8,  set8);
    fit1->FixParameter(9,  set9);
    //fit1->FixParameter(10, set10);
    //fit1->FixParameter(11, set11);
    //fit1->FixParameter(12, 0.01);
    //fit1->FixParameter(13, 0.01);
    hfit->Fit(fit1, "RQMWL", "", 10, 100);
    fit1->SetLineColor(kRed);
    fit1->SetLineWidth(2);
    fit1->Draw("same");
    
    //fit PE Peaks 
    sprintf(ffname, "sff2%s", hfit->GetName());
    TF1* fit2 = new TF1(ffname, background, 0.0, 100000.0, 10);
    fit2->SetParLimits(0,  min0, max0); 
    fit2->SetParLimits(1,  min1, max1); 
    fit2->SetParLimits(2,  min2, max2); 
    fit2->SetParLimits(3,  min3, max3); 
    fit2->SetParLimits(4,  min4, max4); 
    fit2->SetParLimits(5,  min5, max5); 
    fit2->SetParLimits(6,  min6, max6); 
    fit2->SetParLimits(7,  min7, max7); 
    fit2->SetParLimits(8,  min8, max8); 
    fit2->SetParLimits(9,  min9, max9); 
    //fit2->SetParLimits(10, min10, max10); 
    //fit2->SetParLimits(11, min11, max11); 
    //fit2->SetParLimits(12, min12, max12); 
    //fit2->SetParLimits(13, min13, max13); 
    fit2->FixParameter(0,  fit1->GetParameter(0));
    fit2->FixParameter(1,  fit1->GetParameter(1));
    fit2->FixParameter(2,  fit1->GetParameter(2));
    fit2->FixParameter(3,  set3);
    fit2->FixParameter(4,  set4);
    fit2->SetParameter(5,  set5);
    fit2->SetParameter(6,  set6);
    fit2->SetParameter(7,  set7);
    fit2->SetParameter(8,  set8);
    fit2->SetParameter(9,  set9);
    //fit2->FixParameter(10, set10);
    //fit2->FixParameter(11, set11);
    //fit2->FixParameter(12, 0.01);
    //fit2->FixParameter(13, 0.01);
    hfit->Fit(fit2, "RQMWL", "", fit1->GetParameter(2), 400);
    fit2->SetLineColor(kGray);
    fit2->SetLineWidth(2);
    fit2->Draw("same");
    
    //Fine tune PE peaks
    sprintf(ffname, "sff3%s", hfit->GetName());
    TF1* fit3 = new TF1(ffname, background, 0.0, 100000.0, 10);
    double down1 = 0.7;
    double up1   = 1.3;
    fit3->SetParLimits(0,  down1*fit2->GetParameter(0), up1*fit2->GetParameter(0)); 
    fit3->SetParLimits(1,  down1*fit2->GetParameter(1), up1*fit2->GetParameter(1)); 
    fit3->SetParLimits(2,  down1*fit2->GetParameter(2), up1*fit2->GetParameter(2)); 
    fit3->SetParLimits(3,  down1*fit2->GetParameter(3), up1*fit2->GetParameter(3)); 
    fit3->SetParLimits(4,  down1*fit2->GetParameter(4), up1*fit2->GetParameter(4)); 
    fit3->SetParLimits(5,  down1*fit2->GetParameter(5), up1*fit2->GetParameter(5)); 
    fit3->SetParLimits(6,  down1*fit2->GetParameter(6), up1*fit2->GetParameter(6)); 
    fit3->SetParLimits(7,  down1*fit2->GetParameter(7), up1*fit2->GetParameter(7)); 
    fit3->SetParLimits(8,  down1*fit2->GetParameter(8), up1*fit2->GetParameter(8)); 
    fit3->SetParLimits(9,  down1*fit2->GetParameter(9), up1*fit2->GetParameter(9)); 
    //fit3->SetParLimits(10, min10, max10); 
    //fit3->SetParLimits(11, min11, max11); 
    //fit3->SetParLimits(12, min12, max12); 
    //fit3->SetParLimits(13, min13, max13); 
    fit3->SetParameter(0,  fit2->GetParameter(0));
    fit3->SetParameter(1,  fit2->GetParameter(1));
    fit3->SetParameter(2,  fit2->GetParameter(2));
    fit3->FixParameter(3,  fit2->GetParameter(3));
    fit3->FixParameter(4,  fit2->GetParameter(4));
    fit3->SetParameter(5,  fit2->GetParameter(5));
    fit3->SetParameter(6,  fit2->GetParameter(6));
    fit3->SetParameter(7,  fit2->GetParameter(7));
    fit3->SetParameter(8,  fit2->GetParameter(8));
    fit3->SetParameter(9,  fit2->GetParameter(9));
    //fit3->FixParameter(10, set10);
    //fit3->FixParameter(11, set11);
    //fit3->FixParameter(12, 0.01);
    //fit3->FixParameter(13, 0.01);
    hfit->Fit(fit3, "RQMWL", "", fit2->GetParameter(2), 400);
    fit3->SetLineColor(kBlue);
    fit3->SetLineWidth(2);
    fit3->Draw("same");
    
    //Fit MIP peak
    sprintf(ffname, "sff4%s", hfit->GetName());
    TF1* fit4 = new TF1(ffname, background_MIP, 0.0, 100000.0, 14);
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
    fit4->FixParameter(0,  fit3->GetParameter(0));
    fit4->FixParameter(1,  fit3->GetParameter(1));
    fit4->FixParameter(2,  fit3->GetParameter(2));
    fit4->FixParameter(3,  fit3->GetParameter(3));
    fit4->FixParameter(4,  fit3->GetParameter(4));
    fit4->FixParameter(5,  fit3->GetParameter(5));
    fit4->FixParameter(6,  fit3->GetParameter(6));
    fit4->FixParameter(7,  fit3->GetParameter(7));
    fit4->FixParameter(8,  fit3->GetParameter(8));
    fit4->FixParameter(9,  fit3->GetParameter(9));
    fit4->SetParameter(10, set10);
    fit4->SetParameter(11, set11);
    fit4->SetParameter(12, set12);
    fit4->SetParameter(13, set13);
    hfit->Fit(fit4, "RQMWL", "", fit3->GetParameter(2), 4000);
    fit4->SetLineColor(kOrange);
    fit4->SetLineWidth(2);
    fit4->Draw("same");
    
    //Fine Tunning All Fit
    sprintf(ffname, "sff4%s", hfit->GetName());
    TF1* fit5 = new TF1(ffname, background_MIP, 0.0, 100000.0, 14);
    double down = 0.7;
    double up   = 1.3;
    fit5->SetParLimits(0, down*fit4->GetParameter(0),  up*fit4->GetParameter(0));
    fit5->SetParLimits(1, down*fit4->GetParameter(1),  up*fit4->GetParameter(1));
    fit5->SetParLimits(2,      fit4->GetParameter(2),  up*fit4->GetParameter(2));
    fit5->SetParLimits(3, down*fit4->GetParameter(3),  up*fit4->GetParameter(3));
    fit5->SetParLimits(4, down*fit4->GetParameter(4),  up*fit4->GetParameter(4));
    fit5->SetParLimits(5, down*fit4->GetParameter(5),  up*fit4->GetParameter(5));
    fit5->SetParLimits(6, down*fit4->GetParameter(6),  up*fit4->GetParameter(6));
    fit5->SetParLimits(7, down*fit4->GetParameter(7),     fit4->GetParameter(7));
    fit5->SetParLimits(8, down*fit4->GetParameter(8),  up*fit4->GetParameter(8));
    fit5->SetParLimits(9, down*fit4->GetParameter(9),  up*fit4->GetParameter(9));
    fit5->SetParLimits(10,down*fit4->GetParameter(10),    fit4->GetParameter(10));
    fit5->SetParLimits(11,     fit4->GetParameter(11), up*fit4->GetParameter(11));
    fit5->SetParLimits(12,down*fit4->GetParameter(12),    fit4->GetParameter(12));
    fit5->SetParLimits(13,down*fit4->GetParameter(13),    fit4->GetParameter(13));
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
    hfit->Fit(fit5, "RQMWL", "", fit4->GetParameter(2), 4000);
    fit5->SetLineColor(kBlack);
    fit5->SetLineWidth(2);
    fit5->Draw("same");
    
    printf("Chi^2:%10.4f Gain:%10.4f Unc:%10.4f MPV:%10.4f Unc:%10.4f", fit5->GetChisquare(), fit5->GetParameter(7), fit5->GetParError(7), fit5->GetParameter(11), fit5->GetParError(11));
    printf("  Run: %i Hname: %s\n",run, hfit->GetName());

    //std::cout<<"ans: "<<hfit->GetName()<<"  "<<fit5->GetParameter(10)<<"  "<<fit5->GetParError(10)<<"   "<<fit5->GetParameter(6)<<std::endl;
    //printf("Chi^2:%10.4f, P0:%10.4f, P1:%10.4f", fit5->GetChisquare(), fit5->GetParameter(0), fit5->GetParameter(1));
    //printf(" P2:%10.4f, P3:%10.4f, P4:%10.4f, P5:%10.4f", fit5->GetParameter(2), fit5->GetParameter(3), fit5->GetParameter(4), fit5->GetParameter(5));
    //printf(" P6:%10.4f, P7:%10.4f, P8:%10.4f, P9:%10.4f", fit5->GetParameter(6), fit5->GetParameter(7), fit5->GetParameter(8), fit5->GetParameter(9));
    //printf(" P10:%10.4f, P11:%10.4f, P12:%10.4f, P13:%10.4f\n", fit5->GetParameter(10), fit5->GetParameter(11), fit5->GetParameter(12), fit5->GetParameter(13));
    hfit->GetYaxis()->SetTitle("Events /fC");
    hfit->GetXaxis()->SetTitle("Charge [fC]");
    hfit->SetTitleOffset(1,"X");
    hfit->SetTitleOffset(1.2,"Y");
    hfit->SetTitleSize(0.05,"X");
    hfit->SetTitleSize(0.05,"Y");
    
    //TF1* draw5PE = new TF1("BackGround", background, 0.0, 100000.0, 10);
    //draw5PE->FixParameter(0,  fit5->GetParameter(0));
    //draw5PE->FixParameter(1,  fit5->GetParameter(1));
    //draw5PE->FixParameter(2,  fit5->GetParameter(2));
    //draw5PE->FixParameter(3,  fit5->GetParameter(3));
    //draw5PE->FixParameter(4,  fit5->GetParameter(4));
    //draw5PE->FixParameter(5,  fit5->GetParameter(5));
    //draw5PE->FixParameter(6,  fit5->GetParameter(6));
    //draw5PE->FixParameter(7,  fit5->GetParameter(7));
    //draw5PE->FixParameter(8,  fit5->GetParameter(8));
    //draw5PE->FixParameter(9,  fit5->GetParameter(9));
    ////draw5PE->FixParameter(10, fit5->GetParameter(10));
    ////draw5PE->FixParameter(11, fit5->GetParameter(11));
    ////draw5PE->FixParameter(12, fit5->GetParameter(12));
    ////draw5PE->FixParameter(13, fit5->GetParameter(13));
    //draw5PE->SetLineWidth(2);
    //draw5PE->SetLineColor(kBlue);
    //draw5PE->Draw("same");

    TF1* draw5Mip = new TF1("Lang", langaufun, 0.0, 100000.0, 4);
    draw5Mip->FixParameter(0, fit5->GetParameter(12));
    draw5Mip->FixParameter(1, fit5->GetParameter(11));
    draw5Mip->FixParameter(2, fit5->GetParameter(10));
    draw5Mip->FixParameter(3, fit5->GetParameter(13));
    draw5Mip->SetLineWidth(2);
    draw5Mip->SetLineColor(kGreen+2);
    draw5Mip->Draw("same");
    
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
    //sscanf (hfit->GetName(),"adc_nosub_%d_%d_%d", &iEta, &iPhi, &iDepth);
    sscanf (hfit->GetName(),"adc_nosub_constBin_%d_%d_%d", &iEta, &iPhi, &iDepth);
    //sscanf (hfit->GetName(),"ped_adc_%d_%d_%d", &iEta, &iPhi, &iDepth);
    sprintf (chan, "Channel: %d,%d,%d", iEta, iPhi, iDepth);
    
    TLatex* channel = new TLatex(0.93, 0.86, chan);
    channel->SetNDC();
    channel->SetTextFont(42);
    channel->SetTextAlign(32);

    char mpv [100];
    float intmpv = fit5->GetParameter(11);
    sprintf (mpv,"MPV: %0.3f", intmpv);

    TLatex* MPV = new TLatex(0.93, 0.8, mpv);
    MPV->SetNDC();
    MPV->SetTextFont(42);
    MPV->SetTextAlign(31);

    char gain [100];
    float intgain = fit5->GetParameter(7);
    sprintf (gain,"Gain: %0.3f", intgain);

    TLatex* Gain = new TLatex(0.93, 0.75, gain);
    Gain->SetNDC();
    Gain->SetTextFont(42);
    Gain->SetTextAlign(31);

    CMSPrelim1->Draw();
    testbeam->Draw();
    //SiPMTitle->Draw();
    //Muon->Draw();
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
	  int iEta, iPhi, iDepth;
	  sscanf (obj->GetName(),"adc_nosub_constBin_%d_%d_%d", &iEta, &iPhi, &iDepth);

	  //if(strstr(obj->GetName(), "beam_adc_6") != nullptr || strstr(obj->GetName(), "beam_adc_9") != nullptr || strstr(obj->GetName(), "beam_adc_11") != nullptr || strstr(obj->GetName(), "beam_adc_13") != nullptr)
	    //if(strstr(obj->GetName(), "adc_nosub_6") != nullptr || strstr(obj->GetName(), "adc_nosub_9") != nullptr || strstr(obj->GetName(), "adc_nosub_11") != nullptr || strstr(obj->GetName(), "adc_nosub_13") != nullptr)
	    if((strstr(obj->GetName(), "adc_nosub_constBin") != nullptr) && (iPhi == 5) && (iDepth == 2 || iDepth == 3)) 
	    //if(strstr(obj->GetName(), "ped_adc_6") != nullptr || strstr(obj->GetName(), "ped_adc_9") != nullptr || strstr(obj->GetName(), "ped_adc_11") != nullptr || strstr(obj->GetName(), "ped_adc_13") != nullptr)
            {
                TH1* htofit = (TH1*)din->Get(obj->GetName());
                if(htofit->GetMean() > 100)//50
                {
		    //printf("Hname: %s/n", htofit->GetName());
                    //h = htofit;
                    fitSPEMIP(htofit, irun);
                }
            }
        }
    }
}
