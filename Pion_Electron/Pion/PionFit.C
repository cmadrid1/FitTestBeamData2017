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
#include "TH1F.h"

#include <cmath>
#include <cstdio>


// landau-gaussian convolution using numerical integration
Double_t langaufun(Double_t *x, Double_t *par) {

   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

      // Numeric constants
      const Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      const Double_t mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      const Double_t np = 100.0;      // number of convolution steps
      const Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;

      // MP shift correction
      mpc = par[1] - mpshift * par[0];

      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];

      step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);

         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }

      return (par[2] * step * sum * invsq2pi / par[3]);
}

Double_t gaufun(Double_t *x, Double_t *par) {
  //par[0]=mean
  //par[1]=sigma
  //par[2]=amplitude
  
  return par[2]*TMath::Gaus(x[0],par[0],par[1]);
}

//More ugly globals to hold utility functions used in fits
TH1* h;

//This function fits the SiPM MIP distribution
void fitSPEMIP(TH1* hfit, int run)
{
    TCanvas c1("c1","c1",800,800);
    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.12);
    gPad->SetRightMargin(0.05);
    gPad->SetLeftMargin(0.14);
    //c1.SetLogx();
    //c1.SetLogy();
    hfit->GetXaxis()->SetRangeUser(3, 10500);
    
    char ffname[128];

    //Fit Pion Peak
    //sprintf(ffname, "sff%s", hfit->GetName());
    TF1* fit = new TF1("First Try", gaufun, 0.0, 10000.0, 3);
    fit->SetParLimits(0,    0, 2500);
    fit->SetParLimits(1,   80,  250);
    fit->SetParLimits(2,   20,  600);
    fit->SetParameter(0,  2000);
    fit->SetParameter(1,   130);
    fit->SetParameter(2,   400);
    hfit->Fit(fit, "RNQL", "", 100, 4000);
    
    //Fine Tunning Fit
    //sprintf(ffname, "sff%s", hfit->GetName());
    TF1* fit2 = new TF1("Fine Tunning", gaufun, 0.0, 10000.0, 3);
    double down = 0.1;
    double up   = 2.0;
    double me = fit->GetParameter(0);
    fit2->SetParLimits(0, down*fit->GetParameter(0),  up*fit->GetParameter(0));
    fit2->SetParLimits(1, down*fit->GetParameter(1),  up*fit->GetParameter(1));
    fit2->SetParLimits(2, down*fit->GetParameter(2),  up*fit->GetParameter(2));
    fit2->SetParameter(0,  fit->GetParameter(0));
    fit2->SetParameter(1,  fit->GetParameter(1));
    fit2->SetParameter(2,  fit->GetParameter(2));
    hfit->Fit(fit2, "RNQ", "",0.9*me, 1.3*me);

    //Fine Tunning Fit agian
    //sprintf(ffname, "sff%s", hfit->GetName());
    TF1* fit3 = new TF1("Fine Tunning two", gaufun, 0.0, 10000.0, 3);
    double d = 0.1;
    double u = 2.5;
    double m = fit2->GetParameter(0);
    double s = fit2->GetParameter(1)/1.7;
    fit3->SetParLimits(0, d*fit2->GetParameter(0),  u*fit2->GetParameter(0));
    fit3->SetParLimits(1, d*fit2->GetParameter(1),  u*fit2->GetParameter(1));
    fit3->SetParLimits(2, d*fit2->GetParameter(2),  u*fit2->GetParameter(2));
    fit3->SetParameter(0,  fit2->GetParameter(0));
    fit3->SetParameter(1,  fit2->GetParameter(1));
    fit3->SetParameter(2,  fit2->GetParameter(2));
    hfit->Fit(fit3, "RNQ", "", m-s, m+s);

    //Fine Tunning Fit agian again
    //sprintf(ffname, "sff%s", hfit->GetName());
    TF1* fit4  = new TF1("Fine Tunning three", gaufun, 0.0, 10000.0, 3);
    double d4 = 0.8;
    double u4 = 1.2;
    double m4 = fit3->GetParameter(0);
    double s4 = fit3->GetParameter(1);
    fit4->SetParLimits(0, d4*fit3->GetParameter(0),  u4*fit3->GetParameter(0));
    fit4->SetParLimits(1, d4*fit3->GetParameter(1),  u4*fit3->GetParameter(1));
    fit4->SetParLimits(2, d4*fit3->GetParameter(2),  u4*fit3->GetParameter(2));
    fit4->SetParameter(0,  fit3->GetParameter(0));
    fit4->SetParameter(1,  fit3->GetParameter(1));
    fit4->SetParameter(2,  fit3->GetParameter(2));
    hfit->Fit(fit4, "RNQ", "", m4-s4, m4+s4);
    hfit->Fit(fit4, "RNQ", "", m4-s4, m4+s4);
    hfit->Fit(fit4, "RNQ", "", m4-s4, m4+s4);
    hfit->Fit(fit4, "RNQ", "", m4-s4, m4+s4);
    
    printf("Chi^2: %10.4f P0: %10.4f Unc: %10.4f P1: %10.4f Unc: %10.4f P2: %10.4f Unc: %10.4f ", fit4->GetChisquare(), fit4->GetParameter(0),fit4->GetParError(0), fit4->GetParameter(1),fit4->GetParError(1), fit4->GetParameter(2),fit4->GetParError(2));
    printf("Run: %i Hname: %s\n",run, hfit->GetName());
    //std::cout<<"ans: "<<hfit->GetName()<<"  "<<fit4->GetParameter(10)<<"  "<<fit4->GetParError(10)<<std::endl;

    hfit->GetYaxis()->SetTitle("Events /fC");
    hfit->GetXaxis()->SetTitle("Charge [fC]");
    hfit->SetTitleOffset(1,"X");
    hfit->SetTitleOffset(1.2,"Y");
    hfit->SetTitleSize(0.05,"X");
    hfit->SetTitleSize(0.05,"Y");
    //hfit->GetXaxis()->SetRangeUser(0, 10000);
    //hfit->GetXaxis()->SetRange(0, 10000);
    hfit->Draw("hist");
    fit->SetLineWidth(2);
    fit->SetLineColor(kRed);
    fit->Draw("same");
    fit2->SetLineWidth(2);
    fit2->SetLineColor(kGreen+2);
    fit2->Draw("same");
    fit3->SetLineWidth(2);
    fit3->SetLineColor(kBlue);
    fit3->Draw("same");
    fit4->SetLineWidth(2);
    fit4->SetLineColor(kBlack);
    fit4->Draw("same");

    //TF1* testFit  = new TF1("Test Fit", gaufun, m4-s4, m4+s4, 3);
    //std::cout << m4-s4 << "  " << hfit->FindBin(m4-s4) << "  " << m4+s4 << " " <<hfit->FindBin(m4+s4) << std::endl;
    //testFit->FixParameter(0,fit4->GetParameter(0));
    //testFit->FixParameter(1,fit4->GetParameter(1));
    //testFit->FixParameter(2,fit4->GetParameter(2));
    //std::cout<<"Min Chi2:  "<<hfit->Chisquare(testFit,"R")<<std::endl;      
    //
    //for(int i = 0; i<10; i++){
    //  testFit->FixParameter(0, (testFit->GetParameter(0)+0.01*i));
    //  std::cout<<hfit->Chisquare(testFit,"R")<<std::endl;      
    //}
    
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

    TLatex* HPDTitle = new TLatex(0.93, 0.86, "HPD (Hybrid Photodiode)");
    HPDTitle->SetNDC();
    HPDTitle->SetTextFont(42);
    HPDTitle->SetTextAlign(32);

    TLatex* Muon = new TLatex(0.93, 0.81, "150 GeV Muons");
    Muon->SetNDC();
    Muon->SetTextFont(42);
    Muon->SetTextAlign(32);

    char chan [100];
    int iEta, iPhi, iDepth;

    if(sscanf (hfit->GetName(),"cluster_MIP_%d_%d", &iEta, &iPhi) == 2);
    else sscanf (hfit->GetName(),"clusterCal_MIP_%d_%d", &iEta, &iPhi);

    if(sscanf (hfit->GetName(),"iPhi2017_adc_%d", &iPhi) == 2) iEta = 19;
    else sscanf (hfit->GetName(),"iPhi2017Cal_adc_%d", &iPhi), iEta = 19;
    
    sprintf (chan, "Tower: %d,%d", iEta, iPhi);
    
    TLatex* channel = new TLatex(0.93, 0.86, chan);
    channel->SetNDC();
    channel->SetTextFont(42);
    channel->SetTextAlign(32);

    char sigma [100];
    float intsigma = fit4->GetParameter(1);
    sprintf (sigma,"#sigma: %0.3f", intsigma);

    TLatex* Sigma = new TLatex(0.93, 0.8, sigma);
    Sigma->SetNDC();
    Sigma->SetTextFont(42);
    Sigma->SetTextAlign(31);

    char mean [100];    
    float intmean = fit4->GetParameter(0);
    sprintf (mean,"mean: %0.3f", intmean);

    TLatex* Mean = new TLatex(0.93, 0.76, mean);
    Mean->SetNDC();
    Mean->SetTextFont(42);
    Mean->SetTextAlign(31);

    char res [100];
    float intres = intsigma/intmean;
    sprintf (res,"resolution:%0.3f",intres);

    TLatex* Res = new TLatex(0.93, 0.72, res);
    Res->SetNDC();
    Res->SetTextFont(42);
    Res->SetTextAlign(31);
    
    CMSPrelim1->Draw();
    testbeam->Draw();
    //HPDTitle->Draw();
    //Muon->Draw();
    channel->Draw();
    Sigma->Draw();
    Mean->Draw();
    Res->Draw();
    
    char oname[128];
    sprintf(oname, "%s_%i.pdf", hfit->GetName(), run);
    //sprintf(oname, "%s_HPDRuns_330to3475.pdf", hfit->GetName());
    c1.Print(oname);
}

int main()
{
    gROOT->SetStyle("Plain");    

    char fname[] = "analysis_run00%i_EMAP_AC05_01AUG2017_NOMINAL.root";
    
    for(int irun = 3137; irun < 3478; ++irun)
    //for(int irun = 3239; irun < 3240; ++irun)
    {
        char fn[128];
        sprintf(fn, fname, irun);
        TFile *fin = TFile::Open(fn);
    
        if(!fin) continue;
    
        TDirectory* din = (TDirectory*)fin->Get("hcalADCHists");
        
        TIter next(din->GetListOfKeys());
        TObject *obj;
        while((obj = next()))
        {
	  //if(strstr(obj->GetName(), "beam_adc") != nullptr)
	  //if(strstr(obj->GetName(), "cluster_MIP_19") != nullptr || strstr(obj->GetName(), "clusterCal_MIP_19") != nullptr )
	  if(strstr(obj->GetName(), "iPhi2017_adc_4") != nullptr || strstr(obj->GetName(), "iPhi2017Cal_adc_4") != nullptr ||strstr(obj->GetName(), "iPhi2017_adc_5") != nullptr || strstr(obj->GetName(), "iPhi2017Cal_adc_5") != nullptr)
            {
                TH1* htofit = (TH1*)din->Get(obj->GetName());
                if(htofit->GetMean() > 0.0)
                {
		    //printf("Hname: %s\n", htofit->GetName());
                    h = htofit;
                    fitSPEMIP(htofit, irun);
                }
            }
        }
    }
}
