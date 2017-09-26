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


//double crystalBall(double* x, double* p)
//{
//    //p[0]: normalization
//    //p[1]: gaussian mean
//    //p[2]: gaussian sigma
//    //p[3]: switchover point
//    //p[4]: power law slope
//    if((x[0] - p[1])/p[2] > -p[3])
//    {
//        return p[0]*TMath::Gaus(x[0], p[1], p[2]);
//    }
//    else
//    {
//        double A = pow(p[4]/fabs(p[3]), p[4]) * TMath::Exp(-pow(p[3], 2)/2);
//        double B = p[4]/fabs(p[3]) - fabs(p[3]);
//        return A*(B - pow((x[0]-p[1])/p[2], -p[4]));
//    }
//}

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

//Paper on SiPM pixel crosstalk model
//http://arxiv.org/pdf/1302.1455.pdf
//n is the number of neighbor cells in the crosstalk model 
const double n = 4;
//CPn are the constant combinatoric scale factors based on n
double CPn[6];

//More ugly globals to hold utility functions used in fits
TF1* funcs[8];
TF1* funcMIP;

TH1* h;

//Fit function for pedestal and photoelectron spectrum
double ppeFunc(double* x, double* p)
{
    /*
      parameters
      p[0] : pedestal amplitude
      p[1] : pedestal width
      p[2] : background amplidude
      p[3] : background width
      p[4] : overall PE amplitude 
      p[5] : poisson mean number of PE
      p[6] : PE peak spacing 
      p[7] : PE peak width
      p[8] : pixel cross-talk probability 
    */

    double s[7];
    double cp[7];
    double sc[7];
    //s defines the poisson factor for each PE peak based upon the poisson mean PE p[5]
    //in absense of cross-talk this would perfectly describe the PE peak amplitude ratios
    for(int i = 1; i < 6; ++i) s[i] = TMath::Poisson(i, p[5]);
    //cp represents the binomial scale factors for each PE peak based upon the crosstalk probability p[8]
    cp[1] = CPn[1]           * pow(1 - p[8], n);
    cp[2] = CPn[2] * p[8]    * pow(1 - p[8], 2*n - 1);
    cp[3] = CPn[3] * pow(p[8], 2) * pow(1 - p[8], 3*n - 2);
    cp[4] = CPn[4] * pow(p[8], 3) * pow(1 - p[8], 4*n - 3);
    cp[5] = CPn[5] * pow(p[8], 4) * pow(1 - p[8], 5*n - 4);
    //sc holds the convolution of of the poisson PE statistics factors s and binomial cross-talk factors cp 
    sc[1] = s[1] * cp[1];
    sc[2] = s[2] * pow(cp[1], 2) + s[1] * cp[2];
    sc[3] = s[3] * pow(cp[1], 3) + s[1] * cp[3] + s[2] * (2 * cp[1] * cp[2]);
    sc[4] = s[4] * pow(cp[1], 4) + s[1] * cp[4] + s[2] * (2 * cp[1] * cp[3] + pow(cp[2], 2))     + s[3] * (3 * pow(cp[1], 2) * cp[2]);
    sc[5] = s[5] * pow(cp[1], 5) + s[1] * cp[5] + s[2] * (2 * cp[1] * cp[4] + 2 * cp[2] * cp[3]) + s[3] * (3 * cp[1] * pow(cp[2], 2) + 3 * pow(cp[1], 2) * cp[3] + 2 * cp[1] * cp[4]) + s[4] * (3 * pow(cp[1], 3) * cp[2]);
    //Utility functions used parameterize the pedestal, PE peaks, and generic background 
    funcs[1]->SetParameters(p[0],    0.0, p[1]);
    funcs[2]->SetParameters(p[4],   p[6], p[7]);
    funcs[3]->SetParameters(p[4], 2*p[6], p[7]);
    funcs[4]->SetParameters(p[4], 3*p[6], p[7]);
    funcs[5]->SetParameters(p[4], 4*p[6], p[7]);
    funcs[6]->SetParameters(p[4], 5*p[6], p[7]);
    funcs[7]->SetParameters(p[2],    0.0, p[3]);
    //find bin edges 
    int iBin = h->FindBin(x[0]);
    double ll = h->GetBinLowEdge(iBin);
    double ul = ll + h->GetBinWidth(iBin);
    //calculate contribution to bin of interest from each sub-function
    double g1 =         funcs[1]->Integral(ll, ul)/(ul - ll);
    double g2 = sc[1] * funcs[2]->Integral(ll, ul)/(ul - ll);
    double g3 = sc[2] * funcs[3]->Integral(ll, ul)/(ul - ll);
    double g4 = sc[3] * funcs[4]->Integral(ll, ul)/(ul - ll);
    double g5 = sc[4] * funcs[5]->Integral(ll, ul)/(ul - ll);
    double g6 = sc[5] * funcs[6]->Integral(ll, ul)/(ul - ll);
    double g10 =        funcs[7]->Integral(ll, ul)/(ul - ll);
    //return bin value 
    return g1; //+ g2 + g3 + g4 + g5 + g6+ g10;
}

//Functional fit for pedestal plus landgauss MIP peak
double ppeMIPFunc(double* x, double* p)
{
    /*
      parameters
      p[9] : Total MIP peak area
      p[10] : Most probable value of Landau
      p[11] : Landau width parameter
      p[12] : gaussian width 
    */

    //find bin edges
    int iBin = h->FindBin(x[0]);;
    double ll = h->GetBinLowEdge(iBin);
    double ul = ll + h->GetBinWidth(iBin);
    //calculate pedestal, photoelectron, and background contribution
    double ppe = ppeFunc(x, p);
    //set mip parameters
    funcMIP->SetParameters(p[11], p[10], p[9], p[12]);
    //calculate MIP peak contribution
    //double mip = funcMIP->Integral(ll, ul)/(ul - ll);
    double mip = funcMIP->Eval(x[0]);
    //return bin value
    return ppe + mip;
}

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
    hfit->GetXaxis()->SetRangeUser(3, 3000);
    
    char ffname[128];

    //Fit Pedestal peak
    sprintf(ffname, "sff%s", hfit->GetName());
    TF1* fit = new TF1(ffname, ppeFunc, 0.0, 100000.0, 9);
    fit->SetParLimits(0,   0,   50000);
    fit->SetParLimits(1,   0,       5);
    fit->SetParLimits(2,   0,     0.1);
    fit->SetParLimits(3,   0,     0.1);
    fit->SetParLimits(4,   0,     0.1);
    fit->SetParLimits(5,   0,     0.1);
    fit->SetParLimits(6,   0,     0.1);
    fit->SetParLimits(7,   0,     0.1);
    fit->SetParLimits(8,   0,     0.1);
    fit->SetParameter(0,  40000);
    fit->FixParameter(1, 4.01301);
    fit->FixParameter(2,    0.0);
    fit->FixParameter(3,    0.0);
    fit->FixParameter(4,  0.001);
    fit->FixParameter(5,    0.0);
    fit->FixParameter(6, 0.0001);
    fit->FixParameter(7,    0.0);
    fit->FixParameter(8,    0.0);
    hfit->Fit(fit, "RNQL", "", 0, 5);
    
    //Fit MIP peak
    sprintf(ffname, "sff4%s", hfit->GetName());
    TF1* fit4 = new TF1(ffname, ppeMIPFunc, 0.0, 100000.0, 13);
    fit4->SetParLimits(0,    0, 50000);
    fit4->SetParLimits(1,    0,   0.1);
    fit4->SetParLimits(2,    0,   0.1);
    fit4->SetParLimits(3,    0,   0.1);
    fit4->SetParLimits(4,    0,   0.1);
    fit4->SetParLimits(5,    0,   0.1);
    fit4->SetParLimits(6,    0,   0.1);
    fit4->SetParLimits(7,    0,   0.1);
    fit4->SetParLimits(8,    0,   0.1);
    fit4->SetParLimits(9,    0, 10000);
    fit4->SetParLimits(10,   0,    10);
    fit4->SetParLimits(11,   0.1,     4);
    fit4->SetParLimits(12,   0,     5);
    fit4->FixParameter(0, fit->GetParameter(0));
    fit4->FixParameter(1, fit->GetParameter(1));
    fit4->FixParameter(2, fit->GetParameter(2)/2);
    fit4->FixParameter(3, fit->GetParameter(3)/2);
    fit4->FixParameter(4, fit->GetParameter(4));
    fit4->FixParameter(5, fit->GetParameter(5));
    fit4->FixParameter(6, fit->GetParameter(6));
    fit4->FixParameter(7, fit->GetParameter(7));
    fit4->FixParameter(8, fit->GetParameter(8));
    fit4->SetParameter(9, 9000);
    fit4->SetParameter(10,   5);
    fit4->SetParameter(11,   1);
    fit4->SetParameter(12,   4);
    hfit->Fit(fit4, "RNQL", "", 5, 30);

    //Fine Tunning Fit
    sprintf(ffname, "sff4%s", hfit->GetName());
    TF1* fit5 = new TF1(ffname, ppeMIPFunc, 0.0, 100000.0, 13);
    double down = 0.7;
    double up   = 1.3;
    fit5->SetParLimits(0, down*fit4->GetParameter(0),  up*fit4->GetParameter(0));
    fit5->SetParLimits(1, down*fit4->GetParameter(1),  up*fit4->GetParameter(1));
    fit5->SetParLimits(2, down*fit4->GetParameter(2),  up*fit4->GetParameter(2));
    fit5->SetParLimits(3, down*fit4->GetParameter(3),  up*fit4->GetParameter(3));
    fit5->SetParLimits(4, down*fit4->GetParameter(4),  up*fit4->GetParameter(4));
    fit5->SetParLimits(5, down*fit4->GetParameter(5),  up*fit4->GetParameter(5));
    fit5->SetParLimits(6, down*fit4->GetParameter(6),  up*fit4->GetParameter(6));
    fit5->SetParLimits(7, down*fit4->GetParameter(7),  up*fit4->GetParameter(7));
    fit5->SetParLimits(8, down*fit4->GetParameter(8),  up*fit4->GetParameter(8));
    fit5->SetParLimits(9, down*fit4->GetParameter(9),  up*fit4->GetParameter(9));
    fit5->SetParLimits(10,down*fit4->GetParameter(10), up*fit4->GetParameter(10));
    fit5->SetParLimits(11,down*fit4->GetParameter(11), up*fit4->GetParameter(11));
    fit5->SetParLimits(12,down*fit4->GetParameter(12), up*fit4->GetParameter(12));
    fit5->SetParameter(0,  fit4->GetParameter(0));
    fit5->SetParameter(1,  fit4->GetParameter(1));
    fit5->FixParameter(2,  fit4->GetParameter(2));
    fit5->FixParameter(3,  fit4->GetParameter(3));
    fit5->SetParameter(4,  fit4->GetParameter(4));
    fit5->SetParameter(5,  fit4->GetParameter(5));
    fit5->SetParameter(6,  fit4->GetParameter(6));
    fit5->SetParameter(7,  fit4->GetParameter(7));
    fit5->SetParameter(8,  fit4->GetParameter(8));
    fit5->SetParameter(9,  fit4->GetParameter(9));
    fit5->SetParameter(10, fit4->GetParameter(10));
    fit5->SetParameter(11, fit4->GetParameter(11));
    fit5->SetParameter(12, fit4->GetParameter(12));
    hfit->Fit(fit5, "RNQL", "", 5, 30);
    
    printf("Chi^2: %10.4f, P0: %10.4f, P1: %10.4f, P2: %10.4f, P3: %10.4f, P4: %10.4f, P5: %10.4f, P6: %10.4f, P7: %10.4f, P8: %10.4f, P9: %10.4f, P10: %10.4f, P11: %10.4f, P12: %10.4f\n", fit5->GetChisquare(), fit5->GetParameter(0), fit5->GetParameter(1), fit5->GetParameter(2), fit5->GetParameter(3), fit5->GetParameter(4), fit5->GetParameter(5), fit5->GetParameter(6), fit5->GetParameter(7), fit5->GetParameter(8), fit5->GetParameter(9), fit5->GetParameter(10), fit5->GetParameter(11), fit5->GetParameter(12));
    std::cout<<"ans: "<<hfit->GetName()<<"  "<<fit5->GetParameter(10)<<"  "<<fit5->GetParError(10)<<std::endl;
    hfit->Draw("hist");
    hfit->GetYaxis()->SetTitle("Events /fC");
    hfit->GetXaxis()->SetTitle("Charge [fC]");
    hfit->SetTitleOffset(1,"X");
    hfit->SetTitleOffset(1.2,"Y");
    hfit->SetTitleSize(0.05,"X");
    hfit->SetTitleSize(0.05,"Y");
    //fit3->SetLineColor(kRed);
    //fit3->SetLineWidth(2);
    fit5->SetLineWidth(2);

    TF1* fit6 = new TF1("BackGround", ppeFunc, 0.0, 100000.0, 9);
    fit6->FixParameter(0,  fit5->GetParameter(0));
    fit6->FixParameter(1,  fit5->GetParameter(1));
    fit6->FixParameter(2,  fit5->GetParameter(2));
    fit6->FixParameter(3,  fit5->GetParameter(3));
    fit6->FixParameter(4,  fit5->GetParameter(4));
    fit6->FixParameter(5,  fit5->GetParameter(5));
    fit6->FixParameter(6,  fit5->GetParameter(6));
    fit6->FixParameter(7,  fit5->GetParameter(7));
    fit6->FixParameter(8,  fit5->GetParameter(8));
    TF1* fit7 = new TF1("Lang", ppeMIPFunc, 0.0, 100000.0, 13);
    fit7->FixParameter(0,  0.0);
    fit7->FixParameter(1,  0.0);
    fit7->FixParameter(2,  0.0);
    fit7->FixParameter(3,  0.0);
    fit7->FixParameter(4,  0.0);
    fit7->FixParameter(5,  0.0);
    fit7->FixParameter(6,  0.0);
    fit7->FixParameter(7,  0.0);
    fit7->FixParameter(8,  0.0);
    fit7->FixParameter(9,  fit5->GetParameter(9));
    fit7->FixParameter(10, fit5->GetParameter(10));
    fit7->FixParameter(11, fit5->GetParameter(11));
    fit7->FixParameter(12, fit5->GetParameter(12));
    fit6->SetLineColor(kRed);
    fit7->SetLineColor(kGreen+2);
    fit7->Draw("same");
    fit6->Draw("same");
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
    sscanf (hfit->GetName(),"tower_adc_%d_%d", &iEta, &iPhi);
    sprintf (chan, "Tower: %d,%d", iEta, iPhi);
    
    TLatex* channel = new TLatex(0.93, 0.76, chan);
    channel->SetNDC();
    channel->SetTextFont(42);
    channel->SetTextAlign(32);

        char mpv [100];
    float intmpv = fit5->GetParameter(10);
    sprintf (mpv,"MPV: %0.3f", intmpv);

    TLatex* MPV = new TLatex(0.93, 0.7, mpv);
    MPV->SetNDC();
    MPV->SetTextFont(42);
    MPV->SetTextAlign(31);

    CMSPrelim1->Draw();
    testbeam->Draw();
    HPDTitle->Draw();
    Muon->Draw();
    channel->Draw();
    MPV->Draw();
    
    char oname[128];
    sprintf(oname, "%s_%i.pdf", hfit->GetName(), run);
    //sprintf(oname, "%s_HPDRuns_330to3475.pdf", hfit->GetName());
    c1.Print(oname);
}


int main()
{
    gROOT->SetStyle("Plain");    

    //calculation of combinatoric prefactors for cross-talk 
    CPn[1] = 1;
    CPn[2] = n;
    CPn[3] = 0.5*n*(3*n - 1);
    CPn[4] = (1.0/3.0)*n*(8*pow(n, 2) - 6*n + 1);
    CPn[5] = 0.25*n*((125.0/6.0)*pow(n, 3) - 25*pow(n, 2) + (55.0/6.0)*n - 1);
    
    //configure utility functions for fits
    funcs[1] = new TF1("ped", "gaus");
    funcs[2] = new TF1("pe1", "gaus");
    funcs[3] = new TF1("pe2", "gaus");
    funcs[4] = new TF1("pe3", "gaus");
    funcs[5] = new TF1("pe4", "gaus");
    funcs[6] = new TF1("pe5", "gaus");
    funcs[7] = new TF1("bg",  "landau");
    funcMIP = new TF1("mip",  langaufun, 0, 400000, 4);

    char fname[] = "data/analysis_run00%i_EMAP_AC05_01AUG2017_NOMINAL.root";
    
    for(int irun = 3030; irun < 3476; ++irun)
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
	  if(strstr(obj->GetName(), "tower_adc") != nullptr && (strstr(obj->GetName(), "_3") != nullptr || strstr(obj->GetName(), "_4") != nullptr))	      
            {
                TH1* htofit = (TH1*)din->Get(obj->GetName());
                if(htofit->GetMean() > 4.0)
                {
                    printf("Hname: %s\n", htofit->GetName());
                    h = htofit;
                    fitSPEMIP(htofit, irun);
                }
            }
        }
    }
}
