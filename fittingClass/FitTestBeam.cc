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
#include "FitTestBeam.h"

double FitTestBeam::crystalBall(double* x, double* p)
{
  //p[0]: normalization
  //p[1]: gaussian mean
  //p[2]: gaussian sigma
  //p[3]: switchover point
  //p[4]: power law slope
  if((x[0] - p[1])/p[2] > -p[3])
    {
      return p[0]*TMath::Gaus(x[0], p[1], p[2]);
    }
  else
    {
      double A = pow(p[4]/fabs(p[3]), p[4]) * TMath::Exp(-pow(p[3], 2)/2);
      double B = p[4]/fabs(p[3]) - fabs(p[3]);
      return A*(B - pow((x[0]-p[1])/p[2], -p[4]));
    }
}

// landau-gaussian convolution using numerical integration
Double_t FitTestBeam::langaufun(Double_t* x, Double_t* par)
{  
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
  for(i=1.0; i<=np/2; i++){
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);

    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }

  return (par[2] * step * sum * invsq2pi / par[3]);
}

double FitTestBeam::ppeFunc(double* x, double* p)
{
  /*
    parameters
    p[0] : pedestal amplitude
    p[1] : pedestal width
    p[2] : Overall Shift (avg. pedestal)
    p[3] : background amplidude
    p[4] : background width
    p[5] : overall PE amplitude 
    p[6] : poisson mean number of PE
    p[7] : PE peak spacing 
    p[8] : PE peak width
    p[9] : pixel cross-talk probability 
  */

  for(int i = 1; i < nPeaks_+1; ++i) s[i] = TMath::Poisson(i, p[6]);
  //cp represents the binomial scale factors for each PE peak based upon the crosstalk probability p[9]
  cp[1]  = CPn[1]                  * pow(1 - p[9], n);
  cp[2]  = CPn[2]  * p[9]          * pow(1 - p[9], 2*n  -  1);
  cp[3]  = CPn[3]  * pow(p[9], 2)  * pow(1 - p[9], 3*n  -  2);
  cp[4]  = CPn[4]  * pow(p[9], 3)  * pow(1 - p[9], 4*n  -  3);
  cp[5]  = CPn[5]  * pow(p[9], 4)  * pow(1 - p[9], 5*n  -  4);
  cp[6]  = CPn[6]  * pow(p[9], 5)  * pow(1 - p[9], 6*n  -  5);
  cp[7]  = CPn[7]  * pow(p[9], 6)  * pow(1 - p[9], 7*n  -  6);
  cp[8]  = CPn[8]  * pow(p[9], 7)  * pow(1 - p[9], 8*n  -  7);
  cp[9]  = CPn[9]  * pow(p[9], 8)  * pow(1 - p[9], 9*n  -  8);
  cp[10] = CPn[10] * pow(p[9], 9)  * pow(1 - p[9], 10*n -  9);
  cp[11] = CPn[11] * pow(p[9], 10) * pow(1 - p[9], 11*n - 10);
  cp[12] = CPn[12] * pow(p[9], 11) * pow(1 - p[9], 12*n - 11);
  cp[13] = CPn[13] * pow(p[9], 12) * pow(1 - p[9], 13*n - 12);
  cp[14] = CPn[14] * pow(p[9], 13) * pow(1 - p[9], 14*n - 13);
  cp[15] = CPn[15] * pow(p[9], 14) * pow(1 - p[9], 15*n - 14);
    
  //sc holds the convolution of of the poisson PE statistics factors s and binomial cross-talk factors cp 
  sc[15] = 20*s[5]*cp[1]*cp[1]*cp[1]*cp[2]*cp[10] + 2*s[2]*cp[6]*cp[9] + 12*s[4]*cp[1]*cp[1]*cp[6]*cp[7] + 168*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[4]*cp[4] + 4*s[4]*cp[3]*cp[4]*cp[4]*cp[4] + 12*s[4]*cp[2]*cp[3]*cp[3]*cp[7] + 60*s[6]*cp[1]*cp[1]*cp[1]*cp[2]*cp[5]*cp[5] + 12*s[4]*cp[2]*cp[4]*cp[4]*cp[5] + 210*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[3]*cp[6] + 1*s[5]*cp[3]*cp[3]*cp[3]*cp[3]*cp[3] + 4*s[4]*cp[3]*cp[3]*cp[3]*cp[6] + 11*s[11]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[5] + 12*s[4]*cp[1]*cp[2]*cp[2]*cp[10] + 1*s[1]*cp[15] + 6*s[3]*cp[2]*cp[6]*cp[7] + 90*s[10]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[5] + 12*s[4]*cp[1]*cp[2]*cp[6]*cp[6] + 1*s[3]*cp[5]*cp[5]*cp[5] + 10*s[10]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[6] + 30*s[5]*cp[1]*cp[1]*cp[3]*cp[3]*cp[7] + 30*s[6]*cp[1]*cp[1]*cp[1]*cp[1]*cp[5]*cp[6] + 6*s[3]*cp[1]*cp[5]*cp[9] + 12*s[4]*cp[2]*cp[2]*cp[5]*cp[6] + 20*s[5]*cp[2]*cp[2]*cp[2]*cp[4]*cp[5] + 20*s[5]*cp[1]*cp[3]*cp[3]*cp[3]*cp[5] + 84*s[9]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3]*cp[3]*cp[3] + 24*s[4]*cp[1]*cp[3]*cp[4]*cp[7] + 180*s[6]*cp[1]*cp[1]*cp[2]*cp[2]*cp[4]*cp[5] + 4*s[4]*cp[1]*cp[1]*cp[1]*cp[12] + 30*s[6]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3]*cp[8] + 60*s[6]*cp[1]*cp[2]*cp[2]*cp[2]*cp[4]*cp[4] + 30*s[5]*cp[1]*cp[3]*cp[3]*cp[4]*cp[4] + 55*s[11]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3]*cp[3] + 42*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3]*cp[7] + 2*s[2]*cp[3]*cp[12] + 756*s[9]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[3]*cp[3] + 840*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[3]*cp[4] + 6*s[3]*cp[2]*cp[3]*cp[10] + 20*s[6]*cp[1]*cp[1]*cp[1]*cp[4]*cp[4]*cp[4] + 30*s[5]*cp[1]*cp[1]*cp[3]*cp[5]*cp[5] + 140*s[7]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[6] + 12*s[4]*cp[1]*cp[4]*cp[5]*cp[5] + 220*s[12]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2] + 60*s[6]*cp[1]*cp[1]*cp[3]*cp[3]*cp[3]*cp[4] + 30*s[5]*cp[2]*cp[2]*cp[3]*cp[3]*cp[5] + 60*s[6]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[8] + 504*s[9]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[3]*cp[4] + 105*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[7] + 105*s[7]*cp[1]*cp[2]*cp[2]*cp[2]*cp[2]*cp[3]*cp[3] + 12*s[12]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[4] + 36*s[9]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[4]*cp[4] + 105*s[7]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[2]*cp[5] + 30*s[6]*cp[1]*cp[1]*cp[1]*cp[1]*cp[4]*cp[7] + 110*s[11]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[4] + 24*s[4]*cp[1]*cp[2]*cp[3]*cp[9] + 6*s[3]*cp[1]*cp[6]*cp[8] + 30*s[5]*cp[1]*cp[1]*cp[4]*cp[4]*cp[5] + 3*s[3]*cp[2]*cp[2]*cp[11] + 24*s[4]*cp[1]*cp[3]*cp[5]*cp[6] + 60*s[6]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[7] + 3*s[3]*cp[3]*cp[6]*cp[6] + 6*s[3]*cp[1]*cp[2]*cp[12] + 6*s[3]*cp[1]*cp[3]*cp[11] + 12*s[4]*cp[1]*cp[4]*cp[4]*cp[6] + 20*s[5]*cp[1]*cp[2]*cp[2]*cp[2]*cp[8] + 120*s[5]*cp[1]*cp[2]*cp[3]*cp[4]*cp[5] + 12*s[4]*cp[1]*cp[1]*cp[2]*cp[11] + 90*s[10]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3]*cp[4] + 24*s[4]*cp[1]*cp[2]*cp[4]*cp[8] + 5*s[5]*cp[1]*cp[1]*cp[1]*cp[1]*cp[11] + 3*s[3]*cp[3]*cp[3]*cp[9] + 210*s[7]*cp[1]*cp[1]*cp[2]*cp[2]*cp[3]*cp[3]*cp[3] + 105*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3]*cp[3]*cp[5] + 120*s[6]*cp[1]*cp[1]*cp[1]*cp[2]*cp[4]*cp[6] + 42*s[7]*cp[1]*cp[2]*cp[2]*cp[2]*cp[2]*cp[2]*cp[4] + 42*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[4]*cp[6] + 3*s[3]*cp[4]*cp[4]*cp[7] + 6*s[3]*cp[4]*cp[5]*cp[6] + 252*s[9]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[5] + 2*s[2]*cp[4]*cp[11] + 120*s[6]*cp[1]*cp[2]*cp[2]*cp[2]*cp[3]*cp[5] + 2*s[2]*cp[5]*cp[10] + 630*s[9]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[2]*cp[3] + 7*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[9] + 2*s[2]*cp[7]*cp[8] + 280*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[5] + 5*s[5]*cp[2]*cp[2]*cp[2]*cp[2]*cp[7] + 180*s[6]*cp[1]*cp[1]*cp[2]*cp[3]*cp[4]*cp[4] + 12*s[4]*cp[1]*cp[1]*cp[3]*cp[10] + 280*s[8]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[2]*cp[4] + 105*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3]*cp[4]*cp[4] + 12*s[4]*cp[1]*cp[1]*cp[5]*cp[8] + 132*s[12]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[3] + 336*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[3]*cp[5] + 6*s[3]*cp[2]*cp[5]*cp[8] + 9*s[9]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[7] + 60*s[5]*cp[1]*cp[2]*cp[2]*cp[3]*cp[7] + 24*s[4]*cp[2]*cp[3]*cp[4]*cp[6] + 20*s[6]*cp[2]*cp[2]*cp[2]*cp[3]*cp[3]*cp[3] + 180*s[6]*cp[1]*cp[1]*cp[2]*cp[3]*cp[3]*cp[5] + 24*s[4]*cp[1]*cp[2]*cp[5]*cp[7] + 60*s[5]*cp[1]*cp[2]*cp[3]*cp[3]*cp[6] + 2*s[2]*cp[1]*cp[14] + 280*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[3]*cp[3]*cp[3] + 840*s[10]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[3] + 2*s[2]*cp[2]*cp[13] + 30*s[5]*cp[1]*cp[1]*cp[2]*cp[2]*cp[9] + 3*s[3]*cp[1]*cp[7]*cp[7] + 21*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[5]*cp[5] + 30*s[6]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[9] + 495*s[11]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[3] + 420*s[7]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[3]*cp[4] + 84*s[9]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[2]*cp[2]*cp[2] + 78*s[13]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2] + 12*s[4]*cp[2]*cp[2]*cp[4]*cp[7] + 60*s[5]*cp[1]*cp[1]*cp[2]*cp[4]*cp[7] + 30*s[6]*cp[1]*cp[2]*cp[3]*cp[3]*cp[3]*cp[3] + 360*s[10]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[4] + 6*s[3]*cp[2]*cp[4]*cp[9] + 12*s[4]*cp[1]*cp[1]*cp[4]*cp[9] + 168*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3]*cp[3]*cp[4] + 168*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[6] + 210*s[7]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[4]*cp[4] + 20*s[5]*cp[2]*cp[2]*cp[2]*cp[3]*cp[6] + 1*s[15]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1] + 30*s[6]*cp[1]*cp[2]*cp[2]*cp[2]*cp[2]*cp[6] + 72*s[9]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[6] + 30*s[6]*cp[2]*cp[2]*cp[2]*cp[2]*cp[3]*cp[4] + 13*s[13]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3] + 42*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[8] + 12*s[4]*cp[2]*cp[3]*cp[5]*cp[5] + 8*s[8]*cp[1]*cp[2]*cp[2]*cp[2]*cp[2]*cp[2]*cp[2]*cp[2] + 12*s[4]*cp[2]*cp[2]*cp[3]*cp[8] + 360*s[10]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[3]*cp[3] + 120*s[6]*cp[1]*cp[1]*cp[1]*cp[2]*cp[3]*cp[7] + 60*s[6]*cp[1]*cp[1]*cp[1]*cp[3]*cp[3]*cp[6] + 560*s[8]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[3]*cp[3] + 20*s[5]*cp[1]*cp[2]*cp[4]*cp[4]*cp[4] + 20*s[5]*cp[1]*cp[1]*cp[1]*cp[4]*cp[8] + 330*s[11]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[2] + 60*s[5]*cp[1]*cp[1]*cp[3]*cp[4]*cp[6] + 504*s[9]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[4] + 12*s[4]*cp[3]*cp[3]*cp[4]*cp[5] + 168*s[8]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[2]*cp[2]*cp[3] + 56*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[4]*cp[5] + 60*s[5]*cp[1]*cp[2]*cp[2]*cp[4]*cp[6] + 56*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[7] + 72*s[9]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3]*cp[5] + 10*s[5]*cp[1]*cp[1]*cp[1]*cp[6]*cp[6] + 3*s[3]*cp[1]*cp[1]*cp[13] + 180*s[6]*cp[1]*cp[1]*cp[2]*cp[2]*cp[3]*cp[6] + 8*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[8] + 6*s[3]*cp[3]*cp[4]*cp[8] + 252*s[10]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[2]*cp[2] + 30*s[5]*cp[1]*cp[2]*cp[2]*cp[5]*cp[5] + 420*s[7]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[3]*cp[5] + 20*s[5]*cp[1]*cp[1]*cp[1]*cp[3]*cp[9] + 6*s[6]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[10] + 20*s[5]*cp[2]*cp[3]*cp[3]*cp[3]*cp[4] + 56*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3]*cp[6] + 35*s[7]*cp[1]*cp[1]*cp[1]*cp[3]*cp[3]*cp[3]*cp[3] + 14*s[14]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2] + 210*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[4]*cp[5] + 6*s[3]*cp[1]*cp[4]*cp[10] + 30*s[5]*cp[2]*cp[2]*cp[3]*cp[4]*cp[4] + 4*s[4]*cp[2]*cp[2]*cp[2]*cp[9] + 120*s[6]*cp[1]*cp[1]*cp[1]*cp[3]*cp[4]*cp[5] + 6*s[6]*cp[2]*cp[2]*cp[2]*cp[2]*cp[2]*cp[5] + 12*s[4]*cp[1]*cp[3]*cp[3]*cp[8] + 6*s[3]*cp[3]*cp[5]*cp[7] + 60*s[5]*cp[1]*cp[1]*cp[2]*cp[5]*cp[6] + 420*s[7]*cp[1]*cp[1]*cp[1]*cp[2]*cp[3]*cp[3]*cp[4] + 7*s[7]*cp[2]*cp[2]*cp[2]*cp[2]*cp[2]*cp[2]*cp[3] + 60*s[5]*cp[1]*cp[1]*cp[2]*cp[3]*cp[8] + 180*s[6]*cp[1]*cp[2]*cp[2]*cp[3]*cp[3]*cp[4] + 20*s[5]*cp[1]*cp[1]*cp[1]*cp[5]*cp[7];
  sc[14] = 8*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[7] + 110*s[11]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[3] + 42*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[7] + 12*s[4]*cp[1]*cp[1]*cp[4]*cp[8] + 12*s[4]*cp[1]*cp[4]*cp[4]*cp[5] + 4*s[4]*cp[3]*cp[3]*cp[3]*cp[5] + 20*s[5]*cp[1]*cp[1]*cp[1]*cp[3]*cp[8] + 45*s[10]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3]*cp[3] + 420*s[7]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[3]*cp[4] + 12*s[4]*cp[2]*cp[2]*cp[3]*cp[7] + 24*s[4]*cp[1]*cp[2]*cp[4]*cp[7] + 90*s[6]*cp[1]*cp[1]*cp[2]*cp[2]*cp[4]*cp[4] + 6*s[4]*cp[2]*cp[2]*cp[5]*cp[5] + 2*s[2]*cp[4]*cp[10] + 20*s[5]*cp[2]*cp[2]*cp[2]*cp[3]*cp[5] + 1*s[7]*cp[2]*cp[2]*cp[2]*cp[2]*cp[2]*cp[2]*cp[2] + 6*s[4]*cp[3]*cp[3]*cp[4]*cp[4] + 105*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3]*cp[3]*cp[4] + 168*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[5] + 6*s[3]*cp[2]*cp[5]*cp[7] + 2*s[2]*cp[6]*cp[8] + 105*s[7]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[2]*cp[4] + 6*s[3]*cp[2]*cp[4]*cp[8] + 60*s[6]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[7] + 210*s[10]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[2] + 72*s[9]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3]*cp[4] + 6*s[6]*cp[2]*cp[2]*cp[2]*cp[2]*cp[2]*cp[4] + 20*s[5]*cp[1]*cp[2]*cp[2]*cp[2]*cp[7] + 30*s[5]*cp[1]*cp[1]*cp[3]*cp[3]*cp[6] + 7*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[8] + 6*s[3]*cp[3]*cp[4]*cp[7] + 4*s[4]*cp[2]*cp[2]*cp[2]*cp[8] + 60*s[5]*cp[1]*cp[2]*cp[2]*cp[4]*cp[5] + 140*s[7]*cp[1]*cp[1]*cp[1]*cp[2]*cp[3]*cp[3]*cp[3] + 2*s[2]*cp[2]*cp[12] + 15*s[6]*cp[1]*cp[1]*cp[3]*cp[3]*cp[3]*cp[3] + 210*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[3]*cp[5] + 66*s[12]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2] + 12*s[4]*cp[1]*cp[3]*cp[5]*cp[5] + 24*s[4]*cp[2]*cp[3]*cp[4]*cp[5] + 15*s[6]*cp[2]*cp[2]*cp[2]*cp[2]*cp[3]*cp[3] + 140*s[7]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[5] + 12*s[4]*cp[2]*cp[3]*cp[3]*cp[6] + 6*s[4]*cp[1]*cp[1]*cp[6]*cp[6] + 42*s[7]*cp[1]*cp[2]*cp[2]*cp[2]*cp[2]*cp[2]*cp[3] + 60*s[6]*cp[1]*cp[1]*cp[1]*cp[3]*cp[4]*cp[4] + 1*s[2]*cp[7]*cp[7] + 4*s[4]*cp[2]*cp[4]*cp[4]*cp[4] + 10*s[10]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[5] + 12*s[4]*cp[1]*cp[2]*cp[2]*cp[9] + 72*s[9]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[5] + 30*s[5]*cp[2]*cp[2]*cp[3]*cp[3]*cp[4] + 3*s[3]*cp[1]*cp[1]*cp[12] + 120*s[6]*cp[1]*cp[2]*cp[2]*cp[2]*cp[3]*cp[4] + 30*s[5]*cp[1]*cp[1]*cp[2]*cp[5]*cp[5] + 28*s[8]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[2]*cp[2]*cp[2] + 28*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[4]*cp[4] + 9*s[9]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[6] + 56*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[6] + 30*s[6]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3]*cp[7] + 60*s[5]*cp[1]*cp[2]*cp[3]*cp[3]*cp[5] + 60*s[6]*cp[1]*cp[2]*cp[2]*cp[3]*cp[3]*cp[3] + 24*s[4]*cp[1]*cp[2]*cp[5]*cp[6] + 10*s[5]*cp[1]*cp[1]*cp[4]*cp[4]*cp[4] + 336*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[3]*cp[4] + 105*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[6] + 24*s[4]*cp[1]*cp[3]*cp[4]*cp[6] + 6*s[3]*cp[1]*cp[4]*cp[9] + 6*s[3]*cp[3]*cp[5]*cp[6] + 1*s[14]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1] + 20*s[5]*cp[1]*cp[3]*cp[3]*cp[3]*cp[4] + 13*s[13]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2] + 180*s[6]*cp[1]*cp[1]*cp[2]*cp[3]*cp[3]*cp[4] + 5*s[5]*cp[1]*cp[1]*cp[1]*cp[1]*cp[10] + 60*s[5]*cp[1]*cp[2]*cp[3]*cp[4]*cp[4] + 6*s[3]*cp[1]*cp[2]*cp[11] + 30*s[6]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[8] + 165*s[11]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2] + 20*s[5]*cp[1]*cp[1]*cp[1]*cp[2]*cp[9] + 3*s[3]*cp[2]*cp[6]*cp[6] + 60*s[5]*cp[1]*cp[1]*cp[2]*cp[4]*cp[6] + 252*s[9]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[4] + 56*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3]*cp[5] + 5*s[5]*cp[2]*cp[2]*cp[2]*cp[2]*cp[6] + 24*s[4]*cp[1]*cp[2]*cp[3]*cp[8] + 12*s[4]*cp[1]*cp[3]*cp[3]*cp[7] + 6*s[3]*cp[1]*cp[5]*cp[8] + 60*s[6]*cp[1]*cp[1]*cp[1]*cp[3]*cp[3]*cp[5] + 280*s[8]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[2]*cp[3] + 360*s[10]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[3] + 2*s[2]*cp[3]*cp[11] + 20*s[5]*cp[1]*cp[1]*cp[1]*cp[4]*cp[7] + 42*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[4]*cp[5] + 11*s[11]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[4] + 105*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[4]*cp[4] + 6*s[3]*cp[2]*cp[3]*cp[9] + 60*s[5]*cp[1]*cp[2]*cp[2]*cp[3]*cp[6] + 12*s[4]*cp[1]*cp[1]*cp[2]*cp[10] + 120*s[6]*cp[1]*cp[1]*cp[1]*cp[2]*cp[3]*cp[6] + 252*s[9]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[3]*cp[3] + 30*s[5]*cp[1]*cp[1]*cp[2]*cp[2]*cp[8] + 2*s[2]*cp[5]*cp[9] + 504*s[9]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[3] + 120*s[6]*cp[1]*cp[1]*cp[1]*cp[2]*cp[4]*cp[5] + 6*s[6]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[9] + 56*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3]*cp[3]*cp[3] + 30*s[6]*cp[1]*cp[2]*cp[2]*cp[2]*cp[2]*cp[5] + 10*s[5]*cp[2]*cp[2]*cp[2]*cp[4]*cp[4] + 42*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3]*cp[6] + 90*s[10]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[4] + 280*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[4] + 4*s[4]*cp[1]*cp[1]*cp[1]*cp[11] + 3*s[3]*cp[2]*cp[2]*cp[10] + 20*s[5]*cp[1]*cp[1]*cp[1]*cp[5]*cp[6] + 210*s[7]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[3]*cp[3] + 60*s[6]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[6] + 6*s[3]*cp[1]*cp[3]*cp[10] + 15*s[6]*cp[1]*cp[1]*cp[1]*cp[1]*cp[5]*cp[5] + 420*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[3]*cp[3] + 12*s[12]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3] + 2*s[2]*cp[1]*cp[13] + 60*s[5]*cp[1]*cp[1]*cp[3]*cp[4]*cp[5] + 12*s[4]*cp[1]*cp[1]*cp[3]*cp[9] + 3*s[3]*cp[3]*cp[3]*cp[8] + 5*s[5]*cp[2]*cp[3]*cp[3]*cp[3]*cp[3] + 3*s[3]*cp[4]*cp[5]*cp[5] + 6*s[3]*cp[1]*cp[6]*cp[7] + 180*s[6]*cp[1]*cp[1]*cp[2]*cp[2]*cp[3]*cp[5] + 126*s[9]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[2]*cp[2] + 1*s[1]*cp[14] + 3*s[3]*cp[4]*cp[4]*cp[6] + 12*s[4]*cp[1]*cp[1]*cp[5]*cp[7] + 60*s[5]*cp[1]*cp[1]*cp[2]*cp[3]*cp[7] + 12*s[4]*cp[2]*cp[2]*cp[4]*cp[6] + 30*s[6]*cp[1]*cp[1]*cp[1]*cp[1]*cp[4]*cp[6];
  sc[13] = 280*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[3] + 6*s[3]*cp[1]*cp[2]*cp[10] + 12*s[4]*cp[1]*cp[1]*cp[2]*cp[9] + 56*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3]*cp[4] + 20*s[5]*cp[1]*cp[1]*cp[1]*cp[3]*cp[7] + 2*s[2]*cp[5]*cp[8] + 3*s[3]*cp[2]*cp[2]*cp[9] + 2*s[2]*cp[6]*cp[7] + 20*s[5]*cp[1]*cp[2]*cp[2]*cp[2]*cp[6] + 56*s[8]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[2]*cp[2] + 5*s[5]*cp[1]*cp[1]*cp[1]*cp[1]*cp[9] + 42*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3]*cp[5] + 30*s[5]*cp[1]*cp[2]*cp[2]*cp[4]*cp[4] + 72*s[9]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[4] + 30*s[6]*cp[1]*cp[1]*cp[1]*cp[1]*cp[4]*cp[5] + 36*s[9]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3]*cp[3] + 6*s[3]*cp[2]*cp[3]*cp[8] + 105*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[5] + 2*s[2]*cp[1]*cp[12] + 6*s[3]*cp[1]*cp[5]*cp[7] + 2*s[2]*cp[2]*cp[11] + 12*s[4]*cp[2]*cp[3]*cp[4]*cp[4] + 12*s[4]*cp[1]*cp[2]*cp[5]*cp[5] + 42*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[6] + 60*s[5]*cp[1]*cp[2]*cp[2]*cp[3]*cp[5] + 24*s[4]*cp[1]*cp[2]*cp[3]*cp[7] + 30*s[6]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[7] + 11*s[11]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3] + 60*s[5]*cp[1]*cp[2]*cp[3]*cp[3]*cp[4] + 12*s[4]*cp[1]*cp[3]*cp[3]*cp[6] + 12*s[4]*cp[1]*cp[1]*cp[4]*cp[7] + 12*s[4]*cp[2]*cp[2]*cp[4]*cp[5] + 10*s[10]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[4] + 60*s[6]*cp[1]*cp[2]*cp[2]*cp[2]*cp[3]*cp[3] + 12*s[4]*cp[2]*cp[3]*cp[3]*cp[5] + 3*s[3]*cp[3]*cp[3]*cp[7] + 6*s[3]*cp[1]*cp[3]*cp[9] + 3*s[3]*cp[1]*cp[6]*cp[6] + 12*s[4]*cp[1]*cp[2]*cp[2]*cp[8] + 10*s[5]*cp[2]*cp[2]*cp[3]*cp[3]*cp[3] + 8*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[6] + 3*s[3]*cp[4]*cp[4]*cp[5] + 210*s[7]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[3]*cp[3] + 30*s[6]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3]*cp[6] + 12*s[4]*cp[1]*cp[1]*cp[5]*cp[6] + 60*s[6]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[6] + 120*s[10]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2] + 20*s[5]*cp[1]*cp[1]*cp[1]*cp[4]*cp[6] + 126*s[9]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[2] + 5*s[5]*cp[2]*cp[2]*cp[2]*cp[2]*cp[5] + 210*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[3]*cp[4] + 4*s[4]*cp[2]*cp[2]*cp[2]*cp[7] + 140*s[7]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[4] + 9*s[9]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[5] + 12*s[4]*cp[1]*cp[1]*cp[3]*cp[8] + 2*s[2]*cp[3]*cp[10] + 12*s[12]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2] + 7*s[7]*cp[1]*cp[2]*cp[2]*cp[2]*cp[2]*cp[2]*cp[2] + 30*s[5]*cp[1]*cp[1]*cp[2]*cp[2]*cp[7] + 168*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[3]*cp[3] + 5*s[5]*cp[1]*cp[3]*cp[3]*cp[3]*cp[3] + 6*s[6]*cp[2]*cp[2]*cp[2]*cp[2]*cp[2]*cp[3] + 6*s[3]*cp[2]*cp[5]*cp[6] + 60*s[6]*cp[1]*cp[1]*cp[2]*cp[3]*cp[3]*cp[3] + 120*s[6]*cp[1]*cp[1]*cp[1]*cp[2]*cp[3]*cp[5] + 60*s[6]*cp[1]*cp[1]*cp[1]*cp[3]*cp[3]*cp[4] + 105*s[7]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[2]*cp[3] + 4*s[4]*cp[1]*cp[4]*cp[4]*cp[4] + 4*s[4]*cp[3]*cp[3]*cp[3]*cp[4] + 10*s[5]*cp[1]*cp[1]*cp[1]*cp[5]*cp[5] + 20*s[5]*cp[1]*cp[1]*cp[1]*cp[2]*cp[8] + 1*s[13]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1] + 60*s[5]*cp[1]*cp[1]*cp[2]*cp[4]*cp[5] + 56*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[5] + 12*s[4]*cp[2]*cp[2]*cp[3]*cp[6] + 6*s[3]*cp[3]*cp[4]*cp[6] + 30*s[5]*cp[1]*cp[1]*cp[3]*cp[4]*cp[4] + 24*s[4]*cp[1]*cp[2]*cp[4]*cp[6] + 1*s[1]*cp[13] + 60*s[6]*cp[1]*cp[1]*cp[1]*cp[2]*cp[4]*cp[4] + 24*s[4]*cp[1]*cp[3]*cp[4]*cp[5] + 30*s[5]*cp[1]*cp[1]*cp[3]*cp[3]*cp[5] + 168*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[4] + 2*s[2]*cp[4]*cp[9] + 30*s[6]*cp[1]*cp[2]*cp[2]*cp[2]*cp[2]*cp[4] + 60*s[6]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[5] + 60*s[5]*cp[1]*cp[1]*cp[2]*cp[3]*cp[6] + 3*s[3]*cp[1]*cp[1]*cp[11] + 7*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[7] + 6*s[6]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[8] + 6*s[3]*cp[2]*cp[4]*cp[7] + 55*s[11]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2] + 21*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[4]*cp[4] + 35*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3]*cp[3]*cp[3] + 252*s[9]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[3] + 90*s[10]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[3] + 6*s[3]*cp[1]*cp[4]*cp[8] + 180*s[6]*cp[1]*cp[1]*cp[2]*cp[2]*cp[3]*cp[4] + 4*s[4]*cp[1]*cp[1]*cp[1]*cp[10] + 3*s[3]*cp[3]*cp[5]*cp[5] + 20*s[5]*cp[2]*cp[2]*cp[2]*cp[3]*cp[4];
  sc[12] = 42*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3]*cp[4] + 6*s[3]*cp[1]*cp[4]*cp[7] + 12*s[4]*cp[1]*cp[3]*cp[4]*cp[4] + 3*s[3]*cp[2]*cp[5]*cp[5] + 30*s[6]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3]*cp[5] + 168*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[3] + 1*s[2]*cp[6]*cp[6] + 6*s[3]*cp[1]*cp[2]*cp[9] + 72*s[9]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[3] + 60*s[6]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[4] + 6*s[6]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[7] + 21*s[7]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[2]*cp[2] + 6*s[3]*cp[2]*cp[3]*cp[7] + 42*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[5] + 60*s[5]*cp[1]*cp[2]*cp[2]*cp[3]*cp[4] + 3*s[3]*cp[3]*cp[3]*cp[6] + 20*s[6]*cp[1]*cp[1]*cp[1]*cp[3]*cp[3]*cp[3] + 10*s[5]*cp[2]*cp[2]*cp[2]*cp[3]*cp[3] + 60*s[6]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[5] + 1*s[3]*cp[4]*cp[4]*cp[4] + 1*s[4]*cp[3]*cp[3]*cp[3]*cp[3] + 105*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[3]*cp[3] + 12*s[4]*cp[2]*cp[2]*cp[3]*cp[5] + 2*s[2]*cp[1]*cp[11] + 10*s[10]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3] + 70*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[2] + 20*s[5]*cp[1]*cp[1]*cp[1]*cp[4]*cp[5] + 2*s[2]*cp[4]*cp[8] + 30*s[6]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[6] + 140*s[7]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[3] + 4*s[4]*cp[2]*cp[2]*cp[2]*cp[6] + 28*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3]*cp[3] + 15*s[6]*cp[1]*cp[1]*cp[1]*cp[1]*cp[4]*cp[4] + 45*s[10]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2] + 3*s[3]*cp[2]*cp[2]*cp[8] + 56*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[4] + 1*s[6]*cp[2]*cp[2]*cp[2]*cp[2]*cp[2]*cp[2] + 12*s[4]*cp[1]*cp[1]*cp[4]*cp[6] + 12*s[4]*cp[1]*cp[1]*cp[2]*cp[8] + 30*s[5]*cp[1]*cp[1]*cp[2]*cp[4]*cp[4] + 12*s[4]*cp[1]*cp[2]*cp[2]*cp[7] + 6*s[3]*cp[2]*cp[4]*cp[6] + 9*s[9]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[4] + 20*s[5]*cp[1]*cp[2]*cp[2]*cp[2]*cp[5] + 90*s[6]*cp[1]*cp[1]*cp[2]*cp[2]*cp[3]*cp[3] + 6*s[3]*cp[1]*cp[3]*cp[8] + 3*s[3]*cp[1]*cp[1]*cp[10] + 7*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[6] + 2*s[2]*cp[2]*cp[10] + 24*s[4]*cp[1]*cp[2]*cp[3]*cp[6] + 2*s[2]*cp[3]*cp[9] + 4*s[4]*cp[1]*cp[1]*cp[1]*cp[9] + 1*s[1]*cp[12] + 84*s[9]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2] + 12*s[4]*cp[2]*cp[3]*cp[3]*cp[4] + 6*s[4]*cp[2]*cp[2]*cp[4]*cp[4] + 20*s[5]*cp[1]*cp[1]*cp[1]*cp[2]*cp[7] + 1*s[12]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1] + 20*s[5]*cp[1]*cp[1]*cp[1]*cp[3]*cp[6] + 30*s[5]*cp[1]*cp[1]*cp[2]*cp[2]*cp[6] + 8*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[5] + 60*s[5]*cp[1]*cp[1]*cp[2]*cp[3]*cp[5] + 12*s[4]*cp[1]*cp[1]*cp[3]*cp[7] + 20*s[5]*cp[1]*cp[2]*cp[3]*cp[3]*cp[3] + 2*s[2]*cp[5]*cp[7] + 5*s[5]*cp[2]*cp[2]*cp[2]*cp[2]*cp[4] + 6*s[3]*cp[1]*cp[5]*cp[6] + 12*s[4]*cp[1]*cp[3]*cp[3]*cp[5] + 105*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[4] + 11*s[11]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2] + 120*s[6]*cp[1]*cp[1]*cp[1]*cp[2]*cp[3]*cp[4] + 30*s[6]*cp[1]*cp[2]*cp[2]*cp[2]*cp[2]*cp[3] + 30*s[5]*cp[1]*cp[1]*cp[3]*cp[3]*cp[4] + 6*s[3]*cp[3]*cp[4]*cp[5] + 24*s[4]*cp[1]*cp[2]*cp[4]*cp[5] + 5*s[5]*cp[1]*cp[1]*cp[1]*cp[1]*cp[8] + 6*s[4]*cp[1]*cp[1]*cp[5]*cp[5];    
  sc[11] = 2*s[2]*cp[4]*cp[7] + 21*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3]*cp[3] + 1*s[11]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1] + 2*s[2]*cp[5]*cp[6] + 12*s[4]*cp[2]*cp[2]*cp[3]*cp[4] + 42*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[4] + 6*s[3]*cp[1]*cp[2]*cp[8] + 10*s[5]*cp[1]*cp[1]*cp[3]*cp[3]*cp[3] + 1*s[1]*cp[11] + 20*s[5]*cp[1]*cp[2]*cp[2]*cp[2]*cp[4] + 60*s[5]*cp[1]*cp[1]*cp[2]*cp[3]*cp[4] + 8*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[4] + 30*s[5]*cp[1]*cp[1]*cp[2]*cp[2]*cp[5] + 3*s[3]*cp[1]*cp[1]*cp[9] + 60*s[6]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[3] + 2*s[2]*cp[2]*cp[9] + 10*s[10]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2] + 20*s[5]*cp[1]*cp[1]*cp[1]*cp[3]*cp[5] + 4*s[4]*cp[1]*cp[1]*cp[1]*cp[8] + 60*s[6]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[4] + 3*s[3]*cp[3]*cp[3]*cp[5] + 6*s[3]*cp[2]*cp[3]*cp[6] + 5*s[5]*cp[1]*cp[1]*cp[1]*cp[1]*cp[7] + 6*s[3]*cp[2]*cp[4]*cp[5] + 2*s[2]*cp[1]*cp[10] + 24*s[4]*cp[1]*cp[2]*cp[3]*cp[5] + 12*s[4]*cp[1]*cp[1]*cp[3]*cp[6] + 60*s[6]*cp[1]*cp[1]*cp[1]*cp[2]*cp[3]*cp[3] + 30*s[5]*cp[1]*cp[2]*cp[2]*cp[3]*cp[3] + 7*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[5] + 12*s[4]*cp[1]*cp[3]*cp[3]*cp[4] + 5*s[5]*cp[2]*cp[2]*cp[2]*cp[2]*cp[3] + 36*s[9]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2] + 12*s[4]*cp[1]*cp[1]*cp[4]*cp[5] + 105*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[3] + 12*s[4]*cp[1]*cp[1]*cp[2]*cp[7] + 12*s[4]*cp[1]*cp[2]*cp[4]*cp[4] + 6*s[6]*cp[1]*cp[2]*cp[2]*cp[2]*cp[2]*cp[2] + 3*s[3]*cp[1]*cp[5]*cp[5] + 56*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2] + 6*s[3]*cp[1]*cp[4]*cp[6] + 4*s[4]*cp[2]*cp[3]*cp[3]*cp[3] + 30*s[6]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[5] + 4*s[4]*cp[2]*cp[2]*cp[2]*cp[5] + 6*s[3]*cp[1]*cp[3]*cp[7] + 30*s[6]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3]*cp[4] + 6*s[6]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[6] + 56*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[3] + 9*s[9]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3] + 2*s[2]*cp[3]*cp[8] + 3*s[3]*cp[3]*cp[4]*cp[4] + 3*s[3]*cp[2]*cp[2]*cp[7] + 10*s[5]*cp[1]*cp[1]*cp[1]*cp[4]*cp[4] + 20*s[5]*cp[1]*cp[1]*cp[1]*cp[2]*cp[6] + 12*s[4]*cp[1]*cp[2]*cp[2]*cp[6] + 35*s[7]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[2];
  sc[10] = 6*s[3]*cp[2]*cp[3]*cp[5] + 4*s[4]*cp[1]*cp[1]*cp[1]*cp[7] + 9*s[9]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2] + 2*s[2]*cp[2]*cp[8] + 8*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3] + 6*s[3]*cp[1]*cp[4]*cp[5] + 1*s[10]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1] + 60*s[6]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[3] + 1*s[5]*cp[2]*cp[2]*cp[2]*cp[2]*cp[2] + 15*s[6]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2]*cp[2] + 2*s[2]*cp[3]*cp[7] + 28*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2] + 20*s[5]*cp[1]*cp[1]*cp[1]*cp[2]*cp[5] + 4*s[4]*cp[1]*cp[3]*cp[3]*cp[3] + 30*s[5]*cp[1]*cp[1]*cp[2]*cp[3]*cp[3] + 1*s[2]*cp[5]*cp[5] + 3*s[3]*cp[1]*cp[1]*cp[8] + 7*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[4] + 3*s[3]*cp[3]*cp[3]*cp[4] + 12*s[4]*cp[1]*cp[1]*cp[3]*cp[5] + 3*s[3]*cp[2]*cp[4]*cp[4] + 15*s[6]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3]*cp[3] + 30*s[6]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[4] + 6*s[4]*cp[2]*cp[2]*cp[3]*cp[3] + 20*s[5]*cp[1]*cp[2]*cp[2]*cp[2]*cp[3] + 12*s[4]*cp[1]*cp[2]*cp[2]*cp[5] + 2*s[2]*cp[1]*cp[9] + 4*s[4]*cp[2]*cp[2]*cp[2]*cp[4] + 6*s[6]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[5] + 2*s[2]*cp[4]*cp[6] + 35*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2] + 24*s[4]*cp[1]*cp[2]*cp[3]*cp[4] + 1*s[1]*cp[10] + 6*s[3]*cp[1]*cp[3]*cp[6] + 6*s[3]*cp[1]*cp[2]*cp[7] + 3*s[3]*cp[2]*cp[2]*cp[6] + 20*s[5]*cp[1]*cp[1]*cp[1]*cp[3]*cp[4] + 30*s[5]*cp[1]*cp[1]*cp[2]*cp[2]*cp[4] + 6*s[4]*cp[1]*cp[1]*cp[4]*cp[4] + 5*s[5]*cp[1]*cp[1]*cp[1]*cp[1]*cp[6] + 12*s[4]*cp[1]*cp[1]*cp[2]*cp[6] + 42*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[3];    
  sc[9] = 3*s[3]*cp[1]*cp[1]*cp[7] + 5*s[5]*cp[1]*cp[2]*cp[2]*cp[2]*cp[2] + 3*s[3]*cp[1]*cp[4]*cp[4] + 1*s[1]*cp[9] + 12*s[4]*cp[1]*cp[1]*cp[3]*cp[4] + 3*s[3]*cp[2]*cp[2]*cp[5] + 20*s[6]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2] + 1*s[3]*cp[3]*cp[3]*cp[3] + 8*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2] + 12*s[4]*cp[1]*cp[1]*cp[2]*cp[5] + 12*s[4]*cp[1]*cp[2]*cp[2]*cp[4] + 1*s[9]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1] + 30*s[6]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[3] + 12*s[4]*cp[1]*cp[2]*cp[3]*cp[3] + 2*s[2]*cp[3]*cp[6] + 6*s[6]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[4] + 5*s[5]*cp[1]*cp[1]*cp[1]*cp[1]*cp[5] + 2*s[2]*cp[4]*cp[5] + 7*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3] + 30*s[5]*cp[1]*cp[1]*cp[2]*cp[2]*cp[3] + 4*s[4]*cp[1]*cp[1]*cp[1]*cp[6] + 2*s[2]*cp[2]*cp[7] + 6*s[3]*cp[1]*cp[2]*cp[6] + 20*s[5]*cp[1]*cp[1]*cp[1]*cp[2]*cp[4] + 21*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2] + 6*s[3]*cp[1]*cp[3]*cp[5] + 2*s[2]*cp[1]*cp[8] + 10*s[5]*cp[1]*cp[1]*cp[1]*cp[3]*cp[3] + 6*s[3]*cp[2]*cp[3]*cp[4] + 4*s[4]*cp[2]*cp[2]*cp[2]*cp[3];    
  sc[8] = 20*s[5]*cp[1]*cp[1]*cp[1]*cp[2]*cp[3] + 3*s[3]*cp[2]*cp[3]*cp[3] + 3*s[3]*cp[1]*cp[1]*cp[6] + 3*s[3]*cp[2]*cp[2]*cp[4] + 15*s[6]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2] + 1*s[4]*cp[2]*cp[2]*cp[2]*cp[2] + 1*s[2]*cp[4]*cp[4] + 4*s[4]*cp[1]*cp[1]*cp[1]*cp[5] + 2*s[2]*cp[2]*cp[6] + 6*s[4]*cp[1]*cp[1]*cp[3]*cp[3] + 1*s[8]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1] + 1*s[1]*cp[8] + 7*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2] + 12*s[4]*cp[1]*cp[1]*cp[2]*cp[4] + 10*s[5]*cp[1]*cp[1]*cp[2]*cp[2]*cp[2] + 6*s[3]*cp[1]*cp[2]*cp[5] + 2*s[2]*cp[3]*cp[5] + 2*s[2]*cp[1]*cp[7] + 5*s[5]*cp[1]*cp[1]*cp[1]*cp[1]*cp[4] + 6*s[6]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3] + 6*s[3]*cp[1]*cp[3]*cp[4] + 12*s[4]*cp[1]*cp[2]*cp[2]*cp[3];    
  sc[7] = 1*s[7]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1] + 3*s[3]*cp[1]*cp[1]*cp[5] + 4*s[4]*cp[1]*cp[1]*cp[1]*cp[4] + 6*s[6]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2] + 5*s[5]*cp[1]*cp[1]*cp[1]*cp[1]*cp[3] + 4*s[4]*cp[1]*cp[2]*cp[2]*cp[2] + 12*s[4]*cp[1]*cp[1]*cp[2]*cp[3] + 10*s[5]*cp[1]*cp[1]*cp[1]*cp[2]*cp[2] + 2*s[2]*cp[1]*cp[6] + 6*s[3]*cp[1]*cp[2]*cp[4] + 1*s[1]*cp[7] + 2*s[2]*cp[2]*cp[5] + 2*s[2]*cp[3]*cp[4] + 3*s[3]*cp[2]*cp[2]*cp[3] + 3*s[3]*cp[1]*cp[3]*cp[3];    
  sc[6] = 1*s[3]*cp[2]*cp[2]*cp[2] + 1*s[2]*cp[3]*cp[3] + 6*s[4]*cp[1]*cp[1]*cp[2]*cp[2] + 5*s[5]*cp[1]*cp[1]*cp[1]*cp[1]*cp[2] + 6*s[3]*cp[1]*cp[2]*cp[3] + 2*s[2]*cp[1]*cp[5] + 1*s[1]*cp[6] + 2*s[2]*cp[2]*cp[4] + 1*s[6]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1] + 3*s[3]*cp[1]*cp[1]*cp[4] + 4*s[4]*cp[1]*cp[1]*cp[1]*cp[3];
  sc[5] = 3*s[3]*cp[1]*cp[1]*cp[3] + 2*s[2]*cp[1]*cp[4] + 3*s[3]*cp[1]*cp[2]*cp[2] + 1*s[1]*cp[5] + 2*s[2]*cp[2]*cp[3] + 4*s[4]*cp[1]*cp[1]*cp[1]*cp[2] + 1*s[5]*cp[1]*cp[1]*cp[1]*cp[1]*cp[1];    
  sc[4] = 3*s[3]*cp[1]*cp[1]*cp[2] + 2*s[2]*cp[1]*cp[3] + 1*s[4]*cp[1]*cp[1]*cp[1]*cp[1] + 1*s[2]*cp[2]*cp[2] + 1*s[1]*cp[4];    
  sc[3] = 2*s[2]*cp[1]*cp[2] + 1*s[3]*cp[1]*cp[1]*cp[1] + 1*s[1]*cp[3];    
  sc[2] = 1*s[1]*cp[2] + 1*s[2]*cp[1]*cp[1];    
  sc[1] = 1*s[1]*cp[1];
    
  //Utility functions used parameterize the pedestal, PE peaks, and generic background 
  //Real number of photoelectron peaks
  int numPEP = 15;
  int i;
  funcs[0]->SetParameters(p[0],    0.0, p[1]);
  for(i = 1; i < numPEP+1; i++){
    funcs[i]->SetParameters(p[5],  (i-1)*p[7], p[8]);
  }
  funcs[16]->SetParameters(p[3],    0.0, p[4]);
    
  //find bin edges 
  int iBin = h->FindBin(x[0]-p[2]);
  double ll = h->GetBinLowEdge(iBin);
  double ul = ll + h->GetBinWidth(iBin);

  //calculate contribution to bin of interest from each sub-function
  double g = funcs[0]->Integral(ll, ul)/(ul - ll);
  for(i = 1; i < nPeaks_+1; i++){
    //std::string gnum = "g" + std::to_string(i);
    g += sc[i] * funcs[i]->Integral(ll, ul)/(ul - ll);
  }
  g += funcs[i]->Integral(ll, ul)/(ul - ll);
    
  return g;
}

//Functional fit for pedestal plus landgauss MIP peak
double FitTestBeam::mipFunc(double* x, double* p)
{
  /*
    parameters
    p[10] : Total MIP peak area
    p[11] : Most probable value of Landau
    p[12] : Landau width parameter
    p[13] : gaussian width 
  */

  //find bin edges
  int iBin = h->FindBin(x[0]);;
  double ll = h->GetBinLowEdge(iBin);
  double ul = ll + h->GetBinWidth(iBin);

  //calculate pedestal, photoelectron, and background contribution
  double ppe = ppeFunc(x, p);

  //set mip parameters
  funcMIP->SetParameters(p[12], p[11], p[10], p[13]);
    
  //calculate MIP peak contribution
  //double mip = funcMIP->Integral(ll, ul)/(ul - ll);
  double mip = funcMIP->Eval(x[0]);
  //return bin value
  return ppe + mip;
}
