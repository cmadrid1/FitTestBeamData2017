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
    TLegend* leg = new TLegend(0.65, 0.8, 0.99, 0.9);
    leg->SetTextSize(0.03);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hfit,"Testbeam Data","l");
    //c1.SetLogx();
    //c1.SetLogy();
    //hfit->GetXaxis()->SetRangeUser(10, 400);
    //hfit->GetXaxis()->SetRangeUser(10, 4000);
    //hfit->SetMinimum(0.0000001);
    //hfit->SetMaximum(200);
    hfit->SetStats(false);
    hfit->SetTitle("Fit Testbeam Data");
    hfit->SetLineColor(kBlack);
    hfit->GetYaxis()->SetTitle("A.U.");
    hfit->GetXaxis()->SetTitle("Time");
    hfit->SetTitleOffset(1,"X");
    hfit->SetTitleOffset(1.2,"Y");
    hfit->SetTitleSize(0.05,"X");
    hfit->SetTitleSize(0.05,"Y");
    hfit->Draw("hist E");
    leg->Draw();
    
    //////////////////////
    //Fitting info
    //////////////////////
                                                       //////////////Parameter Info////////////////////////////////
    double min1   =     0; double max1   =   10;       //p[1] : Total MIP peak area------------//                //
    double min2   =     0; double max2   =  150;       //p[2] : Most probable value of Landau--//(MPV)           //
    double min3   =     0; double max3   =   10;       //p[3] : Landau width parameter---------//                //
    double min4   =     0; double max4   =   10;       //p[4] : gaussian width-----------------//                //
    double set1   =     6;                             ////////////////////////////////////////////////////////////
    double set2   =   113;
    double set3   =   0.5;
    double set4   = 0.001;
    double fitmin =   100; double fitmax = 200;
    
    TF1* fit1 = new TF1("Lang", langaufun, 0.0, 100000.0, 4);
    fit1->SetParLimits(10, min1, max1); 
    fit1->SetParLimits(11, min2, max2); 
    fit1->SetParLimits(12, min3, max3); 
    fit1->SetParLimits(13, min4, max4); 
    fit1->SetParameter(0, set1);
    fit1->SetParameter(1, set2);
    fit1->SetParameter(2, set3);
    fit1->SetParameter(3, set4);
    fit1->SetLineWidth(2);
    fit1->SetLineColor(kRed);
    hfit->Fit(fit1, "RQML", "", fitmin, fitmax);
    fit1->Draw("same");
    leg->AddEntry(fit1,"Fit","l");

    printf(
           "Chi^2:%10.4f, P1:%10.4f, P2:%10.4f, P3:%10.4f, P4:%10.4f\n",
           fit1->GetChisquare(), fit1->GetParameter(0), fit1->GetParameter(1), fit1->GetParameter(2), fit1->GetParameter(3)
           );

    c1.Print("FittedPlot.pdf");
}

void mipFitsSiPM()
{
    gROOT->SetStyle("Plain");    
    TFile *fin = TFile::Open("10-29.root");
    TH1* htofit = (TH1*)fin->Get("diff");
    fitSPEMIP(htofit, 1);
}

int main()
{
    mipFitsSiPM();
}
