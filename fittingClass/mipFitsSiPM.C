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
    //c1.SetLogx();
    //c1.SetLogy();
    //hfit->GetXaxis()->SetRangeUser(10, 400);
    hfit->GetXaxis()->SetRangeUser(10, 4000);
    //hfit->SetMinimum(0.0000001);
    //hfit->SetMaximum(200);
    
    hfit->Draw("hist");
                                                       //////////////Parameter Info////////////////////////////////
    double min1   =     0; double max1   =   10;       //p[1] : Total MIP peak area------------//                //
    double min2   =     0; double max2   =  150;       //p[2] : Most probable value of Landau--//(MPV)           //
    double min3   =     0; double max3   =   10;       //p[3] : Landau width parameter---------//                //
    double min4   =     0; double max4   =   10;       //p[4] : gaussian width-----------------//                //
    double set1   =     6;                             ////////////////////////////////////////////////////////////
    double set2   =   113;
    double set3   =   0.5;
    double set4   = 0.001;
    double fitmin =   105; double fitmax = 250;
    
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
    fit1->SetLineColor(kBlack);
    hfit->Fit(fit1, "RQML", "", fitmin, fitmax);
    fit1->Draw("same");

    printf(
           "Chi^2:%10.4f, P1:%10.4f, P2:%10.4f, P3:%10.4f, P4:%10.4f\n",
           fit1->GetChisquare(), fit1->GetParameter(0), fit1->GetParameter(1), fit1->GetParameter(2), fit1->GetParameter(3)
           );

    hfit->GetYaxis()->SetTitle("Time [ns]");
    hfit->GetXaxis()->SetTitle("A.U.");
    hfit->SetTitleOffset(1,"X");
    hfit->SetTitleOffset(1.2,"Y");
    hfit->SetTitleSize(0.05,"X");
    hfit->SetTitleSize(0.05,"Y");
        
    //Make Plots Pretty
    //hfit->SetStats(false);
    //hfit->SetTitle("");
    //
    //TLatex* CMSPrelim1 = new TLatex(0.14, 0.91, "CMS #scale[0.9]{#font[52]{Preliminary}}");
    //CMSPrelim1->SetNDC();
    //CMSPrelim1->SetTextFont(62);
    //
    //TLatex* testbeam = new TLatex(0.95, 0.91, "HB Testbeam 2017");
    //testbeam->SetNDC();
    //testbeam->SetTextFont(42);
    //testbeam->SetTextAlign(31);
    //
    //TLatex* SiPMTitle = new TLatex(0.93, 0.86, "SiPM (Silicon Photomultiplier)");
    //SiPMTitle->SetNDC();
    //SiPMTitle->SetTextFont(42);
    //SiPMTitle->SetTextAlign(32);
    //
    //TLatex* Muon = new TLatex(0.93, 0.81, "Parasitic Muons");
    //Muon->SetNDC();
    //Muon->SetTextFont(42);
    //Muon->SetTextAlign(32);
    //
    //char chan [100];
    //int iEta, iPhi, iDepth;
    //sscanf (hfit->GetName(),"adc_nosub_binChris_%d_%d_%d", &iEta, &iPhi, &iDepth);
    //sprintf (chan, "Channel: %d,%d,%d", iEta, iPhi, iDepth);
    //
    //TLatex* channel = new TLatex(0.93, 0.86, chan);
    //channel->SetNDC();
    //channel->SetTextFont(42);
    //channel->SetTextAlign(32);
    //
    //char mpv [100];
    //float intmpv = fit1->GetParameter(11);
    //sprintf (mpv,"MPV: %0.3f", intmpv);
    //
    //TLatex* MPV = new TLatex(0.93, 0.8, mpv);
    //MPV->SetNDC();
    //MPV->SetTextFont(42);
    //MPV->SetTextAlign(31);
    //
    //char gain [100];
    //float intgain = fit1->GetParameter(7);
    //sprintf (gain,"Gain: %0.3f", intgain);
    //
    //TLatex* Gain = new TLatex(0.93, 0.75, gain);
    //Gain->SetNDC();
    //Gain->SetTextFont(42);
    //Gain->SetTextAlign(31);
    //
    //CMSPrelim1->Draw();
    //testbeam->Draw();
    ////SiPMTitle->Draw();
    ////Muon->Draw();
    //channel->Draw();
    //MPV->Draw();
    //Gain->Draw();
    
    char oname[128];
    sprintf(oname, "adcHist_%s_%i.pdf", hfit->GetName(), run);
    c1.Print(oname);
}

void mipFitsSiPM()
{
    gROOT->SetStyle("Plain");    
        
    char fname[] = "10-29.root";
    
    TFile *fin = TFile::Open(fname);
    
    TObject *obj;
    TH1* htofit = (TH1*)fin->Get("diff");
    fitSPEMIP(htofit, 1);
}

int main()
{
    mipFitsSiPM();
}
