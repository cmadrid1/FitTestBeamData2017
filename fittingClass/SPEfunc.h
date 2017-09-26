#ifndef SPEfunc_h
#define SPEfunc_h

#include <memory>

class FitTestBeam;

//Fit function for pedestal and photoelectron spectrum
class PPEFunc
{
  //Paper on SiPM pixel crosstalk model
  //http://arxiv.org/pdf/1302.1455.pdf
  //n is the number of neighbor cells in the crosstalk model 
  const double n = 4;
  bool domip_;
  //CPn are the constant combinatoric scale factors based on n
  std::vector <double> CPn;

  //More ugly globals to hold utility functions used in fits
  std::vector <TF1*> funcs;
  TF1* funcMIP;
  double nPeaks_;
  
  //s defines the poisson factor for each PE peak based upon the poisson mean PE p[5]
  //in absense of cross-talk this would perfectly describe the PE peak amplitude ratios
  std::vector <double> s;
  std::vector <double> cp;
  std::vector <double> sc;
  
public:
  TH1* h;
  //Constructor
  PPEFunc(int nPeak, bool domip) : s(nPeak+1,0.0), cp(nPeak+1,0.0), sc(nPeak+1,0.0)
  {
    nPeaks_ = nPeak;
    domip_ = domip;
    
    //calculation of combinatoric prefactors for cross-talk 
    //CPn.push_back(0);
    //CPn.push_back(1);
    //CPn.push_back(n);
    //CPn.push_back(0.5*n*(3*n - 1));
    //CPn.push_back((1.0/3.0)*n*(8*pow(n, 2) - 6*n + 1));
    //CPn.push_back((0.25*n*((125.0/6.0)*pow(n, 3) - 25*pow(n, 2) + (55.0/6.0)*n - 1)));

    CPn.push_back(0);
    CPn.push_back(1);
    CPn.push_back(n);
    CPn.push_back((n*(-1 + 3*n))/2.);
    CPn.push_back((n*(1 - 6*n + 8*pow(n,2)))/3.);
    CPn.push_back((n*(-3 + 5*n)*(-2 + 5*n)*(-1 + 5*n))/24.);
    CPn.push_back((n*(-1 + 2*n)*(-2 + 3*n)*(-1 + 3*n)*(-1 + 6*n))/10.);
    CPn.push_back((n*(-5 + 7*n)*(-4 + 7*n)*(-3 + 7*n)*(-2 + 7*n)*(-1 + 7*n))/720.);
    CPn.push_back((n*(-1 + 2*n)*(-3 + 4*n)*(-1 + 4*n)*(-5 + 8*n)*(-3 + 8*n)*(-1 + 8*n))/315.);
    CPn.push_back((n*(-2 + 3*n)*(-1 + 3*n)*(-7 + 9*n)*(-5 + 9*n)*(-4 + 9*n)*(-2 + 9*n)*(-1 + 9*n))/4480.);
    CPn.push_back((n*(-1 + 2*n)*(-4 + 5*n)*(-3 + 5*n)*(-2 + 5*n)*(-1 + 5*n)*(-7 + 10*n)*(-3 + 10*n)*(-1 + 10*n))/4536.);
    CPn.push_back((n*(-9 + 11*n)*(-8 + 11*n)*(-7 + 11*n)*(-6 + 11*n)*(-5 + 11*n)*(-4 + 11*n)*(-3 + 11*n)*(-2 + 11*n)*(-1 + 11*n))/3.6288e6);
    CPn.push_back((n*(-1 + 2*n)*(-2 + 3*n)*(-1 + 3*n)*(-3 + 4*n)*(-1 + 4*n)*(-5 + 6*n)*(-1 + 6*n)*(-7 + 12*n)*(-5 + 12*n)*(-1 + 12*n))/11550.);
    CPn.push_back((n*(-11 + 13*n)*(-10 + 13*n)*(-9 + 13*n)*(-8 + 13*n)*(-7 + 13*n)*(-6 + 13*n)*(-5 + 13*n)*(-4 + 13*n)*(-3 + 13*n)*(-2 + 13*n)*(-1 + 13*n))/4.790016e8);
    CPn.push_back((n*(-1 + 2*n)*(-6 + 7*n)*(-5 + 7*n)*(-4 + 7*n)*(-3 + 7*n)*(-2 + 7*n)*(-1 + 7*n)*(-11 + 14*n)*(-9 + 14*n)*(-5 + 14*n)*(-3 + 14*n)*(-1 + 14*n))/1.38996e7);
    CPn.push_back((pow(n,2)*(19802759040 + n*(-425547472896 + n*(4895614496136 + n*(-36007056781696 + n*(180954419071426 + n*(-643524004371248 + n*(1648017660281993 + n*(-3051606955645248 + n*(4048835387589993 + n*(-3751975766595056 + n*(2305159936523171 + n*(-843351196289056 + 139013933454241*n)))))))))))))/6.2270208e9);
    
    //configure utility functions for fits
    funcs.push_back(new TF1("ped", "gaus"));

    for(int i = 1; i <= nPeaks_; i++){
      //fill funcs
      std::string pnum = "pe" + std::to_string(i);
      funcs.push_back(new TF1(pnum.c_str(), "gaus"));
    }

    FitTestBeam* fitFunc = new FitTestBeam();
    
    funcs.push_back(new TF1("bg",  "landau"));
    funcMIP = new TF1("mip",  fitFunc->langaufun, 0, 400000, 4);

  }

  FitTestBeam* fitFunc = new FitTestBeam();
  
  double operator()(double* x, double* p){
    double returnVal = fitFunc->ppeFunc(x,p);
    if(domip_){
      returnVal += fitFunc->mipFunc(x,p);
    }
    return returnVal;
  }
};

#endif
