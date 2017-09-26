#ifndef FitTestBeam_h
#define FitTestBeam_h

#include <memory>

class FitTestBeam
{
  public:
    double crystalBall(double* x, double* p);
    Double_t langaufun(Double_t* x, Double_t* par);
    double ppeFunc(double* x, double* p);
    double mipFunc(double* x, double* p);

};

#endif
