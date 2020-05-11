#ifndef MAXWELL_PARAMETERS_H
#define MAXWELL_PARAMETERS_H

#include <parameters.h>

class MaxwellParameters : public Parameters {
  public:
    MaxwellParameters() : Parameters(1){
      mnpoints = 101;
      mxmin = 0.0;
      mxmax = 1.0;
      mtmax = 1.0;
      mcfl = 1.0;
      moutput_frequency = 1;
      mid_gauss_amp = 1.0;
      mid_gauss_center = 0.0;
      mid_gauss_width = 1.0;
      mKOSigma = 0.0;
    }

    inline void setnpoints(int npoints){
      mnpoints = npoints;
    }

    inline int getnpoints(){
      return mnpoints;
    }

    inline void setxmin(double xmin){
      mxmin = xmin;
    }

    inline double getxmin(){
      return mxmin;
    }

    inline void setxmax(double xmax){
      mxmax = xmax;
    }

    inline double getxmax(){
      return mxmax;
    }

    inline void settmax(double tmax){
      mtmax = tmax;
    }

    inline double gettmax(){
      return mtmax;
    }

    inline void setcfl(double cfl){
      mcfl = cfl;
    }

    inline double getcfl(){
      return mcfl;
    }

    inline void setoutput_frequency(int output_frequency){
      moutput_frequency = output_frequency;
    }

    inline int getoutput_frequency(){
      return moutput_frequency;
    }

    inline void setid_gauss_amp(double id_gauss_amp){
      mid_gauss_amp = id_gauss_amp;
    }

    inline double getid_gauss_amp(){
      return mid_gauss_amp;
    }

    inline void setid_gauss_center(double id_gauss_center){
      mid_gauss_center = id_gauss_center;
    }

    inline double getid_gauss_center(){
      return mid_gauss_center;
    }

    inline void setid_gauss_width(double id_gauss_width){
      mid_gauss_width = id_gauss_width;
    }

    inline double getid_gauss_width(){
      return mid_gauss_width;
    }

    inline void setKOSigma(double KOSigma){
      mKOSigma = KOSigma;
    }

    inline double getKOSigma(){
      return mKOSigma;
    }

  private:
    int mnpoints;
    double mxmin;
    double mxmax;
    double mtmax;
    double mcfl;
    int moutput_frequency;
    double mid_gauss_amp;
    double mid_gauss_center;
    double mid_gauss_width;
    double mKOSigma;
};

#endif
