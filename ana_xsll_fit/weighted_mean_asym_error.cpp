#include <iostream>
#include <iomanip>
#include <math.h>

#include "../Util/Style.h"
#include "../Util/Canvas.h"
#include "../Util/const.h"
#include "../Util/MChain.h"
#include "../Util/MCut.h"
#include "../Util/MCut_array.h"
#include "../Util/FitFunc.h"
#include "../Set/Nominal_cut_selection.h"
#include "../Set/Branch.h"
#include "../Set/makeCut.h"

#include "draws_.h"

#include <vector>
#include <stdlib.h>
#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TString.h>
#include <TF1.h>
#include <TFile.h>
#include <TArrow.h>
#include <TMath.h>
#include <Math/DistFunc.h>
#include <TRandom.h>

Int_t main( Int_t argc, Char_t** argv ){

  TApplication app( "app", &argc, argv );
  Style();
  Double_t x1[][5] = { // {value, stat(+), stat(-), syst(+), syst(-)}
    //{0.006,  0.007,  0.007,  0.0,    0.0}, // PDG check1
    //{0.0027, 0.0015, 0.0015, 0.0037, 0.0017},
    //{0.210, 0.016, 0.016, 0.0,   0.0}, // PDG check2
    //{0.264, 0.009, 0.009, 0.010, 0.026},
    { 0.00, 0.08, 0.07, 0.01, 0.01}, // LHCb (1st bin)
    {-0.20, 0.08, 0.07, 0.01, 0.03}, // LHCb (2nd bin)
    {-0.35, 0.26, 0.23, 0.10, 0.10}, // CDF  (1st bin)
    { 0.29, 0.32, 0.35, 0.15, 0.15}, // CDF  (2nd bin)
    { 0.47, 0.26, 0.32, 0.03, 0.03}, // Belle(1st bin)
    { 0.11, 0.31, 0.36, 0.07, 0.07}, // Belle(2nd bin)
    { 0.14, 0.15, 0.16, 0.20, 0.20}, // BaBar(1st bin)
    { 0.40, 0.18, 0.22, 0.07, 0.07}, // BaBar(2nd bin)
  };
  Int_t n = sizeof(x1)/sizeof(x1[0]);
  Double_t** x2 = new Double_t*[n]; // {value, err(+), err(-), err}
  for( Int_t i=0; i<n; i++ ){
    x2[i] = new Double_t[4];
    x2[i][0] = x1[i][0];
    x2[i][1] = sqrt(x1[i][1]*x1[i][1] + x1[i][3]*x1[i][3]);
    x2[i][2] = sqrt(x1[i][2]*x1[i][2] + x1[i][4]*x1[i][4]);
    x2[i][3] = sqrt(x2[i][1]*x2[i][1]+x2[i][2]*x2[i][2])/2.0;
    std::cout << "var" << i << " : "
	      << std::setw(8) << std::right << x1[i][0] << " + "
	      << std::setw(8) << std::right << x1[i][1] << "(stat) - "
	      << std::setw(8) << std::right << x1[i][2] << "(stat) + "
	      << std::setw(8) << std::right << x1[i][3] << "(syst) - "
	      << std::setw(8) << std::right << x1[i][4] << "(syst) "
	      << " -> "
	      << std::setw(8) << std::right << x2[i][0] << " + "
	      << std::setw(8) << std::right << x2[i][1] << " - "
	      << std::setw(8) << std::right << x2[i][2] << std::endl;
  }

  
  Double_t mean   = 1000000000;
  Double_t d_mean = 100;
  Double_t numerator   = 0; // bunnshi
  Double_t denominator = 0; // bunnbo
  while( d_mean > 0.00000001 ){
    for( Int_t i=0; i<n; i++ ){
      numerator   += x2[i][0]/x2[i][3]/x2[i][3];
      denominator +=      1.0/x2[i][3]/x2[i][3];
    }
    d_mean = mean;
    mean = numerator/denominator;
    d_mean -= mean;
    d_mean = fabs(d_mean);
    for( Int_t i=0; i<n; i++ ){
      if     ( mean > x2[i][0]+x2[i][1] ) x2[i][3] = x2[i][1];
      else if( mean < x2[i][0]-x2[i][2] ) x2[i][3] = x2[i][2];
      else                                x2[i][3] = x2[i][2] + (x2[i][1]-x2[i][2])*( mean - (x2[i][0]-x2[i][2]) )/( x2[i][1]+x2[i][2] );
    }
    //std::cout << "mean = " << mean << " (delta = " << d_mean << ")" << std::endl;
  }

  Double_t err_p = 0;
  Double_t err_m = 0;
  for( Int_t i=0; i<n; i++ ){
    err_p += 1/x2[i][1]/x2[i][1];
    err_m += 1/x2[i][2]/x2[i][2];
  }
  err_p = sqrt(1/err_p);
  err_m = sqrt(1/err_m);

  std::cout << "weighted mean = " << mean << " + " << err_p << " - " << err_m << std::endl;
  
  
  return 0;
}
