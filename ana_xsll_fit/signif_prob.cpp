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

// convert from significance to probability.
// NOTE : range should be less than about six.
//        too large range do not work correctly.
Int_t main( Int_t argc, Char_t** argv ){

  TApplication app( "app", &argc, argv );
  Style();
  if( argc!=2 ) std::cerr << "wrong input" << std::endl, abort();
  const Int_t    sel_fun =   10;
  Double_t range = atof(argv[1]);
  TF1* f = new TF1( "f", make_func(sel_fun), -range, range, n_fitfunc_par(sel_fun) );
  f->SetParNames  ( "area","mean", "sigma" );
  f->SetParameters(    1.0,   0.0,   1.0   );
  std::cout << range << "sigma : " << std::setprecision(15) << f->Integral( -range, range ) << std::endl;
  
  return 0;
}
