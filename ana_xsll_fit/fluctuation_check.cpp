#include <iostream>
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

#include <TApplication.h>
#include <TGraph.h>
#include <TF1.h>

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++;

Int_t main( Int_t argc, Char_t** argv ){
  TApplication app( "app", &argc, argv );
  Style();
  //++++++++++++++++++++++++++++++
  const Int_t   Nq2bin     = 4;
  const Bool_t  fl_message = true;
  const Bool_t  fl_save    = true;

  TGraph** g = new TGraph*[Nq2bin];
  for( Int_t i=0; i<Nq2bin; i++ ){
    g[i] = new TGraph( Form("tmp%d.dat",i), "%*lg %lg %lg" );
    std::cout << "[i = " << i << "] RMS = "<< g[i]->GetRMS() << std::endl;
  }
  return 0;
}


