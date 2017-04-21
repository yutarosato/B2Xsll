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

#include <vector>
#include <stdlib.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TString.h>
#include <TF1.h>
#include <TArrow.h>
#include <TFile.h>

const Double_t sig_mean [3] = { 5.27936, 5.27932, 5.27934 };
const Double_t sig_sigma[3] = { 0.00267, 0.00258, 0.00263 };


Int_t main( Int_t argc, Char_t** argv ){
  TApplication app( "app", &argc, argv );
  Style();

  //TFile file( "pdf/Mbc_peak_cc2_s00.root"          );
  //TFile file( "pdf/Mbc_peak_cc2_s00-5.root"        );
  //TFile file( "pdf/Mbc_peak_swap2_s00.root"        );
  //TFile file( "pdf/Mbc_peak_double2_s00.root"      );
  //TFile file( "pdf_data/Mbc_peak_swap2_s00.root"   );
  TFile file( "pdf_data/Mbc_peak_double2_s00.root" );
  if( file.IsZombie() ) std::cerr << "[ABORT] can not find input-file : " << file.GetName() << std::endl, abort();

  std::cout << "[FILE] " << file.GetName() << std::endl;

  TH1D** hist = new TH1D*[3];
  hist[0] = (TH1D*)file.Get( "hist8"  ); hist[0]->SetLineColor(1); hist[0]->SetMarkerColor(1); // ee
  hist[1] = (TH1D*)file.Get( "hist17" ); hist[1]->SetLineColor(2); hist[1]->SetMarkerColor(2); // mm
  hist[2] = (TH1D*)file.Get( "hist26" ); hist[2]->SetLineColor(3); hist[2]->SetMarkerColor(3); // ll

  const Int_t sel_fun = 15;
  TF1** func = new TF1*[3];
  for( Int_t i=0; i<3; i++ ){
    func[i] = new TF1( Form("func%d",i), make_func(sel_fun), 5.22, 5.30, n_fitfunc_par(sel_fun) );
    func[i]->SetParNames  ( "area",      "mean",       "sigma", "norm","Ebeam","shape" );
    func[i]->SetParameters(    0.5, sig_mean[i],  sig_sigma[i],      3,  5.289,   -10  );
    func[i]->FixParameter( 1, sig_mean [i] );
    func[i]->FixParameter( 2, sig_sigma[i] );
    func[i]->FixParameter( 4, 5.289        );
  }
  
  TCanvas* c1 = Canvas( "c1","c1", 3, 1 );
  c1->Draw();

  for( Int_t i=0; i<3; i++ ){
    c1->cd(i+1);
    hist[i]->Fit( func[i] );

    //hist[i]->Draw();    
  }

  
  for( Int_t i=0; i<3; i++ ) std::cout << Form( "YIELDS%d  : ", i )
				       << func[i]->GetParameter(0)/hist[0]->GetBinWidth (1)
				       << std::endl;
  c1->Update();
  c1->Print( "pic/Mbc_peak_file_fit.eps" );
  if( !gROOT->IsBatch() ) app.Run();
  std::cout << "finish" << std::endl;
  return 0;
}
