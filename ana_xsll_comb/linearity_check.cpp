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
#include <TGraphErrors.h>
#include <TF1.h>

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++;

Int_t main( Int_t argc, Char_t** argv ){
  TApplication app( "app", &argc, argv );
  Style();
  //++++++++++++++++++++++++++++++
  const Char_t* file       = "tmp.dat";
  const Int_t   tmp_i      = 1; // 0(ee), 1(mm), 2(ee+mm)
  const Bool_t  fl_message = true;
  // RESIDUAL
  const Double_t Nsig[3] = { 156.9, 169.3, 326.2 }; // # of signal
  const Double_t Nbkg[3] = {  33.7,  34.3,  68.0 }; // # of peaking B.G.
  //++++++++++++++++++++++++++++++
  
  TCanvas* c1 = Canvas( "c1","c1",2, 2 );

  // READ DAT-FILE (g1)
  TGraphErrors* g1 = new TGraphErrors( file, "%lg %lg %lg" );
  g1->SetTitle( "Linearity Check" );
  g1->Sort();
  Double_t* X  = new Double_t[g1->GetN()];
  Double_t* Y  = new Double_t[g1->GetN()];
  Double_t* EY = new Double_t[g1->GetN()];
  X  = g1->GetX ();
  Y  = g1->GetY ();
  EY = g1->GetEY();

  // FIT FUNC (g1)
  TF1* f1 = new TF1( "f1", "[0]*x+[1]", 0.0, 1.5 );
  f1->SetTitle( "Linearity Check" );
  f1->SetParNames  ( "slope", "offset" );
  f1->SetParameters(    300 ,       10 );

  // FIT and DRAW (g1)
  c1->Draw();
  c1->cd(1);
  g1->Fit("f1");
  g1->Draw("AP");

  // SCALE (g2)
  TGraphErrors* g2 = new TGraphErrors();
  g2->SetTitle("Scale");
  for( Int_t i=0; i<g1->GetN(); i++ ){
    g2->SetPoint( i, X[i], (Y[i]-Nbkg[tmp_i])/Nsig[tmp_i] );
    g2->SetPointError( i, 0, EY[i]/Nsig[tmp_i] );
  }

  // FIT FUNC (g2)
  TF1* f2 = new TF1( "f2", "[0]*x+[1]", 0.0, 1.5 );
  f2->SetTitle( "Linearity Check" );
  f2->SetParNames  ( "slope", "offset" );
  f2->SetParameters(    1 ,       0 );

  // FIT and DRAW (g2)
  c1->cd(2);
  g2->Fit("f2");
  g2->Draw("AP");
  
  // RESIDUAL (g3)
  TGraphErrors* g3 = new TGraphErrors();
  g3->SetTitle("Obs. - True");
  for( Int_t i=0; i<g1->GetN(); i++ ){
    g3->SetPoint( i, X[i], Y[i]-(Nbkg[tmp_i]+X[i]*Nsig[tmp_i]) );
    g3->SetPointError( i, 0, EY[i] );
    if( fl_message ) std::cout << "i = " << i << " : " << X[i]  << ", " << Y[i] << " +- " << EY[i]
			       << " ([TRUE] " << Nbkg[tmp_i]+X[i]*Nsig[tmp_i] << " )"
			       << std::endl;
  }

  c1->cd(3);
  g3->Draw("AP");
  
  // SAVE
  c1->Print( "test.eps" );
  std::cout << "finish" << std::endl;
  app.Run();
  
  delete c1;
  delete g1;
  delete g2;
  delete g3;
  delete f1;
  delete f2;
  
  return 0;
}


