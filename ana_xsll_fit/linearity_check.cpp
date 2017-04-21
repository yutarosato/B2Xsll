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
  if( !(argc==2 || argc==3) ) std::cerr << "wrong input" << std::endl, abort();
  const Int_t   tmp_i     = atoi(argv[1]);
  Int_t         fl_appRun = 1;
  if( argc==3 ) fl_appRun = atoi( argv[2] );
  //++++++++++++++++++++++++++++++
  Char_t* file  = new Char_t[256];
  strcpy( file, Form("tmp%d.dat", tmp_i) );
  const Bool_t  fl_message = true;
  const Bool_t  fl_save    = true;
  // RESIDUAL
  //const Double_t Nsig[3] = { 140.7, 163.7,  304.4 }; // # of signal
  //const Double_t Nbkg[3] = {  29.3,  36.7,   66.0 }; // # of peaking B.G.
  const Double_t Nsig[3] = { 0, 0, 0 }; // # of signal       // tmpppppp
  const Double_t Nbkg[3] = { 0, 0, 0 }; // # of peaking B.G. // tmppppppp
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
    if( Nsig[tmp_i]!=0 ){
      g2->SetPoint( i, X[i], (Y[i]-Nbkg[tmp_i])/Nsig[tmp_i] );
      g2->SetPointError( i, 0, EY[i]/Nsig[tmp_i] );
    }else{
      g2->SetPoint( i, X[i], (Y[i]-Nbkg[tmp_i])/1.0         );
      g2->SetPointError( i, 0, EY[i]/1.0 );
    }

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

  // FIT FUNC (g3)
  TF1* f3 = new TF1( "f3", "[0]", 0.0, 1.5 );
  f3->SetTitle( "Linearity Check" );
  f3->SetParNames ( "offset" );
  f3->SetParameter( 0, 0 );

  // FIT and DRAW(g3)
  c1->cd(3);
  g3->Fit("f3");
  g3->Draw("AP");

  // SAVE
  if( fl_save ) c1->Print( Form("pic/linearity_check_%d.eps",tmp_i) );
  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();
  
  delete c1;
  delete g1;
  delete g2;
  delete g3;
  delete f1;
  delete f2;
  delete f3;
  
  return 0;
}


