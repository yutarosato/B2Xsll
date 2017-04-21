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

#include "form_factor.h"

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
#include <TArrow.h>

const Double_t xmin =  0.0;
//const Double_t xmin =  1e-9;
const Double_t xmax = 25.0;
const Double_t ymax = 0.80;
// mail from Zwicky @20141217 [q2, T3]
const Int_t npoint = 29;
//const Double_t Lis[npoint][2] = {{0,   0.2018841720993929 }, {0.5,  0.20569382259668395},
const Double_t Lis[npoint][2] = {{0.0000001,   0.2018841720993929 }, {0.5,  0.20569382259668395},
				 {1.,  0.20959247349532054}, {1.5,  0.2135809987810585},
				 {2.,  0.21766007330084566}, {2.5,  0.22183003286497263},
				 {3.,  0.22609079700269283}, {3.5,  0.23044177651143646},
				 {4.,  0.23488176283536527}, {4.5,  0.2394087955517761},
				 {5.,  0.24402000368189006}, {5.5,  0.24871141524300153},
				 {6.,  0.25347772885863845}, {6.5,  0.2583120393783522},
				 {7.,  0.26320550727925957}, {7.5,  0.2681469612490331},
				 {8.,  0.2731224166388041 }, {8.5,  0.27811449618597583},
				 {9.,  0.2831017232175272 }, {9.5,  0.2880576631196321},
				 {10., 0.29294987733849265}, {10.5, 0.2977386368479898},
				 {11., 0.3023753399209453 }, {11.5, 0.3068005178824105},
				 {12., 0.3109414316226418 }, {12.5, 0.31470891449528876},
				 {13., 0.3179934561783064 }, {13.5, 0.3206598876443783},
				 {14., 0.32254139232832557}};

Int_t main( Int_t argc, Char_t** argv ){
  TApplication app( "app", &argc, argv );
  TStyle* sty = Style();
  sty->SetTitleOffset(1.2,"y");
  sty->SetFuncStyle(1);
 //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TCanvas* can = Canvas( "can", "can", 2, 1 );
  can->Draw();
  can->cd(1);
  gPad->DrawFrame( xmin, 0.0, xmax, ymax, "T_{3};q^{2} [GeV^{2}]; T_{3}" );
  //gPad->DrawFrame( 1e-9, 0.0, xmax, ymax, "T_{3};q^{2} [GeV^{2}]; T_{3}" );

  // data points from Zwicky
  TGraph* g = new TGraph();
  for( Int_t ip=0; ip<npoint; ip++ ) g->SetPoint( g->GetN(), Lis[ip][0], Lis[ip][1] );
  g->Draw("Psame");


  // function described in paper
  // Ball
  TF1* T3_Kstr_ball = new TF1( "T3_Kstr_ball", func_ball_T3, xmin, xmax, 5 );
  T3_Kstr_ball->SetParNames( "r1", "r2", "mfit2", "r2_", "mfit2_" );
  T3_Kstr_ball->SetParameters( -0.036, 0.368, 48.10, 0.333, 41.41 );
  T3_Kstr_ball->SetLineColor(2);
  T3_Kstr_ball->SetLineStyle(2);
  T3_Kstr_ball->DrawCopy("same");
  // Ali
  TF1* T3_Kstr_ali = new TF1( "T3_Kstr_ali", func_ali, xmin, xmax, 4 );
  T3_Kstr_ali->SetParNames( "F0", "c1", "c2", "c3" );
  T3_Kstr_ali->SetParameters( 0.260, 1.129, 1.128, 0.000 );
  T3_Kstr_ali->SetLineColor(3);
  T3_Kstr_ali->SetLineStyle(2);
  T3_Kstr_ali->DrawCopy("same");

  // fitting with function described in paper
  TF1* T3_Kstr_ball_fit = new TF1( "T3_Kstr_ball_fit", func_ball_T3, xmin, xmax, 5 );
  //TF1* T3_Kstr_ball_fit = new TF1( "T3_Kstr_ball_fit", func_ball_T3, 0.01, xmax, 5 );
  //TF1* T3_Kstr_ball_fit = new TF1( "T3_Kstr_ball_fit", func_ball_T3, 3.0, xmax, 5 );
  T3_Kstr_ball_fit->SetParNames  (   "r_{1,T3tilde}",  "r_{2,T3tilde}", "m_{fit,T3tilde}^{2}", "r_{2,T2}", "m_{fit,T2}^{2}" );
  T3_Kstr_ball_fit->SetParameters( -0.036, 0.368,   48.10, 0.333,   41.41 );
  //T3_Kstr_ball_fit->FixParameter(0, -0.036);
  //T3_Kstr_ball_fit->FixParameter(1, 0.368);
  //T3_Kstr_ball_fit->FixParameter(2, 48.10);
  T3_Kstr_ball_fit->FixParameter(3, 0.333);
  T3_Kstr_ball_fit->FixParameter(4, 41.41);
  T3_Kstr_ball_fit->SetLineColor(4);
  g->Fit("T3_Kstr_ball_fit","NR");
  T3_Kstr_ball_fit->DrawCopy("same");

  // fitting with polynomial
  TF1* func_poly = new TF1( "func_poly", "pol3", xmin, xmax);
  func_poly->SetLineColor(5);
  g->Fit("func_poly","N");
  func_poly->DrawCopy("same");

  TLegend* leg = new TLegend( 0.15,0.45,0.40,0.70 );
  leg->SetFillColor(10);
  leg->AddEntry( T3_Kstr_ball,     "Ball's FFs (paper)", "L" );
  leg->AddEntry( T3_Kstr_ali,      "Ali's FFs (paper)",  "L" );
  leg->AddEntry( T3_Kstr_ball_fit, "Ball's FFs (fit)",   "L" );
  leg->AddEntry( func_poly,        "3rd poly. fit",      "L" );
  leg->Draw();

  ///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  can->cd(2);
  gPad->DrawFrame( xmin, 0.0, xmax, ymax, "T_{3};q^{2} [GeV^{2}]; T_{3}" );
  g->Draw("Psame");

  Double_t threshold;
  for( Int_t i=0; i<2500; i++ ){
    double x = i*(xmax-xmin)/2500.0;
    double yfunc = func_poly   ->Eval( x );
    double yball = T3_Kstr_ball->Eval( x );
    //std::cout << i << " : " << x << " GeV2 : "
    //<< yfunc << " - " << yball << std::endl;
    if( yfunc < yball ){
      threshold = x;
      std::cout << "threshold = " << threshold << " GeV2" << std::endl;
      std::cout << "break" << std::endl;
      break;
    }
  }
  
  T3_Kstr_ball->SetRange( threshold, xmax      );
  func_poly   ->SetRange( xmin,      threshold );
  T3_Kstr_ball->SetLineWidth(2);
  func_poly   ->SetLineWidth(2);
  T3_Kstr_ball->Draw("same");
  func_poly   ->Draw("same");

  TText* tex = new TText();
  tex->SetTextColor(3);
  tex->SetTextSize(0.03);
  tex->DrawTextNDC( 0.2,0.66, "poly. fit + Ball's FFs (paper)"              );
  tex->DrawTextNDC( 0.2,0.62, Form("(threshold = %.2f [GeV2])", threshold ) );

  can->Update();
  can->Print("pic/fit_form_factor.eps");
  std::cout << "finish" << std::endl;
  if( !gROOT->IsBatch() ) app.Run();
  return 0;
}

