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
//const Double_t xmax = 14.0;
const Double_t xmax = 25.0;
const Double_t ymin = -1.0;
const Double_t ymax = 2.0;
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Int_t main( Int_t argc, Char_t** argv ){
  TApplication app( "app", &argc, argv );
  TStyle* sty = Style();
  sty->SetFuncStyle(1);
  sty->SetFuncWidth(2);
  //void plot_form_factor(){

  TF1*  V_Kstr_ball      = new TF1(  "V_Kstr_ball",      func_ball_59, xmin, xmax, 4 );  V_Kstr_ball     ->SetParNames  ( "r1", "mR2", "r2", "mfit2" );
  TF1* A0_Kstr_ball      = new TF1( "A0_Kstr_ball",      func_ball_59, xmin, xmax, 4 ); A0_Kstr_ball     ->SetParNames  ( "r1", "mR2", "r2", "mfit2" );
  TF1* A1_Kstr_ball      = new TF1( "A1_Kstr_ball",      func_ball_61, xmin, xmax, 2 ); A1_Kstr_ball     ->SetParNames  (              "r2", "mfit2" );
  TF1* A2_Kstr_ball      = new TF1( "A2_Kstr_ball",      func_ball_60, xmin, xmax, 3 ); A2_Kstr_ball     ->SetParNames  ( "r1",        "r2", "mfit2" );
  TF1* A3_Kstr_ball      = new TF1( "A3_Kstr_ball",      func_ball_A3, xmin, xmax, 5 ); A3_Kstr_ball     ->SetParNames  (              "r2", "mfit2", "r1_", "r2_", "mfit2_" );
  TF1* T1_Kstr_ball      = new TF1( "T1_Kstr_ball",      func_ball_59, xmin, xmax, 4 ); T1_Kstr_ball     ->SetParNames  ( "r1", "mR2", "r2", "mfit2" );
  TF1* T2_Kstr_ball      = new TF1( "T2_Kstr_ball",      func_ball_61, xmin, xmax, 2 ); T2_Kstr_ball     ->SetParNames  (              "r2", "mfit2" );
  TF1* T3_Kstr_ball      = new TF1( "T3_Kstr_ball",      func_ball_T3, xmin, xmax, 5 ); T3_Kstr_ball     ->SetParNames  ( "r1",        "r2", "mfit2", "r2_", "mfit2_" );
  TF1* T3tilde_Kstr_ball = new TF1( "T3tilde_Kstr_ball", func_ball_60, xmin, xmax, 3 ); T3tilde_Kstr_ball->SetParNames  ( "r1",        "r2", "mfit2" );

  V_Kstr_ball      ->SetParameters(  0.923, pow(5.32,2), -0.511, 49.40 );                        V_Kstr_ball     ->SetLineColor(2);
  A0_Kstr_ball     ->SetParameters(  1.364, pow(5.28,2), -0.990, 36.78 );                       A0_Kstr_ball     ->SetLineColor(2);
  A1_Kstr_ball     ->SetParameters(                       0.290, 40.38 );                       A1_Kstr_ball     ->SetLineColor(2);
  A2_Kstr_ball     ->SetParameters( -0.084,               0.342, 52.00 );                       A2_Kstr_ball     ->SetLineColor(2);
  A3_Kstr_ball     ->SetParameters(                       0.290, 40.38, -0.084, 0.342, 52.00 ); A3_Kstr_ball     ->SetLineColor(2);
  T1_Kstr_ball     ->SetParameters(  0.823, pow(5.32,2), -0.491, 46.31 );                       T1_Kstr_ball     ->SetLineColor(2);
  T2_Kstr_ball     ->SetParameters(                       0.333, 41.41 );                       T2_Kstr_ball     ->SetLineColor(2);
  T3_Kstr_ball     ->SetParameters( -0.036,               0.368, 48.10, 0.333, 41.41 );         T3_Kstr_ball     ->SetLineColor(2);
  T3tilde_Kstr_ball->SetParameters( -0.036,               0.368, 48.10 );                       T3tilde_Kstr_ball->SetLineColor(4);

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  TF1* A1_Kstr_ali = new TF1( "A1_Kstr_ali", func_ali, xmin, xmax, 4 ); A1_Kstr_ali->SetParNames  ( "F0", "c1", "c2", "c3" );
  TF1* A2_Kstr_ali = new TF1( "A2_Kstr_ali", func_ali, xmin, xmax, 4 ); A2_Kstr_ali->SetParNames  ( "F0", "c1", "c2", "c3" );
  TF1* A0_Kstr_ali = new TF1( "A0_Kstr_ali", func_ali, xmin, xmax, 4 ); A0_Kstr_ali->SetParNames  ( "F0", "c1", "c2", "c3" );
  TF1* A3_Kstr_ali = new TF1( "A3_Kstr_ali", func_ali, xmin, xmax, 8 ); A3_Kstr_ali->SetParNames  ( "F0", "c1", "c2", "c3", "F0_", "c1_", "c2_", "c3_" );
  TF1*  V_Kstr_ali = new TF1(  "V_Kstr_ali", func_ali, xmin, xmax, 4 );  V_Kstr_ali->SetParNames  ( "F0", "c1", "c2", "c3" );
  TF1* T1_Kstr_ali = new TF1( "T1_Kstr_ali", func_ali, xmin, xmax, 4 ); T1_Kstr_ali->SetParNames  ( "F0", "c1", "c2", "c3" );
  TF1* T2_Kstr_ali = new TF1( "T2_Kstr_ali", func_ali, xmin, xmax, 4 ); T2_Kstr_ali->SetParNames  ( "F0", "c1", "c2", "c3" );
  TF1* T3_Kstr_ali = new TF1( "T3_Kstr_ali", func_ali, xmin, xmax, 4 ); T3_Kstr_ali->SetParNames  ( "F0", "c1", "c2", "c3" );

  A1_Kstr_ali->SetParameters( 0.337, 0.602, 0.258, 0.000 );                             A1_Kstr_ali->SetLineColor(3);
  A2_Kstr_ali->SetParameters( 0.282, 1.172, 0.567, 0.000 );                             A2_Kstr_ali->SetLineColor(3);
  A0_Kstr_ali->SetParameters( 0.471, 1.505, 0.710, 0.000 );                             A0_Kstr_ali->SetLineColor(3);
  A3_Kstr_ali->SetParameters( 0.337, 0.602, 0.258, 0.000, 0.282, 1.172, 0.567, 0.000 ); A3_Kstr_ali->SetLineColor(3);
  V_Kstr_ali ->SetParameters( 0.457, 1.482, 1.015, 0.000 );                              V_Kstr_ali->SetLineColor(3);
  T1_Kstr_ali->SetParameters( 0.379, 1.519, 1.030, 0.000 );                             T1_Kstr_ali->SetLineColor(3);
  T2_Kstr_ali->SetParameters( 0.379, 0.517, 0.426, 0.000 );                             T2_Kstr_ali->SetLineColor(3);
  T3_Kstr_ali->SetParameters( 0.260, 1.129, 1.128, 0.000 );                             T3_Kstr_ali->SetLineColor(3);


  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TLegend* leg = new TLegend( 0.65,0.15,0.85,0.40 );
  leg->SetFillColor(10);
  leg->AddEntry( V_Kstr_ali,  "Ali's FFs",  "L" );
  leg->AddEntry( V_Kstr_ball, "Ball's FFs", "L" );

  TCanvas* can = new TCanvas( "can", "can", 1400, 1000 );
  can->Draw();
  can->Divide(3,3);
  can->cd(1);
  gPad->DrawFrame(xmin, ymin, xmax, ymax, "V;q^{2} [GeV^{2}]");
  V_Kstr_ali ->Draw("same");
  V_Kstr_ball->Draw("same");
  leg->Draw();

  can->cd(2);
  gPad->DrawFrame(xmin, ymin, xmax, ymax, "A_{0};q^{2} [GeV^{2}]");
  A0_Kstr_ali ->Draw("same");
  A0_Kstr_ball->Draw("same");

  can->cd(3);
  gPad->DrawFrame(xmin, ymin, xmax, ymax, "A_{1};q^{2} [GeV^{2}]");
  A1_Kstr_ali ->Draw("same");
  A1_Kstr_ball->Draw("same");

  can->cd(4);
  gPad->DrawFrame(xmin, ymin, xmax, ymax, "A_{2};q^{2} [GeV^{2}]");
  A2_Kstr_ali ->Draw("same");
  A2_Kstr_ball->Draw("same");

  can->cd(5);
  gPad->DrawFrame(xmin, ymin, xmax, ymax, "A_{3};q^{2} [GeV^{2}]");
  A3_Kstr_ali ->Draw("same");
  A3_Kstr_ball->Draw("same");

  can->cd(6);
  gPad->DrawFrame(xmin, ymin, xmax, ymax, "T_{1};q^{2} [GeV^{2}]");
  T1_Kstr_ali ->Draw("same");
  T1_Kstr_ball->Draw("same");

  can->cd(7);
  gPad->DrawFrame(xmin, ymin, xmax, ymax, "T_{2};q^{2} [GeV^{2}]");
  T2_Kstr_ali ->Draw("same");
  T2_Kstr_ball->Draw("same");

  can->cd(8);
  gPad->DrawFrame(xmin, ymin, xmax, ymax, "T_{3};q^{2} [GeV^{2}]");
  T3_Kstr_ali ->Draw("same");
  T3_Kstr_ball->Draw("same");
  T3tilde_Kstr_ball->Draw("same");

  can->Update();
  can->Print("pic/form_factor_comparison.eps");
  std::cout << "finish" << std::endl;
  if( !gROOT->IsBatch() ) app.Run();

  return 0;
}
