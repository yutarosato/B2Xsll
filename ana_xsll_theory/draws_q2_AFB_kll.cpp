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
#include <TGraph.h>
#include <TFile.h>
#include <TLegend.h>
#include <TLine.h>

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++;

Int_t main( Int_t argc, Char_t** argv ){
  TApplication app( "app", &argc, argv );
  Style();
  //++++++++++++++++++++++++++++++
  Int_t fl_appRun = 1;
  const Bool_t  fl_message = true;
  const Bool_t  fl_save    = true;
  //++++++++++++++++++++++++++++++
  const Int_t ng   = 3; // [ee,mm,ee+mm(sim)]
  const Int_t nafb = 4;
  const Double_t div[ng][4] = {
    { (PDGmass::jpsi-0.325)*(PDGmass::jpsi-0.325), (PDGmass::jpsi+0.125)*(PDGmass::jpsi+0.125), (PDGmass::psi2s-0.20)*(PDGmass::psi2s-0.20), (PDGmass::psi2s+0.10)*(PDGmass::psi2s+0.10) }, // ll
    { (PDGmass::jpsi-0.400)*(PDGmass::jpsi-0.400), (PDGmass::jpsi+0.150)*(PDGmass::jpsi+0.150), (PDGmass::psi2s-0.25)*(PDGmass::psi2s-0.25), (PDGmass::psi2s+0.10)*(PDGmass::psi2s+0.10) }, // ee
    { (PDGmass::jpsi-0.250)*(PDGmass::jpsi-0.250), (PDGmass::jpsi+0.100)*(PDGmass::jpsi+0.100), (PDGmass::psi2s-0.15)*(PDGmass::psi2s-0.15), (PDGmass::psi2s+0.10)*(PDGmass::psi2s+0.10) }, // mm
  };
  const Double_t x[ng][nafb] = {
    {4.3/2.0, (4.3 + div[0][0])/2.0,  (div[0][1] + div[0][2])/2.0, (div[0][3] + 25.0)/2.0}, // ll
    {4.3/2.0, (4.3 + div[1][0])/2.0,  (div[1][1] + div[1][2])/2.0, (div[1][3] + 25.0)/2.0}, // ee
    {4.3/2.0, (4.3 + div[2][0])/2.0,  (div[2][1] + div[2][2])/2.0, (div[2][3] + 25.0)/2.0}, // mm
  };
  const Double_t xE[ng][nafb] = {
    {4.3/2.0, (-4.3 + div[0][0])/2.0,  (-div[0][1] + div[0][2])/2.0, (-div[0][3] + 25.0)/2.0}, // ll
    {4.3/2.0, (-4.3 + div[1][0])/2.0,  (-div[1][1] + div[1][2])/2.0, (-div[1][3] + 25.0)/2.0}, // ee
    {4.3/2.0, (-4.3 + div[2][0])/2.0,  (-div[2][1] + div[2][2])/2.0, (-div[2][3] + 25.0)/2.0}, // mm
  };

  const Double_t afb_true[ng][nafb] = {
    {-0.047,          -0.113,            NULL,               0.12            }, // ee+mm(sim)
    {1.2893*(-0.077), 1.13893*0.01,      1.06319*(-0.350),  1.12141*(-0.107) }, // ee
    {2.08178*0.06,    1.37547*(-0.182),  NULL,              1.08232*0.29     }, // mm
  };
  const Double_t afb_trueE[ng][nafb] = {
    {0.24,         0.29,          NULL,         0.18         },  // ee+mm(sim)
    {1.2893*0.21,  1.13893*0.33,  1.06319*0.38, 1.12141*0.23 }, // ee
    {2.08178*0.23, 1.37547*0.31,  NULL,         1.08232*0.20 }, // mm
  };

  const Double_t xmin =  0.0;
  const Double_t xmax = 25.0;
  const Double_t ymin = -1.0;
  const Double_t ymax =  1.0;
  //++++++++++++++++++++++++++++++ 

  // MAKE TGRAPHERRORS (gt)
  TGraphErrors** gt = new TGraphErrors*[ng];
  for( Int_t i=0; i<ng; i++ ){
    gt[i] = new TGraphErrors();
    for( Int_t j=0; j<nafb; j++ ){
      gt[i]->SetPoint     ( j, x [i][j], afb_true [i][j] );
      gt[i]->SetPointError( j, xE[i][j], afb_trueE[i][j] );
    }
    gt[i]->Sort();
    Deco( gt[i], 0, i+1, i+1 );
    gt[i]->SetMarkerSize(0.8);
  }

  // MAKE LEGEND(leg1)
  TLegend* leg1 = new TLegend( 0.75,0.75,0.99,0.99 );
  leg1->AddEntry( gt[0],"ll",     "PL" );
  leg1->AddEntry( gt[1],"ee",     "PL" );
  leg1->AddEntry( gt[2],"#mu#mu", "PL" );

  // MAKE LINE(l1)
  TLine*** l1 = new TLine**[ng];
  for( Int_t i=0; i<ng; i++ ){
    l1[i] = new TLine*[4];
    for( Int_t j=0; j<4;  j++ ){
      l1[i][j] = new TLine( div[i][j], ymin, div[i][j], ymax );
      l1[i][j]->SetLineStyle(2);
      l1[i][j]->SetLineColor(i+1);
    }
  }
  // MAKE BAND(b1)
  TGraph** b1 = new TGraph*[4];
  b1[0] = new TGraph();
  b1[1] = new TGraph();
  b1[2] = new TGraph();
  b1[3] = new TGraph();
  b1[0]->SetPoint( 0, div[1][0], -1.0 ); b1[0]->SetPoint( 1, div[1][0],  1.0 ); b1[0]->SetPoint( 2, div[1][1], 1.0 ); b1[0]->SetPoint( 3, div[1][1],  -1.0 );
  b1[1]->SetPoint( 0, div[1][2], -1.0 ); b1[1]->SetPoint( 1, div[1][2],  1.0 ); b1[1]->SetPoint( 2, div[1][3], 1.0 ); b1[1]->SetPoint( 3, div[1][3],  -1.0 );
  b1[2]->SetPoint( 0, div[2][0], -1.0 ); b1[2]->SetPoint( 1, div[2][0],  1.0 ); b1[2]->SetPoint( 2, div[2][1], 1.0 ); b1[2]->SetPoint( 3, div[2][1],  -1.0 );
  b1[3]->SetPoint( 0, div[2][2], -1.0 ); b1[3]->SetPoint( 1, div[2][2],  1.0 ); b1[3]->SetPoint( 2, div[2][3], 1.0 ); b1[3]->SetPoint( 3, div[2][3],  -1.0 );

  b1[0]->SetFillColor(kMagenta-10); b1[0]->SetLineColor(kMagenta-10); //b1[0]->SetFillStyle(3003); // ee
  b1[1]->SetFillColor(kMagenta-10); b1[1]->SetLineColor(kMagenta-10); //b1[1]->SetFillStyle(3003); // ee
  b1[2]->SetFillColor(kCyan   -10); b1[2]->SetLineColor(kCyan   -10); //b1[2]->SetFillStyle(3003); // mm
  b1[3]->SetFillColor(kCyan   -10); b1[3]->SetLineColor(kCyan   -10); //b1[3]->SetFillStyle(3003); // mm
  
  // MAKE WAKU(w)
  TH2D* w = new TH2D( "q2_AFB", "q^{2}-A_{FB}", 2, xmin, xmax, 2, ymin, ymax );
  w->GetXaxis()->CenterTitle();
  w->GetYaxis()->CenterTitle();
  w->SetXTitle( "q^{2} [GeV^{2}/c^{4}]" );
  w->SetYTitle( "A_{FB}" );

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TFile file_theory("theory/theory_curve.root");
  TGraph*       g_afb_true          = (TGraph*      )file_theory.Get( "afb_true"          ); 
  TGraphErrors* g_afb_bin_true_lep1 = (TGraphErrors*)file_theory.Get( "afb_bin_true_lep1" ); g_afb_bin_true_lep1->SetMarkerStyle(24); g_afb_bin_true_lep1->SetMarkerColor(2); g_afb_bin_true_lep1->SetLineColor(2);
  TGraphErrors* g_afb_bin_true_lep0 = (TGraphErrors*)file_theory.Get( "afb_bin_true_lep0" ); g_afb_bin_true_lep0->SetMarkerStyle(24); g_afb_bin_true_lep0->SetMarkerColor(3); g_afb_bin_true_lep0->SetLineColor(3);

  TGraphErrors* g_afb_bin_true = new TGraphErrors();
  g_afb_bin_true->SetPoint( 0, x[0][0], ( (2296013.082407-2835784.669695) + (2296013.082407-2835784.669695) )/(5131797.752257 + 5131797.752257) ); g_afb_bin_true->SetPointError( 0, xE[0][0], 0.0 );
  g_afb_bin_true->SetPoint( 1, x[0][1], ( (1190213.457180-952893.844314 ) + (1518830.569079-1165659.312802) )/(2143107.257026 + 2684489.837414) ); g_afb_bin_true->SetPointError( 1, xE[0][1], 0.0 );
  g_afb_bin_true->SetPoint( 2, x[0][2], ( (387902.495273 -199931.420275 ) + (731085.385110 -374995.259962 ) )/(587833.915548  + 1106080.645072) ); g_afb_bin_true->SetPointError( 2, xE[0][2], 0.0 );
  g_afb_bin_true->SetPoint( 3, x[0][3], ( (642843.998691 -276620.678992 ) + (642843.998691 -276620.678992 ) )/(919464.677682  + 919464.677682 ) ); g_afb_bin_true->SetPointError( 3, xE[0][3], 0.0 );
  g_afb_bin_true->SetMarkerStyle(24);
  g_afb_bin_true->SetMarkerColor(1);
  g_afb_bin_true->SetLineColor(1);

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  // DRAW (gt)
  TCanvas* c1 = Canvas( "c1","c1", 1, 1 );
  c1->Draw();
  
  c1->cd(1); // AFB(true)
  w->SetTitle("q^{2} - A_{FB}^{true}");
  w->DrawCopy();
  g_afb_true->Draw("Csame");
  g_afb_bin_true     ->Draw("Psame");
  //g_afb_bin_true_lep1->Draw("Psame");
  //g_afb_bin_true_lep0->Draw("Psame");
  for( Int_t i=0; i<4; i++ ){
    b1[i]->Draw("F");
    //b1[i]->Draw("L");
  }
  for( Int_t i=0; i<ng; i++ ) gt[i]->Draw("P");
  //for( Int_t i=0; i<4; i++ ){
  //l1[1][i]->Draw();
  //l1[2][i]->Draw();
  //}

  leg1->Draw();


  std::cout << "AFB(true,Xsll)" << std::endl;
  g_afb_bin_true->Print();
  
  // SAVE
  if( fl_save ) c1->Print( "pic/q2_AFB_kll.eps" );
  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();
  
  delete c1;
  delete gt;

  return 0;
}
