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
  const Double_t xmin =  0.0;
  //const Double_t xmax = 25.0;
  const Double_t xmax = 22.9;
  const Double_t ymin = -1.0;
  const Double_t ymax =  1.0;
  const Double_t div[ng][4] = {
    //{ (PDGmass::jpsi-0.325)*(PDGmass::jpsi-0.325), (PDGmass::jpsi+0.125)*(PDGmass::jpsi+0.125), (PDGmass::psi2s-0.20)*(PDGmass::psi2s-0.20), (PDGmass::psi2s+0.10)*(PDGmass::psi2s+0.10) }, // ll
    { (PDGmass::jpsi-0.250)*(PDGmass::jpsi-0.250), (PDGmass::jpsi+0.100)*(PDGmass::jpsi+0.100), (PDGmass::psi2s-0.15)*(PDGmass::psi2s-0.15), (PDGmass::psi2s+0.10)*(PDGmass::psi2s+0.10) }, // ll // modified @ 20130812
    { (PDGmass::jpsi-0.400)*(PDGmass::jpsi-0.400), (PDGmass::jpsi+0.150)*(PDGmass::jpsi+0.150), (PDGmass::psi2s-0.25)*(PDGmass::psi2s-0.25), (PDGmass::psi2s+0.10)*(PDGmass::psi2s+0.10) }, // ee
    { (PDGmass::jpsi-0.250)*(PDGmass::jpsi-0.250), (PDGmass::jpsi+0.100)*(PDGmass::jpsi+0.100), (PDGmass::psi2s-0.15)*(PDGmass::psi2s-0.15), (PDGmass::psi2s+0.10)*(PDGmass::psi2s+0.10) }, // mm
  };
  const Double_t x[ng][nafb] = {
    //{4.3/2.0, (4.3 + div[0][0])/2.0,  (div[0][1] + div[0][2])/2.0, (div[0][3] + xmax)/2.0}, // ll
    {4.3/2.0, (4.3 + div[2][0])/2.0,  (div[2][1] + div[2][2])/2.0, (div[2][3] + xmax)/2.0}, // ll // modified @ 20130812
    {4.3/2.0, (4.3 + div[1][0])/2.0,  (div[1][1] + div[1][2])/2.0, (div[1][3] + xmax)/2.0}, // ee
    {4.3/2.0, (4.3 + div[2][0])/2.0,  (div[2][1] + div[2][2])/2.0, (div[2][3] + xmax)/2.0}, // mm
  };
  const Double_t xE[ng][nafb] = {
    //{4.3/2.0, (-4.3 + div[0][0])/2.0,  (-div[0][1] + div[0][2])/2.0, (-div[0][3] + xmax)/2.0}, // ll
    {4.3/2.0, (-4.3 + div[2][0])/2.0,  (-div[2][1] + div[2][2])/2.0, (-div[2][3] + xmax)/2.0}, // ll // modified @ 20130812
    {4.3/2.0, (-4.3 + div[1][0])/2.0,  (-div[1][1] + div[1][2])/2.0, (-div[1][3] + xmax)/2.0}, // ee
    {4.3/2.0, (-4.3 + div[2][0])/2.0,  (-div[2][1] + div[2][2])/2.0, (-div[2][3] + xmax)/2.0}, // mm
  };

  const Double_t afb_true[ng][nafb] = {
    //{0.3422, 0.04095, 0.277,  0.2801}, // ee+mm(sim)
    {0.3422, 0.022,   0.2749, 0.2801}, // ee+mm(sim) // modified @ 20130812
    {0.215,  -0.090,  0.003,  0.089 }, // ee
    {0.633,  0.370,   0.588,  0.355 }, // mm
  };
  const Double_t afb_trueE[ng][nafb] = {
    //{0.2423, 0.3119, 0.2144, 0.1472}, // ee+mm(sim)
    {0.2423, 0.2738, 0.2138, 0.1472}, // ee+mm(sim) // modified @ 20130812
    {0.301,  0.379,  0.301,  0.270 }, // ee
    {0.448,  0.510,  0.297,  0.168 }, // mm
  };

  const Double_t afb_meas[ng][nafb] = {
    {0.239,  0.099,   0.294,  0.229 }, // ee+mm(sum)
    {0.167,  -0.079,  0.002,  0.079 }, // ee
    {0.304,  0.269,   0.569,  0.328 }, // mm
  };
  const Double_t afb_measE[ng][nafb] = {
    {0.156, 0.253, 0.201, 0.132 }, // ee+mm(sum)
    {0.233, 0.333, 0.283, 0.240 }, // ee
    {0.215, 0.371, 0.287, 0.155 }, // mm
  };
  
  //++++++++++++++++++++++++++++++ 

  // MAKE TGRAPHERRORS (gt and gm)
  TGraphErrors** gt = new TGraphErrors*[ng]; // observed A_FB(true) [ee+mm, ee, mm]
  TGraphErrors** gm = new TGraphErrors*[ng]; // raw      A_FB       [ee+mm, ee, mm]
  for( Int_t i=0; i<ng; i++ ){
    gt[i] = new TGraphErrors();
    gm[i] = new TGraphErrors();
    for( Int_t j=0; j<nafb; j++ ){
      gt[i]->SetPoint     ( j, x [i][j], afb_true [i][j] );
      gt[i]->SetPointError( j, xE[i][j], afb_trueE[i][j] );
      gm[i]->SetPoint     ( j, x [i][j], afb_meas [i][j] );
      gm[i]->SetPointError( j, xE[i][j], afb_measE[i][j] );
    }
    gt[i]->Sort();
    gm[i]->Sort();
    Deco( gt[i], 0, i+1, i+1 );
    Deco( gm[i], 0, i+1, i+1 );
    gt[i]->SetMarkerSize(0.8);
    gm[i]->SetMarkerSize(0.8);
  }

  // MAKE LEGEND(leg1)
  TLegend* leg1 = new TLegend( 0.75,0.75,0.99,0.99 );
  leg1->AddEntry( gt[0],"ll",     "PL" );
  leg1->AddEntry( gt[1],"ee",     "PL" );
  leg1->AddEntry( gt[2],"#mu#mu", "PL" );

  // MAKE BAND(b1) for charmonium veto region
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
  TGraph*       g_afb_obs           = (TGraph*      )file_theory.Get( "afb_obs"           ); 
  TGraph*       g_afb_obs_lep1      = (TGraph*      )file_theory.Get( "afb_obs_lep1"      ); 
  TGraph*       g_afb_obs_lep0      = (TGraph*      )file_theory.Get( "afb_obs_lep0"      ); 
  TGraphErrors* g_afb_bin_true_lep1 = (TGraphErrors*)file_theory.Get( "afb_bin_true_lep1" ); g_afb_bin_true_lep1->SetMarkerStyle(24); g_afb_bin_true_lep1->SetMarkerColor(2); g_afb_bin_true_lep1->SetLineColor(2);
  TGraphErrors* g_afb_bin_true_lep0 = (TGraphErrors*)file_theory.Get( "afb_bin_true_lep0" ); g_afb_bin_true_lep0->SetMarkerStyle(24); g_afb_bin_true_lep0->SetMarkerColor(3); g_afb_bin_true_lep0->SetLineColor(3);
  TGraphErrors* g_afb_bin_obs_lep1  = (TGraphErrors*)file_theory.Get( "afb_bin_obs_lep1"  ); g_afb_bin_obs_lep1 ->SetMarkerStyle(24); g_afb_bin_obs_lep1 ->SetMarkerColor(2); g_afb_bin_obs_lep1 ->SetLineColor(2);
  TGraphErrors* g_afb_bin_obs_lep0  = (TGraphErrors*)file_theory.Get( "afb_bin_obs_lep0"  ); g_afb_bin_obs_lep0 ->SetMarkerStyle(24); g_afb_bin_obs_lep0 ->SetMarkerColor(3); g_afb_bin_obs_lep0 ->SetLineColor(3);

  TGraphErrors* g_afb_bin_true = new TGraphErrors(); // obsolete (ee+mm)
  g_afb_bin_true->SetPoint( 0, x[0][0], ( (2296013.082407-2835784.669695) + (2296013.082407-2835784.669695) )/(5131797.752257 + 5131797.752257) ); g_afb_bin_true->SetPointError( 0, xE[0][0], 0.0 );
  g_afb_bin_true->SetPoint( 1, x[0][1], ( (1190213.457180-952893.844314 ) + (1518830.569079-1165659.312802) )/(2143107.257026 + 2684489.837414) ); g_afb_bin_true->SetPointError( 1, xE[0][1], 0.0 );
  g_afb_bin_true->SetPoint( 2, x[0][2], ( (387902.495273 -199931.420275 ) + (731085.385110 -374995.259962 ) )/(587833.915548  + 1106080.645072) ); g_afb_bin_true->SetPointError( 2, xE[0][2], 0.0 );
  g_afb_bin_true->SetPoint( 3, x[0][3], ( (642843.998691 -276620.678992 ) + (642843.998691 -276620.678992 ) )/(919464.677682  + 919464.677682 ) ); g_afb_bin_true->SetPointError( 3, xE[0][3], 0.0 );
  g_afb_bin_true->SetMarkerStyle(24);
  g_afb_bin_true->SetMarkerColor(1);
  g_afb_bin_true->SetLineColor(1);

  TGraphErrors* g_afb_bin_obs = new TGraphErrors();
  g_afb_bin_obs->SetPoint( 0, x[0][0], ( (59169.423526-67689.249445)+(47068.189006-51141.944206) )/( 126858.672971+98210.133212 ) ); g_afb_bin_obs->SetPointError( 0, xE[0][0], 0.0 );
  g_afb_bin_obs->SetPoint( 1, x[0][1], ( (42543.137167-35018.305515)+(50018.023367-40900.989584) )/( 77561.442683 +90919.012951 ) ); g_afb_bin_obs->SetPointError( 1, xE[0][1], 0.0 );
  g_afb_bin_obs->SetPoint( 2, x[0][2], ( (16167.768233-8589.334715 )+(39856.008816-20996.506917) )/( 24757.102948 +60852.515733 ) ); g_afb_bin_obs->SetPointError( 2, xE[0][2], 0.0 );
  g_afb_bin_obs->SetPoint( 3, x[0][3], ( (40573.845255-18875.899619)+(50352.042484-22732.457974) )/( 59449.744874 +73084.500459 ) ); g_afb_bin_obs->SetPointError( 3, xE[0][3], 0.0 );
  g_afb_bin_obs->SetMarkerStyle(24);
  g_afb_bin_obs->SetMarkerColor(1);
  g_afb_bin_obs->SetLineColor(1);

  
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  // DRAW (gt and gm)
  TCanvas* c1 = Canvas( "c1","c1", 2, 2 );
  /*  
  c1->Draw();
  c1->cd(1); // AFB(true)
  w->SetTitle("q^{2} - A_{FB}^{true}");
  w->DrawCopy();
  g_afb_true->Draw("Csame");
  //g_afb_bin_true     ->Draw("Psame");
  //g_afb_bin_true_lep1->Draw("Psame");
  g_afb_bin_true_lep0->Draw("Psame");
  for( Int_t i=0; i<4; i++ ){
    b1[i]->Draw("F");
    //b1[i]->Draw("L");
  }
  for( Int_t i=0; i<ng; i++ ) gt[i]->Draw("P");
  
  c1->cd(2); // AFB(obs,Xsll)
  w->SetTitle("q^{2} - A_{FB}^{obs} (X_{s} l^{+}l^{-})");
  w->DrawCopy();
  g_afb_obs    ->Draw("Csame");
  g_afb_bin_obs->Draw("Psame");
  b1[0]->Draw("F");
  b1[1]->Draw("F");
  b1[2]->Draw("F");
  b1[3]->Draw("F");
  gm[0]->Draw("P");


  c1->cd(3); // AFB(obs,Xsee)
  w->SetTitle("q^{2} - A_{FB}^{obs} (X_{s} e^{+}e^{-})");
  w->DrawCopy();
  g_afb_obs_lep1    ->Draw("Csame");
  g_afb_bin_obs_lep1->Draw("Psame");
  b1[0]->Draw("F");
  b1[1]->Draw("F");
  gm[1]->Draw("P");



  c1->cd(4); // AFB(obs,Xsmm)
  w->SetTitle("q^{2} - A_{FB}^{obs} (X_{s} #mu^{+}#mu^{-})");
  w->DrawCopy();
  g_afb_obs_lep0    ->Draw("Csame");
  g_afb_bin_obs_lep0->Draw("Psame");
  b1[2]->Draw("F");
  b1[3]->Draw("F");
  gm[2]->Draw("P");


  leg1->Draw();


  std::cout << "AFB(true,Xsll)" << std::endl;
  g_afb_bin_true->Print();
  std::cout << "AFB(obs,Xsll)" << std::endl;
  g_afb_bin_obs->Print();
  */

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++;

  TFile file_theory_err("theory_error/theory_curve_error.root");
  TGraph* g_afb_true_err = (TGraph*)file_theory_err.Get( "theory_curve_error_quadsum" );
  g_afb_true_err->SetFillColor(2);

  Double_t tmp_afb_bin[4] = { g_afb_bin_true_lep0->GetY()[0],
			      g_afb_bin_true_lep0->GetY()[1],
			      g_afb_bin_true_lep0->GetY()[3],
			      g_afb_bin_true_lep0->GetY()[5] };
  Double_t tmp_afb_bin_EH[4] = { 0.030, 0.034, 0.036, 0.044 };
  Double_t tmp_afb_bin_EL[4] = { 0.030, 0.034, 0.036, 0.037 };

  TGraph** theory_box = new TGraph*[4];
  const Int_t tmp_ind[4] = {0,1,3,5};
  for( Int_t i=0; i<4; i++ ){
    theory_box[i] = new TGraph();
    theory_box[i]->SetPoint( 0, q2_theta_nonuniform_eff::xbins_afb[0][tmp_ind[i]  ], tmp_afb_bin[i]-tmp_afb_bin_EL[i] ); // bottom left
    theory_box[i]->SetPoint( 1, q2_theta_nonuniform_eff::xbins_afb[0][tmp_ind[i]+1], tmp_afb_bin[i]-tmp_afb_bin_EL[i] ); // bottom right
    theory_box[i]->SetPoint( 2, q2_theta_nonuniform_eff::xbins_afb[0][tmp_ind[i]+1], tmp_afb_bin[i]+tmp_afb_bin_EH[i] ); // top    right
    theory_box[i]->SetPoint( 3, q2_theta_nonuniform_eff::xbins_afb[0][tmp_ind[i]  ], tmp_afb_bin[i]+tmp_afb_bin_EH[i] ); // top    left
    theory_box[i]->SetPoint( 4, q2_theta_nonuniform_eff::xbins_afb[0][tmp_ind[i]  ], tmp_afb_bin[i]-tmp_afb_bin_EL[i] ); // bottom left
    //theory_box[i]->SetFillColor(1);
    if( i==3 ){
      theory_box[i]->SetPoint( 1, xmax, tmp_afb_bin[i]-tmp_afb_bin_EL[i] ); // bottom right
      theory_box[i]->SetPoint( 2, xmax, tmp_afb_bin[i]+tmp_afb_bin_EH[i] ); // top    right
    }
    theory_box[i]->SetLineColor(1);
    theory_box[i]->SetLineStyle(2);
    //theory_box[i]->SetFillStyle(3003);
  }
  TGraph* cp_b1_2 = new TGraph();//(b1[2]);
  TGraph* cp_b1_3 = new TGraph();//(b1[3]);
  for( Int_t i=0; i<b1[2]->GetN(); i++ ) cp_b1_2->SetPoint( i, b1[2]->GetX()[i], b1[2]->GetY()[i] );
  for( Int_t i=0; i<b1[3]->GetN(); i++ ) cp_b1_3->SetPoint( i, b1[3]->GetX()[i], b1[3]->GetY()[i] );
  cp_b1_2->SetFillColor(1);
  cp_b1_3->SetFillColor(1);
  cp_b1_2->SetLineColor(1);
  cp_b1_3->SetLineColor(1);
  cp_b1_2->SetFillStyle(3004);
  cp_b1_3->SetFillStyle(3004);
  
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++;
  // for flipped C7 case
  TFile file_theory_flippedC7("theory/theory_curve_dev_-100_100_100.root");
  TGraph* g_afb_true_flippedC7 = (TGraph*)file_theory_flippedC7.Get( "afb_true" ); 
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++;

  TCanvas* c2 = Canvas( "c2","c2" );
  c2->Draw();
  w->DrawCopy(); // waku
  for( Int_t i=0; i<4; i++ ) b1[i]->Draw("F"); // charmonium veto [jpsi(ee),psi2s(ee),jpsi(mm),psi2s(mm)]
  cp_b1_2->Draw("F");
  cp_b1_3->Draw("F");
  g_afb_true_err      ->Draw("F"    ); // theory_curve_band
  g_afb_true          ->Draw("Csame"); // theory_curve
  g_afb_true_flippedC7->Draw("Csame"); // theory_curve (flipped C7 case)
  for( Int_t i=0; i<4; i++ ) theory_box[i]->Draw("L");
  //g_afb_bin_true_lep0->Draw("Psame"); // binned_theory(mm)
  gt[0]              ->Draw("P"    ); // our results(ee+mm)


  
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++;
  // SAVE
  if( fl_save ){
    c1->Print( "pic/q2_AFB.eps"     );
    c2->Print( "pic/q2_AFB_fin.eps" );
  }
  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();
  
  delete c1;
  delete c2;
  delete gt;
  delete gm;

  return 0;
}
