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
#include <TLine.h>



Int_t main( Int_t argc, Char_t** argv ){
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==1 || argc==2) ) std::cerr << "wrong input" << std::endl, abort();
  Int_t    fl_appRun  = 1;
  if( argc==2 ) fl_appRun = atoi( argv[1] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t sel_fun = 310; // 10(gauss), 310(cb)
  const Double_t threshold = 0.683;
  TChain* c_sigmc = new TChain("h511");
  TChain* c_gmc   = new TChain("h511");
  TChain* c_ccmc  = new TChain("h511");
  std::cout << "[sigMC] " << c_sigmc->Add( "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut2/sig_522/sigMC_*.root" ) << " files" << std::endl;
  std::cout << "[ gMC ] " << c_gmc  ->Add( "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut2/bkg_522/gMC_*.root"   ) << " files" << std::endl;
  std::cout << "[CC-MC] " << c_ccmc ->Add( "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut2/cc_522/CC_*.root"     ) << " files" << std::endl;
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TH1D** hist = new TH1D*[20];
  for( Int_t i=0; i<20; i++ ) hist[i] = new TH1D( Form("hist%d",i), Form("hist%d",i), 100, 0.0, 4.8);
  hist[0]->SetLineColor(1); hist[0]->SetLineStyle(1);
  hist[1]->SetLineColor(2); hist[1]->SetLineStyle(1);
  hist[2]->SetLineColor(3); hist[2]->SetLineStyle(1);
  hist[3]->SetLineColor(4); hist[3]->SetLineStyle(1);
  
  hist[ 4]->SetLineColor(1);  hist[ 4]->SetLineStyle(1);
  hist[ 5]->SetLineColor(2);  hist[ 5]->SetLineStyle(1);
  hist[ 6]->SetLineColor(1);  hist[ 6]->SetLineStyle(1);
  hist[ 7]->SetLineColor(2);  hist[ 7]->SetLineStyle(1);
  hist[ 8]->SetLineColor(1);  hist[ 8]->SetLineStyle(1);
  hist[ 9]->SetLineColor(2);  hist[ 9]->SetLineStyle(1);
  hist[10]->SetLineColor(1);  hist[10]->SetLineStyle(1);
  hist[11]->SetLineColor(2);  hist[11]->SetLineStyle(1);

  hist[12]->SetLineColor(1);  hist[12]->SetLineStyle(1);
  hist[13]->SetLineColor(2);  hist[13]->SetLineStyle(1);
  hist[14]->SetLineColor(1);  hist[14]->SetLineStyle(1);
  hist[15]->SetLineColor(2);  hist[15]->SetLineStyle(1);
  hist[16]->SetLineColor(1);  hist[16]->SetLineStyle(1);
  hist[17]->SetLineColor(2);  hist[17]->SetLineStyle(1);
  hist[18]->SetLineColor(1);  hist[18]->SetLineStyle(1);
  hist[19]->SetLineColor(2);  hist[19]->SetLineStyle(1);
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TH1D** hist2 = new TH1D*[2];
  TF1**  func2 = new TF1* [2];
  for( Int_t i=0; i<2; i++ ){
    hist2[i] = new TH1D( Form("hist_fit%d",i), Form("hist_fit%d",i), 100, 2.5, 3.5);
    func2[i] = new TF1 ( Form("func%d",i), make_func(sel_fun), 2.5, 3.5, n_fitfunc_par(sel_fun) );
    func2[i]->SetLineColor(2);
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  c_sigmc->Project( "hist0",  "cc_m",    "(rm_xs==1001||rm_xs==1010||rm_xs==1101||rm_xs==1110)&&Mbc>5.27&&rm_l==1&&self==1" );
  c_sigmc->Project( "hist1",  "cc_mheg", "(rm_xs==1001||rm_xs==1010||rm_xs==1101||rm_xs==1110)&&Mbc>5.27&&rm_l==1&&self==1" );
  c_sigmc->Project( "hist2",  "cc_mheg", "(rm_xs==1001||rm_xs==1010||rm_xs==1101||rm_xs==1110)&&Mbc>5.27&&rm_l==1&&self==1&&cc_m<sqrt(4.3)"           );
  c_sigmc->Project( "hist3",  "cc_mheg", "(rm_xs==1001||rm_xs==1010||rm_xs==1101||rm_xs==1110)&&Mbc>5.27&&rm_l==1&&self==1&&cc_m>sqrt(4.3)&&cc_m<3.0" );

  c_gmc  ->Project( "hist4",  "cc_m",    "(rm_xs==1001||rm_xs==1010||rm_xs==1101||rm_xs==1110)&&Mbc>5.27&&rm_l==1&&lpmoid==443&&lmmoid==443&&lporg==lmorg&&lporg==korg&&lporg==abs(pi1org)&&lporg==abs(pi2org)&&lporg==abs(pi3org)&&lporg==abs(pi4org)&&rest_sw==0&&pi0org==1"  );
  c_gmc  ->Project( "hist5",  "cc_mheg", "(rm_xs==1001||rm_xs==1010||rm_xs==1101||rm_xs==1110)&&Mbc>5.27&&rm_l==1&&lpmoid==443&&lmmoid==443&&lporg==lmorg&&lporg==korg&&lporg==abs(pi1org)&&lporg==abs(pi2org)&&lporg==abs(pi3org)&&lporg==abs(pi4org)&&rest_sw==0&&pi0org==1"  );
  c_gmc  ->Project( "hist6",  "cc_m",    "(rm_xs==1001||rm_xs==1010||rm_xs==1101||rm_xs==1110)&&Mbc>5.27&&rm_l==1&&lpmoid==443&&lmmoid==443&&lporg==lmorg&&lporg==korg&&lporg==abs(pi1org)&&lporg==abs(pi2org)&&lporg==abs(pi3org)&&lporg==abs(pi4org)&&rest_sw==0&&pi0org==0"  );
  c_gmc  ->Project( "hist7",  "cc_mheg", "(rm_xs==1001||rm_xs==1010||rm_xs==1101||rm_xs==1110)&&Mbc>5.27&&rm_l==1&&lpmoid==443&&lmmoid==443&&lporg==lmorg&&lporg==korg&&lporg==abs(pi1org)&&lporg==abs(pi2org)&&lporg==abs(pi3org)&&lporg==abs(pi4org)&&rest_sw==0&&pi0org==0"  );
  c_gmc  ->Project( "hist8",  "cc_m",    "(rm_xs==1001||rm_xs==1010||rm_xs==1101||rm_xs==1110)&&Mbc>5.27&&rm_l==1&&lpmoid==443&&lmmoid==443&&lporg==lmorg&&lporg==korg&&lporg==abs(pi1org)&&lporg==abs(pi2org)&&lporg==abs(pi3org)&&lporg==abs(pi4org)&&rest_sw==0&&pi0org==-1&& (heg_self==1||leg_self==1)" );
  c_gmc  ->Project( "hist9",  "cc_mheg", "(rm_xs==1001||rm_xs==1010||rm_xs==1101||rm_xs==1110)&&Mbc>5.27&&rm_l==1&&lpmoid==443&&lmmoid==443&&lporg==lmorg&&lporg==korg&&lporg==abs(pi1org)&&lporg==abs(pi2org)&&lporg==abs(pi3org)&&lporg==abs(pi4org)&&rest_sw==0&&pi0org==-1&& (heg_self==1||leg_self==1)" );
  c_gmc  ->Project( "hist10", "cc_m",    "(rm_xs==1001||rm_xs==1010||rm_xs==1101||rm_xs==1110)&&Mbc>5.27&&rm_l==1&&lpmoid==443&&lmmoid==443&&lporg==lmorg&&lporg==korg&&lporg==abs(pi1org)&&lporg==abs(pi2org)&&lporg==abs(pi3org)&&lporg==abs(pi4org)&&rest_sw==0&&pi0org==-1&&!(heg_self==1||leg_self==1)" );
  c_gmc  ->Project( "hist11", "cc_mheg", "(rm_xs==1001||rm_xs==1010||rm_xs==1101||rm_xs==1110)&&Mbc>5.27&&rm_l==1&&lpmoid==443&&lmmoid==443&&lporg==lmorg&&lporg==korg&&lporg==abs(pi1org)&&lporg==abs(pi2org)&&lporg==abs(pi3org)&&lporg==abs(pi4org)&&rest_sw==0&&pi0org==-1&&!(heg_self==1||leg_self==1)" );

  c_ccmc  ->Project( "hist12", "cc_m",    "(rm_xs==1001||rm_xs==1010||rm_xs==1101||rm_xs==1110)&&Mbc>5.27&&rm_l==1&&lpmoid==443&&lmmoid==443&&lporg==lmorg&&lporg==korg&&lporg==abs(pi1org)&&lporg==abs(pi2org)&&lporg==abs(pi3org)&&lporg==abs(pi4org)&&rest_sw==0&&pi0org==1"  );
  c_ccmc  ->Project( "hist13", "cc_mheg", "(rm_xs==1001||rm_xs==1010||rm_xs==1101||rm_xs==1110)&&Mbc>5.27&&rm_l==1&&lpmoid==443&&lmmoid==443&&lporg==lmorg&&lporg==korg&&lporg==abs(pi1org)&&lporg==abs(pi2org)&&lporg==abs(pi3org)&&lporg==abs(pi4org)&&rest_sw==0&&pi0org==1"  );
  c_ccmc  ->Project( "hist14", "cc_m",    "(rm_xs==1001||rm_xs==1010||rm_xs==1101||rm_xs==1110)&&Mbc>5.27&&rm_l==1&&lpmoid==443&&lmmoid==443&&lporg==lmorg&&lporg==korg&&lporg==abs(pi1org)&&lporg==abs(pi2org)&&lporg==abs(pi3org)&&lporg==abs(pi4org)&&rest_sw==0&&pi0org==0"  );
  c_ccmc  ->Project( "hist15", "cc_mheg", "(rm_xs==1001||rm_xs==1010||rm_xs==1101||rm_xs==1110)&&Mbc>5.27&&rm_l==1&&lpmoid==443&&lmmoid==443&&lporg==lmorg&&lporg==korg&&lporg==abs(pi1org)&&lporg==abs(pi2org)&&lporg==abs(pi3org)&&lporg==abs(pi4org)&&rest_sw==0&&pi0org==0"  );
  c_ccmc  ->Project( "hist16", "cc_m",    "(rm_xs==1001||rm_xs==1010||rm_xs==1101||rm_xs==1110)&&Mbc>5.27&&rm_l==1&&lpmoid==443&&lmmoid==443&&lporg==lmorg&&lporg==korg&&lporg==abs(pi1org)&&lporg==abs(pi2org)&&lporg==abs(pi3org)&&lporg==abs(pi4org)&&rest_sw==0&&pi0org==-1&& (heg_self==1||leg_self==1)" );
  c_ccmc  ->Project( "hist17", "cc_mheg", "(rm_xs==1001||rm_xs==1010||rm_xs==1101||rm_xs==1110)&&Mbc>5.27&&rm_l==1&&lpmoid==443&&lmmoid==443&&lporg==lmorg&&lporg==korg&&lporg==abs(pi1org)&&lporg==abs(pi2org)&&lporg==abs(pi3org)&&lporg==abs(pi4org)&&rest_sw==0&&pi0org==-1&& (heg_self==1||leg_self==1)" );
  c_ccmc  ->Project( "hist18", "cc_m",    "(rm_xs==1001||rm_xs==1010||rm_xs==1101||rm_xs==1110)&&Mbc>5.27&&rm_l==1&&lpmoid==443&&lmmoid==443&&lporg==lmorg&&lporg==korg&&lporg==abs(pi1org)&&lporg==abs(pi2org)&&lporg==abs(pi3org)&&lporg==abs(pi4org)&&rest_sw==0&&pi0org==-1&&!(heg_self==1||leg_self==1)" );
  c_ccmc  ->Project( "hist19", "cc_mheg", "(rm_xs==1001||rm_xs==1010||rm_xs==1101||rm_xs==1110)&&Mbc>5.27&&rm_l==1&&lpmoid==443&&lmmoid==443&&lporg==lmorg&&lporg==korg&&lporg==abs(pi1org)&&lporg==abs(pi2org)&&lporg==abs(pi3org)&&lporg==abs(pi4org)&&rest_sw==0&&pi0org==-1&&!(heg_self==1||leg_self==1)" );

  c_gmc  ->Project( "hist_fit0", "cc_mheg", "(rm_xs==1001||rm_xs==1010||rm_xs==1101||rm_xs==1110)&&Mbc>5.27&&rm_l==1&&lpmoid==443&&lmmoid==443&&lporg==lmorg&&lporg==korg&&lporg==abs(pi1org)&&lporg==abs(pi2org)&&lporg==abs(pi3org)&&lporg==abs(pi4org)&&rest_sw==0&&pi0org==-1&& (heg_self==1||leg_self==1)" );
  c_ccmc ->Project( "hist_fit1", "cc_mheg", "(rm_xs==1001||rm_xs==1010||rm_xs==1101||rm_xs==1110)&&Mbc>5.27&&rm_l==1&&lpmoid==443&&lmmoid==443&&lporg==lmorg&&lporg==korg&&lporg==abs(pi1org)&&lporg==abs(pi2org)&&lporg==abs(pi3org)&&lporg==abs(pi4org)&&rest_sw==0&&pi0org==-1&& (heg_self==1||leg_self==1)" );
  //(x1.2, -2.97),+0.02

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  TCanvas* c1 = Canvas("c1","c1",4,3);
  c1->Draw();

  c1->cd(1);
  hist[4]->SetTitle("[gMC] pi0org==1");
  hist[4]->Draw();
  hist[5]->Draw("same");
  c1->cd(2);
  hist[6]->SetTitle("[gMC] pi0org==0");
  hist[6]->Draw();
  hist[7]->Draw("same");
  c1->cd(3);
  hist[9]->SetTitle("[gMC] pi0org==-1&&[hl]eg_self==1");
  hist[9]->Draw();
  hist[8]->Draw("same");
  c1->cd(4);
  hist[10]->SetTitle("[gMC] pi0org==-1&&[hl]eg_self!=1");
  hist[10]->Draw();
  hist[11]->Draw("same");

  c1->cd(5);
  hist[12]->SetTitle("[ccMC] pi0org==1");
  hist[12]->Draw();
  hist[13]->Draw("same");
  c1->cd(6);
  hist[14]->SetTitle("[ccMC] pi0org==0");
  hist[14]->Draw();
  hist[15]->Draw("same");
  c1->cd(7);
  hist[17]->SetTitle("[ccMC] pi0org==-1&&[hl]eg_self==1");
  hist[17]->Draw();
  hist[16]->Draw("same");
  c1->cd(8);
  hist[18]->SetTitle("[ccMC] pi0org==-1&&[hl]eg_self!=1");
  hist[18]->Draw();
  hist[19]->Draw("same");

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  c1->cd(9);
  hist[0]->SetTitle("[sigMC]");
  hist[0]->Draw();
  hist[1]->Draw("same");
  hist[2]->Draw("same");
  hist[3]->Draw("same");
  
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  c1->cd(10);
  func_set_parameters(sel_fun, func2[0], hist2[0], 100, 2.5, 3.5);
  hist2[0]->SetTitle("[gMC]");
  hist2[0]->Fit("func0");
  hist2[0]->Draw();

  c1->cd(11);
  func_set_parameters(sel_fun, func2[1], hist2[1], 100, 2.5, 3.5);
  hist2[1]->SetTitle("[ccMC]");
  hist2[1]->Fit("func1");


  Double_t mean         = func2[1]->GetParameter(3);
  Double_t entry_high   = func2[1]->Integral(mean, 3.5 );
  Double_t entry_low    = func2[1]->Integral(2.5,  mean);
  Double_t step = 0.01;
  Double_t low_edge, high_edge;

  TLine** l = new TLine*[3];
  l[0] = new TLine(mean, 0.0, mean, 0.8*hist2[1]->GetMaximum() );
  // high side scan
  for( Int_t m=1; m<=100; m++ ){
    if( func2[1]->Integral(mean, mean+m*step)/entry_high > threshold ){
      high_edge = mean+m*step;
      l[1] = new TLine(mean+m*step, 0.0, mean+m*step, 0.8*hist2[1]->GetMaximum() );
      break;
    }
  }
  // low side scan
  for( Int_t m=1; m<=100; m++ ){
    if( func2[1]->Integral(mean-m*step, mean)/entry_low > threshold ){
      low_edge = mean-m*step;
      l[2] = new TLine(mean-m*step, 0.0, mean-m*step, 0.8*hist2[1]->GetMaximum() );
      break;
    }
  }

  l[0]->SetLineColor(3); l[0]->Draw();
  l[1]->SetLineColor(3); l[1]->Draw();
  l[2]->SetLineColor(3); l[2]->Draw();

  TPaveText* pave = new TPaveText( 0.05, 0.75, 0.40, 0.90,"BRNDC" );
  
  pave->AddText( Form("Low  = %f (%f%)", low_edge,  threshold ) );
  pave->AddText( Form("High = %f (%f%)", high_edge, threshold ) );
  pave->Draw();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  
  c1->Update();
  c1->Print("pic/M_llpi0.eps");

  
  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();
  
  return 0;
}
