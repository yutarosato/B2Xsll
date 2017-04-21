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
#include <TLine.h>
#include <TF1.h>
#include <TLegend.h>

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++;

Int_t main( Int_t argc, Char_t** argv ){
  TApplication app( "app", &argc, argv );
  Style();
  //++++++++++++++++++++++++++++++
  if( !(argc==1 || argc==2) ) std::cerr << "wrong input" << std::endl, abort();
  Int_t         fl_appRun = 1;
  if( argc==2 ) fl_appRun = atoi( argv[1] );
  //++++++++++++++++++++++++++++++
  const Char_t*  file_gen       = "tmp_gen.dat";
  const Char_t*  file_rec       = "tmp_rec.dat";
  const Bool_t   fl_save    = true;
  const Bool_t   fl_message = true;
  const Double_t scale[2]   = {6.0,100.0}; // [gmc, cc-mc]
  //++++++++++++++++++++++++++++++
  TGraphErrors*** graph_gen  = new TGraphErrors**[2]; // [mc][mode]
  TGraphErrors*** graph_rec  = new TGraphErrors**[2]; // [mc][mode]
  TGraphErrors*** graph_eff  = new TGraphErrors**[2]; // [mc][mode]
  for( Int_t i=0; i<2; i++ ){
    graph_gen[i] = new TGraphErrors*[6];
    graph_rec[i] = new TGraphErrors*[6];
    graph_eff[i] = new TGraphErrors*[6];
    for( Int_t j=0; j<6; j++ ){
      graph_gen[i][j] = new TGraphErrors();
      graph_rec[i][j] = new TGraphErrors();
      graph_eff[i][j] = new TGraphErrors();
    }
  }
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Char_t buf[1024];
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // [GEN]
  Int_t cnt_gen[2][6][nexpno_gmc]={0};
  ifstream fin_gen;
  fin_gen.open(file_gen);
  
  while(!fin_gen.eof()){
    fin_gen.getline(buf,1024);
    if(buf[0]=='#') continue; // comment
    if(buf[0]=='*') break;    // finish
    std::istringstream sTmp(buf);
    Int_t mc, mode, exp, nevt;
    sTmp >> mc >> mode >> exp >> nevt;
    Int_t tmp_mc  = mc ? 1 : 0;
    Int_t tmp_exp = -1;
    for( Int_t e=0; e<nexpno_gmc; e++ ){
      if( exp == expno_gmc_int[e] ) tmp_exp = e;
    }
    if( tmp_exp == -1 ) std::cerr << "[ABORT] Wrong expNo : " << exp << std::endl, abort();
    std::cout << std::setw(3) << std::right << mc      << " -> "
	      << std::setw(3) << std::right << tmp_mc  << ", "
	      << std::setw(3) << std::right << mode    << ", "
	      << std::setw(3) << std::right << exp     << " -> "
	      << std::setw(3) << std::right << tmp_exp << ", "
	      << std::setw(8) << std::right << nevt    << ", "
	      << std::setw(8) << std::right << nbb_gmc[tmp_exp] << ", "
	      << std::setw(8) << std::right << (Double_t)nevt << " -> "
	      << std::setw(8) << std::right << (Double_t)nevt/scale[tmp_mc] << " -> "
	      << std::setw(8) << std::right << (Double_t)nevt/scale[tmp_mc]/nbb_gmc[tmp_exp] << ", "
	      << std::endl;
    
    cnt_gen[tmp_mc][mode][tmp_exp] += nevt;
  }
  fin_gen.close();
  
  for( Int_t i=0; i<2; i++ ){
    for( Int_t j=0; j<6; j++ ){
      for( Int_t e=0; e<nexpno_gmc; e++ ){
	graph_gen[i][j]->SetPoint     ( e, expno_gmc_int[e],      (Double_t)cnt_gen[i][j][e] /scale[i]/nbb_ccmc[e] );
	graph_gen[i][j]->SetPointError( e,                0, sqrt((Double_t)cnt_gen[i][j][e])/scale[i]/nbb_ccmc[e] );
      }
    }
  }
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // [REC]
  Int_t cnt_rec[2][6][nexpno_gmc]={0};
  ifstream fin_rec;
  fin_rec.open(file_rec);
  
  while(!fin_rec.eof()){
    fin_rec.getline(buf,1024);
    if(buf[0]=='#') continue; // comment
    if(buf[0]=='*') break;    // finish
    std::istringstream sTmp(buf);
    Int_t mc, mode, exp, nevt;
    sTmp >> mc >> mode >> exp >> nevt;
    Int_t tmp_mc  = mc ? 1 : 0;
    Int_t tmp_exp = -1;
    for( Int_t e=0; e<nexpno_gmc; e++ ){
      if( exp == expno_gmc_int[e] ) tmp_exp = e;
    }
    if( tmp_exp == -1 ) std::cerr << "[ABORT] Wrong expNo : " << exp << std::endl, abort();
    std::cout << std::setw(3) << std::right << mc      << " -> "
	      << std::setw(3) << std::right << tmp_mc  << ", "
	      << std::setw(3) << std::right << mode    << ", "
	      << std::setw(3) << std::right << exp     << " -> "
	      << std::setw(3) << std::right << tmp_exp << ", "
	      << std::setw(8) << std::right << nevt    << ", "
	      << std::setw(8) << std::right << nbb_gmc[tmp_exp] << ", "
	      << std::setw(8) << std::right << (Double_t)nevt << " -> "
	      << std::setw(8) << std::right << (Double_t)nevt/scale[tmp_mc] << " -> "
	      << std::setw(8) << std::right << (Double_t)nevt/scale[tmp_mc]/nbb_gmc[tmp_exp] << ", "
	      << std::endl;
    
    cnt_rec[tmp_mc][mode][tmp_exp] += nevt;
  }
  fin_rec.close();
  
  for( Int_t i=0; i<2; i++ ){
    for( Int_t j=0; j<6; j++ ){
      for( Int_t e=0; e<nexpno_gmc; e++ ){
	graph_rec[i][j]->SetPoint     ( e, expno_gmc_int[e],      (Double_t)cnt_rec[i][j][e] /scale[i]/nbb_ccmc[e] );
	graph_rec[i][j]->SetPointError( e,                0, sqrt((Double_t)cnt_rec[i][j][e])/scale[i]/nbb_ccmc[e] );
      }
    }
  }
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // EFF and total
  for( Int_t i=0; i<2; i++ ){
    for( Int_t j=0; j<6; j++ ){
      Double_t tot_gen = 0;
      Double_t tot_rec = 0;
      for( Int_t e=0; e<nexpno_gmc; e++ ){
	Double_t p = (Double_t)cnt_rec[i][j][e]/(Double_t)cnt_gen[i][j][e];
	graph_eff[i][j]->SetPoint     ( e, expno_gmc_int[e], p                                        );
	graph_eff[i][j]->SetPointError( e,                0, p*sqrt((1-p)/(Double_t)cnt_rec[i][j][e]) );
	tot_gen += cnt_gen[i][j][e];
	tot_rec += cnt_rec[i][j][e];
      }
      Double_t p = tot_rec/tot_gen;
      graph_eff[i][j]->SetPoint     ( nexpno_gmc, 75, p                     );
      graph_eff[i][j]->SetPointError( nexpno_gmc,  0, p*sqrt((1-p)/tot_rec) );
      Double_t tmp_scale =0;
      for( Int_t e=0; e<nexpno_gmc; e++ ){
	tmp_scale += scale[i]*nbb_ccmc[e];
      }
      graph_gen[i][j]->SetPoint     ( nexpno_gmc, 75, tot_gen/tmp_scale       );
      graph_gen[i][j]->SetPointError( nexpno_gmc,  0, sqrt(tot_gen)/tmp_scale );
      graph_rec[i][j]->SetPoint     ( nexpno_gmc, 75, tot_rec/tmp_scale       );
      graph_rec[i][j]->SetPointError( nexpno_gmc,  0, sqrt(tot_rec)/tmp_scale );
    }
  }

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TCanvas* c1 = Canvas( "c1","c1",6, 3 );
  for( Int_t j=0; j<6; j++ ){
    for( Int_t i=0; i<2; i++ ){
      graph_gen[i][j]->Sort();
      graph_gen[i][j]->SetMarkerStyle(20);
      graph_gen[i][j]->SetMarkerColor(i+2);
      graph_gen[i][j]->SetLineColor  (i+2);
      graph_gen[i][j]->SetMinimum(0);
      graph_gen[i][j]->SetTitle( Form("[Gen] mode%d",j) );
      graph_rec[i][j]->Sort();
      graph_rec[i][j]->SetMarkerStyle(20);
      graph_rec[i][j]->SetMarkerColor(i+2);
      graph_rec[i][j]->SetLineColor  (i+2);
      graph_rec[i][j]->SetMinimum(0);
      graph_rec[i][j]->SetTitle( Form("[Rec] mode%d",j) );
      graph_eff[i][j]->Sort();
      graph_eff[i][j]->SetMarkerStyle(20);
      graph_eff[i][j]->SetMarkerColor(i+2);
      graph_eff[i][j]->SetLineColor  (i+2);
      graph_eff[i][j]->SetMinimum(0);
      graph_eff[i][j]->SetTitle( Form("[Eff] mode%d",j) );
    }
    c1->cd(j+1);
    graph_gen[0][j]->Draw("AP");
    graph_gen[1][j]->Draw("Psame");
    c1->cd(j+1+6);
    graph_rec[0][j]->Draw("AP");
    graph_rec[1][j]->Draw("Psame");
    c1->cd(j+1+12);
    graph_eff[0][j]->Draw("AP");
    graph_eff[1][j]->Draw("Psame");
  }
  // Tlegend
  c1->cd(6);
  TLegend* legend1 = new TLegend( 0.75,0.75,0.99,0.99 );
  legend1->AddEntry( graph_gen[0][0], "gMC",   "PL" );
  legend1->AddEntry( graph_gen[1][0], "cc-MC", "PL" );
  legend1->Draw();

  // SAVE
  if( fl_save ) c1->Print( "pic/entry_eff.eps");
  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();
  return 0;
}

