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
#include <TLatex.h>
#include <TPaveText.h>

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++;

Int_t main( Int_t argc, Char_t** argv ){
  TApplication app( "app", &argc, argv );
  Style();
  //++++++++++++++++++++++++++++++
  if( !(argc==2 || argc==3) ) std::cerr << "wrong input" << std::endl, abort();
  const Int_t   tmp_i     = atoi(argv[1]); // 0(ee), 1(mm), 2(ee+mm)
  Int_t         fl_appRun = 1;
  if( argc==3 ) fl_appRun = atoi( argv[2] );
  //++++++++++++++++++++++++++++++
  Char_t* file_ntot  = new Char_t[256];
  Char_t* file_nq2   = new Char_t[256];
  Char_t* file_afb   = new Char_t[256];
  strcpy( file_ntot, Form("tmp%d_ntot.dat", tmp_i) );
  strcpy( file_nq2,  Form("tmp%d_nq2.dat",  tmp_i) );
  strcpy( file_afb,  Form("tmp%d_afb.dat",  2    ) );
  const Bool_t   fl_message = true;
  const Bool_t   fl_save    = true;
  const Int_t    nbin_q2    = 4; // 4 or 7 // # of pad should be adjusted.
  const Double_t xmin       = -0.5;
  const Double_t xmax       =  9.5;
  //++++++++++++++++++++++++++++++
  TGraphErrors*  g_ntot = new TGraphErrors();
  TGraphErrors** g_nq2  = new TGraphErrors*[nbin_q2];
  TGraphErrors** g_afb  = new TGraphErrors*[nbin_q2];
  TPaveText**    t_afb  = new TPaveText*   [nbin_q2];

  TF1*  f_ntot      = new TF1( "f_ntot", "[0]", xmin, xmax );
  TF1** f_nq2       = new TF1*[nbin_q2];
  TF1** f_afb       = new TF1*[nbin_q2];
  TGraph*  f_ntot_e = new TGraph();
  TGraph** f_nq2_e  = new TGraph*();
  TGraph** f_afb_e  = new TGraph*();

  f_ntot->SetTitle( "N_{sig}(total)" );
  f_ntot->SetParNames ( "offset" );
  f_ntot->SetParameter( 0, 300);
  f_ntot->SetLineColor(2);
  f_ntot->SetLineStyle(7);

  for( Int_t i=0; i<nbin_q2; i++ ){
    g_nq2  [i] = new TGraphErrors();
    g_afb  [i] = new TGraphErrors();
    f_nq2_e[i] = new TGraph();
    f_afb_e[i] = new TGraph();

    f_nq2[i] = new TF1( Form("f_nq2%d",i), "[0]", xmin, xmax );
    f_afb[i] = new TF1( Form("f_afb%d",i), "[0]", xmin, xmax );
    f_nq2[i]->SetTitle( Form( "N_{q2}, (%d q^{2}-bin)", i+1) );
    f_nq2[i]->SetParNames ( "offset" );
    f_nq2[i]->SetParameter( 0, 70 );
    f_nq2[i]->SetLineColor(2);
    f_nq2[i]->SetLineStyle(7);
    f_afb[i]->SetTitle( Form( "A_{FB}, (%d q^{2}-bin)", i+1) );
    f_afb[i]->SetParNames  ( "offset" );
    f_afb[i]->SetParameter( 0, 0 );
    f_afb[i]->SetLineColor(2);
    f_afb[i]->SetLineStyle(7);
  }
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TLine** l_afb_gen  = new TLine*[nbin_q2];
  TLine** l_afb_eff  = new TLine*[nbin_q2];
  TLine** l_nq2_exp  = new TLine*[nbin_q2];
  if( tmp_i==0 ){ // ee
    l_nq2_exp[0] = new TLine( xmin, 58.4,            xmax, 58.4            );
    l_nq2_exp[1] = new TLine( xmin, 34.2,            xmax, 34.2            );
    l_nq2_exp[2] = new TLine( xmin, 13.5,            xmax, 13.5            );
    l_nq2_exp[3] = new TLine( xmin, 34.6+5.667+6.33, xmax, 34.6+5.667+6.33 ); // include higher resonance
  }else if( tmp_i==1 ){ // mm
    l_nq2_exp[0] = new TLine( xmin, 43.7,             xmax, 43.7             );
    l_nq2_exp[1] = new TLine( xmin, 42.5,             xmax, 42.5             );
    l_nq2_exp[2] = new TLine( xmin, 32.3,             xmax, 32.3             );
    l_nq2_exp[3] = new TLine( xmin, 45.2+7.333+7.167, xmax, 45.2+7.333+7.167 ); // include higher resonance
  }

  for( Int_t i=0; i<nbin_q2; i++ ){
    l_afb_gen[i] = new TLine(xmin, AFB_gen_6qbin[2][i], xmax, AFB_gen_6qbin[2][i]);
    l_afb_gen[i]->SetLineColor(3);
    l_afb_gen[i]->SetLineStyle(7);
    l_afb_eff[i] = new TLine(xmin, AFB_eff_6qbin[2][i], xmax, AFB_eff_6qbin[2][i]);
    l_afb_eff[i]->SetLineColor(4);
    l_afb_eff[i]->SetLineStyle(7);
    l_nq2_exp[i]->SetLineColor(3);
    l_nq2_exp[i]->SetLineStyle(7);
  }
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Char_t buf[1024];
  // [Ntot]
  ifstream fin_ntot;
  fin_ntot.open(file_ntot);
  Int_t cnt_ntot=0;
  if( fl_message ) std::cout << "[READ DATA(Ntot)]" << std::endl;
  while(!fin_ntot.eof()){
    fin_ntot.getline(buf,1024);
    if(buf[0]=='#') continue; // comment
    if(buf[0]=='*') break;    // finish
    std::istringstream sTmp(buf);
    Double_t stream, x, xE;
    Int_t    q2bin;
    sTmp >> stream >> x >> xE;
    if( fl_message ) std::cout << std::setw(5) << std::right << Form("cnt = %d, ", cnt_ntot)
			       << std::setw(8) << std::right << Form("stream = %f, ", stream)
			       << std::setw(8) << std::right << Form("x =%f ", x )
			       << std::setw(8) << std::right << Form("+- %f",  xE)
			       << std::endl;
    g_ntot->SetPoint     ( cnt_ntot, stream, x  );
    g_ntot->SetPointError( cnt_ntot,      0, xE );
    cnt_ntot++;
  }
  fin_ntot.close();
  // [Nq2]
  ifstream fin_nq2;
  fin_nq2.open(file_nq2);
  Int_t cnt_nq2[nbin_q2] = {0};
  if( fl_message ) std::cout << "[READ DATA(Nq2)]" << std::endl;
  while(!fin_nq2.eof()){
    fin_nq2.getline(buf,1024);
    if(buf[0]=='#') continue; // comment
    if(buf[0]=='*') break;    // finish
    std::istringstream sTmp(buf);
    Double_t stream, x, xE;
    Int_t    q2bin;
    sTmp >> stream >> q2bin >> x >> xE;
    if( fl_message ) std::cout << std::setw(5) << std::right << Form("cnt = %d, ", cnt_nq2[q2bin])
			       << std::setw(8) << std::right << Form("stream = %f, ", stream)
			       << std::setw(8) << std::right << Form("q2 bin = %d, ", q2bin)
			       << std::setw(8) << std::right << Form("x =%f ", x )
			       << std::setw(8) << std::right << Form("+- %f",  xE)
			       << std::endl;
    g_nq2[q2bin]->SetPoint     ( cnt_nq2[q2bin], stream, x  );
    g_nq2[q2bin]->SetPointError( cnt_nq2[q2bin],      0, xE );
    cnt_nq2[q2bin]++;
  }
  fin_nq2.close();

  // [AFB]
  ifstream fin_afb;
  fin_afb.open(file_afb);
  Int_t cnt_afb[nbin_q2] = {0};
  if( fl_message ) std::cout << "[READ DATA(AFB)]" << std::endl;
  while(!fin_afb.eof()){
    fin_afb.getline(buf,1024);
    if(buf[0]=='#') continue; // comment
    if(buf[0]=='*') break;    // finish
    std::istringstream sTmp(buf);
    Double_t stream, x, xE;
    Int_t    q2bin;
    sTmp >> stream >> q2bin >> x >> xE;
    if( fl_message ) std::cout << std::setw(5) << std::right << Form("cnt = %d, ", cnt_afb[q2bin])
			       << std::setw(8) << std::right << Form("stream = %f, ", stream)
			       << std::setw(8) << std::right << Form("q2 bin = %d, ", q2bin)
			       << std::setw(8) << std::right << Form("x =%f ", x )
			       << std::setw(8) << std::right << Form("+- %f",  xE)
			       << std::endl;
    g_afb[q2bin]->SetPoint     ( cnt_afb[q2bin], stream, x  );
    g_afb[q2bin]->SetPointError( cnt_afb[q2bin],      0, xE );
    cnt_afb[q2bin]++;
  }
  fin_afb.close();
  
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  TCanvas* c1 = Canvas( "c1","c1",4, 3 );
  for( Int_t i=0; i<nbin_q2; i++ ){
    c1->cd(i+1);
    g_nq2[i]->Sort();
    g_nq2[i]->SetTitle( Form("N_{sig} (%d q^{2} bin)",i+1) );
    g_nq2[i]->GetXaxis()->SetTitle("stream No." );
    g_nq2[i]->GetYaxis()->SetTitle("N_{sig} (output)");
    g_nq2[i]->GetXaxis()->CenterTitle();
    g_nq2[i]->GetYaxis()->CenterTitle();
    g_nq2[i]->GetXaxis()->SetLimits( xmin, xmax);
    g_nq2[i]->SetMarkerStyle(20);
    g_nq2[i]->SetMarkerSize(1.2);
    g_nq2[i]->Draw("AP");
    g_nq2[i]->Fit( Form("f_nq2%d",i) );
    f_nq2_e[i]->SetPoint(0, xmin, f_nq2[i]->GetParameter(0)-f_nq2[i]->GetParError(0) );
    f_nq2_e[i]->SetPoint(1, xmin, f_nq2[i]->GetParameter(0)+f_nq2[i]->GetParError(0) );
    f_nq2_e[i]->SetPoint(2, xmax, f_nq2[i]->GetParameter(0)+f_nq2[i]->GetParError(0) );
    f_nq2_e[i]->SetPoint(3, xmax, f_nq2[i]->GetParameter(0)-f_nq2[i]->GetParError(0) );
    f_nq2_e[i]->SetFillColor(9);
    //f_nq2_e[i]->SetFillStyle(3003);
    f_nq2_e[i]->SetLineColor(0);
    f_nq2_e[i]->Draw("F"); // !!!
    f_nq2_e[i]->Draw("L"); // !!!
    l_nq2_exp[i]->Draw("same");

    c1->cd(i+1+nbin_q2);
    g_afb[i]->Sort();
    g_afb[i]->SetTitle( Form("A_{FB} (%d q^{2} bin)", i+1) );
    g_afb[i]->GetXaxis()->SetTitle("stream No." );
    g_afb[i]->GetYaxis()->SetTitle("A_{FB}(output)" );
    g_afb[i]->GetXaxis()->CenterTitle();
    g_afb[i]->GetYaxis()->CenterTitle();
    g_afb[i]->GetXaxis()->SetLimits( xmin, xmax);
    g_afb[i]->SetMarkerStyle(20);
    g_afb[i]->SetMarkerSize(1.2);
    g_afb[i]->Draw("AP");
    g_afb[i]->Fit( Form("f_afb%d",i) );
    f_afb_e[i]->SetPoint(0, xmin, f_afb[i]->GetParameter(0)-f_afb[i]->GetParError(0) );
    f_afb_e[i]->SetPoint(1, xmin, f_afb[i]->GetParameter(0)+f_afb[i]->GetParError(0) );
    f_afb_e[i]->SetPoint(2, xmax, f_afb[i]->GetParameter(0)+f_afb[i]->GetParError(0) );
    f_afb_e[i]->SetPoint(3, xmax, f_afb[i]->GetParameter(0)-f_afb[i]->GetParError(0) );
    f_afb_e[i]->SetFillColor(9);
    //f_afb_e[i]->SetFillStyle(3003);
    f_afb_e[i]->SetLineColor(0);
    f_afb_e[i]->Draw("F"); // !!!
    f_afb_e[i]->Draw("L"); // !!!
    l_afb_gen[i]->Draw("same");
    l_afb_eff[i]->Draw("same");
    t_afb[i] = new TPaveText( 0.05, 0.75, 0.40, 0.90,"BRNDC" );
    if( nbin_q2 == 4 ){
      t_afb[i]->AddText( Form("A_{FB}^{gen} = %f", AFB_gen_6qbin[2][i]) );
      t_afb[i]->AddText( Form("A_{FB}^{eff} = %f", AFB_eff_6qbin[2][i]) );
    }else if( nbin_q2 == 7 ){
      t_afb[i]->AddText( Form("A_{FB}^{gen} = %f", AFB_gen_9qbin[2][i]) );
      t_afb[i]->AddText( Form("A_{FB}^{eff} = %f", AFB_eff_9qbin[2][i]) );
    }
    t_afb[i]->Draw();
  }

  c1->cd(2*nbin_q2+1);
  g_ntot->Sort();
  g_ntot->SetTitle( "N_{sig}(total)" );
  g_ntot->GetXaxis()->SetTitle("stream No." );
  g_ntot->GetYaxis()->SetTitle("N_{sig}(total,output)");
  g_ntot->GetXaxis()->CenterTitle();
  g_ntot->GetYaxis()->CenterTitle();
  g_ntot->GetXaxis()->SetLimits( xmin, xmax);
  g_ntot->SetMarkerStyle(20);
  g_ntot->SetMarkerSize(1.2);
  g_ntot->Draw("AP");
  g_ntot->Fit("f_ntot");
  f_ntot_e->SetPoint(0, xmin, f_ntot->GetParameter(0)-f_ntot->GetParError(0) );
  f_ntot_e->SetPoint(1, xmin, f_ntot->GetParameter(0)+f_ntot->GetParError(0) );
  f_ntot_e->SetPoint(2, xmax, f_ntot->GetParameter(0)+f_ntot->GetParError(0) );
  f_ntot_e->SetPoint(3, xmax, f_ntot->GetParameter(0)-f_ntot->GetParError(0) );
  f_ntot_e->SetFillColor(9);
  //f_ntot_e->SetFillStyle(3003);
  f_ntot_e->SetLineColor(0);
  f_ntot_e->Draw("F"); // !!!
  f_ntot_e->Draw("L"); // !!!
    
  // Tlegend
  c1->cd(2*nbin_q2+2);
  TLegend* legend1 = new TLegend( 0.75,0.75,0.99,0.99 );
  legend1->AddEntry( l_afb_gen[0], "Gen",           "L" );
  legend1->AddEntry( l_afb_eff[0], "Eff",           "L" );
  legend1->AddEntry( f_afb[0],     "Rec(Fit->Gen)", "L" );
  legend1->Draw();

  // SAVE
  if( fl_save ){
    c1->Print( Form("pic/linearity_check_bin2_ensemble_%d_sim.eps",tmp_i) );
    c1->Print( Form("pic/linearity_check_bin2_ensemble_%d_sim.pdf",tmp_i) );
  }
  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();

  return 0;
}

