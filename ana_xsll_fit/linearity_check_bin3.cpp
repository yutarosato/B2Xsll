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
  const Bool_t  fl_message = true;
  const Bool_t  fl_save    = true;
  const Int_t   nbin_q2    = 4; // 4 or 7 // # of pad should be adjusted.
  const Double_t incl_asym_min =  0.0;
  const Double_t incl_asym_max = 10.0;
  //++++++++++++++++++++++++++++++
  TGraphErrors*  g_ntot = new TGraphErrors();
  TGraphErrors** g_nq2  = new TGraphErrors*[nbin_q2];
  TGraphErrors** g_afb  = new TGraphErrors*[nbin_q2];

  TF1*  f_ntot = new TF1( "f_ntot", "[0]", incl_asym_min, incl_asym_max );
  TF1** f_nq2  = new TF1*[nbin_q2];
  TF1** f_afb  = new TF1*[nbin_q2];

  f_ntot->SetTitle( "N_{sig}(total)" );
  f_ntot->SetParNames  ( "offset" );
  f_ntot->SetParameter( 0, 300 );
  f_ntot->SetLineColor(2);
  f_ntot->SetLineStyle(7);

  for( Int_t i=0; i<nbin_q2; i++ ){
    g_nq2[i] = new TGraphErrors();
    g_afb[i] = new TGraphErrors();

    f_nq2[i] = new TF1( Form("f_nq2%d",i), "[0]",       incl_asym_min, incl_asym_max );
    f_afb[i] = new TF1( Form("f_afb%d",i), "[0]*x+[1]", incl_asym_min, incl_asym_max );
    f_nq2[i]->SetTitle( Form( "N_{q2}, (%d q^{2}-bin)", i+1) );
    f_nq2[i]->SetParNames  ( "offset" );
    f_nq2[i]->SetParameter( 0, 300 );
    f_nq2[i]->SetLineColor(2);
    f_nq2[i]->SetLineStyle(7);
    f_afb[i]->SetTitle( Form( "A_{FB}, (%d q^{2}-bin)", i+1) );
    f_afb[i]->SetParNames  ( "slope", "offset" );
    f_afb[i]->SetParameters(   1,        0 );
    f_afb[i]->SetLineColor(2);
    f_afb[i]->SetLineStyle(7);
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
    Double_t incl_ratio, x, xE;
    Int_t    q2bin;
    sTmp >> incl_ratio >> x >> xE;
    if( fl_message ) std::cout << std::setw(5) << std::right << Form("cnt = %d, ", cnt_ntot)
			       << std::setw(8) << std::right << Form("incl_ratio = %f, ", incl_ratio)
			       << std::setw(8) << std::right << Form("x =%f ", x )
			       << std::setw(8) << std::right << Form("+- %f",  xE)
			       << std::endl;
    g_ntot->SetPoint     ( cnt_ntot, incl_ratio, x  );
    g_ntot->SetPointError( cnt_ntot,          0, xE );
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
    Double_t incl_ratio, x, xE;
    Int_t    q2bin;
    sTmp >> incl_ratio >> q2bin >> x >> xE;
    if( fl_message ) std::cout << std::setw(5) << std::right << Form("cnt = %d, ", cnt_nq2[q2bin])
			       << std::setw(8) << std::right << Form("incl_ratio = %f, ", incl_ratio)
			       << std::setw(8) << std::right << Form("q2 bin = %d, ", q2bin)
			       << std::setw(8) << std::right << Form("x =%f ", x )
			       << std::setw(8) << std::right << Form("+- %f",  xE)
			       << std::endl;
    g_nq2[q2bin]->SetPoint     ( cnt_nq2[q2bin], incl_ratio, x  );
    g_nq2[q2bin]->SetPointError( cnt_nq2[q2bin],          0, xE );
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
    Double_t incl_ratio, x, xE;
    Int_t    q2bin;
    sTmp >> incl_ratio >> q2bin >> x >> xE;
    if( fl_message ) std::cout << std::setw(5) << std::right << Form("cnt = %d, ", cnt_afb[q2bin])
			       << std::setw(8) << std::right << Form("incl_ratio = %f, ", incl_ratio)
			       << std::setw(8) << std::right << Form("q2 bin = %d, ", q2bin)
			       << std::setw(8) << std::right << Form("x =%f ", x )
			       << std::setw(8) << std::right << Form("+- %f",  xE)
			       << std::endl;
    g_afb[q2bin]->SetPoint     ( cnt_afb[q2bin], incl_ratio, x  );
    g_afb[q2bin]->SetPointError( cnt_afb[q2bin],          0, xE );
    cnt_afb[q2bin]++;
  }
  fin_afb.close();
  
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  TCanvas* c1 = Canvas( "c1","c1",4, 3 );
  for( Int_t i=0; i<nbin_q2; i++ ){
    c1->cd(i+1);
    g_nq2[i]->Sort();
    g_nq2[i]->SetTitle( Form("N_{sig} (%d q^{2} bin)",i+1) );
    g_nq2[i]->GetXaxis()->SetTitle("A_{FB}^{true} (input)" );
    g_nq2[i]->GetYaxis()->SetTitle("N_{sig}");
    g_nq2[i]->GetXaxis()->CenterTitle();
    g_nq2[i]->GetYaxis()->CenterTitle();
    g_nq2[i]->Draw("AP");
    g_nq2[i]->Fit( Form("f_nq2%d",i) );

    c1->cd(i+1+nbin_q2);
    g_afb[i]->Sort();
    g_afb[i]->SetTitle( Form("A_{FB} (%d q^{2} bin)", i+1) );
    g_afb[i]->GetXaxis()->SetTitle("A_{FB}^{true} (input)"  );
    g_afb[i]->GetYaxis()->SetTitle("A_{FB}^{true} (output)" );
    g_afb[i]->GetXaxis()->CenterTitle();
    g_afb[i]->GetYaxis()->CenterTitle();
    g_afb[i]->Draw("AP");
    g_afb[i]->Fit( Form("f_afb%d",i) );
  }

  c1->cd(2*nbin_q2+1);
  g_ntot->Sort();
  g_ntot->SetTitle( "N_{sig}(total)" );
  g_ntot->GetXaxis()->SetTitle("A_{FB}^{true} (input)");
  g_ntot->GetYaxis()->SetTitle("N_{sig}(total)");
  g_ntot->GetXaxis()->CenterTitle();
  g_ntot->GetYaxis()->CenterTitle();
  g_ntot->Draw("AP");
  g_ntot->Fit("f_ntot");

  // SAVE
  if( fl_save ) c1->Print( Form("pic/linearity_check_bin3_%d_sim.eps",tmp_i) );
  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();

  return 0;
}

