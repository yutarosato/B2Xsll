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
  const Double_t incl_ratio_min =  0.0;
  const Double_t incl_ratio_max = 1000.0;
  //++++++++++++++++++++++++++++++
  TGraphErrors*  g_ntot = new TGraphErrors();
  TGraphErrors** g_nq2  = new TGraphErrors*[nbin_q2];
  TGraphErrors** g_afb  = new TGraphErrors*[nbin_q2];
  TPaveText**    t_afb  = new TPaveText*   [nbin_q2];

  TF1*  f_ntot = new TF1( "f_ntot", "[0]*x+[1]", incl_ratio_min, incl_ratio_max );
  TF1** f_nq2  = new TF1*[nbin_q2];
  TF1** f_afb  = new TF1*[nbin_q2];

  f_ntot->SetTitle( "N_{sig}(total)" );
  f_ntot->SetParNames  ( "slope", "offset" );
  f_ntot->SetParameters(    300,       10 );
  f_ntot->SetLineColor(2);
  f_ntot->SetLineStyle(7);

  for( Int_t i=0; i<nbin_q2; i++ ){
    g_nq2[i] = new TGraphErrors();
    g_afb[i] = new TGraphErrors();

    f_nq2[i] = new TF1( Form("f_nq2%d",i), "[0]*x+[1]", incl_ratio_min, incl_ratio_max );
    f_afb[i] = new TF1( Form("f_afb%d",i), "[0]",       incl_ratio_min, incl_ratio_max );
    f_nq2[i]->SetTitle( Form( "N_{q2}, (%d q^{2}-bin)", i+1) );
    f_nq2[i]->SetParNames  ( "slope", "offset" );
    f_nq2[i]->SetParameters(    1 ,       0 );
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
  for( Int_t i=0; i<nbin_q2; i++ ){
    l_afb_gen[i] = new TLine(incl_ratio_min, AFB_gen_6qbin[2][i], incl_ratio_max, AFB_gen_6qbin[2][i]);
    l_afb_gen[i]->SetLineColor(3);
    l_afb_gen[i]->SetLineStyle(7);
    l_afb_eff[i] = new TLine(incl_ratio_min, AFB_eff_6qbin[2][i], incl_ratio_max, AFB_eff_6qbin[2][i]);
    l_afb_eff[i]->SetLineColor(4);
    l_afb_eff[i]->SetLineStyle(7);
  }
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Char_t buf[1024];

  Double_t pre_incl_ratio= 0;
  Int_t    pre_q2bin     = 0;
  Int_t    tmp_cnt       = 0;
  Int_t    tmp_cnt_line  = 0;
  Double_t calx          = 0;
  Double_t calxE         = 0;
  //======================================================================================  
  // [Ntot]
  ifstream fin_ntot;
  fin_ntot.open(file_ntot);
  Int_t    cnt_ntot  = 0;
  Double_t var_ntot  = 0;
  Double_t var_ntotE = 0;
  if( fl_message ) std::cout << "[READ DATA(Ntot)]" << std::endl;
  pre_incl_ratio = 0;
  pre_q2bin      = 0;
  tmp_cnt        = 0;
  tmp_cnt_line   = 0;
  calx           = 0;
  calxE          = 0;
  while(!fin_ntot.eof()){
    fin_ntot.getline(buf,1024);
    if(buf[0]=='#') continue; // comment
    if(buf[0]=='*'){
      std::cout << "[DUMP] " << pre_incl_ratio << " : " << calx/tmp_cnt << " +- " << sqrt(calxE)/tmp_cnt << std::endl;
      g_ntot->SetPoint     ( cnt_ntot, pre_incl_ratio,        calx/tmp_cnt );
      g_ntot->SetPointError( cnt_ntot,              0, sqrt(calxE)/tmp_cnt );
      cnt_ntot++;
      tmp_cnt = 0;
      pre_incl_ratio = 0;
      break;    // finish
    }
    std::istringstream sTmp(buf);
    Double_t incl_ratio, x, xE;
    Int_t    q2bin;
    sTmp >> incl_ratio >> x >> xE;
    if( fl_message ) std::cout << std::setw(5) << std::right << Form("cnt = %d, ", cnt_ntot)
			       << std::setw(8) << std::right << Form("incl_ratio = %f, ", incl_ratio)
			       << std::setw(8) << std::right << Form("x =%f ", x )
			       << std::setw(8) << std::right << Form("+- %f",  xE)
			       << std::endl;
    if     ( tmp_i==0 ) incl_ratio *= 139.90; // tmpppp
    else if( tmp_i==1 ) incl_ratio *= 163.71; // tmpppp
    else if( tmp_i==2 ) incl_ratio *= 303.61; // tmpppp
    if( tmp_cnt_line==0 ) pre_incl_ratio = incl_ratio;
    tmp_cnt_line++;
    if( incl_ratio != pre_incl_ratio ){
      std::cout << "[DUMP] " << pre_incl_ratio << " : " << calx/tmp_cnt << " +- " << sqrt(calxE)/tmp_cnt << std::endl;
      g_ntot->SetPoint     ( cnt_ntot, pre_incl_ratio,       calx/tmp_cnt  );
      g_ntot->SetPointError( cnt_ntot,              0, sqrt(calxE)/tmp_cnt );
      cnt_ntot++;
      calx    = x;
      calxE   = xE*xE;
      tmp_cnt = 1;
      pre_incl_ratio = incl_ratio;
    }else{
      calx  += x;
      calxE += xE*xE;
      tmp_cnt++;
    }
  }
  fin_ntot.close();

  //======================================================================================  
  // [Nq2]
  ifstream fin_nq2;
  fin_nq2.open(file_nq2);
  Int_t    cnt_nq2 [nbin_q2] = {0};
  Double_t var_nq2 [nbin_q2] = {0};
  Double_t var_nq2E[nbin_q2] = {0};
  if( fl_message ) std::cout << "[READ DATA(Nq2)]" << std::endl;
  pre_incl_ratio = 0;
  pre_q2bin      = 0;
  tmp_cnt        = 0;
  tmp_cnt_line   = 0;
  calx           = 0;
  calxE          = 0;
  while(!fin_nq2.eof()){
    fin_nq2.getline(buf,1024);
    if(buf[0]=='#') continue; // comment
    if(buf[0]=='*'){
      std::cout << "[DUMP] " << pre_incl_ratio << " : " << calx/tmp_cnt << " +- " << sqrt(calxE)/tmp_cnt << std::endl;
      g_nq2[pre_q2bin]->SetPoint     ( cnt_nq2[pre_q2bin], pre_incl_ratio,        calx/tmp_cnt );
      g_nq2[pre_q2bin]->SetPointError( cnt_nq2[pre_q2bin],              0, sqrt(calxE)/tmp_cnt );
      cnt_nq2[pre_q2bin]++;
      tmp_cnt = 0;
      pre_incl_ratio = 0;
      pre_q2bin      = 0;
      break;    // finish
    }
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
    if     ( tmp_i==0 && nbin_q2==4 && q2bin==0 ) incl_ratio *= 27.65+30.25; // tmpppp
    else if( tmp_i==0 && nbin_q2==4 && q2bin==1 ) incl_ratio *= 18.46+15.44; // tmpppp
    else if( tmp_i==0 && nbin_q2==4 && q2bin==2 ) incl_ratio *=  8.58+ 4.90; // tmpppp
    else if( tmp_i==0 && nbin_q2==4 && q2bin==3 ) incl_ratio *= 21.25+13.37; // tmpppp
    else if( tmp_i==1 && nbin_q2==4 && q2bin==0 ) incl_ratio *= 21.36+22.29; // tmpppp
    else if( tmp_i==1 && nbin_q2==4 && q2bin==1 ) incl_ratio *= 23.12+19.37; // tmpppp
    else if( tmp_i==1 && nbin_q2==4 && q2bin==2 ) incl_ratio *= 20.48+11.83; // tmpppp
    else if( tmp_i==1 && nbin_q2==4 && q2bin==3 ) incl_ratio *= 27.89+17.37; // tmpppp
    else if( tmp_i==2 && nbin_q2==4 && q2bin==0 ) incl_ratio *= 49.01+52.54; // tmpppp
    else if( tmp_i==2 && nbin_q2==4 && q2bin==1 ) incl_ratio *= 41.58+34.81; // tmpppp
    else if( tmp_i==2 && nbin_q2==4 && q2bin==2 ) incl_ratio *= 29.06+16.73; // tmpppp
    else if( tmp_i==2 && nbin_q2==4 && q2bin==3 ) incl_ratio *= 49.14+30.74; // tmpppp
    if( tmp_cnt_line==0 ){
      pre_incl_ratio = incl_ratio;
      pre_q2bin      = q2bin;
    }
    tmp_cnt_line++;
    if( incl_ratio != pre_incl_ratio || q2bin != pre_q2bin ){
      std::cout << "[DUMP] " << pre_incl_ratio << " : " << calx/tmp_cnt << " +- " << sqrt(calxE)/tmp_cnt << std::endl;
      g_nq2[pre_q2bin]->SetPoint     ( cnt_nq2[pre_q2bin], pre_incl_ratio,        calx/tmp_cnt );
      g_nq2[pre_q2bin]->SetPointError( cnt_nq2[pre_q2bin],              0, sqrt(calxE)/tmp_cnt );
      cnt_nq2[pre_q2bin]++;
      calx    = x;
      calxE   = xE*xE;
      tmp_cnt = 1;
      pre_incl_ratio = incl_ratio;
      pre_q2bin      = q2bin;
    }else{
      calx  += x;
      calxE += xE*xE;
      tmp_cnt++;
    }
  }
  fin_nq2.close();

  //======================================================================================  
  // [AFB]
  ifstream fin_afb;
  fin_afb.open(file_afb);
  Int_t    cnt_afb [nbin_q2] = {0};
  Double_t var_afb [nbin_q2] = {0};
  Double_t var_afbE[nbin_q2] = {0};

  if( fl_message ) std::cout << "[READ DATA(AFB)]" << std::endl;
  pre_incl_ratio = 0;
  pre_q2bin      = 0;
  tmp_cnt        = 0;
  tmp_cnt_line   = 0;
  calx           = 0;
  calxE          = 0;
  while(!fin_afb.eof()){
    fin_afb.getline(buf,1024);
    if(buf[0]=='#') continue; // comment
    if(buf[0]=='*'){
      std::cout << "[DUMP] " << pre_incl_ratio << " : " << calx/tmp_cnt << " +- " << sqrt(calxE)/tmp_cnt << std::endl;
      g_afb[pre_q2bin]->SetPoint     ( cnt_nq2[pre_q2bin], pre_incl_ratio,        calx/tmp_cnt );
      g_afb[pre_q2bin]->SetPointError( cnt_nq2[pre_q2bin],              0, sqrt(calxE)/tmp_cnt );
      cnt_afb[pre_q2bin]++;
      tmp_cnt = 0;
      pre_incl_ratio = 0;
      pre_q2bin      = 0;
      break;    // finish
    }
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
    if( tmp_cnt_line==0 ){
      pre_incl_ratio = incl_ratio;
      pre_q2bin      = q2bin;
    }
    tmp_cnt_line++;
    if( incl_ratio != pre_incl_ratio || q2bin != pre_q2bin ){
      std::cout << "[DUMP] " << pre_incl_ratio << " : " << calx/tmp_cnt << " +- " << sqrt(calxE)/tmp_cnt << std::endl;
      g_afb[pre_q2bin]->SetPoint     ( cnt_afb[pre_q2bin], pre_incl_ratio,        calx/tmp_cnt );
      g_afb[pre_q2bin]->SetPointError( cnt_afb[pre_q2bin],              0, sqrt(calxE)/tmp_cnt );
      cnt_afb[pre_q2bin]++;
      calx    = x;
      calxE   = xE*xE;
      tmp_cnt = 1;
      pre_incl_ratio = incl_ratio;
      pre_q2bin      = q2bin;
    }else{
      calx  += x;
      calxE += xE*xE;
      tmp_cnt++;
    }
  }
  fin_afb.close();
  
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  TCanvas* c1 = Canvas( "c1","c1",4, 3 );
  for( Int_t i=0; i<nbin_q2; i++ ){
    c1->cd(i+1);
    g_nq2[i]->Sort();
    g_nq2[i]->SetTitle( Form("N_{sig} (%d q^{2} bin)",i+1) );
    g_nq2[i]->GetXaxis()->SetTitle("N_{sig} (input)" );
    g_nq2[i]->GetYaxis()->SetTitle("N_{sig} (output)");
    g_nq2[i]->GetXaxis()->CenterTitle();
    g_nq2[i]->GetYaxis()->CenterTitle();
    g_nq2[i]->Draw("AP");
    g_nq2[i]->Fit( Form("f_nq2%d",i) );

    c1->cd(i+1+nbin_q2);
    g_afb[i]->Sort();
    g_afb[i]->SetTitle( Form("A_{FB} (%d q^{2} bin)", i+1) );
    g_afb[i]->GetXaxis()->SetTitle("incl ratio" );
    g_afb[i]->GetYaxis()->SetTitle("A_{FB}" );
    g_afb[i]->GetXaxis()->CenterTitle();
    g_afb[i]->GetYaxis()->CenterTitle();
    g_afb[i]->Draw("AP");
    g_afb[i]->Fit( Form("f_afb%d",i) );
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
  g_ntot->GetXaxis()->SetTitle("N_{sig}(total,input)" );
  g_ntot->GetYaxis()->SetTitle("N_{sig}(total,output)");
  g_ntot->GetXaxis()->CenterTitle();
  g_ntot->GetYaxis()->CenterTitle();
  g_ntot->Draw("AP");
  g_ntot->Fit("f_ntot");

  // Tlegend
  c1->cd(2*nbin_q2+2);
  TLegend* legend1 = new TLegend( 0.75,0.75,0.99,0.99 );
  legend1->AddEntry( l_afb_gen[0], "Gen",           "L" );
  legend1->AddEntry( l_afb_eff[0], "Eff",           "L" );
  legend1->AddEntry( f_afb[0],     "Rec(Fit->Gen)", "L" );
  legend1->Draw();

  // SAVE
  if( fl_save ) c1->Print( Form("pic/linearity_check_bin2_toy_%d_sim.eps",tmp_i) );
  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();

  return 0;
}

