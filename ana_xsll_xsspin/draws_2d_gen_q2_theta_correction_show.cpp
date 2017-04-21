#include <iostream>
#include <math.h>

#include "../Util/Style.h"
#include "../Util/Canvas.h"
#include "../Util/const.h"
#include "../Util/MChain.h"
#include "../Util/MCut.h"
#include "../Util/MCut_array.h"
#include "../Util/FitFunc.h"
#include "../Util/Manip.h"
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
#include <TH2D.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TString.h>
#include <TArrow.h>

using namespace q2_theta_uniform_gen;
//using namespace q2_theta_nonuniform_gen;
const Bool_t  flag_save = true;

Int_t main( Int_t argc, Char_t** argv ){
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==4 || argc==5) ) std::cerr << "wrong input" << std::endl
					<< " Usage : ./draws_2d (int)fl_mode_ll (int)fl_xsid (char*)setname [(int)fl_appRun]" << std::endl
					<< "[fl_mode_ll] : 1(e), 0(mu)"                << std::endl
					<< "[ fl_xsid  ] : 0(all), 1(K), 2(K*), 3(Xs)" << std::endl, abort();
  
  Int_t   fl_mode_ll = atoi( argv[1] ); // 1(e), 0(mu)
  Int_t   fl_xsid    = atoi( argv[2] ); // 0(all), 1(K), 2(K*), 3(Xs)
  Char_t* setname    = argv[3]; 
  Int_t   fl_appRun  = 1;
  if( argc==5 ) fl_appRun = atoi( argv[4] );
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t  Nchain        = 3;
  const Int_t  nfile[Nchain] = {0};
  
  Char_t** infile = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ) infile[i] = new Char_t[1024]; 
  
  std::stringstream sTmp;
  for( Int_t i=0; i<Nchain; i++ ){
    if     ( i==0 ) sTmp << "~/ewp/ana/data/sigmc/hbk6/right/hbk_org/sigMC_*_set["           << setname << "]";
    else if( i==1 ) sTmp << "~/ewp/ana/data/sigmc_xsspin/hbk6/right/hbk_org/sigMC_*_set["    << setname << "]";
    else if( i==2 ) sTmp << "~/ewp/ana/data/sigmc_xsspin/hbk6/right/hbk_org/sigMC_*_set["    << setname << "]";
    //else if( i==1 ) sTmp << "~/ewp/modules/xsll_sigmc_xsspin/data_noweight/hbk/sigMC_*_set[" << setname << "]";
    //else if( i==2 ) sTmp << "~/ewp/modules/xsll_sigmc_xsspin/data_noweight/hbk/sigMC_*_set[" << setname << "]";
    strcpy( infile[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    add_cut[i] = new Char_t[4096];
    if     ( fl_xsid==0 ) sTmp << "1"; // h12 (all)
    else if( fl_xsid==1 ) sTmp << "1&&gm_fl_xs<0&&Xs_m<0.50"; // h12 (K),  gm_fl_xs : 1(xs), -1(k,k*), gm_xs, gm_bg, Xs_m
    else if( fl_xsid==2 ) sTmp << "1&&gm_fl_xs<0&&Xs_m>0.50"; // h12 (K*),
    else if( fl_xsid==3 ) sTmp << "1&&gm_fl_xs>0&&Xs_m>1.10"; // h12 (Xs),  

    if     ( i==1 ) sTmp << " && (abs(Xs_id)==30343 || abs(Xs_id)==30353) "; // xs-spin0
    else if( i==2 ) sTmp << " && (abs(Xs_id)==30344 || abs(Xs_id)==30354) "; // xs-spin1

    //sTmp << " && llg_m*llg_m<4.3"; // tmpppppp for checking cos-theta distribution
    strcpy( add_cut[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain      = new MChain*[Nchain];
  TH2D**    tmphist    = new TH2D*  [Nchain];
  TH2D**    hist       = new TH2D*  [Nchain+1];
  TH1D**    histx_gen  = new TH1D*  [Nchain];
  TH1D**    histy_gen  = new TH1D*  [Nchain];
  TH1D**    hist_afb   = new TH1D*  [Nchain+1];
  TCanvas*  c1         = Canvas( "c1","c1",3, (int)ceil(Nchain/3.0)+2 );
  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++

  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form(" ************************ make tmphist ( %s, %s ) *************************************",tname,axis) << std::endl;
    std::cout << Form( "<infile %d > ", j );
    chain[j] = new MChain( infile[j], tname, branch_table(), nfile[j], "*.root" );
    nominal_cut_selection( chain[j], fl_mode_ll )( chain[j]->GetCut(), tname );
  }
  
  // ++++++++++++++++++++++++
  // cut change
  for( Int_t j=0; j<Nchain; j++ ) chain[j]->GetCut()->Set( "gm_fl_xs*Xs_m", 1, 0.0, 0.0, 1.10 );
  //chain[0]->GetCut()->Set(    441, 0 );
  //chain[0]->GetCut()->Set(    443, 0 );
  //chain[0]->GetCut()->Set( 100443, 0 );
  //
  // ++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ changed cut *************************************" << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<infile %d > ", j ) << std::endl;
    chain[j]->GetCut()->Display(0);
    chain[j]->MakeTree();
    std::cout << "add_cut : " << add_cut[j] << std::endl;
  }
  
  // +++++++ make tmphist ++++++++++++++++++++++++++++++++++
  for( Int_t j=0; j<Nchain; j++ ){
    tmphist[j] = new TH2D( Form("tmphist%d",j), Form("%s",chain[j]->GetChange()), xbin,xmin,xmax, ybin,ymin,ymax );
    chain[j]->GetTree()->Project( Form("tmphist%d",j), axis, add_cut[j] );
  }
  
  // +++++++ make hist  ++++++++++++++++++++++++++++++++++
  for( Int_t i=0; i<Nchain+1; i++ ){
    hist[i] = new TH2D( Form("hist%d",i),Form("hist%d",i), xbin,xmin,xmax,ybin,ymin,ymax );
    hist[i]->GetXaxis()->CenterTitle();
    hist[i]->GetYaxis()->CenterTitle();
    hist[i]->SetXTitle( xlabel );
    hist[i]->SetYTitle( ylabel );
    if( i<Nchain ) hist[i]->Add( tmphist[i] );
    else{
      hist[i]->Add( hist[1] );
      hist[i]->Add( hist[2] );
    }
    // +++++++ display ++++++++++++++++++++++++++++++++++
    std::cout << std::setw(10) << Form("<hist%d>",i);
    std::cout << std::setw(9)  << hist[i]->GetEntries() << " events"
	      << std::endl; 
  }

  hist[0]->SetTitle("defalut X_{s}" );
  hist[1]->SetTitle("X_{s} (spin-0)");
  hist[2]->SetTitle("X_{s} (spin-1)");
  hist[3]->SetTitle("X_{s} (spin)"  );
  
  // +++++++ projection hist ++++++++++++++++++
  for( Int_t i=0; i<Nchain; i++ ){
    histx_gen[i] = new TH1D(*(hist[i]->ProjectionX( Form("_px%d",i), 0, -1, "e")) ); Deco( histx_gen[i], 0, i+1, i+1 );
    histy_gen[i] = new TH1D(*(hist[i]->ProjectionY( Form("_py%d",i), 0, -1, "e")) ); Deco( histy_gen[i], 0, i+1, i+1 );
  }

  // +++++++ make hist(afb) ++++++++++++++++++
  for( Int_t i=0; i<Nchain+1; i++ ){
    hist_afb[i] = new TH1D( Form("afb_rec%d",i), Form("afb_rec%d",i), xbin, xmin, xmax);
    hist_afb[i]->SetXTitle( xlabel   );
    hist_afb[i]->SetYTitle( "A_{FB}" );
    Deco( hist_afb[i], 0, i+1, i+1 );
  }
  
  for( Int_t m=0; m<xbin; m++ ){
    for( Int_t i=0; i<Nchain; i++ ){
      Double_t theta_p = hist[i]->Integral( m+1,m+1,ybin/2+1, ybin   );
      Double_t theta_m = hist[i]->Integral( m+1,m+1,       1, ybin/2 );
      if( theta_p+theta_m==0 ) continue;
      Double_t afb     = (theta_p - theta_m ) / (theta_p + theta_m );
      Double_t afbE    = 2/(theta_p+theta_m)/(theta_p+theta_m)*sqrt(theta_p*theta_m*(theta_p+theta_m));
      hist_afb[i]->SetBinContent( m+1, afb  );
      hist_afb[i]->SetBinError  ( m+1, afbE );
    }

    // for (spin0)+(spin1)
    Double_t tmp_theta_p = hist[1]->Integral( m+1,m+1,ybin/2+1, ybin   ) + hist[2]->Integral( m+1,m+1,ybin/2+1, ybin   );
    Double_t tmp_theta_m = hist[1]->Integral( m+1,m+1,       1, ybin/2 ) + hist[2]->Integral( m+1,m+1,       1, ybin/2 );
    if( tmp_theta_p+tmp_theta_m==0 ) continue;
    Double_t tmp_afb     = (tmp_theta_p - tmp_theta_m ) / (tmp_theta_p + tmp_theta_m );
    Double_t tmp_afbE    = 2/(tmp_theta_p+tmp_theta_m)/(tmp_theta_p+tmp_theta_m)*sqrt(tmp_theta_p*tmp_theta_m*(tmp_theta_p+tmp_theta_m));
    hist_afb[Nchain]->SetBinContent( m+1, tmp_afb  );
    hist_afb[Nchain]->SetBinError  ( m+1, tmp_afbE );
  }

  // +++++++ draw ++++++++++++++++++++++++++++++++++
  TLegend* legend = new TLegend( 0.75,0.75,0.99,0.99 );
  legend->SetHeader( Form("lep%d, xsid%d", fl_mode_ll, fl_xsid) );
  for( Int_t i=0; i<Nchain+1; i++ ) legend->AddEntry( hist_afb[i], hist[i]->GetTitle(), "P" );
  
  c1->Draw();
  c1->cd(1);
  TH2D* waku = new TH2D("AFB","A_{FB}",2,0,25,2,-0.4,0.4);
  waku->SetXTitle( xlabel   );
  waku->SetYTitle( "A_{FB}" );
  waku->GetXaxis()->CenterTitle();
  waku->GetYaxis()->CenterTitle();
  waku->Draw();
  for( Int_t i=0; i<Nchain+1; i++ ){
    hist_afb[i]->Draw("PE0same");
  }

  for( Int_t i=0; i<Nchain; i++ ) histx_gen[i]->Scale(1/histx_gen[i]->Integral());
  for( Int_t i=0; i<Nchain; i++ ) histy_gen[i]->Scale(1/histy_gen[i]->Integral());

  c1->cd(2);
  TH2D* wakugenx = Waku( Nchain, histx_gen, xlabel, Form("%s (gen)", xlabel), Form("%s (gen)", xlabel) );
  wakugenx->Draw();
  legend->Draw();
  for( Int_t i=0; i<Nchain; i++ ) histx_gen[i]->Draw("Psame");

  c1->cd(3);
  TH2D* wakugeny = Waku( Nchain, histy_gen, ylabel, Form("%s (gen)", ylabel), Form("%s (gen)", ylabel) );
  wakugeny->Draw();
  for( Int_t i=0; i<Nchain; i++ ) histy_gen[i]->Draw("Psame");

  for( Int_t i=0; i<Nchain; i++ ){
    c1->cd(i+7);
    hist[i]->Draw("COLZ");
  }
  
  
  // +++++++ make hist (Xs, ll) ++++++++++++++++++++++++++++++++++
  TH1D**    hist_xs    = new TH1D*[Nchain];
  TH1D**    hist_ll    = new TH1D*[Nchain];
  TH1D**    hist_q2    = new TH1D*[Nchain];
  TF1**     func       = new TF1* [Nchain];
  const Int_t sel_fun  = 3; // cubic

  for( Int_t j=0; j<Nchain; j++ ){
    hist_xs[j] = new TH1D( Form("hist_xs_%d",j), Form("%s",chain[j]->GetChange()), 100, 0,  5.0 );
    hist_ll[j] = new TH1D( Form("hist_ll_%d",j), Form("%s",chain[j]->GetChange()), 100, 0,  5.0 );
    hist_q2[j] = new TH1D( Form("hist_q2_%d",j), Form("%s",chain[j]->GetChange()), 100, 0, 25.0 );
    chain[j]->GetTree()->Project( Form("hist_xs_%d",j), "Xs_m",        add_cut[j] );
    chain[j]->GetTree()->Project( Form("hist_ll_%d",j), "llg_m",       add_cut[j] );
    chain[j]->GetTree()->Project( Form("hist_q2_%d",j), "llg_m*llg_m", add_cut[j] );
    hist_xs[j]->Sumw2();
    hist_ll[j]->Sumw2();
    hist_q2[j]->Sumw2();
    if( j!=0 ) hist_xs[j]->Scale( 2.0 );
    Deco( hist_xs[j], 0, j+1, j+1 );
    Deco( hist_ll[j], 0, j+1, j+1 );
    Deco( hist_q2[j], 0, j+1, j+1 );
    //c1->cd(4); hist_xs[j]->DrawCopy( j==Nchain-1 ? "P" : "Psame" );
    c1->cd(4); hist_xs[j]->DrawCopy( j==0 ? "P" : "Psame" );
  }
  c1->cd(5);
  /*
  hist_ll[2]->DrawCopy( "P"     );
  hist_ll[0]->DrawCopy( "Psame" );
  hist_ll[1]->DrawCopy( "Psame" );
  */
  ///*
  hist_ll[0]->DrawCopy( "P"     );
  hist_ll[1]->DrawCopy( "Psame" );
  hist_ll[2]->DrawCopy( "Psame" );
  hist_ll[2]->Add( hist_ll[1] );
  Deco( hist_ll[2], 0, Nchain+1, Nchain+1 );
  hist_ll[2]->SetLineStyle(7);
  hist_ll[2]->DrawCopy( "Psame" );
  //*/
  for( Int_t j=1; j<Nchain; j++ ){
    c1->cd(6); hist_xs[j]->Divide(hist_xs[0],hist_xs[j]); hist_xs[j]->DrawCopy( j==1 ? "P" : "Psame" );
    func[j] = new TF1( Form("func%d",j), "[0]*exp(-[1]*x*x*x*x+[2]*x*x*x+[3]*x*x+[4]*x+[5])", 1.1, 5);
    func[j]->SetParameters(0.1, 0.1, 0.1, 1.0);
    Deco( func[j], 0, j+1, j+1 );
    hist_xs[j]->Fit( func[j],"R", j==1 ? "PE0" : "PE0same" );
  }
  
  c1->Update();
  if( flag_save ) c1->Print( Form("pic/2d_gen_q2_theta_correction_show_lep%d_xsid%d_set%s.eps", fl_mode_ll, fl_xsid, setname) );
  std::cout << "finish" << std::flush;  
  
  if( fl_appRun ){
    //gROOT->Idle( 10, ".q");
    app.Run();
    //std::cout << "continueeeee" << std::endl;
    //c1->clear();
  }


  delete[] chain;
  delete[] tmphist;
  delete[] hist;
  delete   c1;
    
  return 0;
}
