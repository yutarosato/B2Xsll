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
#include <TFile.h>
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

Int_t main( Int_t argc, Char_t** argv ){
  using namespace sig_bkg;
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==4 || argc==5) ) std::cerr << "wrong input" << std::endl
					<< " Usage : ./draws_1d_diagonal (char*)stream (char*)setname (int)fl_mode_ll [(int)fl_appRun]" << std::endl
					<< "[  stream  ] 0, 1, ..., 5, 01, 0-2"     	  << std::endl
    					<< "[  setname ] A, B, ..., U, AB, A-U"     	  << std::endl, abort();
  Char_t* stream     = argv[1];
  Char_t* setname    = argv[2];
  Int_t   fl_mode_ll = atoi( argv[3] ); // 1-> e, 0-> mu
  Int_t   fl_appRun  = 1;
  if( argc==5 ) fl_appRun = atoi( argv[4] );

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t Nchain        = 6;     // [sig(w/o pi0), sig(w/  pi0), uds, charm, mixed, charged]          // fixed !!
  const Int_t Nhist         = 5;     // [sig(w/o pi0), sig(w/  pi0), bkg(qq), bkg(bb), bkg(total)] // fixed !!
  const Int_t fl_sb[Nchain] = {1,2}; // 1(sig[w/o pi0]), 2(sig[w/ pi0]), 0(bkg)
  const Int_t nfile[Nchain] = {0};
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  Char_t** infile = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    infile[i] = new Char_t[1024];
    if( fl_sb[i] ) sTmp << indir[1] << "/sigMC_*_m9999m_caseB_set[" << setname << "]";         // sig
    else           sTmp << indir[0] << "/gMC_" << bkgtype[i-2] << "_*_s0[" << stream << "]";  // bkg
    strcpy( infile[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t add[Nhist][Nchain] ={
    {1},           // sig(w/o pi0)
    {0,1},         // sig(w/  pi0)
    {0,0,1,1},     // bkg(qq)
    {0,0,0,0,1,1}, // bkg(bb)
    {0,0,1,1,1,1}, // bkg(total)
  };

  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++) add_cut[i] = new Char_t[4096];
  strcpy( add_cut[0], makeCut_pi0(0,1).c_str() ); // sig(w/o pi0)
  strcpy( add_cut[1], makeCut_pi0(1,1).c_str() ); // sig(w/  pi0)
  add_cut[2] = "";
  add_cut[3] = "";
  add_cut[4] = "lpgt!=3 || lmgt!=3";
  add_cut[5] = "lpgt!=3 || lmgt!=3";

  const Char_t*  tname    = "h511";
  const Char_t*  axis     = "de";
  const Double_t offset   = 0.0;
  const Int_t    xbin     = 120; 
  const Double_t xmin     = -0.12;
  const Double_t xmax     =  0.12;
  const Double_t xmin_fit = -0.05-(Double_t)fl_mode_ll*0.05;
  const Double_t xmax_fit =  0.05;
  const Char_t*  xlabel   = "#DeltaE [GeV]";

  const Bool_t  flag_fit   = true;
  const Int_t   sel_fun[3] = {2, 30000,410}; // bkg,sig(w/o pi0),sig(w/ pi0)
  
  const Bool_t flag_save  = true; // outfile.eps and outfile.root
  Bool_t       flag_scale = true;

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain   = new MChain*[Nchain];
  TH1D**    tmphist = new TH1D*  [Nchain];
  TH1D**    hist    = new TH1D*  [Nhist];
  TF1**     func    = new TF1*   [Nhist];
  TCanvas*  c1      = Canvas( "c1","c1",2 );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++
  std::cout << Form(" ************************ make tmphist ( %s, %s ) *************************************",tname,axis) << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<infile %d > ", j );
    chain[j] = new MChain( infile[j], tname, branch_table(), nfile[j], tail );
    nominal_cut_selection( chain[j], fl_mode_ll )( chain[j]->GetCut(), tname );
  }

  // ++++++++++++++++++++++++
  // cut change
  for( Int_t j=0; j<Nchain; j++ ){
    //chain[j]->GetCut()->Set( "Mbc", 1,  5.22, 0.0, 5.30 );
    chain[j]->GetCut()->Set( "de",  1, -0.20, 0.0, 0.20 );
  }
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
    tmphist[j] = new TH1D( Form("tmphist%d",j), Form("%s",chain[j]->GetChange()), xbin,offset+xmin,offset+xmax );
    chain[j]->GetTree()->Project( Form("tmphist%d",j), axis, add_cut[j] );
  }
  
  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make hist *************************************" << std::endl;
  hist[0] = new TH1D( Form("pdf_hist_de_lep%d_sig",     fl_mode_ll), Form("pdf_hist_de_lep%d_set%s_sig",    fl_mode_ll,setname), xbin,offset+xmin,offset+xmax );
  hist[1] = new TH1D( Form("pdf_hist_de_lep%d_pi0_sig", fl_mode_ll), Form("pdf_hist_de_lep%d_pi0_set%s_sig",fl_mode_ll,setname), xbin,offset+xmin,offset+xmax );
  hist[2] = new TH1D( Form("pdf_hist_de_lep%d_qq_bkg",  fl_mode_ll), Form("pdf_hist_de_lep%d_qq_s0%s_bkg",  fl_mode_ll,stream),  xbin,offset+xmin,offset+xmax );
  hist[3] = new TH1D( Form("pdf_hist_de_lep%d_bb_bkg",  fl_mode_ll), Form("pdf_hist_de_lep%d_bb_s0%s_bkg",  fl_mode_ll,stream),  xbin,offset+xmin,offset+xmax );
  hist[4] = new TH1D( Form("pdf_hist_de_lep%d_tot_bkg", fl_mode_ll), Form("pdf_hist_de_lep%d_tot_s0%s_bkg", fl_mode_ll,stream),  xbin,offset+xmin,offset+xmax );
  for( Int_t i=0; i<Nhist; i++ ){
    //hist[i] = new TH1D( Form("pdf_hist_de_lep%d_s0%s_bkg",fl_mode_ll,stream),Form("pdf_hist_de_lep%d_s0%s_bkg",fl_mode_ll,stream), xbin,offset+xmin,offset+xmax );
    Deco( hist[i], 2, i+1, i+1 );
    for( Int_t j=0; j<Nchain; j++ ){
      if( add[i][j] ) hist[i]->Add( tmphist[j] );
    }
  }
  // +++++++ display ++++++++++++++++++++++++++++++++++
  for( Int_t i=0; i<Nhist; i++ ){
    std::cout << Form("<hist %d> ",i);
    Double_t entry_all    = hist[i]->GetEntries();
    Double_t entry_canvas = hist[i]->Integral();
    Double_t entry_under  = hist[i]->GetBinContent(0);
    Double_t entry_over   = hist[i]->GetBinContent(xbin+1);
    TString sTmp = "Added-Files( ";
    for( Int_t j=0; j<Nchain; j++ ) if( add[i][j] ) sTmp += j, sTmp += ",";
    sTmp += ")[";
    sTmp += entry_all;
    sTmp += " events ( canvas : ";
    sTmp += entry_canvas;
    sTmp += " / under : ";
    sTmp += entry_under;
    sTmp += " / over  : ";
    sTmp += entry_over;
    sTmp += "]";
    
    hist[i]->SetTitle( sTmp.Data() );
    std::cout << sTmp.Data() << std::endl;
    if( flag_scale ){
      hist[i]->Sumw2();
      hist[i]->Scale( 1 / hist[i]->Integral(hist[i]->FindBin(xmin_fit+0.0000001),hist[i]->FindBin(xmax_fit-0.0000001)) );
    }
  }
  
  // +++++++ make fit-func ++++++++++++++++++++++++++++++++++
  if( flag_fit ){
    func[0] = new TF1( Form("pdf_func_de_lep%d_sig",     fl_mode_ll), make_func(sel_fun[fl_sb[0]]),offset+xmin_fit,offset+xmax_fit,n_fitfunc_par(sel_fun[fl_sb[0]]) );
    func[1] = new TF1( Form("pdf_func_de_lep%d_pi0_sig", fl_mode_ll), make_func(sel_fun[fl_sb[1]]),offset+xmin_fit,offset+xmax_fit,n_fitfunc_par(sel_fun[fl_sb[1]]) );
    func[2] = new TF1( Form("pdf_func_de_lep%d_qq_bkg",  fl_mode_ll), make_func(sel_fun[fl_sb[2]]),offset+xmin_fit,offset+xmax_fit,n_fitfunc_par(sel_fun[fl_sb[2]]) );
    func[3] = new TF1( Form("pdf_func_de_lep%d_bb_bkg",  fl_mode_ll), make_func(sel_fun[fl_sb[3]]),offset+xmin_fit,offset+xmax_fit,n_fitfunc_par(sel_fun[fl_sb[3]]) );
    func[4] = new TF1( Form("pdf_func_de_lep%d_tot_bkg", fl_mode_ll), make_func(sel_fun[fl_sb[4]]),offset+xmin_fit,offset+xmax_fit,n_fitfunc_par(sel_fun[fl_sb[4]]) );
    for( Int_t i=0; i<Nhist; i++ ){
      //func[i] = new TF1( Form("pdf_func_de_lep%d_s0%s_bkg",fl_mode_ll,stream), make_func(sel_fun[fl_sb[i]]),offset+xmin_fit,offset+xmax_fit,n_fitfunc_par(sel_fun[fl_sb[i]]) );
      func[i]->SetLineColor(i+1);
      if( sel_fun[fl_sb[i]] == 30000 ){
	func[i]->SetParNames  ("area","alpha_l","n_l","mean","sigma_l","alpha_h","n_h","sigma_h");
	func[i]->SetParLimits(1,  0.00,  5.0);
	func[i]->SetParLimits(2,  0.00,100.0);
	func[i]->SetParLimits(3, -0.01, 0.01);
	func[i]->SetParLimits(4,  0.00, 0.50);
	func[i]->SetParLimits(5,  0.00,  5.0);
	if( i==0 ) func[i]->SetParLimits(6,  0.00,100.0);
	else if( i==1 ) func[i]->SetParLimits(6,  1000000.0, 10000000000.0);
	func[i]->SetParLimits(7,  0.00, 0.50);
	if( fl_mode_ll==1 && i==0 ){ // ee w/o pi0
	  Double_t tmp_sigma = 0.015;
	  Double_t tmp_mu    = hist[i]->GetBinCenter(hist[i]->GetMaximumBin());
	  Double_t tmp_area  = hist[i]->GetMaximum()*sqrt(TMath::TwoPi())*tmp_sigma;
	  func[i]->SetParameters(tmp_area,  1.2,  1.2, tmp_mu, tmp_sigma,  1.4,   8.3,  tmp_sigma);
	}else if( fl_mode_ll==1 && i==1 ){ // ee w/  pi0
	  Double_t tmp_sigma = 0.015;
	  Double_t tmp_mu    = hist[i]->GetBinCenter(hist[i]->GetMaximumBin());
	  Double_t tmp_area  = hist[i]->GetMaximum()*sqrt(TMath::TwoPi())*tmp_sigma;
	  func[i]->SetParameters(tmp_area,  1.0,  1.7, tmp_mu, 1.8*tmp_sigma,  1.2,   10000000,   tmp_sigma);
	  //func[i]->SetParameters(tmp_area,  1.0,  1.7, tmp_mu, tmp_sigma,  1.2,   10000000,   tmp_sigma);
	}else if( fl_mode_ll==0 && i==0 ){ // mu w/o pi0
	  Double_t tmp_sigma = 0.012;
	  Double_t tmp_mu    = hist[i]->GetBinCenter(hist[i]->GetMaximumBin());
	  Double_t tmp_area  = hist[i]->GetMaximum()*sqrt(TMath::TwoPi())*tmp_sigma;
	  func[i]->SetParameters(tmp_area,  1.7,  2.1, tmp_mu, tmp_sigma,  1.6,  10,   tmp_sigma);
	}else if( fl_mode_ll==0 && i==1 ){ // mu w/  pi0
	  Double_t tmp_sigma = 0.015;
	  Double_t tmp_mu    = hist[i]->GetBinCenter(hist[i]->GetMaximumBin());
	  Double_t tmp_area  = hist[i]->GetMaximum()*sqrt(TMath::TwoPi())*tmp_sigma;
	  func[i]->SetParameters(tmp_area,  0.7,  21, tmp_mu, tmp_sigma,  1.2,    1000000,  tmp_sigma);
	}
      }else if( sel_fun[fl_sb[i]] == 410 ){
	func[i]->SetParNames( "area","alpha","n","mean","sigma_l", "sigma_h" );
	func[i]->SetParLimits(1,  0.00,  5.0);
	func[i]->SetParLimits(2,  0.00,100.0);
	func[i]->SetParLimits(3, -0.01, 0.01);
	func[i]->SetParLimits(4,  0.00, 0.50);
	func[i]->SetParLimits(5,  0.00, 0.50);
	if( fl_mode_ll==1 && i==0 ){ // ee w/o pi0
	  Double_t tmp_sigma = 0.015;
	  Double_t tmp_mu    = hist[i]->GetBinCenter(hist[i]->GetMaximumBin());
	  Double_t tmp_area  = hist[i]->GetMaximum()*sqrt(TMath::TwoPi())*tmp_sigma;
	  func[i]->SetParameters(tmp_area,  1.0,  1.3, tmp_mu, tmp_sigma, tmp_sigma);
	  func[i]->SetParLimits(3, -0.01, 0.01);
	}else if( fl_mode_ll==1 && i==1 ){ // ee w/  pi0
	  Double_t tmp_sigma = 0.015;
	  Double_t tmp_mu    = hist[i]->GetBinCenter(hist[i]->GetMaximumBin());
	  Double_t tmp_area  = hist[i]->GetMaximum()*sqrt(TMath::TwoPi())*tmp_sigma;
	  func[i]->SetParameters(tmp_area,  0.8,  2.0, tmp_mu, tmp_sigma, tmp_sigma);
	  func[i]->SetParLimits(3, -0.01, 0.01);
	}else if( fl_mode_ll==0 && i==0 ){ // mu w/o pi0
	  Double_t tmp_sigma = 0.012;
	  Double_t tmp_mu    = hist[i]->GetBinCenter(hist[i]->GetMaximumBin());
	  Double_t tmp_area  = hist[i]->GetMaximum()*sqrt(TMath::TwoPi())*tmp_sigma;
	  func[i]->SetParameters(tmp_area,  1.5,  2.8, tmp_mu, tmp_sigma, tmp_sigma);
	  func[i]->SetParLimits(3, -0.01, 0.01);
	}else if( fl_mode_ll==0 && i==1 ){ // mu w/  pi0
	  Double_t tmp_sigma = 0.015;
	  Double_t tmp_mu    = hist[i]->GetBinCenter(hist[i]->GetMaximumBin());
	  Double_t tmp_area  = hist[i]->GetMaximum()*sqrt(TMath::TwoPi())*tmp_sigma;
	  func[i]->SetParameters(tmp_area,  0.7,  21, tmp_mu, tmp_sigma, tmp_sigma);
	}
      }else{
	func_set_parameters( sel_fun[fl_sb[i]], func[i], hist[i], xbin, offset+xmin, offset+xmax );
      }
      
    }
  }
  
  
  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();
  c1->cd(1);
  TH2D* waku = Waku( Nhist, hist, xlabel );
  ((TGaxis*)waku->GetXaxis())->SetMaxDigits(3);
  ((TGaxis*)waku->GetYaxis())->SetMaxDigits(3);
  waku->SetTitle( Form("%s (fl_mode_ll=%d)", xlabel, fl_mode_ll) );
  waku->Draw();
  
  for(Int_t i=0; i<Nhist; i++ ){
    if( flag_fit ){
      hist[i]->Fit(func[i],"R","PE0same"); // PE0same // RRRRR
      std::cout << "chi2/NDF = " << func[i]->GetChisquare() << " / " << func[i]->GetNDF()
		<< " = " << func[i]->GetChisquare()/func[i]->GetNDF()
		<< std::endl;
    }else hist[i]->Draw( "same" );
  }
  // +++++++ tlegend1 ++++++++++++++++++++++++++++++++++
  TLegend* legend1 = new TLegend( 0.75,0.75,0.99,0.99 );
  legend1->AddEntry( hist[0],"sig(w/o pi0)", "L" );
  legend1->AddEntry( hist[1],"sig(w/  pi0)", "L" );
  legend1->AddEntry( hist[2],"bkg(  qq   )", "L" );
  legend1->AddEntry( hist[3],"bkg(  bb   )", "L" );
  legend1->AddEntry( hist[4],"bkg( tot   )", "L" );
  legend1->Draw();
  //=================================================================
  Double_t de_range[2][2] = {
    { -0.05, 0.05 }, // for mm
    { -0.10, 0.05 }, // for ee
  };
  TArrow* ar1 = new TArrow( de_range[fl_mode_ll][0], 1.02*((TAxis*)waku->GetYaxis())->GetXmax(),
			    de_range[fl_mode_ll][1], 1.02*((TAxis*)waku->GetYaxis())->GetXmax(),
			    0.005,"<|>" );
  ar1->SetLineColor(2);
  ar1->SetFillColor(2);
  ar1->SetLineWidth(2);
  ar1->Draw();
  //=================================================================

  // +++++++ tlegend ++++++++++++++++++++++++++++++++++
  c1->cd(2);
  TPaveText* box    = new TPaveText( 0.0,0.4,1.0,1.0 );
  box->SetTextSize(0.015);
  TLegend*   legend = new TLegend  ( 0.0,0.0,1.0,0.4 );
  legend->SetTextSize(0.015);
  for( Int_t j=0; j<Nchain; j++ ){
    box->AddText( Form("tmphist%d : Mode(%s : %s : %d files)", j,chain[j]->GetModeName(),gSystem->BaseName(chain[j]->GetFileName()),chain[j]->GetNfile()) );
    //box->AddText( Form("   %s", gSystem->DirName(chain[j]->GetFileName())) );
    //box->AddText( Form("   %s", tmphist[j]->GetTitle()) );
    box->AddText( Form("   %s", add_cut[j]) );
  }

  for( Int_t i=0; i<Nhist; i++ ){
    legend->AddEntry( hist[i],"","PL" );
  }
  box->Draw();
  legend->Draw();

  c1->Update();
  if( flag_save ){
    c1->Print(     Form("pic/pdf_%s_lep%d_s0%s_set%s.eps",  axis, fl_mode_ll,stream,setname) );
    TFile outfile( Form("pic/pdf_%s_lep%d_s0%s_set%s.root", axis, fl_mode_ll,stream,setname), "RECREATE" );
    c1->Write();
      for( Int_t i=0; i<Nhist; i++ ) hist[i]->Write();
    if( flag_fit ){
      for( Int_t i=0; i<Nhist; i++ ) func[i]->Write();
    }
    outfile.Close();
  }

  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();

  delete[] chain;
  delete[] tmphist;
  delete[] hist;
  delete[] func;
  delete[] add_cut;
  delete   c1;

  return 0;
}
