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
#include <TGraph.h>
#include <TArrow.h>
#include <TMath.h>
#include <Math/DistFunc.h>

Int_t main( Int_t argc, Char_t** argv ){
  using namespace sigmc;
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==2 || argc==3) ) std::cerr << "wrong input" << std::endl
					<< " Usage : ./draws_1d_diagonal (char*)setname [(int)fl_appRun]"
					<< std::endl, abort();
  Char_t* setname = argv[1];
  Int_t fl_appRun = 1;
  if( argc==3 ) fl_appRun = atoi( argv[2] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t Nchain        = 12; // 0 tracks[e,mu],1 tracks[e,mu],2 tracks[e,mu],3 tracks[e,mu],4 tracks[e,mu],5 tracks[e,mu]
  const Int_t Nhist         =  6; // 0 tracks, 1tracks, 2tracks, 3tracks, [45]tracks, all tracks
  const Int_t nfile[Nchain] = {0};
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Char_t** infile = new Char_t*[Nchain];
  std::stringstream sTmp;
  for( Int_t i=0; i<Nchain; i++ ){
    infile[i] = new Char_t[1024];
    sTmp << indir << "sigMC_*_m9999m_caseB_set[" << setname << "]";
    strcpy( infile[i], (Char_t*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t add[Nhist][Nchain] ={
    {1,1,0,0,0,0,0,0,0,0,0,0}, // 0tracks
    {0,0,1,1,0,0,0,0,0,0,0,0}, // 1tracks
    {0,0,0,0,1,1,0,0,0,0,0,0}, // 2tracks
    {0,0,0,0,0,0,1,1,0,0,0,0}, // 3tracks
    {0,0,0,0,0,0,0,0,1,1,1,1}, // [45]tracks
    {1,1,1,1,1,1,1,1,1,1,1,1}, // all tracks
  };
  Int_t fl_mode_ll[Nchain] = {1,0,1,0,1,0,1,0,1,0,1,0};

  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    add_cut[i] = new Char_t[4096];
    strcpy( add_cut[i], makeCut_track(i/2,1).c_str() );
  }
  
  const Char_t*  tname    = "h511";
  const Char_t*  axis     = "kfbchi/kfbdgf";
  const Double_t offset   = 0.0;
  const Int_t    xbin     = 200;
  const Double_t xmin     = 0;
  const Double_t xmax     = 10.0;
  const Double_t xmin_fit = xmin;
  const Double_t xmax_fit = xmax;
  const Char_t*  xlabel   = "C.L.";

  const Bool_t  flag_fit   = true;
  const Int_t   sel_fun    = 1092;
  
  const Bool_t flag_save  = true; // outfile.eps and outfile.root
  Bool_t       flag_scale = true;

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain   = new MChain*[Nchain];
  TH1D**    tmphist = new TH1D*  [Nchain];
  TH1D**    hist    = new TH1D*  [Nhist];
  TF1**     func    = new TF1*   [Nhist];
  TCanvas*  c1      = Canvas( "c1","c1",2       );
  TCanvas*  c2      = Canvas( "c2","c2",Nhist+1 );
  TCanvas*  c3      = Canvas( "c3","c3",Nhist+1 );

  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++
  std::cout << Form(" ************************ make tmphist ( %s, %s ) *************************************",tname,axis) << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<infile %d > ", j );
    chain[j] = new MChain( infile[j], tname, branch_table(), nfile[j], tail );
    nominal_cut_selection( chain[j], fl_mode_ll[j] )( chain[j]->GetCut(), tname );
  }
  // ++++++++++++++++++++++++
  // cut change
  for( Int_t j=0; j<Nchain; j++ ){
    //chain[j]->GetCut()->Set( "Mbc", 1,  5.22, 0.0, 5.30 );
    //chain[j]->GetCut()->Set( "de",  1, -0.20, 0.0, 0.20 );
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
  for( Int_t i=0; i<Nhist; i++ ){
    if( i==Nhist-1 ) hist[i] = new TH1D(      "pdf_hist_kfbchi2_all_sig",        Form("pdf_hist_kfbchi2_all_set%s_sig",       setname), xbin,offset+xmin,offset+xmax );
    else             hist[i] = new TH1D( Form("pdf_hist_kfbchi2_%dtracks_sig",i),Form("pdf_hist_kfbchi2_%dtracks_set%s_sig",i,setname), xbin,offset+xmin,offset+xmax );
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
      hist[i]->Scale( 1/(hist[i]->GetEntries()) );
    }
  }

  // +++++++ make fit-func ++++++++++++++++++++++++++++++++++
  if( flag_fit ){
    for( Int_t i=0; i<Nhist; i++ ){
      if( i==Nhist-1 ) func[i] = new TF1(      "pdf_func_kfbchi2_all_sig",         make_func(sel_fun),offset+xmin_fit,offset+xmax_fit,n_fitfunc_par(sel_fun) );
      else             func[i] = new TF1( Form("pdf_func_kfbchi2_%dtracks_sig",i), make_func(sel_fun),offset+xmin_fit,offset+xmax_fit,n_fitfunc_par(sel_fun) );
      //func_set_parameters(sel_fun, func[i], hist[i], xbin, offset+xmin, offset+xmax);
      func[i]->SetLineColor(i+2);
      func[i]->SetLineWidth(2);
      if( sel_fun==1092){
	func[i]->SetParNames( "area1","alpha1","n1","area2","alpha2", "n2" );
	if( i==0 ){ // 0-track
	  func[i]->SetParameters(  0.07,   1.6,   1.5, 0.004,   0.5,  1.3  );
	}else if( i==1 ){ // 1-track
	  func[i]->SetParameters(  0.09,   2.2,   2.3, 0.005,   0.5,  1.6  );
	}else if( i==2 ){ // 2-track
	  func[i]->SetParameters(  0.11,   2.9,   3.0, 0.005,  0.50,   1.9 );
	}else if( i==3 ){ // 3-track
	  func[i]->SetParameters(  0.007,  0.53,  2.0, 0.13,    3.5,   3.7  );
	}else if( i==4 ){ // [45]-track
	  func[i]->SetParameters(  0.009,  0.52,  2.1, 0.12,    3.7,   4.2  );
	}else if( i==5 ){ // all-track
	  func[i]->SetParameters(  0.09,   2.2,   2.3, 0.005,   0.5,  1.6  );	  
	}
      }else{
	  func_set_parameters(sel_fun, func[i], hist[i], xbin, offset+xmin, offset+xmax);
      }
    }
  }
  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();
  c2->Draw();
  c3->Draw();
  c1->cd(1);
  gPad->SetLogy();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++

  for(Int_t i=0; i<Nhist; i++ ){
    if( flag_fit ){
      if( i==0 ) hist[i]->Fit(func[i],"LR","PE0"); // PE0same // RRRRR
      else hist[i]->Fit(func[i],"LR","PE0same");   // PE0same // RRRRR
      std::cout << "chi2/NDF = " << func[i]->GetChisquare() << " / " << func[i]->GetNDF()
		<< " = " << func[i]->GetChisquare()/func[i]->GetNDF()
		<< std::endl;
      c2->cd(i+1);
      gPad->SetLogy();
      hist[i]->Draw();
      c2->cd(Nhist+1);
      gPad->SetLogy();
      if( i==0 ) func[i]->Draw();
      else       func[i]->Draw("same");
      c3->cd(i+1);
      hist[i]->Draw();
      c3->cd(Nhist+1);
      if( i==0 ) func[i]->Draw();
      else       func[i]->Draw("same");
      c1->cd(1);
    }else hist[i]->Draw( "same" );
  }

  // +++++++ tlegend ++++++++++++++++++++++++++++++++++

  c1->cd(2);
  TPaveText* box    = new TPaveText( 0.0,0.4,1.0,1.0 );
  box->SetTextSize(0.015);
  TLegend*   legend = new TLegend  ( 0.0,0.0,1.0,0.4 );
  legend->SetTextSize(0.015);
  for( Int_t j=0; j<Nchain; j++ ){
    box->AddText( Form("tmphist%d : Mode(%s : %s : %d files)", j,chain[j]->GetModeName(),gSystem->BaseName(chain[j]->GetFileName()),chain[j]->GetNfile()) );
    box->AddText( Form("   %s", add_cut[j]) );
  }

  for( Int_t i=0; i<Nhist; i++ ){
    legend->AddEntry( hist[i],"","PL" );
  }
  box->Draw();
  legend->Draw();

  c1->Update();
  c2->Update();
  c3->Update();
  if( flag_save ){
    c1->Print(     Form("pic/pdf_kfbchi2_1_set%s_sig.eps",setname) );
    c2->Print(     Form("pic/pdf_kfbchi2_2_set%s_sig.eps",setname) );
    c3->Print(     Form("pic/pdf_kfbchi2_3_set%s_sig.eps",setname) );
    TFile outfile( Form("pic/pdf_kfbchi2_set%s_sig.root", setname), "RECREATE" );
    //c1->Write();
    //c2->Write();
    //c3->Write();
    for( Int_t i=0; i<Nhist; i++ ){
      hist[i]->Write();
      if( flag_fit ) func[i]->Write(), func[i]->Print();
    }
    outfile.Close();
  }

  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();

  delete[] infile;
  delete[] add_cut;
  delete[] chain;
  delete[] tmphist;
  delete[] hist;
  delete[] func;
  delete   box;
  delete   legend;
  delete   c1;
  delete   c2;
  delete   c3;
    
  return 0;
}
