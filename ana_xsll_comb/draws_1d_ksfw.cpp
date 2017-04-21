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
					<< " Usage : ./draws_1d_ksfw (char*)stream (char*)setname (int)fl_axis [(int)fl_appRun]" << std::endl
					<< "[  stream  ] 0, 1, ..., 5, 01, 0-2"     	  << std::endl
    					<< "[  setname ] A, B, ..., U, AB, A-U"     	  << std::endl, abort();

  Char_t* stream     = argv[1];
  Char_t* setname    = argv[2];
  Int_t   fl_axis    = atoi(argv[3]);
  Int_t   fl_appRun  = 1;
  if( argc==5 ) fl_appRun = atoi( argv[4] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t Nchain             = 10; // [sig, uds, charm, mixed, charged](ee,mm)
  const Int_t Nhist              =  6; // [sig, qq, bb](ee,mm)
  const Int_t fl_sb[Nchain]      = {1,0,0,0,0,
				    1,0,0,0,0}; // 1(sig) 0(bkg)
  const Int_t fl_mode_ll[Nchain] = {1,1,1,1,1,
				    0,0,0,0,0}; // 1(e), 0(mu)
  const Int_t nfile[Nchain] = {0};
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  Char_t** infile = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    infile[i] = new Char_t[1024];
    if( fl_sb[i] ) sTmp << indir[1] << "/sigMC_*_m9999m_caseB_set[" << setname << "]";         // sig
    else           sTmp << indir[0] << "/gMC_" << bkgtype[i%5-1] << "_*_s0[" << stream << "]"; // bkg
    strcpy( infile[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t add[Nhist][Nchain] ={
    {1},                   // sig(ee)
    {0,1,1},               // bkg(qq,ee)
    {0,0,0,1,1},           // bkg(bb,ee)
    {0,0,0,0,0,1},         // sig(mm)
    {0,0,0,0,0,0,1,1},     // bkg(qq,mm)
    {0,0,0,0,0,0,0,0,1,1}, // bkg(bb,mm)
  };
  
  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++) add_cut[i] = new Char_t[4096];
  add_cut[0] = "self==1";            // sig(ee)
  add_cut[1] = "";                   // uds(ee)
  add_cut[2] = "";                   // charm(ee)
  add_cut[3] = "lpgt!=3 || lmgt!=3"; // mixed(ee)
  add_cut[4] = "lpgt!=3 || lmgt!=3"; // charged(ee)
  add_cut[5] = "self==1";            // sig(mm)
  add_cut[6] = "";                   // uds(mm)
  add_cut[7] = "";                   // charm(mm)
  add_cut[8] = "lpgt!=3 || lmgt!=3"; // mixed(mm)
  add_cut[9] = "lpgt!=3 || lmgt!=3"; // charged(mm)

  //using namespace ksfw;
  using namespace ksfw_fsp;

  const Bool_t flag_save  = true; // outfile.eps and outfile.root
  Bool_t       flag_scale = true;

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain   = new MChain*[Nchain];
  TH1D**    tmphist = new TH1D*  [Nchain];
  TH1D**    hist    = new TH1D*  [Nhist];
  TF1**     func    = new TF1*   [Nhist];
  TCanvas*  c1      = Canvas( "c1","c1",4 );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++
  std::cout << Form(" ************************ make tmphist ( %s, %s ) *************************************",tname, axis[fl_axis]) << std::endl;
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
    tmphist[j] = new TH1D( Form("tmphist%d",j), Form("%s",chain[j]->GetChange()), xbin[fl_axis],offset[fl_axis]+xmin[fl_axis],offset[fl_axis]+xmax[fl_axis] );
    chain[j]->GetTree()->Project( Form("tmphist%d",j), axis[fl_axis], add_cut[j] );
  }
  
  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make hist *************************************" << std::endl;
  hist[0] = new TH1D( "ksfw_sig_lep1", "ksfw_sig_lep1", xbin[fl_axis],offset[fl_axis]+xmin[fl_axis],offset[fl_axis]+xmax[fl_axis] );
  hist[1] = new TH1D( "ksfw_qq_lep1",  "ksfw_qq_lep1",  xbin[fl_axis],offset[fl_axis]+xmin[fl_axis],offset[fl_axis]+xmax[fl_axis] );
  hist[2] = new TH1D( "ksfw_bb_lep1",  "ksfw_bb_lep1",  xbin[fl_axis],offset[fl_axis]+xmin[fl_axis],offset[fl_axis]+xmax[fl_axis] );
  hist[3] = new TH1D( "ksfw_sig_lep0", "ksfw_sig_lep0", xbin[fl_axis],offset[fl_axis]+xmin[fl_axis],offset[fl_axis]+xmax[fl_axis] );
  hist[4] = new TH1D( "ksfw_qq_lep0",  "ksfw_qq_lep0",  xbin[fl_axis],offset[fl_axis]+xmin[fl_axis],offset[fl_axis]+xmax[fl_axis] );
  hist[5] = new TH1D( "ksfw_bb_lep0",  "ksfw_bb_lep0",  xbin[fl_axis],offset[fl_axis]+xmin[fl_axis],offset[fl_axis]+xmax[fl_axis] );


  for( Int_t i=0; i<Nhist; i++ ){
    //hist[i] = new TH1D( Form("pdf_hist_de_lep%d_s0%s_bkg",fl_mode_ll,stream),Form("pdf_hist_de_lep%d_s0%s_bkg",fl_mode_ll,stream), xbin[fl_axis],offset[fl_axis]+xmin[fl_axis],offset[fl_axis]+xmax[fl_axis] );
    Deco( hist[i], 2, i/3+1, i/3+1 );
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
    Double_t entry_over   = hist[i]->GetBinContent(xbin[fl_axis]+1);
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
  
  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();
  TH2D* waku = Waku( Nhist, hist, xlabel[fl_axis] );
  ((TGaxis*)waku->GetXaxis())->SetMaxDigits(3);
  ((TGaxis*)waku->GetYaxis())->SetMaxDigits(3);
  waku->SetTitle( Form("%s(set%s,stream%s)", axis[fl_axis],setname,stream) );

  // +++++++ tlegend1 ++++++++++++++++++++++++++++++++++
  TLegend** legend1 = new TLegend*[Nhist/2];
  for( Int_t i=0; i<Nhist/2; i++ ){
    legend1[i] = new TLegend( 0.75,0.75,0.99,0.99 );
    legend1[i]->AddEntry( hist[0],"(ee)",     "L" );
    legend1[i]->AddEntry( hist[3],"(#mu#mu)", "L" );
  }
  legend1[0]->SetHeader( "sig" );
  legend1[1]->SetHeader( "qq"  );
  legend1[2]->SetHeader( "bb"  );
  
  //for(Int_t i=0; i<Nhist; i++ ) hist[i]->Draw( "same" );
  for( Int_t i=0; i<Nhist/2; i++ ){
  c1->cd(i+1);
  waku->Draw();
  hist[i]    ->Draw( "same" );
  hist[i+3]  ->Draw( "same" );
  legend1[i] ->Draw();
  }
  // +++++++ tlegend ++++++++++++++++++++++++++++++++++
  c1->cd(4);
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
  if( flag_save ){
    c1->Print(     Form("pic/ksfw_%s_s0%s_set%s.eps",  fname[fl_axis], stream,setname) );
    TFile outfile( Form("pic/ksfw_%s_s0%s_set%s.root", fname[fl_axis], stream,setname), "RECREATE" );
    c1->Write();
    for( Int_t i=0; i<Nhist; i++ ) hist[i]->Write();
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
