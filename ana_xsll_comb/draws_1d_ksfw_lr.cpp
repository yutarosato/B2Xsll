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

Int_t main( Int_t argc, Char_t** argv ){
  using namespace sig_bkg;  
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==5 || argc==6) ) std::cerr << "wrong input" << std::endl
    					<< " Usage : ./draws_1d_ksfw_lr (char*)stream (char*) setname (int)fl_mode_ll (int)fl_xs_region [(int)fl_appRun]" << std::endl
					<< " fl_xs_region : 0(total) -1(low-Xs) 1(high-Xs)"
					<< std::endl, abort();
  Char_t* stream        = argv[1];
  Char_t* setname       = argv[2];
  Int_t   fl_mode_ll    = atoi(argv[3]);
  Int_t   fl_xs         = atoi(argv[4]);
  Int_t   fl_appRun     = 1;
  if( argc==6 ) fl_appRun = atoi( argv[5] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t    Nchain          = nbkgtype+3; // [uds,charm,mixed,charged] [false except (t,t) , false with (t,t), true]
  const Int_t    Nhist           = 5;          // [qq,bb] [false except (t,t) , false with (t,t), true]
  const Int_t    nfile[Nchain]   = {0};
  const Int_t    fl_sb[Nchain]   = {0,0,0,0,1,1,1}; // 0(bkg), 1(sig)
  const Bool_t   flag_scale      = true;
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  Char_t** infile = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    infile[i] = new Char_t[1024];
    if( fl_sb[i] ) sTmp << indir[fl_sb[i]] << "sigMC_*_m9999m_caseB_set[" << setname << "]";      // sig
    else           sTmp << indir[fl_sb[i]] << "gMC_" << bkgtype[i] << "_*_s0[" << stream << "]";  // bkg
    strcpy( infile[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t add[Nhist][Nchain] ={
    {1,1},           // bkg(qq)
    {0,0,1,1},       // bkg(bb)
    {0,0,0,0,1},     // false except(t,t)
    {0,0,0,0,0,1},   // false with(t,t)
    {0,0,0,0,0,0,1}, // true
  };

  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ) add_cut[i] = new Char_t[4096];
  //for( Int_t i=0; i<nbkgtype; i++ ) add_cut[i] = "lpgt!=3 || lmgt!=3";
  for( Int_t i=0; i<nbkgtype; i++ ){
    if     ( fl_xs>0 ) add_cut[i] = "xs_m>1.1";
    else if( fl_xs<0 ) add_cut[i] = "xs_m<1.1";
    else               add_cut[i] = "1";
  }

  sTmp << "self!=1 &&";
  sTmp << makeCut_q2fl( 1,1,1 );
  if     ( fl_xs>0 ) sTmp << " && xs_m>1.1";
  else if( fl_xs<0 ) sTmp << " && xs_m<1.1";
  strcpy( add_cut[4], (char*)sTmp.str().c_str() ); // false except q2=true, fl=true
  sTmp.str("");
  sTmp.clear();

  sTmp << "self!=1 &&";
  sTmp << makeCut_q2fl( 1,1 );
  if     ( fl_xs>0 ) sTmp << " && xs_m>1.1";
  else if( fl_xs<0 ) sTmp << " && xs_m<1.1";
  strcpy( add_cut[5], (char*)sTmp.str().c_str() ); // false with q2=true, fl=true
  sTmp.str("");
  sTmp.clear();

  sTmp << "self==1";
  if     ( fl_xs>0 ) sTmp << " && xs_m>1.1";
  else if( fl_xs<0 ) sTmp << " && xs_m<1.1";
  strcpy( add_cut[6], (char*)sTmp.str().c_str() ); // true
  sTmp.str("");
  sTmp.clear();

  using namespace ksfw_lr_tot;
  
  const Bool_t flag_save  = true; // outfile.eps and outfile.root

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain   = new MChain*[Nchain];
  TH1D**    tmphist = new TH1D*  [Nchain];
  TH1D**    hist    = new TH1D*  [Nhist];
  TF1**     func    = new TF1*   [Nhist];
  TCanvas*  c1      = Canvas( "c1","c1",2 );
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
    //chain[j]->GetCut()->Set( "Mbc", 1,  5.20, 0.0, 5.30 );
    //chain[j]->GetCut()->Set( "de",  1, -0.50, 0.0, 0.50 );
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
    hist[i] = new TH1D( Form("hist%d",i),Form("hist%d",i), xbin,offset+xmin,offset+xmax );
    Deco( hist[i], 1, i+1, i+1 );
    for( Int_t j=0; j<Nchain; j++ ){
      if( add[i][j] ) hist[i]->Add( tmphist[j] );
    }
  }
  
  // +++++++ display ++++++++++++++++++++++++++++++++++
  Double_t entry_sig[Nhist] = {0}; // # of events in signal box region
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
  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();
  c1->cd(1);
  TH2D* waku = Waku( Nhist, hist, xlabel );
  ((TGaxis*)waku->GetXaxis())->SetMaxDigits(3);
  ((TGaxis*)waku->GetYaxis())->SetMaxDigits(3);
  waku->SetTitle( Form("%s (lep=%d,xs=%d,stream=%s,set=%s)", axis,fl_mode_ll,fl_xs,stream,setname) );
  waku->Draw();
  for(Int_t i=0; i<Nhist; i++ ) hist[i]->Draw( "same" );
  // +++++++ tlegend1 ++++++++++++++++++++++++++++++++++
  TLegend* legend1 = new TLegend( 0.75,0.75,0.99,0.99 );
  legend1->AddEntry( hist[4],"Sig.(true)",                                 "P" );
  legend1->AddEntry( hist[3],"Sig.(false with correct q^{2} and flavor",   "P" );
  legend1->AddEntry( hist[2],"Sig.(false except correct q^{2} and flavor", "P" );
  legend1->AddEntry( hist[1],"Bkg.(bb)",                                   "P" );
  legend1->AddEntry( hist[0],"Bkg.(qq)",                                   "P" );
  legend1->Draw();
  // +++++++ tlegend2 ++++++++++++++++++++++++++++++++++
  c1->cd(2);
  TPaveText* box     = new TPaveText( 0.0,0.4,1.0,1.0 );
  box->SetTextSize(0.015);
  TLegend*   legend2 = new TLegend  ( 0.0,0.0,1.0,0.4 );
  legend2->SetTextSize(0.015);
  for( Int_t j=0; j<Nchain; j++ ){
    box->AddText( Form("tmphist%d : Mode(%s : %s : %d files)", j,chain[j]->GetModeName(),gSystem->BaseName(chain[j]->GetFileName()),chain[j]->GetNfile()) );
    box->AddText( Form("   %s", add_cut[j]) );
  }

  for( Int_t i=0; i<Nhist; i++ ){
    legend2->AddEntry( hist[i],"","PL" );
  }
  box->Draw();
  legend2->Draw();

  c1->Update();
  if( flag_save ){
    c1->Print( Form("pic/%s_lep%d_xs_%d_s0%s_set%s.eps",  fname, fl_mode_ll, fl_xs, stream, setname) );
    c1->Print( Form("pic/%s_lep%d_xs_%d_s0%s_set%s.root", fname, fl_mode_ll, fl_xs, stream, setname) );
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

