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
    					<< " Usage : ./draws_1d_ksfw_var (char*)stream (char*) setname (int)fl_mode_ll (int)fl_ksfw_region [(int)fl_appRun]" << std::endl
					<< std::endl, abort();
  Char_t* stream        = argv[1];
  Char_t* setname       = argv[2];
  Int_t   fl_mode_ll    = atoi(argv[3]);
  Int_t   fl_ksfw       = atoi(argv[4]);
  Int_t   fl_appRun     = 1;
  if( argc==6 ) fl_appRun = atoi( argv[5] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t Ncategory       = nbkgtype+1;  // [uds,charm,mixed,charged,sig(true)]
  const Int_t Nplot           = 3;           // [qq,bb,sig(true)]
  const Int_t Nchain          = Ncategory*5; // x [5body]
  const Int_t Nhist           = Nplot    *5; // x [5body]
  const Int_t nfile[Nchain]   = {0};
  const Int_t fl_sb[Nchain]   = {0,0,0,0,1,
				 0,0,0,0,1,
				 0,0,0,0,1,
				 0,0,0,0,1,
				 0,0,0,0,1,
  }; // 0(bkg), 1(sig)
  const Bool_t   flag_scale      = true;
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  Char_t** infile = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    infile[i] = new Char_t[1024];
    if( fl_sb[i] ) sTmp << indir[fl_sb[i]] << "sigMC_*_m9999m_caseB_set[" << setname << "]";      // sig
    else           sTmp << indir[fl_sb[i]] << "gMC_" << bkgtype[i%5] << "_*_s0[" << stream << "]";  // bkg
    strcpy( infile[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t add[Nhist][Nchain] ={
    {1,1},           // bkg(qq), 1body
    {0,0,1,1},       // bkg(bb), 1body
    {0,0,0,0,1},     // true,    1body
    {0,0,0,0,0,1,1},           // bkg(qq), 2body
    {0,0,0,0,0,0,0,1,1},       // bkg(bb), 2body
    {0,0,0,0,0,0,0,0,0,1},     // true,    2body
    {0,0,0,0,0,0,0,0,0,0,1,1},           // bkg(qq), 3body
    {0,0,0,0,0,0,0,0,0,0,0,0,1,1},       // bkg(bb), 3body
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},     // true,    3body
    {0,0,0,0,0,0,0,0,0,0,0,0,0,1,1},           // bkg(qq), 4body
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1},       // bkg(bb), 4body
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},     // true,    4body
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1},           // bkg(qq), 5body
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1},       // bkg(bb), 5body
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},     // true,    5body
  };

  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ) add_cut[i] = new Char_t[4096];
  for( Int_t i=0; i<Nchain; i++ ){
  if     ( fl_ksfw == 0 ) sTmp << "( k0mm2 <-0.5 )";
  else if( fl_ksfw == 1 ) sTmp << "( k0mm2 >=-0.5 && k0mm2 < 0.3 )";
  else if( fl_ksfw == 2 ) sTmp << "( k0mm2 >= 0.3 && k0mm2 < 1.0 )";
  else if( fl_ksfw == 3 ) sTmp << "( k0mm2 >= 1.0 && k0mm2 < 2.0 )";
  else if( fl_ksfw == 4 ) sTmp << "( k0mm2 >= 2.0 && k0mm2 < 3.5 )";
  else if( fl_ksfw == 5 ) sTmp << "( k0mm2 >= 3.5 && k0mm2 < 6.0 )";
  else if( fl_ksfw == 6 ) sTmp << "( k0mm2 >= 6.0 )";
  else std::cerr << "[ABORT] invalid ksfw region : " << fl_ksfw << std::endl, abort();

  if     ( i/5==0 ) sTmp << " && " << makeCut_body(1);
  else if( i/5==1 ) sTmp << " && " << makeCut_body(2);
  else if( i/5==2 ) sTmp << " && " << makeCut_body(3);
  else if( i/5==3 ) sTmp << " && " << makeCut_body(4);
  else if( i/5==4 ) sTmp << " && " << makeCut_body(5);
  
  if( i%5==4 ) sTmp << " && self==1";
  
  strcpy( add_cut[i], (char*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();
  }
  using namespace ksfw_var_tot;
  
  const Bool_t flag_save  = true; // outfile.eps and outfile.root

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain   = new MChain*[Nchain];
  TH1D**    tmphist = new TH1D*  [Nchain];
  TH1D**    hist    = new TH1D*  [Nhist];
  TF1**     func    = new TF1*   [Nhist];
  TCanvas*  c1      = Canvas( "c1","c1", 3, 1 );
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
    Deco( hist[i], 1, i/Nplot+1, i/Nplot+1 );
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
  TH2D* waku = Waku( Nhist, hist, xlabel );
  ((TGaxis*)waku->GetXaxis())->SetMaxDigits(3);
  ((TGaxis*)waku->GetYaxis())->SetMaxDigits(3);
  waku->SetTitle( Form("%s (lep=%d,ksfw=%d,stream=%s,set=%s)", axis,fl_mode_ll,fl_ksfw,stream,setname) );
  // +++++++ tlegend1 ++++++++++++++++++++++++++++++++++
  TLegend** legend1 = new TLegend*[Nplot];

  for(Int_t i=0; i<Nplot; i++ ){
    c1->cd(i+1);
    waku->Draw();

    legend1[i] = new TLegend( 0.75,0.75,0.99,0.99 );
    if     ( i== 0 ) legend1[i]->SetHeader("qq"       );
    else if( i== 1 ) legend1[i]->SetHeader("bb"       );
    else if( i== 2 ) legend1[i]->SetHeader("sig(true)");
    legend1[i]->AddEntry( hist[i        ],"1body", "P" );
    legend1[i]->AddEntry( hist[i+1*Nplot],"2body", "P" );
    legend1[i]->AddEntry( hist[i+2*Nplot],"3body", "P" );
    legend1[i]->AddEntry( hist[i+3*Nplot],"4body", "P" );
    legend1[i]->AddEntry( hist[i+4*Nplot],"5body", "P" );

    hist[i        ]->Draw( "same" );
    hist[i+1*Nplot]->Draw( "same" );
    hist[i+2*Nplot]->Draw( "same" );
    hist[i+3*Nplot]->Draw( "same" );
    hist[i+4*Nplot]->Draw( "same" );
    legend1[i]->Draw();
  }


  c1->Update();
  if( flag_save ){
    c1->Print( Form("pic/%s_lep%d_ksfw_%d_mode_s0%s_set%s.eps",  fname, fl_mode_ll,fl_ksfw, stream, setname) );
    c1->Print( Form("pic/%s_lep%d_ksfw_%d_mode_s0%s_set%s.root", fname, fl_mode_ll,fl_ksfw, stream, setname) );
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

