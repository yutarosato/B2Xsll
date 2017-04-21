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
  using namespace sigmc;
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==7 || argc==8) ) std::cerr << "wrong input" << std::endl
					<< " Usage : ./draws_1d_diagonal (int)fl_mode_ll (int)4*fl_mode_xs (char*)setname [(int)fl_appRun]" << std::endl
					<< "[fl_mode_xs] : four digits"           << std::endl
					<< "[fl_mode_ll] : 1(e), 0(mu)"           << std::endl,abort();
  Int_t   fl_mode_ll = atoi( argv[1] ); // 1(e), 0(mu)
  Char_t* setname    = argv[6];
  Int_t   fl_appRun  = 1;
  if( argc==8 ) fl_appRun = atoi( argv[7] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t Nchain        = 4;
  const Int_t Nhist         = 4;
  const Int_t nfile[Nchain] = {0};
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Int_t fl_mode_xs[Nchain];
  fl_mode_xs[0] = atoi( argv[2] );
  fl_mode_xs[1] = atoi( argv[3] );
  fl_mode_xs[2] = atoi( argv[4] );
  fl_mode_xs[3] = atoi( argv[5] );
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
    {1},
    {0,1},
    {0,0,1},
    {0,0,0,1},
  };
  
  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    add_cut[i] = new Char_t[4096];
    strcpy( add_cut[i], makeCut_mode_category(fl_mode_xs[i], 1).c_str() ); // true
  }

  const Char_t*  tname    = "h511";
  const Char_t*  axis     = "kfbchi/kfbdgf";
  const Double_t offset   = 0.0;
  const Int_t    xbin     = 50;
  const Double_t xmin     = 0;
  const Double_t xmax     = 10.0;
  const Double_t xmin_fit = xmin;
  const Double_t xmax_fit = xmax;
  const Char_t*  xlabel   = "C.L.";
  
  //using namespace Mbc;
  //using namespace deltaE;
  //using namespace M_Xs;
  //using namespace M_ll;
  
  const Bool_t flag_save  = !true; // outfile.eps and outfile.root
  const Bool_t flag_scale = true;
  
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
    chain[j]->GetCut()->Set( "Mbc", 1,  5.22, 0.0, 5.30 );
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
  for( Int_t i=0; i<Nhist; i++ ){
    hist[i] = new TH1D( Form("hist%d",i),Form("hist%d",i), xbin,offset+xmin,offset+xmax );
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

  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();
  c1->cd(1);
  gPad->SetLogy();
  //TH2D* waku = Waku( Nhist, hist, xlabel );
  //((TGaxis*)waku->GetXaxis())->SetMaxDigits(3);
  //((TGaxis*)waku->GetYaxis())->SetMaxDigits(3);
  //waku->SetTitle( Form("%s (fl_mode_ll=%d, set=%s)", xlabel, fl_mode_ll, setname) );
  //waku->Draw();
  hist[0]->Draw();
  for(Int_t i=Nhist-1; i>=0; i-- ) hist[i]->Draw( "same" );
  
  // +++++++ tlegend1 ++++++++++++++++++++++++++++++++++
  TLegend* legend1 = new TLegend( 0.75,0.75,0.99,0.99 );
  legend1->AddEntry( hist[0],Form("m%dm",fl_mode_xs[0]), "L" );
  legend1->AddEntry( hist[1],Form("m%dm",fl_mode_xs[1]), "L" );
  legend1->AddEntry( hist[2],Form("m%dm",fl_mode_xs[2]), "L" );
  legend1->AddEntry( hist[3],Form("m%dm",fl_mode_xs[3]), "L" );
  legend1->Draw();
  // +++++++ tlegend ++++++++++++++++++++++++++++++++++
  c1->cd(2);
  TPaveText* box    = new TPaveText( 0.0,0.4,1.0,1.0 );
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
    c1->Print( Form("pic/diag_%s_true_m%dm_lep%d_set%s.eps",  axis, fl_mode_xs[0], fl_mode_ll, setname) );
    c1->Print( Form("pic/diag_%s_true_m%dm_lep%d_set%s.root", axis, fl_mode_xs[0], fl_mode_ll, setname) );
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

