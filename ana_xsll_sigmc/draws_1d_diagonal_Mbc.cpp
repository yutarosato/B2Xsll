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
  if( !(argc==4 || argc==5) ) std::cerr << "wrong input" << std::endl
					<< " Usage : ./draws_1d_diagonal (int)fl_mode_xs (int)fl_mode_ll (char*)setname [(int)fl_appRun]" << std::endl
    					<< "[fl_mode_xs] : four digits"           << std::endl
					<< "[fl_mode_ll] : 1(e), 0(mu)"           << std::endl, abort();

  Int_t   fl_mode_xs = atoi( argv[1] );
  Int_t   fl_mode_ll = atoi( argv[2] ); // 1(e), 0(mu)
  Char_t* setname    = argv[3]; 
  Int_t   fl_appRun  = 1;
  if( argc==5 ) fl_appRun = atoi( argv[4] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  sTmp << indir << "sigMC_*_m9999m_caseB_set[" << setname << "]";
  Char_t tmp_infile[255];
  strcpy( tmp_infile, (Char_t*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();
  const Char_t* infile = (const Char_t*)tmp_infile;  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t Nchain        = 5;
  const Int_t Nhist         = 5;
  const Int_t nfile[Nchain] = {0};
  
  const Int_t add[Nhist][Nchain] ={
    {1},
    {1,1},
    {1,1,1},
    {1,1,1,1},
    {0,0,0,0,1},
  };

  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ) add_cut[i] = new Char_t[4096];

  strcpy( add_cut[0], makeCut_mode_category(fl_mode_xs, 4).c_str() ); // other
  strcpy( add_cut[1], makeCut_mode_category(fl_mode_xs, 3).c_str() ); // off-diagonal
  strcpy( add_cut[2], makeCut_mode_category(fl_mode_xs, 0).c_str() ); // false
  strcpy( add_cut[3], makeCut_mode_category(fl_mode_xs, 1).c_str() ); // true
  strcpy( add_cut[4], makeCut_mode_category(fl_mode_xs   ).c_str() ); // total for check !

  using namespace Mbc;
  //using namespace deltaE;
  //using namespace M_Xs;
  //using namespace M_ll;
  
  const Bool_t flag_save = true; // outfile.eps and outfile.root

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
    chain[j] = new MChain( infile, tname, branch_table(), nfile[j], tail );
    nominal_cut_selection( chain[j], fl_mode_ll )( chain[j]->GetCut(), tname );
  }

  // ++++++++++++++++++++++++
  // cut change
  for( Int_t j=0; j<Nchain; j++ ){
    chain[j]->GetCut()->Set( "Mbc", 1,  5.22, 0.0, 5.30 );
    //chain[j]->GetCut()->Set( "de",  1, -0.20, 0.0, 0.20 );
  }
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
    Deco( hist[i], 3, col_fil[i], col_fil[i] );
    for( Int_t j=0; j<Nchain; j++ ){
      if( add[i][j] ) hist[i]->Add( tmphist[j] );
    }
  }

  // +++++++ display +++++++++++++++++++++++++++++++++++++++++++++++++
  Double_t entry_all[Nhist] = {0};
  Double_t entry_sig[Nhist] = {0}; // # of events in signal box region
  for( Int_t i=0; i<Nhist; i++ ){
    std::cout << Form("<hist %d> ",i);
    entry_all[i] = hist[i]->GetEntries();
    Double_t entry_canvas = hist[i]->Integral();
    Double_t entry_under  = hist[i]->GetBinContent(0);
    Double_t entry_over   = hist[i]->GetBinContent(xbin+1);
    for( Int_t j=hist[i]->FindBin(5.27+0.0000001); j<=xbin; j++ ) entry_sig[i] += hist[i]->GetBinContent(j);
    TString sTmp = "Added-Files( ";
    for( Int_t j=0; j<Nchain; j++ ) if( add[i][j] ) sTmp += j, sTmp += ",";
    sTmp += ")[";
    sTmp += entry_all[i];
    sTmp += " events ( canvas : ";
    sTmp += entry_canvas;
    sTmp += " / under : ";
    sTmp += entry_under;
    sTmp += " / over  : ";
    sTmp += entry_over;
    sTmp += " / sig  : ";
    sTmp += entry_sig[i];
    sTmp += "]";
    
    hist[i]->SetTitle( sTmp.Data() );
    std::cout << sTmp.Data() << std::endl;
  }
  if( entry_all[Nhist-1]!=entry_all[Nhist-2] ) std::cerr << "mismatch # of events -> abort()" << std::endl, abort();

  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();
  c1->cd(1);
  TH2D* waku = Waku( Nhist, hist, xlabel );
  ((TGaxis*)waku->GetXaxis())->SetMaxDigits(3);
  ((TGaxis*)waku->GetYaxis())->SetMaxDigits(3);
  waku->SetTitle( Form("%s (fl_mode_xs=%d, fl_mode_ll=%d, set=%s)", xlabel, fl_mode_xs, fl_mode_ll, setname) );
  waku->Draw();

  for(Int_t i=Nhist-2; i>=0; i-- ) hist[i]->Draw( "same" );

  //=================================================================
  TArrow* ar1 = new TArrow( 5.27, 1.02*((TAxis*)waku->GetYaxis())->GetXmax(),
			    5.29, 1.02*((TAxis*)waku->GetYaxis())->GetXmax(),
			    0.01,"<|>" );
  ar1->SetLineColor(2);
  ar1->SetFillColor(2);
  ar1->SetLineWidth(2);
  ar1->Draw();
  //=================================================================
  // +++++++ tlegend1 ++++++++++++++++++++++++++++++++++
  TLegend* legend1 = new TLegend( 0.75,0.75,0.99,0.99 );
  legend1->AddEntry( hist[3],"True",    "P" );
  legend1->AddEntry( hist[2],"False",   "P" );
  legend1->AddEntry( hist[1],"Off-diag","P" );
  legend1->AddEntry( hist[0],"Other",   "P" );
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
    legend2->AddEntry( hist[i],"","P" );
  }
  box->Draw();
  legend2->Draw();

  c1->Update();
  if( flag_save ){
    c1->Print( Form("pic/diag_%s_m%dm_lep%d_set%s.eps", axis,fl_mode_xs, fl_mode_ll, setname) );
    c1->Print( Form("pic/diag_%s_m%dm_lep%d_set%s.root",axis,fl_mode_xs, fl_mode_ll, setname) );
  }
  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();

  delete[] chain;
  delete[] tmphist;
  delete[] hist;
  delete[] func;
  delete[] add_cut;
  delete   c1;
  delete   legend1;
  delete   legend2;
  delete   ar1;
  delete   box;
    
  return 0;
}

