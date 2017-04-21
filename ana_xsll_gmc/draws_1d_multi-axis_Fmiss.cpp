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
  using namespace gmc;  
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==4 || argc==5) ) std::cerr << "wrong input" << std::endl
    					<< " Usage : ./draws_1d_Mbc (char*)event_type (char*)stream (int)fl_mode_ll [(int)fl_appRun]" << std::endl
    					<< "[event_type] all, uds, charm, mixed, charged" << std::endl, abort();
  Char_t* event_type = argv[1];
  Char_t* stream     = argv[2];
  Int_t   fl_mode_ll = atoi(argv[3]);
  Int_t   fl_appRun  = 1;
  if( argc==5 ) fl_appRun = atoi( argv[4] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t Nchain        = 6;
  const Int_t Nhist         = 6;
  const Int_t nfile[Nchain] = {0};
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  if( strcmp((const char*)event_type,"all") ) sTmp << indir << "/gMC_" << event_type << "_*_s0[" << stream << "]";
  else sTmp << indir << "/gMC_*_s0[" << stream << "]";  
  Char_t tmp_infile[512];
  strcpy( tmp_infile, (Char_t*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();
  const Char_t* infile = (const Char_t*)tmp_infile;
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t add[Nhist][Nchain] ={
    {1},
    {0,1},
    {0,0,1},
    {0,0,0,1},
    {0,0,0,0,1},
    {0,0,0,0,0,1},
  };
  
  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ) add_cut[i] = new Char_t[4096];
  add_cut[0] = "";
  add_cut[1] = "";
  add_cut[2] = "";
  add_cut[3] = "";
  add_cut[4] = "";
  add_cut[5] = "";

  //using namespace Evis_multi;
  using namespace Mmiss_multi;
  //using namespace Mmiss_multi_sign;
  //using namespace Deroe_multi;

  const Bool_t flag_save = true; // outfile.eps and outfile.root

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain   = new MChain*[Nchain];
  TH1D**    tmphist = new TH1D*  [Nchain];
  TH1D**    hist    = new TH1D*  [Nhist];
  TF1**     func    = new TF1*   [Nhist];
  TCanvas*  c1      = Canvas( "c1","c1",2 );
  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++
  std::cout << Form(" ************************ make tmphist ( %s, %s ) *************************************",tname,axis[0]) << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<infile %d > ", j );
    chain[j] = new MChain( infile, tname, branch_table(), nfile[j], tail );
    nominal_cut_selection( chain[j], fl_mode_ll )( chain[j]->GetCut(), tname );
  }

  // ++++++++++++++++++++++++
  // cut change
  for( Int_t j=0; j<Nchain; j++ ){
    chain[j]->GetCut()->Set( "Mbc", 1,  5.20, 0.0, 5.30 );
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
    chain[j]->GetTree()->Project( Form("tmphist%d",j), axis[j], add_cut[j] );
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
  Double_t entry_all[Nhist] = {0};
  Double_t entry_sig[Nhist] = {0}; // # of events in signal box region
  for( Int_t i=0; i<Nhist; i++ ){
    std::cout << Form("<hist %d> ",i);
    entry_all[i] = hist[i]->GetEntries();
    Double_t entry_canvas = hist[i]->Integral();
    Double_t entry_under  = hist[i]->GetBinContent(0);
    Double_t entry_over   = hist[i]->GetBinContent(xbin+1);
    for( Int_t j=hist[i]->FindBin(5.27+0.00000001); j<=xbin; j++ ) entry_sig[i] += hist[i]->GetBinContent(j);
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
  waku->SetTitle( Form("%s (%s, lep=%d,stream=%s)", axis[0], event_type, fl_mode_ll,stream) );
  waku->Draw();

  for(Int_t i=0; i<Nhist; i++ ) hist[i]->Draw( "same" );
  // +++++++ tlegend1 ++++++++++++++++++++++++++++++++++
  TLegend* legend1 = new TLegend( 0.75,0.75,0.99,0.99 );
  for(Int_t i=0; i<Nhist; i++ ) legend1->AddEntry( hist[i],Form("%s",axis[i]), "L" );
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
    legend2->AddEntry( hist[i],"","L" );
  }
  box->Draw();
  legend2->Draw();

  c1->Update();
  if( flag_save ){
    c1->Print( Form("pic/%s_%s_lep%d_s0%s.eps",  axis[0], event_type, fl_mode_ll, stream) );
    c1->Print( Form("pic/%s_%s_lep%d_s0%s.root", axis[0], event_type, fl_mode_ll, stream) );
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
