#include <iostream>
#include <math.h>

#include "../Util/Style.h"
#include "../Util/Canvas.h"
#include "../Util/const.h"
#include "../Util/MChain.h"
#include "../Util/MCut.h"
#include "../Util/MCut_array.h"
#include "../Util/FitFunc.h"
#include "../Set_dstr_fl/Nominal_cut_selection.h"
#include "../Set_dstr_fl/Branch.h"

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
  if( !(argc==3 || argc==4) ) std::cerr << "wrong input" << std::endl
					<< " Usage : ./draws_1d (int)fl_pid (int)fl_par [(int)fl_appRun]" << std::endl
					<< "[fl_pid] ; 1(e), 0(mu), 2(K)"       << std::endl
					<< "[fl_par] : 0(K), 1(pi), 2(slow-pi)" << std::endl, abort();

  const Int_t fl_pid    = atoi( argv[1] );
  const Int_t fl_par    = atoi( argv[2] );
  Int_t       fl_appRun = 1;
  if( argc==4 ) fl_appRun = atoi( argv[3] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  sTmp << indir << "*-e55r*"; // small sample for test
  Char_t tmp_infile[255];
  strcpy( tmp_infile, (Char_t*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();
  const Char_t* infile = (const Char_t*)tmp_infile;  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t Nchain        = 10;
  const Int_t Nhist         = 5;
  const Int_t nfile[Nchain] = {0};
  
  const Int_t add[Nhist][Nchain] ={
    {1,1,1,1},             // fake-D
    {1,1,1,1,1,1},         // fake-D(switch K-pi)
    {1,1,1,1,1,1,1,1},     // true-D + fake-slow-pion
    {1,1,1,1,1,1,1,1,1},   // true-D*
    {0,0,0,0,0,0,0,0,0,1}, // total
  };

  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ) add_cut[i] = new Char_t[4096];
  add_cut[0] = "( !(mc_0==1 && mc_1==1 && idmo_0==idmo_1 && abs(idmo_0)==421) && !(mc_0==0 && mc_1==0 && idmo_0==idmo_1 && abs(idmo_0)==421) && mc_0==0 && mc_1==0 )";                                   // fake-D (fake-K & fake-pi)
  add_cut[1] = "( !(mc_0==1 && mc_1==1 && idmo_0==idmo_1 && abs(idmo_0)==421) && !(mc_0==0 && mc_1==0 && idmo_0==idmo_1 && abs(idmo_0)==421) && mc_0==1 && mc_1==0 )";                                   // fake-D (true-K & fake-pi)
  add_cut[2] = "( !(mc_0==1 && mc_1==1 && idmo_0==idmo_1 && abs(idmo_0)==421) && !(mc_0==0 && mc_1==0 && idmo_0==idmo_1 && abs(idmo_0)==421) && mc_0==0 && mc_1==1 )";                                   // fake-D (fake-K & true-pi)
  add_cut[3] = "( !(mc_0==1 && mc_1==1 && idmo_0==idmo_1 && abs(idmo_0)==421) && !(mc_0==0 && mc_1==0 && idmo_0==idmo_1 && abs(idmo_0)==421) && mc_0==1 && mc_1==1 )";                                   // fake-D (true-K & true-pi)
  add_cut[4] = "( !(mc_0==1 && mc_1==1 && idmo_0==idmo_1 && abs(idmo_0)==421) &&  (mc_0==0 && mc_1==0 && idmo_0==idmo_1 && abs(idmo_0)==421) && !(mc_2==1 && abs(idmo_2)==413) )"; // fake-D(switch K-pi) + fake-slow-pion
  add_cut[5] = "( !(mc_0==1 && mc_1==1 && idmo_0==idmo_1 && abs(idmo_0)==421) &&  (mc_0==0 && mc_1==0 && idmo_0==idmo_1 && abs(idmo_0)==421) &&  (mc_2==1 && abs(idmo_2)==413) )"; // fake-D(switch K-pi) + true-slow-pion
  add_cut[6] = "(  (mc_0==1 && mc_1==1 && idmo_0==idmo_1 && abs(idmo_0)==421) &&  (mc_2==0) )";                     // true-D, fake-slow-pion(not pion)
  add_cut[7] = "(  (mc_0==1 && mc_1==1 && idmo_0==idmo_1 && abs(idmo_0)==421) &&  (mc_2==1 && abs(idmo_2)!=413) )"; // true-D, fake-slow-pion
  add_cut[8] = "(  (mc_0==1 && mc_1==1 && idmo_0==idmo_1 && abs(idmo_0)==421) &&  (mc_2==1 && abs(idmo_2)==413) )"; // true-D*
  add_cut[9] = "1";

  using namespace pid;
  
  const Bool_t flag_save = true; // outfile.eps and outfile.root

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain   = new MChain*[Nchain];
  TH1D**    tmphist = new TH1D*  [Nchain];
  TH1D**    hist    = new TH1D*  [Nhist];
  TF1**     func    = new TF1*   [Nhist];
  TCanvas*  c1      = Canvas( "c1","c1",1 );
  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++
  std::cout << Form(" ************************ make tmphist ( %s, %s, %d ) *************************************",tname,axis[fl_pid], fl_par) << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<infile %d > ", j );
    chain[j] = new MChain( infile, tname, branch_table(), nfile[j], tail );
    nominal_cut_selection( chain[j] )( chain[j]->GetCut(), tname );
  }

  // ++++++++++++++++++++++++
  // cut change
  for( Int_t j=0; j<Nchain; j++ ){
    //chain[j]->GetCut()->Set( "Mbc", 1,  5.22, 0.0, 5.30 );
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
    tmphist[j] = new TH1D( Form("tmphist%d",j), Form("%s",chain[j]->GetChange()), xbin[fl_pid],offset[fl_pid]+xmin[fl_pid],offset[fl_pid]+xmax[fl_pid] );
    chain[j]->GetTree()->Project( Form("tmphist%d",j), Form(axis[fl_pid], fl_par), add_cut[j] );
  }
  
  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make hist *************************************" << std::endl;
  for( Int_t i=0; i<Nhist; i++ ){
    hist[i] = new TH1D( Form("hist%d",i),Form("hist%d",i), xbin[fl_pid],offset[fl_pid]+xmin[fl_pid],offset[fl_pid]+xmax[fl_pid] );
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
    Double_t entry_over   = hist[i]->GetBinContent(xbin[fl_pid]+1);
    if(       strcmp(axis[fl_pid],"eprob_%d")==0 ) {
      for( Int_t k=hist[i]->FindBin(0.80+0.0000001); k<=xbin[fl_pid]; k++ ) entry_sig[i] += hist[i]->GetBinContent(k); // e-ID
    }else if( strcmp(axis[fl_pid],"mulmu_%d")==0 ) {
      for( Int_t k=hist[i]->FindBin(0.97+0.0000001); k<=xbin[fl_pid]; k++ ) entry_sig[i] += hist[i]->GetBinContent(k); // mu-ID
    }
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

  for( Int_t j=0; j<Nchain; j++ ){
    double tmp_entry = 0;

    if(       strcmp(axis[fl_pid],"eprob_%d")==0 ){
      for( Int_t k=tmphist[j]->FindBin(0.80+0.0000001); k<=xbin[fl_pid]; k++ ) tmp_entry += tmphist[j]->GetBinContent(k); // e-ID
    }else if( strcmp(axis[fl_pid],"mulmu_%d")==0 ) {
      for( Int_t k=tmphist[j]->FindBin(0.97+0.0000001); k<=xbin[fl_pid]; k++ ) tmp_entry += tmphist[j]->GetBinContent(k); // mu-ID
    }
    std::cout << "tmphist" << std::setw(3) << std::right << j << " : " << std::setw(12) << std::right << tmp_entry << " ( " << std::setw(12) << std::right << 100 * tmp_entry/entry_sig[Nhist-1] << " % )" << std::endl;
    }

  if( entry_sig[Nhist-1] != entry_sig[Nhist-2] ) std::cerr << "[ABORT] : wrong entry : " << entry_sig[Nhist-1] << " ?=? " << entry_sig[Nhist-2] << std::endl, abort();

  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();
  c1->cd(1);
  TH2D* waku = Waku( Nhist, hist, xlabel[fl_pid] );
  ((TGaxis*)waku->GetXaxis())->SetMaxDigits(3);
  ((TGaxis*)waku->GetYaxis())->SetMaxDigits(3);
  waku->SetTitle( Form("%s_%d", xlabel[fl_pid], fl_par) );
  waku->Draw();

  for(Int_t i=Nhist-2; i>=0; i-- ) hist[i]->Draw( "same" );

  // +++++++ tlegend1 ++++++++++++++++++++++++++++++++++
  TLegend* legend1 = new TLegend( 0.75,0.75,0.99,0.99 );
  legend1->AddEntry( hist[3],"true-D*",                "P" );
  legend1->AddEntry( hist[2],"true-D + fake-slow-#pi", "P" );
  legend1->AddEntry( hist[1],"false-D(switch K-#pi)",  "P" );
  legend1->AddEntry( hist[0],"false-D",                "P" );
  legend1->Draw();

  c1->Update();
  if( flag_save ){
    c1->Print( Form("pic/pid_%s_%d.eps", fname[fl_pid], fl_par) );
    c1->Print( Form("pic/pid_%s_%d.root",fname[fl_pid], fl_par) );
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
    
  return 0;
}

