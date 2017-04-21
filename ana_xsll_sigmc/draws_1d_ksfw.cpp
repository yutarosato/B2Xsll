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
  using namespace ksfw_comb;
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==4 || argc==5) ) std::cerr << "wrong input" << std::endl
					<< " Usage : ./draws_1d_ksfw (int)fl_mode_xs (int)fl_mode_ll (char*)setname [(int)fl_appRun]" << std::endl
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
  const Int_t Nhist = 2*n_axis;
  const Int_t nfile = 0;
  
  Char_t* add_cut = new Char_t[1024];
  strcpy( add_cut, makeCut_mode_category(fl_mode_xs, 1).c_str() ); // true
 
  const Bool_t flag_save = true; // outfile.eps and outfile.root

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  TH1D**    hist    = new TH1D* [Nhist];
  TCanvas*  c1      = Canvas( "c1","c1", 18 );
  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++
  std::cout << Form(" ************************ make chain ( %s, ksfw's ) *************************************",tname) << std::endl;
  std::cout << "<infile 0 > ";
  MChain*  chain = new MChain( infile, tname, branch_table(), nfile, tail );
  nominal_cut_selection( chain, fl_mode_ll )( chain->GetCut(), tname );

  // ++++++++++++++++++++++++
  // cut change
  //chain->GetCut()->Set( "Mbc", 1,  5.22, 0.0, 5.30 );
  //chain->GetCut()->Set( "de",  1, -0.20, 0.0, 0.20 );
  // ++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ changed cut *************************************" << std::endl;
  std::cout << "<infile 0 > " << std::endl;
  chain->GetCut()->Display(0);
  chain->MakeTree();
  std::cout << "add_cut : " << add_cut << std::endl;

  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make hist *************************************" << std::endl;
  for( Int_t i=0; i<Nhist; i++ ){
    hist[i] = new TH1D( Form("hist%d",i), Form("%s",chain->GetChange()), xbin[i%n_axis],offset[i%n_axis]+xmin[i%n_axis],offset[i%n_axis]+xmax[i%n_axis] );
    //Deco( hist[i], 3, col_fil[i], col_fil[i] );
    hist[i]->SetXTitle( Form(axis[i%n_axis],i/n_axis) );
    Deco( hist[i], 1, i/n_axis+1, i/n_axis+1 );
    chain->GetTree()->Project( Form("hist%d",i), Form(axis[i%n_axis],i/n_axis), add_cut );
  }

  // +++++++ display +++++++++++++++++++++++++++++++++++++++++++++++++
  Double_t entry_all[Nhist] = {0};
  Double_t entry_sig[Nhist] = {0}; // # of events in signal box region
  for( Int_t i=0; i<Nhist; i++ ){
    std::cout << Form("<hist %d> ",i);
    std::cout << std::setw(7) << std::right << Form(axis[i%n_axis],i/n_axis);
    entry_all[i] = hist[i]->GetEntries();
    Double_t entry_canvas = hist[i]->Integral();
    Double_t entry_under  = hist[i]->GetBinContent(0);
    Double_t entry_over   = hist[i]->GetBinContent(xbin[i%n_axis]+1);
    TString sTmp = "[";
    sTmp += entry_all[i];
    sTmp += " events ( canvas : ";
    sTmp += entry_canvas;
    sTmp += " / under : ";
    sTmp += entry_under;
    sTmp += " / over  : ";
    sTmp += entry_over;
    sTmp += "]";
    
    hist[i]->SetTitle( sTmp.Data() );
    std::cout << sTmp.Data() << std::endl;
  }

  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();
  for( Int_t i=0; i<n_axis; i++ ){
    c1->cd(i+1);
    hist[i]->Draw( );  
    hist[i+n_axis]->Draw("same");  
  }

  // +++++++ tlegend1 ++++++++++++++++++++++++++++++++++
  TLegend* legend1 = new TLegend( 0.75,0.75,0.99,0.99 );
  legend1->AddEntry( hist[0],       "k0", "P" );
  legend1->AddEntry( hist[0+n_axis],"k1", "P" );
  legend1->Draw();
  c1->Update();
  if( flag_save ){
    c1->Print( Form("pic/ksfw_m%dm_lep%d_set%s.eps", fl_mode_xs, fl_mode_ll, setname) );
    c1->Print( Form("pic/ksfw_m%dm_lep%d_set%s.root",fl_mode_xs, fl_mode_ll, setname) );
  }

  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();

  delete   chain;
  delete[] hist;
  delete   add_cut;
  delete   c1;
  delete   legend1;
    
  return 0;
}

