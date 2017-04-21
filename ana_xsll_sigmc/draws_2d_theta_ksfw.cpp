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
#include <TH2D.h>
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
  if( !(argc==3 || argc==4) ) std::cerr << "wrong input" << std::endl
					<< " Usage : ./draws_2d_theta_ksfw (char*)setname (int)fl_axis [(int)fl_appRun]" << std::endl
					<< "[  stream  ] 0, 1, ..., 5, 01, 0-2"     	  << std::endl, abort();

  Char_t* setname    = argv[1];
  Int_t   fl_axis    = atoi(argv[2]);
  Int_t   fl_appRun  = 1;
  if( argc==4 ) fl_appRun = atoi( argv[3] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t Nchain             = 2; // e, mu
  const Int_t Nhist              = 3; // e, mu, e+mu
  const Int_t fl_mode_ll[Nchain] = {1,0}; // 1(e), 0(mu)
  const Int_t nfile[Nchain] = {0};
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  Char_t** infile = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    infile[i] = new Char_t[1024];
    sTmp << indir << "/sigMC_*_m9999m_caseB_set[" << setname << "]";
    strcpy( infile[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t add[Nhist][Nchain] ={
    {1},   // e
    {0,1}, // mu
    {1,1}, // e+mu
  };
  
  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++) add_cut[i] = new Char_t[4096];
  add_cut[0] = "self==1";
  add_cut[1] = "self==1";

  using namespace ksfw;
  //using namespace ksfw_fsp;

  const Bool_t flag_save  = true; // outfile.eps and outfile.root
  Bool_t       flag_scale = !true;
  
  const Char_t*  yaxis  = "coslp";
  const Int_t    ybin   =    50;
  const Double_t ymin   = -1.00;
  const Double_t ymax   =  1.00;
  const Char_t*  ylabel = "cos#theta";

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain   = new MChain*[Nchain];
  TH2D**    tmphist = new TH2D*  [Nchain];
  TH2D**    hist    = new TH2D*  [Nhist];
  TF1**     func    = new TF1*   [Nhist];
  TCanvas*  c1      = Canvas( "c1","c1",4 );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++
  std::cout << Form(" ************************ make tmphist ( %s, %s v.s. %s ) *************************************",tname, axis[fl_axis],yaxis) << std::endl;
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
    tmphist[j] = new TH2D( Form("tmphist%d",j), Form("%s",chain[j]->GetChange()), xbin[fl_axis],offset[fl_axis]+xmin[fl_axis],offset[fl_axis]+xmax[fl_axis], ybin,ymin,ymax );
    chain[j]->GetTree()->Project( Form("tmphist%d",j), Form("%s:%s",yaxis,axis[fl_axis]), add_cut[j] );
  }
  
  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make hist *************************************" << std::endl;
  for( Int_t i=0; i<Nhist; i++ ){
    hist[i] = new TH2D( Form("hist%d",i),Form("hist%d",i), xbin[fl_axis],offset[fl_axis]+xmin[fl_axis],offset[fl_axis]+xmax[fl_axis], ybin,ymin,ymax );
    ((TGaxis*)hist[i]->GetXaxis())->SetMaxDigits(3);
    ((TGaxis*)hist[i]->GetYaxis())->SetMaxDigits(3);
    hist[i]->GetXaxis()->CenterTitle();
    hist[i]->GetYaxis()->CenterTitle();
    hist[i]->SetXTitle( xlabel[fl_axis] );
    hist[i]->SetYTitle( ylabel );
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
  for( Int_t i=0; i<Nhist; i++ ){
  c1->cd(i+1);
  hist[i]->Draw("COLZ");
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
    c1->Print(     Form("pic/theta_ksfw_%s_set%s.eps",  fname[fl_axis], setname) );
    TFile outfile( Form("pic/theta_ksfw_%s_set%s.root", fname[fl_axis], setname), "RECREATE" );
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
