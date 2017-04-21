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
  using namespace sigmc;
  using namespace bgsup;
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==2 || argc==3) ) std::cout << "wrong input, please input axis-name" << std::endl
					<< " Usage : ./draws_1d_overwrite (int)fl_axis [(int)fl_appRun]" << std::endl
					<< " [fl_axis] 0(evis) 1(mmiss) 2(deroe) 3(kfbchi) 4(dzll) 5(cos-theta_B)"
					<< std::endl, abort();
  Int_t fl_axis   = atoi(argv[1]);
  Int_t fl_appRun = 1;
  if( argc==3 ) fl_appRun = atoi( argv[2] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t Nhist          = 2; // can not changed
  const Int_t number         = 4;
  const Char_t* event_type   = "all";
  const Char_t* stream       = "0";
  const Char_t* indir[Nhist] = { "./pic/", // sig
				 "../ana_xsll_gmc/pic/pic_bgsup/" }; // bkg
  const Bool_t flag_save   = true; // outfile.eps and outfile.root
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  std::stringstream sTmp;
  TH1D** hist  = new TH1D*[Nhist];

  for( Int_t i=0; i<Nhist; i++ ){
    Char_t tmp_infile[1024];
    if( i==0 ){
      sTmp << "bgsup_body_" << fname[fl_axis] << "_sigmc";
      strcpy( tmp_infile, (char*)sTmp.str().c_str() );
      sTmp.str("");
      sTmp.clear();
    }else{
      sTmp << "bgsup_body_" << fname[fl_axis] << "_" << event_type << "_s" << stream << "_gmc";
      strcpy( tmp_infile, (char*)sTmp.str().c_str() );
      sTmp.str("");
      sTmp.clear();
    }
    
    TFile* fin = new TFile( Form("%s/%s.root", indir[i], tmp_infile) );
    hist[i] = (TH1D*)gDirectory->Get( Form("%s%d", tmp_infile, number) );
    Deco( hist[i], 2, i+2, i+2 );
    if( !hist[i] ) std::cout << "check hist-name : " << indir[i] << tmp_infile << " (" << number << ")" << std::endl, abort();
    //if( i==0 ) hist[i]->Draw();
    //else hist[i]->Draw("same");
    std::cout << "[" << i  << "] " << indir[i] << tmp_infile << std::endl;
  }

  const Int_t xbin = hist[0]->GetNbinsX();
  TCanvas*    c1   = Canvas( "c1","c1",2 );

  // +++++++ display ++++++++++++++++++++++++++++++++++
 
  for( Int_t i=0; i<Nhist; i++ ){
    Double_t entry_all    = hist[i]->GetEntries();
    Double_t entry_canvas = hist[i]->Integral();
    Double_t entry_under  = hist[i]->GetBinContent(0);
    Double_t entry_over   = hist[i]->GetBinContent(xbin+1);
    std::cout << "[" << i << "] " 
	      << entry_all     <<  " events ( canvas : "
	      << entry_canvas  << " / under : "
	      << entry_under   << " / over  : "
	      << entry_over    << "]" << std::endl;
  }

  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();
  c1->cd(1);

  TH2D* waku = Waku( Nhist, hist, xlabel[fl_axis] );
  ((TGaxis*)waku->GetXaxis())->SetMaxDigits(3);
  ((TGaxis*)waku->GetYaxis())->SetMaxDigits(3);
  waku->SetTitle( xlabel[fl_axis] );
  waku->Draw();

  for( Int_t i=0; i<Nhist; i++ ) hist[i]->Draw( "same" );

  // +++++++ tlegend ++++++++++++++++++++++++++++++++++

  c1->cd(2);
  TLegend* legend = new TLegend  ( 0.3,0.3,0.6,0.6 );
  //legend->SetTextSize(0.015);

  legend->SetHeader( Form("%s",xlabel[fl_axis]) );
  legend->AddEntry( hist[0],Form("sig-MC"),"PL" );
  legend->AddEntry( hist[1],Form("bkg-MC"),"PL" );

  legend->Draw();

  c1->Update();
  if( flag_save ){
    c1->Print(     Form("pic/bgsup_overwrite_%s.eps",  fname[fl_axis]) );
    TFile outfile( Form("pic/bgsup_overwrite_%s.root", fname[fl_axis]), "RECREATE" );
    c1->Write();
    for( Int_t i=0; i<Nhist; i++ ) hist[i]->Write();
    outfile.Close();
  }

  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();


  delete[] hist;
  delete   legend;
  delete   c1;

  return 0;
}

