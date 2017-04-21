#include <iostream>

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
#include <TArrow.h>


Int_t main( Int_t argc, Char_t** argv ){
  //using namespace sig_gmc_rd_cut2_beforebgsup;
  using namespace sig_gmc_rd_cut3_beforebgsup;
  //using namespace sig_gmc_rd_emu_beforebgsup;
  using namespace nb_lep;
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==4 || argc==5) ) std::cerr << "wrong input" << std::endl
					<< " Usage : ./draws_2d (char*)stream (double)used_nstream (int)fl_mode_ll [(int)fl_appRun]" << std::endl
					<< "[  stream  ] 0,1,2,3,4,5,01,0-5" << std::endl, abort();
  
  const Char_t*  stream       = argv[1];
  const Double_t used_nstream = atof(argv[2]); 
  const Int_t    fl_mode_ll   = atoi(argv[3]); // 1(e), 0(mu)
  Int_t          fl_appRun    = 1;
  if( argc==5 ) fl_appRun = atoi( argv[4] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t    Nchain          = 2; // gmc, rd
  const Int_t    Nhist           = 2; // gmc, rd
  const Int_t    fl_sbr[Nchain]  = {0,2}; // 0(gmc), 1(sigmc), 2(rd)
  const Int_t    nfile[Nchain]   = {0};
  const Double_t scale_event_bkg = used_nstream; // gmc : N -> N/alpha
  const Bool_t   flag_save       = true;
  const Bool_t   flag_k4pi       = true; // 1(veto  K4pi    modes)
  const Bool_t   flag_unflavor   = true; // 1(veto unflavor modes)
  const Bool_t   flag_ccpi0      = true; // 1(veto ccpi0 peak    )
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  Char_t** infile = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    infile[i] = new Char_t[1024];
    if     ( fl_sbr[i]==0 ) sTmp << indir[fl_sbr[i]] << "/gMC_*_s0[" << stream << "]"; // gmc
    else if( fl_sbr[i]==2 ) sTmp << indir[fl_sbr[i]] << "/RD_";                        // rd
    //if     ( fl_sbr[i]==0 ) sTmp << indir[fl_sbr[i]] << "/gMC_*_e03?_s0[" << stream << "]"; // gmc // small sample
    //else if( fl_sbr[i]==2 ) sTmp << indir[fl_sbr[i]] << "/RD_e03?_";                        // rd // small sample
    else                    std::cerr << "[ABORT] Wrong fl_sbr : " << fl_sbr[i] << std::endl, abort();
    strcpy( infile[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t add[Nhist][Nchain] = {
    {1},   // gmc
    {0,1}, // rd
  };
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t fl_q2 = 0;
  Char_t* LRNB_cut = new Char_t[2048];
  LRNB_cut = "1";
  std::cout << "LRNB cut : " << LRNB_cut << std::endl;

  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t j=0; j<Nchain; j++ ){
    add_cut[j] = new Char_t[4096];
    sTmp << " 1 ";
    if( flag_k4pi     ) sTmp << " && " << makeCut_5body_veto().c_str();
    if( flag_unflavor ) sTmp << " && " << makeCut_unflavor_veto().c_str();
    if( flag_ccpi0    ) sTmp << " && " << cut_ccpi0;
    sTmp << " && " << LRNB_cut;
    ///*
    sTmp << "&& ( " // Kll
	 << Form( " (rm_xs==  1) || " ) // K+
	 << Form( " (rm_xs== 10) || " ) // Ks
	 << Form( " (rm_xs==101 && abs(xs_m-%f)<0.050) || ", PDGmass::kstr0 ) // K+pi-
	 << Form( " (rm_xs==110 && abs(xs_m-%f)<0.050) ",    PDGmass::kstrp ) // Kspi-
	 << " )";
    //*/
    strcpy( add_cut[j], (Char_t*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  MChain**  chain   = new MChain*[Nchain];
  TH2D**    tmphist = new TH2D*  [Nchain];
  TH2D**    hist2   = new TH2D*  [Nhist];
  TH1D**    histx   = new TH1D*  [Nhist];
  TH1D**    histy   = new TH1D*  [Nhist];
  TCanvas*  c1      = Canvas( "c1","c1",2, 2 );
  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form(" ************************ make chain-tree ( %s, %s ) *************************************",tname,axis) << std::endl;
    std::cout << Form( "<infile %d > ", j );
    chain[j] = new MChain( infile[j], tname, branch_table(), nfile[j], tail );
    nominal_cut_selection( chain[j], fl_mode_ll )( chain[j]->GetCut(), tname );
  }

  // ++++++++++++++++++++++++
  // cut change
  for( Int_t j=0; j<Nchain; j++ ){
    //chain[j]->GetCut()->Set( "Mbc", 1,  5.20, 0.0, 5.26 ); // Mbc sideband
    //chain[j]->GetCut()->Set( "Mbc", 1,  5.20, 0.0, 5.30 ); // Mbc (sideband+signalbox)
  }
  // ++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ changed cut *************************************" << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<infile %d > ", j ) << std::endl;
    chain[j]->GetCut()->Display(0);
    chain[j]->MakeTree();
  }

  // +++++++ make tmphist ++++++++++++++++++++++++++++++++++
  std::cout << " ************************ make tmphist *************************************" << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    tmphist[j] = new TH2D( Form("tmphist%d",j), Form("%s",chain[j]->GetChange()), xbin,xmin,xmax, ybin,ymin,ymax );
    chain[j]->GetTree()->Project( Form("tmphist%d",j), Form(axis, fl_mode_ll, fl_mode_ll), add_cut[j] ); 
    tmphist[j]->Sumw2();
    if( j==0 ) tmphist[j]->Scale( 1/scale_event_bkg );
    // +++++++ display ++++++++++++++++++++++++++++++++++
    std::cout << std::setw(10) << Form("<infile%d -> tmphist%d>",j,j);
    std::cout << std::setw(12) << tmphist[j]->GetEntries() << " events : " 
	      << add_cut[j]  << std::endl;
  }

  // +++++++ make hist2  ++++++++++++++++++++++++++++++++++
  std::cout << " ************************ make hist *************************************" << std::endl;
  for( Int_t i=0; i<Nhist; i++ ){
    if     ( i==0 ) hist2[i] = new TH2D( Form("hist%d",i),Form("hist%d (gmc)", i), xbin,xmin,xmax, ybin,ymin,ymax );
    else if( i==1 ) hist2[i] = new TH2D( Form("hist%d",i),Form("hist%d (rd)",  i), xbin,xmin,xmax, ybin,ymin,ymax );
    ((TGaxis*)hist2[i]->GetXaxis())->SetMaxDigits(3);
    ((TGaxis*)hist2[i]->GetYaxis())->SetMaxDigits(3);
    hist2[i]->GetXaxis()->CenterTitle();
    hist2[i]->GetYaxis()->CenterTitle();
    hist2[i]->SetXTitle( xlabel );
    hist2[i]->SetYTitle( ylabel );
    for( Int_t j=0; j<Nchain; j++ ){
      if( add[i][j] ) hist2[i]->Add( tmphist[j], (Double_t) add[i][j] );
    }
    // +++++++ display ++++++++++++++++++++++++++++++++++
    std::cout << std::setw(10) << Form("<hist%d>",i);
    //std::cout << std::setw(12) << hist2[i]->GetEntries() << " events, Added-Files( ";
    std::cout << std::setw(12) << hist2[i]->Integral() << " events, Added-Files( ";
    for( Int_t j=0; j<Nchain; j++ ) if( add[i][j] ) std::cout << add[i][j]*j << ",";
    std::cout << " )" << std::endl; 
  }

  // +++++++ make hist[xy] ++++++++++++++++++++++++++++++++++
  for( Int_t i=0; i<Nhist; i++ ){
    histx[i] = hist2[i]->ProjectionX();
    histy[i] = hist2[i]->ProjectionY();
    histx[i]->Scale( 1/histx[i]->Integral() );
    histy[i]->Scale( 1/histy[i]->Integral() );
    //if( i==0 ) histx[i]->Scale( 1/scale_event_bkg );
    //if( i==0 ) histy[i]->Scale( 1/scale_event_bkg );
    Deco( histx[i], 2, i+1, i+1 );
    Deco( histy[i], 2, i+1, i+1 );
  }

  // +++++++ draw ++++++++++++++++++++++++++++++++++
  Double_t allymax_x = GetYMax( Nhist, histx);
  Double_t allymax_y = GetYMax( Nhist, histy);
  for( Int_t i=0; i<Nhist; i++ ){
    histx[i]->SetMaximum( 1.1*allymax_x );
    histy[i]->SetMaximum( 1.1*allymax_y );
    histx[i]->SetTitle( Form("%s projection", xlabel) );
    histy[i]->SetTitle( Form("%s projection", ylabel) );
  }
  c1->Draw();
  for(Int_t i=0; i<Nhist; i++ ){
    c1->cd(i+1);
    //hist2[i]->Draw( "COLZ" );
    hist2[i]->DrawNormalized( "COLZ" );
  }
  c1->cd(3);
  histx[0]->Draw("hist");
  histx[1]->Draw("same");
  c1->cd(4);
  histy[0]->Draw("hist");
  histy[1]->Draw("same");


  // +++++++ tlegend ++++++++++++++++++++++++++++++++++
  TLegend* legend1 = new TLegend( 0.75,0.75,0.99,0.99 );
  legend1->AddEntry( histx[0], "gmc", "L" );
  legend1->AddEntry( histx[1], "rd",  "L" );

  legend1->Draw();
  legend1->Draw();

  // +++++++ save +++++++++++++++++++++++++++++++++++++++
  c1->Update();

  if( flag_save ){
    c1->Print( Form("pic/2d_nb_lep%d_s0%s.eps", fl_mode_ll, stream) );

    TFile outfile( Form("pic/2d_nb_lep%d_s0%s.root", fl_mode_ll, stream), "RECREATE" );
    c1->Write();
    for( Int_t i=0; i<Nhist; i++ ) hist2[i]->Write();
    for( Int_t i=0; i<Nhist; i++ ) histx[i]->Write();
    for( Int_t i=0; i<Nhist; i++ ) histy[i]->Write();
    outfile.Close();
  }
  std::cout << "finish" << std::flush;  
  if( fl_appRun ) app.Run();
  
  delete[] chain;
  delete[] tmphist;
  delete[] hist2;
  delete[] histx;
  delete[] histy;
  delete   c1;
  delete   legend1;
  
  return 0;
}
