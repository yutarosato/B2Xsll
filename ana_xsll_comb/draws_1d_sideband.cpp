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
  using namespace sig_bkg;
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==5 || argc==6) ) std::cerr << "wrong input" << std::endl
					<< " Usage : ./draws_1d_sideband (char*)stream (char*)setname (int)fl_axis (int)fl_mode_ll [(int)fl_appRun]" << std::endl
					<< "[  stream  ] 0, 1, ..., 5, 01, 0-2"     	  << std::endl
    					<< "[  setname ] A, B, ..., U, AB, A-U"     	  << std::endl, abort();

  Char_t* stream     = argv[1];
  Char_t* setname    = argv[2];
  Int_t   fl_axis    = atoi(argv[3]);
  Int_t   fl_mode_ll = atoi(argv[4]);
  Int_t   fl_appRun  = 1;
  if( argc==6 ) fl_appRun = atoi( argv[5] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t Nchain             = 13; // [sig] [(uds, charm, mixed, charged)x(low-band,signal-box,high-band)]
  const Int_t Nhist              =  7; // [sig] [(qq, bb)x(low-band,signal-box,high-band)]
  const Int_t fl_sb[Nchain]      = {1}; // 1(sig) 0(bkg)
  const Int_t nfile[Nchain] = {0};
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  Char_t** infile = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    infile[i] = new Char_t[1024];
    if( fl_sb[i] ) sTmp << indir[1] << "/sigMC_*_m9999m_caseB_set[" << setname << "]";         // sig
    else           sTmp << indir[0] << "/gMC_" << bkgtype[(i-1)%4] << "_*_s0[" << stream << "]"; // bkg
    strcpy( infile[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t add[Nhist][Nchain] ={
    {1},                         // sig
    {0,1,1},                     // qq-bkg(signal-box)
    {0,0,0,1,1},                 // bb-bkg(signal-box)
    {0,0,0,0,0,1,1},             // qq-bkg( low-band )
    {0,0,0,0,0,0,0,1,1},         // bb-bkg( low-band )
    {0,0,0,0,0,0,0,0,0,1,1},     // qq-bkg(high-band )
    {0,0,0,0,0,0,0,0,0,0,0,1,1}, // bb-bkg(high-band )
  };
  
  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++) add_cut[i] = new Char_t[4096];
  for( Int_t i=0; i<Nchain; i++){
    if     ( i==0 )         sTmp << "self==1";              // sig
    else if( (i-1)%4 < 1.5) sTmp << "1";                    // qq
    else                    sTmp << "(lpgt!=3 || lmgt!=3)"; // bb

    if(              i<= 4 && fl_mode_ll==1 ) sTmp << " && (de>-0.10 && de< 0.05)"; // signal-box(ee)
    else if( i>=5 && i<= 8 && fl_mode_ll==1 ) sTmp << " && (de>-0.20 && de<-0.10)"; //  low-band(ee)
    else if( i>=9 && i<=12 && fl_mode_ll==1 ) sTmp << " && (de> 0.05 && de< 0.20)"; // high-band(ee)
    else if(         i<= 4 && fl_mode_ll==0 ) sTmp << " && (de>-0.05 && de< 0.05)"; // signal-box(mm)
    else if( i>=5 && i<= 8 && fl_mode_ll==0 ) sTmp << " && (de>-0.15 && de<-0.05)"; //  low-band(mm)
    else if( i>=9 && i<=12 && fl_mode_ll==0 ) sTmp << " && (de> 0.05 && de< 0.20)"; // high-band(mm)

    strcpy( add_cut[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }


  using namespace bgsup;

  const Bool_t flag_save  = true; // outfile.eps and outfile.root
  Bool_t       flag_scale = true;

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain   = new MChain*[Nchain];
  TH1D**    tmphist = new TH1D*  [Nchain];
  TH1D**    hist    = new TH1D*  [Nhist];
  TF1**     func    = new TF1*   [Nhist];
  TCanvas*  c1      = Canvas( "c1","c1",2 );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++
  std::cout << Form(" ************************ make tmphist ( %s, %s ) *************************************",tname, axis[fl_axis]) << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<infile %d > ", j );
    chain[j] = new MChain( infile[j], tname, branch_table(), nfile[j], tail );
    nominal_cut_selection( chain[j], fl_mode_ll )( chain[j]->GetCut(), tname );
  }

  // ++++++++++++++++++++++++
  // cut change
  for( Int_t j=0; j<Nchain; j++ ){
    //chain[j]->GetCut()->Set( "Mbc", 1,  5.22, 0.0, 5.30 );
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
    tmphist[j] = new TH1D( Form("tmphist%d",j), Form("%s",chain[j]->GetChange()), xbin[fl_axis],offset[fl_axis]+xmin[fl_axis],offset[fl_axis]+xmax[fl_axis] );
    chain[j]->GetTree()->Project( Form("tmphist%d",j), axis[fl_axis], add_cut[j] );
  }
  
  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make hist *************************************" << std::endl;
  hist[0] = new TH1D( Form("bgsup_sideband_sig_lep%d",fl_mode_ll), Form("bgsup_sideband_sig_lep%d",fl_mode_ll), xbin[fl_axis],offset[fl_axis]+xmin[fl_axis],offset[fl_axis]+xmax[fl_axis] );
  hist[1] = new TH1D( Form("bgsup_sideband_qq_lep%d_sigbox",   fl_mode_ll), Form("bgsup_sideband_qq_lep%d_sigbox",   fl_mode_ll), xbin[fl_axis],offset[fl_axis]+xmin[fl_axis],offset[fl_axis]+xmax[fl_axis] );
  hist[2] = new TH1D( Form("bgsup_sideband_bb_lep%d_sigbox",   fl_mode_ll), Form("bgsup_sideband_bb_lep%d_sigbox",   fl_mode_ll), xbin[fl_axis],offset[fl_axis]+xmin[fl_axis],offset[fl_axis]+xmax[fl_axis] );
  hist[3] = new TH1D( Form("bgsup_sideband_qq_lep%d_lowband",  fl_mode_ll), Form("bgsup_sideband_qq_lep%d_lowband",  fl_mode_ll), xbin[fl_axis],offset[fl_axis]+xmin[fl_axis],offset[fl_axis]+xmax[fl_axis] );
  hist[4] = new TH1D( Form("bgsup_sideband_bb_lep%d_lowband",  fl_mode_ll), Form("bgsup_sideband_bb_lep%d_lowband",  fl_mode_ll), xbin[fl_axis],offset[fl_axis]+xmin[fl_axis],offset[fl_axis]+xmax[fl_axis] );
  hist[5] = new TH1D( Form("bgsup_sideband_qq_lep%d_highband", fl_mode_ll), Form("bgsup_sideband_qq_lep%d_highband", fl_mode_ll), xbin[fl_axis],offset[fl_axis]+xmin[fl_axis],offset[fl_axis]+xmax[fl_axis] );
  hist[6] = new TH1D( Form("bgsup_sideband_bb_lep%d_highband", fl_mode_ll), Form("bgsup_sideband_bb_lep%d_highband", fl_mode_ll), xbin[fl_axis],offset[fl_axis]+xmin[fl_axis],offset[fl_axis]+xmax[fl_axis] );



  for( Int_t i=0; i<Nhist; i++ ){
    //hist[i] = new TH1D( Form("pdf_hist_de_lep%d_s0%s_bkg",fl_mode_ll,stream),Form("pdf_hist_de_lep%d_s0%s_bkg",fl_mode_ll,stream), xbin[fl_axis],offset[fl_axis]+xmin[fl_axis],offset[fl_axis]+xmax[fl_axis] );
    if( i==0 ) Deco( hist[i], 2, i+1,       i+1 );
    else       Deco( hist[i], 2, (i-1)/2+2, (i-1)/2+2 );
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
      //hist[i]->Sumw2();
      hist[i]->Scale( 1/(hist[i]->GetEntries()) );
    }
  }
  
  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();
  TH2D* waku = Waku( Nhist, hist, xlabel[fl_axis] );
  ((TGaxis*)waku->GetXaxis())->SetMaxDigits(3);
  ((TGaxis*)waku->GetYaxis())->SetMaxDigits(3);
  waku->SetTitle( Form("%s(set%s,stream%s)", axis[fl_axis],setname,stream) );

  // +++++++ tlegend1 ++++++++++++++++++++++++++++++++++
  TLegend** legend1 = new TLegend*[2];
  for( Int_t i=0; i<4; i++ ) legend1[i] = new TLegend( 0.75,0.75,0.99,0.99 );
  for( Int_t i=0; i<2; i++ ){
    legend1[i]->AddEntry( hist[0  ],"sig",             "L" );
    legend1[i]->AddEntry( hist[i+1],"bkg(signal-box)", "L" );
    legend1[i]->AddEntry( hist[i+3],"bkg( low-band )", "L" );
    legend1[i]->AddEntry( hist[i+5],"bkg(high-band )", "L" );
  }

  legend1[0]->SetHeader( "qq"  );
  legend1[1]->SetHeader( "bb"  );

  c1->cd(1);
  waku->Draw();
  hist[0]->Draw( "same" );
  hist[1]->Draw( "same" );
  hist[3]->Draw( "same" );
  hist[5]->Draw( "same" );
  legend1[0] ->Draw();

  c1->cd(2);
  waku->Draw();
  hist[0]->Draw( "same" );
  hist[2]->Draw( "same" );
  hist[4]->Draw( "same" );
  hist[6]->Draw( "same" );
  legend1[1] ->Draw();


  c1->Update();
  if( flag_save ){
    c1->Print(     Form("pic/bgsup_sideband_%s_lep%d_s0%s_set%s.eps",  fname[fl_axis], fl_mode_ll,stream,setname) );
    TFile outfile( Form("pic/bgsup_sideband_%s_lep%d_s0%s_set%s.root", fname[fl_axis], fl_mode_ll,stream,setname), "RECREATE" );
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
