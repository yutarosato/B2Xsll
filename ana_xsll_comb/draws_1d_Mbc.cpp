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
  using namespace sig_gmc_rd;
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==7 || argc==8) ) std::cerr << "wrong input" << std::endl
    					<< " Usage : ./draws_1d_Mbc (char*)stream (char*) setname (double)used_nstream (double)used_nset (int)fl_mode_ll (int)fl_q2 [(int)fl_appRun]" << std::endl
					<< std::endl, abort();
  Char_t* stream        = argv[1];
  Char_t* setname       = argv[2];
  Double_t used_nstream = atof(argv[3]); 
  Double_t used_nset    = atof(argv[4]); 
  Int_t   fl_mode_ll    = atoi(argv[5]);
  Int_t   fl_q2         = atoi(argv[6]);
  Int_t   fl_appRun     = 1;
  if( argc==8 ) fl_appRun = atoi( argv[7] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t    fl_5body_veto   = !true;
  const Int_t    Nchain          = nbkgtype+3; // [uds,charm,mixed,charged] [false except (t,t) , false with (t,t), true]
  const Int_t    Nhist           = 5;          // [qq,bb] [false except (t,t) , false with (t,t), true]
  const Int_t    nfile[Nchain]   = {0};
  const Int_t    fl_sb[Nchain]   = {0,0,0,0,1,1,1}; // 0(bkg), 1(sig)
  const Double_t scale_event_sig = sigmc_amount*(used_nset/(Double_t)nset); // sigmc : N -> N/alpha
  const Double_t scale_event_bkg = used_nstream;                            //   gmc : N -> N/alpha
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Char_t* LRNB_cut = new Char_t[4096];
  LRNB_cut = "1";
  //strcpy( LRNB_cut, (char*)makeCut_LRNB_1d    ("lr1_%s",                             fl_q2, 0.85, 0.94                                    ).c_str() ); // 1d LR
  //strcpy( LRNB_cut, (char*)makeCut_LRNB_2d    ("lr1_%s",                             fl_q2, 0.83, 0.56, 0.86, 0.93                        ).c_str() ); // 2d LR    (bcs=bb)
  //strcpy( LRNB_cut, (char*)makeCut_LRNB_1d    ("lr1_%s",                             fl_q2, 0.87, 0.95                                    ).c_str() ); // 1d LR             (remove deroe)
  //strcpy( LRNB_cut, (char*)makeCut_LRNB_2d    ("lr1_%s",                             fl_q2, 0.78, 0.58, 0.83, 0.93                        ).c_str() ); // 2d LR    (bcs=bb) (remove deroe)

  //strcpy( LRNB_cut, (char*)makeCut_LRNB_1d    ("nb_orgksfw_vtxcl_fmiss1_%s",         fl_q2, 0.68, 0.89                                    ).c_str() ); // 1d NB
  //strcpy( LRNB_cut, (char*)makeCut_LRNB_2d    ("nb_orgksfw_vtxcl_fmiss1_%s",         fl_q2, 0.92, 0.45, 0.94, 0.87                        ).c_str() ); // 2d NB    (bcs=bb)
  //strcpy( LRNB_cut, (char*)makeCut_LRNB_2d    ("nb_orgksfw_vtxcl_fmiss1_%s",         fl_q2, 0.88, 0.45, 0.94, 0.91                        ).c_str() ); // 2d NB    (bcs=bb)
  //strcpy( LRNB_cut, (char*)makeCut_LRNB_1d_lep("nb_lep%d_orgksfw_vtxcl_fmiss1_%s",   fl_q2, 0.72, 0.89, 0.59, 0.88                        ).c_str() ); // 1d NB_lep
  //strcpy( LRNB_cut, (char*)makeCut_LRNB_2d_lep("nb_lep%d_orgksfw_vtxcl_fmiss1_%s",   fl_q2, 0.91, 0.36, 0.94, 0.92, 0.87, 0.57, 0.92, 0.89).c_str() ); // 2d NB_lep(bcs=bb)
  
  //strcpy( LRNB_cut, (char*)makeCut_LRNB_1d    ("nb_madeksfw_vtxcl_fmiss1_%s",        fl_q2, 0.67, 0.87                                    ).c_str() ); // 1d NB madeksfw
  //strcpy( LRNB_cut, (char*)makeCut_LRNB_2d    ("nb_madeksfw_vtxcl_fmiss1_%s",        fl_q2, 0.90, 0.56, 0.93, 0.88                        ).c_str() ); // 2d NB madeksfw   (bcs=bb)
  //strcpy( LRNB_cut, (char*)makeCut_LRNB_1d_lep("nb_lep%d_madeksfw_vtxcl_fmiss1_%s",  fl_q2, 0.66, 0.90, 0.76, 0.92                        ).c_str() ); // 1d NB_lep(madeksfw)
  //strcpy( LRNB_cut, (char*)makeCut_LRNB_2d_lep("nb_lep%d_madeksfw_vtxcl_fmiss1_%s",  fl_q2, 0.86, 0.42, 0.93, 0.89, 0.91, 0.67, 0.94, 0.90).c_str() ); // 2d NB_lep(madeksfw,bcs=bb)

  //strcpy( LRNB_cut, (char*)makeCut_LRNB_1d_lep("nb_lep%d_madeksfw_vtxcl_fmiss1_%s",  fl_q2, 0.91, 0.94, 0.88, 0.93                         ).c_str() ); // 1d NB_lep2 madeksfw
  //strcpy( LRNB_cut, (char*)makeCut_LRNB_2d_lep("nb_lep%d_madeksfw_vtxcl_fmiss1_%s",   fl_q2, 0.95, 0.93, 0.95, 0.95, 0.96, 0.92, 0.96, 0.94).c_str() ); // 2d NB_lep2 madeksfw (bcs=bb)
  
  //strcpy( LRNB_cut, (char*)makeCut_LRNB_2d_lep("nb_lep%d_orgksfw_vtxcl_fmiss1_%s",   fl_q2, 0.91, 0.37, 0.94, 0.92, 0.85, 0.55, 0.93, 0.88).c_str() ); // 2d NB_lep(bcs=bb)
  //strcpy( LRNB_cut, (char*)makeCut_LRNB_2d_lep("nb_lep%d_orgksfw_vtxcl_fmiss1_%s",   fl_q2, 0.89, 0.39, 0.93, 0.91, 0.86, 0.56, 0.92, 0.87).c_str() ); // 2d NB_lep(bcs=bb) // after calibration(dzll<150, Mbc>5.20)
  //strcpy( LRNB_cut, (char*)makeCut_LRNB_2d_lep("nb_lep%d_orgksfw_vtxcl_fmiss1_%s",   fl_q2, 0.76, 0.36, 0.93, 0.91, 0.86, 0.56, 0.92, 0.87).c_str() ); // 2d NB_lep(bcs=bb) // after calibration(dzll<190, Mbc>5.20)
  //strcpy( LRNB_cut, (char*)makeCut_LRNB_2d_lep("nb_lep%d_orgksfw_vtxcl_fmiss1_%s",   fl_q2, 0.91, 0.39, 0.93, 0.91, 0.86, 0.56, 0.92, 0.87).c_str() ); // 2d NB_lep(bcs=bb) // after calibration(dzll<150, Mbc>5.22)
  //strcpy( LRNB_cut, (char*)makeCut_LRNB_2d_lep("nb_lep%d_orgksfw_vtxcl_fmiss1_%s",   fl_q2, 0.91, 0.39, 0.93, 0.91, 0.86, 0.56, 0.92, 0.87).c_str() ); // 2d NB_lep(bcs=bb) // after calibration(dzll<190, Mbc>5.22)

  std::cout << "LRNB cut : " << LRNB_cut << std::endl;
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  Char_t** infile = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    infile[i] = new Char_t[1024];
    if( fl_sb[i] ) sTmp << indir[fl_sb[i]] << "sigMC_*_m9999m_caseB_set[" << setname << "]";      // sig
    else           sTmp << indir[fl_sb[i]] << "gMC_" << bkgtype[i] << "_*_s0[" << stream << "]";  // bkg
    strcpy( infile[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t add[Nhist][Nchain] ={
    {1,1},           // bkg(qq)
    {1,1,1,1},       // bkg(bb)
    {1,1,1,1,1},     // false except(t,t)
    {1,1,1,1,1,1},   // false with(t,t)
    {1,1,1,1,1,1,1}, // true
  };

  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ) add_cut[i] = new Char_t[4096];
  for( Int_t i=0; i<nbkgtype; i++ ){
    sTmp << LRNB_cut;
    if( fl_5body_veto ) sTmp << " && ( rm_xs!=401 && rm_xs!=410 && rm_xs!=1301 && rm_xs!=1310 )"; // remove 5-body modes
    strcpy( add_cut[i], (char*)sTmp.str().c_str() ); // bkg
    sTmp.str("");
    sTmp.clear();
  }
  sTmp << LRNB_cut << " && self!=1 &&";
  sTmp << makeCut_q2fl( 1,1,1 );
  if( fl_5body_veto ) sTmp << " && ( rm_xs!=401 && rm_xs!=410 && rm_xs!=1301 && rm_xs!=1310 )"; // remove 5-body modes
  strcpy( add_cut[4], (char*)sTmp.str().c_str() ); // false except q2=true, fl=true
  sTmp.str("");
  sTmp.clear();
  
  sTmp << LRNB_cut << " && self!=1 &&";
  sTmp << makeCut_q2fl( 1,1 );
  if( fl_5body_veto ) sTmp << " && ( rm_xs!=401 && rm_xs!=410 && rm_xs!=1301 && rm_xs!=1310 )"; // remove 5-body modes
  strcpy( add_cut[5], (char*)sTmp.str().c_str() ); // false with q2=true, fl=true
  sTmp.str("");
  sTmp.clear();
  
  sTmp << LRNB_cut << " && self==1";
  if( fl_5body_veto ) sTmp << " && ( rm_xs!=401 && rm_xs!=410 && rm_xs!=1301 && rm_xs!=1310 )"; // remove 5-body modes
  strcpy( add_cut[6], (char*)sTmp.str().c_str() ); // true
  sTmp.str("");
  sTmp.clear();

  using namespace Mbc_comb;
  
  const Bool_t flag_save  = true; // outfile.eps and outfile.root

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
    chain[j]->GetTree()->Project( Form("tmphist%d",j), axis, add_cut[j] );
    //tmphist[j]->Sumw2();
    if( fl_sb[j] ) tmphist[j]->Scale( 1/scale_event_sig );
    else           tmphist[j]->Scale( 1/scale_event_bkg );
  }
  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make hist *************************************" << std::endl;
  for( Int_t i=0; i<Nhist; i++ ){
    hist[i] = new TH1D( Form("hist%d",i),Form("hist%d",i), xbin,offset+xmin,offset+xmax );
    //Deco( hist[i], 3, col_fil[i], col_fil[i] );
    for( Int_t j=0; j<Nchain; j++ ){
      if( add[i][j] ) hist[i]->Add( tmphist[j] );
    }
  }
  
  if( fl_mode_ll==1 ){
    Deco( hist[0], 3, 16, 16 ); // light green
    Deco( hist[1], 3,  4,  4 ); // green
    Deco( hist[2], 3,  9,  9 ); // yellow
    Deco( hist[3], 3, 15, 15 ); // light blue
    Deco( hist[4], 3,  3,  3 ); // blue
  }else if( fl_mode_ll==0 ){
    Deco( hist[0], 3, 16, 16 ); // light green
    Deco( hist[1], 3,  4,  4 ); // green
    Deco( hist[2], 3,  9,  9 ); // yelloe
    Deco( hist[3], 3, 13, 13 ); // light red
    Deco( hist[4], 3,  2,  2 ); // red
  }
  // +++++++ display ++++++++++++++++++++++++++++++++++
  std::cout << "=================================" << std::endl
	    << "[Scale Factor]"                    << std::endl
    	    << "sig : " << scale_event_sig
	    << "("      << used_nset << " set)"    << std::endl
	    << "bkg : " << scale_event_bkg         << std::endl
	    << "=================================" << std::endl;
  Double_t entry_sig[Nhist] = {0}; // # of events in signal box region
  for( Int_t i=0; i<Nhist; i++ ){
    std::cout << Form("<hist %d> ",i);
    Double_t entry_all    = hist[i]->GetEntries();
    Double_t entry_canvas = hist[i]->Integral();
    Double_t entry_under  = hist[i]->GetBinContent(0);
    Double_t entry_over   = hist[i]->GetBinContent(xbin+1);
    for( Int_t j=hist[i]->FindBin(5.27+0.000001); j<=xbin; j++ ) entry_sig[i] += hist[i]->GetBinContent(j);
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
    sTmp += " / sig  : ";
    sTmp += entry_sig[i];
    sTmp += "]";
    
    hist[i]->SetTitle( sTmp.Data() );
    std::cout << sTmp.Data() << std::endl;
  }
  std::cout << "S = " << std::setw(8) << std::right << entry_sig[4] - entry_sig[2] << ", "
	    << "B = " << std::setw(8) << std::right << entry_sig[2]
	    << " (= " << std::setw(8) << std::right << entry_sig[0]                << " [qq], "
	    << " + "  << std::setw(8) << std::right << entry_sig[1] - entry_sig[0] << " [bb], "
	    << " + "  << std::setw(8) << std::right << entry_sig[2] - entry_sig[1] << " [false except(t,t)] ) "
	    << std::endl
	    << "Significance = " << (entry_sig[4]-entry_sig[2])/sqrt(entry_sig[4]) << std::endl;

  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();
  c1->cd(1);
  TH2D* waku = Waku( Nhist, hist, xlabel );
  ((TGaxis*)waku->GetXaxis())->SetMaxDigits(3);
  ((TGaxis*)waku->GetYaxis())->SetMaxDigits(3);
  waku->SetTitle( Form("%s (lep=%d,q2=%d,stream=%s,set=%s)", axis,fl_mode_ll,fl_q2,stream,setname) );
  waku->Draw();
  for(Int_t i=Nhist-1; i>=0; i-- ) hist[i]->Draw( "same" );
  // +++++++ tlegend1 ++++++++++++++++++++++++++++++++++
  TLegend* legend1 = new TLegend( 0.75,0.75,0.99,0.99 );
  legend1->AddEntry( hist[4],"Sig.(true)",                                 "P" );
  legend1->AddEntry( hist[3],"Sig.(false with correct q^{2} and flavor",   "P" );
  legend1->AddEntry( hist[2],"Sig.(false except correct q^{2} and flavor", "P" );
  legend1->AddEntry( hist[1],"Bkg.(bb)",                                   "P" );
  legend1->AddEntry( hist[0],"Bkg.(qq)",                                   "P" );
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
    legend2->AddEntry( hist[i],"","PL" );
  }
  box->Draw();
  legend2->Draw();

  c1->Update();
  if( flag_save ){
    c1->Print( Form("pic/%s_comb_lep%d_q2_%d_s0%s_set%s.eps",  axis, fl_mode_ll, fl_q2, stream, setname) );
    c1->Print( Form("pic/%s_comb_lep%d_q2_%d_s0%s_set%s.root", axis, fl_mode_ll, fl_q2, stream, setname) );
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
