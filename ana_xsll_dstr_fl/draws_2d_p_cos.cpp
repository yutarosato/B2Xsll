#include <iostream>
#include <math.h>
#include <string.h>

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
  if( !(argc==4 || argc==5) ) std::cerr << "wrong input" << std::endl
					<< " Usage : ./draws_2d (int)fl_pid (int)fl_par (int)fl_adef [(int)fl_appRun]" << std::endl
					<< "[fl_pid ] : 1(e-id), 0(mu-id)"         << std::endl
					<< "[fl_par ] : 0(K), 1(pi), 2(slow-pi)"   << std::endl
					<< "[fl_adef] : 0(0.4GeV), 1(0.6GeV) "     << std::endl, abort();

  
  const Int_t fl_pid    = atoi( argv[1] );
  const Int_t fl_par    = atoi( argv[2] );
  const Int_t fl_adef   = atoi( argv[3] ); // definition of momentum-axis (0) 0.4 GeV, (1) 0.6 GeV
  Int_t       fl_appRun = 1;
  if( argc==5 ) fl_appRun = atoi( argv[4] );

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  sTmp << indir << "*";
  //sTmp << indir << "*-e55r*"; // small sample for test
  //sTmp << indir << "*-e21r*"; // small sample for test
  Char_t tmp_infile[255];
  strcpy( tmp_infile, (Char_t*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();
  const Char_t* infile = (const Char_t*)tmp_infile;  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t Nchain = 2;
  const Int_t Nhist  = 2;
  const Int_t nfile  = 0;
  
  const Int_t add[Nhist][Nchain] ={
    {1},   // no-lepton-ID
    {0,1}, //    lepton-ID
  };

  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ) add_cut[i] = new Char_t[4096];
  add_cut[0] = "1";
  
  if     ( fl_pid ==0 ) sTmp << Form("mulmu_%d && much2_%d>0", fl_par, fl_par ); // mu-id;
  else if( fl_pid ==1 ) sTmp << Form("eprob_%d",               fl_par         ); // e-id;
  else                  std::cerr << "[ABORT] Invalid fl_pid : " << fl_pid << std::endl, abort();
  strcpy( add_cut[1], (Char_t*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();

  using namespace p_cos;
  
  const Bool_t flag_save = true; // outfile.eps and outfile.root

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TH2D**    tmphist = new TH2D*  [Nchain];
  TH2D**    hist    = new TH2D*  [Nhist];
  TF1**     func    = new TF1*   [Nhist];
  TCanvas*  c1      = Canvas( "c1","c1", Nhist+1, 1 );
  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++

  std::cout << "<infile> ";
  MChain* chain = new MChain( infile, tname, branch_table(), nfile, tail );
  nominal_cut_selection( chain )( chain->GetCut(), tname );
  std::cout << Form(" ************************ make tmphist ( %s, %s, %d, %d ) *************************************",tname,axis,fl_pid,fl_par) << std::endl;
  // ++++++++++++++++++++++++
  // cut change
  for( Int_t j=0; j<Nchain; j++ ){   
    //if( strcmp(axis[fl_pid],"m_d_vf")==0 ) chain[j]->GetCut()->Set( "m_d_vf", 0 );
    //if( strcmp(axis[fl_pid],"dm_vf" )==0 ) chain[j]->GetCut()->Set(  "dm_vf", 0 );
    //chain->GetCut()->Set( "xdst_vf", 0 );
  }
  // ++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ changed cut *************************************" << std::endl;
    chain->GetCut()->Display(0);
    chain->MakeTree();
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<tmphist %d > ", j ) << std::endl;
    std::cout << "add_cut : " << add_cut[j] << std::endl;
  }
  // +++++++ make tmphist ++++++++++++++++++++++++++++++++++
  for( Int_t j=0; j<Nchain; j++ ){
    tmphist[j] = new TH2D( Form("tmphist%d",j), Form("%s",chain->GetChange()), xbin,offset+xmin,offset+xmax, ybin[fl_adef],ymin[fl_adef],ymax[fl_adef] );
    chain->GetTree()->Project( Form("tmphist%d",j), Form(axis,fl_par,fl_par), add_cut[j] );
  }
  
  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make hist *************************************" << std::endl;
  for( Int_t i=0; i<Nhist; i++ ){
    if     ( i==0 ) hist[i] = new TH2D( Form("hist%d",i), Form("No-Lepton-ID, par%d, adef%d",         fl_par, fl_adef), xbin,offset+xmin,offset+xmax, ybin[fl_adef],ymin[fl_adef],ymax[fl_adef] );
    else if( i==1 ) hist[i] = new TH2D( Form("hist%d",i), Form("pid%d, par%d, adef%d",        fl_pid, fl_par, fl_adef), xbin,offset+xmin,offset+xmax, ybin[fl_adef],ymin[fl_adef],ymax[fl_adef] );
    hist[i]->GetXaxis()->CenterTitle();
    hist[i]->GetYaxis()->CenterTitle();
    hist[i]->SetXTitle( xlabel );
    hist[i]->SetYTitle( ylabel );
    for( Int_t j=0; j<Nchain; j++ ){
      if( add[i][j] ) hist[i]->Add( tmphist[j] );
    }
  }

  // +++++++ display +++++++++++++++++++++++++++++++++++++++++++++++++
  Double_t entry_all[Nhist] = {0};
  for( Int_t i=0; i<Nhist; i++ ){
    std::cout << Form("<hist %d> ",i);
    entry_all[i] = hist[i]->GetEntries();
    Double_t entry_canvas = hist[i]->Integral();
    Double_t entry_under  = hist[i]->GetBinContent(0);
    Double_t entry_over   = hist[i]->GetBinContent(xbin+1);
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
    sTmp += "]";
    
    //hist[i]->SetTitle( sTmp.Data() );
    std::cout << sTmp.Data() << std::endl;
  }

  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();
  for(Int_t i=0; i<Nhist; i++ ){
    c1->cd(i+1);
    hist[i]->SetLabelSize( 0.03, "Z" );
    hist[i]->GetZaxis()->SetTitleOffset(-1.3);
    hist[i]->Draw( "COLZ" );
  }


 TH2D* hist_ratio = new TH2D( "fake_rate", "fake_rate", xbin,offset+xmin,offset+xmax, ybin[fl_adef],ymin[fl_adef],ymax[fl_adef] );
 hist_ratio->GetXaxis()->CenterTitle();
 hist_ratio->GetYaxis()->CenterTitle();
 hist_ratio->SetXTitle( xlabel );
 hist_ratio->SetYTitle( ylabel );
 hist_ratio->SetLabelSize( 0.03, "Z" );
 hist_ratio->GetZaxis()->SetTitleOffset(-1.3);
 for( Int_t jx=0; jx<xbin; jx++ ){
   for( Int_t jy=0; jy<ybin[fl_adef]; jy++ ){
     if( hist[0]->GetBinContent(jx,jy) ) hist_ratio->SetBinContent( jx,jy,hist[1]->GetBinContent(jx,jy)/hist[0]->GetBinContent(jx,jy) );
     else hist_ratio->SetBinContent( jx,jy,0 );
   }
 }

 c1->cd(3);
 hist_ratio->Draw("COLZ");
  
  c1->Update();
  if( flag_save ){
    c1->Print( Form("pic/p_cos_pid%d_par%d_adef%d.eps", fl_pid,fl_par,fl_adef) );
    c1->Print( Form("pic/p_cos_pid%d_par%d_adef%d.root",fl_pid,fl_par,fl_adef) );
  }
  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();

  delete   chain;
  delete[] tmphist;
  delete[] hist;
  delete[] func;
  delete[] add_cut;
  delete   c1;
    
  return 0;
}
