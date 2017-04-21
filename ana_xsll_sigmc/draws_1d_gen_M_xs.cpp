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

using namespace sigmc;
using namespace M_Xs_gen;
const Bool_t  flag_save  = true; // outfile.eps and outfile.root
const Bool_t  flag_scale = true;

Int_t main( Int_t argc, Char_t** argv ){
  TApplication app( "app", &argc, argv );
  Style();

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==5 || argc==6) ) std::cerr << "wrong input" << std::endl
					<< " Usage : ./draws_1d_gen_M_xs (int)fl_mode_ll (int)fl_klong (char*)setname (int) used_nset [(int)fl_appRun]" << std::endl
					<< std::endl, abort();
  Int_t   fl_mode_ll = atoi( argv[1] ); // 1(e) 0(mu) 2(all)
  Int_t   fl_klong   = atoi( argv[2] ); // 1(include K_L) 0(not include K_L)
  Char_t* setname    = argv[3];
  Int_t   used_nset  = atoi(argv[4]);
  Int_t   fl_appRun  = 1;
  if( argc==6 ) fl_appRun = atoi( argv[5] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  //sTmp << indir_gen << "sigMC_*_m1m_caseB_set[" << setname << "]";
  sTmp << indir_gen << "sigMC_*_m9999m_caseB_set[" << setname << "]";
  Char_t tmp_infile[255];
  strcpy( tmp_infile, (Char_t*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();
  const Char_t* infile = (const Char_t*)tmp_infile;
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  const Int_t Nchain        = 7;
  const Int_t Nhist         = 7;
  const Int_t nfile[Nchain] = {0};
  
  const Int_t add[Nhist][Nchain] ={
    {1},
    {1,1},
    {1,1,1},
    {1,1,1,1},
    {1,1,1,1,1},
    {1,1,1,1,1,1},
    {0,0,0,0,0,0,1},
  };
  const Double_t scale_event_sig = sigmc_amount*(used_nset/(Double_t)nset); // sigmc : N -> N/alpha
  
  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    add_cut[i] = new Char_t[4096];

    // Lepton selection
    if( fl_mode_ll==0 || fl_mode_ll==1 ) sTmp << "gm_l==" << fl_mode_ll;
    else sTmp << 1;
    
    if( i==0 ){ // other mode
      sTmp << "&& (";
      if( fl_klong ){
	sTmp << "(gm_bg!=0 && gm_bg!=1000)";
	sTmp << "||";
	sTmp << "(gm_bg==1000 && gm_xs!=0 && gm_xs!=100 && gm_xs!=200 && gm_xs!=300 && gm_xs!=400 && gm_xs!=1000 && gm_xs!=1100 && gm_xs!=1200 && gm_xs!=1300)";
      }else{
	sTmp << "(gm_bg!=0)";
      }
      sTmp << "||";
      sTmp << "( gm_bg==0 && ";
      for( Int_t k=0; k<nmode; k++ ){
	sTmp << " gm_xs!=" << mode[k];
	if( k!= nmode-1 ) sTmp << " && ";
      }
      sTmp << ")";
      sTmp << ")";
      strcpy( add_cut[i], (char*)sTmp.str().c_str() );
      sTmp.str("");
      sTmp.clear();
    }else if( i==1 ){
      sTmp << "&& (";
      sTmp << "( (gm_xs==1 || gm_xs==10) && gm_bg==0    )"; // 1 body
      if( fl_klong ){
	sTmp << " || ";
	sTmp << "( gm_xs==0 && gm_bg==1000 )"; // 1 body (K_L)
      }
      sTmp << " )";
      strcpy( add_cut[i], (char*)sTmp.str().c_str() );
      sTmp.str("");
      sTmp.clear();
    }else if( i==2 ){
      sTmp << "&& (";
      sTmp << "( (gm_xs==101 || gm_xs==110 || gm_xs==1001 || gm_xs==1010) && gm_bg==0    )"; // 2 body
      if( fl_klong ){
	sTmp << " || ";
	sTmp << "( (              gm_xs==100 ||                gm_xs==1000) && gm_bg==1000 )"; // 2 body (K_L)
      }
      sTmp << " )";
      strcpy( add_cut[i], (char*)sTmp.str().c_str() );
      sTmp.str("");
      sTmp.clear();
    }else if( i==3 ){
      sTmp << "&& (";
      sTmp << "( (gm_xs==201 || gm_xs==210 || gm_xs==1101 || gm_xs==1110) && gm_bg==0    )"; // 3 body
      if( fl_klong ){
      sTmp << " || ";
      sTmp << "( (              gm_xs==200 ||                gm_xs==1100) && gm_bg==1000 )"; // 3 body (K_L)
      }
      sTmp << " )";
      strcpy( add_cut[i], (char*)sTmp.str().c_str() );
      sTmp.str("");
      sTmp.clear();
    }else if( i==4 ){
      sTmp << "&& (";
      sTmp << "( (gm_xs==301 || gm_xs==310 || gm_xs==1201 || gm_xs==1210) && gm_bg==0   )"; // 4 body
      if( fl_klong ){
	sTmp << " || ";
	sTmp << "( (              gm_xs==300 ||                gm_xs==1200) && gm_bg==1000 )"; // 4 body (K_L)
      }
      sTmp << " )";
      strcpy( add_cut[i], (char*)sTmp.str().c_str() );
      sTmp.str("");
      sTmp.clear();
    }else if( i==5 ){
      sTmp << "&& (";
      sTmp << "( (gm_xs==401 || gm_xs==410 || gm_xs==1301 || gm_xs==1310) && gm_bg==0    )"; // 5 body
      if( fl_klong ){
	sTmp << " || ";
	sTmp << "( (              gm_xs==400 ||                gm_xs==1300) && gm_bg==1000 )"; // 5 body (K_L)
      }
      sTmp << " )";
      strcpy( add_cut[i], (char*)sTmp.str().c_str() );
      sTmp.str("");
      sTmp.clear();
    }else if( i==6 ){
      strcpy( add_cut[i], (char*)sTmp.str().c_str() );
      sTmp.str("");
      sTmp.clear();
    }
  }

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
    //chain[j]->GetCut()->Set( "gm_llg*llg_m", 0 );
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
    if( flag_scale ) tmphist[j]->Scale( 1/scale_event_sig );
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

  // +++++++ display ++++++++++++++++++++++++++++++++++
  Double_t entry_all[Nhist] = {0};
  for( Int_t i=0; i<Nhist; i++ ){
    std::cout << Form("<hist %d> ",i);
    entry_all[i] = hist[i]->GetEntries();
    Double_t entry_canvas = hist[i]->Integral();
    Double_t entry_under  = hist[i]->GetBinContent(0);
    Double_t entry_over   = hist[i]->GetBinContent(xbin+1);
    TString sTmp2 = "Added-Files( ";
    for( Int_t j=0; j<Nchain; j++ ) if( add[i][j] ) sTmp2 += j, sTmp2 += ",";
    sTmp2 += ")[";
    sTmp2 += entry_all[i];
    sTmp2 += " events ( canvas : ";
    sTmp2 += entry_canvas;
    sTmp2 += " / under : ";
    sTmp2 += entry_under;
    sTmp2 += " / over  : ";
    sTmp2 += entry_over;
    sTmp2 += "]";
    
    hist[i]->SetTitle( sTmp2.Data() );
    std::cout << sTmp2.Data() << std::endl;
  }
  if( entry_all[Nhist-1]!=entry_all[Nhist-2] ) std::cerr << "mismatch # of events -> abort()" << std::endl, abort();
  
  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();
  c1->cd(1);
 
  TH2D* waku = Waku( Nhist, hist, xlabel );
  ((TGaxis*)waku->GetXaxis())->SetMaxDigits(3);
  ((TGaxis*)waku->GetYaxis())->SetMaxDigits(3);
  waku->SetTitle( Form("%s (fl_mode_ll=%d,fl_klong==%d,set==%s)", xlabel, fl_mode_ll, fl_klong, setname) );
  waku->Draw();
  
  for(Int_t i=Nhist-2; i>=0; i-- ) hist[i]->Draw( "same" );

  // +++++++ tlegend1 ++++++++++++++++++++++++++++++++++
  TLegend* legend1 = new TLegend( 0.75,0.75,0.99,0.99 );
  legend1->AddEntry( hist[0], "Other", "P" );
  legend1->AddEntry( hist[1], "1 body", "P" );
  legend1->AddEntry( hist[2], "2 body", "P" );
  legend1->AddEntry( hist[3], "3 body", "P" );
  legend1->AddEntry( hist[4], "4 body", "P" );
  legend1->AddEntry( hist[5], "5 body", "P" );
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
    c1->Print( Form("pic/%s_lep%d_klong%d_set%s.eps",  axis, fl_mode_ll,fl_klong,setname) );
    c1->Print( Form("pic/%s_lep%d_klong%d_set%s.root", axis, fl_mode_ll,fl_klong,setname) );
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

