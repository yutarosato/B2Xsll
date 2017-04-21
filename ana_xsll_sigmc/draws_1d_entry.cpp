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

Int_t main( Int_t argc, Char_t** argv ){
  using namespace sigmc;
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==4 || argc==5) ) std::cerr << "wrong input" << std::endl
					<< " Usage : ./draws_1d_entry (int)fl_mode_ll (int)fl_xsid (char*)fl_setname [(int)fl_appRun]" << std::endl
					<< "[fl_mode_ll] 1(e), 0(mu)"          << std::endl
					<< "[fl_xsid] 0(all) 1(K) 2(K*) 3(Xs)" << std::endl, abort();
  Int_t   fl_mode_ll   = atoi( argv[1] ); // 1-> e, 0-> mu, 2-> e and mu
  Int_t   fl_xsid      = atoi( argv[2] );
  Char_t* fl_setname   = argv[3];
  Int_t   fl_appRun    = 1;
  if( argc==5 ) fl_appRun = atoi( argv[4] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t Nchain          = nmode+2;
  const Int_t   Nhist         = nmode+2;
  const Int_t   nfile[Nchain] = {0};
  std::stringstream sTmp;
  sTmp << indir_gen << "sigMC_*_m9999m_";
  Char_t tmp_infile[255];
  strcpy( tmp_infile, (Char_t*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();
  const Char_t* infile = (const Char_t*)tmp_infile;

  sTmp << "*_set[" << fl_setname << "]*.root";
  Char_t tmp_tail[255];
  strcpy( tmp_tail, (Char_t*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();
  const Char_t* tail_entry = (const Char_t*)tmp_tail;
  
  const Int_t add[Nhist][Nchain] ={
    {1},
    {0,1},
    {0,0,1},
    {0,0,0,1},
    {0,0,0,0,1},
    {0,0,0,0,0,1},
    {0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
  };
  
  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ) add_cut[i] = new Char_t[4096];

  for( Int_t k=0; k<nmode; k++ ){
    if( fl_xsid == 1 )      sTmp << "( abs(Xs_id)==  311 || abs(Xs_id)==  321 ) && "; // (K )
    else if( fl_xsid == 2 ) sTmp << "( abs(Xs_id)==  313 || abs(Xs_id)==  323 ) && "; // (K*)
    else if( fl_xsid == 3 ) sTmp << "( abs(Xs_id)==30343 || abs(Xs_id)==30353 ) && "; // (Xs)
    sTmp << "gm_xs==" << mode[k] << " && gm_bg==0"; // [each mode]
    if( fl_mode_ll==0 || fl_mode_ll==1 ) sTmp << " && gm_l==" << fl_mode_ll;
    strcpy( add_cut[k], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }

  if( fl_xsid == 1 )      sTmp << "( abs(Xs_id)==  311 || abs(Xs_id)==  321 ) && "; // (K )
  else if( fl_xsid == 2 ) sTmp << "( abs(Xs_id)==  313 || abs(Xs_id)==  323 ) && "; // (K*)
  else if( fl_xsid == 3 ) sTmp << "( abs(Xs_id)==30343 || abs(Xs_id)==30353 ) && "; // (Xs)
  sTmp << "( gm_bg!=0 || ( gm_bg==0 && "; // [other modes]
  for( Int_t k=0; k<nmode; k++ ){
    sTmp << " gm_xs!=" << mode[k];
    if( k!= nmode-1 ) sTmp << " && ";
  }
  sTmp << ") )";
  if( fl_mode_ll==0 || fl_mode_ll==1 ) sTmp << " && gm_l==" << fl_mode_ll;
  strcpy( add_cut[nmode], (char*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();
  
  
  if( fl_xsid == 1 )      sTmp << "( abs(Xs_id)==  311 || abs(Xs_id)==  321 ) && "; // (K )
  else if( fl_xsid == 2 ) sTmp << "( abs(Xs_id)==  313 || abs(Xs_id)==  323 ) && "; // (K*)
  else if( fl_xsid == 3 ) sTmp << "( abs(Xs_id)==30343 || abs(Xs_id)==30353 ) && "; // (Xs)
  if( fl_mode_ll==0 || fl_mode_ll==1 ) sTmp << "gm_l==" << fl_mode_ll; // [total]
  strcpy( add_cut[nmode+1], (char*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();
 
  using namespace M_Xs_gen;
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain   = new MChain*[Nchain];
  TH1D**    tmphist = new TH1D*  [Nchain];
  TH1D**    hist    = new TH1D*  [Nhist];
  TF1**     func    = new TF1*   [Nhist];
  //TCanvas*  c1      = Canvas( "c1","c1",2 );
  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++
  std::cout << Form(" ************************ make tmphist ( %s, %s ) *************************************",tname,axis) << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<infile %d > ", j );
    chain[j] = new MChain( infile, tname, branch_table(), nfile[j], tail_entry );
    nominal_cut_selection( chain[j], fl_mode_ll )( chain[j]->GetCut(), tname );
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
    tmphist[j] = new TH1D( Form("tmphist%d",j), Form("%s",chain[j]->GetChange()), xbin,offset+xmin,offset+xmax );
    chain[j]->GetTree()->Project( Form("tmphist%d",j), axis, add_cut[j] );
    std::cout << "[test1] " << chain[j]->GetTree()->GetEntries( add_cut[j] ) << std::endl;
    std::cout << "[test2] " << chain[j]->GetTree()->GetEntries( Form( "%s&&%s", add_cut[j], "Xs_m>2.0" ) ) << std::endl;
    std::cout << "[test3] " << chain[j]->GetTree()->GetEntries( Form( "%s&&%s", add_cut[j], "Xs_m<2.0" ) ) << std::endl;
  }
  
  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make hist *************************************" << std::endl;
  for( Int_t i=0; i<Nhist; i++ ){
    hist[i] = new TH1D( Form("hist%d",i),Form("hist%d",i), xbin,offset+xmin,offset+xmax );
    Deco( hist[i], 0, i+1, i+1 );
    for( Int_t j=0; j<Nchain; j++ ){
      if( add[i][j] ) hist[i]->Add( tmphist[j] );
    }
  }

  // +++++++ display ++++++++++++++++++++++++++++++++++
  std::ofstream fout( Form("log/log_gen_entry_lep%d_gmxs%d_set%s.dat", fl_mode_ll, fl_xsid, fl_setname) );
  
  for( Int_t i=0; i<Nhist; i++ ){
    std::cout << Form("<hist %d> ",i);
    Double_t entry_all    = hist[i]->GetEntries();
    fout << std::setw(5)  << fl_mode_ll << " "
	 << std::setw(5)  << i          << " "
	 << std::setw(15) << std::setprecision(15)
	 << entry_all     << std::endl;
    Double_t entry_canvas = hist[i]->Integral();
    Double_t entry_under  = hist[i]->GetBinContent(0);
    Double_t entry_over   = hist[i]->GetBinContent(xbin+1);
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
  }
  fout << "*" << std::endl;
  fout.close();
  /*
  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();
  c1->cd(1);
  TH2D* waku = Waku( Nhist, hist, xlabel );
  ((TGaxis*)waku->GetXaxis())->SetMaxDigits(3);
  ((TGaxis*)waku->GetYaxis())->SetMaxDigits(3);
  waku->Draw();
  // (fit?)  
  for( Int_t i=0; i<Nhist; i++ ) hist[i]->Draw( "same" );

  // +++++++ tlegend ++++++++++++++++++++++++++++++++++

  c1->cd(2);
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
  */
  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();

  delete[] chain;
  delete[] tmphist;
  delete[] hist;
  delete[] func;
  delete[] add_cut;
  //delete   c1;
    
  return 0;
}
