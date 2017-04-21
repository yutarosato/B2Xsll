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
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==5 || argc==6) ) std::cerr << "wrong input" << std::endl
					<< " Usage : ./draws_1d (int)fl_mode_xs (int)fl_mode_ll (int)fl_xsid (char*)setname [(int)fl_appRun]" << std::endl
					<< "[fl_mode_xs] : four digits"                << std::endl
					<< "[fl_mode_ll] : 1(e), 0(mu)"                << std::endl
					<< "[ fl_xsid  ] : 0(all), 1(K), 2(K*), 3(Xs)" << std::endl, abort();
  Int_t   fl_mode_xs = atoi( argv[1] );
  Int_t   fl_mode_ll = atoi( argv[2] ); // 1(e), 0(mu)
  Int_t   fl_xsid    = atoi( argv[3] ); // 0(all), 1(K), 2(K*), 3(Xs)
  Char_t* setname    = argv[4]; 
  Int_t fl_appRun    = 1;
  if( argc==6 ) fl_appRun = atoi( argv[5] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  sTmp << indir << "sigMC_*_m9999m_caseB_set[" << setname << "]";
  Char_t tmp_infile[255];
  strcpy( tmp_infile, (Char_t*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();
  const Char_t* infile = (const Char_t*)tmp_infile;
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  Char_t*     LRNB_cut = new Char_t[4096];
  const Int_t fl_q2    = 0;
  //strcpy( LRNB_cut, "1" );
  strcpy( LRNB_cut, Form("%s",     cut_ccpi0                         ) );

  //strcpy( LRNB_cut, Form("%s&&(%s||%s)",     cut_ccpi0, makeCut_q2(fl_mode_ll,1).c_str(), makeCut_q2(fl_mode_ll,2).c_str()                                  ) ); // tmpppp
  //strcpy( LRNB_cut, Form("%s&& %s",          cut_ccpi0, makeCut_q2(fl_mode_ll,3).c_str()                                                                    ) ); // tmpppp
  //strcpy( LRNB_cut, Form("%s&& %s",          cut_ccpi0, makeCut_q2(fl_mode_ll,5).c_str()                                                                    ) ); // tmpppp
  //strcpy( LRNB_cut, Form("%s&&(%s||%s||%s)", cut_ccpi0, makeCut_q2(fl_mode_ll,7).c_str(), makeCut_q2(fl_mode_ll,8).c_str(), makeCut_q2(fl_mode_ll,9).c_str()) ); // tmpppp
  
  //strcpy( LRNB_cut, (char*)makeCut_LRNB_2d_lep("nb_lep%d_orgksfw_vtxcl_fmiss1_%s",   fl_q2, 0.91, 0.36, 0.94, 0.92, 0.87, 0.57, 0.92, 0.89).c_str() ); // 2d NB_lep(bcs=bb)
  //strcpy( LRNB_cut, (char*)makeCut_LRNB_1d    ("lr1_%s",                             fl_q2, 0.85, 0.94                                    ).c_str() ); // 1d LR
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t Nchain           = nmode+9;
  const Int_t   Nhist          = nmode+9;
  const Int_t   nfile[Nchain]  = {0};
  
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
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},                   // each moede(only false)
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},                 // other modes
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},               // total
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},             // true
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},           // false( q2=true,  fl=true    )
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},         // false( q2=true,  fl=unknown )
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},       // false( q2=true,  fl=false   )
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},     // false( q2=false, fl=true    )
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},   // false( q2=false, fl=unknown )
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}, // false( q2=false, fl=false   )
  };
  
  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ) add_cut[i] = new Char_t[4096];

  for( Int_t k=0; k<nmode; k++ ){
    sTmp << makeCut_Xs(fl_xsid) << " && ";
    sTmp << LRNB_cut;
    sTmp << " && rm_xs==" <<  fl_mode_xs << "&& self!=1 && gm_bg==0 && gm_xs==" << mode[k]; // [each mode]
    strcpy( add_cut[k], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }

  sTmp << makeCut_Xs(fl_xsid) << " && ";
  sTmp << LRNB_cut;
  sTmp << " && ";
  sTmp << makeCut_mode_category( fl_mode_xs, 4 ); // [other modes]
  strcpy( add_cut[nmode], (char*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();

  sTmp << makeCut_Xs(fl_xsid) << " && ";
  sTmp << LRNB_cut;
  sTmp << " && ";
  sTmp << makeCut_mode_category( fl_mode_xs ); // [total]
  strcpy( add_cut[nmode+1], (char*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();

  sTmp << makeCut_Xs(fl_xsid) << " && ";
  sTmp << LRNB_cut;
  sTmp << " && ";
  sTmp << makeCut_mode_category( fl_mode_xs, 1 ); // [true]
  strcpy( add_cut[nmode+2], (char*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();

  sTmp << makeCut_Xs(fl_xsid) << " && ";
  sTmp << LRNB_cut;
  sTmp << " && self!=1 && ";
  sTmp << makeCut_mode_q2fl( fl_mode_xs,1,1 ); // false( q2=true, fl=true )
  strcpy( add_cut[nmode+3], (char*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();

  sTmp << makeCut_Xs(fl_xsid) << " && ";
  sTmp << LRNB_cut;
  sTmp << " && self!=1 && ";
  sTmp << makeCut_mode_q2fl( fl_mode_xs,1,2 ); // false( q2=true, fl=unknown )
  strcpy( add_cut[nmode+4], (char*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();

  sTmp << makeCut_Xs(fl_xsid) << " && ";
  sTmp << LRNB_cut;
  sTmp << " && self!=1 && ";
  sTmp << makeCut_mode_q2fl( fl_mode_xs,1,0 ); // false( q2=true, fl=false )
  strcpy( add_cut[nmode+5], (char*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();
  
  sTmp << makeCut_Xs(fl_xsid) << " && ";
  sTmp << LRNB_cut;
  sTmp << " && self!=1 && ";
  sTmp << makeCut_mode_q2fl( fl_mode_xs,0,1 ); // false( q2=false, fl=true )
  strcpy( add_cut[nmode+6], (char*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();

  sTmp << makeCut_Xs(fl_xsid) << " && ";
  sTmp << LRNB_cut;
  sTmp << " && self!=1 && ";
  sTmp << makeCut_mode_q2fl( fl_mode_xs,0,2 ); // false( q2=false, fl=unknown )
  strcpy( add_cut[nmode+7], (char*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();

  sTmp << makeCut_Xs(fl_xsid) << " && ";
  sTmp << LRNB_cut;
  sTmp << " && self!=1 && ";
  sTmp << makeCut_mode_q2fl( fl_mode_xs,0,0 ); // false( q2=false, fl=false )
  strcpy( add_cut[nmode+8], (char*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();

  
  using namespace Mbc;
  //using namespace deltaE;
  //using namespace M_Xs;
  //using namespace M_ll;

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
    chain[j] = new MChain( infile, tname, branch_table(), nfile[j], tail );
    nominal_cut_selection( chain[j], fl_mode_ll )( chain[j]->GetCut(), tname );
  }

  // ++++++++++++++++++++++++
  // cut change
  for( Int_t j=0; j<Nchain; j++ ){
    //chain[j]->GetCut()->Set( "dzll3d", 1,  0.015 );
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
  std::ofstream fout( Form("log/log_rec_entry_%d_lep%d_gmxs%d_set%s.dat",fl_mode_xs, fl_mode_ll, fl_xsid, setname) );
  Double_t entry_all[Nhist] = {0};
  Double_t total_event1     = 0; // sum of each mode
  Double_t total_event2     = 0; // sum of q2fl
  for( Int_t i=0; i<Nhist; i++ ){
    std::cout << Form("<hist %d> ",i);
    entry_all[i] = hist[i]->GetEntries();
    if( i!=nmode+1 && i<=nmode+2 ) total_event1 += entry_all[i];
    if( i!=nmode+1 && i>=nmode+2 ) total_event2 += entry_all[i];
    fout << std::setw(5)  << fl_mode_xs << " "
	 << std::setw(5)  << fl_mode_ll << " "
	 << std::setw(5)  << i          << " "
      	 << std::setw(15) << std::setprecision(15)
	 << entry_all[i]  << std::endl;
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
    
    hist[i]->SetTitle( sTmp.Data() );
    std::cout << sTmp.Data() << std::endl;
  }
  fout << "*" << std::endl;
  fout.close();
  if( total_event1 != entry_all[nmode+1] ) std::cerr << "mismatch # of events(sum of each mode) -> abort()" << std::endl, abort();
  if( total_event2 != entry_all[nmode+1] ) std::cerr << "mismatch # of events(sum of   q2-fl  ) -> abort()" << std::endl, abort();
  /*
  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();
  c1->cd(1);
  TH2D* waku = Waku( Nhist, hist, xlabel );
  ((TGaxis*)waku->GetXaxis())->SetMaxDigits(3);
  ((TGaxis*)waku->GetYaxis())->SetMaxDigits(3);
  waku->Draw();

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

