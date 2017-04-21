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
#include <TFile.h>
#include <TArrow.h>

Int_t main( Int_t argc, Char_t** argv ){
  using namespace rd;
  using namespace d0mass;
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==6 || argc==7) ) std::cerr << "wrong input" << std::endl
					<< " Usage : ./draws_1d (int)fl_pid (int)fl_par (int)fl_mom (int)fl_cos (int)fl_adef [(int)fl_appRun]" << std::endl
					<< "[fl_pid ] : 1(e-id), 0(mu-id), 2(no-pid)" << std::endl
					<< "[fl_par ] : 0(K), 1(pi)"                  << std::endl
					<< "[fl_adef] : 0(0.4GeV), 1(0.6GeV) "        << std::endl, abort();

  const Int_t fl_pid  = atoi( argv[1] );
  const Int_t fl_par  = atoi( argv[2] );
  const Int_t fl_mom  = atoi( argv[3] );
  const Int_t fl_cos  = atoi( argv[4] );
  const Int_t fl_adef = atoi( argv[5] );
  Int_t       fl_appRun = 1;
  if( argc==7 ) fl_appRun = atoi( argv[6] );
  if( fl_mom > n_mom[fl_adef] ) std::cerr << "[ABORT] Invalid fl_mom : " << fl_mom << std::endl, abort();
  if( fl_cos > n_cos          ) std::cerr << "[ABORT] Invalid fl_cos : " << fl_cos << std::endl, abort();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  sTmp << indir << "*"; // small sample for test
  //sTmp << indir << "*-e55r*"; // small sample for test
  //sTmp << indir << "*-e21r*"; // small sample for test
  Char_t tmp_infile[255];
  strcpy( tmp_infile, (Char_t*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();
  const Char_t* infile = (const Char_t*)tmp_infile;  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t Nchain = 1;
  const Int_t Nhist  = 1;
  const Int_t nfile  = 0;
  
  const Int_t add[Nhist][Nchain] ={
    {1}, // total
  };

  const Bool_t flag_save = true; // outfile.eps and outfile.root
  
  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ) add_cut[i] = new Char_t[4096];
  for( Int_t i=0; i<Nchain; i++ ){

    // [PID-cut]
    if     ( fl_pid==1 ) sTmp << Form( "(eprob_%d>0.8) ",                 fl_par );         //  e-id
    else if( fl_pid==0 ) sTmp << Form( "(mulmu_%d>0.97 && much2_%d>0 ) ", fl_par, fl_par ); // mu-id
    else if( fl_pid==2)  sTmp << "1";                                                       // no-pid
    else                 std::cerr << "[ABORT] invalid pid flag : " << fl_pid << std::endl, abort();
    // [Category]
    sTmp << " && 1";

    // [Momentum]
    if( fl_mom>0 ) sTmp << " && (pmag_" << fl_par << " > " << mom_min[fl_adef] << " + " << del_mom[fl_adef]*(fl_mom-1) << ")"
			<< " && (pmag_" << fl_par << " < " << mom_min[fl_adef] << " + " << del_mom[fl_adef]* fl_mom    << ")";
    // [Direction]
    if( fl_cos>0 ) sTmp << " && (cos_" << fl_par << " > " << cos_min << " + " << del_cos*(fl_cos-1) << ")"
			<< " && (cos_" << fl_par << " < " << cos_min << " + " << del_cos* fl_cos    << ")";

    strcpy( add_cut[i], (Char_t*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TH1D**    tmphist = new TH1D*  [Nchain];
  TH1D**    hist    = new TH1D*  [Nhist];
  TCanvas*  c1      = Canvas( "c1","c1", 2, 1 );
  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++
  std::cout << Form(" ************************ make tmphist ( %s, %s ) *************************************",tname,axis) << std::endl;
  std::cout << "<infile %> ";
  MChain* chain = new MChain( infile, tname, branch_table(), nfile, tail );
  nominal_cut_selection( chain )( chain->GetCut(), tname );

  // ++++++++++++++++++++++++
  // cut change
  chain->GetCut()->Set( "m_d_vf",  0 );
  //chain->GetCut()->Set( "xdst_vf", 0 );

  // ++++++++++++++++++++++++
  chain->GetCut()->Display(0);
  chain->MakeTree();
  
  std::cout << std::endl
	    << " ************************ changed cut *************************************" << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<tmphist %d > ", j ) << std::endl;
    std::cout << "add_cut : " << add_cut[j] << std::endl;
  }
  // +++++++ make tmphist ++++++++++++++++++++++++++++++++++
  for( Int_t j=0; j<Nchain; j++ ){
    tmphist[j] = new TH1D( Form("tmphist%d",j), Form("%s",chain->GetChange()), xbin,offset+xmin,offset+xmax );
    chain->GetTree()->Project( Form("tmphist%d",j), axis, add_cut[j] );
  }
  
  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make hist *************************************" << std::endl;
  for( Int_t i=0; i<Nhist; i++ ){
    hist[i] = new TH1D( Form("hist%d",i),Form("hist%d",i), xbin,offset+xmin,offset+xmax );
    Deco( hist[i], 3, col_fil[i], 1 );
    for( Int_t j=0; j<Nchain; j++ ){
      if( add[i][j] ) hist[i]->Add( tmphist[j] );
    }
  }

  // +++++++ display +++++++++++++++++++++++++++++++++++++++++++++++++
  Double_t entry_all[Nhist] = {0};
  Double_t entry_sig[Nhist] = {0}; // # of events in signal box region
  for( Int_t i=0; i<Nhist; i++ ){
    std::cout << Form("<hist %d> ",i);
    entry_all[i] = hist[i]->GetEntries();
    Double_t entry_canvas = hist[i]->Integral();
    Double_t entry_under  = hist[i]->GetBinContent(0);
    Double_t entry_over   = hist[i]->GetBinContent(xbin+1);

    for( Int_t j=hist[i]->FindBin(PDGmass::d0-0.030+0.0000001); j<=hist[i]->FindBin(PDGmass::d0+0.030-0.000000001); j++ ) entry_sig[i] += hist[i]->GetBinContent(j);
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
    sTmp += " / sig  : ";
    sTmp += entry_sig[i];
    sTmp += "]";
    
    hist[i]->SetTitle( sTmp.Data() );
    std::cout << sTmp.Data() << std::endl;
  }

  for( Int_t j=0; j<Nchain; j++ ){
    double tmp_entry = 0;
    tmp_entry = tmphist[j]->Integral( tmphist[j]->FindBin(PDGmass::d0               -0.030 +0.00000001),tmphist[j]->FindBin(PDGmass::d0               +0.030 -0.000000001) );
    std::cout << "tmphist" << std::setw(3) << std::right << j << " : " << std::setw(12) << std::right << tmp_entry << " ( " << std::setw(12) << std::right << 100 * tmp_entry/entry_sig[Nhist-1] << " % )" << std::endl;
    }
  
  if( hist[Nhist-1]->Integral()==0 ){
    std::cout << Form( "[EXIT] No events in canvas : pid%d, par%d, mom%d, cos%d", fl_pid, fl_par, fl_mom, fl_cos ) << std::endl;
    return 0;
  }

  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();

  c1->cd(1);
  TH2D* waku = Waku( Nhist, hist, xlabel );
  ((TGaxis*)waku->GetXaxis())->SetMaxDigits(3);
  ((TGaxis*)waku->GetYaxis())->SetMaxDigits(3);
  waku->SetTitle( Form("%s", xlabel) );
  waku->Draw();

  // ++++++++ decoration ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TArrow* ar_d = new TArrow( PDGmass::d0-0.030, 1.02*((TAxis*)waku->GetYaxis())->GetXmax(),
			     PDGmass::d0+0.030, 1.02*((TAxis*)waku->GetYaxis())->GetXmax(),
			     0.01,"<|>" );
  ar_d->SetLineColor(2);
  ar_d->SetFillColor(2);
  ar_d->SetLineWidth(2);
  ar_d->Draw();

  // +++++++ tlegend1 ++++++++++++++++++++++++++++++++++
  TLegend* legend1 = new TLegend( 0.75,0.75,0.99,0.99 );
  legend1->AddEntry( hist[0],"Total", "F" );

  
  for(Int_t i=Nhist-1; i>=0; i-- ) hist[i]->Draw( "same" );

  // +++++++ draw ++++++++++++++++++++++++++++++++++

  c1->cd(2);
  gPad->SetLogy();
  ar_d->Draw();
  for(Int_t i=Nhist-1; i>=0; i-- ){
    if( i==Nhist-1 ) hist[i]->Draw();
    else             hist[i]->Draw( "same" );
  }

  c1->cd(2);
  legend1->Draw();
  
  c1->Update();
  if( flag_save ){
    c1->Print    ( Form( "pic/d0_plot_data_pid%d_par%d_mom%d_cos%d_adef%d.eps",  fl_pid, fl_par, fl_mom, fl_cos, fl_adef) );
    TFile outfile( Form( "pic/d0_plot_data_pid%d_par%d_mom%d_cos%d_adef%d.root", fl_pid, fl_par, fl_mom, fl_cos, fl_adef), "RECREATE" );
    c1->Write();
    for( Int_t i=0; i<Nhist; i++ ) hist[i]->Write();
    outfile.Close();
  }
  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();

  delete   chain;
  delete[] tmphist;
  delete[] hist;
  delete[] add_cut;
  delete   c1;
  delete   legend1;
  delete   ar_d;
    
  return 0;
}

