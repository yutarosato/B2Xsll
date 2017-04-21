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
  if( !(argc==2 || argc==3) ) std::cerr << "wrong input" << std::endl
					<< " Usage : ./draws_1d (int)fl_axis [(int)fl_appRun]" << std::endl, abort();

  const Int_t fl_axis   = atoi( argv[1] );
  Int_t       fl_appRun = 1;
  if( argc==3 ) fl_appRun = atoi( argv[2] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  sTmp << indir << "*-e55r*"; // small sample for test
  Char_t tmp_infile[255];
  strcpy( tmp_infile, (Char_t*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();
  const Char_t* infile = (const Char_t*)tmp_infile;  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t Nchain        = 10;
  const Int_t Nhist         = 5;
  const Int_t nfile[Nchain] = {0};
  
  const Int_t add[Nhist][Nchain] ={
    {1,1,1,1},             // fake-D
    {1,1,1,1,1,1},         // fake-D(switch K-pi)
    {1,1,1,1,1,1,1,1},     // true-D + fake-slow-pion
    {1,1,1,1,1,1,1,1,1},   // true-D*
    {0,0,0,0,0,0,0,0,0,1}, // total
  };

  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ) add_cut[i] = new Char_t[4096];
  
  add_cut[0] = "( !(mc_0==1 && mc_1==1 && idmo_0==idmo_1 && abs(idmo_0)==421) && !(mc_0==0 && mc_1==0 && idmo_0==idmo_1 && abs(idmo_0)==421) && mc_0==0 && mc_1==0 )";                                   // fake-D (fake-K & fake-pi)
  add_cut[1] = "( !(mc_0==1 && mc_1==1 && idmo_0==idmo_1 && abs(idmo_0)==421) && !(mc_0==0 && mc_1==0 && idmo_0==idmo_1 && abs(idmo_0)==421) && mc_0==1 && mc_1==0 )";                                   // fake-D (true-K & fake-pi)
  add_cut[2] = "( !(mc_0==1 && mc_1==1 && idmo_0==idmo_1 && abs(idmo_0)==421) && !(mc_0==0 && mc_1==0 && idmo_0==idmo_1 && abs(idmo_0)==421) && mc_0==0 && mc_1==1 )";                                   // fake-D (fake-K & true-pi)
  add_cut[3] = "( !(mc_0==1 && mc_1==1 && idmo_0==idmo_1 && abs(idmo_0)==421) && !(mc_0==0 && mc_1==0 && idmo_0==idmo_1 && abs(idmo_0)==421) && mc_0==1 && mc_1==1 )";                                   // fake-D (true-K & true-pi)
  add_cut[4] = "( !(mc_0==1 && mc_1==1 && idmo_0==idmo_1 && abs(idmo_0)==421) &&  (mc_0==0 && mc_1==0 && idmo_0==idmo_1 && abs(idmo_0)==421) && !(mc_2==1 && abs(idmo_2)==413) )"; // fake-D(switch K-pi) + fake-slow-pion
  add_cut[5] = "( !(mc_0==1 && mc_1==1 && idmo_0==idmo_1 && abs(idmo_0)==421) &&  (mc_0==0 && mc_1==0 && idmo_0==idmo_1 && abs(idmo_0)==421) &&  (mc_2==1 && abs(idmo_2)==413) )"; // fake-D(switch K-pi) + true-slow-pion
  add_cut[6] = "(  (mc_0==1 && mc_1==1 && idmo_0==idmo_1 && abs(idmo_0)==421) &&  (mc_2==0) )";                     // true-D, fake-slow-pion(not pion)
  add_cut[7] = "(  (mc_0==1 && mc_1==1 && idmo_0==idmo_1 && abs(idmo_0)==421) &&  (mc_2==1 && abs(idmo_2)!=413) )"; // true-D, fake-slow-pion
  add_cut[8] = "(  (mc_0==1 && mc_1==1 && idmo_0==idmo_1 && abs(idmo_0)==421) &&  (mc_2==1 && abs(idmo_2)==413) )"; // true-D*
  add_cut[9] = "1";


  using namespace var;
  
  const Bool_t flag_save = true; // outfile.eps and outfile.root

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain   = new MChain*[Nchain];
  TH1D**    tmphist = new TH1D*  [Nchain];
  TH1D**    hist    = new TH1D*  [Nhist];
  TF1**     func    = new TF1*   [Nhist];
  TCanvas*  c1      = Canvas( "c1","c1",2 );
  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++
  std::cout << Form(" ************************ make tmphist ( %s, %s ) *************************************",tname,axis[fl_axis]) << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<infile %d > ", j );
    chain[j] = new MChain( infile, tname, branch_table(), nfile[j], tail );
    nominal_cut_selection( chain[j] )( chain[j]->GetCut(), tname );
  }

  // ++++++++++++++++++++++++
  // cut change
  for( Int_t j=0; j<Nchain; j++ ){
    if( strcmp(axis[fl_axis],"m_d_vf" )==0 ) chain[j]->GetCut()->Set( "m_d_vf", 0 );
    if( strcmp(axis[fl_axis], "dm_vf" )==0 ) chain[j]->GetCut()->Set(  "dm_vf", 0 );
    if( strcmp(axis[fl_axis],"m_drev" )==0 ) chain[j]->GetCut()->Set( "m_drev", 0 );
  }
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
  for( Int_t i=0; i<Nhist; i++ ){
    hist[i] = new TH1D( Form("hist%d",i),Form("hist%d",i), xbin[fl_axis],offset[fl_axis]+xmin[fl_axis],offset[fl_axis]+xmax[fl_axis] );
    Deco( hist[i], 3, col_fil[i], col_fil[i] );
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
    Double_t entry_over   = hist[i]->GetBinContent(xbin[fl_axis]+1);

    if( strcmp(axis[fl_axis],"m_d_vf")==0 ) {
      for( Int_t j=hist[i]->FindBin(PDGmass::d0-0.030+0.0000001); j<=hist[i]->FindBin(PDGmass::d0+0.030-0.000000001); j++ ) entry_sig[i] += hist[i]->GetBinContent(j);
    }else if( strcmp(axis[fl_axis],"dm_vf")==0 ) {
      for( Int_t j=hist[i]->FindBin(PDGmass::dstrp-PDGmass::d0-0.0015+0.0000001); j<=hist[i]->FindBin(PDGmass::dstrp-PDGmass::d0+0.0015-0.0000001); j++ ) entry_sig[i] += hist[i]->GetBinContent(j);
    }
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
    if( strcmp(axis[fl_axis],"m_d_vf")==0 ) tmp_entry = tmphist[j]->Integral( tmphist[j]->FindBin(PDGmass::d0               -0.030 +0.00000001),tmphist[j]->FindBin(PDGmass::d0               +0.030 -0.000000001) );
    if( strcmp(axis[fl_axis],"dm_vf" )==0 ) tmp_entry = tmphist[j]->Integral( tmphist[j]->FindBin(PDGmass::dstrp-PDGmass::d0-0.0015+0.00000001),tmphist[j]->FindBin(PDGmass::dstrp-PDGmass::d0+0.0015-0.000000001) );
    std::cout << "tmphist" << std::setw(3) << std::right << j << " : " << std::setw(12) << std::right << tmp_entry << " ( " << std::setw(12) << std::right << 100 * tmp_entry/entry_sig[Nhist-1] << " % )" << std::endl;
    }
  
  if( entry_sig[Nhist-1] != entry_sig[Nhist-2] ) std::cerr << "[ABORT] : wrong entry : " << entry_sig[Nhist-1] << " ?=? " << entry_sig[Nhist-2] << std::endl, abort();

  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();

  c1->cd(1);
  TH2D* waku = Waku( Nhist, hist, xlabel[fl_axis] );
  ((TGaxis*)waku->GetXaxis())->SetMaxDigits(3);
  ((TGaxis*)waku->GetYaxis())->SetMaxDigits(3);
  waku->SetTitle( Form("%s", xlabel[fl_axis]) );
  waku->Draw();

  // ++++++++ decoration ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TArrow* ar_d = new TArrow( PDGmass::d0-0.030, 1.02*((TAxis*)waku->GetYaxis())->GetXmax(),
			     PDGmass::d0+0.030, 1.02*((TAxis*)waku->GetYaxis())->GetXmax(),
			     0.01,"<|>" );
  ar_d->SetLineColor(2);
  ar_d->SetFillColor(2);
  ar_d->SetLineWidth(2);
  if( strcmp(axis[fl_axis],"m_d_vf")==0 ) ar_d->Draw();
  if( strcmp(axis[fl_axis],"m_drev")==0 ) ar_d->Draw();

  TArrow* ar_deltam = new TArrow( PDGmass::dstrp-PDGmass::d0-0.0015, 1.02*((TAxis*)waku->GetYaxis())->GetXmax(),
				  PDGmass::dstrp-PDGmass::d0+0.0015, 1.02*((TAxis*)waku->GetYaxis())->GetXmax(),
				  0.01,"<|>" );
  ar_deltam->SetLineColor(2);
  ar_deltam->SetFillColor(2);
  ar_deltam->SetLineWidth(2);
  if( strcmp(axis[fl_axis],"dm_vf")==0 ) ar_deltam->Draw();
  
  for(Int_t i=Nhist-2; i>=0; i-- ) hist[i]->Draw( "same" );

  // +++++++ tlegend1 ++++++++++++++++++++++++++++++++++
  TLegend* legend1 = new TLegend( 0.75,0.75,0.99,0.99 );
  legend1->AddEntry( hist[3],"true-D*",                "P" );
  legend1->AddEntry( hist[2],"true-D + fake-slow-#pi", "P" );
  legend1->AddEntry( hist[1],"false-D(switch K-#pi)",  "P" );
  legend1->AddEntry( hist[0],"false-D",                "P" );
  legend1->Draw();

  // +++++++ draw ++++++++++++++++++++++++++++++++++

  c1->cd(2);
  gPad->SetLogy();
  if( strcmp(axis[fl_axis],"m_d_vf")==0 ) ar_d->Draw();
  if( strcmp(axis[fl_axis],"m_drev")==0 ) ar_d->Draw();
  if( strcmp(axis[fl_axis],"dm_vf")==0 ) ar_deltam->Draw();  
  for(Int_t i=Nhist-2; i>=0; i-- ){
    if( i==Nhist-2 ) hist[i]->Draw();
    else             hist[i]->Draw( "same" );
  }
  
  c1->Update();
  if( flag_save ){
    c1->Print( Form("pic/var_%s.eps", fname[fl_axis]) );
    c1->Print( Form("pic/var_%s.root",fname[fl_axis]) );
  }
  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();

  delete[] chain;
  delete[] tmphist;
  delete[] hist;
  delete[] func;
  delete[] add_cut;
  delete   c1;
  delete   legend1;
  delete ar_d;
  delete ar_deltam;
    
  return 0;
}

