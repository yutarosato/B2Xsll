#include <iostream>
#include <math.h>

#include "../Util/Style.h"
#include "../Util/Canvas.h"
#include "../Util/const.h"
#include "../Util/MChain.h"
#include "../Util/MCut.h"
#include "../Util/MCut_array.h"
#include "../Util/FitFunc.h"
#include "../Set_fake_lepton/Nominal_cut_selection.h"
#include "../Set_fake_lepton/Branch.h"
#include "../Set_fake_lepton/makeCut.h"

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
#include <TLeaf.h>

Int_t main( Int_t argc, Char_t** argv ){
  using namespace gmc_fl;
  TApplication app( "app", &argc, argv );
  Style();

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  const Int_t Nchain = 1;
  std::stringstream sTmp;
  Char_t** infile = new Char_t*[Nchain];
  infile[0] = new Char_t[1024];
  sTmp << "~/ewp/ana/data/gmc_fl/hbk5/right/hbk_single_org/mixed/9999/gMC_*_e03*";
  strcpy( infile[0], (char*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  const Int_t xbin = 7;

  const Double_t xbins_ee[xbin+1] = {
    -1.0, //cos(2*TMath::Pi()*151/360), // -0.875
    cos(2*TMath::Pi()*132/360), // -0.669
    cos(2*TMath::Pi()*125/360), // -0.574
    cos(2*TMath::Pi()*60/360),  // +0.50
    cos(2*TMath::Pi()*40/360),  // +0.766
    cos(2*TMath::Pi()*35/360),  // +0.819
    cos(2*TMath::Pi()*25/360),  // +0.906
    1.0, //cos(2*TMath::Pi()*18/360), // +0.951
  };

  const Double_t xbins_mm[xbin+1] = {
    -1.0, //cos(2*TMath::Pi()*150/360), // -0.866
    cos(2*TMath::Pi()*145/360), // -0.819
    cos(2*TMath::Pi()*130/360), // -0.643
    cos(2*TMath::Pi()*117/360), // -0.454
    cos(2*TMath::Pi()*51/360), // 0.629
    cos(2*TMath::Pi()*37/360), // 0.799
    cos(2*TMath::Pi()*25/360), // 0.906
    1.0, //cos(2*TMath::Pi()*17/360), // 0.956
  };
  
  const Int_t    ybin = 10;
  const Double_t ymin = 0.0;
  const Double_t ymax = 5.0;

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  TChain* c = new TChain("h511");
  std::cout << c->Add( infile[0] ) << " files" << std::endl;
  TH2D** h_ee_eff = new TH2D*[3];
  h_ee_eff[0] = new TH2D( "eff_e_pre", "eff_e_pre", xbin, xbins_ee, ybin, ymin, ymax );
  h_ee_eff[1] = new TH2D( "eff_e_pos", "eff_e_pos", xbin, xbins_ee, ybin, ymin, ymax );
  h_ee_eff[2] = new TH2D( "eff_e",     "eff_e",     xbin, xbins_ee, ybin, ymin, ymax );
  h_ee_eff[2]->SetAxisRange( 0.0, 1.0, "Z" );
  TH2D** h_ee_fake = new TH2D*[3];
  h_ee_fake[0] = new TH2D( "fake_e_pre", "fake_e_pre", xbin, xbins_ee, ybin, ymin, ymax );
  h_ee_fake[1] = new TH2D( "fake_e_pos", "fake_e_pos", xbin, xbins_ee, ybin, ymin, ymax );
  h_ee_fake[2] = new TH2D( "fake_e",     "fake_e",     xbin, xbins_ee, ybin, ymin, ymax );
  h_ee_fake[2]->SetAxisRange( 0.0, 1.0, "Z" );
  TH2D** h_mm_eff = new TH2D*[3];
  h_mm_eff[0] = new TH2D( "eff_mu_pre", "eff_mu_pre", xbin, xbins_mm, ybin, ymin, ymax );
  h_mm_eff[1] = new TH2D( "eff_mu_pos", "eff_mu_pos", xbin, xbins_mm, ybin, ymin, ymax );
  h_mm_eff[2] = new TH2D( "eff_mu",     "eff_mu",     xbin, xbins_mm, ybin, ymin, ymax );
  h_mm_eff[2]->SetAxisRange( 0.0, 1.0, "Z" );
  TH2D** h_mm_fake = new TH2D*[3];
  h_mm_fake[0] = new TH2D( "fake_mu_pre", "fake_mu_pre", xbin, xbins_mm, ybin, ymin, ymax );
  h_mm_fake[1] = new TH2D( "fake_mu_pos", "fake_mu_pos", xbin, xbins_mm, ybin, ymin, ymax );
  h_mm_fake[2] = new TH2D( "fake_mu",     "fake_mu",     xbin, xbins_mm, ybin, ymin, ymax );
  h_mm_fake[2]->SetAxisRange( 0.0, 1.0, "Z" );

  ///* LEPTON TRAKCKS
  c->Project("eff_e_pre",   "pi1p:pi1c", "abs(pi1selfi)==11");
  c->Project("eff_e_pos",   "pi1p:pi1c", "abs(pi1selfi)==11&&pi1_eid>0.80");
  c->Project("fake_e_pre",  "pi1p:pi1c", "abs(pi1selfi)==11");
  c->Project("fake_e_pos",  "pi1p:pi1c", "abs(pi1selfi)==11&&pi1_eid<0.80 && pi1_kid<0.6");
  c->Project("eff_mu_pre",  "pi1p:pi1c", "abs(pi1selfi)==13");
  c->Project("eff_mu_pos",  "pi1p:pi1c", "abs(pi1selfi)==13&&pi1_muid>0.97");
  c->Project("fake_mu_pre", "pi1p:pi1c", "abs(pi1selfi)==13");
  c->Project("fake_mu_pos", "pi1p:pi1c", "abs(pi1selfi)==13&&pi1_muid<0.97&& pi1_kid<0.6");
  //*/

  /* PION TRACKS
  c->Project("eff_e_pre",   "pi1p:pi1c", "rm_l==1&&pi1self==1");
  c->Project("eff_e_pos",   "pi1p:pi1c", "rm_l==1&&pi1self==1&&pi1_kid<0.60");
  c->Project("fake_e_pre",  "pi1p:pi1c", "rm_l==1&&pi1self==1");
  c->Project("fake_e_pos",  "pi1p:pi1c", "rm_l==1&&pi1self==1&&pi1_kid>0.60");
  c->Project("eff_mu_pre",  "pi1p:pi1c", "rm_l==0&&pi1self==1");
  c->Project("eff_mu_pos",  "pi1p:pi1c", "rm_l==0&&pi1self==1&&pi1_kid<0.60");
  c->Project("fake_mu_pre", "pi1p:pi1c", "rm_l==0&&pi1self==1");
  c->Project("fake_mu_pos", "pi1p:pi1c", "rm_l==0&&pi1self==1&&pi1_kid>0.60");
  */

  for( Int_t tx=1; tx<=xbin; tx++ ){
    for( Int_t ty=1; ty<=ybin; ty++ ){
      if( h_ee_eff[0] ->GetBinContent(tx,ty) ) h_ee_eff[2] ->SetBinContent(tx,ty, h_ee_eff[1] ->GetBinContent(tx,ty)/h_ee_eff[0] ->GetBinContent(tx,ty) );
      if( h_ee_fake[0]->GetBinContent(tx,ty) ) h_ee_fake[2]->SetBinContent(tx,ty, h_ee_fake[1]->GetBinContent(tx,ty)/h_ee_fake[0]->GetBinContent(tx,ty) );
      if( h_mm_eff[0] ->GetBinContent(tx,ty) ) h_mm_eff[2] ->SetBinContent(tx,ty, h_mm_eff[1] ->GetBinContent(tx,ty)/h_mm_eff[0] ->GetBinContent(tx,ty) );
      if( h_mm_fake[0]->GetBinContent(tx,ty) ) h_mm_fake[2]->SetBinContent(tx,ty, h_mm_fake[1]->GetBinContent(tx,ty)/h_mm_fake[0]->GetBinContent(tx,ty) );
    }
  }

  TCanvas* c1 = Canvas( "c1","c1",2, 2 );
  c1->cd(1); h_ee_eff[2] ->Draw("COLZ");
  c1->cd(2); h_ee_fake[2]->Draw("COLZ");
  c1->cd(3); h_mm_eff[2] ->Draw("COLZ");
  c1->cd(4); h_mm_fake[2]->Draw("COLZ");
  
  std::cout << "finish" << std::endl;
  app.Run();
    
  return 0;
}
