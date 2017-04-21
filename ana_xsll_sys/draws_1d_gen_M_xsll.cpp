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

  const Int_t   fl_appRun     = 1;
  const Bool_t  flag_save     = true; // outfile.eps and outfile.root
  const Char_t* tname         = "h12";
  const Char_t* setname       =   "A-U";
  const Int_t   Nchain        =     3;
  const Char_t* axis      [2] = {            "Xs_m",        "llg_m" };
  const Char_t* axislabel [2] = { "M_{X_{s}} [GeV]", "M_{ll} [GeV]" };
  const Int_t   fl_mode_ll[2] = {                 1,              0 };
  const Int_t   fl_par        = 1; // 0(fm), 1(mb)

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
  TChain** chain = new TChain*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ) chain[i] = new TChain( tname );
  if( fl_par==0 ){
    std::cout << "[0 : fm410(default)] " << chain[0]->Add( Form("~/ewp/ana/data/sigmc/hbk6/right/hbk_org/sigMC_*_caseB_set[%s]*.root",         setname) ) << " files" << std::endl;
    std::cout << "[1 : fm200         ] " << chain[1]->Add( Form("~/ewp/ana/data/sigmc_fm/hbk6/right/hbk_org/sigMC_*_caseB_fm200_set[%s]*.root",setname) ) << " files" << std::endl;
    std::cout << "[2 : fm480         ] " << chain[2]->Add( Form("~/ewp/ana/data/sigmc_fm/hbk6/right/hbk_org/sigMC_*_caseB_fm480_set[%s]*.root",setname) ) << " files" << std::endl;
  }else if( fl_par==1 ){
    std::cout << "[0 : mb480(default)] " << chain[0]->Add( Form("~/ewp/ana/data/sigmc/hbk6/right/hbk_org/sigMC_*_caseB_set[%s]*.root",         setname) ) << " files" << std::endl;
    std::cout << "[1 : mb465         ] " << chain[1]->Add( Form("~/ewp/ana/data/sigmc_mb/hbk6/right/hbk_org/sigMC_*_caseB_mb465_set[%s]*.root",setname) ) << " files" << std::endl;
    std::cout << "[2 : mb495         ] " << chain[2]->Add( Form("~/ewp/ana/data/sigmc_mb/hbk6/right/hbk_org/sigMC_*_caseB_mb495_set[%s]*.root",setname) ) << " files" << std::endl;
  }
  
  TH1D**** hist = new TH1D***[2]; // [ee,mm][xs,cc][fm/mb]
  for( Int_t j=0; j<2; j++ ){
    hist[j] = new TH1D**[2];
    for( Int_t k=0; k<2; k++ ){
      hist[j][k] = new TH1D*[2];
      for( Int_t i=0; i<Nchain; i++ ){
	hist[j][k][i] = new TH1D( Form("hist%d_%d_%d",j,k,i), Form("hist%d_%d_%d",j,k,i), 100, 0.0, 5.0);
	if( k==0 ) chain[i]->Project ( Form("hist%d_%d_%d",j,k,i), axis[k], Form("gm_l==%d && llg_m>0.2 && gm_fl_xs>0", fl_mode_ll[j]) );
	else       chain[i]->Project ( Form("hist%d_%d_%d",j,k,i), axis[k], Form("gm_l==%d && llg_m>0.2 && gm_fl_xs>0 && Xs_m>1.1 && Xs_m<2.0", fl_mode_ll[j]) );
	Deco( hist[j][k][i], 1, i+1, i+1 );
      }
    }
  }

  Double_t scale_factor[Nchain] = {0};
  for( Int_t j=0; j<2; j++ ){
    for( Int_t i=0; i<Nchain; i++ ){
      scale_factor[i] = hist[j][0][i]->Integral( hist[j][0][i]->FindBin(1.1+0.000001), hist[j][0][i]->FindBin(2.0-0.000001) );
    }
  }
  //std::cout << hist[0][0][0]->FindBin(1.1+0.000001) << " : " << hist[0][0][0]->FindBin(2.0-0.000001) << std::endl;
  for( Int_t j=0; j<2; j++ ){
    for( Int_t k=0; k<2; k++ ){
      for( Int_t i=0; i<Nchain; i++ ){
	hist[j][k][i]->Scale(1/scale_factor[i]);
	//if( k==0 ) hist[j][k][i]->Scale(1/scale_factor[i]);
	//else       hist[j][k][i]->Scale(1/hist[j][k][i]->Integral());
      }
    }
  }
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
  TCanvas* c1 = Canvas( "c1","c1",2,2 );
  c1->Draw();
  TH2D***   w  = new TH2D**[2];
  for( Int_t j=0; j<2; j++ ){
    w[j] = new TH2D*[2];
    for( Int_t k=0; k<2; k++ ){
      w[j][k] = Waku( Nchain, &hist[j][k][0], axislabel[k], Form("w%d_%d",j,k), Form("%s lep%d",axis[k],fl_mode_ll[j]) );
      c1->cd(2*j+k+1);
      w[j][k]->Draw();
      for( Int_t i=0; i<Nchain; i++ ) hist[j][k][i]->Draw("same");
      
    }
  }

  TLegend* leg = new TLegend( 0.70,0.75,0.99,0.99 );
  if( fl_par==0 ){
    leg->AddEntry( hist[0][0][0], "p_{f} = 0.41 GeV", "L" );
    leg->AddEntry( hist[0][0][1], "p_{f} = 0.20 GeV", "L" );
    leg->AddEntry( hist[0][0][2], "p_{f} = 0.48 GeV", "L" );
  }else if( fl_par==1 ){
    leg->AddEntry( hist[0][0][0], "M_{b} = 4.80 GeV", "L" );
    leg->AddEntry( hist[0][0][1], "M_{b} = 4.65 GeV", "L" );
    leg->AddEntry( hist[0][0][2], "M_{b} = 4.95 GeV", "L" );
  }
  leg->Draw();
  
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
  c1->Update();
  if( flag_save ){
    if( fl_par==0 ){
      c1->Print( Form("pic/M_xsll_fm_set%s.eps", setname) );
      c1->Print( Form("pic/M_xsll_fm_set%s.root",setname) );
    }else if( fl_par==1 ){
      c1->Print( Form("pic/M_xsll_mb_set%s.eps", setname) );
      c1->Print( Form("pic/M_xsll_mb_set%s.root",setname) );
    }
  }
  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();
  return 0;
}

