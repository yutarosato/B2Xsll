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
#include <TFile.h>

void manip_func( TF1* func ){
  func->FixParameter( 4, 5.289 );
}

using namespace Mbc;
const Bool_t flag_save     = true; // outfile.eps and outfile.root
const Bool_t flag_fit      = true;
const Bool_t flag_file     = !true; // read hist from file
const Int_t  sel_fun       = 15; // 10(gauss), 15(gauss+argus)
const Bool_t flag_k4pi     = !true; // 1(veto  K4pi    modes)
const Bool_t flag_unflavor = !true; // 1(veto unflavor modes)
const Bool_t flag_ccpi0    = true; // 1(veto ccpi0 peak    )

Int_t main( Int_t argc, Char_t** argv ){
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //const Char_t* infile    = "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut3/cc_522/hbk_orgksfw_vtxcl_fmiss1_bb_bcs_522_lrnb_before_combine/CC_mixedjpsi_f???_*_"; // cc-MC (small sample)
  const Char_t* infile    = "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut3/cc_522/CC_*_"; // cc-MC
  //const Char_t* infile    = "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut3/bkg_522/gMC_*_"; // gMC
  //const Char_t* infile    = "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut3/rd_522/RD_*_"; // RD
  Int_t         fl_appRun = 1;
  fl_appRun=0;
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t Nchain        = 2; // ee, mm
  const Int_t Nhist         = 3; // ee, mm, ee+mm
  const Int_t nfile[Nchain] = {0};

  const Int_t fl_mode_ll[Nchain] = {1,0}; // 1(e), 0(mu)
  const Int_t add[Nhist][Nchain] ={
    {1},   // ee
    {0,1}, // mm
    {1,1}, // ee+mm
  };

  std::stringstream sTmp;
  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    add_cut[i] = new Char_t[4096];
    if     ( i==0 ) sTmp << "3.096916-0.40<cc_m && 3.096916+0.15>cc_m && 3.096916-0.40<cc_morg && 3.096916+0.15>cc_morg"; // ee for J/psi samples
    else if( i==1 ) sTmp << "3.096916-0.25<cc_m && 3.096916+0.10>cc_m && 3.096916-0.25<cc_morg && 3.096916+0.10>cc_morg"; // mm for J/psi samples
    if( flag_k4pi     ) sTmp << " && " << makeCut_5body_veto().c_str();
    if( flag_unflavor ) sTmp << " && " << makeCut_unflavor_veto().c_str();
    if( flag_ccpi0    ) sTmp << " && " << cut_ccpi0;
    strcpy( add_cut[i], (Char_t*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain   = new MChain*[Nchain];
  TH1D**    tmphist = new TH1D*  [Nchain];
  TH1D**    hist    = new TH1D*  [Nhist];
  TF1**     func    = new TF1*   [Nhist];
  TCanvas*  c1      = Canvas( "c1","c1",Nhist, 1 );
  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++
  TFile file_hist("pic/Mbc_cc_fit_hz_all.root");
  if( !flag_file ) file_hist.Close();
  if( flag_file ){
    if( file_hist.IsZombie() ) std::cerr << "[ABORT] can not find input-file : " << file_hist.GetName() << std::endl, abort();
    std::cout << Form(" ************************ read hist from file ( %s ) *************************************", file_hist.GetName()) << std::endl;
    for( Int_t i=0; i<Nhist; i++ ){
      hist[i] = (TH1D*)file_hist.Get( Form("hist%d",i) );
      Deco( hist[i], 1, 1, 1 );
    }
    //file_hist.Close();
  }else{
    std::cout << Form(" ************************ make chain-tree ( %s, %s ) *************************************",tname,axis) << std::endl;
    for( Int_t j=0; j<Nchain; j++ ){
      std::cout << Form( "<infile %d > ", j );
      chain[j] = new MChain( infile, tname, branch_table(), nfile[j], "*.root" );
      nominal_cut_selection( chain[j], fl_mode_ll[j] )( chain[j]->GetCut(), tname );
    }
    
    // ++++++++++++++++++++++++
    // cut change
    for( Int_t j=0; j<Nchain; j++ ){
      chain[j]->GetCut()->Set( "Mbc", 1,  5.22, 0.0, 5.30 );
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
      
      Deco( hist[i], 1, 1, 1 );
      for( Int_t j=0; j<Nchain; j++ ){
	if( add[i][j] ) hist[i]->Add( tmphist[j] );
      }
    }
  }  
  // +++++++ display ++++++++++++++++++++++++++++++++++
  Double_t entry_sig[Nhist] = {0}; // # of events in signal box region
  for( Int_t i=0; i<Nhist; i++ ){
    std::cout << Form("<hist %d> ",i);
    Double_t entry_all    = hist[i]->GetEntries();
    Double_t entry_canvas = hist[i]->Integral();
    Double_t entry_under  = hist[i]->GetBinContent(0);
    Double_t entry_over   = hist[i]->GetBinContent(xbin+1);
    for( Int_t j=hist[i]->FindBin(5.27+0.000000001); j<=xbin; j++ ) entry_sig[i] += hist[i]->GetBinContent(j);
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

  // +++++++ make fit-func ++++++++++++++++++++++++++++++++++
  if( flag_fit ){
    for( Int_t i=0; i<Nhist; i++ ){
      func[i] = new TF1( Form("func%d",i), make_func(sel_fun),offset+xmin_fit,offset+xmax_fit,n_fitfunc_par(sel_fun) );
      if( sel_fun == 15 ){
	Double_t a     = ( hist[i]->GetBinContent(xbin)-hist[i]->GetBinContent(1))/(hist[i]->GetBinCenter(xbin)-hist[i]->GetBinCenter(1) );
	Double_t b     = hist[i]->GetBinContent(1)-a*hist[i]->GetBinCenter(1);
	Double_t mu    = hist[i]->GetBinCenter( hist[i]->GetMaximumBin() );
	Double_t sigma = 3*hist[i]->GetBinWidth(1);
	Double_t area  = ( hist[i]->GetMaximum() - (a*mu+b) )*sqrt( TMath::TwoPi() )*sigma;
	Double_t norm  = 5*hist[i]->GetBinContent(1);
	func[i]->SetParNames  ( "area","mean","sigma","norm","Ebeam","shape" );
	std::cout << area << ", " << mu << ", " << sigma << ", " << norm << std::endl
		  << hist[i]->GetMaximum() << ", " << a << ", " << b << std::endl;
	
	func[i]->SetParameters( area,   mu,    sigma,  norm,   5.289,   -60 );
	func[i]->FixParameter( 4, 5.289 );
      }else func_set_parameters(sel_fun, func[i], hist[i], xbin, offset+xmin, offset+xmax);
      func[i]->SetLineColor(2);
    }
  }
  

  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();
  TH2D** waku = new TH2D*[Nhist];
  waku[0] = Waku( Nhist, hist, Form("%s (ee)",       xlabel) );
  waku[1] = Waku( Nhist, hist, Form("%s (#mu#mu)",   xlabel) );
  waku[2] = Waku( Nhist, hist, Form("%s (ee+#mu#mu)",xlabel) );
    
  for(Int_t i=0; i<Nhist; i++ ){
    c1->cd(i+1);
    if     ( i==0 ) hist[i]->SetTitle( Form("ee"       ) );
    else if( i==1 ) hist[i]->SetTitle( Form("#mu#mu"   ) );
    else if( i==2 ) hist[i]->SetTitle( Form("ee+#mu#mu") );
    hist[i]->SetXTitle(xlabel);
    hist[i]->SetYTitle(waku[i]->GetYaxis()->GetTitle());
    hist[i]->GetXaxis()->CenterTitle();
    hist[i]->GetYaxis()->CenterTitle();
    //waku[i]->Draw();
    if( flag_fit ){
      std::cout << Form( "================================== FUNC%d =================================", i ) << std::endl;
      Double_t init_var[n_fitfunc_par(sel_fun)];
      for( Int_t m=0; m<n_fitfunc_par(sel_fun); m++ ) init_var[m] = func[i]->GetParameter(m);
      iterative_fit( hist[i], func[i], n_fitfunc_par(sel_fun), init_var, "PE0", manip_func );
      //hist[i]->Fit(func[i],"R","PE0");
      //std::cout << "chi2/NDF = " << func[i]->GetChisquare() << " / " << func[i]->GetNDF()
      //<< " = " << func[i]->GetChisquare()/func[i]->GetNDF()
      //<< std::endl;
    }else hist[i]->Draw( "same" );
  }

  c1->Update();
  if( flag_save ){
    c1->Print( Form("pic/%s_cc_fit_c1.eps", axis) );
    TFile outfile( Form("pic/%s_cc_fit.root",axis), "RECREATE" );
    c1->Write();
    for( Int_t i=0; i<Nhist; i++ ) hist[i]->Write();
    outfile.Close();
  }
  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();
  return 0;
}

