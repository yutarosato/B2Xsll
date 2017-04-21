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

void manip_func( TF1* func ){
  func->FixParameter( 2, 0.0     ); // mean
  func->SetParLimits( 1, 0.0, 1.0); // ratio(area)
  func->SetParLimits( 4, 0.0, 1.0); // ratio(sigma)
}

Int_t main( Int_t argc, Char_t** argv ){
  //using namespace sig_gmc_rd_calib;
  //using namespace sig_gmc_rd_cut2_beforebgsup;
  using namespace sig_gmc_rd_cut3_beforebgsup;
  //using namespace sig_gmc_rd_emu_beforebgsup;
  //using namespace sig_gmc_rd_emu_beforebgsup_novtxinfo;
  using namespace dzll3dorg;
  //using namespace dzll3d;
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==3 || argc==4) ) std::cerr << "wrong input" << std::endl
    					<< " Usage : ./draws_1d_var (char*)stream (double)used_nstream [(int)fl_appRun]" << std::endl
					<< std::endl, abort();
  Char_t*  stream        = argv[1];
  Double_t used_nstream  = atof(argv[2]); 
  Int_t   fl_appRun      = 1;
  if( argc==4 ) fl_appRun = atoi( argv[3] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Bool_t flag_save  = true; // outfile.eps and outfile.root
  const Bool_t flag_fit   = true;
  const Int_t  sel_fun    = 20;
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t    Nchain             = 2*2; // [gmc,rd] x [ee,mm]
  const Int_t    Nhist              = 2*3; // [gmc,rd] x [ee,mm,ee+mm]
  const Int_t    nfile[Nchain]      = {0};
  const Int_t    fl_sb[Nchain]      = {0,2,0,2}; // 0(bkg), 1(sig), 2(rd)
  //const Int_t    fl_mode_ll[Nchain] = {31,31,  // 1(ee)
  //30,30}; // 0(mm)
  const Int_t    fl_mode_ll[Nchain] = {1,1,  // 1(ee)
				       0,0}; // 0(mm)
  const Double_t scale_event_bkg = used_nstream; // gmc : N -> N/alpha
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Char_t* LRNB_cut = new Char_t[4096];
  LRNB_cut = "1";
  std::cout << "LRNB cut : " << LRNB_cut << std::endl;
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  Char_t** infile = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    infile[i] = new Char_t[1024];
    if     ( fl_sb[i]==0 ) sTmp << indir[fl_sb[i]] << "gMC_*s0[" << stream << "]_";     // bkg
    else if( fl_sb[i]==2 ) sTmp << indir[fl_sb[i]] << "RD_";                            // rd
    else                   std::cerr << "[ABORT] Wrong fl_sb : " << fl_sb[i] << std::endl, abort();
    strcpy( infile[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t add[Nhist][Nchain] ={
    {1},       // gmc, ee
    {0,1},     // rd,  ee
    {0,0,1},   // gmc, mm
    {0,0,0,1}, // rd,  mm
    {1,0,1},   // gmc, ee+mm
    {0,1,0,1}, // rd,  ee+mm
  };

  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    add_cut[i] = new Char_t[4096];
    sTmp << LRNB_cut;
    ///*
    sTmp << "&& ( " // Kll
	 << Form( " (rm_xs==  1) || " ) // K+
	 << Form( " (rm_xs== 10) || " ) // Ks
	 << Form( " (rm_xs==101 && abs(xs_m-%f)<0.050) || ", PDGmass::kstr0 ) // K+pi-
	 << Form( " (rm_xs==110 && abs(xs_m-%f)<0.050) ",    PDGmass::kstrp ) // Kspi-
	 << " )";
    //*/
    //sTmp << " && kfbclcalib>1.0e-18";
    if     ( fl_sb[i]==0 ) " && kfbclcalib>1.0e-18"; // mc
    else if( fl_sb[i]==2 ) " && kfbcl>1.0e-18"; // rd
    strcpy( add_cut[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain   = new MChain*[Nchain];
  TH1D**    tmphist = new TH1D*  [Nchain];
  TH1D**    hist    = new TH1D*  [Nhist];
  TF1**     func    = new TF1*   [Nhist];
  TCanvas*  c1      = Canvas( "c1","c1",2, 3 );
  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++
  std::cout << Form(" ************************ make tmphist ( %s, %s ) *************************************",tname,axis) << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<infile %d > ", j );
    chain[j] = new MChain( infile[j], tname, branch_table(), nfile[j], tail );
    nominal_cut_selection( chain[j], fl_mode_ll[j] )( chain[j]->GetCut(), tname );
  }

  // ++++++++++++++++++++++++
  // cut change
  for( Int_t j=0; j<Nchain; j++ ){
    chain[j]->GetCut()->Set( "dzll3dcalib", 0 );
    //chain[j]->GetCut()->Set( "Mbc", 1,  5.20, 0.0, 5.26 );
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
    if     ( fl_sb[j]==0 ) chain[j]->GetTree()->Project( Form("tmphist%d",j), axis,   add_cut[j] ); // gmc
    else if( fl_sb[j]==2 ) chain[j]->GetTree()->Project( Form("tmphist%d",j), axisrd, add_cut[j] ); // rd
    //chain[j]->GetTree()->Project( Form("tmphist%d",j), axis, add_cut[j] );
    tmphist[j]->Sumw2();
    //if( fl_sb[j]==0 ) tmphist[j]->Scale( 1/scale_event_bkg );
  }
  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make hist *************************************" << std::endl;
  for( Int_t i=0; i<Nhist; i++ ){
    hist[i] = new TH1D( Form("hist%d",i),Form("hist%d",i), xbin,offset+xmin,offset+xmax );
    //Deco( hist[i], 2, i%2+1, i%2+1 );
    Deco( hist[i], 2, 1, 1 );
    for( Int_t j=0; j<Nchain; j++ ){
      if( add[i][j] ) hist[i]->Add( tmphist[j] );
    }
    hist[i]->Scale( 1/hist[i]->GetEntries() );
    hist[i]->SetXTitle( xlabel );
    hist[i]->GetXaxis()->CenterTitle();
  }

  hist[0]->SetTitle( "ee (gmc)"        );
  hist[1]->SetTitle( "ee (rd)"         );
  hist[2]->SetTitle( "#mu#mu (gmc)"    );
  hist[3]->SetTitle( "#mu#mu (rd)"     );
  hist[4]->SetTitle( "ee+#mu#mu (gmc)" );
  hist[5]->SetTitle( "ee+#mu#mu (rd)"  );

  // +++++++ display ++++++++++++++++++++++++++++++++++
  std::cout << "=================================" << std::endl
	    << "[Scale Factor]"                    << std::endl
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

    for( Int_t j=0; j<Nchain; j++ ) if( add[i][j] ) std::cout << j << ",";
    std::cout <<  ")[ "
	      << std::setw(12) << std::right << entry_all
	      <<  " events ( canvas : "
	      << std::setw(12) << std::right << std::setprecision(3) << entry_canvas
	      << " / under : "
	      << std::setw(12) << std::right << std::setprecision(3) << entry_under
	      << " / over  : "
	      << std::setw(12) << std::right << std::setprecision(3) << entry_over
	      << " / sig  : "
	      << std::setw(12) << std::right << std::setprecision(3) << entry_sig[i]
	      << "]"
	      << std::endl;
  }
  // +++++++ make fit-func ++++++++++++++++++++++++++++++++++
  if( flag_fit ){
    for( Int_t i=0; i<Nhist; i++ ){
      func[i] = new TF1( Form("func%d", i), make_func(sel_fun),offset+xmin_fit,offset+xmax_fit,n_fitfunc_par(sel_fun) );
      func[i]->SetLineColor(i+2);
      if( sel_fun == 10 ){
	Double_t tmp_mu     = hist[i]->GetBinCenter(hist[i]->GetMaximumBin());
	Double_t tmp_sigma  = 35*hist[i]->GetBinWidth(1);
	Double_t tmp_area   = hist[i]->GetMaximum()*sqrt(TMath::TwoPi())*tmp_sigma;
	func[i]->SetParNames  (  "area", "mean", "sigma" );
	func[i]->SetParameters( tmp_area, tmp_mu,tmp_sigma );
      }else if( sel_fun == 11 ){
	Double_t tmp_slope  = (hist[i]->GetBinContent(xbin)-hist[i]->GetBinContent(1))/(hist[i]->GetBinCenter(xbin)-hist[i]->GetBinCenter(1));
	Double_t tmp_offset = hist[i]->GetBinContent(1)-tmp_slope*hist[i]->GetBinCenter(1);
	Double_t tmp_mu     = hist[i]->GetBinCenter(hist[i]->GetMaximumBin());
	Double_t tmp_sigma  = 20*hist[i]->GetBinWidth(1);
	Double_t tmp_area   = (hist[i]->GetMaximum() - (tmp_slope*tmp_mu+tmp_offset))*sqrt(TMath::TwoPi())*tmp_sigma;
	func[i]->SetParNames  ("area",    "mean", "sigma",   "slope",   "offset"    );
	func[i]->SetParameters(tmp_area,  tmp_mu, tmp_sigma, tmp_slope,  tmp_offset );
      }else if( sel_fun == 20 ){
	Double_t tmp_mu     = hist[i]->GetBinCenter(hist[i]->GetMaximumBin());
	Double_t tmp_sigma  = 35*hist[i]->GetBinWidth(1);
	Double_t tmp_area   = hist[i]->GetMaximum()*sqrt(TMath::TwoPi())*tmp_sigma;
	func[i]->SetParNames  (  "area", "area_ratio", "mean", "sigma", "sigma_ratio" );
	func[i]->SetParameters( tmp_area,   0.6,      tmp_mu,  tmp_sigma,   0.35      );
      }else if( sel_fun == 21 ){
	Double_t tmp_slope  = (hist[i]->GetBinContent(xbin)-hist[i]->GetBinContent(1))/(hist[i]->GetBinCenter(xbin)-hist[i]->GetBinCenter(1));
	Double_t tmp_offset = hist[i]->GetBinContent(1)-tmp_slope*hist[i]->GetBinCenter(1);
	Double_t tmp_mu     = hist[i]->GetBinCenter(hist[i]->GetMaximumBin());
	Double_t tmp_sigma  = 20*hist[i]->GetBinWidth(1);
	Double_t tmp_area   = (hist[i]->GetMaximum() - (tmp_slope*tmp_mu+tmp_offset))*sqrt(TMath::TwoPi())*tmp_sigma;
	func[i]->SetParNames  ("area", "area_ratio", "mean", "sigma", "sigma_ratio", "slope", "offset"     );
	func[i]->SetParameters(tmp_area,   0.6,      tmp_mu, tmp_sigma,   0.4,         tmp_slope,  tmp_offset );
      }else{
	func_set_parameters( sel_fun, func[i], hist[i], xbin, offset+xmin, offset+xmax );
      }
    }
  }

  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();

  for(Int_t i=0; i<Nhist; i++ ){
    c1->cd(i+1);
    if( flag_fit ){
      std::cout << Form( "================================== FUNC%d =================================", i ) << std::endl;
      Double_t* init_var = new Double_t[n_fitfunc_par(sel_fun)];
      for( Int_t m=0; m<n_fitfunc_par(sel_fun); m++ ) init_var[m] = func[i]->GetParameter(m);
      iterative_fit( hist[i], func[i], n_fitfunc_par(sel_fun), init_var, "PE0", manip_func );
    }else hist[i]->Draw();
  }

  c1->Update();
  if( flag_save ){
    c1->Print( Form("pic/%s_fun%d_s0%s.eps",  fname, sel_fun, stream) );
    c1->Print( Form("pic/%s_fun%d_s0%s.root", fname, sel_fun, stream) );
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

