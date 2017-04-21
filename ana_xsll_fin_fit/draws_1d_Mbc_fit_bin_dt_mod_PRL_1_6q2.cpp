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
#include <TFile.h>
#include <TArrow.h>
#include <TMath.h>
#include <Math/DistFunc.h>
#include <TRandom.h>


#include <RooFit.h>
#include <RooGlobalFunc.h>
#include <RooRealVar.h>
#include <RooCategory.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooArgSet.h>
#include <RooArgList.h>
#include <RooGaussian.h>
#include <RooArgusBG.h>
#include <RooChebychev.h>
#include <RooFormulaVar.h>
#include <RooAddPdf.h>
#include <RooHistPdf.h>
#include <RooSimultaneous.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooAbsReal.h>
#include <RooGenericPdf.h>
#include <RooCFunction1Binding.h>
#include <RooCFunction2Binding.h>
#include <RooCFunction3Binding.h>
#include <RooCFunction4Binding.h>
#include <RooTFnBinding.h>

using namespace RooFit;
using namespace sig_gmc_rd_cut2;

using namespace Mbc_bkg;
//using namespace Mbc_bkg_wide;

const Bool_t flag_save    = true; // outfile.eps and outfile.root
const Bool_t flag_fit     = true;
const Bool_t flag_roof    = true; // flase(skip roofit)
const Bool_t flag_message = !true;

const Bool_t   flag_fix_end = true; // true(fix endpoint parameter)
const Double_t endpoint     = 5.289;

const Bool_t   flag_scale    = true;
const Bool_t   flag_k4pi     = true; // 1(veto  K4pi    modes)
const Bool_t   flag_unflavor = true; // 1(veto unflavor modes)
const Bool_t   flag_ccpi0    = true; // 1(veto ccpi0 peak    )

void manip_func( TF1* func ){
  if( func->GetNpar()==3 || func->GetNpar()==4 ){ // (modified-)argus
    if( flag_fix_end ) func->FixParameter( 1, endpoint );
    else               func->SetParLimits( 1, 5.280, 5.295 );
  }else if( func->GetNpar()==6 || func->GetNpar()==7 ){ // gaussian + (modified-)argus
    func->FixParameter( 1, 5.2794  );
    func->FixParameter( 2, 0.00273 );
    if( flag_fix_end ) func->FixParameter( 4, endpoint );
    else               func->SetParLimits( 4, 5.280, 5.295 );
  }
}

Int_t main( Int_t argc, Char_t** argv ){
  gROOT->SetBatch(true); // tmpppppp
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==3 || argc==4) ) std::cerr << "wrong input" << std::endl
					<< " Usage : ./draws_1d_Mbc_fit (int)sel_fun (int)fl_q2bin[(int)fl_appRun]" << std::endl
					<< " [sel_fun ] : 50(argus), 51(modified-argus), 15(gauss+argus), 151(gauss+modified-argus)" << std::endl
					<< " [Nbin_afb] : 8 or 14"                                                                   << std::endl
					<< " [fl_q2bin] ; 8[0-3], 14[0-6]"
					<< std::endl, abort();
  Int_t    sel_fun      = atoi(argv[1]);
  Int_t    fl_q2bin     = atoi(argv[2]); // 8[0-3], 14[0-6]
  Int_t    fl_appRun    = 1;
  if( argc==4 ) fl_appRun = atoi( argv[3] );

  const Int_t Nbin_afb = 8; // 8 or 14
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t Nchain        = 2*3; // [gmc,sigmc,rd] x [ee,mm]
  const Int_t Ncategory     = 1;
  const Int_t Nplot         = 1;
  const Int_t Ntmp          = Ncategory*2; // x[ee,mm      ]
  const Int_t Nhist         = Nplot    *3; // x[ee,mm,ee+mm]
  const Int_t fl_sb[Ntmp]      = {2,5}; // 0(bkg,ee), 1(sig,ee), 2(rd,ee), 3(bkg,mm), 4(sig,mm), 5(rd,mm), 
  const Int_t fl_mode_ll[Nchain] = {1,1,1,0,0,0};   // 1(e), 0(mu)

  const Int_t add[Nhist][Ntmp] ={
    {1,0}, // ee
    {0,1}, // mm
    {1,1}, // ee+mm
  };

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  Char_t** infile = new Char_t*[3];
  for( Int_t i=0; i<3; i++ ){
    infile[i] = new Char_t[1024];
    if     ( i==0 ) sTmp << indir[i] << "gMC_*_s0[0-5]";    // bkg (dummy)
    else if( i==1 ) sTmp << indir[i] << "sigMC_*_set[A-U]"; // sig (dummy)
    else if( i==2 ) sTmp << indir[i] << "RD_";              // rd
    strcpy( infile[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Char_t* LRNB_cut = new Char_t[2048];
  strcpy( LRNB_cut, "1" ); // lrnb cut should be already applied.
  std::cout << "LRNB cut : " << LRNB_cut << std::endl;

  Char_t*** add_cut = new Char_t**[Nbin_afb];
  for( Int_t i=0; i<Nbin_afb; i++ ){
    add_cut[i] = new Char_t*[Ntmp];
    for( Int_t j=0; j<Ntmp; j++ ){
      Int_t k = j+i*Ntmp;
      add_cut[i][j] = new Char_t[4096];
      sTmp << LRNB_cut;
      Int_t tmp_fl_mode_ll = 1;
      if( fl_sb[j]>2 ) tmp_fl_mode_ll = 0;
      if( Nbin_afb == 14 ){
	if     ( i ==  0 ) sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 1,  1 ).c_str();
	else if( i ==  1 ) sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 1, -1 ).c_str();
	else if( i ==  2 ) sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 2,  1 ).c_str();
	else if( i ==  3 ) sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 2, -1 ).c_str();
	else if( i ==  4 ) sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 3,  1 ).c_str();
	else if( i ==  5 ) sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 3, -1 ).c_str();
	else if( i ==  6 ) sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 5,  1 ).c_str();
	else if( i ==  7 ) sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 5, -1 ).c_str();
	else if( i ==  8 ) sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 7,  1 ).c_str();
	else if( i ==  9 ) sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 7, -1 ).c_str();
	else if( i == 10 ) sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 8,  1 ).c_str();
	else if( i == 11 ) sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 8, -1 ).c_str();
	else if( i == 12 ) sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 9,  1 ).c_str();
	else if( i == 13 ) sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 9, -1 ).c_str();
      }else if( Nbin_afb == 8 ){
	if     ( i ==  0 ){
	  sTmp << " && ( " << makeCut_q2( tmp_fl_mode_ll, 1,  1 ).c_str()
	       << " || "
	       << makeCut_q2( tmp_fl_mode_ll, 2,  1 ).c_str()
	       << " )";
	}else if( i ==  1 ){
	  sTmp << " && ( " << makeCut_q2( tmp_fl_mode_ll, 1, -1 ).c_str()
	       << " || "
	       << makeCut_q2( tmp_fl_mode_ll, 2, -1 ).c_str()
	       << " )";
	}else if( i ==  2 ){
	  sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 3,  1 ).c_str();
	}else if( i ==  3 ){
	  sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 3, -1 ).c_str();
	}else if( i ==  4 ){
	  sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 5,  1 ).c_str();
	}else if( i ==  5 ){
	  sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 5, -1 ).c_str();
	}else if( i ==  6 ){
	  sTmp << " && ( " << makeCut_q2( tmp_fl_mode_ll, 7,  1 ).c_str()
	       << " || "   << makeCut_q2( tmp_fl_mode_ll, 8,  1 ).c_str()
	       << " || "   << makeCut_q2( tmp_fl_mode_ll, 9,  1 ).c_str()
	       << " )";
	}else if( i ==  7 ){
	  sTmp << " && ( " << makeCut_q2( tmp_fl_mode_ll, 7, -1 ).c_str()
	       << " || "   << makeCut_q2( tmp_fl_mode_ll, 8, -1 ).c_str()
	       << " || "   << makeCut_q2( tmp_fl_mode_ll, 9, -1 ).c_str()
	       << " )";
	}
      }else std::cerr << "[Abort] Wrong Nbin_afb : " << Nbin_afb << std::endl, abort();
      if( flag_k4pi     ) sTmp << " && " << makeCut_5body_veto().c_str();
      if( flag_unflavor ) sTmp << " && " << makeCut_unflavor_veto().c_str();
      if( flag_ccpi0    ) sTmp << " && " << cut_ccpi0;
      strcpy( add_cut[i][j], (Char_t*)sTmp.str().c_str() );
      sTmp.str("");
      sTmp.clear();
    }
  }

  if( flag_k4pi != flag_unflavor ) std::cerr << "[ABORT] Wrong setting flag_k4pi and flag_unflavor" << std::endl, abort();

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**   chain   = new MChain*[Nchain];  // [gmc, sigmc, rd] x [ee, mm]
  TH1D***    tmphist = new TH1D**[Nbin_afb]; // [Ntmp ]
  TH1D***    hist    = new TH1D**[Nbin_afb]; // [Nhist]
  TF1***     func    = new TF1** [Nbin_afb]; // [Nhist]
  Int_t      N0bin[Nbin_afb][Nhist] = {0};   // # of zero-bins for chi2 calculation in RooFit
  for( Int_t i=0; i<Nbin_afb; i++ ){
    tmphist[i] = new TH1D*[Ntmp];
    hist   [i] = new TH1D*[Nhist];
    func   [i] = new TF1* [Nhist];
  }
  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++
  std::cout << Form(" ************************ make chain-tree ( %s, %s ) *************************************",tname,axis) << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<infile %d > ", j );
    if( !(j==2 || j==5) ) continue;
    chain[j] = new MChain( infile[j%3], tname, branch_table(), 0, tail );
    nominal_cut_selection( chain[j], fl_mode_ll[j] )( chain[j]->GetCut(), tname );
  }

  // ++++++++++++++++++++++++
  // cut change
  for( Int_t j=0; j<Nchain; j++ ){
    if( !(j==2 || j==5) ) continue;
    chain[j]->GetCut()->Set( "Mbc", 1,  5.20, 0.0, 5.30 );
  }

  // ++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ changed cut *************************************" << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<infile %d > ", j ) << std::endl;
    if( !(j==2 || j==5) ) continue;
    chain[j]->GetCut()->Display(0);
    chain[j]->MakeTree();
  }

  // +++++++ make tmphist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make tmp hist *************************************" << std::endl;
  for( Int_t i=0; i<Nbin_afb; i++ ){
    for( Int_t j=0; j<Ntmp; j++ ){
      Int_t k = j+i*Ntmp;
      tmphist[i][j] = new TH1D( Form("tmphist%d",k), Form("%s",chain[fl_sb[j]]->GetChange()), xbin,offset+xmin,offset+xmax );
      chain[fl_sb[j]]->GetTree()->Project( Form("tmphist%d",k), axis, add_cut[i][j] );
      if( flag_message ) std::cout << Form( "<tmphist %d > ", k )
				   << "add_cut : " << add_cut[i][j] << std::endl;
      tmphist[i][j]->Sumw2();
    }
  }
  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make hist *************************************" << std::endl;
  for( Int_t i=0; i<Nbin_afb; i++ ){
    for( Int_t j=0; j<Nhist; j++ ){
      Int_t k = j+i*Nhist;
      hist[i][j] = new TH1D( Form("hist%d",k),Form("hist%d",k), xbin,offset+xmin,offset+xmax );
      Deco( hist[i][j], 3, col_fil[j%Nplot], 1 );
      for( Int_t m=0; m<Ntmp; m++ ){
	if( add[j][m] ) hist[i][j]->Add( tmphist[i][m] );
      }
    }
  }
  // +++++++ display ++++++++++++++++++++++++++++++++++
  Double_t entry_sig[Nbin_afb][Nhist] = {0}; // # of events in signal box region
  Double_t total_evt_entry  = 0;
  Double_t total_evt_canvas = 0;
  Double_t total_evt_under  = 0;
  Double_t total_evt_over   = 0;
  Double_t total_evt_sig    = 0;
  for( Int_t i=0; i<Nbin_afb; i++ ){
    for( Int_t j=0; j<Nhist; j++ ){
      Int_t k = j+i*Nhist;
      std::cout << Form("<hist %d> ",k);
      Double_t entry_all    = hist[i][j]->GetEntries();
      Double_t entry_canvas = hist[i][j]->Integral();
      Double_t entry_under  = hist[i][j]->GetBinContent(0);
      Double_t entry_over   = hist[i][j]->GetBinContent(xbin+1);
      for( Int_t m=hist[i][j]->FindBin(5.27+0.000000001); m<=xbin; m++ ) entry_sig[i][j] += hist[i][j]->GetBinContent(m);
      std::cout << "Added-Files( ";
      for( Int_t m=0; m<Ntmp; m++ ) if( add[j][m] ) std::cout << i*Ntmp+m <<  ",";
      std::cout <<  ")["
		<< entry_all
		<< " events ( canvas : "
		<< entry_canvas
		<< " / under : "
		<< entry_under
		<< " / over  : "
		<< entry_over
		<< " / sig  : "
		<< entry_sig[i][j]
		<< "]"
		<< std::endl;
      total_evt_entry  += entry_all;
      total_evt_canvas += entry_canvas;
      total_evt_under  += entry_under;
    total_evt_over   += entry_over;
    total_evt_sig    += entry_sig[i][j];
    }
  }
  std::cout <<  "TOTAL["
	    << total_evt_entry
	    << " events ( canvas : "
	    << total_evt_canvas
	    << " / under : "
	    << total_evt_under
	    << " / over  : "
	    << total_evt_over
	    << " / sig  : "
	    << total_evt_sig
	    << "]"
	    << std::endl << std::endl;
  

  Double_t entry_sig_each[Nbin_afb][Ntmp] = {0}; // # of events in signal box region (tmphist)
  Double_t total_tmp_evt_canvas = 0;
  Double_t total_tmp_evt_sig    = 0;
  for( Int_t i=0; i<Nbin_afb; i++ ){
    for( Int_t j=0; j<Ntmp; j++ ){
      Int_t k = j+i*Ntmp;
      std::cout << Form("<tmphist %d> ",k);
      for( Int_t m=tmphist[i][j]->FindBin(5.27+0.000000001); m<=xbin; m++ ) entry_sig_each[i][j] += tmphist[i][j]->GetBinContent(m);
      std::cout << entry_sig_each[i][j]      << " events(signal-box), "
		<< tmphist[i][j]->Integral() << " events(canvas)"
		<< std::endl;
      total_tmp_evt_canvas += tmphist[i][j]->Integral();
      total_tmp_evt_sig    += entry_sig_each[i][j];
    }
  }
  std::cout << "TOTAL "
	    << total_tmp_evt_sig    << " events(signal-box), "
	    << total_tmp_evt_canvas << " events(canvas)"
	    << std::endl << std::endl;
  
  // +++++++ make fit-func ++++++++++++++++++++++++++++++++++
  if( flag_fit ){
    for( Int_t i=0; i<Nbin_afb; i++ ){
      for( Int_t j=0; j<Nhist; j++ ){
	Int_t k = j+i*Nhist;
	func[i][j] = new TF1( Form("func%d",k), make_func(sel_fun),offset+xmin_fit,offset+xmax_fit,n_fitfunc_par(sel_fun) );
	if( sel_fun==15 ){
	  Double_t mu    = 5.279;
	  Double_t sigma = 0.00255;
	  Double_t area  = ( hist[i][j]->GetMaximum() - hist[i][j]->GetBinContent(1) )*sqrt( TMath::TwoPi() )*sigma;
	  Double_t norm  = 1.5*hist[i][j]->GetBinContent(1);
	  Double_t shape = -15;
	  if( Nbin_afb==8 ){
	    if     ( i==0 && j== 0 ){ area = 0.06979; norm = 21.3;  shape = -13.9;  } // ee
	    else if( i==1 && j== 0 ){ area = 0.07733; norm = 9.23;  shape = -10.33; }
	    else if( i==2 && j== 0 ){ area = 0.06016; norm = 22.1;  shape = -15.07; }
	    else if( i==3 && j== 0 ){ area = 0.04288; norm = 11.46; shape = -24.25; }
	    else if( i==4 && j== 0 ){ area = 0.02805; norm = 6.68;  shape = -12.77; }
	    else if( i==5 && j== 0 ){ area = 0.00985; norm = 4.164; shape = -29.87; }
	    else if( i==6 && j== 0 ){ area = 0.06818; norm = 8.411; shape = -21.13; }
	    else if( i==7 && j== 0 ){ area = 0.05043; norm = 5.097; shape = -27.21; }
	    else if( i==0 && j== 1 ){ area = 0.06825; norm = 14.8;  shape = 0.4128; }// mm
	    else if( i==1 && j== 1 ){ area = 0.0640 ; norm = 9.763; shape = 0.162;  }
	    else if( i==2 && j== 1 ){ area = 0.08339; norm = 36.48; shape = -16.29; }
	    else if( i==3 && j== 1 ){ area = 0.07181; norm = 11.12; shape = -0.3936;}
	    else if( i==4 && j== 1 ){ area = 0.0618;  norm = 22.41; shape = -20.73; }
	    else if( i==5 && j== 1 ){ area = 0.03317; norm = 9.657; shape = -12.15; }
	    else if( i==6 && j== 1 ){ area = 0.09878; norm = 14.06; shape = -26.82; }
	    else if( i==7 && j== 1 ){ area = 0.06771; norm = 7.673; shape = -24.88; }
	    else if( i==0 && j== 2 ){ area = 0.01369; norm = 36.7;  shape = -8.148; } // ee+mm
	    else if( i==1 && j== 2 ){ area = 0.01409; norm = 19.5;  shape = -5.83;  }
	    else if( i==2 && j== 2 ){ area = 0.14;    norm = 60.66; shape = -17.5;  }
	    else if( i==3 && j== 2 ){ area = 0.1138;  norm = 22.97; shape = -11.46; }
	    else if( i==4 && j== 2 ){ area = 0.08928; norm = 29.76; shape = -18.87; }
	    else if( i==5 && j== 2 ){ area = 0.04263; norm = 14.08; shape = -17.05; }
	    else if( i==6 && j== 2 ){ area = 0.1721;  norm = 23.34; shape = -25.98; }
	    else if( i==7 && j== 2 ){ area = 0.1176;  norm = 13.28; shape = -26.5;  }
	  }

	  std::cout << "i = " << i << ", j = " << j << ", area = " << area << ", norm = " << norm << std::endl;
	  func[i][j]->SetParNames  ( "area","mean","sigma","norm",  "Ebeam",  "shape" );
	  func[i][j]->SetParameters( area,   mu,    sigma,  norm,    endpoint, shape  );
	  func[i][j]->FixParameter( 1, 5.2794  );
	  func[i][j]->FixParameter( 2, 0.00273 );
	  if( flag_fix_end ) func[i][j]->FixParameter( 4, endpoint );
	  else               func[i][j]->SetParLimits( 4, 5.280, 5.295 );
	}else if( sel_fun==151 ){
	  Double_t mu    = 5.279;
	  Double_t sigma = 0.00256;
	  Double_t area  = ( hist[i][j]->GetMaximum() - hist[i][j]->GetBinContent(1) )*sqrt( TMath::TwoPi() )*sigma;
	  func[i][j]->SetParNames  ( "area","mean","sigma","norm",                            "Ebeam", "shape",  "new" );
	  func[i][j]->SetParameters( area,   mu,    sigma,  2*hist[i][j]->GetBinContent(1),  endpoint,   -15,     0.46 );
	  func[i][j]->FixParameter( 1, 5.2794 );
	  func[i][j]->FixParameter( 2, 0.00273 );
	  if( flag_fix_end ) func[i][j]->FixParameter( 4, endpoint );
	  else               func[i][j]->SetParLimits( 4, 5.280, 5.295 );
	}else func_set_parameters(sel_fun, func[i][j], hist[i][j], xbin, offset+xmin, offset+xmax);
	func[i][j]->SetLineColor(2);
      }
    }
  }
  
  // +++++++ draw ++++++++++++++++++++++++++++++++++
  TCanvas*** c1 = new TCanvas**[3];
  for( Int_t m=0; m<3; m++ ){ // [ee,mm,ee+mm]
    c1[m] = new TCanvas*[2];
    for( Int_t n=0; n<2; n++ ){ // [hist, fitted hist]
      c1[m][n] = Canvas( Form("c1_%d_%d",m,n), Form("c1_%d_%d",m,n), 4, 4 );
      c1[m][n]->Draw();
    }
  }
  
  for( Int_t i=0; i<Nbin_afb; i++ ){
    // NO FIT
    for(Int_t j=Nhist-1; j>=0; j-- ){
      Int_t k = j+i*Nhist;
      Int_t n;

      sTmp << i/2+1 << " q^{2} bin, ";
      if     ( i%2==0 ) sTmp << "Forward, " ;
      else if( i%2==1 ) sTmp << "Backward, ";
      if     ( 0*Nplot <= j && j < 1*Nplot ) sTmp << "ee",        n = 0; // ee
      else if( 1*Nplot <= j && j < 2*Nplot ) sTmp << "#mu#mu",    n = 1; // mm
      else if( 2*Nplot <= j && j < 3*Nplot ) sTmp << "ee+#mu#mu", n = 2; // ee+mm
      hist[i][j]->SetTitle( sTmp.str().c_str() );
      sTmp.str("");
      sTmp.clear();
      hist[i][j]->SetXTitle(xlabel);
      sTmp << "Events/" <<  Double_t((hist[i][j]->GetXaxis()->GetXmax()-hist[i][j]->GetXaxis()->GetXmin())/xbin);
      hist[i][j]->SetYTitle( sTmp.str().c_str() );
      sTmp.str("");
      sTmp.clear();
      hist[i][j]->GetXaxis()->CenterTitle();
      hist[i][j]->GetYaxis()->CenterTitle();
      c1[n][0]->cd(i+1);
      if( j%Nplot==Nplot-1 ) hist[i][j]->Draw("hist");
      else                   hist[i][j]->Draw("hist same");
    }
    
    
    // FIT
    for(Int_t j=0; j<Nhist; j++ ){
      Int_t k = j+i*Nhist;
      Int_t n;
      if     ( 0*Nplot <= j && j < 1*Nplot ) n = 0;
      else if( 1*Nplot <= j && j < 2*Nplot ) n = 1;
      else if( 2*Nplot <= j && j < 3*Nplot ) n = 2;
      if( flag_fit && j%Nplot==Nplot-1 ){
	c1[n][1]->cd(i+1);
	std::cout << Form( "================================== FUNC%d =================================", k ) << std::endl;
	Double_t init_var[n_fitfunc_par(sel_fun)];
	for( Int_t m=0; m<n_fitfunc_par(sel_fun); m++ ) init_var[m] = func[i][j]->GetParameter(m);
	iterative_fit( hist[i][j], func[i][j], n_fitfunc_par(sel_fun), init_var, "PE0", manip_func );
      }
    }
  }

  if( !flag_roof ){
    for( Int_t m=0; m<3; m++ ){ // [ee,mm,ee+mm]
      for( Int_t n=0; n<2; n++ ){ // [hist, fitted hist]
	Int_t k = n + 2*m;
	c1[m][n]->Update();
	if( fl_q2bin==0 ) c1[m][n]->Print( Form("pic/%s_fit_bin_dt_mod_func%d_c1_%d.eps",  axis, sel_fun, k+1) );
      }
    }
    std::cout << "finish" << std::endl;
    if( fl_appRun ) app.Run();
    return 0;
  }

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++ RooFit ++++++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  for( Int_t i=0; i<Nbin_afb; i++ ){
    for( Int_t j=0; j<Nhist; j++ ){
      Int_t k = j+i*Nhist;
      N0bin[i][j] = 0;
      for( Int_t m=0; m<xbin; m++ ){
	if( hist[i][j]->GetBinContent(m+1)==0 ) N0bin[i][j]++;
      }
      if( flag_message ) std::cout << Form("<hist%d> ", k) << "# of zero-bin : " << N0bin[i][j] << std::endl;
    }
  }
  
  std::cout << std::endl
	    << " *************************************************************" << std::endl
	    << " *********************** RooFit Start ************************" << std::endl
	    << " *************************************************************" << std::endl;

  // --- Observable ---
  const Int_t    Nroohist                 = 3; // ee, mm , ee+mm
  const Int_t    rooad   [Nroohist]       = {0, 1, 2};
  Int_t          evt_cnt [Ntmp]           = {0}; // event counter for each q2-cos bin
  Int_t          evt_cnt2[Nbin_afb][Ntmp] = {0}; // event counter for total bin

  RooRealVar**  obs = new RooRealVar* [Nroohist];
  RooCategory** tag = new RooCategory*[Nroohist];
  for(Int_t i=0; i<Nroohist; i++ ){
    tag[i] = new RooCategory( Form("bin%d",i), Form("bin%d",i) );
    if( Nbin_afb==14 ){
      tag[i]->defineType ( "1st_p",  1 ); tag[i]->defineType ( "1st_m",  2 );
      tag[i]->defineType ( "2nd_p",  3 ); tag[i]->defineType ( "2nd_m",  4 );
      tag[i]->defineType ( "3rd_p",  5 ); tag[i]->defineType ( "3rd_m",  6 );
      tag[i]->defineType ( "5th_p",  7 ); tag[i]->defineType ( "5th_m",  8 );
      tag[i]->defineType ( "7th_p",  9 ); tag[i]->defineType ( "7th_m", 10 );
      tag[i]->defineType ( "8th_p", 11 ); tag[i]->defineType ( "8th_m", 12 );
      tag[i]->defineType ( "9th_p", 13 ); tag[i]->defineType ( "9th_m", 14 );
    }else if( Nbin_afb==8 ){
      tag[i]->defineType ( "12_p",  1 ); tag[i]->defineType ( "12_m",  2 );
      tag[i]->defineType ( "3_p",   3 ); tag[i]->defineType ( "3_m",   4 );
      tag[i]->defineType ( "5_p",   5 ); tag[i]->defineType ( "5_m",   6 );
      tag[i]->defineType ( "789_p", 7 ); tag[i]->defineType ( "789_m", 8 );
    }
  }

 
  for( Int_t i=0; i<Nroohist; i++ ) obs[i] = new RooRealVar( Form("%s%d",axis,i), xlabel, offset+xmin_fit, offset+xmax_fit );

  RooRealVar*  obs_sim = new RooRealVar( Form("%s_sim",axis), xlabel, offset+xmin_fit, offset+xmax_fit );
  RooCategory* tag_sim = new RooCategory( "bin_sim", "bin_sim" );
  if( Nbin_afb==14 ){
    tag_sim->defineType ( "1st_p_e", 101 ); tag_sim->defineType ( "1st_m_e", 102 );
    tag_sim->defineType ( "2nd_p_e", 103 ); tag_sim->defineType ( "2nd_m_e", 104 );
    tag_sim->defineType ( "3rd_p_e", 105 ); tag_sim->defineType ( "3rd_m_e", 106 );
    tag_sim->defineType ( "5th_p_e", 107 ); tag_sim->defineType ( "5th_m_e", 108 );
    tag_sim->defineType ( "7th_p_e", 109 ); tag_sim->defineType ( "7th_m_e", 110 );
    tag_sim->defineType ( "8th_p_e", 111 ); tag_sim->defineType ( "8th_m_e", 112 );
    tag_sim->defineType ( "9th_p_e", 113 ); tag_sim->defineType ( "9th_m_e", 114 );
    tag_sim->defineType ( "1st_p_mu",  1 ); tag_sim->defineType ( "1st_m_mu",  2 );
    tag_sim->defineType ( "2nd_p_mu",  3 ); tag_sim->defineType ( "2nd_m_mu",  4 );
    tag_sim->defineType ( "3rd_p_mu",  5 ); tag_sim->defineType ( "3rd_m_mu",  6 );
    tag_sim->defineType ( "5th_p_mu",  7 ); tag_sim->defineType ( "5th_m_mu",  8 );
    tag_sim->defineType ( "7th_p_mu",  9 ); tag_sim->defineType ( "7th_m_mu", 10 );
    tag_sim->defineType ( "8th_p_mu", 11 ); tag_sim->defineType ( "8th_m_mu", 12 );
    tag_sim->defineType ( "9th_p_mu", 13 ); tag_sim->defineType ( "9th_m_mu", 14 );
  }else if( Nbin_afb==8 ){
    tag_sim->defineType ( "12_p_e",  101 ); tag_sim->defineType ( "12_m_e",  102 );
    tag_sim->defineType ( "3_p_e",   103 ); tag_sim->defineType ( "3_m_e",   104 );
    tag_sim->defineType ( "5_p_e",   105 ); tag_sim->defineType ( "5_m_e",   106 );
    tag_sim->defineType ( "789_p_e", 107 ); tag_sim->defineType ( "789_m_e", 108 );
    tag_sim->defineType ( "12_p_mu",   1 ); tag_sim->defineType ( "12_m_mu",   2 );
    tag_sim->defineType ( "3_p_mu",    3 ); tag_sim->defineType ( "3_m_mu",    4 );
    tag_sim->defineType ( "5_p_mu",    5 ); tag_sim->defineType ( "5_m_mu",    6 );
    tag_sim->defineType ( "789_p_mu",  7 ); tag_sim->defineType ( "789_m_mu",  8 );
  }
  
  Float_t x_obs, cc_obs, cos_obs, rm_l;
  // --- Data Set ---
  RooDataSet** data   = new RooDataSet*[Nroohist];
  TTree***     chain2 = new TTree**    [Nbin_afb]; //[Ntmp]; // apply add_cut[]

  for( Int_t i=0; i<Nroohist;      i++ ) data[i] = new RooDataSet( Form("data%d",i), Form("data%d",i), RooArgSet(*obs[i], *tag[i]) );
  RooDataSet* data_sim = new RooDataSet( "data_sim", "data_sim", RooArgSet(*obs_sim, *tag_sim) );
  for( Int_t i=0; i<Nbin_afb; i++ ){
    chain2[i] = new TTree*[Ntmp];
    for( Int_t j=0; j<Ntmp; j++ ){
      chain2[i][j] = chain[fl_sb[j]]->GetTree()->CopyTree( add_cut[i][j] ); // [Nbin_afb][Nroohist]
    }
  }
  
  for( Int_t i=0; i<Nbin_afb; i++ ){
    for( Int_t j=0; j<Ntmp; j++ ){
      Int_t k = j+i*Ntmp;
      Int_t nevt = chain2[i][j]->GetTree()->GetEntries();
      std::cout << Form("<chain%d> nevt = %d", k, nevt) << std::endl;
      chain2[i][j]->GetTree()->SetBranchAddress( axis,    &x_obs      );
      chain2[i][j]->GetTree()->SetBranchAddress( "cc_m",  &cc_obs     );
      chain2[i][j]->GetTree()->SetBranchAddress( "coslp", &cos_obs    );
      chain2[i][j]->GetTree()->SetBranchAddress( "rm_l",  &rm_l       );
      for( Int_t m=0; m<nevt; m++ ){
	chain2[i][j]->GetTree()->GetEntry(m);
	if( x_obs>=offset+xmin && x_obs<=offset+xmax ){
	  Double_t tmp_rnd = gRandom->Rndm();
	  for( Int_t n=0; n<Nroohist; n++ ){
	    obs[n] ->setVal( x_obs );
	    obs_sim->setVal( x_obs );
	    if( add[rooad[n]][j] ){
	      if(      Nbin_afb == 14 && !((tag_q2_cos_region ((Int_t)rm_l, cc_obs*cc_obs, cos_obs) == 2*fl_q2bin+1) || (tag_q2_cos_region ((Int_t)rm_l, cc_obs*cc_obs, cos_obs) == 2*fl_q2bin+2))  ) continue;
	      else if( Nbin_afb ==  8 && !((tag_q2_cos_region2((Int_t)rm_l, cc_obs*cc_obs, cos_obs) == 2*fl_q2bin+1) || (tag_q2_cos_region2((Int_t)rm_l, cc_obs*cc_obs, cos_obs) == 2*fl_q2bin+2))  ) continue;
	      if( fl_sb[j]==2 || fl_sb[j]==5 ){ // RD
		if     ( Nbin_afb == 14 ) tag[n]->setIndex( tag_q2_cos_region ( (Int_t)rm_l, cc_obs*cc_obs, cos_obs) );
		else if( Nbin_afb ==  8 ) tag[n]->setIndex( tag_q2_cos_region2( (Int_t)rm_l, cc_obs*cc_obs, cos_obs) );
		data[n]->add( RooArgSet(*obs[n],*tag[n]) );
		if( n==Nroohist-1 ){
		  if     ( Nbin_afb == 14 ) tag_sim->setIndex( 100*(Int_t)rm_l + tag_q2_cos_region ( (Int_t)rm_l, cc_obs*cc_obs, cos_obs) );
		  else if( Nbin_afb ==  8 ) tag_sim->setIndex( 100*(Int_t)rm_l + tag_q2_cos_region2( (Int_t)rm_l, cc_obs*cc_obs, cos_obs) );
		  data_sim->add( RooArgSet(*obs_sim,*tag_sim) );
		  evt_cnt[j]++;
		  evt_cnt2[i][j]++;
		}
	      }
	    }
	  }
	}
      }
    }
  }

  for( Int_t i=0; i<Ntmp;     i++ ) std::cout << Form("< evt_cnt %d> ", i )    << evt_cnt[i]
					      << std::endl;
  for( Int_t i=0; i<Nroohist; i++ ) std::cout << Form("< data%d> nevt = ", i ) << data[i]->numEntries()
					      << "( " << data[i]->sumEntries() << " )"
					      << std::endl; // sum of weight
  std::cout << "< data_sim> nevt = " << data_sim->numEntries()
	    << "( " << data_sim->sumEntries() << " )"
	    << std::endl; // sum of weight

  // ------------------------- Fit Fucntion ----------------------------
  const Int_t nparam = 5; // # of fit parameters

  // Sig-PDF(gauss)
  RooRealVar**  sig_mean      = new RooRealVar* [Nroohist];
  RooRealVar**  sig_mean_sim  = new RooRealVar* [Nroohist];
  RooRealVar**  sig_sigma     = new RooRealVar* [Nroohist];
  RooRealVar**  sig_sigma_sim = new RooRealVar* [Nroohist];
  RooGaussian** gauss         = new RooGaussian*[Nroohist];
  RooGaussian** gauss_sim     = new RooGaussian*[Nroohist];
  // Fixed Parameters // 20121016 [gMC(6st) and RD, AFB-modes]
  //sig_mean     [0] = new RooRealVar ( "#mu0",          "#mu",          5.27925); // MC
  //sig_mean     [1] = new RooRealVar ( "#mu1",          "#mu",          5.27922); // MC
  //sig_mean     [2] = new RooRealVar ( "#mu2",          "#mu",          5.27924); // MC
  //sig_sigma    [0] = new RooRealVar ( "#sigma0",       "#sigma",       0.00261); // MC
  //sig_sigma    [1] = new RooRealVar ( "#sigma1",       "#sigma",       0.00256); // MC
  //sig_sigma    [2] = new RooRealVar ( "#sigma2",       "#sigma",       0.00258); // MC
  //sig_mean_sim [0] = new RooRealVar ( "#mu0^{sim}",    "#mu^{sim}",    5.27925); // MC
  //sig_mean_sim [1] = new RooRealVar ( "#mu1^{sim}",    "#mu^{sim}",    5.27922); // MC
  //sig_mean_sim [2] = new RooRealVar ( "#mu2^{sim}",    "#mu^{sim}",    5.27924); // MC
  //sig_sigma_sim[0] = new RooRealVar ( "#sigma0^{sim}", "#sigma^{sim}", 0.00261); // MC
  //sig_sigma_sim[1] = new RooRealVar ( "#sigma1^{sim}", "#sigma^{sim}", 0.00256); // MC
  //sig_sigma_sim[2] = new RooRealVar ( "#sigma2^{sim}", "#sigma^{sim}", 0.00258); // MC
  sig_mean     [0] = new RooRealVar ( "#mu0",          "#mu",          5.27936); // RD
  sig_mean     [1] = new RooRealVar ( "#mu1",          "#mu",          5.27932); // RD
  sig_mean     [2] = new RooRealVar ( "#mu2",          "#mu",          5.27934); // RD
  sig_sigma    [0] = new RooRealVar ( "#sigma0",       "#sigma",       0.00267); // RD
  sig_sigma    [1] = new RooRealVar ( "#sigma1",       "#sigma",       0.00258); // RD
  sig_sigma    [2] = new RooRealVar ( "#sigma2",       "#sigma",       0.00263); // RD
  sig_mean_sim [0] = new RooRealVar ( "#mu0^{sim}",    "#mu^{sim}",    5.27936); // RD
  sig_mean_sim [1] = new RooRealVar ( "#mu1^{sim}",    "#mu^{sim}",    5.27932); // RD
  sig_mean_sim [2] = new RooRealVar ( "#mu2^{sim}",    "#mu^{sim}",    5.27934); // RD
  sig_sigma_sim[0] = new RooRealVar ( "#sigma0^{sim}", "#sigma^{sim}", 0.00267); // RD
  sig_sigma_sim[1] = new RooRealVar ( "#sigma1^{sim}", "#sigma^{sim}", 0.00258); // RD
  sig_sigma_sim[2] = new RooRealVar ( "#sigma2^{sim}", "#sigma^{sim}", 0.00263); // RD
  for( Int_t i=0; i<Nroohist; i++ ) gauss    [i] = new RooGaussian( Form("gauss%d",      i), "gauss(x,mean,sigma)", *obs[i],  *sig_mean    [i], *sig_sigma    [i] );
  for( Int_t i=0; i<Nroohist; i++ ) gauss_sim[i] = new RooGaussian( Form("gauss^{sim}%d",i), "gauss(x,mean,sigma)", *obs_sim, *sig_mean_sim[i], *sig_sigma_sim[i] );

  // Bkg-PDF
  RooRealVar*** arg_end       = new RooRealVar**[Nroohist]; // [Nbin_afb]
  RooRealVar*** arg_end_sim   = new RooRealVar**[Nroohist]; // [Nbin_afb]
  RooRealVar*** arg_shape     = new RooRealVar**[Nroohist]; // [Nbin_afb]
  RooRealVar*** arg_shape_sim = new RooRealVar**[Nroohist]; // [Nbin_afb]
  RooRealVar*** arg_new       = new RooRealVar**[Nroohist]; // [Nbin_afb]
  RooRealVar*** arg_new_sim   = new RooRealVar**[Nroohist]; // [Nbin_afb]
  RooAbsPdf***  modargus      = new RooAbsPdf** [Nroohist]; // [Nbin_afb]
  RooAbsPdf***  modargus_sim  = new RooAbsPdf** [Nroohist]; // [Nbin_afb]
  for( Int_t i=0; i<Nroohist; i++ ){
    arg_end      [i] = new RooRealVar*[Nbin_afb];
    arg_end_sim  [i] = new RooRealVar*[Nbin_afb];
    arg_shape    [i] = new RooRealVar*[Nbin_afb];
    arg_shape_sim[i] = new RooRealVar*[Nbin_afb];
    arg_new      [i] = new RooRealVar*[Nbin_afb];
    arg_new_sim  [i] = new RooRealVar*[Nbin_afb];
    modargus     [i] = new RooAbsPdf* [Nbin_afb];
    modargus_sim [i] = new RooAbsPdf* [Nbin_afb];
    for( Int_t j=0; j<Nbin_afb; j++ ){
      Int_t k = j+i*Nbin_afb;
      arg_end    [i][j] = new RooRealVar( Form("end%d",k),   "argus endpoint", endpoint, 5.285, 5.290 );
      arg_end_sim[i][j] = new RooRealVar( Form("end%d",k),   "argus endpoint", endpoint, 5.285, 5.290 );
      if( flag_fix_end ) arg_end    [i][j]->setConstant(kTRUE);
      if( flag_fix_end ) arg_end_sim[i][j]->setConstant(kTRUE);
      if( sel_fun== 15 ){
	Double_t tmp_shape = func[j][rooad[i]]->GetParameter(5);
	arg_shape    [i][j] = new RooRealVar( Form("shape%d_%d",      i,j), "argus shape parameter", tmp_shape, tmp_shape-100.0, tmp_shape+100.0 );
	arg_shape_sim[i][j] = new RooRealVar( Form("shape^{sim}%d_%d",i,j), "argus shape parameter", tmp_shape, tmp_shape-100.0, tmp_shape+100.0 );
	//arg_shape    [i][j]->setConstant(kTRUE); // tmpppppppp
	//arg_shape_sim[i][j]->setConstant(kTRUE); // tmpppppppp
	modargus     [i][j] = bindPdf( Form("modargus%d",      k), func_roo_argus, *obs[i],  *arg_end    [i][j], *arg_shape    [i][j] );
	modargus_sim [i][j] = bindPdf( Form("modargus^{sim}%d",k), func_roo_argus, *obs_sim, *arg_end_sim[i][j], *arg_shape_sim[i][j] );
      }else if( sel_fun==151 ){
	Double_t tmp_shape = func[j][rooad[i]]->GetParameter(5);
	Double_t tmp_new   = func[j][rooad[i]]->GetParameter(6);
	arg_shape    [i][j] = new RooRealVar( Form("shape%d_%d",      i,j), "argus shape parameter", tmp_shape, tmp_shape-20.0, tmp_shape+20.0 );
	arg_shape_sim[i][j] = new RooRealVar( Form("shape^{sim}%d_%d",i,j), "argus shape parameter", tmp_shape, tmp_shape-20.0, tmp_shape+20.0 );
	arg_new      [i][j] = new RooRealVar( Form("new%d_%d",        i,j), "argus sqrt parameter",  tmp_new,   tmp_new  -0.15, tmp_new  +0.15 );
	arg_new_sim  [i][j] = new RooRealVar( Form("new^{sim}%d_%d",  i,j), "argus sqrt parameter",  tmp_new,   tmp_new  -0.15, tmp_new  +0.15 );
	modargus     [i][j] = bindPdf( Form("modargus%d",      k), func_roo_modargus,   *obs[i],  *arg_end    [i][j], *arg_shape    [i][j], *arg_new    [i][j] );
	modargus_sim [i][j] = bindPdf( Form("modargus^{sim}%d",k), func_roo_modargus,   *obs_sim, *arg_end_sim[i][j], *arg_shape_sim[i][j], *arg_new_sim[i][j] );
      }
    }
  }


  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // READ PDF (self cross-feed)
  //TFile file_scf("pdf/Mbc_self_cf_setA-U.root");
  //TFile file_scf("pdf/Mbc_self_cf2_setA-U.root");
  TFile file_scf("pdf_for_1_6q2/Mbc_self_cf2_setA-U.root");
  if( file_scf.IsZombie() ) std::cerr << "[ABORT] can not find input-file : " << file_scf.GetName() << std::endl, abort();

  TH1D**** hist_scf = new TH1D***[Nroohist]; // [Nbin_afb][4]

  for( Int_t i=0; i<Nroohist; i++ ){
    hist_scf[i] = new TH1D**[Nbin_afb];
    for( Int_t j=0; j<Nbin_afb; j++ ){
      hist_scf[i][j] = new TH1D*[4]; // true,(t,t), (t,f), (f,?)
      for( Int_t k=0; k<4; k++ ){
	Int_t m = k+4*(j+i*Nbin_afb);
	hist_scf[i][j][k] = (TH1D*)file_scf.Get( Form("hist%d",m) ); 
	if( hist_scf[i][j][k] == NULL ) std::cerr << "[ABORT] can not find histgram (scf)" << k << std::endl, abort();
      }
    }
  }
  RooDataHist*** scf_tt         = new RooDataHist**[Nroohist]; // [Nbin_afb]
  RooDataHist*** scf_tt_sim     = new RooDataHist**[Nroohist]; // [Nbin_afb]
  RooDataHist*** scf_tf         = new RooDataHist**[Nroohist]; // [Nbin_afb]
  RooDataHist*** scf_tf_sim     = new RooDataHist**[Nroohist]; // [Nbin_afb]
  RooDataHist*** scf_f          = new RooDataHist**[Nroohist]; // [Nbin_afb]
  RooDataHist*** scf_f_sim      = new RooDataHist**[Nroohist]; // [Nbin_afb]
  RooHistPdf***  pdf_scf_tt     = new RooHistPdf** [Nroohist]; // [Nbin_afb]
  RooHistPdf***  pdf_scf_tt_sim = new RooHistPdf** [Nroohist]; // [Nbin_afb]
  RooHistPdf***  pdf_scf_tf     = new RooHistPdf** [Nroohist]; // [Nbin_afb]
  RooHistPdf***  pdf_scf_tf_sim = new RooHistPdf** [Nroohist]; // [Nbin_afb]
  RooHistPdf***  pdf_scf_f      = new RooHistPdf** [Nroohist]; // [Nbin_afb]
  RooHistPdf***  pdf_scf_f_sim  = new RooHistPdf** [Nroohist]; // [Nbin_afb]

  for( Int_t i=0; i<Nroohist; i++ ){
    scf_tt        [i] = new RooDataHist*[Nbin_afb];
    scf_tt_sim    [i] = new RooDataHist*[Nbin_afb];
    scf_tf        [i] = new RooDataHist*[Nbin_afb];
    scf_tf_sim    [i] = new RooDataHist*[Nbin_afb];
    scf_f         [i] = new RooDataHist*[Nbin_afb];
    scf_f_sim     [i] = new RooDataHist*[Nbin_afb];
    pdf_scf_tt    [i] = new RooHistPdf* [Nbin_afb];
    pdf_scf_tt_sim[i] = new RooHistPdf* [Nbin_afb];
    pdf_scf_tf    [i] = new RooHistPdf* [Nbin_afb];
    pdf_scf_tf_sim[i] = new RooHistPdf* [Nbin_afb];
    pdf_scf_f     [i] = new RooHistPdf* [Nbin_afb];
    pdf_scf_f_sim [i] = new RooHistPdf* [Nbin_afb];

    for( Int_t j=0; j<Nbin_afb; j++ ){
      Int_t k = j+i*Nbin_afb;
      scf_tt        [i][j] = new RooDataHist( Form("scf_tt%d",          k), Form("scf_tt%d",          k), *obs[i],  hist_scf[i][j][1]  );
      scf_tt_sim    [i][j] = new RooDataHist( Form("scf_tt^{sim}%d",    k), Form("scf_tt^{sim}%d",    k), *obs_sim, hist_scf[i][j][1]  );
      scf_tf        [i][j] = new RooDataHist( Form("scf_tf%d",          k), Form("scf_tf%d",          k), *obs[i],  hist_scf[i][j][2]  );
      scf_tf_sim    [i][j] = new RooDataHist( Form("scf_tf^{sim}%d",    k), Form("scf_tf^{sim}%d",    k), *obs_sim, hist_scf[i][j][2]  );
      scf_f         [i][j] = new RooDataHist( Form("scf_f%d",           k), Form("scf_f%d",           k), *obs[i],  hist_scf[i][j][3]  );
      scf_f_sim     [i][j] = new RooDataHist( Form("scf_f^{sim}%d",     k), Form("scf_f^{sim}%d",     k), *obs_sim, hist_scf[i][j][3]  );
      pdf_scf_tt    [i][j] = new RooHistPdf ( Form("pdf_scf_tt%d",      k), Form("pdf_scf_tt%d",      k), *obs[i],  *scf_tt [i][j]     );
      pdf_scf_tt_sim[i][j] = new RooHistPdf ( Form("pdf_scf_tt^{sim}%d",k), Form("pdf_scf_tt^{sim}%d",k), *obs_sim, *scf_tt_sim [i][j] );
      pdf_scf_tf    [i][j] = new RooHistPdf ( Form("pdf_scf_tf%d",      k), Form("pdf_scf_tf%d",      k), *obs[i],  *scf_tf     [i][j] );
      pdf_scf_tf_sim[i][j] = new RooHistPdf ( Form("pdf_scf_tf^{sim}%d",k), Form("pdf_scf_tf^{sim}%d",k), *obs_sim, *scf_tf_sim [i][j] );
      pdf_scf_f     [i][j] = new RooHistPdf ( Form("pdf_scf_f%d",       k), Form("pdf_scf_f%d",       k), *obs[i],  *scf_f      [i][j] );
      pdf_scf_f_sim [i][j] = new RooHistPdf ( Form("pdf_scf_f^{sim}%d", k), Form("pdf_scf_f^{sim}%d", k), *obs_sim, *scf_f_sim  [i][j] );
    }
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // SETTING of COEFFICIENT (signal and self cross-feed)
  RooFormulaVar*** nsig        = new RooFormulaVar**[Nroohist]; //[Nbin_afb]
  RooRealVar***    nbkg        = new RooRealVar   **[Nroohist]; //[Nbin_afb]
  RooFormulaVar*** nsig_tt     = new RooFormulaVar**[Nroohist]; //[Nbin_afb]
  RooFormulaVar*** nsig_tf     = new RooFormulaVar**[Nroohist]; //[Nbin_afb]
  RooRealVar***    r_tt        = new RooRealVar   **[Nroohist]; //[Nbin_afb]
  RooRealVar***    r_tf        = new RooRealVar   **[Nroohist]; //[Nbin_afb]
  RooRealVar***    r_f         = new RooRealVar   **[Nroohist]; //[Nbin_afb]
  RooRealVar***    nsig_q2     = new RooRealVar   **[Nroohist]; //[Nbin_afb/2]
  RooRealVar***    AFB         = new RooRealVar   **[Nroohist]; //[Nbin_afb/2]

  RooFormulaVar*** nsig_sim     = new RooFormulaVar**[Nroohist]; //[Nbin_afb]
  RooRealVar***    nbkg_sim     = new RooRealVar   **[Nroohist]; //[Nbin_afb]
  RooFormulaVar*** nsig_tt_sim  = new RooFormulaVar**[Nroohist]; //[Nbin_afb]
  RooFormulaVar*** nsig_tf_sim  = new RooFormulaVar**[Nroohist]; //[Nbin_afb]
  RooRealVar***    nsig_q2_sim  = new RooRealVar   **[Nroohist]; //[Nbin_afb/2]
  RooFormulaVar*** AFB_meas     = new RooFormulaVar**[Nroohist]; //[Nbin_afb]
  RooRealVar**     AFB_true     = new RooRealVar   * [Nbin_afb/2];
  
  for( Int_t i=0; i<Nroohist; i++ ){
    r_tt[i] = new RooRealVar*[Nbin_afb];
    r_tf[i] = new RooRealVar*[Nbin_afb];
    r_f [i] = new RooRealVar*[Nbin_afb];
  }

  if( Nbin_afb == 14 ){
  }else if( Nbin_afb == 8 ){ // for 1-6 q2 measurement
    r_tt[ 0][ 0] = new RooRealVar( "r_tt0",    "r_tt0",     0.08758);   r_tf[ 0][ 0] = new RooRealVar( "r_tf0",    "r_tf0",     0.01574);   r_f [ 0][ 0] = new RooRealVar( "r_f0",     "r_f0",    0.0001533); 
    r_tt[ 0][ 1] = new RooRealVar( "r_tt1",    "r_tt1",     0.08964);   r_tf[ 0][ 1] = new RooRealVar( "r_tf1",    "r_tf1",     0.01483);   r_f [ 0][ 1] = new RooRealVar( "r_f1",     "r_f1",    0.0001631); 
    r_tt[ 0][ 2] = new RooRealVar( "r_tt2",    "r_tt2",     0.09527);   r_tf[ 0][ 2] = new RooRealVar( "r_tf2",    "r_tf2",      0.0199);   r_f [ 0][ 2] = new RooRealVar( "r_f2",     "r_f2",    0.0004972); 
    r_tt[ 0][ 3] = new RooRealVar( "r_tt3",    "r_tt3",     0.09526);   r_tf[ 0][ 3] = new RooRealVar( "r_tf3",    "r_tf3",     0.02027);   r_f [ 0][ 3] = new RooRealVar( "r_f3",     "r_f3",    0.0005379); 
    r_tt[ 0][ 4] = new RooRealVar( "r_tt4",    "r_tt4",       0.154);   r_tf[ 0][ 4] = new RooRealVar( "r_tf4",    "r_tf4",     0.04828);   r_f [ 0][ 4] = new RooRealVar( "r_f4",     "r_f4",    0.0002512); 
    r_tt[ 0][ 5] = new RooRealVar( "r_tt5",    "r_tt5",      0.1461);   r_tf[ 0][ 5] = new RooRealVar( "r_tf5",    "r_tf5",     0.04598);   r_f [ 0][ 5] = new RooRealVar( "r_f5",     "r_f5",    0.0002309); 
    r_tt[ 0][ 6] = new RooRealVar( "r_tt6",    "r_tt6",      0.2157);   r_tf[ 0][ 6] = new RooRealVar( "r_tf6",    "r_tf6",      0.1213);   r_f [ 0][ 6] = new RooRealVar( "r_f6",     "r_f6",     0.001751); 
    r_tt[ 0][ 7] = new RooRealVar( "r_tt7",    "r_tt7",      0.1928);   r_tf[ 0][ 7] = new RooRealVar( "r_tf7",    "r_tf7",      0.1107);   r_f [ 0][ 7] = new RooRealVar( "r_f7",     "r_f7",     0.002015); 
    r_tt[ 1][ 0] = new RooRealVar( "r_tt8",    "r_tt0",     0.06484);   r_tf[ 1][ 0] = new RooRealVar( "r_tf8",    "r_tf0",     0.01191);   r_f [ 1][ 0] = new RooRealVar( "r_f8",     "r_f0",    7.134e-05); 
    r_tt[ 1][ 1] = new RooRealVar( "r_tt9",    "r_tt1",     0.06478);   r_tf[ 1][ 1] = new RooRealVar( "r_tf9",    "r_tf1",     0.01277);   r_f [ 1][ 1] = new RooRealVar( "r_f9",     "r_f1",    7.001e-05); 
    r_tt[ 1][ 2] = new RooRealVar( "r_tt10",   "r_tt2",     0.07951);   r_tf[ 1][ 2] = new RooRealVar( "r_tf10",   "r_tf2",     0.01881);   r_f [ 1][ 2] = new RooRealVar( "r_f10",    "r_f2",    0.0003771); 
    r_tt[ 1][ 3] = new RooRealVar( "r_tt11",   "r_tt3",      0.0804);   r_tf[ 1][ 3] = new RooRealVar( "r_tf11",   "r_tf3",     0.01766);   r_f [ 1][ 3] = new RooRealVar( "r_f11",    "r_f3",     0.000358); 
    r_tt[ 1][ 4] = new RooRealVar( "r_tt12",   "r_tt4",      0.1244);   r_tf[ 1][ 4] = new RooRealVar( "r_tf12",   "r_tf4",     0.03945);   r_f [ 1][ 4] = new RooRealVar( "r_f12",    "r_f4",    0.0005924); 
    r_tt[ 1][ 5] = new RooRealVar( "r_tt13",   "r_tt5",      0.1182);   r_tf[ 1][ 5] = new RooRealVar( "r_tf13",   "r_tf5",     0.04148);   r_f [ 1][ 5] = new RooRealVar( "r_f13",    "r_f5",    0.0005955); 
    r_tt[ 1][ 6] = new RooRealVar( "r_tt14",   "r_tt6",      0.1399);   r_tf[ 1][ 6] = new RooRealVar( "r_tf14",   "r_tf6",      0.0773);   r_f [ 1][ 6] = new RooRealVar( "r_f14",    "r_f6",     0.001496); 
    r_tt[ 1][ 7] = new RooRealVar( "r_tt15",   "r_tt7",      0.1249);   r_tf[ 1][ 7] = new RooRealVar( "r_tf15",   "r_tf7",     0.07171);   r_f [ 1][ 7] = new RooRealVar( "r_f15",    "r_f7",     0.001631); 
    r_tt[ 2][ 0] = new RooRealVar( "r_tt16",   "r_tt0",     0.07811);   r_tf[ 2][ 0] = new RooRealVar( "r_tf16",   "r_tf0",     0.01417);   r_f [ 2][ 0] = new RooRealVar( "r_f16",    "r_f0",    0.0002247); 
    r_tt[ 2][ 1] = new RooRealVar( "r_tt17",   "r_tt1",     0.07947);   r_tf[ 2][ 1] = new RooRealVar( "r_tf17",   "r_tf1",     0.01398);   r_f [ 2][ 1] = new RooRealVar( "r_f17",    "r_f1",    0.0002331); 
    r_tt[ 2][ 2] = new RooRealVar( "r_tt18",   "r_tt2",     0.08802);   r_tf[ 2][ 2] = new RooRealVar( "r_tf18",   "r_tf2",     0.01941);   r_f [ 2][ 2] = new RooRealVar( "r_f18",    "r_f2",    0.0008743); 
    r_tt[ 2][ 3] = new RooRealVar( "r_tt19",   "r_tt3",     0.08854);   r_tf[ 2][ 3] = new RooRealVar( "r_tf19",   "r_tf3",     0.01907);   r_f [ 2][ 3] = new RooRealVar( "r_f19",    "r_f3",     0.000896); 
    r_tt[ 2][ 4] = new RooRealVar( "r_tt20",   "r_tt4",      0.1331);   r_tf[ 2][ 4] = new RooRealVar( "r_tf20",   "r_tf4",     0.04207);   r_f [ 2][ 4] = new RooRealVar( "r_f20",    "r_f4",    0.0008437); 
    r_tt[ 2][ 5] = new RooRealVar( "r_tt21",   "r_tt5",      0.1265);   r_tf[ 2][ 5] = new RooRealVar( "r_tf21",   "r_tf5",      0.0428);   r_f [ 2][ 5] = new RooRealVar( "r_f21",    "r_f5",    0.0008264); 
    r_tt[ 2][ 6] = new RooRealVar( "r_tt22",   "r_tt6",      0.1728);   r_tf[ 2][ 6] = new RooRealVar( "r_tf22",   "r_tf6",      0.0966);   r_f [ 2][ 6] = new RooRealVar( "r_f22",    "r_f6",     0.003247); 
    r_tt[ 2][ 7] = new RooRealVar( "r_tt23",   "r_tt7",      0.1547);   r_tf[ 2][ 7] = new RooRealVar( "r_tf23",   "r_tf7",     0.08863);   r_f [ 2][ 7] = new RooRealVar( "r_f23",    "r_f7",     0.003646); 
  }

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // SETTING of INITIAL VALUE
  for( Int_t i=0; i<Nroohist; i++ ){
    nsig_q2    [i] = new RooRealVar*[Nbin_afb/2];
    nsig_q2_sim[i] = new RooRealVar*[Nbin_afb/2];
    AFB        [i] = new RooRealVar*[Nbin_afb/2];
    for( Int_t j=0; j<Nbin_afb/2; j++ ){
      Int_t k = j + i*Nbin_afb/2;
      Double_t tmp_nsig_q2 = 30;
      Double_t tmp_AFB     =  0;
      nsig_q2    [i][j] = new RooRealVar( Form("N_{sig}^{q^{2}}%d_%d",    i,j), Form("N_{sig}^{q^{2}}%d_%d",    i,j), tmp_nsig_q2, tmp_nsig_q2-100, tmp_nsig_q2+100 );
      nsig_q2_sim[i][j] = new RooRealVar( Form("N_{sig}^{q^{2},sim}%d_%d",i,j), Form("N_{sig}^{q^{2},sim}%d_%d",i,j), tmp_nsig_q2, tmp_nsig_q2-100, tmp_nsig_q2+100 );
      AFB        [i][j] = new RooRealVar( Form("A_{FB}%d_%d",             i,j), Form("A_{FB}%d_%d",             i,j), tmp_AFB,     -1.5,           1.5            );
      AFB_true      [j] = new RooRealVar( Form("A_{FB}^{true}_%d",          j), Form("A_{FB}^{true}_%d",          j), tmp_AFB,     -1.5,           1.5            );
      std::cout << Form("[i=%d, j=%d] nsig_q2 = %f, tmp_AFB = %f", i,j,tmp_nsig_q2,tmp_AFB) << std::endl;
    }
  }
  
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  RooRealVar*** cf_slope = new RooRealVar**[Nroohist];
  for( Int_t i=0; i<Nroohist;   i++ ) cf_slope[i]    = new RooRealVar*[Nbin_afb/2];
  for( Int_t j=0; j<Nbin_afb/2; j++ ) cf_slope[2][j] = new RooRealVar( Form("cf_2_%d",j), Form("cf_2_%d",j), 0.0 ); // dummy

  if( Nbin_afb == 14 ){
  }else if( Nbin_afb == 8 ){
    ///* 6qbin with correction for RD for 1-6 q2 measurement
    cf_slope[ 0][ 0] = new RooRealVar( "cf_0_0", "cf_0_0",    1.53866 );
    cf_slope[ 0][ 1] = new RooRealVar( "cf_0_1", "cf_0_1",    1.25511 );
    //cf_slope[ 0][ 1] = new RooRealVar( "cf_0_1", "cf_0_1",    1.24664 ); // for systr from flip
    cf_slope[ 0][ 2] = new RooRealVar( "cf_0_2", "cf_0_2",    1.06319 );
    cf_slope[ 0][ 3] = new RooRealVar( "cf_0_3", "cf_0_3",    1.12141 );
    cf_slope[ 1][ 0] = new RooRealVar( "cf_1_0", "cf_1_0",    2.16957 );
    cf_slope[ 1][ 1] = new RooRealVar( "cf_1_1", "cf_1_1",    1.86334 );
    //cf_slope[ 1][ 1] = new RooRealVar( "cf_1_1", "cf_1_1",    1.84830 ); // for syst from flip
    cf_slope[ 1][ 2] = new RooRealVar( "cf_1_2", "cf_1_2",    1.03293 );
    cf_slope[ 1][ 3] = new RooRealVar( "cf_1_3", "cf_1_3",    1.08232 );
    //*/
  }
  /*
  RooRealVar*** cf_offset = new RooRealVar**[Nroohist];
  for( Int_t i=0; i<Nroohist;   i++ ) cf_offset[i]    = new RooRealVar*[Nbin_afb/2];
  for( Int_t j=0; j<Nbin_afb/2; j++ ) cf_offset[2][j] = new RooRealVar( Form("cf_2_%d",j), Form("cf_2_%d",j), 0.0 ); // dummy
  cf_offset[ 0][ 0] = new RooRealVar( "cfo_0_0", "cfo_0_0",    0.007238  ); // dummy
  cf_offset[ 0][ 1] = new RooRealVar( "cfo_0_1", "cfo_0_1",    0.001633  ); // dummy
  cf_offset[ 0][ 2] = new RooRealVar( "cfo_0_2", "cfo_0_2",    0.0063933 ); // dummy
  cf_offset[ 0][ 3] = new RooRealVar( "cfo_0_3", "cfo_0_3",    0.003193  ); // dummy
  cf_offset[ 1][ 0] = new RooRealVar( "cfo_1_0", "cfo_1_0",    0.0065723 ); // dummy
  cf_offset[ 1][ 1] = new RooRealVar( "cfo_1_1", "cfo_1_1",    0.007965  ); // dummy
  cf_offset[ 1][ 2] = new RooRealVar( "cfo_1_2", "cfo_1_2",    0.00012404); // dummy
  cf_offset[ 1][ 3] = new RooRealVar( "cfo_1_3", "cfo_1_3",    0.0012754 ); // dummy

  cf_offset[ 0][ 1] = new RooRealVar( "cfo_0_1", "cfo_0_1",    0.008917 );
  cf_offset[ 1][ 1] = new RooRealVar( "cfo_1_1", "cfo_1_1",    0.001063 );
  */
  for( Int_t i=0; i<Nroohist; i++ ){
    AFB_meas[i] = new RooFormulaVar*[Nbin_afb];
    for( Int_t j=0; j<Nbin_afb/2; j++ ){
      //AFB_meas[i][j] = new RooFormulaVar( Form("A_{FB}^{meas.}_%d_%d",i,j),  "1/@0*@1", RooArgSet(*cf_slope[i][j], *AFB_true[j]) );
      //if(      i==0 && j==1 ) AFB_meas[i][j] = new RooFormulaVar( Form("A_{FB}^{meas.}_%d_%d",i,j),  "0.131560/0.110736/@0*@1", RooArgSet(*cf_slope[i][j], *AFB_true[j]) ); // modified @ 20130813
      //else if( i==0 && j==2 ) AFB_meas[i][j] = new RooFormulaVar( Form("A_{FB}^{meas.}_%d_%d",i,j),  "0.321939/0.319769/@0*@1", RooArgSet(*cf_slope[i][j], *AFB_true[j]) ); // modified @ 20130813
      //if(      i==0 && j==1 ) AFB_meas[i][j] = new RooFormulaVar( Form("A_{FB}^{meas.}_%d_%d",i,j),  "1.01868/@0*@1", RooArgSet(*cf_slope[i][j], *AFB_true[j]) ); // modified @ 20130924
      //else if( i==0 && j==2 ) AFB_meas[i][j] = new RooFormulaVar( Form("A_{FB}^{meas.}_%d_%d",i,j),  "1.0034/@0*@1", RooArgSet(*cf_slope[i][j], *AFB_true[j]) ); // modified @ 20130924
      //else                    AFB_meas[i][j] = new RooFormulaVar( Form("A_{FB}^{meas.}_%d_%d",i,j),  "1/@0*@1", RooArgSet(*cf_slope[i][j], *AFB_true[j]) ); // modified @ 20130813
      
      if( i==0 && j==2 ) AFB_meas[i][j] = new RooFormulaVar( Form("A_{FB}^{meas.}_%d_%d",i,j),  "1.0034/@0*@1", RooArgSet(*cf_slope[i][j], *AFB_true[j]) ); // for 1-6 q2 measurement
      else               AFB_meas[i][j] = new RooFormulaVar( Form("A_{FB}^{meas.}_%d_%d",i,j),  "1/@0*@1", RooArgSet(*cf_slope[i][j], *AFB_true[j]) );
      //AFB_meas[i][j] = new RooFormulaVar( Form("A_{FB}^{meas.}_%d_%d",i,j),  "1/@0*(@1-@2*sqrt(1+@0*@0))", RooArgSet(*cf_slope[i][j], *AFB_true[j], *cf_offset[i][j]) );
    }
  }
  
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  

  for( Int_t i=0; i<Nroohist; i++ ){
    nsig    [i] = new RooFormulaVar*[Nbin_afb];
    nsig_sim[i] = new RooFormulaVar*[Nbin_afb];
    for( Int_t j=0; j<Nbin_afb/2; j++ ){
      nsig    [i][2*j+0] = new RooFormulaVar( Form("nsig_%d_%d",      i,2*j+0),  "@0*(1+@1)/2", RooArgSet(*nsig_q2    [i][j], *AFB     [i][j]) );
      nsig    [i][2*j+1] = new RooFormulaVar( Form("nsig_%d_%d",      i,2*j+1),  "@0*(1-@1)/2", RooArgSet(*nsig_q2    [i][j], *AFB     [i][j]) );
      nsig_sim[i][2*j+0] = new RooFormulaVar( Form("nsig^{sim}_%d_%d",i,2*j+0),  "@0*(1+@1)/2", RooArgSet(*nsig_q2_sim[i][j], *AFB_meas[i][j]) );
      nsig_sim[i][2*j+1] = new RooFormulaVar( Form("nsig^{sim}_%d_%d",i,2*j+1),  "@0*(1-@1)/2", RooArgSet(*nsig_q2_sim[i][j], *AFB_meas[i][j]) );
      std::cout << Form("[Formula, i=%d, j=%d] nsig = %f, nsig = %f", i,j,nsig[i][2*j+0]->getVal(),nsig[i][2*j+1]->getVal()) << std::endl;
    }
  }

  for( Int_t i=0; i<Nroohist; i++ ){
    nsig_tt    [i] = new RooFormulaVar*[Nbin_afb];
    nsig_tf    [i] = new RooFormulaVar*[Nbin_afb];
    nsig_tt_sim[i] = new RooFormulaVar*[Nbin_afb];
    nsig_tf_sim[i] = new RooFormulaVar*[Nbin_afb];
    for( Int_t j=0; j<Nbin_afb/2; j++ ){
      nsig_tt    [i][2*j+0] = new RooFormulaVar( Form("nsig_tt_%d_%d",    i,2*j+0), "@0*@1", RooArgSet(*r_tt[i][2*j+0], *nsig        [i][2*j+0]) );
      nsig_tt    [i][2*j+1] = new RooFormulaVar( Form("nsig_tt_%d_%d",    i,2*j+1), "@0*@1", RooArgSet(*r_tt[i][2*j+1], *nsig        [i][2*j+1]) );
      nsig_tf    [i][2*j+0] = new RooFormulaVar( Form("nsig_tf_%d_%d",    i,2*j+0), "@0*@1", RooArgSet(*r_tf[i][2*j+0], *nsig        [i][2*j+1]) );
      nsig_tf    [i][2*j+1] = new RooFormulaVar( Form("nsig_tf_%d_%d",    i,2*j+1), "@0*@1", RooArgSet(*r_tf[i][2*j+1], *nsig        [i][2*j+0]) );
      nsig_tt_sim[i][2*j+0] = new RooFormulaVar( Form("nsig_tt_sim_%d_%d",i,2*j+0), "@0*@1", RooArgSet(*r_tt[i][2*j+0], *nsig_sim    [i][2*j+0]) );
      nsig_tt_sim[i][2*j+1] = new RooFormulaVar( Form("nsig_tt_sim_%d_%d",i,2*j+1), "@0*@1", RooArgSet(*r_tt[i][2*j+1], *nsig_sim    [i][2*j+1]) );
      nsig_tf_sim[i][2*j+0] = new RooFormulaVar( Form("nsig_tf_sim_%d_%d",i,2*j+0), "@0*@1", RooArgSet(*r_tf[i][2*j+0], *nsig_sim    [i][2*j+1]) );
      nsig_tf_sim[i][2*j+1] = new RooFormulaVar( Form("nsig_tf_sim_%d_%d",i,2*j+1), "@0*@1", RooArgSet(*r_tf[i][2*j+1], *nsig_sim    [i][2*j+0]) );
    }
  }

  for( Int_t i=0; i<Nroohist; i++ ){
    nbkg    [i] = new RooRealVar*[Nbin_afb];
    nbkg_sim[i] = new RooRealVar*[Nbin_afb];
    for( Int_t j=0; j<Nbin_afb; j++ ){
      Double_t tmp_nbkg = 300;
      nbkg    [i][j] = new RooRealVar ( Form("N_{bkg}%d_%d",      i,j), Form("N_{bkg}%d_%d",      i,j), tmp_nbkg, tmp_nbkg-300, tmp_nbkg+450 );
      nbkg_sim[i][j] = new RooRealVar ( Form("N_{bkg}^{sim}%d_%d",i,j), Form("N_{bkg}^{sim}%d_%d",i,j), tmp_nbkg, tmp_nbkg-300, tmp_nbkg+450 );
    }
  }
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // READ PDF (charmonium)
  //TFile file_cc("pdf/Mbc_peak_cc_s00-5.root"); // gMC
  //TFile file_cc("pdf/Mbc_peak_cc_s00.root"); // CC-MC
  //TFile file_cc("pdf/Mbc_peak_cc2_s00-5.root"); // gMC
  //TFile file_cc("pdf/Mbc_peak_cc2_s00.root"); // CC-MC
  TFile file_cc("pdf_for_1_6q2/Mbc_peak_cc2_s00-5.root"); // gMC
  //TFile file_cc("pdf_for_1_6q2/Mbc_peak_cc2_s00.root"); // CC-MC
  if( file_cc.IsZombie() ) std::cerr << "[ABORT] can not find input-file : " << file_cc.GetName() << std::endl, abort();
  
  TH1D*** hist_cc = new TH1D**[Nroohist];
  for( Int_t i=0; i<Nroohist; i++ ){
    hist_cc[i] = new TH1D*[Nbin_afb+1];
    for( Int_t j=0; j<Nbin_afb+1; j++ ){
      Int_t k = j + i*(Nbin_afb+1);
      hist_cc[i][j] = (TH1D*)file_cc.Get( Form("hist%d",k) );
      if( hist_cc[i][j] == NULL ) std::cerr << "[ABORT] can not find histgram (charmonium)" << k << std::endl, abort();
    }
  }

  RooDataHist*** peak_cc         = new RooDataHist**[Nroohist];
  RooDataHist*** peak_cc_sim     = new RooDataHist**[Nroohist];
  RooHistPdf***  pdf_peak_cc     = new RooHistPdf** [Nroohist];
  RooHistPdf***  pdf_peak_cc_sim = new RooHistPdf** [Nroohist];
  
  for( Int_t i=0; i<Nroohist; i++ ){
    peak_cc        [i] = new RooDataHist*[Nbin_afb+1];
    peak_cc_sim    [i] = new RooDataHist*[Nbin_afb+1];
    pdf_peak_cc    [i] = new RooHistPdf* [Nbin_afb+1];
    pdf_peak_cc_sim[i] = new RooHistPdf* [Nbin_afb+1];
    for( Int_t j=0; j<Nbin_afb+1; j++ ){
      Int_t k = j + i*(Nbin_afb+1);
      peak_cc        [i][j] = new RooDataHist( Form("peak_cc%d",           k), Form("peak_cc%d",       k), *obs[i],  hist_cc     [i][j] );
      peak_cc_sim    [i][j] = new RooDataHist( Form("peak_cc^{sim}%d",     k), Form("peak_cc^{sim}%d", k), *obs_sim, hist_cc     [i][j] );
      pdf_peak_cc    [i][j] = new RooHistPdf ( Form("pdf_peak_cc%d",       k), Form("peak_cc%d",       k), *obs[i],  *peak_cc    [i][j] );
      pdf_peak_cc_sim[i][j] = new RooHistPdf ( Form("pdf_peak_cc^{sim}%d", k), Form("peak_cc^{sim}%d", k), *obs_sim, *peak_cc_sim[i][j] );
    }
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // SETTING of COEFFICIENT (charmonium)
  RooRealVar*** r_peak_cc = new RooRealVar**[Nroohist];
  Double_t scale_cc[3] = {0.823, 0.741, 0.779 };
  //scale_cc[0] += scale_cc[0];
  //scale_cc[1] += scale_cc[1];
  //scale_cc[2] += scale_cc[2];
  for( Int_t i=0; i<Nroohist; i++ ) r_peak_cc[i] = new RooRealVar*[Nbin_afb+1];
  if( Nbin_afb == 14 ){
  }else if( Nbin_afb == 8 ){
    // CC-MC for 1-6 q2 measurement
    /*
    r_peak_cc[ 0][ 0] = new RooRealVar( "r_peak_cc0", "r_peak_cc0", 0*scale_cc[0] );
    r_peak_cc[ 0][ 1] = new RooRealVar( "r_peak_cc1", "r_peak_cc1", 0*scale_cc[0] );
    r_peak_cc[ 0][ 2] = new RooRealVar( "r_peak_cc2", "r_peak_cc2", 0.85*scale_cc[0] );
    r_peak_cc[ 0][ 3] = new RooRealVar( "r_peak_cc3", "r_peak_cc3", 0.97*scale_cc[0] );
    r_peak_cc[ 0][ 4] = new RooRealVar( "r_peak_cc4", "r_peak_cc4", 0.71*scale_cc[0] );
    r_peak_cc[ 0][ 5] = new RooRealVar( "r_peak_cc5", "r_peak_cc5", 0.61*scale_cc[0] );
    r_peak_cc[ 0][ 6] = new RooRealVar( "r_peak_cc6", "r_peak_cc6", 0*scale_cc[0] );
    r_peak_cc[ 0][ 7] = new RooRealVar( "r_peak_cc7", "r_peak_cc7", 0*scale_cc[0] );
    r_peak_cc[ 0][ 8] = new RooRealVar( "r_peak_cc8", "r_peak_cc8", 3.14*scale_cc[0] );
    r_peak_cc[ 1][ 0] = new RooRealVar( "r_peak_cc9", "r_peak_cc9", 0*scale_cc[1] );
    r_peak_cc[ 1][ 1] = new RooRealVar( "r_peak_cc10", "r_peak_cc10", 0.01*scale_cc[1] );
    r_peak_cc[ 1][ 2] = new RooRealVar( "r_peak_cc11", "r_peak_cc11", 0.42*scale_cc[1] );
    r_peak_cc[ 1][ 3] = new RooRealVar( "r_peak_cc12", "r_peak_cc12", 0.39*scale_cc[1] );
    r_peak_cc[ 1][ 4] = new RooRealVar( "r_peak_cc13", "r_peak_cc13", 0.45*scale_cc[1] );
    r_peak_cc[ 1][ 5] = new RooRealVar( "r_peak_cc14", "r_peak_cc14", 0.44*scale_cc[1] );
    r_peak_cc[ 1][ 6] = new RooRealVar( "r_peak_cc15", "r_peak_cc15", 0.01*scale_cc[1] );
    r_peak_cc[ 1][ 7] = new RooRealVar( "r_peak_cc16", "r_peak_cc16", 0.01*scale_cc[1] );
    r_peak_cc[ 1][ 8] = new RooRealVar( "r_peak_cc17", "r_peak_cc17", 1.73*scale_cc[1] );
    r_peak_cc[ 2][ 0] = new RooRealVar( "r_peak_cc18", "r_peak_cc18", 0*scale_cc[2] );
    r_peak_cc[ 2][ 1] = new RooRealVar( "r_peak_cc19", "r_peak_cc19", 0.01*scale_cc[2] );
    r_peak_cc[ 2][ 2] = new RooRealVar( "r_peak_cc20", "r_peak_cc20", 1.27*scale_cc[2] );
    r_peak_cc[ 2][ 3] = new RooRealVar( "r_peak_cc21", "r_peak_cc21", 1.36*scale_cc[2] );
    r_peak_cc[ 2][ 4] = new RooRealVar( "r_peak_cc22", "r_peak_cc22", 1.16*scale_cc[2] );
    r_peak_cc[ 2][ 5] = new RooRealVar( "r_peak_cc23", "r_peak_cc23", 1.05*scale_cc[2] );
    r_peak_cc[ 2][ 6] = new RooRealVar( "r_peak_cc24", "r_peak_cc24", 0.01*scale_cc[2] );
    r_peak_cc[ 2][ 7] = new RooRealVar( "r_peak_cc25", "r_peak_cc25", 0.01*scale_cc[2] );
    r_peak_cc[ 2][ 8] = new RooRealVar( "r_peak_cc26", "r_peak_cc26", 4.87*scale_cc[2] );
    */
    // gMC for 1-6 q2 measurement
    ///*
      r_peak_cc[ 0][ 0] = new RooRealVar( "r_peak_cc0", "r_peak_cc0", 0*scale_cc[0] );
      r_peak_cc[ 0][ 1] = new RooRealVar( "r_peak_cc1", "r_peak_cc1", 0.1667*scale_cc[0] );
      r_peak_cc[ 0][ 2] = new RooRealVar( "r_peak_cc2", "r_peak_cc2", 3*scale_cc[0] );
      r_peak_cc[ 0][ 3] = new RooRealVar( "r_peak_cc3", "r_peak_cc3", 3.5*scale_cc[0] );
      r_peak_cc[ 0][ 4] = new RooRealVar( "r_peak_cc4", "r_peak_cc4", 0.8333*scale_cc[0] );
      r_peak_cc[ 0][ 5] = new RooRealVar( "r_peak_cc5", "r_peak_cc5", 0.6667*scale_cc[0] );
      r_peak_cc[ 0][ 6] = new RooRealVar( "r_peak_cc6", "r_peak_cc6", 5.667*scale_cc[0] );
      r_peak_cc[ 0][ 7] = new RooRealVar( "r_peak_cc7", "r_peak_cc7", 6.333*scale_cc[0] );
      r_peak_cc[ 0][ 8] = new RooRealVar( "r_peak_cc8", "r_peak_cc8", 20.17*scale_cc[0] );
      r_peak_cc[ 1][ 0] = new RooRealVar( "r_peak_cc9", "r_peak_cc9", 0*scale_cc[1] );
      r_peak_cc[ 1][ 1] = new RooRealVar( "r_peak_cc10", "r_peak_cc10", 0*scale_cc[1] );
      r_peak_cc[ 1][ 2] = new RooRealVar( "r_peak_cc11", "r_peak_cc11", 0.3333*scale_cc[1] );
      r_peak_cc[ 1][ 3] = new RooRealVar( "r_peak_cc12", "r_peak_cc12", 0.6667*scale_cc[1] );
      r_peak_cc[ 1][ 4] = new RooRealVar( "r_peak_cc13", "r_peak_cc13", 0*scale_cc[1] );
      r_peak_cc[ 1][ 5] = new RooRealVar( "r_peak_cc14", "r_peak_cc14", 1*scale_cc[1] );
      r_peak_cc[ 1][ 6] = new RooRealVar( "r_peak_cc15", "r_peak_cc15", 7.333*scale_cc[1] );
      r_peak_cc[ 1][ 7] = new RooRealVar( "r_peak_cc16", "r_peak_cc16", 7.167*scale_cc[1] );
      r_peak_cc[ 1][ 8] = new RooRealVar( "r_peak_cc17", "r_peak_cc17", 16.5*scale_cc[1] );
      r_peak_cc[ 2][ 0] = new RooRealVar( "r_peak_cc18", "r_peak_cc18", 0*scale_cc[2] );
      r_peak_cc[ 2][ 1] = new RooRealVar( "r_peak_cc19", "r_peak_cc19", 0.1667*scale_cc[2] );
      r_peak_cc[ 2][ 2] = new RooRealVar( "r_peak_cc20", "r_peak_cc20", 3.333*scale_cc[2] );
      r_peak_cc[ 2][ 3] = new RooRealVar( "r_peak_cc21", "r_peak_cc21", 4.167*scale_cc[2] );
      r_peak_cc[ 2][ 4] = new RooRealVar( "r_peak_cc22", "r_peak_cc22", 0.8333*scale_cc[2] );
      r_peak_cc[ 2][ 5] = new RooRealVar( "r_peak_cc23", "r_peak_cc23", 1.667*scale_cc[2] );
      r_peak_cc[ 2][ 6] = new RooRealVar( "r_peak_cc24", "r_peak_cc24", 13*scale_cc[2] );
      r_peak_cc[ 2][ 7] = new RooRealVar( "r_peak_cc25", "r_peak_cc25", 13.5*scale_cc[2] );
      r_peak_cc[ 2][ 8] = new RooRealVar( "r_peak_cc26", "r_peak_cc26", 36.67*scale_cc[2] );
      //*/
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // READ PDF (double miss-PID)
  //TFile file_double("pdf/Mbc_peak_double_s00.root"); // MC
  //TFile file_double("pdf/Mbc_peak_double2_s00.root"); // MC
  //TFile file_double("pdf_data/Mbc_peak_double_s00.root"); // RD
  //TFile file_double("pdf_data/Mbc_peak_double2_s00.root"); // RD
  TFile file_double("pdf_for_1_6q2/Mbc_peak_double2_s00.root"); // RD
  if( file_double.IsZombie() ) std::cerr << "[ABORT] can not find input-file : " << file_double.GetName() << std::endl, abort();
  
  TH1D*** hist_double = new TH1D**[Nroohist];
  for( Int_t i=0; i<Nroohist; i++ ){
    hist_double[i] = new TH1D*[Nbin_afb+1];
    for( Int_t j=0; j<Nbin_afb+1; j++ ){
      Int_t k = j + i*(Nbin_afb+1);
      hist_double[i][j] = (TH1D*)file_double.Get( Form("hist%d",k) );
      if( hist_double[i][j] == NULL ) std::cerr << "[ABORT] can not find histgram (double)" << k << std::endl, abort();
    }
  }

  RooDataHist*** peak_double         = new RooDataHist**[Nroohist];
  RooDataHist*** peak_double_sim     = new RooDataHist**[Nroohist];
  RooHistPdf***  pdf_peak_double     = new RooHistPdf** [Nroohist];
  RooHistPdf***  pdf_peak_double_sim = new RooHistPdf** [Nroohist];
  
  for( Int_t i=0; i<Nroohist; i++ ){
    peak_double        [i] = new RooDataHist*[Nbin_afb+1];
    peak_double_sim    [i] = new RooDataHist*[Nbin_afb+1];
    pdf_peak_double    [i] = new RooHistPdf* [Nbin_afb+1];
    pdf_peak_double_sim[i] = new RooHistPdf* [Nbin_afb+1];
    for( Int_t j=0; j<Nbin_afb+1; j++ ){
      Int_t k = j + i*(Nbin_afb+1);
      peak_double        [i][j] = new RooDataHist( Form("peak_double%d",           k), Form("peak_double%d",       k), *obs[i],  hist_double     [i][j] );
      peak_double_sim    [i][j] = new RooDataHist( Form("peak_double^{sim}%d",     k), Form("peak_double^{sim}%d", k), *obs_sim, hist_double     [i][j] );
      pdf_peak_double    [i][j] = new RooHistPdf ( Form("pdf_peak_double%d",       k), Form("peak_double%d",       k), *obs[i],  *peak_double    [i][j] );
      pdf_peak_double_sim[i][j] = new RooHistPdf ( Form("pdf_peak_double^{sim}%d", k), Form("peak_double^{sim}%d", k), *obs_sim, *peak_double_sim[i][j] );
    }
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // SETTING of COEFFICIENT (double miss-PID)
  RooRealVar*** r_peak_double = new RooRealVar**[Nroohist];
  Double_t scale_double[3] = {1.0, 1.0, 1.0};
  //scale_double[0] += 0.55*scale_double[0];
  //scale_double[1] += 0.07*scale_double[1];
  //scale_double[2] += 0.55*scale_double[2];
  for( Int_t i=0; i<Nroohist; i++ ) r_peak_double[i] = new RooRealVar*[Nbin_afb+1];
  if( Nbin_afb == 14 ){
  }else if( Nbin_afb == 8 ){
    // MC
    /*
    r_peak_double[ 0][ 0] = new RooRealVar( "r_peak_double0", "r_peak_double0", 0.2983/0.273 );
    r_peak_double[ 0][ 1] = new RooRealVar( "r_peak_double1", "r_peak_double1", 0.3003/0.273 );
    r_peak_double[ 0][ 2] = new RooRealVar( "r_peak_double2", "r_peak_double2", 0.1071/0.273 );
    r_peak_double[ 0][ 3] = new RooRealVar( "r_peak_double3", "r_peak_double3", 0.07031/0.273 );
    r_peak_double[ 0][ 4] = new RooRealVar( "r_peak_double4", "r_peak_double4", 0.05812/0.273 );
    r_peak_double[ 0][ 5] = new RooRealVar( "r_peak_double5", "r_peak_double5", 0.03668/0.273 );
    r_peak_double[ 0][ 6] = new RooRealVar( "r_peak_double6", "r_peak_double6", 0.138/0.273 );
    r_peak_double[ 0][ 7] = new RooRealVar( "r_peak_double7", "r_peak_double7", 0.08922/0.273 );
    r_peak_double[ 0][ 8] = new RooRealVar( "r_peak_double8", "r_peak_double8", 1.098/0.273 );
    r_peak_double[ 1][ 0] = new RooRealVar( "r_peak_double9", "r_peak_double9", 6.632/0.273 );
    r_peak_double[ 1][ 1] = new RooRealVar( "r_peak_double10", "r_peak_double10", 6.358/0.273 );
    r_peak_double[ 1][ 2] = new RooRealVar( "r_peak_double11", "r_peak_double11", 4.3/0.273 );
    r_peak_double[ 1][ 3] = new RooRealVar( "r_peak_double12", "r_peak_double12", 2.403/0.273 );
    r_peak_double[ 1][ 4] = new RooRealVar( "r_peak_double13", "r_peak_double13", 2.103/0.273 );
    r_peak_double[ 1][ 5] = new RooRealVar( "r_peak_double14", "r_peak_double14", 1.289/0.273 );
    r_peak_double[ 1][ 6] = new RooRealVar( "r_peak_double15", "r_peak_double15", 1.599/0.273 );
    r_peak_double[ 1][ 7] = new RooRealVar( "r_peak_double16", "r_peak_double16", 0.9708/0.273 );
    r_peak_double[ 1][ 8] = new RooRealVar( "r_peak_double17", "r_peak_double17", 25.65/0.273 );
    r_peak_double[ 2][ 0] = new RooRealVar( "r_peak_double18", "r_peak_double18", 6.93/0.273 );
    r_peak_double[ 2][ 1] = new RooRealVar( "r_peak_double19", "r_peak_double19", 6.658/0.273 );
    r_peak_double[ 2][ 2] = new RooRealVar( "r_peak_double20", "r_peak_double20", 4.407/0.273 );
    r_peak_double[ 2][ 3] = new RooRealVar( "r_peak_double21", "r_peak_double21", 2.474/0.273 );
    r_peak_double[ 2][ 4] = new RooRealVar( "r_peak_double22", "r_peak_double22", 2.161/0.273 );
    r_peak_double[ 2][ 5] = new RooRealVar( "r_peak_double23", "r_peak_double23", 1.325/0.273 );
    r_peak_double[ 2][ 6] = new RooRealVar( "r_peak_double24", "r_peak_double24", 1.737/0.273 );
    r_peak_double[ 2][ 7] = new RooRealVar( "r_peak_double25", "r_peak_double25", 1.06/0.273 );
    r_peak_double[ 2][ 8] = new RooRealVar( "r_peak_double26", "r_peak_double26", 26.75/0.273 );
    */
    // RD for 1-6 q2 measurement
    ///*
    r_peak_double[ 0][ 0] = new RooRealVar( "r_peak_double0", "r_peak_double0", 0.01693/0.273*scale_double[0] );
    r_peak_double[ 0][ 1] = new RooRealVar( "r_peak_double1", "r_peak_double1", 0.01886/0.273*scale_double[0] );
    r_peak_double[ 0][ 2] = new RooRealVar( "r_peak_double2", "r_peak_double2", 0.04762/0.273*scale_double[0] );
    r_peak_double[ 0][ 3] = new RooRealVar( "r_peak_double3", "r_peak_double3", 0.04082/0.273*scale_double[0] );
    r_peak_double[ 0][ 4] = new RooRealVar( "r_peak_double4", "r_peak_double4", 0.01268/0.273*scale_double[0] );
    r_peak_double[ 0][ 5] = new RooRealVar( "r_peak_double5", "r_peak_double5", 0.007482/0.273*scale_double[0] );
    r_peak_double[ 0][ 6] = new RooRealVar( "r_peak_double6", "r_peak_double6", 0.03242/0.273*scale_double[0] );
    r_peak_double[ 0][ 7] = new RooRealVar( "r_peak_double7", "r_peak_double7", 0.02005/0.273*scale_double[0] );
    r_peak_double[ 0][ 8] = new RooRealVar( "r_peak_double8", "r_peak_double8", 0.1969/0.273*scale_double[0] );
    r_peak_double[ 1][ 0] = new RooRealVar( "r_peak_double9", "r_peak_double9", 1.386/0.273*scale_double[1] );
    r_peak_double[ 1][ 1] = new RooRealVar( "r_peak_double10", "r_peak_double10", 1.415/0.273*scale_double[1] );
    r_peak_double[ 1][ 2] = new RooRealVar( "r_peak_double11", "r_peak_double11", 2.259/0.273*scale_double[1] );
    r_peak_double[ 1][ 3] = new RooRealVar( "r_peak_double12", "r_peak_double12", 1.766/0.273*scale_double[1] );
    r_peak_double[ 1][ 4] = new RooRealVar( "r_peak_double13", "r_peak_double13", 1.12/0.273*scale_double[1] );
    r_peak_double[ 1][ 5] = new RooRealVar( "r_peak_double14", "r_peak_double14", 0.6334/0.273*scale_double[1] );
    r_peak_double[ 1][ 6] = new RooRealVar( "r_peak_double15", "r_peak_double15", 0.863/0.273*scale_double[1] );
    r_peak_double[ 1][ 7] = new RooRealVar( "r_peak_double16", "r_peak_double16", 0.5026/0.273*scale_double[1] );
    r_peak_double[ 1][ 8] = new RooRealVar( "r_peak_double17", "r_peak_double17", 9.945/0.273*scale_double[1] );
    r_peak_double[ 2][ 0] = new RooRealVar( "r_peak_double18", "r_peak_double18", 1.403/0.273*scale_double[2] );
    r_peak_double[ 2][ 1] = new RooRealVar( "r_peak_double19", "r_peak_double19", 1.434/0.273*scale_double[2] );
    r_peak_double[ 2][ 2] = new RooRealVar( "r_peak_double20", "r_peak_double20", 2.306/0.273*scale_double[2] );
    r_peak_double[ 2][ 3] = new RooRealVar( "r_peak_double21", "r_peak_double21", 1.807/0.273*scale_double[2] );
    r_peak_double[ 2][ 4] = new RooRealVar( "r_peak_double22", "r_peak_double22", 1.132/0.273*scale_double[2] );
    r_peak_double[ 2][ 5] = new RooRealVar( "r_peak_double23", "r_peak_double23", 0.6409/0.273*scale_double[2] );
    r_peak_double[ 2][ 6] = new RooRealVar( "r_peak_double24", "r_peak_double24", 0.8954/0.273*scale_double[2] );
    r_peak_double[ 2][ 7] = new RooRealVar( "r_peak_double25", "r_peak_double25", 0.5226/0.273*scale_double[2] );
    r_peak_double[ 2][ 8] = new RooRealVar( "r_peak_double26", "r_peak_double26", 10.14/0.273*scale_double[2] );
    //*/
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // READ PDF (swapped miss-PID)
  //TFile file_swap("pdf/Mbc_peak_swap_s00.root"); // MC
  //TFile file_swap("pdf/Mbc_peak_swap2_s00.root"); // MC
  //TFile file_swap("pdf_data/Mbc_peak_swap_s00.root"); // RD
  //TFile file_swap("pdf_data/Mbc_peak_swap2_s00.root"); // RD
  TFile file_swap("pdf_for_1_6q2/Mbc_peak_swap2_s00.root"); // RD
  if( file_swap.IsZombie() ) std::cerr << "[ABORT] can not find input-file : " << file_swap.GetName() << std::endl, abort();
  
  TH1D*** hist_swap = new TH1D**[Nroohist];
  for( Int_t i=0; i<Nroohist; i++ ){
    hist_swap[i] = new TH1D*[Nbin_afb+1];
    for( Int_t j=0; j<Nbin_afb+1; j++ ){
      Int_t k = j + i*(Nbin_afb+1);
      hist_swap[i][j] = (TH1D*)file_swap.Get( Form("hist%d",k) );
      if( hist_swap[i][j] == NULL ) std::cerr << "[ABORT] can not find histgram (swap)" << k << std::endl, abort();
    }
  }

  RooDataHist*** peak_swap         = new RooDataHist**[Nroohist];
  RooDataHist*** peak_swap_sim     = new RooDataHist**[Nroohist];
  RooHistPdf***  pdf_peak_swap     = new RooHistPdf** [Nroohist];
  RooHistPdf***  pdf_peak_swap_sim = new RooHistPdf** [Nroohist];
  
  for( Int_t i=0; i<Nroohist; i++ ){
    peak_swap        [i] = new RooDataHist*[Nbin_afb+1];
    peak_swap_sim    [i] = new RooDataHist*[Nbin_afb+1];
    pdf_peak_swap    [i] = new RooHistPdf* [Nbin_afb+1];
    pdf_peak_swap_sim[i] = new RooHistPdf* [Nbin_afb+1];
    for( Int_t j=0; j<Nbin_afb+1; j++ ){
      Int_t k = j + i*(Nbin_afb+1);
      peak_swap        [i][j] = new RooDataHist( Form("peak_swap%d",           k), Form("peak_swap%d",           k), *obs[i],  hist_swap     [i][j] );
      peak_swap_sim    [i][j] = new RooDataHist( Form("peak_swap^{sim}%d",     k), Form("peak_swap^{sim}%d",     k), *obs_sim, hist_swap     [i][j] );
      pdf_peak_swap    [i][j] = new RooHistPdf ( Form("pdf_peak_swap%d",       k), Form("pdf_peak_swap%d",       k), *obs[i],  *peak_swap    [i][j] );
      pdf_peak_swap_sim[i][j] = new RooHistPdf ( Form("pdf_peak_swap^{sim}%d", k), Form("pdf_peak_swap^{sim}%d", k), *obs_sim, *peak_swap_sim[i][j] );
    }
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // SETTING of COEFFICIENT (swapped miss-PID)
  RooRealVar*** r_peak_swap = new RooRealVar**[Nroohist];
  Double_t scale_swap[3] = {1.0, 1.0, 1.0};
  //scale_swap[0] += 0.28*scale_swap[0];
  //scale_swap[1] += 0.04*scale_swap[1];
  //scale_swap[2] += 0.28*scale_swap[2];
  for( Int_t i=0; i<Nroohist; i++ ) r_peak_swap[i] = new RooRealVar*[Nbin_afb+1];
  if( Nbin_afb == 14 ){
  }else if( Nbin_afb == 8 ){
    // MC
    /*
    r_peak_swap[ 0][ 0] = new RooRealVar( "r_peak_swap0", "r_peak_swap0", 0.05571 );
    r_peak_swap[ 0][ 1] = new RooRealVar( "r_peak_swap1", "r_peak_swap1", 0.06249 );
    r_peak_swap[ 0][ 2] = new RooRealVar( "r_peak_swap2", "r_peak_swap2", 0.1116 );
    r_peak_swap[ 0][ 3] = new RooRealVar( "r_peak_swap3", "r_peak_swap3", 0.04579 );
    r_peak_swap[ 0][ 4] = new RooRealVar( "r_peak_swap4", "r_peak_swap4", 0.0036 );
    r_peak_swap[ 0][ 5] = new RooRealVar( "r_peak_swap5", "r_peak_swap5", 0.0004204 );
    r_peak_swap[ 0][ 6] = new RooRealVar( "r_peak_swap6", "r_peak_swap6", 0 );
    r_peak_swap[ 0][ 7] = new RooRealVar( "r_peak_swap7", "r_peak_swap7", 0 );
    r_peak_swap[ 0][ 8] = new RooRealVar( "r_peak_swap8", "r_peak_swap8", 0.2796 );
    r_peak_swap[ 1][ 0] = new RooRealVar( "r_peak_swap9", "r_peak_swap9", 1.105 );
    r_peak_swap[ 1][ 1] = new RooRealVar( "r_peak_swap10", "r_peak_swap10", 1.208 );
    r_peak_swap[ 1][ 2] = new RooRealVar( "r_peak_swap11", "r_peak_swap11", 3.947 );
    r_peak_swap[ 1][ 3] = new RooRealVar( "r_peak_swap12", "r_peak_swap12", 1.976 );
    r_peak_swap[ 1][ 4] = new RooRealVar( "r_peak_swap13", "r_peak_swap13", 0.2017 );
    r_peak_swap[ 1][ 5] = new RooRealVar( "r_peak_swap14", "r_peak_swap14", 0.02547 );
    r_peak_swap[ 1][ 6] = new RooRealVar( "r_peak_swap15", "r_peak_swap15", 0 );
    r_peak_swap[ 1][ 7] = new RooRealVar( "r_peak_swap16", "r_peak_swap16", 0 );
    r_peak_swap[ 1][ 8] = new RooRealVar( "r_peak_swap17", "r_peak_swap17", 8.463 );
    r_peak_swap[ 2][ 0] = new RooRealVar( "r_peak_swap18", "r_peak_swap18", 1.161 );
    r_peak_swap[ 2][ 1] = new RooRealVar( "r_peak_swap19", "r_peak_swap19", 1.27 );
    r_peak_swap[ 2][ 2] = new RooRealVar( "r_peak_swap20", "r_peak_swap20", 4.058 );
    r_peak_swap[ 2][ 3] = new RooRealVar( "r_peak_swap21", "r_peak_swap21", 2.022 );
    r_peak_swap[ 2][ 4] = new RooRealVar( "r_peak_swap22", "r_peak_swap22", 0.2053 );
    r_peak_swap[ 2][ 5] = new RooRealVar( "r_peak_swap23", "r_peak_swap23", 0.02589 );
    r_peak_swap[ 2][ 6] = new RooRealVar( "r_peak_swap24", "r_peak_swap24", 0 );
    r_peak_swap[ 2][ 7] = new RooRealVar( "r_peak_swap25", "r_peak_swap25", 0 );
    r_peak_swap[ 2][ 8] = new RooRealVar( "r_peak_swap26", "r_peak_swap26", 8.742 );
    */
    // RD for 1-6 q2 measurement
    ///*
    r_peak_swap[ 0][ 0] = new RooRealVar( "r_peak_swap0", "r_peak_swap0", 0.001771*scale_swap[0] );
    r_peak_swap[ 0][ 1] = new RooRealVar( "r_peak_swap1", "r_peak_swap1", 0.001702*scale_swap[0] );
    r_peak_swap[ 0][ 2] = new RooRealVar( "r_peak_swap2", "r_peak_swap2", 0.04234*scale_swap[0] );
    r_peak_swap[ 0][ 3] = new RooRealVar( "r_peak_swap3", "r_peak_swap3", 0.02299*scale_swap[0] );
    r_peak_swap[ 0][ 4] = new RooRealVar( "r_peak_swap4", "r_peak_swap4", 0.002827*scale_swap[0] );
    r_peak_swap[ 0][ 5] = new RooRealVar( "r_peak_swap5", "r_peak_swap5", 0.0005247*scale_swap[0] );
    r_peak_swap[ 0][ 6] = new RooRealVar( "r_peak_swap6", "r_peak_swap6", 0*scale_swap[0] );
    r_peak_swap[ 0][ 7] = new RooRealVar( "r_peak_swap7", "r_peak_swap7", 0*scale_swap[0] );
    r_peak_swap[ 0][ 8] = new RooRealVar( "r_peak_swap8", "r_peak_swap8", 0.07216*scale_swap[0] );
    r_peak_swap[ 1][ 0] = new RooRealVar( "r_peak_swap9", "r_peak_swap9", 0.05881*scale_swap[1] );
    r_peak_swap[ 1][ 1] = new RooRealVar( "r_peak_swap10", "r_peak_swap10", 0.07246*scale_swap[1] );
    r_peak_swap[ 1][ 2] = new RooRealVar( "r_peak_swap11", "r_peak_swap11", 0.9639*scale_swap[1] );
    r_peak_swap[ 1][ 3] = new RooRealVar( "r_peak_swap12", "r_peak_swap12", 0.9059*scale_swap[1] );
    r_peak_swap[ 1][ 4] = new RooRealVar( "r_peak_swap13", "r_peak_swap13", 0.09681*scale_swap[1] );
    r_peak_swap[ 1][ 5] = new RooRealVar( "r_peak_swap14", "r_peak_swap14", 0.06818*scale_swap[1] );
    r_peak_swap[ 1][ 6] = new RooRealVar( "r_peak_swap15", "r_peak_swap15", 0*scale_swap[1] );
    r_peak_swap[ 1][ 7] = new RooRealVar( "r_peak_swap16", "r_peak_swap16", 0*scale_swap[1] );
    r_peak_swap[ 1][ 8] = new RooRealVar( "r_peak_swap17", "r_peak_swap17", 2.166*scale_swap[1] );
    r_peak_swap[ 2][ 0] = new RooRealVar( "r_peak_swap18", "r_peak_swap18", 0.06059*scale_swap[2] );
    r_peak_swap[ 2][ 1] = new RooRealVar( "r_peak_swap19", "r_peak_swap19", 0.07416*scale_swap[2] );
    r_peak_swap[ 2][ 2] = new RooRealVar( "r_peak_swap20", "r_peak_swap20", 1.006*scale_swap[2] );
    r_peak_swap[ 2][ 3] = new RooRealVar( "r_peak_swap21", "r_peak_swap21", 0.9289*scale_swap[2] );
    r_peak_swap[ 2][ 4] = new RooRealVar( "r_peak_swap22", "r_peak_swap22", 0.09964*scale_swap[2] );
    r_peak_swap[ 2][ 5] = new RooRealVar( "r_peak_swap23", "r_peak_swap23", 0.0687*scale_swap[2] );
    r_peak_swap[ 2][ 6] = new RooRealVar( "r_peak_swap24", "r_peak_swap24", 0*scale_swap[2] );
    r_peak_swap[ 2][ 7] = new RooRealVar( "r_peak_swap25", "r_peak_swap25", 0*scale_swap[2] );
    r_peak_swap[ 2][ 8] = new RooRealVar( "r_peak_swap26", "r_peak_swap26", 2.238*scale_swap[2] );
    //*/
  }
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Total-PDF
  RooAddPdf***      pdf            = new RooAddPdf**     [Nroohist]; // [Nbin_afb];
  RooAddPdf***      pdf_sim        = new RooAddPdf**     [Nroohist]; // [Nbin_afb];
  RooSimultaneous** tot_pdf        = new RooSimultaneous*[Nroohist];
  RooFitResult**    fit_result     = new RooFitResult*   [Nroohist];
  RooSimultaneous*  tot_pdf_sim    = new RooSimultaneous( "tot_pdf_sim", "tot_pdf_sim", *tag_sim );

  for( Int_t i=0; i<Nroohist; i++ ){
    tot_pdf[i] = new RooSimultaneous( Form("tot_pdf%d",i), "tot_pdf", *tag[i] );
    pdf    [i] = new RooAddPdf* [Nbin_afb];
    pdf_sim[i] = new RooAddPdf* [Nbin_afb];
    for( Int_t j=0; j<Nbin_afb; j++ ){
      Int_t k = j + i*Nbin_afb;
      if( fl_q2bin!=3 ){
	pdf    [i][j] = new RooAddPdf( Form("pdf%d",    k), Form("pdf%d",    k), RooArgList(*gauss    [i], *modargus    [i][j], *pdf_scf_tt    [i][j], *pdf_scf_tf    [i][j], *pdf_peak_cc    [i][j], *pdf_peak_double    [i][j], *pdf_peak_swap    [i][j]), RooArgList(*nsig     [i][j], *nbkg    [i][j], *nsig_tt    [i][j], *nsig_tf    [i][j], *r_peak_cc[i][j], *r_peak_double[i][j], *r_peak_swap[i][j] ) ); // signal + scf + peak
	pdf_sim[i][j] = new RooAddPdf( Form("pdf_sim%d",k), Form("pdf_sim%d",k), RooArgList(*gauss_sim[i], *modargus_sim[i][j], *pdf_scf_tt_sim[i][j], *pdf_scf_tf_sim[i][j], *pdf_peak_cc_sim[i][j], *pdf_peak_double_sim[i][j], *pdf_peak_swap_sim[i][j]), RooArgList(*nsig_sim [i][j], *nbkg_sim[i][j], *nsig_tt_sim[i][j], *nsig_tf_sim[i][j], *r_peak_cc[i][j], *r_peak_double[i][j], *r_peak_swap[i][j] ) ); // signal + scf + peak
      }else{
	pdf    [i][j] = new RooAddPdf( Form("pdf%d",    k), Form("pdf%d",    k), RooArgList(*gauss    [i], *modargus    [i][j], *pdf_scf_tt    [i][j], *pdf_scf_tf    [i][j], *pdf_peak_double    [i][j]), RooArgList(*nsig     [i][j], *nbkg    [i][j], *nsig_tt    [i][j], *nsig_tf    [i][j], *r_peak_double[i][j]) ); // signal + scf + peak [for 6th q2 bin]
	pdf_sim[i][j] = new RooAddPdf( Form("pdf_sim%d",k), Form("pdf_sim%d",k), RooArgList(*gauss_sim[i], *modargus_sim[i][j], *pdf_scf_tt_sim[i][j], *pdf_scf_tf_sim[i][j], *pdf_peak_double_sim[i][j]), RooArgList(*nsig_sim [i][j], *nbkg_sim[i][j], *nsig_tt_sim[i][j], *nsig_tf_sim[i][j], *r_peak_double[i][j]) ); // signal + scf + peak [for 6th q2 bin]
      }
    }
  }


  // FIT (ee, mm, ee+mm)
  for( Int_t i=0; i<Nroohist; i++ ){
    if( Nbin_afb == 14 ){
      if( fl_q2bin == 0 ) tot_pdf[i]->addPdf( *pdf[i][ 0],  "1st_p" );
      if( fl_q2bin == 0 ) tot_pdf[i]->addPdf( *pdf[i][ 1],  "1st_m" );
      if( fl_q2bin == 1 ) tot_pdf[i]->addPdf( *pdf[i][ 2],  "2nd_p" );
      if( fl_q2bin == 1 ) tot_pdf[i]->addPdf( *pdf[i][ 3],  "2nd_m" );
      if( fl_q2bin == 2 ) tot_pdf[i]->addPdf( *pdf[i][ 4],  "3rd_p" );
      if( fl_q2bin == 2 ) tot_pdf[i]->addPdf( *pdf[i][ 5],  "3rd_m" );
      if( fl_q2bin == 3 ) tot_pdf[i]->addPdf( *pdf[i][ 6],  "5th_p" );
      if( fl_q2bin == 3 ) tot_pdf[i]->addPdf( *pdf[i][ 7],  "5th_m" );
      if( fl_q2bin == 4 ) tot_pdf[i]->addPdf( *pdf[i][ 8],  "7th_p" );
      if( fl_q2bin == 4 ) tot_pdf[i]->addPdf( *pdf[i][ 9],  "7th_m" );
      if( fl_q2bin == 5 ) tot_pdf[i]->addPdf( *pdf[i][10],  "8th_p" );
      if( fl_q2bin == 5 ) tot_pdf[i]->addPdf( *pdf[i][11],  "8th_m" );
      if( fl_q2bin == 6 ) tot_pdf[i]->addPdf( *pdf[i][12],  "9th_p" );
      if( fl_q2bin == 6 ) tot_pdf[i]->addPdf( *pdf[i][13],  "9th_m" );
    }else if( Nbin_afb == 8 ){
      if( fl_q2bin == 0 ) tot_pdf[i]->addPdf( *pdf[i][ 0],  "12_p"  );
      if( fl_q2bin == 0 ) tot_pdf[i]->addPdf( *pdf[i][ 1],  "12_m"  );
      if( fl_q2bin == 1 ) tot_pdf[i]->addPdf( *pdf[i][ 2],  "3_p"   );
      if( fl_q2bin == 1 ) tot_pdf[i]->addPdf( *pdf[i][ 3],  "3_m"   );
      if( fl_q2bin == 2 ) tot_pdf[i]->addPdf( *pdf[i][ 4],  "5_p"   );
      if( fl_q2bin == 2 ) tot_pdf[i]->addPdf( *pdf[i][ 5],  "5_m"   );
      if( fl_q2bin == 3 ) tot_pdf[i]->addPdf( *pdf[i][ 6],  "789_p" );
      if( fl_q2bin == 3 ) tot_pdf[i]->addPdf( *pdf[i][ 7],  "789_m" );
    }
    fit_result[i] = tot_pdf[i]->fitTo( *data[i], Extended(), Save(true) );
  }

  // SIMULTANEOUS FIT( ee and mm )
  if( Nbin_afb == 14 ){
    if( fl_q2bin == 0 ) tot_pdf_sim->addPdf( *pdf_sim[0][ 0],  "1st_p_e"  );
    if( fl_q2bin == 0 ) tot_pdf_sim->addPdf( *pdf_sim[0][ 1],  "1st_m_e"  );
    if( fl_q2bin == 1 ) tot_pdf_sim->addPdf( *pdf_sim[0][ 2],  "2nd_p_e"  );
    if( fl_q2bin == 1 ) tot_pdf_sim->addPdf( *pdf_sim[0][ 3],  "2nd_m_e"  );
    if( fl_q2bin == 2 ) tot_pdf_sim->addPdf( *pdf_sim[0][ 4],  "3rd_p_e"  );
    if( fl_q2bin == 2 ) tot_pdf_sim->addPdf( *pdf_sim[0][ 5],  "3rd_m_e"  );
    if( fl_q2bin == 3 ) tot_pdf_sim->addPdf( *pdf_sim[0][ 6],  "5th_p_e"  );
    if( fl_q2bin == 3 ) tot_pdf_sim->addPdf( *pdf_sim[0][ 7],  "5th_m_e"  );
    if( fl_q2bin == 4 ) tot_pdf_sim->addPdf( *pdf_sim[0][ 8],  "7th_p_e"  );
    if( fl_q2bin == 4 ) tot_pdf_sim->addPdf( *pdf_sim[0][ 9],  "7th_m_e"  );
    if( fl_q2bin == 5 ) tot_pdf_sim->addPdf( *pdf_sim[0][10],  "8th_p_e"  );
    if( fl_q2bin == 5 ) tot_pdf_sim->addPdf( *pdf_sim[0][11],  "8th_m_e"  );
    if( fl_q2bin == 6 ) tot_pdf_sim->addPdf( *pdf_sim[0][12],  "9th_p_e"  );
    if( fl_q2bin == 6 ) tot_pdf_sim->addPdf( *pdf_sim[0][13],  "9th_m_e"  );
    if( fl_q2bin == 0 ) tot_pdf_sim->addPdf( *pdf_sim[1][ 0],  "1st_p_mu" );
    if( fl_q2bin == 0 ) tot_pdf_sim->addPdf( *pdf_sim[1][ 1],  "1st_m_mu" );
    if( fl_q2bin == 1 ) tot_pdf_sim->addPdf( *pdf_sim[1][ 2],  "2nd_p_mu" );
    if( fl_q2bin == 1 ) tot_pdf_sim->addPdf( *pdf_sim[1][ 3],  "2nd_m_mu" );
    if( fl_q2bin == 2 ) tot_pdf_sim->addPdf( *pdf_sim[1][ 4],  "3rd_p_mu" );
    if( fl_q2bin == 2 ) tot_pdf_sim->addPdf( *pdf_sim[1][ 5],  "3rd_m_mu" );
    if( fl_q2bin == 3 ) tot_pdf_sim->addPdf( *pdf_sim[1][ 6],  "5th_p_mu" );
    if( fl_q2bin == 3 ) tot_pdf_sim->addPdf( *pdf_sim[1][ 7],  "5th_m_mu" );
    if( fl_q2bin == 4 ) tot_pdf_sim->addPdf( *pdf_sim[1][ 8],  "7th_p_mu" );
    if( fl_q2bin == 4 ) tot_pdf_sim->addPdf( *pdf_sim[1][ 9],  "7th_m_mu" );
    if( fl_q2bin == 5 ) tot_pdf_sim->addPdf( *pdf_sim[1][10],  "8th_p_mu" );
    if( fl_q2bin == 5 ) tot_pdf_sim->addPdf( *pdf_sim[1][11],  "8th_m_mu" );
    if( fl_q2bin == 6 ) tot_pdf_sim->addPdf( *pdf_sim[1][12],  "9th_p_mu" );
    if( fl_q2bin == 6 ) tot_pdf_sim->addPdf( *pdf_sim[1][13],  "9th_m_mu" );
  }else if( Nbin_afb == 8 ){
    if( fl_q2bin == 0 ) tot_pdf_sim->addPdf( *pdf_sim[0][ 0],  "12_p_e"   );
    if( fl_q2bin == 0 ) tot_pdf_sim->addPdf( *pdf_sim[0][ 1],  "12_m_e"   );
    if( fl_q2bin == 1 ) tot_pdf_sim->addPdf( *pdf_sim[0][ 2],  "3_p_e"    );
    if( fl_q2bin == 1 ) tot_pdf_sim->addPdf( *pdf_sim[0][ 3],  "3_m_e"    );
    if( fl_q2bin == 2 ) tot_pdf_sim->addPdf( *pdf_sim[0][ 4],  "5_p_e"    );
    if( fl_q2bin == 2 ) tot_pdf_sim->addPdf( *pdf_sim[0][ 5],  "5_m_e"    );
    if( fl_q2bin == 3 ) tot_pdf_sim->addPdf( *pdf_sim[0][ 6],  "789_p_e"  );
    if( fl_q2bin == 3 ) tot_pdf_sim->addPdf( *pdf_sim[0][ 7],  "789_m_e"  );
    if( fl_q2bin == 0 ) tot_pdf_sim->addPdf( *pdf_sim[1][ 0],  "12_p_mu"  );
    if( fl_q2bin == 0 ) tot_pdf_sim->addPdf( *pdf_sim[1][ 1],  "12_m_mu"  );
    if( fl_q2bin == 1 ) tot_pdf_sim->addPdf( *pdf_sim[1][ 2],  "3_p_mu"   );
    if( fl_q2bin == 1 ) tot_pdf_sim->addPdf( *pdf_sim[1][ 3],  "3_m_mu"   );
    if( fl_q2bin == 2 ) tot_pdf_sim->addPdf( *pdf_sim[1][ 4],  "5_p_mu"   );
    if( fl_q2bin == 2 ) tot_pdf_sim->addPdf( *pdf_sim[1][ 5],  "5_m_mu"   );
    if( fl_q2bin == 3 ) tot_pdf_sim->addPdf( *pdf_sim[1][ 6],  "789_p_mu" );
    if( fl_q2bin == 3 ) tot_pdf_sim->addPdf( *pdf_sim[1][ 7],  "789_m_mu" );
  }
  RooFitResult* fit_result_sim = tot_pdf_sim->fitTo( *data_sim, Extended(), Save(true) );
  
  // +++++++++++++++++++++++++ Draw(RooPlot) +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  RooPlot**    frame_param = new RooPlot*   [Nroohist];
  RooPlot***   frame       = new RooPlot**  [Nroohist]; // [Nbin_afb];
  TPaveText*** box         = new TPaveText**[Nroohist]; // [Nbin_afb];

  for( Int_t i=0; i<Nroohist; i++ ){
    frame[i] = new RooPlot*  [Nbin_afb];
    box  [i] = new TPaveText*[Nbin_afb];
    for( Int_t j=2*fl_q2bin; j<2*(fl_q2bin+1); j++ ){
      Int_t k = j+i*Nbin_afb;
      frame[i][j] = obs[i]->frame();
      frame[i][j]->GetXaxis()->CenterTitle();
      frame[i][j]->GetYaxis()->CenterTitle();
      frame[i][j]->SetTitleOffset( 1.00,"x" );
      frame[i][j]->SetTitleOffset( 1.30,"y" );
      if     ( i==0 ) frame[i][j]->SetTitle( Form("bin=%d, ee       ",j+1) );
      else if( i==1 ) frame[i][j]->SetTitle( Form("bin=%d, #mu#mu   ",j+1) );
      else if( i==2 ) frame[i][j]->SetTitle( Form("bin=%d, ee+#mu#mu",j+1) );
      if( Nbin_afb == 14 ){
	if     ( j== 0 ) data[i]->plotOn ( frame[i][j], Binning(xbin), LineWidth(1), Cut( Form("bin%d==bin%d::1st_p",i,i) ) );
	else if( j== 1 ) data[i]->plotOn ( frame[i][j], Binning(xbin), LineWidth(1), Cut( Form("bin%d==bin%d::1st_m",i,i) ) );
	else if( j== 2 ) data[i]->plotOn ( frame[i][j], Binning(xbin), LineWidth(1), Cut( Form("bin%d==bin%d::2nd_p",i,i) ) );
	else if( j== 3 ) data[i]->plotOn ( frame[i][j], Binning(xbin), LineWidth(1), Cut( Form("bin%d==bin%d::2nd_m",i,i) ) );
	else if( j== 4 ) data[i]->plotOn ( frame[i][j], Binning(xbin), LineWidth(1), Cut( Form("bin%d==bin%d::3rd_p",i,i) ) );
	else if( j== 5 ) data[i]->plotOn ( frame[i][j], Binning(xbin), LineWidth(1), Cut( Form("bin%d==bin%d::3rd_m",i,i) ) );
	else if( j== 6 ) data[i]->plotOn ( frame[i][j], Binning(xbin), LineWidth(1), Cut( Form("bin%d==bin%d::5th_p",i,i) ) );
	else if( j== 7 ) data[i]->plotOn ( frame[i][j], Binning(xbin), LineWidth(1), Cut( Form("bin%d==bin%d::5th_m",i,i) ) );
	else if( j== 8 ) data[i]->plotOn ( frame[i][j], Binning(xbin), LineWidth(1), Cut( Form("bin%d==bin%d::7th_p",i,i) ) );
	else if( j== 9 ) data[i]->plotOn ( frame[i][j], Binning(xbin), LineWidth(1), Cut( Form("bin%d==bin%d::7th_m",i,i) ) );
	else if( j==10 ) data[i]->plotOn ( frame[i][j], Binning(xbin), LineWidth(1), Cut( Form("bin%d==bin%d::8th_p",i,i) ) );
	else if( j==11 ) data[i]->plotOn ( frame[i][j], Binning(xbin), LineWidth(1), Cut( Form("bin%d==bin%d::8th_m",i,i) ) );
	else if( j==12 ) data[i]->plotOn ( frame[i][j], Binning(xbin), LineWidth(1), Cut( Form("bin%d==bin%d::9th_p",i,i) ) );
	else if( j==13 ) data[i]->plotOn ( frame[i][j], Binning(xbin), LineWidth(1), Cut( Form("bin%d==bin%d::9th_m",i,i) ) );
      }else if( Nbin_afb == 8 ){
	if     ( j== 0 ) data[i]->plotOn ( frame[i][j], Binning(xbin), LineWidth(1), Cut( Form("bin%d==bin%d::12_p" ,i,i) ) );
	else if( j== 1 ) data[i]->plotOn ( frame[i][j], Binning(xbin), LineWidth(1), Cut( Form("bin%d==bin%d::12_m" ,i,i) ) );
	else if( j== 2 ) data[i]->plotOn ( frame[i][j], Binning(xbin), LineWidth(1), Cut( Form("bin%d==bin%d::3_p"  ,i,i) ) );
	else if( j== 3 ) data[i]->plotOn ( frame[i][j], Binning(xbin), LineWidth(1), Cut( Form("bin%d==bin%d::3_m"  ,i,i) ) );
	else if( j== 4 ) data[i]->plotOn ( frame[i][j], Binning(xbin), LineWidth(1), Cut( Form("bin%d==bin%d::5_p"  ,i,i) ) );
	else if( j== 5 ) data[i]->plotOn ( frame[i][j], Binning(xbin), LineWidth(1), Cut( Form("bin%d==bin%d::5_m"  ,i,i) ) );
	else if( j== 6 ) data[i]->plotOn ( frame[i][j], Binning(xbin), LineWidth(1), Cut( Form("bin%d==bin%d::789_p",i,i) ) );
	else if( j== 7 ) data[i]->plotOn ( frame[i][j], Binning(xbin), LineWidth(1), Cut( Form("bin%d==bin%d::789_m",i,i) ) );
      }
      tag[i]->setIndex(j+1); tot_pdf[i] ->plotOn ( frame[i][j], Slice(*tag[i]), ProjWData(*tag[i],*data[i]), LineWidth(1) );
      box[i][j] = new TPaveText( 0.05, 0.87, 0.35, 0.93,"BRNDC" );
      box[i][j]->AddText( Form("#chi^{2}/ndf = %f", frame[i][j]->chiSquare(nparam+N0bin[j][rooad[i]])) );
      tag[i]->setIndex(j+1); tot_pdf[i]->plotOn ( frame[i][j], Slice(*tag[i]), ProjWData(*tag[i],*data[i]), Components( *gauss   [i]                                                                ), LineStyle(7), LineColor(2), LineWidth(1) ); // sig (gauss)
      tag[i]->setIndex(j+1); tot_pdf[i]->plotOn ( frame[i][j], Slice(*tag[i]), ProjWData(*tag[i],*data[i]), Components( *modargus[i][j]                                                             ), LineStyle(7), LineColor(4), LineWidth(1) ); // bkg (arguss)
      tag[i]->setIndex(j+1); tot_pdf[i]->plotOn ( frame[i][j], Slice(*tag[i]), ProjWData(*tag[i],*data[i]), Components( RooArgSet(*pdf_scf_tt [i][j], *pdf_scf_tf[i][j])                            ), LineStyle(7), LineColor(5), LineWidth(1) ); // self cros-feed
      tag[i]->setIndex(j+1); tot_pdf[i]->plotOn ( frame[i][j], Slice(*tag[i]), ProjWData(*tag[i],*data[i]), Components( RooArgSet(*pdf_peak_cc[i][j], *pdf_peak_double[i][j], *pdf_peak_swap[i][j]) ), LineStyle(7), LineColor(6), LineWidth(1) ); // peak
      frame[i][j]->addObject( box[i][j] );
    }
    frame_param[i] = obs[i]->frame();
    tot_pdf[i] ->paramOn( frame_param[i], Format("NELU", AutoPrecision(2)), Layout(0.40, 0.99, 0.99), ShowConstants(kFALSE) );
    frame_param[i]->getAttText()->SetTextSize(0.045);
  }

  // +++++++++++++++++++++++++ Draw(RooPlot, simultaneous fit) +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  RooPlot**    frame_param_sim = new RooPlot*   [2];
  RooPlot***   frame_sim       = new RooPlot**  [2]; // [Nbin_afb];
  TPaveText*** box_sim         = new TPaveText**[2]; // [Nbin_afb];

  for( Int_t i=0; i<2; i++ ){
    frame_sim[i] = new RooPlot*  [Nbin_afb];
    box_sim  [i] = new TPaveText*[Nbin_afb];
    for( Int_t j=2*fl_q2bin; j<2*(fl_q2bin+1); j++ ){
      Int_t k = j+i*Nbin_afb;
      frame_sim[i][j] = obs_sim->frame();
      frame_sim[i][j]->GetXaxis()->CenterTitle();
      frame_sim[i][j]->GetYaxis()->CenterTitle();
      frame_sim[i][j]->SetTitleOffset( 1.00,"x" );
      frame_sim[i][j]->SetTitleOffset( 1.30,"y" );
      if     ( i==0 ) frame_sim[i][j]->SetTitle( Form("bin=%d, ee",    j+1) );
      else if( i==1 ) frame_sim[i][j]->SetTitle( Form("bin=%d, #mu#mu",j+1) );
      if( Nbin_afb == 14 ){
	if( i==0 ){
	  if     ( j== 0 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::1st_p_e" ) );
	  else if( j== 1 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::1st_m_e" ) );
	  else if( j== 2 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::2nd_p_e" ) );
	  else if( j== 3 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::2nd_m_e" ) );
	  else if( j== 4 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::3rd_p_e" ) );
	  else if( j== 5 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::3rd_m_e" ) );
	  else if( j== 6 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::5th_p_e" ) );
	  else if( j== 7 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::5th_m_e" ) );
	  else if( j== 8 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::7th_p_e" ) );
	  else if( j== 9 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::7th_m_e" ) );
	  else if( j==10 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::8th_p_e" ) );
	  else if( j==11 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::8th_m_e" ) );
	  else if( j==12 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::9th_p_e" ) );
	  else if( j==13 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::9th_m_e" ) );
	}else{
	  if     ( j== 0 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::1st_p_mu" ) );
	  else if( j== 1 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::1st_m_mu" ) );
	  else if( j== 2 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::2nd_p_mu" ) );
	  else if( j== 3 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::2nd_m_mu" ) );
	  else if( j== 4 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::3rd_p_mu" ) );
	  else if( j== 5 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::3rd_m_mu" ) );
	  else if( j== 6 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::5th_p_mu" ) );
	  else if( j== 7 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::5th_m_mu" ) );
	  else if( j== 8 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::7th_p_mu" ) );
	  else if( j== 9 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::7th_m_mu" ) );
	  else if( j==10 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::8th_p_mu" ) );
	  else if( j==11 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::8th_m_mu" ) );
	  else if( j==12 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::9th_p_mu" ) );
	  else if( j==13 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::9th_m_mu" ) );
	}
      }else if( Nbin_afb == 8 ){
	if( i==0 ){
	  if     ( j== 0 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::12_p_e"  ) );
	  else if( j== 1 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::12_m_e"  ) );
	  else if( j== 2 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::3_p_e"   ) );
	  else if( j== 3 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::3_m_e"   ) );
	  else if( j== 4 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::5_p_e"   ) );
	  else if( j== 5 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::5_m_e"   ) );
	  else if( j== 6 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::789_p_e" ) );
	  else if( j== 7 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::789_m_e" ) );
	}else{
	  if     ( j== 0 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::12_p_mu"  ) );
	  else if( j== 1 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::12_m_mu"  ) );
	  else if( j== 2 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::3_p_mu"   ) );
	  else if( j== 3 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::3_m_mu"   ) );
	  else if( j== 4 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::5_p_mu"   ) );
	  else if( j== 5 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::5_m_mu"   ) );
	  else if( j== 6 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::789_p_mu" ) );
	  else if( j== 7 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::789_m_mu" ) );
	}
      }
      tag_sim->setIndex(i==0 ? 100+j+1 : j+1 ); tot_pdf_sim->plotOn( frame_sim[i][j], Slice(*tag_sim), ProjWData(*tag_sim,*data_sim), LineWidth(1) );
      box_sim[i][j] = new TPaveText( 0.05, 0.87, 0.35, 0.93,"BRNDC" );
      box_sim[i][j]->AddText( Form("#chi^{2}/ndf = %f", frame_sim[i][j]->chiSquare(nparam+N0bin[j][rooad[i]])) );
      tag_sim->setIndex(i==0 ? 100+j+1 : j+1); tot_pdf_sim->plotOn ( frame_sim[i][j], Slice(*tag_sim), ProjWData(*tag_sim,*data_sim), Components( RooArgSet(*gauss_sim[i],*pdf_scf_tt_sim[i][j],*pdf_scf_tf_sim[i][j]) ), LineStyle(7), LineColor(2), LineWidth(1) ); // sig (gauss) + self cross-feed
      //tag_sim->setIndex(i==0 ? 100+j+1 : j+1); tot_pdf_sim->plotOn ( frame_sim[i][j], Slice(*tag_sim), ProjWData(*tag_sim,*data_sim), Components( *gauss_sim   [i]                                                                      ), LineStyle(7), LineColor(2), LineWidth(1) ); // sig (gauss)
      //tag_sim->setIndex(i==0 ? 100+j+1 : j+1); tot_pdf_sim->plotOn ( frame_sim[i][j], Slice(*tag_sim), ProjWData(*tag_sim,*data_sim), Components( RooArgSet(*pdf_scf_tt_sim [i][j],*pdf_scf_tf_sim     [i][j])                          ), LineStyle(7), LineColor(5), LineWidth(1) ); // self cross-feed
      tag_sim->setIndex(i==0 ? 100+j+1 : j+1); tot_pdf_sim->plotOn ( frame_sim[i][j], Slice(*tag_sim), ProjWData(*tag_sim,*data_sim), Components( *modargus_sim[i][j]                                                                   ), LineStyle(7), LineColor(4), LineWidth(1) ); // bkg (arguss)

      if( fl_q2bin!=3 ){
	//tag_sim->setIndex(i==0 ? 100+j+1 : j+1); tot_pdf_sim->plotOn ( frame_sim[i][j], Slice(*tag_sim), ProjWData(*tag_sim,*data_sim), Components( RooArgSet(*pdf_peak_cc_sim[i][j],*pdf_peak_double_sim[i][j],*pdf_peak_swap_sim[i][j]) ), LineStyle(7), LineColor(6), LineWidth(1) ); // peak
	tag_sim->setIndex(i==0 ? 100+j+1 : j+1); tot_pdf_sim->plotOn ( frame_sim[i][j], Slice(*tag_sim), ProjWData(*tag_sim,*data_sim), Components( RooArgSet(*pdf_peak_cc_sim[i][j],*pdf_peak_double_sim[i][j],*pdf_peak_swap_sim[i][j]) ), FillStyle(1001), FillColor(11), LineWidth(1), VLines(), DrawOption("F") ); // peak // Fill style
      }else{
	//tag_sim->setIndex(i==0 ? 100+j+1 : j+1); tot_pdf_sim->plotOn ( frame_sim[i][j], Slice(*tag_sim), ProjWData(*tag_sim,*data_sim), Components( *pdf_peak_double_sim[i][j]                                                            ), LineStyle(7), LineColor(6), LineWidth(1) ); // peak for 6th q2 bin
	tag_sim->setIndex(i==0 ? 100+j+1 : j+1); tot_pdf_sim->plotOn ( frame_sim[i][j], Slice(*tag_sim), ProjWData(*tag_sim,*data_sim), Components( *pdf_peak_double_sim[i][j]                                                            ), FillStyle(1001), FillColor(11), LineWidth(1), VLines(), DrawOption("F") ); // peak for 6th q2 bin // Fill style
      }
      frame_sim[i][j]->addObject( box_sim[i][j] );
    }
    frame_param_sim[i] = obs_sim->frame();
    tot_pdf_sim->paramOn( frame_param_sim[i], Format("NELU", AutoPrecision(2)), Layout(0.40, 0.99, 0.99), ShowConstants(kFALSE) );
    frame_param_sim[i]->getAttText()->SetTextSize(0.045);
  }
  // +++++++++++++++++++++++++ Draw(TCanvas) +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TCanvas* c2 = Canvas( "c2", "c2", 3, 3 );
  ///*
  for( Int_t i=0; i<Nroohist; i++ ){
    frame[i][2*fl_q2bin+0]->SetAxisRange( 5.22, 5.30, "X" );
    frame[i][2*fl_q2bin+1]->SetAxisRange( 5.22, 5.30, "X" );
    if( i==0 ){
      if( fl_q2bin==0 ){ // 1st q2 bin
	frame[i][0]->SetAxisRange( 0, 30, "Y" ); // ee
	frame[i][1]->SetAxisRange( 0, 30, "Y" ); // ee
      }else if( fl_q2bin==1 ){ // 2nd q2 bin
	frame[i][2]->SetAxisRange( 0, 24, "Y" ); // ee
	frame[i][3]->SetAxisRange( 0, 24, "Y" ); // ee
      }else if( fl_q2bin==2 ){ // 4th q2 bin
	frame[i][4]->SetAxisRange( 0, 14, "Y" ); // ee
	frame[i][5]->SetAxisRange( 0, 14, "Y" ); // ee
      }else if( fl_q2bin==3 ){ // 6th q2 bin
	frame[i][6]->SetAxisRange( 0, 23, "Y" ); // ee
	frame[i][7]->SetAxisRange( 0, 23, "Y" ); // ee
      }
    }else if( i==1 ){
      if( fl_q2bin==0 ){ // 1st q2 bin
	frame[i][0]->SetAxisRange( 0, 22, "Y" ); // mm
	frame[i][1]->SetAxisRange( 0, 22, "Y" ); // mm
      }else if( fl_q2bin==1 ){ // 2nd q2 bin
	frame[i][2]->SetAxisRange( 0, 32, "Y" ); // mm
	frame[i][3]->SetAxisRange( 0, 32, "Y" ); // mm
      }else if( fl_q2bin==2 ){ // 4th q2 bin
	frame[i][4]->SetAxisRange( 0, 25, "Y" ); // mm
	frame[i][5]->SetAxisRange( 0, 25, "Y" ); // mm
      }else if( fl_q2bin==3 ){ // 6th q2 bin
	frame[i][6]->SetAxisRange( 0, 29, "Y" ); // mm
	frame[i][7]->SetAxisRange( 0, 29, "Y" ); // mm
      }
    }else if( i==2 ){
      if( fl_q2bin==0 ){ // 1st q2 bin
	frame[i][0]->SetAxisRange( 0,  45, "Y" ); // ll
	frame[i][1]->SetAxisRange( 0,  45, "Y" ); // ll
      }else if( fl_q2bin==1 ){ // 2nd q2 bin
	frame[i][2]->SetAxisRange( 0,  50, "Y" ); // ll
	frame[i][3]->SetAxisRange( 0,  50, "Y" ); // ll
      }else if( fl_q2bin==2 ){ // 4th q2 bin
	frame[i][4]->SetAxisRange( 0,  34, "Y" ); // ll
	frame[i][5]->SetAxisRange( 0,  34, "Y" ); // ll
      }else if( fl_q2bin==3 ){ // 6th q2 bin
	frame[i][6]->SetAxisRange( 0,  50, "Y" ); // ll
	frame[i][7]->SetAxisRange( 0,  50, "Y" ); // ll
      }
    }
  }
  //*/
  
  for( Int_t i=0; i<Nroohist; i++ ){
    c2->cd(3*i+1);
    frame[i][2*fl_q2bin+0]->Draw();
    c2->cd(3*i+2);
    frame[i][2*fl_q2bin+1]->Draw();
    c2->cd(3*i+3);
    frame_param[i]->Draw();
  }

  // +++++++++++++++++++++++++ Draw(TCanvas, simultaneous fit) +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TCanvas* c3 = Canvas( "c3", "c3", 3, 2 );

  for( Int_t i=0; i<2; i++ ){
    frame_sim[i][2*fl_q2bin+0]->SetAxisRange( 5.22, 5.30, "X" );
    frame_sim[i][2*fl_q2bin+1]->SetAxisRange( 5.22, 5.30, "X" );
    if( i==0 ){
      if( fl_q2bin==0 ){ // 1st q2 bin
	frame_sim[i][0]->SetAxisRange( 0, 30, "Y" ); // ee
	frame_sim[i][1]->SetAxisRange( 0, 30, "Y" ); // ee
      }else if( fl_q2bin==1 ){ // 2nd q2 bin
	frame_sim[i][2]->SetAxisRange( 0, 24, "Y" ); // ee
	frame_sim[i][3]->SetAxisRange( 0, 24, "Y" ); // ee
      }else if( fl_q2bin==2 ){ // 4th q2 bin
	frame_sim[i][4]->SetAxisRange( 0, 14, "Y" ); // ee
	frame_sim[i][5]->SetAxisRange( 0, 14, "Y" ); // ee
      }else if( fl_q2bin==3 ){ // 6th q2 bin
	frame_sim[i][6]->SetAxisRange( 0, 23, "Y" ); // ee
	frame_sim[i][7]->SetAxisRange( 0, 23, "Y" ); // ee
      }
    }else if( i==1 ){
      if( fl_q2bin==0 ){ // 1st q2 bin
	frame_sim[i][0]->SetAxisRange( 0, 22, "Y" ); // mm
	frame_sim[i][1]->SetAxisRange( 0, 22, "Y" ); // mm
      }else if( fl_q2bin==1 ){ // 2nd q2 bin
	frame_sim[i][2]->SetAxisRange( 0, 32, "Y" ); // mm
	frame_sim[i][3]->SetAxisRange( 0, 32, "Y" ); // mm
      }else if( fl_q2bin==2 ){ // 4th q2 bin
	frame_sim[i][4]->SetAxisRange( 0, 25, "Y" ); // mm
	frame_sim[i][5]->SetAxisRange( 0, 25, "Y" ); // mm
      }else if( fl_q2bin==3 ){ // 6th q2 bin
	frame_sim[i][6]->SetAxisRange( 0, 29, "Y" ); // mm
	frame_sim[i][7]->SetAxisRange( 0, 29, "Y" ); // mm
      }
    }
  }


  for( Int_t i=0; i<2; i++ ){
    c3->cd(3*i+1);
    //frame_sim[0][2*fl_q2bin+0]->SetMaximum(35);
    //frame_sim[1][2*fl_q2bin+0]->SetMaximum(30);
    frame_sim[i][2*fl_q2bin+0]->Draw();
    c3->cd(3*i+2);
    //frame_sim[0][2*fl_q2bin+1]->SetMaximum(35);
    //frame_sim[1][2*fl_q2bin+1]->SetMaximum(30);
    frame_sim[i][2*fl_q2bin+1]->Draw();
    c3->cd(3*i+3);
    frame_param_sim[i]->Draw();
  }
  
  // +++++++++++++++++++++++++ Likelihood Scan for simultaneous fitting+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TCanvas* c4 = Canvas( "c4","c4", 1, 1 );
  c4->Draw();
  RooPlot* frame_scan = AFB_true[fl_q2bin]->frame( Bins(50),Range(-1.0, 1.0), Title("profile") );
  RooAbsReal* nll     = tot_pdf_sim->createNLL( *data_sim );
  nll->plotOn( frame_scan, LineColor(kBlue), ShiftToZero() );
  RooAbsReal* pnll = nll->createProfile( *AFB_true[fl_q2bin] ) ;
  pnll->plotOn( frame_scan,LineColor(kRed) ) ;
  
  frame_scan->GetXaxis()->CenterTitle();
  frame_scan->GetYaxis()->CenterTitle();
  frame_scan->SetTitleOffset( 1.00, "x" );
  frame_scan->SetTitleOffset( 1.30, "y" );
  frame_scan->SetTitle( "Likelihood Scan" );
  frame_scan->SetMinimum(0);
  frame_scan->SetMaximum(20);
  frame_scan->Draw();

  //+++++++++++++++++++++++++++++++++ LOG +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::cout << std::endl;
  //std::cout << "correlation matrix" << std::endl; fit_result_sim->correlationMatrix().Print();
  //std::cout << "covariance matrix"  << std::endl; fit_result_sim->covarianceMatrix ().Print();

  std::cout << "AFB(true)    = "    << std::setw(8) << std::right << AFB_true[fl_q2bin]->getVal()    << " +- " << std::setw(7) << std::right << AFB_true[fl_q2bin]->getError()                               << std::endl;
  std::cout << "AFB(meas,ee) = "    << std::setw(8) << std::right << AFB_meas[0][fl_q2bin]->getVal() << " +- " << std::setw(7) << std::right << AFB_meas[0][fl_q2bin]->getPropagatedError( *fit_result_sim );
  std::cout << "  [cf_slope(ee) = " << std::setw(7) << std::right << cf_slope[0][fl_q2bin]->getVal() << " ]"    << std::endl;
  std::cout << "AFB(meas,mm) = "    << std::setw(8) << std::right << AFB_meas[1][fl_q2bin]->getVal() << " +- " << std::setw(7) << std::right << AFB_meas[1][fl_q2bin]->getPropagatedError( *fit_result_sim );
  std::cout << "  [cf_slope(mm) = " << std::setw(7) << std::right << cf_slope[1][fl_q2bin]->getVal() << " ]"    << std::endl;
  
  std::cout << "Nsig(ee) = "  << std::setw(8) << std::right << nsig_q2_sim[0][fl_q2bin]->getVal()   << " = " << std::setw(7) << std::right << nsig_sim[0][2*fl_q2bin+0]->getVal() << " + " << std::setw(7) << std::right << nsig_sim[0][2*fl_q2bin+1]->getVal() << std::endl;
  std::cout << "Nsig(mm) = "  << std::setw(8) << std::right << nsig_q2_sim[1][fl_q2bin]->getVal()   << " = " << std::setw(7) << std::right << nsig_sim[1][2*fl_q2bin+0]->getVal() << " + " << std::setw(7) << std::right << nsig_sim[1][2*fl_q2bin+1]->getVal() << std::endl;
  std::cout << "Nsig(ee)E = " << std::setw(8) << std::right << nsig_q2_sim[0][fl_q2bin]->getError()
	    << " = "          << std::setw(8) << std::right << nsig_sim[0][2*fl_q2bin+0]->getPropagatedError( *fit_result_sim ) << " ++ " << std::setw(7) << std::right << nsig_sim[0][2*fl_q2bin+1]->getPropagatedError( *fit_result_sim ) << std::endl;
  std::cout << "Nsig(mm)E = " << std::setw(8) << std::right << nsig_q2_sim[1][fl_q2bin]->getError()
	    << " = "          << std::setw(8) << std::right << nsig_sim[1][2*fl_q2bin+0]->getPropagatedError( *fit_result_sim ) << " ++ " << std::setw(7) << std::right << nsig_sim[1][2*fl_q2bin+1]->getPropagatedError( *fit_result_sim ) << std::endl;

  // [Nsig for each bin]
  for( Int_t i=0; i<2; i++ ){
    Int_t j = fl_q2bin;
    std::cout << std::setw( 3) << std::right << j                << " ";
    std::cout << std::setw(10) << std::right << nsig_q2_sim[i][j]->getVal()                              << " "
	      << std::setw(10) << std::right << nsig_q2_sim[i][j]->getPropagatedError( *fit_result_sim ) << " ";
    std::cout << "  hooge" << i 
	      << "  func"  << sel_fun
	      << std::endl;
  }
  
  // [AFB]
  {
    Int_t j = fl_q2bin;
    std::cout << std::setw( 3) << std::right << j                                                  << " "
	      << std::setw(10) << std::right << AFB_true[j]->getVal()                              << " "
	      << std::setw(10) << std::right << AFB_true[j]->getPropagatedError( *fit_result_sim ) << " "
	      << "  hoooge" << 2 
	      << "  func"   << sel_fun
	      << std::endl;
  }
  
  //+++++++++++++++++++++++++++++++++ SAVE +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  for( Int_t i=0; i<Nroohist; i++ ){
    for( Int_t j=0; j<2; j++ ) c1[i][j]->Update();  // [hist, fitted hist]
  }
  c2->Update();
  c3->Update();
  c4->Update();

  if( flag_save ){
    TFile outfile( Form("pic/%s_fit_bin_dt_func%d_%dq2.root",  axis, sel_fun, fl_q2bin), "RECREATE" );
    for( Int_t i=0; i<Nroohist; i++ ){
      for( Int_t j=0; j<2; j++ ){ // [hist, fitted hist]
	Int_t k = j + 2*i;
	if( fl_q2bin==0 ) c1[i][j]->Print( Form("pic/%s_fit_bin_dt_mod_func%d_c1_%d.eps",  axis, sel_fun, k+1) );
	if( fl_q2bin==0 ) c1[i][j]->Write();
      }
    }

    c2->Print( Form("pic/%s_fit_bin_dt_mod_func%d_%dq2_c2.eps",  axis, sel_fun, fl_q2bin) );
    c2->Write();

    c3->Print( Form("pic/%s_fit_bin_dt_mod_func%d_%dq2_c3.eps",  axis, sel_fun, fl_q2bin) );
    c3->Write();
    c4->Print( Form("pic/%s_fit_bin_dt_mod_func%d_%dq2_c4.eps",  axis, sel_fun, fl_q2bin) );
    c4->Write();

    outfile.Close();
  }
  
  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();
  return 0;
}
