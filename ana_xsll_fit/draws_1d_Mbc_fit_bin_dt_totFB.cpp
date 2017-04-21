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

  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==2 || argc==3) ) std::cerr << "wrong input" << std::endl
					<< " Usage : ./draws_1d_Mbc_fit (int)sel_fun [(int)fl_appRun]" << std::endl
					<< " [sel_fun ] : 50(argus), 51(modified-argus), 15(gauss+argus), 151(gauss+modified-argus)" << std::endl
					<< std::endl, abort();
  Int_t    sel_fun      = atoi(argv[1]);
  Int_t    fl_appRun    = 1;
  if( argc==3 ) fl_appRun = atoi( argv[2] );

  const Int_t Nbin_afb = 2;
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

      
      if( i ==  0 ){
	sTmp << " && ( "
	     << makeCut_q2( tmp_fl_mode_ll, 1,  1 ).c_str()
	     << " || "
	     << makeCut_q2( tmp_fl_mode_ll, 2,  1 ).c_str()
	     << " || "
	     << makeCut_q2( tmp_fl_mode_ll, 3,  1 ).c_str()
	     << " || "
	     << makeCut_q2( tmp_fl_mode_ll, 5,  1 ).c_str()
	     << " || "
	     << makeCut_q2( tmp_fl_mode_ll, 7,  1 ).c_str()
	     << " || "
	     << makeCut_q2( tmp_fl_mode_ll, 8,  1 ).c_str()
	     << " || "
	     << makeCut_q2( tmp_fl_mode_ll, 9,  1 ).c_str()
	     << " )";
      }else if( i ==  1 ){
	sTmp << " && ( "
	     << makeCut_q2( tmp_fl_mode_ll, 1, -1 ).c_str()
	     << " || "
	     << makeCut_q2( tmp_fl_mode_ll, 2, -1 ).c_str()
	     << " || "
	     << makeCut_q2( tmp_fl_mode_ll, 3, -1 ).c_str()
	     << " || "
	     << makeCut_q2( tmp_fl_mode_ll, 5, -1 ).c_str()
	     << " || "
	     << makeCut_q2( tmp_fl_mode_ll, 7, -1 ).c_str()
	     << " || "
	     << makeCut_q2( tmp_fl_mode_ll, 8, -1 ).c_str()
	     << " || "
	     << makeCut_q2( tmp_fl_mode_ll, 9, -1 ).c_str()
	     << " )";
      }
      
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
	  if     ( i==0 && j== 0 ){ area = 0.20; norm = 60.0; shape = -20.0; } // ee
	  else if( i==1 && j== 0 ){ area = 0.20; norm = 30.0; shape = -20.0; }
	  else if( i==0 && j== 1 ){ area = 0.30; norm = 85.0; shape = -13.0; }// mm
	  else if( i==1 && j== 1 ){ area = 0.15; norm = 40.0; shape = -20.0; }
	  else if( i==0 && j== 2 ){ area = 0.50; norm =  145; shape = -17.0; } // ee+mm
	  else if( i==1 && j== 2 ){ area = 0.35; norm =   70; shape = -19.0; }

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
    tag[i]->defineType ( "1235789_p",  1 ); tag[i]->defineType ( "1235789_m",  2 );

  }

 
  for( Int_t i=0; i<Nroohist; i++ ) obs[i] = new RooRealVar( Form("%s%d",axis,i), xlabel, offset+xmin_fit, offset+xmax_fit );

  RooRealVar*  obs_sim = new RooRealVar( Form("%s_sim",axis), xlabel, offset+xmin_fit, offset+xmax_fit );
  RooCategory* tag_sim = new RooCategory( "bin_sim", "bin_sim" );
  tag_sim->defineType ( "1235789_p_e",  101 ); tag_sim->defineType ( "1235789_m_e",  102 );
  tag_sim->defineType ( "1235789_p_mu",   1 ); tag_sim->defineType ( "1235789_m_mu",   2 );

  
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
	      if( Nbin_afb == 2 && !((tag_q2_cos_region_totFB((Int_t)rm_l, cc_obs*cc_obs, cos_obs)== 1) || (tag_q2_cos_region_totFB((Int_t)rm_l, cc_obs*cc_obs, cos_obs)== 2))  ) continue;
	      if( fl_sb[j]==2 || fl_sb[j]==5 ){ // RD
		tag [n]->setIndex( tag_q2_cos_region_totFB ( (Int_t)rm_l, cc_obs*cc_obs, cos_obs) );
		data[n]->add( RooArgSet(*obs[n],*tag[n]) );
		if( n==Nroohist-1 ){
		  tag_sim ->setIndex( 100*(Int_t)rm_l + tag_q2_cos_region_totFB( (Int_t)rm_l, cc_obs*cc_obs, cos_obs) );
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
  TFile file_scf("pdf_totFB/Mbc_self_cf2_totFB_setA-U.root");
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

  r_tt[ 0][ 0] = new RooRealVar( "r_tt0",    "r_tt0",      0.1363);   r_tf[ 0][ 0] = new RooRealVar( "r_tf0",    "r_tf0",      0.0435);   r_f [ 0][ 0] = new RooRealVar( "r_f0",     "r_f0",     0.002444); 
  r_tt[ 0][ 1] = new RooRealVar( "r_tt1",    "r_tt1",      0.1199);   r_tf[ 0][ 1] = new RooRealVar( "r_tf1",    "r_tf1",     0.04886);   r_f [ 0][ 1] = new RooRealVar( "r_f1",     "r_f1",     0.002682); 
  r_tt[ 1][ 0] = new RooRealVar( "r_tt2",    "r_tt0",      0.1091);   r_tf[ 1][ 0] = new RooRealVar( "r_tf2",    "r_tf0",     0.03689);   r_f [ 1][ 0] = new RooRealVar( "r_f2",     "r_f0",     0.002515); 
  r_tt[ 1][ 1] = new RooRealVar( "r_tt3",    "r_tt1",     0.09791);   r_tf[ 1][ 1] = new RooRealVar( "r_tf3",    "r_tf1",     0.04008);   r_f [ 1][ 1] = new RooRealVar( "r_f3",     "r_f1",     0.002543); 
  r_tt[ 2][ 0] = new RooRealVar( "r_tt4",    "r_tt0",      0.1213);   r_tf[ 2][ 0] = new RooRealVar( "r_tf4",    "r_tf0",     0.04002);   r_f [ 2][ 0] = new RooRealVar( "r_f4",     "r_f0",     0.004959); 
  r_tt[ 2][ 1] = new RooRealVar( "r_tt5",    "r_tt1",      0.1083);   r_tf[ 2][ 1] = new RooRealVar( "r_tf5",    "r_tf1",     0.04402);   r_f [ 2][ 1] = new RooRealVar( "r_f5",     "r_f1",     0.005226);
  
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // SETTING of INITIAL VALUE
  for( Int_t i=0; i<Nroohist; i++ ){
    nsig_q2    [i] = new RooRealVar*[Nbin_afb/2];
    nsig_q2_sim[i] = new RooRealVar*[Nbin_afb/2];
    AFB        [i] = new RooRealVar*[Nbin_afb/2];
    for( Int_t j=0; j<Nbin_afb/2; j++ ){
      Int_t k = j + i*Nbin_afb/2;
      Double_t tmp_nsig_q2 = 220;
      Double_t tmp_AFB     =  0;
      nsig_q2    [i][j] = new RooRealVar( Form("N_{sig}^{q^{2}}%d_%d",    i,j), Form("N_{sig}^{q^{2}}%d_%d",    i,j), tmp_nsig_q2, tmp_nsig_q2-200, tmp_nsig_q2+200 );
      nsig_q2_sim[i][j] = new RooRealVar( Form("N_{sig}^{q^{2},sim}%d_%d",i,j), Form("N_{sig}^{q^{2},sim}%d_%d",i,j), tmp_nsig_q2, tmp_nsig_q2-200, tmp_nsig_q2+200 );
      AFB        [i][j] = new RooRealVar( Form("A_{FB}%d_%d",             i,j), Form("A_{FB}%d_%d",             i,j), tmp_AFB,     -1.5,           1.5            );
      AFB_true      [j] = new RooRealVar( Form("A_{FB}^{true}_%d",          j), Form("A_{FB}^{true}_%d",          j), tmp_AFB,     -1.5,           1.5            );
      std::cout << Form("[i=%d, j=%d] nsig_q2 = %f, tmp_AFB = %f", i,j,tmp_nsig_q2,tmp_AFB) << std::endl;
    }
  }
  
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  RooRealVar*** cf_slope = new RooRealVar**[Nroohist];
  for( Int_t i=0; i<Nroohist;   i++ ) cf_slope[i]    = new RooRealVar*[Nbin_afb/2];
  for( Int_t j=0; j<Nbin_afb/2; j++ ) cf_slope[2][j] = new RooRealVar( Form("cf_2_%d",j), Form("cf_2_%d",j), 0.0 ); // dummy
  
  cf_slope[ 0][ 0] = new RooRealVar( "cf_0_0", "cf_0_0",    1.000 ); // dummy
  cf_slope[ 0][ 1] = new RooRealVar( "cf_0_1", "cf_0_1",    1.000 ); // dummy
  cf_slope[ 0][ 2] = new RooRealVar( "cf_0_2", "cf_0_2",    1.000 ); // dummy
  cf_slope[ 0][ 3] = new RooRealVar( "cf_0_3", "cf_0_3",    1.000 ); // dummy
  cf_slope[ 1][ 0] = new RooRealVar( "cf_1_0", "cf_1_0",    1.000 ); // dummy
  cf_slope[ 1][ 1] = new RooRealVar( "cf_1_1", "cf_1_1",    1.000 ); // dummy
  cf_slope[ 1][ 2] = new RooRealVar( "cf_1_2", "cf_1_2",    1.000 ); // dummy
  cf_slope[ 1][ 3] = new RooRealVar( "cf_1_3", "cf_1_3",    1.000 ); // dummy

  for( Int_t i=0; i<Nroohist; i++ ){
    AFB_meas[i] = new RooFormulaVar*[Nbin_afb];
    for( Int_t j=0; j<Nbin_afb/2; j++ ){
      AFB_meas[i][j] = new RooFormulaVar( Form("A_{FB}^{meas.}_%d_%d",i,j),  "1/@0*@1", RooArgSet(*cf_slope[i][j], *AFB_true[j]) );
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
      Double_t tmp_nbkg = 600;
      if( i==Nroohist-1 ) tmp_nbkg *= 2;
      nbkg    [i][j] = new RooRealVar ( Form("N_{bkg}%d_%d",      i,j), Form("N_{bkg}%d_%d",      i,j), tmp_nbkg, tmp_nbkg-1000, tmp_nbkg+1500 );
      nbkg_sim[i][j] = new RooRealVar ( Form("N_{bkg}^{sim}%d_%d",i,j), Form("N_{bkg}^{sim}%d_%d",i,j), tmp_nbkg, tmp_nbkg-1000, tmp_nbkg+1500 );
    }
  }
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // READ PDF (charmonium)
  TFile file_cc("pdf_totFB/Mbc_peak_cc2_totFB_s00.root"); // CC-MC
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
  const Double_t scale_cc[3] = {0.823, 0.741, 0.779 };
  for( Int_t i=0; i<Nroohist; i++ ) r_peak_cc[i] = new RooRealVar*[Nbin_afb+1];
  // CC-MC
  r_peak_cc[ 0][ 0] = new RooRealVar( "r_peak_cc0", "r_peak_cc0", scale_cc[0]*4.22 );
  r_peak_cc[ 0][ 1] = new RooRealVar( "r_peak_cc1", "r_peak_cc1", scale_cc[0]*3.96 );
  r_peak_cc[ 0][ 2] = new RooRealVar( "r_peak_cc2", "r_peak_cc2", scale_cc[0]*8.18 );
  r_peak_cc[ 1][ 0] = new RooRealVar( "r_peak_cc3", "r_peak_cc3", scale_cc[1]*5.37 );
  r_peak_cc[ 1][ 1] = new RooRealVar( "r_peak_cc4", "r_peak_cc4", scale_cc[1]*5.43 );
  r_peak_cc[ 1][ 2] = new RooRealVar( "r_peak_cc5", "r_peak_cc5", scale_cc[1]*10.82 );
  r_peak_cc[ 2][ 0] = new RooRealVar( "r_peak_cc6", "r_peak_cc6", scale_cc[2]*9.59 );
  r_peak_cc[ 2][ 1] = new RooRealVar( "r_peak_cc7", "r_peak_cc7", scale_cc[2]*9.39 );
  r_peak_cc[ 2][ 2] = new RooRealVar( "r_peak_cc8", "r_peak_cc8", scale_cc[2]*19 );
  


  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // READ PDF (double miss-PID)
  TFile file_double("pdf_data_totFB/Mbc_peak_double2_totFB_s00.root"); // RD
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
  for( Int_t i=0; i<Nroohist; i++ ) r_peak_double[i] = new RooRealVar*[Nbin_afb+1];
  r_peak_double[ 0][ 0] = new RooRealVar( "r_peak_double0", "r_peak_double0", 0.1205/0.273 );
  r_peak_double[ 0][ 1] = new RooRealVar( "r_peak_double1", "r_peak_double1", 0.09409/0.273 );
  r_peak_double[ 0][ 2] = new RooRealVar( "r_peak_double2", "r_peak_double2", 0.2146/0.273 );
  r_peak_double[ 1][ 0] = new RooRealVar( "r_peak_double3", "r_peak_double3", 7.189/0.273 );
  r_peak_double[ 1][ 1] = new RooRealVar( "r_peak_double4", "r_peak_double4", 5.147/0.273 );
  r_peak_double[ 1][ 2] = new RooRealVar( "r_peak_double5", "r_peak_double5", 12.34/0.273 );
  r_peak_double[ 2][ 0] = new RooRealVar( "r_peak_double6", "r_peak_double6", 7.309/0.273 );
  r_peak_double[ 2][ 1] = new RooRealVar( "r_peak_double7", "r_peak_double7", 5.241/0.273 );
  r_peak_double[ 2][ 2] = new RooRealVar( "r_peak_double8", "r_peak_double8", 12.55/0.273 );
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // READ PDF (swapped miss-PID)
  TFile file_swap("pdf_data_totFB/Mbc_peak_swap2_totFB_s00.root"); // RD
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
  for( Int_t i=0; i<Nroohist; i++ ) r_peak_swap[i] = new RooRealVar*[Nbin_afb+1];
  r_peak_swap[ 0][ 0] = new RooRealVar( "r_peak_swap0", "r_peak_swap0", 0.06054 );
  r_peak_swap[ 0][ 1] = new RooRealVar( "r_peak_swap1", "r_peak_swap1", 0.03178 );
  r_peak_swap[ 0][ 2] = new RooRealVar( "r_peak_swap2", "r_peak_swap2", 0.09232 );
  r_peak_swap[ 1][ 0] = new RooRealVar( "r_peak_swap3", "r_peak_swap3", 2.351 );
  r_peak_swap[ 1][ 1] = new RooRealVar( "r_peak_swap4", "r_peak_swap4", 1.579 );
  r_peak_swap[ 1][ 2] = new RooRealVar( "r_peak_swap5", "r_peak_swap5", 3.93 );
  r_peak_swap[ 2][ 0] = new RooRealVar( "r_peak_swap6", "r_peak_swap6", 2.412 );
  r_peak_swap[ 2][ 1] = new RooRealVar( "r_peak_swap7", "r_peak_swap7", 1.61 );
  r_peak_swap[ 2][ 2] = new RooRealVar( "r_peak_swap8", "r_peak_swap8", 4.022 );
    
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
      pdf    [i][j] = new RooAddPdf( Form("pdf%d",    k), Form("pdf%d",    k), RooArgList(*gauss    [i], *modargus    [i][j], *pdf_scf_tt    [i][j], *pdf_scf_tf    [i][j], *pdf_peak_cc    [i][j], *pdf_peak_double    [i][j], *pdf_peak_swap    [i][j]), RooArgList(*nsig     [i][j], *nbkg    [i][j], *nsig_tt    [i][j], *nsig_tf    [i][j], *r_peak_cc[i][j], *r_peak_double[i][j], *r_peak_swap[i][j] ) ); // signal + scf + peak
      pdf_sim[i][j] = new RooAddPdf( Form("pdf_sim%d",k), Form("pdf_sim%d",k), RooArgList(*gauss_sim[i], *modargus_sim[i][j], *pdf_scf_tt_sim[i][j], *pdf_scf_tf_sim[i][j], *pdf_peak_cc_sim[i][j], *pdf_peak_double_sim[i][j], *pdf_peak_swap_sim[i][j]), RooArgList(*nsig_sim [i][j], *nbkg_sim[i][j], *nsig_tt_sim[i][j], *nsig_tf_sim[i][j], *r_peak_cc[i][j], *r_peak_double[i][j], *r_peak_swap[i][j] ) ); // signal + scf + peak
    }
  }
  
  
  // FIT (ee, mm, ee+mm)
  for( Int_t i=0; i<Nroohist; i++ ){
    tot_pdf[i]->addPdf( *pdf[i][ 0],  "1235789_p"  );
    tot_pdf[i]->addPdf( *pdf[i][ 1],  "1235789_m"  );
    fit_result[i] = tot_pdf[i]->fitTo( *data[i], Extended(), Save(true) );
  }

  // SIMULTANEOUS FIT( ee and mm )
    tot_pdf_sim->addPdf( *pdf_sim[0][ 0],  "1235789_p_e"   );
    tot_pdf_sim->addPdf( *pdf_sim[0][ 1],  "1235789_m_e"   );
    tot_pdf_sim->addPdf( *pdf_sim[1][ 0],  "1235789_p_mu"  );
    tot_pdf_sim->addPdf( *pdf_sim[1][ 1],  "1235789_m_mu"  );

  RooFitResult* fit_result_sim = tot_pdf_sim->fitTo( *data_sim, Extended(), Save(true) );
  
  // +++++++++++++++++++++++++ Draw(RooPlot) +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  RooPlot**    frame_param = new RooPlot*   [Nroohist];
  RooPlot***   frame       = new RooPlot**  [Nroohist]; // [Nbin_afb];
  TPaveText*** box         = new TPaveText**[Nroohist]; // [Nbin_afb];

  for( Int_t i=0; i<Nroohist; i++ ){
    frame[i] = new RooPlot*  [Nbin_afb];
    box  [i] = new TPaveText*[Nbin_afb];
    for( Int_t j=0; j<2; j++ ){
      Int_t k = j+i*Nbin_afb;
      frame[i][j] = obs[i]->frame();
      frame[i][j]->GetXaxis()->CenterTitle();
      frame[i][j]->GetYaxis()->CenterTitle();
      frame[i][j]->SetTitleOffset( 1.00,"x" );
      frame[i][j]->SetTitleOffset( 1.30,"y" );
      if     ( i==0 ) frame[i][j]->SetTitle( Form("bin=%d, ee       ",j+1) );
      else if( i==1 ) frame[i][j]->SetTitle( Form("bin=%d, #mu#mu   ",j+1) );
      else if( i==2 ) frame[i][j]->SetTitle( Form("bin=%d, ee+#mu#mu",j+1) );

      if     ( j== 0 ) data[i]->plotOn ( frame[i][j], Binning(xbin), LineWidth(1), Cut( Form("bin%d==bin%d::1235789_p" ,i,i) ) );
      else if( j== 1 ) data[i]->plotOn ( frame[i][j], Binning(xbin), LineWidth(1), Cut( Form("bin%d==bin%d::1235789_m" ,i,i) ) );
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
    for( Int_t j=0; j<2; j++ ){
      Int_t k = j+i*Nbin_afb;
      frame_sim[i][j] = obs_sim->frame();
      frame_sim[i][j]->GetXaxis()->CenterTitle();
      frame_sim[i][j]->GetYaxis()->CenterTitle();
      frame_sim[i][j]->SetTitleOffset( 1.00,"x" );
      frame_sim[i][j]->SetTitleOffset( 1.30,"y" );
      if     ( i==0 ) frame_sim[i][j]->SetTitle( Form("bin=%d, ee",    j+1) );
      else if( i==1 ) frame_sim[i][j]->SetTitle( Form("bin=%d, #mu#mu",j+1) );

      if( i==0 ){
	if     ( j== 0 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::1235789_p_e"  ) );
	else if( j== 1 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::1235789_m_e"  ) );
      }else{
	if     ( j== 0 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::1235789_p_mu"  ) );
	else if( j== 1 ) data_sim->plotOn ( frame_sim[i][j], Binning(xbin), LineWidth(1), Cut( "bin_sim==bin_sim::1235789_m_mu"  ) );
      }

      tag_sim->setIndex(i==0 ? 100+j+1 : j+1 ); tot_pdf_sim->plotOn( frame_sim[i][j], Slice(*tag_sim), ProjWData(*tag_sim,*data_sim), LineWidth(1) );
      box_sim[i][j] = new TPaveText( 0.05, 0.87, 0.35, 0.93,"BRNDC" );
      box_sim[i][j]->AddText( Form("#chi^{2}/ndf = %f", frame_sim[i][j]->chiSquare(nparam+N0bin[j][rooad[i]])) );
      tag_sim->setIndex(i==0 ? 100+j+1 : j+1); tot_pdf_sim->plotOn ( frame_sim[i][j], Slice(*tag_sim), ProjWData(*tag_sim,*data_sim), Components( RooArgSet(*gauss_sim[i],*pdf_scf_tt_sim [i][j],*pdf_scf_tf_sim[i][j])                 ), LineStyle(7), LineColor(5), LineWidth(1) ); // sig + self cros-feed
      //tag_sim->setIndex(i==0 ? 100+j+1 : j+1); tot_pdf_sim->plotOn ( frame_sim[i][j], Slice(*tag_sim), ProjWData(*tag_sim,*data_sim), Components( *gauss_sim   [i]                                                                      ), LineStyle(7), LineColor(2), LineWidth(1) ); // sig (gauss)
      tag_sim->setIndex(i==0 ? 100+j+1 : j+1); tot_pdf_sim->plotOn ( frame_sim[i][j], Slice(*tag_sim), ProjWData(*tag_sim,*data_sim), Components( *modargus_sim[i][j]                                                                   ), LineStyle(7), LineColor(4), LineWidth(1) ); // bkg (arguss)
      //tag_sim->setIndex(i==0 ? 100+j+1 : j+1); tot_pdf_sim->plotOn ( frame_sim[i][j], Slice(*tag_sim), ProjWData(*tag_sim,*data_sim), Components( RooArgSet(*pdf_scf_tt_sim [i][j],*pdf_scf_tf_sim     [i][j])                          ), LineStyle(7), LineColor(5), LineWidth(1) ); // self cros-feed
      //tag_sim->setIndex(i==0 ? 100+j+1 : j+1); tot_pdf_sim->plotOn ( frame_sim[i][j], Slice(*tag_sim), ProjWData(*tag_sim,*data_sim), Components( RooArgSet(*pdf_peak_cc_sim[i][j],*pdf_peak_double_sim[i][j],*pdf_peak_swap_sim[i][j]) ), LineStyle(7), LineColor(6), LineWidth(1) ); // peak
      tag_sim->setIndex(i==0 ? 100+j+1 : j+1); tot_pdf_sim->plotOn ( frame_sim[i][j], Slice(*tag_sim), ProjWData(*tag_sim,*data_sim), Components( RooArgSet(*pdf_peak_cc_sim[i][j],*pdf_peak_double_sim[i][j],*pdf_peak_swap_sim[i][j]) ), FillStyle(1001), FillColor(11), LineWidth(1), VLines(), DrawOption("F") ); // peak // Fill Style
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
    frame[i][0]->SetAxisRange( 5.22, 5.30, "X" );
    frame[i][1]->SetAxisRange( 5.22, 5.30, "X" );
    if( i==0 ){
      frame[i][0]->SetAxisRange( 0, 70, "Y" ); // ee
      frame[i][1]->SetAxisRange( 0, 70, "Y" ); // ee
    }else if( i==1 ){
      frame[i][0]->SetAxisRange( 0, 95, "Y" ); // mm
      frame[i][1]->SetAxisRange( 0, 95, "Y" ); // mm
    }else if( i==2 ){
      frame[i][0]->SetAxisRange( 0, 160, "Y" ); // ll
      frame[i][1]->SetAxisRange( 0, 160, "Y" ); // ll
    }
  }
  //*/
  for( Int_t i=0; i<Nroohist; i++ ){
    c2->cd(3*i+1);
    frame[i][0]->Draw();
    c2->cd(3*i+2);
    frame[i][1]->Draw();
    c2->cd(3*i+3);
    frame_param[i]->Draw();
  }

  // +++++++++++++++++++++++++ Draw(TCanvas, simultaneous fit) +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TCanvas* c3 = Canvas( "c3", "c3", 3, 2 );
  ///*
  for( Int_t i=0; i<2; i++ ){
    frame_sim[i][0]->SetAxisRange( 5.22, 5.30, "X" );
    frame_sim[i][1]->SetAxisRange( 5.22, 5.30, "X" );
    if( i==0 ){
      frame_sim[i][0]->SetAxisRange( 0, 70, "Y" ); // ee
      frame_sim[i][1]->SetAxisRange( 0, 70, "Y" ); // ee
    }else if( i==1 ){
      frame_sim[i][0]->SetAxisRange( 0, 95, "Y" ); // mm
      frame_sim[i][1]->SetAxisRange( 0, 95, "Y" ); // mm
    }
  }
  //*/
  for( Int_t i=0; i<2; i++ ){
    c3->cd(3*i+1);
    frame_sim[i][0]->Draw();
    c3->cd(3*i+2);
    frame_sim[i][1]->Draw();
    c3->cd(3*i+3);
    frame_param_sim[i]->Draw();
  }
  
  //+++++++++++++++++++++++++++++++++ LOG +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::cout << std::endl;
  //std::cout << "correlation matrix" << std::endl; fit_result_sim->correlationMatrix().Print();
  //std::cout << "covariance matrix"  << std::endl; fit_result_sim->covarianceMatrix ().Print();

  std::cout << "AFB(true)    = "    << std::setw(8) << std::right << AFB_true[0]->getVal()    << " +- " << std::setw(7) << std::right << AFB_true[0]->getError()                               << std::endl;
  std::cout << "AFB(meas,ee) = "    << std::setw(8) << std::right << AFB_meas[0][0]->getVal() << " +- " << std::setw(7) << std::right << AFB_meas[0][0]->getPropagatedError( *fit_result_sim );
  std::cout << "  [cf_slope(ee) = " << std::setw(7) << std::right << cf_slope[0][0]->getVal() << " ]"    << std::endl;
  std::cout << "AFB(meas,mm) = "    << std::setw(8) << std::right << AFB_meas[1][0]->getVal() << " +- " << std::setw(7) << std::right << AFB_meas[1][0]->getPropagatedError( *fit_result_sim );
  std::cout << "  [cf_slope(mm) = " << std::setw(7) << std::right << cf_slope[1][0]->getVal() << " ]"    << std::endl;
  
  std::cout << "Nsig(ee) = "  << std::setw(8) << std::right << nsig_q2_sim[0][0]->getVal()   << " = " << std::setw(7) << std::right << nsig_sim[0][0]->getVal() << " + " << std::setw(7) << std::right << nsig_sim[0][1]->getVal() << std::endl;
  std::cout << "Nsig(mm) = "  << std::setw(8) << std::right << nsig_q2_sim[1][0]->getVal()   << " = " << std::setw(7) << std::right << nsig_sim[1][0]->getVal() << " + " << std::setw(7) << std::right << nsig_sim[1][1]->getVal() << std::endl;
  std::cout << "Nsig(ee)E = " << std::setw(8) << std::right << nsig_q2_sim[0][0]->getError()
	    << " = "          << std::setw(8) << std::right << nsig_sim[0][0]->getPropagatedError( *fit_result_sim ) << " ++ " << std::setw(7) << std::right << nsig_sim[0][1]->getPropagatedError( *fit_result_sim ) << std::endl;
  std::cout << "Nsig(mm)E = " << std::setw(8) << std::right << nsig_q2_sim[1][0]->getError()
	    << " = "          << std::setw(8) << std::right << nsig_sim[1][0]->getPropagatedError( *fit_result_sim ) << " ++ " << std::setw(7) << std::right << nsig_sim[1][1]->getPropagatedError( *fit_result_sim ) << std::endl;

  // [Nsig for each bin]
  for( Int_t i=0; i<2; i++ ){
    Int_t j = 0;
    std::cout << std::setw( 3) << std::right << j                << " ";
    std::cout << std::setw(10) << std::right << nsig_q2_sim[i][j]->getVal()                              << " "
	      << std::setw(10) << std::right << nsig_q2_sim[i][j]->getPropagatedError( *fit_result_sim ) << " ";
    std::cout << "  hooge" << i 
	      << "  func"  << sel_fun
	      << std::endl;
  }
  
  // [AFB]
  {
    Int_t j = 0;
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

  if( flag_save ){
    TFile outfile( Form("pic/%s_fit_bin_dt_func%d_totq2.root",  axis, sel_fun), "RECREATE" );
    for( Int_t i=0; i<Nroohist; i++ ){
      for( Int_t j=0; j<2; j++ ){ // [hist, fitted hist]
	Int_t k = j + 2*i;
	c1[i][j]->Print( Form("pic/%s_fit_bin_dt_func%d_totq2_c1_%d.eps",  axis, sel_fun, k+1) );
	c1[i][j]->Write();
      }
    }

    c2->Print( Form("pic/%s_fit_bin_dt_func%d_totq2_c2.eps",  axis, sel_fun) );
    c2->Write();

    c3->Print( Form("pic/%s_fit_bin_dt_func%d_totq2_c3.eps",  axis, sel_fun) );
    c3->Write();

    outfile.Close();
  }
  
  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();
  return 0;
}
