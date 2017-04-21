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
#include <RooSimultaneous.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooAbsReal.h>

void manip_func( TF1* func ){
  func->SetParLimits( 1, 5.276, 5.282 );
  func->SetParLimits( 2, 0.001, 0.005 );
  func->FixParameter( 4, 5.289 );
}

Int_t main( Int_t argc, Char_t** argv ){
  //using namespace sig_gmc_rd_calib;
  //using namespace sig_gmc_rd_cut3_beforebgsup;
  using namespace sig_gmc_rd_cut3;
  
  using namespace Mbc_bkg;
  using namespace RooFit;
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==3 || argc==4) ) std::cerr << "wrong input" << std::endl
    					<< " Usage : ./draws_1d_Mbc (char*)stream (double)used_nstream [(int)fl_appRun]" << std::endl
					<< std::endl, abort();
  Char_t*  stream        = argv[1];
  Double_t used_nstream = atof(argv[2]); 
  Int_t   fl_appRun     = 1;
  if( argc==4 ) fl_appRun = atoi( argv[3] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Bool_t flag_k4pi      = true; // 1(veto  K4pi    modes)
  const Bool_t flag_unflavor  = true; // 1(veto unflavor modes)
  const Bool_t flag_ccpi0     = true; // 1(veto ccpi0 peak    )
  const Bool_t flag_save      = true; // outfile.eps and outfile.root
  const Bool_t flag_fit       = true;
  const Int_t  sel_fun        = 15;   // 15(gauss+argus)
  const Int_t  Nfun           = 6;    // [MC,RD] x [ee,mm,ee+mm]
  const Int_t  sel_hist[Nfun] = {1,4,7,2,5,8}; // MC(ee), MC(mm), MC(ee+mm), RD(ee), RD(mm), RD(ee+mm)
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t    Nchain             = 3*2; // [gmc(qq,bb),rd] x [ee,mm]
  const Int_t    Nhist              = 3*3; // [gmc(qq,bb),rd] x [ee,mm,ee+mm]
  const Int_t    nfile[Nchain]      = {0};
  const Int_t    fl_sb[Nchain]      = {0,0,2,0,0,2}; // 0(bkg), 1(sig), 2(rd)
  const Int_t    fl_mode_ll[Nchain] = {101,101,101,  // 1(ee)
				       100,100,100}; // 0(mm)
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
    //if     ( fl_sb[i]==0 ) sTmp << indir[fl_sb[i]] << "gMC_*_e031*s0[" << stream << "]_";   // bkg // tmppppppppppppppp
    //else if( fl_sb[i]==2 ) sTmp << indir[fl_sb[i]] << "RD_e031";                            // rd // tmpppppppppppppp
    if     ( fl_sb[i]==0 ) sTmp << indir[fl_sb[i]] << "gMC_*_s0[" << stream << "]";     // bkg
    else if( fl_sb[i]==2 ) sTmp << indir[fl_sb[i]] << "RD_";                            // rd
    strcpy( infile[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t add[Nhist][Nchain] ={
    {1},     // gmc(qq), ee
    {1,1},   // gmc(bb), ee
    {0,0,1}, // rd,      ee
    {0,0,0,1},     // gmc(qq), mm
    {0,0,0,1,1},   // gmc(bb), mm
    {0,0,0,0,0,1}, // rd,      mm
    {1,0,0,1},        // gmc(qq), ee+mm
    {1,1,0,1,1},      // gmc(bb), ee+mm
    {0,0,1,0,0,1},    // rd,      ee+mm
  };

  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    add_cut[i] = new Char_t[4096];
    sTmp << LRNB_cut;
    if( flag_k4pi     ) sTmp << " && " << makeCut_5body_veto().c_str();
    if( flag_unflavor ) sTmp << " && " << makeCut_unflavor_veto().c_str();
    if( flag_ccpi0    ) sTmp << " && " << cut_ccpi0;

    if     ( i==0 || i==3 ) sTmp << " && genbfl==0"; // gmc(qq)
    else if( i==1 || i==4 ) sTmp << " && genbfl!=0"; // gmc(bb)
    
    /*
    sTmp << "&& ( "//  with   pi0 modes
	 << Form( "(rm_xs==1001 && abs(xs_m-%f)<0.050) || ", PDGmass::kstrp ) // K+pi0
	 << Form( "(rm_xs==1010 && abs(xs_m-%f)<0.050)",     PDGmass::kstr0 ) // Kspi0
	 << " ) ";
    */
    /*
      sTmp << "&& ( " // Kll
      << Form( " (rm_xs==  1) || " ) // K+
      << Form( " (rm_xs== 10) || " ) // Ks
      << Form( " (rm_xs==101 && abs(xs_m-%f)<0.050) || ", PDGmass::kstr0 ) // K+pi-
      << Form( " (rm_xs==110 && abs(xs_m-%f)<0.050) ",    PDGmass::kstrp ) // Kspi-
      << " )";
    */
    /*
    sTmp << "&& ( "//  with K/K* ll
	 << Form( " (rm_xs==  1) || " ) // K+
	 << Form( " (rm_xs== 10) || " ) // Ks
	 << Form( " (rm_xs==101 && abs(xs_m-%f)<0.050) || ", PDGmass::kstr0 ) // K+pi-
	 << Form( " (rm_xs==110 && abs(xs_m-%f)<0.050) || ", PDGmass::kstrp ) // Kspi-
	 << Form( "(rm_xs==1001 && abs(xs_m-%f)<0.050) || ", PDGmass::kstrp ) // K+pi0
	 << Form( "(rm_xs==1010 && abs(xs_m-%f)<0.050) ",    PDGmass::kstr0 ) // Kspi0
	 << " ) ";
    */
    /*
    sTmp << "&& ( "//  with kll
	 << Form( " (rm_xs==  1) || " ) // K+
	 << Form( " (rm_xs== 10) "    ) // Ks
	 << " ) ";
    */

    strcpy( add_cut[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain   = new MChain*[Nchain];
  TH1D**    tmphist = new TH1D*  [Nchain];
  TH1D**    hist    = new TH1D*  [Nhist];
  TF1**     func    = new TF1*   [Nfun];
  TCanvas*  c1      = Canvas( "c1","c1",3, 3 );
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
    chain[j]->GetCut()->Set( "Mbc", 1,  5.20, 0.0, 5.30 );
  }
  //
  // ++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ changed cut *************************************" << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<infile %d > ", j ) << std::endl;
    chain[j]->GetCut()->Display(0);
    if( fl_sb[j]==2 && fl_mode_ll[j]==101 ){ // for e
      chain[j]->GetCut()->Set( "decalib",  0                   );
      chain[j]->GetCut()->Set( "de",       1, -0.10, 0.0, 0.05 );
    }else if( fl_sb[j]==2 && fl_mode_ll[j]==100 ){ // for mu
      chain[j]->GetCut()->Set( "decalib",  0                   );
      chain[j]->GetCut()->Set( "de",       1, -0.05, 0.0, 0.05 );
    }
    chain[j]->MakeTree();
  }
  // +++++++ make tmphist ++++++++++++++++++++++++++++++++++
  for( Int_t j=0; j<Nchain; j++ ){
    tmphist[j] = new TH1D( Form("tmphist%d",j), Form("%s",chain[j]->GetChange()), xbin,offset+xmin,offset+xmax );
    chain[j]->GetTree()->Project( Form("tmphist%d",j), axis, add_cut[j] );
    tmphist[j]->Sumw2();
    if( fl_sb[j]==0 ) tmphist[j]->Scale( 1/scale_event_bkg );
  }
  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make hist *************************************" << std::endl;
  for( Int_t i=0; i<Nhist; i++ ){
    hist[i] = new TH1D( Form("hist%d",i),Form("hist%d",i), xbin,offset+xmin,offset+xmax );
    //if( i!=3 && i!=7 && i!=11 ) Deco( hist[i], 3, col_fil[i%4], col_fil[i%4] );
    Deco( hist[i], 3, col_fil[i%3], col_fil[i%3] );
    for( Int_t j=0; j<Nchain; j++ ){
      if( add[i][j] ) hist[i]->Add( tmphist[j] );
    }
  }
  hist[0]->SetTitle("gMC(qq), ee"   );
  hist[1]->SetTitle("gMC(bb), ee"   );
  hist[2]->SetTitle("RD,      ee"   );
  hist[3]->SetTitle("gMC(qq), mm"   );
  hist[4]->SetTitle("gMC(bb), mm"   );
  hist[5]->SetTitle("RD,      mm"   );
  hist[6]->SetTitle("gMC(qq), ee+mm");
  hist[7]->SetTitle("gMC(bb), ee+mm");
  hist[8]->SetTitle("RD,      ee+mm");

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
    for( Int_t i=0; i<Nfun; i++ ){
      func[i] = new TF1( Form("func%d",i), make_func(sel_fun),offset+xmin_fit,offset+xmax_fit,n_fitfunc_par(sel_fun) );
      if( sel_fun == 15 ){
	Double_t mu    = 5.279;
	Double_t sigma = 0.003;
	Double_t area  = ( hist[i]->GetMaximum() - hist[i]->GetBinContent(1) )*sqrt( TMath::TwoPi() )*sigma;
	func[i]->SetParNames  ( "area","mean","sigma","norm",                          "Ebeam","shape" );
	func[i]->SetParameters( area,   mu,    sigma,  100*hist[i]->GetBinContent(1),  5.289,   -30 );
	func[i]->SetParLimits( 1, 5.276, 5.282 );
	func[i]->SetParLimits( 2, 0.001, 0.005 );
	func[i]->FixParameter( 4, 5.289 );
      }else func_set_parameters(sel_fun, func[i], hist[i], xbin, offset+xmin, offset+xmax);
      func[i]->SetLineColor(2);
    }
  }
  

  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();
  c1->cd(1);

  for(Int_t i=Nhist/3-2; i>=0; i-- ){
    c1->cd(3);
    if( i==Nhist/3-2 ) hist[i+2*Nhist/3]->Draw( "hist");
    else               hist[i+2*Nhist/3]->Draw( "hist same" );
    c1->cd(2);
    if( i==Nhist/3-2 ) hist[i+1*Nhist/3]->Draw( "hist" );
    else               hist[i+1*Nhist/3]->Draw( "hist same" );
    c1->cd(1);
    if( i==Nhist/3-2 ) hist[i+0*Nhist/3]->Draw( "hist" );
    else               hist[i+0*Nhist/3]->Draw( "hist same" );
  }

  // rd
  c1->cd(1); hist[2]->SetLineColor(2); hist[2]->SetMarkerColor(2); hist[2]->Draw("PE0same");
  c1->cd(2); hist[5]->SetLineColor(2); hist[5]->SetMarkerColor(2); hist[5]->Draw("PE0same");
  c1->cd(3); hist[8]->SetLineColor(2); hist[8]->SetMarkerColor(2); hist[8]->Draw("PE0same");

  // +++++++ tlegend1 ++++++++++++++++++++++++++++++++++
  TLegend* legend1 = new TLegend( 0.75,0.75,0.99,0.99 );
  legend1->AddEntry( hist[2],"rd",       "P" );
  legend1->AddEntry( hist[1],"gmc.(bb)", "P" );
  legend1->AddEntry( hist[0],"gmc.(qq)", "P" );
  legend1->Draw();


  TH1D** hist_fit = new TH1D*[Nfun];
  for(Int_t k=0; k<Nfun; k++ ){
    hist_fit[k] = new TH1D(*hist[sel_hist[k]]);
    Deco( hist_fit[k], 1, 1, 1 );
  }
  hist_fit[0]->SetTitle( "MC, ee"        );
  hist_fit[1]->SetTitle( "MC, #mu#mu"    );
  hist_fit[2]->SetTitle( "MC, ee+#mu#mu" );
  hist_fit[3]->SetTitle( "RD, ee"        );
  hist_fit[4]->SetTitle( "RD, #mu#mu"    );
  hist_fit[5]->SetTitle( "RD, ee+#mu#mu" );

  for(Int_t k=0; k<Nfun; k++ ){
    c1->cd(k+4);
    if( flag_fit ){
      std::cout << Form( "================================== FUNC%d =================================", k ) << std::endl;
      Double_t* init_var = new Double_t[n_fitfunc_par(sel_fun)];
      for( Int_t m=0; m<n_fitfunc_par(sel_fun); m++ ) init_var[m] = func[k]->GetParameter(m);
      iterative_fit( hist_fit[k], func[k], n_fitfunc_par(sel_fun), init_var, "PE0", manip_func );
    }else hist_fit[k]->Draw();
  }

  c1->Update();
  if( flag_save ) c1->Print( Form("pic/%s_comb_jpsi_s0%s_c1.eps",  axis, stream) );

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++ RooFit ++++++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  Int_t* N0bin = new Int_t[Nfun]; // # of zero-bins for chi2 calculation in RooFit
  
  for( Int_t i=0; i<Nfun; i++ ){
    N0bin[i] = 0;
    for( Int_t k=0; k<xbin; k++ ){
      if( hist_fit[i]->GetBinContent(k+1)==0 ) N0bin[i]++;
    }
    std::cout << Form("<hist%d> ", i) << "# of zero-bin : " << N0bin[i] << std::endl;
  }
  
  std::cout << std::endl
	    << " *************************************************************" << std::endl
	    << " *********************** RooFit Start ************************" << std::endl
	    << " *************************************************************" << std::endl;
  TCanvas* c2 = Canvas( "c2","c2",Nfun/2, 2 );

  // --- Observable ---
  RooRealVar** obs = new RooRealVar*[Nfun];
  for( Int_t i=0; i<Nfun; i++ ) obs[i] = new RooRealVar( axis, xlabel, offset+xmin_fit, offset+xmax_fit );
  Float_t    x_obs;

  // --- Data Set ---
  RooDataSet** data   = new RooDataSet*[Nfun];
  TTree**      chain2 = new TTree*[Nchain]; // apply add_cut[]

  for( Int_t i=0; i<Nfun;   i++ ) data[i]   = new RooDataSet( Form("data%d",i), Form("data%d",i), RooArgSet(*obs[i]) );
  for( Int_t j=0; j<Nchain; j++ ) chain2[j] = chain[j]->GetTree()->CopyTree( add_cut[j] );

    for( Int_t j=0; j<Nchain; j++ ){
    Int_t nevt = chain2[j]->GetTree()->GetEntries();
    std::cout << Form("<chain%d> nevt = %d", j, nevt) << std::endl;
    chain2[j]->GetTree()->SetBranchAddress( axis, &x_obs );
    for( Int_t k=0; k<nevt; k++ ){
      chain2[j]->GetTree()->GetEntry(k);
      if( x_obs>=offset+xmin && x_obs<=offset+xmax ){
	for( Int_t i=0; i<Nfun; i++ ){
	  obs[i]->setVal( x_obs );
	  if( add[sel_hist[i]][j] ) data[i]->add( RooArgSet(*obs[i]) );
	}
      }
    }
  }

  for( Int_t i=0; i<Nfun; i++ ) std::cout << Form("< data%d> nevt = ", i ) << data[i]->numEntries() << std::endl;
  // ------------------------- Fit Fucntion ----------------------------

  const Int_t nparam = 5; // # of fit parameters

  // Sig-PDF(gauss)
  RooRealVar**  sig_mean  = new RooRealVar*[Nfun];
  RooRealVar**  sig_sigma = new RooRealVar*[Nfun];
  RooGaussian** gauss     = new RooGaussian*[Nfun];
  for( Int_t i=0; i<Nfun; i++ ){
    sig_mean [i] = new RooRealVar ( "#mu",    "#mu",    5.279, 5.276, 5.282 );
    sig_sigma[i] = new RooRealVar ( "#sigma", "#sigma", 0.002, 0.001, 0.004 );
    gauss    [i] = new RooGaussian( "gauss", "gauss(x,mean,sigma)", *obs[i], *sig_mean[i], *sig_sigma[i] );
  }

  // Bkg-PDF
  const Double_t endpoint = 5.289;
  RooRealVar** arg_end   = new RooRealVar*[Nfun];
  RooRealVar** arg_shape = new RooRealVar*[Nfun];
  RooArgusBG** argus     = new RooArgusBG*[Nfun];
  for( Int_t i=0; i<Nfun; i++ ){
    arg_end  [i] = new RooRealVar( "end",   "argus endpoint",        endpoint );
    arg_shape[i] = new RooRealVar( "shape", "argus shape parameter", -30.0, -200.0, 0.0 );
    arg_end[i]->setConstant(kTRUE);
    argus[i] = new RooArgusBG( "argus","Argus PDF", *obs[i], *arg_end[i], *arg_shape[i] );
  }

  // Total-PDF
  RooAddPdf**    pdf        = new RooAddPdf*[Nfun];
  RooRealVar**   nsig       = new RooRealVar*[Nfun];
  RooRealVar**   nbkg       = new RooRealVar*[Nfun];
  RooFitResult** fit_result = new RooFitResult*[Nfun];
  for( Int_t i=0; i<Nfun; i++ ){
    nsig[i]       = new RooRealVar ( "N_{sig}", "N_{sig}", 0.50*data[i]->numEntries(), 0.00*data[i]->numEntries(), 1.00*data[i]->numEntries() );
    nbkg[i]       = new RooRealVar ( "N_{bkg}", "N_{bkg}", 0.50*data[i]->numEntries(), 0.00*data[i]->numEntries(), 1.00*data[i]->numEntries() );
    pdf[i]        = new RooAddPdf( Form("pdf%d",i), Form("pdf%d",i), RooArgList(*gauss[i],*argus[i]), RooArgList(*nsig[i], *nbkg[i]) );
    fit_result[i] = pdf[i]->fitTo( *data[i], Extended() );  
  } 
  
  // ------------------------- Draw ----------------------------

  RooPlot**    frame = new RooPlot*  [Nfun];
  TPaveText**  box   = new TPaveText*[Nfun];

  for( Int_t i=0; i<Nfun; i++ ){
    c2->cd(i+1);
    frame[i] = obs[i]->frame();
    frame[i]->GetXaxis()->CenterTitle();
    frame[i]->GetYaxis()->CenterTitle();
    frame[i]->SetTitleOffset( 1.00,"x" );
    frame[i]->SetTitleOffset( 1.30,"y" );
    frame[i]->SetTitle( hist_fit[i]->GetTitle() );

    Double_t tmp_scale = 1;
    if( i==0 || i==1 || i==2 ) tmp_scale /= used_nstream;
    data[i]->plotOn ( frame[i], Binning(xbin),       LineWidth(1), Rescale( tmp_scale ) );
    pdf[i] ->plotOn ( frame[i], ProjWData(*data[i]), LineWidth(1), Normalization( tmp_scale, 1 ) );
    box[i] = new TPaveText( 0.05, 0.87, 0.35, 0.93,"BRNDC" );
    box[i]->AddText( Form("#chi^{2}/ndf = %f", frame[i]->chiSquare(nparam+N0bin[i])) );
    pdf[i] ->plotOn ( frame[i], ProjWData(*data[i]), Components(*gauss[i]), LineStyle(7), LineColor(2), LineWidth(1), Normalization( tmp_scale, 1 ) );
    pdf[i] ->plotOn ( frame[i], ProjWData(*data[i]), Components(*argus[i]), LineStyle(7), LineColor(4), LineWidth(1), Normalization( tmp_scale, 1 ) );
    pdf[i] ->paramOn( frame[i], Format("NELU", AutoPrecision(2)), Layout(0.40, 0.99, 0.99) );

    frame[i]->addObject( box[i] );
    frame[i]->SetAxisRange(0,1.4*hist_fit[i]->GetMaximum(),"Y");
    frame[i]->Draw();
  }


  c2->Update();
  if( flag_save ) c2->Print( Form("pic/%s_comb_jpsi_s0%s_c2.eps",  axis, stream) );
  
  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();

  delete[] chain;
  delete[] tmphist;
  delete[] hist;
  delete[] hist_fit;
  delete[] func;
  delete[] add_cut;
  delete   c1;

  delete[] chain2;
  delete[] obs;
  delete[] data;
  delete[] sig_mean;
  delete[] sig_sigma;
  delete[] arg_end;
  delete[] arg_shape;
  delete[] gauss;
  delete[] argus;
  delete[] nsig;
  delete[] nbkg;
  delete[] pdf;
  delete[] fit_result;
  delete[] frame;
  delete[] box;
  delete   c2;

  return 0;
}

