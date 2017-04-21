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
#include <TMath.h>
#include <Math/DistFunc.h>

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
#include <RooGenericPdf.h>
#include "RooCFunction1Binding.h"
#include "RooCFunction2Binding.h"
#include "RooCFunction3Binding.h"
#include "RooCFunction4Binding.h" 
#include "RooTFnBinding.h" 

void manip_func( TF1* func ){
  if( func->GetNpar()==3 || func->GetNpar()==4 ){ // (modified-)argus
    func->FixParameter( 1, 5.289 ); 
  }else if( func->GetNpar()==6 || func->GetNpar()==7 ){ // gaussian + (modified-)argus
    func->FixParameter( 1, 5.2794  );
    func->FixParameter( 2, 0.00273 );
    func->FixParameter( 4, 5.289   );
  }
}

Int_t main( Int_t argc, Char_t** argv ){
  //using namespace sig_gmc_rd_emu_beforebgsup;
  using namespace sig_gmc_rd_emu;
  
  using namespace Mbc_bkg_wide;
  using namespace RooFit;
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==5 || argc==6) ) std::cerr << "wrong input" << std::endl
    					<< " Usage : ./draws_1d_Mbc (char*)stream (double)used_nstream (int)fl_mode_ll (int)sel_fun [(int)fl_appRun]" << std::endl
					<< std::endl, abort();
  Char_t*  stream       = argv[1];
  Double_t used_nstream = atof(argv[2]);
  Int_t    fl_mode_ll   = atoi(argv[3]);
  Int_t    sel_fun      = atoi(argv[4]); // 15(gauss+argus), 151(gauss+modified-argus), 50(argus), 51(modified-argus)
  Int_t   fl_appRun     = 1;
  if( argc==6 ) fl_appRun = atoi( argv[5] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Bool_t flag_k4pi      = true; // 1(veto  K4pi    modes)
  const Bool_t flag_unflavor  = true; // 1(veto unflavor modes)
  const Bool_t flag_save      = true; // outfile.eps and outfile.root
  const Bool_t flag_fit       = !true;

  const Int_t  Nfun           = 2;    // [MC,RD]
  const Int_t  sel_hist[Nfun] = {1,2}; // MC, RD
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t    Nchain             = 3; // [gmc(qq,bb),rd]
  const Int_t    Nhist              = 3; // [gmc(qq,bb),rd]
  const Int_t    nfile[Nchain]      = {0};
  const Int_t    fl_sb[Nchain]      = {0,0,2}; // 0(bkg), 1(sig), 2(rd)
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
    //if     ( fl_sb[i]==0 ) sTmp << indir[fl_sb[i]] << "gMC_*_s0[" << stream << "]"; // bkg(before B.G. suppression)
    //else if( fl_sb[i]==2 ) sTmp << indir[fl_sb[i]] << "RD_";                        // rd (before B.G. suppression)
    if     ( fl_sb[i]==0 ) sTmp << indir[fl_sb[i]] << "gMC_*_s0[" << stream << Form("]_*_emu%d", fl_mode_ll); // bkg(after B.G. suppression)
    else if( fl_sb[i]==2 ) sTmp << indir[fl_sb[i]] << Form("RD_*emu%d", fl_mode_ll);                          // rd (after B.G. suppression)
    strcpy( infile[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t add[Nhist][Nchain] ={
    {1},     // gmc(qq)
    {1,1},   // gmc(bb)
    {0,0,1}, // rd
  };

  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    add_cut[i] = new Char_t[4096];
    sTmp << LRNB_cut;
    if( flag_k4pi     ) sTmp << " && " << makeCut_5body_veto().c_str();
    if( flag_unflavor ) sTmp << " && " << makeCut_unflavor_veto().c_str();

    if     ( i==0 ) sTmp << " && genbfl==0"; // gmc(qq)
    else if( i==1 ) sTmp << " && genbfl!=0"; // gmc(bb)
    
    strcpy( add_cut[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain   = new MChain*[Nchain];
  TH1D**    tmphist = new TH1D*  [Nchain];
  TH1D**    hist    = new TH1D*  [Nhist];
  TF1**     func    = new TF1*   [Nfun];
  TCanvas*  c1      = Canvas( "c1","c1",3, 1 );
  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++
  std::cout << Form(" ************************ make tmphist ( %s, %s ) *************************************",tname,axis) << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<infile %d > ", j );
    chain[j] = new MChain( infile[j], tname, branch_table(), nfile[j], tail );
    //nominal_cut_selection( chain[j], fl_mode_ll )( chain[j]->GetCut(), tname );
    nominal_cut_selection( chain[j], 1 )( chain[j]->GetCut(), tname );
  }

  // ++++++++++++++++++++++++
  // cut change
  for( Int_t j=0; j<Nchain; j++ ){
    chain[j]->GetCut()->Set( "Mbc", 1,  5.20, 0.0, 5.30 );
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
    chain[j]->GetTree()->Project( Form("tmphist%d",j), axis, add_cut[j] );
    tmphist[j]->Sumw2();
    if( fl_sb[j]==0 ) tmphist[j]->Scale( 1/scale_event_bkg );
  }
  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make hist *************************************" << std::endl;
  for( Int_t i=0; i<Nhist; i++ ){
    hist[i] = new TH1D( Form("hist%d",i),Form("hist%d",i), xbin,offset+xmin,offset+xmax );
    Deco( hist[i], 3, col_fil[i%3], col_fil[i%3] );
    for( Int_t j=0; j<Nchain; j++ ){
      if( add[i][j] ) hist[i]->Add( tmphist[j] );
    }
  }
  hist[0]->SetTitle("gMC(qq)"   );
  hist[1]->SetTitle("gMC(bb)"   );
  hist[2]->SetTitle("RD"        );

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
      if( sel_fun==50 ){
	func[i]->SetParNames  ("norm","Ebeam","shape");
	func[i]->SetParameters( 100*hist[i]->GetBinContent(1), 5.289, -30);
	func[i]->FixParameter( 1, 5.289 );
      }else if( sel_fun==51 ){
	func[i]->SetParNames  ("norm","Ebeam","shape","new");
	func[i]->SetParameters( 100*hist[i]->GetBinContent(1), 5.289, -30, 0.5);
	func[i]->FixParameter( 1, 5.289 );
      }else if( sel_fun==15 ){
	Double_t mu    = 5.279;
	Double_t sigma = 0.003;
	Double_t area  = ( hist[i]->GetMaximum() - hist[i]->GetBinContent(1) )*sqrt( TMath::TwoPi() )*sigma;
	func[i]->SetParNames  ( "area","mean","sigma","norm",                          "Ebeam","shape" );
	func[i]->SetParameters( area,   mu,    sigma,  100*hist[i]->GetBinContent(1),  5.289,   -30 );
	//func[i]->SetParLimits( 1, 5.276, 5.282 );
	//func[i]->SetParLimits( 2, 0.001, 0.005 );
	func[i]->FixParameter( 1, 5.2794 );
	func[i]->FixParameter( 2, 0.00273 );
	func[i]->FixParameter( 4, 5.289 );
      }else if( sel_fun==151 ){
	Double_t mu    = 5.279;
	Double_t sigma = 0.003;
	Double_t area  = ( hist[i]->GetMaximum() - hist[i]->GetBinContent(1) )*sqrt( TMath::TwoPi() )*sigma;
	func[i]->SetParNames  ( "area","mean","sigma","norm",                          "Ebeam","shape", "new" );
	func[i]->SetParameters( area,   mu,    sigma,  100*hist[i]->GetBinContent(1),  5.289,   -30,     0.5 );
	//func[i]->SetParLimits( 1, 5.276, 5.282 );
	//func[i]->SetParLimits( 2, 0.001, 0.005 );
	func[i]->FixParameter( 1, 5.2794 );
	func[i]->FixParameter( 2, 0.00273 );
	func[i]->FixParameter( 4, 5.289 );
      }else func_set_parameters(sel_fun, func[i], hist[i], xbin, offset+xmin, offset+xmax);
      func[i]->SetLineColor(2);
    }
  }
  
  
  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();
  c1->cd(1);

  hist[1]->Draw( "hist");
  hist[0]->Draw( "hist same" );
  // rd
  c1->cd(1); hist[2]->SetLineColor(2); hist[2]->SetMarkerColor(2); hist[2]->Draw("PE0same");


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
  hist_fit[0]->SetTitle( "MC" );
  hist_fit[1]->SetTitle( "RD" );
  
  for(Int_t k=0; k<Nfun; k++ ){
    c1->cd(k+2);
    if( flag_fit ){
      std::cout << Form( "================================== FUNC%d =================================", k ) << std::endl;
      Double_t init_var[n_fitfunc_par(sel_fun)];
      for( Int_t m=0; m<n_fitfunc_par(sel_fun); m++ ) init_var[m] = func[k]->GetParameter(m);
      iterative_fit( hist_fit[k], func[k], n_fitfunc_par(sel_fun), init_var, "PE0", manip_func );
    }else hist_fit[k]->Draw();
  }

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

  c1->Update();
  if( flag_save ) c1->Print( Form("pic/%s_comb_emu_lep%d_func%d_s0%s_c1.eps",  axis, fl_mode_ll, sel_fun, stream) );

  if( flag_fit ){
    std::cout << std::endl
	      << " *************************************************************" << std::endl
	      << " *********************** RooFit Start ************************" << std::endl
	      << " *************************************************************" << std::endl;
    TCanvas* c2 = Canvas( "c2","c2",2, 1 );
    
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
      sig_mean [i] = new RooRealVar ( "#mu",    "#mu",    5.2794,  5.276, 5.282 );
      sig_sigma[i] = new RooRealVar ( "#sigma", "#sigma", 0.00273, 0.001, 0.004 );
      gauss    [i] = new RooGaussian( "gauss", "gauss(x,mean,sigma)", *obs[i], *sig_mean[i], *sig_sigma[i] );
      sig_mean [i]->setConstant(kTRUE);
      sig_sigma[i]->setConstant(kTRUE);
    }
    
    // Bkg-PDF
    const Double_t endpoint = 5.289;
    RooRealVar**    arg_end   = new RooRealVar*[Nfun];
    RooRealVar**    arg_shape = new RooRealVar*[Nfun];
    RooRealVar**    arg_new   = new RooRealVar*[Nfun];
    RooAbsPdf**     modargus  = new RooAbsPdf*[Nfun];

    for( Int_t i=0; i<Nfun; i++ ){
      arg_end  [i] = new RooRealVar( "end",   "argus endpoint",        endpoint, 5.280, 5.300 );
      arg_shape[i] = new RooRealVar( "shape", "argus shape parameter", -30.0,   -200.0, 0.0   );
      arg_new  [i] = new RooRealVar( "new",   "argus sqrt parameter",    0.5,      0.0, 0.6   );
      arg_end[i]->setConstant(kTRUE);
      if     ( sel_fun==50 || sel_fun== 15 ) modargus[i] = bindPdf( Form("modargus%d",i), func_roo_argus,    *obs[i], *arg_end[i], *arg_shape[i]);
      else if( sel_fun==51 || sel_fun==151 ) modargus[i] = bindPdf( Form("modargus%d",i), func_roo_modargus, *obs[i], *arg_end[i], *arg_shape[i], *arg_new[i]  );
    }
    
    // Total-PDF
    RooAddPdf**    pdf        = new RooAddPdf*[Nfun];
    RooRealVar**   nsig       = new RooRealVar*[Nfun];
    RooRealVar**   nbkg       = new RooRealVar*[Nfun];
    RooFitResult** fit_result = new RooFitResult*[Nfun];
    for( Int_t i=0; i<Nfun; i++ ){
      nsig[i] = new RooRealVar ( "N_{sig}", "N_{sig}", 0.50*data[i]->numEntries(), 0.00*data[i]->numEntries(), 1.00*data[i]->numEntries() );
      nbkg[i] = new RooRealVar ( "N_{bkg}", "N_{bkg}", 0.50*data[i]->numEntries(), 0.00*data[i]->numEntries(), 1.00*data[i]->numEntries() );
      pdf[i]  = new RooAddPdf( Form("pdf%d",i), Form("pdf%d",i), RooArgList(*gauss[i],*modargus[i]), RooArgList(*nsig[i], *nbkg[i]) );

      if     ( sel_fun==50 || sel_fun== 51 ) fit_result[i] = modargus[i]->fitTo( *data[i]             );
      else if( sel_fun==15 || sel_fun==151 ) fit_result[i] = pdf[i]     ->fitTo( *data[i], Extended() );
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
      if( i==0 ) tmp_scale /= used_nstream;
      data[i]->plotOn ( frame[i], Binning(xbin),       LineWidth(1), Rescale( tmp_scale ) );

      if     ( sel_fun==15 || sel_fun==151 ) pdf[i]     ->plotOn( frame[i], ProjWData(*data[i]), LineWidth(1), Normalization( tmp_scale, 1 ) );
      else if( sel_fun==50 || sel_fun== 51 ) modargus[i]->plotOn( frame[i], ProjWData(*data[i]), LineWidth(1), Normalization( tmp_scale, 1 ) );
      box[i] = new TPaveText( 0.05, 0.87, 0.35, 0.93,"BRNDC" );
      box[i]->AddText( Form("#chi^{2}/ndf = %f", frame[i]->chiSquare(nparam+N0bin[i])) );
      
      if( sel_fun==15 || sel_fun==151 ){
	pdf[i] ->plotOn( frame[i], ProjWData(*data[i]), Components(*gauss[i]), LineStyle(7), LineColor(2), LineWidth(1), Normalization( tmp_scale, 1 ) );
	pdf[i] ->plotOn( frame[i], ProjWData(*data[i]), Components(*modargus[i]), LineStyle(7), LineColor(4), LineWidth(1), Normalization( tmp_scale, 1 ) );
	pdf[i] ->paramOn( frame[i], Format("NELU", AutoPrecision(2)), Layout(0.40, 0.99, 0.99) );
      }else if( sel_fun==50 || sel_fun==51 ) modargus[i]->paramOn( frame[i], Format("NELU", AutoPrecision(2)), Layout(0.40, 0.99, 0.99) );
      
      frame[i]->addObject( box[i] );
      frame[i]->SetAxisRange(0,1.4*hist_fit[i]->GetMaximum(),"Y");
      frame[i]->Draw();
    }
    
    c2->Update();
    if( flag_save ) c2->Print( Form("pic/%s_comb_emu_lep%d_func%d_s0%s_c2.eps",  axis, fl_mode_ll, sel_fun, stream) );
  }
  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();
  
  delete[] chain;
  delete[] tmphist;
  delete[] hist;
  delete[] hist_fit;
  delete[] func;
  delete[] add_cut;
  delete   c1;
  
  return 0;
}

