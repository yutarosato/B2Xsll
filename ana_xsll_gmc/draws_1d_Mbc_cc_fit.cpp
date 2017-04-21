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

Int_t main( Int_t argc, Char_t** argv ){
  using namespace RooFit;
  using namespace gmc;
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==2 || argc==3) ) std::cerr << "wrong input" << std::endl
					<< " Usage : ./draws_1d_Mbc_fit (char*)stream [(int)fl_appRun]"
					<< std::endl, abort();

  Char_t* stream    = argv[1];
  Int_t   fl_appRun = 1;
  if( argc==3 ) fl_appRun = atoi( argv[2] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  sTmp << indir << "/gMC_*_s0[" << stream << "]";
  Char_t tmp_infile[255];
  strcpy( tmp_infile, (Char_t*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();
  const Char_t* infile = (const Char_t*)tmp_infile;

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
  
  Char_t** add_cut = new Char_t*[Nchain];
  add_cut[0] = "3.096916-0.40<cc_m && 3.096916+0.15>cc_m && 3.096916-0.40<cc_morg && 3.096916+0.15>cc_morg"; // ee for J/psi samples
  add_cut[1] = "3.096916-0.25<cc_m && 3.096916+0.10>cc_m && 3.096916-0.25<cc_morg && 3.096916+0.10>cc_morg"; // mm for J/psi samples
  
  using namespace Mbc;
  //using namespace deltaE;
  //using namespace M_Xs;
  //using namespace M_ll;
  //using namespace M_Xs_gen;

  const Bool_t flag_save = true; // outfile.eps and outfile.root
  const Bool_t  flag_fit = true;
  const Int_t   sel_fun  = 15; // 10(gauss), 15(gauss+argus)
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain   = new MChain*[Nchain];
  TH1D**    tmphist = new TH1D*  [Nchain];
  TH1D**    hist    = new TH1D*  [Nhist];
  TF1**     func    = new TF1*   [Nhist];
  TCanvas*  c1      = Canvas( "c1","c1",Nhist, 1 );

  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++
  std::cout << Form(" ************************ make chain-tree ( %s, %s ) *************************************",tname,axis) << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<infile %d > ", j );
    chain[j] = new MChain( infile, tname, branch_table(), nfile[j], tail );
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
	Double_t sigma = 4*hist[i]->GetBinWidth(1);
	Double_t area  = ( hist[i]->GetMaximum() - (a*mu+b) )*sqrt( TMath::TwoPi() )*sigma;
	func[i]->SetParNames  ( "area","mean","sigma","norm","Ebeam","shape" );
	func[i]->SetParameters( area,   mu,    sigma,  500,   5.289,   -60 );
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
      hist[i]->Fit(func[i],"R","PE0");
      std::cout << "chi2/NDF = " << func[i]->GetChisquare() << " / " << func[i]->GetNDF()
		<< " = " << func[i]->GetChisquare()/func[i]->GetNDF()
		<< std::endl;
    }else hist[i]->Draw( "same" );
  }

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++ RooFit ++++++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  Int_t* N0bin = new Int_t[Nhist]; // # of zero-bins for chi2 calculation
  
  for( Int_t i=0; i<Nhist; i++ ){
    N0bin[i] = 0;
    for( Int_t k=0; k<xbin; k++ ){
      if( hist[i]->GetBinContent(k+1)==0 ) N0bin[i]++;
    }
    std::cout << Form("<hist%d> ", i) << "# of zero-bin : " << N0bin[i] << std::endl;
  }
  
  std::cout << std::endl
	    << " *************************************************************" << std::endl
	    << " *********************** RooFit Start ************************" << std::endl
	    << " *************************************************************" << std::endl;
  TCanvas* c2 = Canvas( "c2","c2",Nhist, 1 );
  // --- Observable ---
  RooRealVar** obs = new RooRealVar*[Nhist];
  for( Int_t i=0; i<Nhist; i++ ) obs[i] = new RooRealVar( axis, xlabel, offset+xmin_fit, offset+xmax_fit );
  Float_t    x_obs;
  Float_t    cc_m, cc_morg, rm_l;
  // --- Data Set ---
  RooDataSet** data   = new RooDataSet*[Nhist];
  TTree**      chain2 = new TTree*[Nhist]; // apply add_cut[]

  for( Int_t i=0; i<Nhist;  i++ ) data[i]   = new RooDataSet( Form("data%d",i), Form("data%d",i), RooArgSet(*obs[i]) );
  for( Int_t j=0; j<Nchain; j++ ) chain2[j] = chain[j]->GetTree()->CopyTree( add_cut[j] );

    for( Int_t j=0; j<Nchain; j++ ){
    Int_t nevt = chain2[j]->GetTree()->GetEntries();
    std::cout << Form("<chain%d> nevt = %d", j, nevt) << std::endl;
    chain2[j]->GetTree()->SetBranchAddress( axis,      &x_obs   );
    for( Int_t k=0; k<nevt; k++ ){
      chain2[j]->GetTree()->GetEntry(k);
      if( x_obs>=offset+xmin && x_obs<=offset+xmax ){
	for( Int_t i=0; i<Nhist; i++ ){
	  obs[i]->setVal( x_obs );
	  if( add[i][j] ) data[i]->add( RooArgSet(*obs[i]) );
	}
      }
    }
  }

  for( Int_t i=0; i<Nhist; i++ ) std::cout << Form("< data%d> nevt = ", i ) << data[i]->numEntries() << std::endl;
  // ------------------------- Fit Fucntion ----------------------------

  const Int_t nparam = 5; // # of fit parameters

  // Sig-PDF(gauss)
  RooRealVar**  sig_mean  = new RooRealVar*[Nhist];
  RooRealVar**  sig_sigma = new RooRealVar*[Nhist];
  RooGaussian** gauss     = new RooGaussian*[Nhist];
  for( Int_t i=0; i<Nhist; i++ ){
    sig_mean [i] = new RooRealVar ( "#mu",    "#mu",    5.280, 5.275, 5.285 );
    sig_sigma[i] = new RooRealVar ( "#sigma", "#sigma", 0.001, 0.000, 0.005 );
    gauss    [i] = new RooGaussian( "gauss", "gauss(x,mean,sigma)", *obs[i], *sig_mean[i], *sig_sigma[i] );
  }

  // Bkg-PDF
  const Double_t endpoint = 5.289;
  RooRealVar** arg_end   = new RooRealVar*[Nhist];
  RooRealVar** arg_shape = new RooRealVar*[Nhist];
  RooArgusBG** argus     = new RooArgusBG*[Nhist];
  for( Int_t i=0; i<Nhist; i++ ){
    arg_end  [i] = new RooRealVar( "end",   "argus endpoint",        endpoint,  5.41, 5.417 );
    arg_shape[i] = new RooRealVar( "shape", "argus shape parameter", -50.0,   -200.0, 0.0   );
    arg_end[i]->setConstant(kTRUE);
    argus[i] = new RooArgusBG( "argus","Argus PDF", *obs[i], *arg_end[i], *arg_shape[i] );
  }

  // Total-PDF
  RooAddPdf**    pdf        = new RooAddPdf*[Nhist];
  RooRealVar**   nsig       = new RooRealVar*[Nhist];
  RooRealVar**   nbkg       = new RooRealVar*[Nhist];
  RooFitResult** fit_result = new RooFitResult*[Nhist];
  for( Int_t i=0; i<Nhist; i++ ){
    nsig[i]       = new RooRealVar ( "N_{sig}", "N_{sig}", 0.70*data[i]->numEntries(), 0.60*data[i]->numEntries(), 0.80*data[i]->numEntries() );
    nbkg[i]       = new RooRealVar ( "N_{bkg}", "N_{bkg}", 0.30*data[i]->numEntries(), 0.20*data[i]->numEntries(), 0.40*data[i]->numEntries() );
    pdf[i]        = new RooAddPdf( Form("pdf%d",i), Form("pdf%d",i), RooArgList(*gauss[i],*argus[i]), RooArgList(*nsig[i], *nbkg[i]) );
    fit_result[i] = pdf[i]->fitTo( *data[i], Extended() );  
  } 

  
  // ------------------------- Draw ----------------------------

  RooPlot**    frame       = new RooPlot*  [Nhist];
  TPaveText**  box         = new TPaveText*[Nhist];
  //RooAbsReal** int_pdf_bin = new RooAbsReal*[Nhist];
  //Double_t*    cnt_chi2    = new Double_t[Nhist];
  //Int_t*       cnt_ndf     = new Int_t[Nhist];

  for( Int_t i=0; i<Nhist; i++ ){
    c2->cd(i+1);
    frame[i] = obs[i]->frame();
    frame[i]->GetXaxis()->CenterTitle();
    frame[i]->GetYaxis()->CenterTitle();
    frame[i]->SetTitleOffset( 1.00,"x" );
    frame[i]->SetTitleOffset( 1.30,"y" );
    if     ( i==0 ) frame[i]->SetTitle( "ee"        );
    else if( i==1 ) frame[i]->SetTitle( "#mu#mu   " );
    else if( i==2 ) frame[i]->SetTitle( "ee+#mu#mu" );
    data[i]->plotOn ( frame[i], Binning(xbin), LineWidth(1) );
    pdf[i] ->plotOn ( frame[i], ProjWData(*data[i]), LineWidth(1) );
    box[i] = new TPaveText( 0.05, 0.87, 0.35, 0.93,"BRNDC" );
    box[i]->AddText( Form("#chi^{2}/ndf = %f", frame[i]->chiSquare(nparam+N0bin[i])) );
    pdf[i] ->plotOn ( frame[i], ProjWData(*data[i]), Components(*gauss[i]), LineStyle(7), LineColor(2), LineWidth(1) );
    pdf[i] ->plotOn ( frame[i], ProjWData(*data[i]), Components(*argus[i]), LineStyle(7), LineColor(4), LineWidth(1) );
    pdf[i] ->paramOn( frame[i], Format("NELU", AutoPrecision(2)), Layout(0.40, 0.99, 0.99) );


    /* -------------- calculation of chi2 -------------
    cnt_chi2[i] = 0;
    cnt_ndf [i] = 0;
    int_pdf_bin[i] = pdf[i]->createIntegral(*obs[i], NormSet(*obs[i]), Range("bin_int") );
    for( Int_t k=0; k<xbin; k++ ){
      obs[i]->setRange("bin_int", hist[i]->GetBinLowEdge(k+1), hist[i]->GetBinLowEdge(k+1) + hist[i]->GetBinWidth(k+1) );
      if( hist[i]->GetBinContent(k+1)==0 ) continue;
      Double_t tmp_chi = hist[i]->GetBinContent(k+1) - (nsig[i]->getVal()+nbkg[i]->getVal()) * int_pdf_bin[i]->getVal();
      tmp_chi = tmp_chi*tmp_chi;
      tmp_chi /= hist[i]->GetBinError(k+1)*hist[i]->GetBinError(k+1);
      cnt_chi2[i] += tmp_chi;
      cnt_ndf[i]++;
    }
    */

    //box[i]->AddText( Form("#chi^{2}/ndf = %f", cnt_chi2[i]/(Double_t)(cnt_ndf[i]-nparam)) );
    frame[i]->addObject( box[i] );
    frame[i]->Draw();
  }

  c1->Update();
  c2->Update();

  if( flag_save ){
    c1->Print( Form("pic/%s_cc_fit_s0%s_c1.eps", axis, stream) );
    c1->Print( Form("pic/%s_cc_fit_s0%s_c1.root",axis, stream) );
    c2->Print( Form("pic/%s_cc_fit_s0%s_c2.eps", axis, stream) );
    c2->Print( Form("pic/%s_cc_fit_s0%s_c2.root",axis, stream) );
  }
  
  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();

  delete[] chain;
  delete[] tmphist;
  delete[] hist;
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

