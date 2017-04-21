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

using namespace Mbc_bkg;

const Bool_t flag_save = true; // outfile.eps and outfile.root

const Double_t endpoint     = 5.289;

const Bool_t   flag_k4pi     = !true; // 1(veto  K4pi    modes)
const Bool_t   flag_unflavor = !true; // 1(veto unflavor modes)
const Bool_t   flag_ccpi0    = true; // 1(veto ccpi0 peak    )

void manip_func( TF1* func ){
  func->FixParameter( 1, 5.2794  );
  func->FixParameter( 2, 0.00273 );
  func->FixParameter( 4, endpoint );
}

Int_t main( Int_t argc, Char_t** argv ){

  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==4 || argc==5) ) std::cerr << "wrong input" << std::endl
					<< " Usage : ./draws_1d_Mbc_fit (int)fl_mode (int)fl_mc (char*)setname"
					<< std::endl, abort();
  const Int_t   fl_mode            = atoi(argv[1]);
  const Int_t   fl_mc              = atoi(argv[2]);
  const Char_t* setname            = argv[3];
  Int_t         fl_appRun          = 1;
  if( argc==5 ) fl_appRun = atoi( argv[4] );
  if ( abs(fl_mode)!=10 && abs(fl_mode)!=1000 ) std::cerr << "wrong fl_mc : " << fl_mc << std::endl, abort();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t   Nchain             = 2; // [ee,mm]
  const Int_t   Nhist              = 3; // [ee,mm,ee+mm]
  //const Char_t* infile[2]          = {"~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut3/rd_522/RD_",
  //"~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut3/xsjpsi_522/sigMC_*_set[%s]"};
  const Char_t* infile[2]          = {"~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut2/rd_522/RD_",
				      "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut2/sig_522/sigMC_*_set[%s]"};
  const Int_t   fl_mode_ll[Nchain] = {1,0};
  const Int_t   sel_fun            = 15;
  // +++++++ add cut +++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  Char_t* add_cut = new Char_t[2048];
  strcpy( add_cut, (char*)sTmp.str().c_str() );
  if     ( fl_mode ==    10 ) sTmp << " (rm_xs==10 || rm_xs==110 || rm_xs==1010 || rm_xs==210 || rm_xs==1110 || rm_xs==310 || rm_xs==1210 || rm_xs==410 || rm_xs==1310)"; // mode     including Ks
  else if( fl_mode ==   -10 ) sTmp << "!(rm_xs==10 || rm_xs==110 || rm_xs==1010 || rm_xs==210 || rm_xs==1110 || rm_xs==310 || rm_xs==1210 || rm_xs==410 || rm_xs==1310)"; // mode not including Ks
  else if( fl_mode ==  1000 ) sTmp << "rm_xs > 999"; // mode     including pi0
  else if( fl_mode == -1000 ) sTmp << "rm_xs < 999"; // mode not including pi0

  sTmp << " && xs_m >1.1";
  if( flag_k4pi     ) sTmp << " && " << makeCut_5body_veto().c_str();
  if( flag_unflavor ) sTmp << " && " << makeCut_unflavor_veto().c_str();
  if( flag_ccpi0    ) sTmp << " && " << cut_ccpi0;
  strcpy( add_cut, (char*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();
  std::cout << "add_cut : " << add_cut << std::endl;
  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++
  MChain** chain   = new MChain*[Nchain];
  std::cout << Form(" ************************ make chain-tree ( %s, %s ) *************************************",tname,axis) << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<infile %d > ", j );
    if( fl_mc==0 ) chain[j] = new MChain( infile[fl_mc],               tname, branch_table(), 0, "*.root" );
    else           chain[j] = new MChain( Form(infile[fl_mc],setname), tname, branch_table(), 0, "*.root" );
    nominal_cut_selection( chain[j], fl_mode_ll[j] )( chain[j]->GetCut(), tname );
  }
  
  // ++++++++++++++++++++++++
  // cut change
  for( Int_t j=0; j<Nchain; j++ ){
    chain[j]->GetCut()->Set( "Mbc", 1,  5.20, 0.0, 5.30 );
    chain[j]->GetCut()->Display(0);
    chain[j]->MakeTree();
    std::cout << "[hist] " << j << " : " << chain[j]->GetTree()->GetEntries() << std::endl;
  }

  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  TH1D** hist = new TH1D* [Nhist];
  std::cout << std::endl
	    << " ************************ make hist *************************************" << std::endl;
  for( Int_t i=0; i<Nhist; i++ ){
    if(      i==0 ) hist[i] = new TH1D( Form("hist%d",i), Form("fl_mode=%d, ee",        fl_mode), xbin,offset+xmin,offset+xmax );
    else if( i==1 ) hist[i] = new TH1D( Form("hist%d",i), Form("fl_mode=%d, #mu#mu",    fl_mode), xbin,offset+xmin,offset+xmax );
    else if( i==2 ) hist[i] = new TH1D( Form("hist%d",i), Form("fl_mode=%d, ee+#mu#mu", fl_mode), xbin,offset+xmin,offset+xmax );
    if( i!=2 ) chain[i]->GetTree()->Project( Form("hist%d",i), axis, add_cut );
    else{
      hist[i]-> Add( hist[0] );
      hist[i]-> Add( hist[1] );
    }
    if( i!=2 ) hist[i]->Sumw2();
    hist[i]->GetXaxis()->CenterTitle();
    hist[i]->GetYaxis()->CenterTitle();
    hist[i]->SetXTitle( xlabel );
    std::cout << "[test] " << i << " : " << hist[i]->GetEntries() << std::endl;
  }

  // +++++++ make fit-func ++++++++++++++++++++++++++++++++++
  TF1** func = new TF1*   [Nhist];
  for( Int_t i=0; i<Nhist; i++ ){
    func[i] = new TF1( Form("func%d",i), make_func(sel_fun),offset+xmin_fit,offset+xmax_fit,n_fitfunc_par(sel_fun) );
    Double_t mu    = 5.2793;
    Double_t sigma = 0.00262;
    Double_t area  = ( hist[i]->GetMaximum() - hist[i]->GetBinContent(1) )*sqrt( TMath::TwoPi() )*sigma;
    func[i]->SetParNames  ( "area","mean","sigma","norm",                         "Ebeam", "shape" );
    func[i]->SetParameters( area,   mu,    sigma,  2*hist[i]->GetBinContent(1),  endpoint,   -15   );
    func[i]->FixParameter( 1, mu       );
    func[i]->FixParameter( 2, sigma    );
    func[i]->FixParameter( 4, endpoint );
    func[i]->SetLineColor(2);
  }

  // +++++++ draw ++++++++++++++++++++++++++++++++++
  TCanvas*  c1 = Canvas( "c1","c1", 3, 1 );
  c1->Draw();
  
  for(Int_t i=0; i<Nhist; i++ ){
    c1->cd(i+1);
    std::cout << Form( "================================== FUNC%d =================================", i ) << std::endl;
    Double_t init_var[n_fitfunc_par(sel_fun)];
    for( Int_t m=0; m<n_fitfunc_par(sel_fun); m++ ) init_var[m] = func[i]->GetParameter(m);
    iterative_fit( hist[i], func[i], n_fitfunc_par(sel_fun), init_var, "PE0", manip_func );
    std::cout << std::setw( 3) << std::right << i
	      << std::setw( 3) << std::right << fl_mc
	      << std::setw( 6) << std::right << fl_mode
	      << std::setw(12) << std::right << func[i]->GetParameter(0)
	      << std::setw(12) << std::right << func[i]->GetParError (0)
	      << "   HOGE"
	      << std::endl;
      
  }
  
  c1->Update();
  c1->Print( Form("pic/Mbc_fit_mode_mode%d_mc%d_set%s.eps", fl_mode, fl_mc, setname) );
  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();
  
  return 0;
}

