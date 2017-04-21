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

const Bool_t flag_save = true; // outfile.eps and outfile.root
const Bool_t flag_fit  = true;
const Bool_t flag_roof = true; // flase(skip roofit)

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
  if( !(argc==7 || argc==8) ) std::cerr << "wrong input" << std::endl
					<< " Usage : ./draws_1d_Mbc_fit (char*)stream (char*)setname (double)nstream (double)nset (int)sel_fun (int)incl_ratio [(int)fl_appRun]" << std::endl
					<< " [sel_fun] : 50(argus), 51(modified-argus), 15(gauss+argus), 151(gauss+modified-argus)"
					<< std::endl, abort();
  Char_t*  stream       = argv[1];
  Char_t*  setname      = argv[2];
  Double_t used_nstream = atof(argv[3]);
  Double_t used_nset    = atof(argv[4]);
  Int_t    sel_fun      = atoi(argv[5]);
  Int_t    incl_ratio   = atoi(argv[6]); // [%]
  Int_t    fl_appRun    = 1;
  if( argc==8 ) fl_appRun = atoi( argv[7] );

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ///*
  unsigned int seed = time(NULL);
  seed += (1001*incl_ratio);
  seed = (seed>>16 | seed<<16 );
  seed += (2003*incl_ratio);
  gRandom->SetSeed( seed );
  //*/
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t Nchain        = 2*3; // [gmc,sigmc,rd] x [ee,mm]
  const Int_t Ncategory     = 8; //[gmc(qq,non-peak,peak),sigmc((f,f), (f,t), (t,f), (t,t), true)]
  const Int_t Nplot         = 4; //[bkg(qq,non-peak,peak),sigmc]
  const Int_t Ntmp          = Ncategory*2; // x[ee,mm      ]
  const Int_t Nhist         = Nplot    *3; // x[ee,mm,ee+mm]
  const Int_t fl_sb[Ntmp]      = {0,0,0,1,1,1,1,1,
				  3,3,3,4,4,4,4,4}; // 0(bkg,ee), 1(sig,ee), 2(rd,ee), 3(bkg,mm), 4(sig,mm), 5(rd,mm), 
  const Int_t fl_mode_ll[Nchain] = {1,1,1,0,0,0};   // 1(e), 0(mu)
  const Double_t scale_event_sig = sigmc_amount*(used_nset/(Double_t)nset)/((Double_t)incl_ratio/100.0); // sigmc : N -> N/alpha
  const Double_t scale_event_bkg = used_nstream;                                                         //   gmc : N -> N/alpha

  const Int_t add[Nhist][Ntmp] ={
    {1,0,0,0,0,0}, // ee(qq)
    {1,1,0,1,1,1}, // ee(non-peak)
    {1,1,1,1,1,1}, // ee(peak)
    //{1,1,1,1,1,1,1,1}, // ee(sig)
    {1,1,0,0,0,0,0,0}, // ee(sig) // tmppp
    {0,0,0,0,0,0,0,0,1,0,0,0,0,0},     // mm
    {0,0,0,0,0,0,0,0,1,1,0,1,1,1},     // mm
    {0,0,0,0,0,0,0,0,1,1,1,1,1,1},     // mm
    //{0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1}, // mm
    {0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0}, // mm // tmpppp
    {1,0,0,0,0,0,0,0,1,0,0,0,0,0},      // ee+mm
    {1,1,0,1,1,1,0,0,1,1,0,1,1,1},      // ee+mm
    {1,1,1,1,1,1,0,0,1,1,1,1,1,1},      // ee+mm
    //{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},  // ee+mm
    {1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0},  // ee+mm // tmppppp
  };

  const Int_t Nbin_afb = 14;

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  Char_t** infile = new Char_t*[3];
  for( Int_t i=0; i<3; i++ ){
    infile[i] = new Char_t[1024];
    if     ( i==0 ) sTmp << indir[i] << "gMC_*_s0["    << stream << "]";  // bkg
    else if( i==1 ) sTmp << indir[i] << "sigMC_*_set[" << setname << "]"; // sig
    else if( i==2 ) sTmp << indir[i] << "RD_";                            // rd
    strcpy( infile[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Char_t* LRNB_cut = new Char_t[2048];
  strcpy( LRNB_cut, "1" ); // lrnb cut should be already applied.
  std::cout << "LRNB cut : " << LRNB_cut << std::endl;

  Char_t** add_cut = new Char_t*[Ntmp];
  for( Int_t i=0; i<Ntmp; i++ ){
    add_cut[i] = new Char_t[4096];
    sTmp << LRNB_cut;
    if     ( i%Ncategory==0 ){ // gmc(qq)
      sTmp << " && genbfl==0";
      
    }else if( i%Ncategory==1 ){ // gmc(bb-non-peak)
      ///*
      sTmp << " && genbfl!=0"
	   << " && !( (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && ( pi0org==1 || (pi0org<0.5&&rm_xs>999) || (pi0org<0.5&&rm_xs<999)) )" // charmonium events(right-pi0)
	   << " && !(!(lpgt==3 && lmgt==3 ) && (lpself==0 && lmself==0) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3) "  // double miss-id
           << " && !(!(lpgt==3 && lmgt==3 ) && ( (lpself==1&&lmself==0) || (lpself==0&&lmself==1) ) && ( lpgt==3 || lmgt==3 ) && rest_sw==0 && gm_bg1<3              && ( korg==lporg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg ) && lporg==lmorg)" // single miss-id (cc)
	   << " && !(!(lpgt==3 && lmgt==3 ) && ( (lpself==1&&lmself==0) || (lpself==0&&lmself==1) ) &&   lpgt!=3 && lmgt!=3   && rest_sw==0 && gm_bg1<3 && rm_xs>999 && ( korg==lporg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg ) && lporg==lmorg)"; // single miss-id with neutrino
      //*/
      //sTmp << " && genbfl!=0 && !(lpgt==3 && lmgt==3) && lpself==1 && lmself==1"; // gmc(semi-leptonic b) // tmpppppp
    }else if( i%Ncategory==2 ){ // gmc(bb-peak)
      ///*
      sTmp << " && genbfl!=0";
      sTmp << " && ("
	   << "    ( (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && ( pi0org==1 || (pi0org<0.5&&rm_xs>999) || (pi0org<0.5&&rm_xs<999)) )" // charmonium events(right-pi0)
	   << " || (!(lpgt==3 && lmgt==3 ) && (lpself==0 && lmself==0) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3) "  // double miss-id
	   << " || (!(lpgt==3 && lmgt==3 ) && ( (lpself==1&&lmself==0) || (lpself==0&&lmself==1) ) && ( lpgt==3 || lmgt==3 ) && rest_sw==0 && gm_bg1<3              && ( korg==lporg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg ) && lporg==lmorg)" // single miss-id (cc)
	   << " || (!(lpgt==3 && lmgt==3 ) && ( (lpself==1&&lmself==0) || (lpself==0&&lmself==1) ) &&   lpgt!=3 && lmgt!=3   && rest_sw==0 && gm_bg1<3 && rm_xs>999 && ( korg==lporg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg ) && lporg==lmorg)"; // single miss-id wit
      sTmp << " )";
      //*/
      //sTmp << " && genbfl!=0 && !(!(lpgt==3 && lmgt==3) && lpself==1 && lmself==1)"; // gmc(semi-leptonic b) // tmpppppp
    }else if( i%Ncategory==3 ){ // false with (f,f)
      sTmp << " && self!=1 && ";
      sTmp << makeCut_q2fl( 0,0 );
    }else if( i%Ncategory==4 ){ // false with (f,t)
      sTmp << " && self!=1 && ";
      sTmp << makeCut_q2fl( 0,1 );
    }else if( i%Ncategory==5 ){ // false with (t,f)
      sTmp << " && self!=1 && ";
      sTmp << makeCut_q2fl( 1,0 );
    }else if( i%Ncategory==6 ){ // false with (t,t)
      sTmp << " && self!=1 && ";
      sTmp << makeCut_q2fl( 1,1 );
    }else if( i%Ncategory==7 ){ // true
      sTmp << " && self==1 ";
    }
    
    if( flag_k4pi     ) sTmp << " && " << makeCut_5body_veto().c_str();
    if( flag_unflavor ) sTmp << " && " << makeCut_unflavor_veto().c_str();
    if( flag_ccpi0    ) sTmp << " && " << cut_ccpi0;
    strcpy( add_cut[i], (Char_t*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain   = new MChain*[Nchain]; // [gmc, sigmc, rd] x [ee, mm]
  TH1D**    tmphist = new TH1D*  [Ntmp];
  TH1D**    hist    = new TH1D*  [Nhist];
  TF1**     func    = new TF1*   [Nhist];
  Int_t*    N0bin   = new Int_t  [Nhist]; // # of zero-bins for chi2 calculation in RooFit
  TCanvas*  c1      = Canvas( "c1","c1", 3, 2 );

  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++
  std::cout << Form(" ************************ make chain-tree ( %s, %s ) *************************************",tname,axis) << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<infile %d > ", j );
    chain[j] = new MChain( infile[j%3], tname, branch_table(), 0, tail );
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
    chain[j]->MakeTree();
  }
  // +++++++ make tmphist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make tmp hist *************************************" << std::endl;
  for( Int_t j=0; j<Ntmp; j++ ){
    tmphist[j] = new TH1D( Form("tmphist%d",j), Form("%s",chain[fl_sb[j]]->GetChange()), xbin,offset+xmin,offset+xmax );
    chain[fl_sb[j]]->GetTree()->Project( Form("tmphist%d",j), axis, add_cut[j] );
    std::cout << Form( "<tmphist %d > ", j )
	      << "add_cut : " << add_cut[j] << std::endl;
    if( flag_scale ){
      tmphist[j]->Sumw2();
      if     ( fl_sb[j]==0 || fl_sb[j]==3) tmphist[j]->Scale( 1/scale_event_bkg );
      else if( fl_sb[j]==1 || fl_sb[j]==4) tmphist[j]->Scale( 1/scale_event_sig );
    }
  }
  
  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make hist *************************************" << std::endl;
  for( Int_t i=0; i<Nhist; i++ ){
    hist[i] = new TH1D( Form("hist%d",i),Form("hist%d",i), xbin,offset+xmin,offset+xmax );
    Deco( hist[i], 3, col_fil[i%Nplot], col_fil[i%Nplot] );
    for( Int_t j=0; j<Ntmp; j++ ){
      if( add[i][j] ) hist[i]->Add( tmphist[j] );
    }
  }

  // +++++++ display ++++++++++++++++++++++++++++++++++
  Double_t entry_sig[Nhist] = {0}; // # of events in signal box region
  Double_t total_evt_entry  = 0;
  Double_t total_evt_canvas = 0;
  Double_t total_evt_under  = 0;
  Double_t total_evt_over   = 0;
  Double_t total_evt_sig    = 0;
  for( Int_t i=0; i<Nhist; i++ ){
    std::cout << Form("<hist %d> ",i);
    Double_t entry_all    = hist[i]->GetEntries();
    Double_t entry_canvas = hist[i]->Integral();
    Double_t entry_under  = hist[i]->GetBinContent(0);
    Double_t entry_over   = hist[i]->GetBinContent(xbin+1);
    for( Int_t j=hist[i]->FindBin(5.27+0.000000001); j<=xbin; j++ ) entry_sig[i] += hist[i]->GetBinContent(j);
    std::cout << "Added-Files( ";
    for( Int_t j=0; j<Ntmp; j++ ) if( add[i][j] ) std::cout << j <<  ",";
    std::cout <<  ")["
	      << entry_all
	      << " events ( canvas : "
	      << entry_canvas
	      << " / under : "
	      << entry_under
	      << " / over  : "
	      << entry_over
	      << " / sig  : "
	      << entry_sig[i]
	      << "]"
	      << std::endl;
    total_evt_entry  += entry_all;
    total_evt_canvas += entry_canvas;
    total_evt_under  += entry_under;
    total_evt_over   += entry_over;
    total_evt_sig    += entry_sig[i];
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

  Double_t entry_sig_each[Ntmp] = {0}; // # of events in signal box region (tmphist)
  Double_t total_tmp_evt_canvas = 0;
  Double_t total_tmp_evt_sig    = 0;
  for( Int_t i=0; i<Ntmp; i++ ){
    std::cout << Form("<tmphist %d> ",i);
    for( Int_t j=tmphist[i]->FindBin(5.27+0.000000001); j<=xbin; j++ ) entry_sig_each[i] += tmphist[i]->GetBinContent(j);
    std::cout << entry_sig_each[i]      << " events(signal-box), "
	      << tmphist[i]->Integral() << " events(canvas)" << std::endl;
    total_tmp_evt_canvas += tmphist[i]->Integral();
    total_tmp_evt_sig    += entry_sig_each[i];
  }
  std::cout << "TOTAL "
	    << total_tmp_evt_sig    << " events(signal-box), "
	    << total_tmp_evt_canvas << " events(canvas)"
	    << std::endl << std::endl;

  // +++++++ make fit-func ++++++++++++++++++++++++++++++++++
  if( flag_fit ){
    for( Int_t i=0; i<Nhist; i++ ){
      func[i] = new TF1( Form("func%d",i), make_func(sel_fun),offset+xmin_fit,offset+xmax_fit,n_fitfunc_par(sel_fun) );
      if( sel_fun==15 ){
	Double_t mu    = 5.279;
	Double_t sigma = 0.003;
	Double_t area  = ( hist[i]->GetMaximum() - hist[i]->GetBinContent(1) )*sqrt( TMath::TwoPi() )*sigma;
	func[i]->SetParNames  ( "area","mean","sigma","norm",                         "Ebeam", "shape" );
	func[i]->SetParameters( area,   mu,    sigma,  2*hist[i]->GetBinContent(1),  endpoint,   -15   );
	func[i]->FixParameter( 1, 5.2794  );
	func[i]->FixParameter( 2, 0.00273 );
	if( flag_fix_end ) func[i]->FixParameter( 4, endpoint );
	else               func[i]->SetParLimits( 4, 5.280, 5.295 );
      }else if( sel_fun==151 ){
	Double_t mu    = 5.279;
	Double_t sigma = 0.003;
	Double_t area  = ( hist[i]->GetMaximum() - hist[i]->GetBinContent(1) )*sqrt( TMath::TwoPi() )*sigma;
	func[i]->SetParNames  ( "area","mean","sigma","norm",                          "Ebeam","shape", "new" );
	func[i]->SetParameters( area,   mu,    sigma,  2*hist[i]->GetBinContent(1),  endpoint,   -15,     0.46 );
	func[i]->FixParameter( 1, 5.2794 );
	func[i]->FixParameter( 2, 0.00273 );
	if( flag_fix_end ) func[i]->FixParameter( 4, endpoint );
	else               func[i]->SetParLimits( 4, 5.280, 5.295 );
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
  
  for(Int_t i=Nhist-1; i>=0; i-- ){
    c1->cd(i/Nplot+1);
    if     ( 0*Nplot <= i && i < 1*Nplot ) hist[i]->SetTitle( Form("ee"       ) );
    else if( 1*Nplot <= i && i < 2*Nplot ) hist[i]->SetTitle( Form("#mu#mu"   ) );
    else if( 2*Nplot <= i && i < 3*Nplot ) hist[i]->SetTitle( Form("ee+#mu#mu") );

    hist[i]->SetXTitle(xlabel);
    hist[i]->SetYTitle(waku[i/Nplot]->GetYaxis()->GetTitle());
    hist[i]->GetXaxis()->CenterTitle();
    hist[i]->GetYaxis()->CenterTitle();
    //waku[i]->Draw();

    if( i%Nplot==Nplot-1 ) hist[i]->Draw("hist");
    else                   hist[i]->Draw("hist same");
  }

  for(Int_t i=0; i<Nhist; i++ ){
    if( flag_fit && i%Nplot==Nplot-1 ){
      c1->cd(i/Nplot+4);
      std::cout << Form( "================================== FUNC%d =================================", i ) << std::endl;
      Double_t init_var[n_fitfunc_par(sel_fun)];
      for( Int_t m=0; m<n_fitfunc_par(sel_fun); m++ ) init_var[m] = func[i]->GetParameter(m);
      iterative_fit( hist[i], func[i], n_fitfunc_par(sel_fun), init_var, "PE0", manip_func );
    }
  }

  if( !flag_roof ){
    c1->Update();
    c1->Print( Form("pic/%s_fit_func%d_ratio%d_s0%s_set%s_c1.eps", axis, sel_fun, incl_ratio, stream, setname) );
    std::cout << "finish" << std::endl;
    if( fl_appRun ) app.Run();
    return 0;
  }
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++ RooFit ++++++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
  TCanvas* c2 = Canvas( "c2","c2", 3, 1 );
  // --- Observable ---
  const Int_t    Nroohist        = 3; // ee, mm , ee+mm
  const Int_t    rooad[Nroohist] = {3, 7, 11};
  Int_t          evt_cnt[Ntmp] = {0};
  

  RooRealVar** obs = new RooRealVar*[Nroohist];
  RooRealVar** wgt = new RooRealVar*[Nroohist];
  for( Int_t i=0; i<Nroohist; i++ ) obs[i] = new RooRealVar( axis, xlabel, offset+xmin_fit, offset+xmax_fit );
  for( Int_t i=0; i<Nroohist; i++ ) wgt[i] = new RooRealVar( Form("wgt%d",i), Form("wg%d",i), 0, 1.0 );
  Float_t    x_obs;
  // --- Data Set ---
  RooDataSet** data   = new RooDataSet*[Nroohist];
  TTree**      chain2 = new TTree*[Ntmp]; // apply add_cut[]
  for( Int_t i=0; i<Nroohist; i++ ) data  [i] = new RooDataSet( Form("data%d",i), Form("data%d",i), RooArgSet(*obs[i]) );
  for( Int_t j=0; j<Ntmp;     j++ ) chain2[j] = chain[fl_sb[j]]->GetTree()->CopyTree( add_cut[j] );

  for( Int_t j=0; j<Ntmp; j++ ){
    Int_t nevt = chain2[j]->GetTree()->GetEntries();
    std::cout << Form("<chain%d> nevt = %d", j, nevt) << std::endl;
    chain2[j]->GetTree()->SetBranchAddress( axis, &x_obs );
    for( Int_t k=0; k<nevt; k++ ){
      chain2[j]->GetTree()->GetEntry(k);
      if( x_obs>=offset+xmin && x_obs<=offset+xmax ){
	Double_t tmp_rnd = gRandom->Rndm();

	for( Int_t i=0; i<Nroohist; i++ ){
	  obs[i]->setVal( x_obs );
	  if( add[rooad[i]][j] ){
	    if( flag_scale && (fl_sb[j]==1 || fl_sb[j]==4) ){
	      if( tmp_rnd > 1-1/scale_event_sig ){
		data[i]->add( RooArgSet(*obs[i]) );
		if( i==Nroohist-1 ) evt_cnt[j]++;
	      }
	    }else if( flag_scale && (fl_sb[j]==0 || fl_sb[j]==3) ){
	      if( tmp_rnd > 1-1/scale_event_bkg ){
		data[i]->add( RooArgSet(*obs[i]) );
		if( i==Nroohist-1 ) evt_cnt[j]++;
	      }
	    }else{
	      data[i]->add( RooArgSet(*obs[i]) );
	      if( i==Nroohist-1 ) evt_cnt[j]++;
	    }
	  }
	}
      }
    }
  }
  
  const Double_t evt_max[Ntmp] = {197.8, 893.5, 45.67, 1.84, 1.41, 6.35, 18.25, 140.8,  //[qq, non-peak, peak, (f,f), (f,t), (t,f), (t,t), true]x[ee,mm] in Mbc>5.22, afb-modes
				  332.5, 1290,  45.67, 1.76, 1.41, 6.38, 17.06, 163.7}; // <- used only for initial values and display
  //const Double_t evt_max[Ntmp] = {197.8, 675.2, 264.0, 1.84, 1.41, 6.35, 18.25, 140.8,  // [qq, semi-leptonic b, except semi-leptonic b, (f,f), (f,t), (t,f), (t,t), true]x[ee,mm] in Mbc>5.22, afb-modes
  //332.5, 867.8, 468.3, 1.76, 1.41, 6.38, 17.06, 163.7}; // tmppppp
  
  for( Int_t i=0; i<Ntmp;     i++ ) std::cout << Form("< evt_cnt %d> ", i )    << evt_cnt[i]
					      << ", Norg = " << evt_max    [i] << " -> " << (Double_t)incl_ratio/100.0*evt_max    [i]
					      << std::endl;
  for( Int_t i=0; i<Nroohist; i++ ) std::cout << Form("< data%d> nevt = ", i ) << data[i]->numEntries()
					      << "( " << data[i]->sumEntries() << " )"
					      << std::endl; // sum of weight
  // ------------------------- Fit Fucntion ----------------------------

  const Int_t nparam = 5; // # of fit parameters

  // Sig-PDF(gauss)
  RooRealVar**  sig_mean  = new RooRealVar*[Nroohist];
  RooRealVar**  sig_sigma = new RooRealVar*[Nroohist];
  RooGaussian** gauss     = new RooGaussian*[Nroohist];
  // Fixed Parameters // 20121016 [gMC(6st) and RD, AFB-modes]
  sig_mean [0] = new RooRealVar ( "#mu",    "#mu",    5.27925); // MC
  sig_mean [1] = new RooRealVar ( "#mu",    "#mu",    5.27922); // MC
  sig_mean [2] = new RooRealVar ( "#mu",    "#mu",    5.27924); // MC
  sig_sigma[0] = new RooRealVar ( "#sigma", "#sigma", 0.00261); // MC
  sig_sigma[1] = new RooRealVar ( "#sigma", "#sigma", 0.00256); // MC
  sig_sigma[2] = new RooRealVar ( "#sigma", "#sigma", 0.00258); // MC
  //sig_mean [0] = new RooRealVar ( "#mu",    "#mu",    5.27936); // RD
  //sig_mean [1] = new RooRealVar ( "#mu",    "#mu",    5.27932); // RD
  //sig_mean [2] = new RooRealVar ( "#mu",    "#mu",    5.27934); // RD
  //sig_sigma[0] = new RooRealVar ( "#sigma", "#sigma", 0.00267); // RD
  //sig_sigma[1] = new RooRealVar ( "#sigma", "#sigma", 0.00258); // RD
  //sig_sigma[2] = new RooRealVar ( "#sigma", "#sigma", 0.00263); // RD
  for( Int_t i=0; i<Nroohist; i++ ) gauss[i] = new RooGaussian( "gauss", "gauss(x,mean,sigma)", *obs[i], *sig_mean[i], *sig_sigma[i] );

  // Bkg-PDF
  RooRealVar** arg_end   = new RooRealVar*[Nroohist];
  RooRealVar** arg_shape = new RooRealVar*[Nroohist];
  RooRealVar** arg_new   = new RooRealVar*[Nroohist];
  RooAbsPdf**  modargus  = new RooAbsPdf* [Nroohist];

  for( Int_t i=0; i<Nroohist; i++ ){
    arg_end  [i] = new RooRealVar( "end",   "argus endpoint",        endpoint, 5.285, 5.290 );
    if( flag_fix_end ) arg_end[i]->setConstant(kTRUE);
    if( sel_fun== 15 ){
      Double_t tmp_shape = func[rooad[i]]->GetParameter(5);
      arg_shape[i] = new RooRealVar( "shape", "argus shape parameter", tmp_shape, tmp_shape-20.0, tmp_shape+20.0 );
      modargus [i] = bindPdf( Form("modargus%d",i), func_roo_argus,      *obs[i],   *arg_end[i], *arg_shape[i] );
    }else if( sel_fun==151 ){
      Double_t tmp_shape = func[rooad[i]]->GetParameter(5);
      Double_t tmp_new   = func[rooad[i]]->GetParameter(6);
      arg_shape[i] = new RooRealVar( "shape", "argus shape parameter", tmp_shape, tmp_shape-20.0, tmp_shape+20.0 );
      arg_new  [i] = new RooRealVar( "new",   "argus sqrt parameter",  tmp_new,   tmp_new  -0.15, tmp_new  +0.15 );
      modargus[i] = bindPdf( Form("modargus%d",i), func_roo_modargus,   *obs[i], *arg_end[i], *arg_shape[i], *arg_new[i] );
    }
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // READ PDF (self cross-feed)
  TFile file_scf("pdf/Mbc_self_cf_tot_setA-U.root");
  if( file_scf.IsZombie() ) std::cerr << "[ABORT] can not find input-file : " << file_scf.GetName() << std::endl, abort();

  TH1D** hist_scf = new TH1D*[12];
  for( Int_t i=0; i<12; i++ ){
    hist_scf[i] = (TH1D*)file_scf.Get( Form("hist%d",i) ); 
    if( hist_scf[i] == NULL ) std::cerr << "[ABORT] can not find histgram (scf)" << i << std::endl, abort();
  }
  RooDataHist** scf_tt     = new RooDataHist*[Nroohist];
  RooDataHist** scf_tf     = new RooDataHist*[Nroohist];
  RooDataHist** scf_f      = new RooDataHist*[Nroohist];
  RooHistPdf**  pdf_scf_tt = new RooHistPdf* [Nroohist];
  RooHistPdf**  pdf_scf_tf = new RooHistPdf* [Nroohist];
  RooHistPdf**  pdf_scf_f  = new RooHistPdf* [Nroohist];
  
  for( Int_t i=0; i<3; i++ ){
    scf_tt[i] = new RooDataHist( Form("scf_tt%d", i), Form("scf_tt%d", i), *obs[i], hist_scf[4*i+1] );
    scf_tf[i] = new RooDataHist( Form("scf_tf%d", i), Form("scf_tf%d", i), *obs[i], hist_scf[4*i+2] );
    scf_f [i] = new RooDataHist( Form("scf_f%d",  i), Form("scf_f%d",  i), *obs[i], hist_scf[4*i+3] );
    pdf_scf_tt[i] = new RooHistPdf( Form("pdf_scf_tt%d",i), Form("pdf_scf_tt%d",i), *obs[i], *scf_tt[i] );
    pdf_scf_tf[i] = new RooHistPdf( Form("pdf_scf_tf%d",i), Form("pdf_scf_tf%d",i), *obs[i], *scf_tf[i] );
    pdf_scf_f [i] = new RooHistPdf( Form("pdf_scf_f%d", i), Form("pdf_scf_f%d", i), *obs[i], *scf_f [i] );
  }
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // SETTING of COEFFICIENT (signal and self cross-feed)
  RooRealVar**    nsig    = new RooRealVar*   [Nroohist];
  RooRealVar**    nbkg    = new RooRealVar*   [Nroohist];
  RooFormulaVar** nsig_tt = new RooFormulaVar*[Nroohist];
  RooFormulaVar** nsig_tf = new RooFormulaVar*[Nroohist];
  RooFormulaVar** nsig_f  = new RooFormulaVar*[Nroohist];
  RooRealVar**   r_tt     = new RooRealVar*   [Nroohist];
  RooRealVar**   r_tf     = new RooRealVar*   [Nroohist];
  RooRealVar**   r_f      = new RooRealVar*   [Nroohist];
  // AFB-modes
  r_tt[0] = new RooRealVar( "r_tt0", "r_tt0",     0.1299 );   r_tf[0] = new RooRealVar( "r_tf0", "r_tf0",    0.04642 );   r_f [0] = new RooRealVar( "r_f0",  "r_f0",     0.02237 );
  r_tt[1] = new RooRealVar( "r_tt1", "r_tt1",     0.1042 );   r_tf[1] = new RooRealVar( "r_tf1", "r_tf1",    0.03869 );   r_f [1] = new RooRealVar( "r_f1",  "r_f1",     0.01873 );
  r_tt[2] = new RooRealVar( "r_tt2", "r_tt2",     0.1161 );   r_tf[2] = new RooRealVar( "r_tf2", "r_tf2",    0.04226 );   r_f [2] = new RooRealVar( "r_f2",  "r_f2",     0.02041 );
  // All-Modes
  //r_tt[0] = new RooRealVar( "r_tt0", "r_tt0",     0.1174 );   r_tf[0] = new RooRealVar( "r_tf0", "r_tf0",    0.04127 );   r_f [0] = new RooRealVar( "r_f0",  "r_f0",     0.02007 );
  //r_tt[1] = new RooRealVar( "r_tt1", "r_tt1",    0.09531 );   r_tf[1] = new RooRealVar( "r_tf1", "r_tf1",    0.03481 );   r_f [1] = new RooRealVar( "r_f1",  "r_f1",     0.01704 );
  //r_tt[2] = new RooRealVar( "r_tt2", "r_tt2",     0.1055 );   r_tf[2] = new RooRealVar( "r_tf2", "r_tf2",     0.0378 );   r_f [2] = new RooRealVar( "r_f2",  "r_f2",     0.01844 );

  for( Int_t i=0; i<Nroohist; i++ ){
    Double_t tmp_nsig;
    if     ( i==0 ) tmp_nsig =  (evt_max[ 7]             *incl_ratio/100);
    else if( i==1 ) tmp_nsig =  (evt_max[15]             *incl_ratio/100);
    else if( i==2 ) tmp_nsig = ((evt_max[ 7]+evt_max[15])*incl_ratio/100);

    //const Double_t evt_peak[Nroohist]   = {47.8, 126.4, 174.2}; // ee, mm, ee+mm <- used only for initail values
    const Double_t evt_peak[Nroohist]   = {0.0, 0.0, 0.0}; // ee, mm, ee+mm // tmppppp

    nbkg[i] = new RooRealVar ( "N_{bkg}", "N_{bkg}", data[i]->sumEntries()-evt_peak[i], data[i]->sumEntries()-evt_peak[i]-200, data[i]->sumEntries()-evt_peak[i]+200 ); // tmppppppp
    //nbkg[i] = new RooRealVar ( "N_{bkg}", "N_{bkg}", data[i]->sumEntries()-tmp_nsig, data[i]->sumEntries()-tmp_nsig-200, data[i]->sumEntries()-tmp_nsig+200 ); // tmpppp
    //nbkg[i] = new RooRealVar ( "N_{bkg}", "N_{bkg}", data[i]->sumEntries()-tmp_nsig, data[i]->sumEntries()-evt_peak[i]-tmp_nsig-200, data[i]->sumEntries()-evt_peak[i]-tmp_nsig+200 );
    nsig[i] = new RooRealVar ( "N_{sig}", "N_{sig}", 10,                     -200,                       200 ); // tmppppp
    //nsig[i] = new RooRealVar ( "N_{sig}", "N_{sig}", tmp_nsig  ,                     tmp_nsig-200,                       tmp_nsig+200 );
    nsig_tt[i] = new RooFormulaVar( Form("nsig_tt%d",i),  "@0*@1", RooArgSet(*r_tt[i], *nsig[i]) );
    nsig_tf[i] = new RooFormulaVar( Form("nsig_tf%d",i),  "@0*@1", RooArgSet(*r_tf[i], *nsig[i]) );
    nsig_f [i] = new RooFormulaVar( Form("nsig_f%d", i),  "@0*@1", RooArgSet(*r_f [i], *nsig[i]) );
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // READ PDF (charmonium)
  TFile file_cc("pdf/Mbc_peak_cc_s00-5.root"); // gMC
  //TFile file_cc("pdf/Mbc_peak_cc_s00.root"  ); // CC-MC
  if( file_cc.IsZombie() ) std::cerr << "[ABORT] can not find input-file : " << file_cc.GetName() << std::endl, abort();
  
  TH1D** hist_cc = new TH1D*[(Nbin_afb+1)*3];
  for( Int_t i=0; i<(Nbin_afb+1)*3; i++ ){
    hist_cc[i] = (TH1D*)file_cc.Get( Form("hist%d",i) ); 
    if( hist_cc[i] == NULL ) std::cerr << "[ABORT] can not find histgram (charmonium)" << i << std::endl, abort();
  }
  RooDataHist** peak_cc     = new RooDataHist*[Nroohist];
  RooHistPdf**  pdf_peak_cc = new RooHistPdf* [Nroohist];
  
  for( Int_t i=0; i<3; i++ ){
    std::cout << "test" << i << std::endl;
    peak_cc[i]     = new RooDataHist( Form("peak_cc%d",     i), Form("peak_cc%d",     i), *obs[i], hist_cc[(Nbin_afb+1)*i+Nbin_afb] );
    std::cout << "test" << i << std::endl;
    pdf_peak_cc[i] = new RooHistPdf ( Form("pdf_peak_cc%d", i), Form("pdf_peak_cc%d", i), *obs[i], *peak_cc[i]   );
    std::cout << "test" << i << std::endl;
  }
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // SETTING of COEFFICIENT (charmonium)
  RooRealVar** r_peak_cc = new RooRealVar*[Nroohist];
  const Double_t scale_cc[3] = {1.0, 1.0, 1.0 }; // for MC
  std::cerr << "[WARNING] cc-normalization is obsolete !! should be updated" << std::endl;
  // CC-MC, AFB-modes
  r_peak_cc[0] = new RooRealVar( "r_peak_cc0", "r_peak_cc0", scale_cc[0]*14.25 );
  r_peak_cc[1] = new RooRealVar( "r_peak_cc1", "r_peak_cc1", scale_cc[1]*10.82 );
  r_peak_cc[2] = new RooRealVar( "r_peak_cc2", "r_peak_cc2", scale_cc[2]*25.07 );
  // CC-MC, All-modes
  //const Double_t scale_cc[3] = {0.815, 0.738, 0.774};
  //r_peak_cc[0] = new RooRealVar( "r_peak_cc0", "r_peak_cc0", scale_cc[0]*17.16 );
  //r_peak_cc[1] = new RooRealVar( "r_peak_cc1", "r_peak_cc1", scale_cc[1]*12.76 );
  //r_peak_cc[2] = new RooRealVar( "r_peak_cc2", "r_peak_cc2", scale_cc[2]*29.92 );
  // gMC, AFB-modes  
  //const Double_t scale_cc[3] = {0.882, 0.801, 0.838};
  //r_peak_cc[0] = new RooRealVar( "r_peak_cc", "r_peak_cc", scale_cc[0]*43.67 );
  //r_peak_cc[1] = new RooRealVar( "r_peak_cc", "r_peak_cc", scale_cc[1]*25.33 );
  //r_peak_cc[2] = new RooRealVar( "r_peak_cc", "r_peak_cc", scale_cc[2]*69.00 );
  // gMC, All-modes  
  //const Double_t scale_cc[3] = {0.864, 0.787, 0.823};
  //r_peak_cc[0] = new RooRealVar( "r_peak_cc", "r_peak_cc", scale_cc[0]*50.50 );
  //r_peak_cc[1] = new RooRealVar( "r_peak_cc", "r_peak_cc", scale_cc[1]*27.50 );
  //r_peak_cc[2] = new RooRealVar( "r_peak_cc", "r_peak_cc", scale_cc[2]*78.00 );
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // READ PDF (double miss-PID)
  TFile file_double("pdf/Mbc_peak_double_s00.root");
  //TFile file_double("pdf_data/Mbc_peak_double_s00.root");
  if( file_double.IsZombie() ) std::cerr << "[ABORT] can not find input-file : " << file_double.GetName() << std::endl, abort();

  TH1D** hist_double = new TH1D*[(Nbin_afb+1)*3];
  for( Int_t i=0; i<(Nbin_afb+1)*3; i++ ){
    hist_double[i] = (TH1D*)file_double.Get( Form("hist%d",i) ); 
    if( hist_double[i] == NULL ) std::cerr << "[ABORT] can not find histgram (double)" << i << std::endl, abort();
  }
  RooDataHist** peak_double     = new RooDataHist*[Nroohist];
  RooHistPdf**  pdf_peak_double = new RooHistPdf* [Nroohist];
  
  for( Int_t i=0; i<3; i++ ){
    peak_double[i]     = new RooDataHist( Form("peak_double%d",     i), Form("peak_double%d",     i), *obs[i], hist_double[(Nbin_afb+1)*i+Nbin_afb] );
    pdf_peak_double[i] = new RooHistPdf ( Form("pdf_peak_double%d", i), Form("pdf_peak_double%d", i), *obs[i], *peak_double[i]   );
  }
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // SETTING of COEFFICIENT (double miss-PID)
  RooRealVar** r_peak_double = new RooRealVar*[Nroohist];
  // AFB-modes(MC)
  r_peak_double[0] = new RooRealVar( "r_peak_double0", "r_peak_double0",  1.098/0.273 ); //  ~4.0
  r_peak_double[1] = new RooRealVar( "r_peak_double1", "r_peak_double1", 25.65 /0.273 ); // ~94.0
  r_peak_double[2] = new RooRealVar( "r_peak_double2", "r_peak_double2", 26.75 /0.273 ); // ~98.0
  // All-modes(MC)
  //r_peak_double[0] = new RooRealVar( "r_peak_double0", "r_peak_double0",  1.268/0.273 );
  //r_peak_double[1] = new RooRealVar( "r_peak_double1", "r_peak_double1", 30.45 /0.273 );
  //r_peak_double[2] = new RooRealVar( "r_peak_double2", "r_peak_double2", 31.72 /0.273 );
  // AFB-modes(RD)
  //r_peak_double[0] = new RooRealVar( "r_peak_double0", "r_peak_double0",  0.2146/0.273 );
  //r_peak_double[1] = new RooRealVar( "r_peak_double1", "r_peak_double1", 12.34  /0.273 );
  //r_peak_double[2] = new RooRealVar( "r_peak_double2", "r_peak_double2", 12.55  /0.273 );
  // All-modes(RD)
  //r_peak_double[0] = new RooRealVar( "r_peak_double0", "r_peak_double0",  0.2473/0.273 );
  //r_peak_double[1] = new RooRealVar( "r_peak_double1", "r_peak_double1", 14.63  /0.273 );
  //r_peak_double[2] = new RooRealVar( "r_peak_double2", "r_peak_double2", 14.87  /0.273 );

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // READ PDF (swapped miss-PID)
  TFile file_swap("pdf/Mbc_peak_swap_s00.root");
  //TFile file_swap("pdf_data/Mbc_peak_swap_s00.root");
  if( file_swap.IsZombie() ) std::cerr << "[ABORT] can not find input-file : " << file_swap.GetName() << std::endl, abort();

  TH1D** hist_swap = new TH1D*[(Nbin_afb+1)*3];
  for( Int_t i=0; i<(Nbin_afb+1)*3; i++ ){
    hist_swap[i] = (TH1D*)file_swap.Get( Form("hist%d",i) ); 
    if( hist_swap[i] == NULL ) std::cerr << "[ABORT] can not find histgram (swap)" << i << std::endl, abort();
  }
  RooDataHist** peak_swap     = new RooDataHist*[Nroohist];
  RooHistPdf**  pdf_peak_swap = new RooHistPdf* [Nroohist];
  
  for( Int_t i=0; i<3; i++ ){
    peak_swap[i]     = new RooDataHist( Form("peak_swap%d",     i), Form("peak_swap%d",     i), *obs[i], hist_swap[(Nbin_afb+1)*i+Nbin_afb] );
    pdf_peak_swap[i] = new RooHistPdf ( Form("pdf_peak_swap%d", i), Form("pdf_peak_swap%d", i), *obs[i], *peak_swap[i]   );
  }
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // SETTING of COEFFICIENT (swapped miss-PID)
  RooRealVar** r_peak_swap = new RooRealVar*[Nroohist];
  // AFB-modes(MC)
  r_peak_swap[0] = new RooRealVar( "r_peak_swap0", "r_peak_swap0", 0.2796 );
  r_peak_swap[1] = new RooRealVar( "r_peak_swap1", "r_peak_swap1", 8.463  );
  r_peak_swap[2] = new RooRealVar( "r_peak_swap2", "r_peak_swap2", 8.742  );
  // All-modes(MC)
  //r_peak_swap[0] = new RooRealVar( "r_peak_swap0", "r_peak_swap0", 0.317   );
  //r_peak_swap[1] = new RooRealVar( "r_peak_swap1", "r_peak_swap1", 8.9.337 );
  //r_peak_swap[2] = new RooRealVar( "r_peak_swap2", "r_peak_swap2", 9.654   );
  // AFB-modes(RD)
  //r_peak_swap[0] = new RooRealVar( "r_peak_swap0", "r_peak_swap0", 0.09232 );
  //r_peak_swap[1] = new RooRealVar( "r_peak_swap1", "r_peak_swap1", 3.93    );
  //r_peak_swap[2] = new RooRealVar( "r_peak_swap2", "r_peak_swap2", 4.022   );
  // All-modes(RD)
  //r_peak_swap[0] = new RooRealVar( "r_peak_swap0", "r_peak_swap0", 0.09828 );
  //r_peak_swap[1] = new RooRealVar( "r_peak_swap1", "r_peak_swap1", 4.103   );
  //r_peak_swap[2] = new RooRealVar( "r_peak_swap2", "r_peak_swap2", 4.202   );

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Total-PDF
  RooAddPdf**    pdf        = new RooAddPdf*[Nroohist];
  RooFitResult** fit_result = new RooFitResult*[Nroohist];
  for( Int_t i=0; i<Nroohist; i++ ){
    pdf[i]  = new RooAddPdf( Form("pdf%d",i), Form("pdf%d",i), RooArgList(*gauss[i], *modargus[i]                                               ),                                                          RooArgList(*nsig[i], *nbkg[i]                                                                                         ) ); // signal
    //pdf[i]  = new RooAddPdf( Form("pdf%d",i), Form("pdf%d",i), RooArgList(*gauss[i], *modargus[i], *pdf_scf_tt[i], *pdf_scf_tf[i], *pdf_scf_f[i]),                                                          RooArgList(*nsig[i], *nbkg[i], *nsig_tt[i], *nsig_tf[i], *nsig_f[i]                                                   ) ); // signal + scf
    //pdf[i]  = new RooAddPdf( Form("pdf%d",i), Form("pdf%d",i), RooArgList(*gauss[i], *modargus[i], *pdf_peak_cc[i], *pdf_peak_double[i], *pdf_peak_swap[i]),                                                RooArgList(*nsig[i], *nbkg[i], *r_peak_cc[i], *r_peak_double[i], *r_peak_swap[i]                                      ) ); // signal + peak
    //pdf[i]  = new RooAddPdf( Form("pdf%d",i), Form("pdf%d",i), RooArgList(*gauss[i], *modargus[i], *pdf_peak_cc[i]), RooArgList(*nsig[i], *nbkg[i], *r_peak_cc[i] ) ); // signal + peak(cc)
    //pdf[i]  = new RooAddPdf( Form("pdf%d",i), Form("pdf%d",i), RooArgList(*gauss[i], *modargus[i], *pdf_peak_double[i]), RooArgList(*nsig[i], *nbkg[i], *r_peak_double[i] ) ); // signal + peak(double)
    //pdf[i]  = new RooAddPdf( Form("pdf%d",i), Form("pdf%d",i), RooArgList(*gauss[i], *modargus[i], *pdf_scf_tt[i], *pdf_scf_tf[i], *pdf_scf_f[i], *pdf_peak_cc[i], *pdf_peak_double[i], *pdf_peak_swap[i]), RooArgList(*nsig[i], *nbkg[i], *nsig_tt[i], *nsig_tf[i], *nsig_f[i], *r_peak_cc[i], *r_peak_double[i], *r_peak_swap[i]) ); // signal + scf + peak

    if( sel_fun==15 || sel_fun==151 ){
      fit_result[i] = pdf[i]->fitTo( *data[i], Extended(), Save(true) );
      /*
      arg_end  [i]->setConstant(kTRUE );
      arg_shape[i]->setConstant(kFALSE);
      nsig     [i]->setConstant(kFALSE);
      nbkg     [i]->setConstant(kFALSE);
      fit_result[i] = pdf[i]->fitTo( *data[i], Extended(), Save(true) );
      arg_end  [i]->setConstant(kFALSE);
      arg_shape[i]->setConstant(kTRUE );
      nsig     [i]->setConstant(kTRUE );
      nbkg     [i]->setConstant(kTRUE );
      fit_result[i] = pdf[i]->fitTo( *data[i], Extended(), Save(true) );
      */
    }
  }
  /*
  TCanvas* c3 = Canvas( "c3","c3", 3, 3 );
  c3->Draw();
  std::cout << "[TESTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT]" << std::endl;
  RooAbsReal* nll1 = pdf[0]->createNLL( *data[0] );
  RooPlot* frame1 = arg_end[0]->frame( Bins(50),Range(5.285, 5.295), Title("profile") );
  nll1->plotOn( frame1,ShiftToZero() );
  frame1->SetMinimum(0);
  frame1->SetMaximum(20);
  c3->cd(1);
  frame1->Draw();

  RooPlot* frame2 = arg_shape[0]->frame( Bins(50),Range(-30, 0), Title("profile") );
  nll1->plotOn( frame2,ShiftToZero() );
  frame2->SetMinimum(0);
  frame2->SetMaximum(3);
  c3->cd(2);
  frame2->Draw();

  RooPlot* frame3 = nsig[0]->frame( Bins(50),Range(-50, 50), Title("profile") );
  nll1->plotOn( frame3,ShiftToZero() );
  frame3->SetMinimum(0);
  frame3->SetMaximum(3);
  c3->cd(3);
  frame3->Draw();

  ////////////////////////////////
  RooAbsReal* nll2 = pdf[1]->createNLL( *data[1] );
  RooPlot* frame4 = arg_end[1]->frame( Bins(50),Range(5.285, 5.295), Title("profile") );
  nll2->plotOn( frame4,ShiftToZero() );
  frame4->SetMinimum(0);
  frame4->SetMaximum(20);
  c3->cd(4);
  frame4->Draw();

  RooPlot* frame5 = arg_shape[1]->frame( Bins(50),Range(-30, 0), Title("profile") );
  nll2->plotOn( frame5,ShiftToZero() );
  frame5->SetMinimum(0);
  frame5->SetMaximum(3);
  c3->cd(5);
  frame5->Draw();

  RooPlot* frame6 = nsig[1]->frame( Bins(50),Range(-50, 50), Title("profile") );
  nll2->plotOn( frame6,ShiftToZero() );
  frame6->SetMinimum(0);
  frame6->SetMaximum(3);
  c3->cd(6);
  frame6->Draw();

  ////////////////////////////////
  RooAbsReal* nll3 = pdf[2]->createNLL( *data[2] );
  RooPlot* frame7 = arg_end[2]->frame( Bins(50),Range(5.285, 5.295), Title("profile") );
  nll3->plotOn( frame7,ShiftToZero() );
  frame7->SetMinimum(0);
  frame7->SetMaximum(20);
  c3->cd(7);
  frame7->Draw();

  RooPlot* frame8 = arg_shape[2]->frame( Bins(50),Range(-30, 0), Title("profile") );
  nll3->plotOn( frame8,ShiftToZero() );
  frame8->SetMinimum(0);
  frame8->SetMaximum(3);
  c3->cd(8);
  frame8->Draw();

  RooPlot* frame9 = nsig[2]->frame( Bins(50),Range(-50, 50), Title("profile") );
  nll3->plotOn( frame9,ShiftToZero() );
  frame9->SetMinimum(0);
  frame9->SetMaximum(3);
  c3->cd(9);
  frame9->Draw();

  c3->Print( Form("pic/%s_fit_func%d_ratio%d_s0%s_set%s_c3.eps",  axis, sel_fun, incl_ratio, stream, setname) );
  */
  // ------------------------- Draw ----------------------------

  RooPlot**    frame = new RooPlot*  [Nroohist];
  TPaveText**  box   = new TPaveText*[Nroohist];

  for( Int_t i=0; i<Nroohist; i++ ){
    frame[i] = obs[i]->frame();
    frame[i]->GetXaxis()->CenterTitle();
    frame[i]->GetYaxis()->CenterTitle();
    frame[i]->SetTitleOffset( 1.00,"x" );
    frame[i]->SetTitleOffset( 1.30,"y" );
    if     ( i==0 ) frame[i]->SetTitle( "ee"        );
    else if( i==1 ) frame[i]->SetTitle( "#mu#mu   " );
    else if( i==2 ) frame[i]->SetTitle( "ee+#mu#mu" );
    data[i]->plotOn ( frame[i], Binning(xbin), LineWidth(1) );
    if( sel_fun==50 || sel_fun==51 ){
      modargus[i]->plotOn( frame[i], ProjWData(*data[i]), LineWidth(1) );      
      box[i] = new TPaveText( 0.05, 0.87, 0.35, 0.93,"BRNDC" );
      box[i]->AddText( Form("#chi^{2}/ndf = %f", frame[i]->chiSquare(nparam+N0bin[i])) );
      modargus[i] ->paramOn( frame[i], Format("NELU", AutoPrecision(2)), Layout(0.40, 0.99, 0.99) );
    }else if( sel_fun==15 || sel_fun==151 ){
      pdf[i] ->plotOn ( frame[i], ProjWData(*data[i]), LineWidth(1) );
      box[i] = new TPaveText( 0.05, 0.87, 0.35, 0.93,"BRNDC" );
      box[i]->AddText( Form("#chi^{2}/ndf = %f", frame[i]->chiSquare(nparam+N0bin[i])) );
      pdf[i] ->plotOn ( frame[i], ProjWData(*data[i]), Components(*gauss          [i]), LineStyle(7), LineColor(2), LineWidth(1) ); // sig (gauss)
      pdf[i] ->plotOn ( frame[i], ProjWData(*data[i]), Components(*modargus       [i]), LineStyle(7), LineColor(4), LineWidth(1) ); // bkg (arguss)
      //pdf[i] ->plotOn ( frame[i], ProjWData(*data[i]), Components(*pdf_scf_tt     [i]), LineStyle(1), LineColor(5), LineWidth(1) ); // self cros-feed (t,t)
      //pdf[i] ->plotOn ( frame[i], ProjWData(*data[i]), Components(*pdf_scf_tf     [i]), LineStyle(1), LineColor(5), LineWidth(1) ); // self cros-feed (t,f)
      //pdf[i] ->plotOn ( frame[i], ProjWData(*data[i]), Components(*pdf_scf_f      [i]), LineStyle(1), LineColor(5), LineWidth(1) ); // self cros-feed (f,?)
      //pdf[i] ->plotOn ( frame[i], ProjWData(*data[i]), Components(*pdf_peak_cc    [i]), LineStyle(1), LineColor(5), LineWidth(1) ); // peak(cc)
      //pdf[i] ->plotOn ( frame[i], ProjWData(*data[i]), Components(*pdf_peak_double[i]), LineStyle(1), LineColor(5), LineWidth(1) ); // peak(double)
      //pdf[i] ->plotOn ( frame[i], ProjWData(*data[i]), Components(*pdf_peak_swap  [i]), LineStyle(1), LineColor(5), LineWidth(1) ); // peak(swap)
      pdf[i] ->paramOn( frame[i], Format("NELU", AutoPrecision(2)), Layout(0.40, 0.99, 0.99) );
    }

    frame[i]->addObject( box[i] );
    c2->cd(i+1);
    frame[i]->Draw();
  }

  c1->Update();
  c2->Update();

  if( flag_save ){
    c1->Print( Form("pic/%s_fit_func%d_ratio%d_s0%s_set%s_c1.eps",  axis, sel_fun, incl_ratio, stream, setname) );
    c2->Print( Form("pic/%s_fit_func%d_ratio%d_s0%s_set%s_c2.eps",  axis, sel_fun, incl_ratio, stream, setname) );
  }

  for( Int_t i=0; i<Nroohist; i++ ){
    
    
    std::cout << std::setw( 5) << std::right << incl_ratio/100.0                              << " "
	      << std::setw(10) << std::right << nsig[i]->getVal()                             << " "
	      << std::setw(10) << std::right << nsig[i]->getPropagatedError( *fit_result[i] ) << " "
	      << "  HOGE" << i 
	      << "  FUNC" << sel_fun
      //<< std::setw(10) << std::right << fit_result[i]->status()  << " "
	      << std::endl;
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
  delete[] modargus;
  delete[] nsig;
  delete[] nbkg;
  delete[] pdf;
  delete[] fit_result;
  delete[] frame;
  delete[] box;
  delete   c2;

  return 0;
}

