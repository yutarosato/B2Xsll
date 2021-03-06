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
#include <RooCFunction1Binding.h>
#include <RooCFunction2Binding.h>
#include <RooCFunction3Binding.h>
#include <RooCFunction4Binding.h>
#include <RooTFnBinding.h>

const Bool_t   flag_fix_end = !true; // true(fix endpoint parameter)
const Double_t endpoint     = 5.289;

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
  using namespace RooFit;
  //using namespace sig_gmc_rd_cut2;
  using namespace sig_gmc_rd;
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==6 || argc==7) ) std::cerr << "wrong input" << std::endl
					<< " Usage : ./draws_1d_Mbc_fit (char*)stream (char*)setname (double)nstream (double)nset (int)sel_fun [(int)fl_appRun]" << std::endl
					<< " [sel_fun] : 50(argus), 51(modified-argus), 15(gauss+argus), 151(gauss+modified-argus)"
					<< std::endl, abort();

  Char_t*  stream       = argv[1];
  Char_t*  setname      = argv[2];
  Double_t used_nstream = atof(argv[3]);
  Double_t used_nset    = atof(argv[4]);
  Int_t    sel_fun      = atoi(argv[5]);
  Int_t    fl_appRun    = 1;
  if( argc==7 ) fl_appRun = atoi( argv[6] );

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t Ncategory     = 6; //[gmc(qq,non-peak,peak),sigmc(false except(t,t), false with (t,t), true)]
  const Int_t Nplot         = 4; //[bkg(qq,non-peak,peak),sigmc]
  const Int_t Nchain        = Ncategory*2; // x[ee,mm      ]
  const Int_t Nhist         = Nplot    *3; // x[ee,mm,ee+mm]
  const Int_t fl_sb[Nchain] = {0,0,0,1,1,1,
			       0,0,0,1,1,1}; // 0(bkg), 1(sig), 2(rd)
  const Int_t fl_mode_ll[Nchain] = {1,1,1,1,1,1,
				    0,0,0,0,0,0}; // 1(e), 0(mu)
  const Int_t nfile[Nchain] = {0};
  const Double_t scale_event_sig = sigmc_amount*(used_nset/(Double_t)nset); // sigmc : N -> N/alpha
  const Double_t scale_event_bkg = used_nstream;                            //   gmc : N -> N/alpha

  const Int_t add[Nhist][Nchain] ={
    {1,0,0,0,0,0}, // ee(qq)
    {1,1,0,1,0,0}, // ee(non-peak)
    {1,1,1,1,0,0}, // ee(peak)
    {1,1,1,1,1,1}, // ee(sig)
    {0,0,0,0,0,0,1,0,0,0,0,0}, // mm
    {0,0,0,0,0,0,1,1,0,1,0,0}, // mm
    {0,0,0,0,0,0,1,1,1,1,0,0}, // mm
    {0,0,0,0,0,0,1,1,1,1,1,1}, // mm
    {1,0,0,0,0,0,1,0,0,0,0,0}, // ee+mm
    {1,1,0,1,0,0,1,1,0,1,0,0}, // ee+mm
    {1,1,1,1,0,0,1,1,1,1,0,0}, // ee+mm
    {1,1,1,1,1,1,1,1,1,1,1,1}, // ee+mm
  };

  const Bool_t   flag_scale    = true;
  const Bool_t   flag_k4pi     = true; // 1(veto  K4pi    modes)
  const Bool_t   flag_unflavor = true; // 1(veto unflavor modes)
  const Bool_t   flag_ccpi0    = true; // 1(veto ccpi0 peak    )

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  Char_t** infile = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    infile[i] = new Char_t[1024];
    if     ( fl_sb[i]==0 ) sTmp << indir[fl_sb[i]] << "gMC_*_s0[" << stream << "]";     // bkg
    else if( fl_sb[i]==1 ) sTmp << indir[fl_sb[i]] << "sigMC_*_set[" << setname << "]"; // sig
    else if( fl_sb[i]==2 ) sTmp << indir[fl_sb[i]] << "RD_";                            // rd
    strcpy( infile[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Char_t* LRNB_cut = new Char_t[2048];
  strcpy( LRNB_cut, "1" ); // lrnb cut should be already applied.
  //strcpy( LRNB_cut, (char*)makeCut_LRNB_1d    ("lr1_%s",0, 0.85, 0.94 ).c_str() ); // 1d LR
  std::cout << "LRNB cut : " << LRNB_cut << std::endl;

  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    add_cut[i] = new Char_t[4096];
    sTmp << LRNB_cut;
    if     ( i%Ncategory==0 ){ // gmc(qq)
      sTmp << " && genbfl==0";
    }else if( i%Ncategory==1 ){ // gmc(bb-non-peak)
      sTmp << " && genbfl!=0"
	   << " && !( (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && ( pi0org==1 || (pi0org<0.5&&rm_xs>999) || (pi0org<0.5&&rm_xs<999)) )" // charmonium events(right-pi0)
	   << " && !(!(lpgt==3 && lmgt==3 ) && (lpself==0 && lmself==0) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3) "  // double miss-id
           << " && !(!(lpgt==3 && lmgt==3 ) && ( (lpself==1&&lmself==0) || (lpself==0&&lmself==1) ) && ( lpgt==3 || lmgt==3 ) && rest_sw==0 && gm_bg1<3              && ( korg==lporg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg ) && lporg==lmorg)" // single miss-id (cc)
	   << " && !(!(lpgt==3 && lmgt==3 ) && ( (lpself==1&&lmself==0) || (lpself==0&&lmself==1) ) &&   lpgt!=3 && lmgt!=3   && rest_sw==0 && gm_bg1<3 && rm_xs>999 && ( korg==lporg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg ) && lporg==lmorg)"; // single miss-id with neutrino
    }else if( i%Ncategory==2 ){ // gmc(bb-peak)
      sTmp << " && genbfl!=0";
      sTmp << " && ("
	   << "    ( (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && ( pi0org==1 || (pi0org<0.5&&rm_xs>999) || (pi0org<0.5&&rm_xs<999)) )" // charmonium events(right-pi0)
	   << " || (!(lpgt==3 && lmgt==3 ) && (lpself==0 && lmself==0) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3) "  // double miss-id
	   << " || (!(lpgt==3 && lmgt==3 ) && ( (lpself==1&&lmself==0) || (lpself==0&&lmself==1) ) && ( lpgt==3 || lmgt==3 ) && rest_sw==0 && gm_bg1<3              && ( korg==lporg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg ) && lporg==lmorg)" // single miss-id (cc)
	   << " || (!(lpgt==3 && lmgt==3 ) && ( (lpself==1&&lmself==0) || (lpself==0&&lmself==1) ) &&   lpgt!=3 && lmgt!=3   && rest_sw==0 && gm_bg1<3 && rm_xs>999 && ( korg==lporg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg ) && lporg==lmorg)"; // single miss-id wit
      sTmp << " )";
    }else if( i%Ncategory==3 ){ // false except (t,t)
      sTmp << LRNB_cut << " && self!=1 && ";
      sTmp << makeCut_q2fl( 1,1,1 );
    }else if( i%Ncategory==4 ){ // false with (t,t)
      sTmp << LRNB_cut << " && self!=1 && ";
      sTmp << makeCut_q2fl( 1,1 );
    }else if( i%Ncategory==5 ){ // true
      sTmp << LRNB_cut << " && self==1 ";
    }
    
    if( flag_k4pi     ) sTmp << " && " << makeCut_5body_veto().c_str();
    if( flag_unflavor ) sTmp << " && " << makeCut_unflavor_veto().c_str();
    if( flag_ccpi0    ) sTmp << " && " << cut_ccpi0;
    strcpy( add_cut[i], (Char_t*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
    
  using namespace Mbc_bkg;
  //using namespace Mbc_bkg_wide;

  const Bool_t flag_save = true; // outfile.eps and outfile.root
  const Bool_t flag_fit  = true;
  const Bool_t flag_roof = false; // flase(skip roofit)
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain   = new MChain*[Nchain];
  TH1D**    tmphist = new TH1D*  [Nchain];
  TH1D**    hist    = new TH1D*  [Nhist];
  TF1**     func    = new TF1*   [Nhist];
  Int_t*    N0bin   = new Int_t  [Nhist]; // # of zero-bins for chi2 calculation in RooFit
  TCanvas*  c1      = Canvas( "c1","c1", 3, 2 );

  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++
  std::cout << Form(" ************************ make chain-tree ( %s, %s ) *************************************",tname,axis) << std::endl;
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
    chain[j]->MakeTree();
    std::cout << "add_cut : " << add_cut[j] << std::endl;
  }
  // +++++++ make tmphist ++++++++++++++++++++++++++++++++++
  for( Int_t j=0; j<Nchain; j++ ){
    tmphist[j] = new TH1D( Form("tmphist%d",j), Form("%s",chain[j]->GetChange()), xbin,offset+xmin,offset+xmax );
    chain[j]->GetTree()->Project( Form("tmphist%d",j), axis, add_cut[j] );
    if( flag_scale ){
      /*
      tmphist[j]->Sumw2();
      if     ( fl_sb[j]==1 ) tmphist[j]->Scale( 1/scale_event_sig );
      else if( fl_sb[j]==0 ) tmphist[j]->Scale( 1/scale_event_bkg );
      */
      if( fl_sb[j]==1 ) tmphist[j]->Scale( scale_event_bkg/scale_event_sig );
      tmphist[j]->Sumw2();
      tmphist[j]->Scale( 1/scale_event_bkg );
    }
  }
  
  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make hist *************************************" << std::endl;
  for( Int_t i=0; i<Nhist; i++ ){
    hist[i] = new TH1D( Form("hist%d",i),Form("hist%d",i), xbin,offset+xmin,offset+xmax );
    Deco( hist[i], 3, col_fil[i%Nplot], col_fil[i%Nplot] );
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
    std::cout << "Added-Files( ";
    for( Int_t j=0; j<Nchain; j++ ) if( add[i][j] ) std::cout << j <<  ",";
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
  }
  std::cout << std::endl;

  Double_t entry_sig_each[Nchain] = {0}; // # of events in signal box region (tmphist)
  for( Int_t i=0; i<Nchain; i++ ){
    std::cout << Form("<tmphist %d> ",i);
    for( Int_t j=tmphist[i]->FindBin(5.27+0.000000001); j<=xbin; j++ ) entry_sig_each[i] += tmphist[i]->GetBinContent(j);
    std::cout << entry_sig_each[i] << " events" << std::endl;
  }

  // +++++++ make fit-func ++++++++++++++++++++++++++++++++++
  if( flag_fit ){
    for( Int_t i=0; i<Nhist; i++ ){
      func[i] = new TF1( Form("func%d",i), make_func(sel_fun),offset+xmin_fit,offset+xmax_fit,n_fitfunc_par(sel_fun) );
      if( sel_fun==15 ){
	Double_t mu    = 5.279;
	Double_t sigma = 0.003;
	Double_t area  = ( hist[i]->GetMaximum() - hist[i]->GetBinContent(1) )*sqrt( TMath::TwoPi() )*sigma;
	func[i]->SetParNames  ( "area","mean","sigma","norm",                          "Ebeam","shape" );
	func[i]->SetParameters( area,   mu,    sigma,  2*hist[i]->GetBinContent(1),  endpoint,   -15 );
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
    else                           hist[i]->Draw("hist same");
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
    c1->Print( Form("pic/%s_peak_fit_func%d_s0%s_set%s_c1.eps", axis, sel_fun, stream, setname) );
    /*
    TCanvas* c3 = Canvas( "c3","c3", 2, 2 );
    c3->cd(1); tmphist[ 3]->Draw();
    c3->cd(2); tmphist[ 4]->Draw();
    c3->cd(3); tmphist[ 9]->Draw();
    c3->cd(4); tmphist[10]->Draw();
    c3->Update();
    c3->Print( Form("pic/%s_peak_fit_func%d_s0%s_set%s_c3.eps", axis, sel_fun, stream, setname) );
    */
  
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
  const Int_t Nroohist = 3;
  const Int_t rooad[Nroohist] = {3, 7, 11};
  RooRealVar** obs = new RooRealVar*[Nroohist];
  RooRealVar** wgt = new RooRealVar*[Nroohist];
  for( Int_t i=0; i<Nroohist; i++ ) obs[i] = new RooRealVar( axis, xlabel, offset+xmin_fit, offset+xmax_fit );
  for( Int_t i=0; i<Nroohist; i++ ) wgt[i] = new RooRealVar( Form("wgt%d",i), Form("wg%d",i), 0, 1.0 );
  Float_t    x_obs;
  Float_t    cc_m, cc_morg, rm_l;
  // --- Data Set ---
  RooDataSet** data   = new RooDataSet*[Nroohist];
  TTree**      chain2 = new TTree*[Nchain]; // apply add_cut[]
  for( Int_t i=0; i<Nroohist;  i++ ) data[i]   = new RooDataSet( Form("data%d",i), Form("data%d",i), RooArgSet(*obs[i],*wgt[i]), WeightVar(*wgt[i]) );
  for( Int_t j=0; j<Nchain; j++ ) chain2[j] = chain[j]->GetTree()->CopyTree( add_cut[j] );

  for( Int_t j=0; j<Nchain; j++ ){
    Int_t nevt = chain2[j]->GetTree()->GetEntries();
    std::cout << Form("<chain%d> nevt = %d", j, nevt) << std::endl;
    chain2[j]->GetTree()->SetBranchAddress( axis, &x_obs );
    for( Int_t k=0; k<nevt; k++ ){
      chain2[j]->GetTree()->GetEntry(k);
      if( x_obs>=offset+xmin && x_obs<=offset+xmax ){
	for( Int_t i=0; i<Nroohist; i++ ){
	  obs[i]->setVal( x_obs );
	  if( add[rooad[i]][j] ){
	    if     ( flag_scale && fl_sb[j]==1 ) data[i]->add( RooArgSet(*obs[i]), 1/scale_event_sig );
	    else if( flag_scale && fl_sb[j]==0 ) data[i]->add( RooArgSet(*obs[i]), 1/scale_event_bkg );
	    else                                 data[i]->add( RooArgSet(*obs[i]), 1.0 );
	  }
	}
      }
    }
  }

  for( Int_t i=0; i<Nroohist; i++ ) std::cout << Form("< data%d> nevt = ", i ) << data[i]->numEntries()
					      << "( " << data[i]->sumEntries() << " )"
					      << std::endl; // sum of weight
  // ------------------------- Fit Fucntion ----------------------------

  const Int_t nparam = 5; // # of fit parameters

  // Sig-PDF(gauss)
  RooRealVar**  sig_mean  = new RooRealVar*[Nroohist];
  RooRealVar**  sig_sigma = new RooRealVar*[Nroohist];
  RooGaussian** gauss     = new RooGaussian*[Nroohist];
  // Fixed Parameters
  sig_mean [0] = new RooRealVar ( "#mu",    "#mu",    5.27927);
  sig_mean [1] = new RooRealVar ( "#mu",    "#mu",    5.27925);
  sig_mean [2] = new RooRealVar ( "#mu",    "#mu",    5.27926);
  sig_sigma[0] = new RooRealVar ( "#sigma", "#sigma", 0.00259);
  sig_sigma[1] = new RooRealVar ( "#sigma", "#sigma", 0.00254);
  sig_sigma[2] = new RooRealVar ( "#sigma", "#sigma", 0.00256);
  for( Int_t i=0; i<Nroohist; i++ ){
    //sig_mean [i] = new RooRealVar ( "#mu",    "#mu",    5.280, 5.275, 5.285 );
    //sig_sigma[i] = new RooRealVar ( "#sigma", "#sigma", 0.001, 0.000, 0.005 );
    gauss    [i] = new RooGaussian( "gauss", "gauss(x,mean,sigma)", *obs[i], *sig_mean[i], *sig_sigma[i] );
  }
  // fix

  // Bkg-PDF
  RooRealVar** arg_end   = new RooRealVar*[Nroohist];
  RooRealVar** arg_shape = new RooRealVar*[Nroohist];
  RooRealVar** arg_new   = new RooRealVar*[Nroohist];
  RooAbsPdf**  modargus  = new RooAbsPdf* [Nroohist];

  for( Int_t i=0; i<Nroohist; i++ ){
    arg_end  [i] = new RooRealVar( "end",   "argus endpoint",        endpoint, 5.285, 5.290 );
    if( flag_fix_end ) arg_end[i]->setConstant(kTRUE);
    if( sel_fun== 15 ){
      arg_shape[0] = new RooRealVar( "shape", "argus shape parameter", -13.14, -25.0, -8.0 );
      arg_shape[1] = new RooRealVar( "shape", "argus shape parameter", -12.22, -17.0, -8.0 );
      arg_shape[2] = new RooRealVar( "shape", "argus shape parameter", -13.60, -20.0, -8.0 );
      modargus [i] = bindPdf( Form("modargus%d",i), func_roo_argus,      *obs[i], *arg_end[i], *arg_shape[i] );
    }else if( sel_fun==151 ){
      arg_shape[0] = new RooRealVar( "shape", "argus shape parameter", -23.75, -30.0, -10.0 );
      arg_shape[1] = new RooRealVar( "shape", "argus shape parameter",  -8.36, -20.0,  -4.0 );
      arg_shape[2] = new RooRealVar( "shape", "argus shape parameter", -14.35, -20.0, -10.0 );
      arg_new  [0] = new RooRealVar( "new",   "argus sqrt parameter",    0.50,  0.35,  0.60 );
      arg_new  [1] = new RooRealVar( "new",   "argus sqrt parameter",    0.39,  0.30,  0.60 );
      arg_new  [2] = new RooRealVar( "new",   "argus sqrt parameter",    0.46,  0.35,  0.60 );
      modargus[i] = bindPdf( Form("modargus%d",i), func_roo_modargus,   *obs[i], *arg_end[i], *arg_shape[i], *arg_new[i] );
    }
  }

  // Total-PDF
  RooAddPdf**    pdf        = new RooAddPdf*[Nroohist];
  RooRealVar**   nsig       = new RooRealVar*[Nroohist];
  RooRealVar**   nbkg       = new RooRealVar*[Nroohist];
  RooFitResult** fit_result = new RooFitResult*[Nroohist];
  for( Int_t i=0; i<Nroohist; i++ ){
    nsig[i] = new RooRealVar ( "N_{sig}", "N_{sig}", 0.10*data[i]->sumEntries(),  0.05*data[i]->sumEntries(), 0.15*data[i]->sumEntries() );
    nbkg[i] = new RooRealVar ( "N_{bkg}", "N_{bkg}", 0.90*data[i]->sumEntries(),  0.85*data[i]->sumEntries(), 0.95*data[i]->sumEntries() );
    pdf[i]  = new RooAddPdf( Form("pdf%d",i), Form("pdf%d",i), RooArgList(*gauss[i],*modargus[i]), RooArgList(*nsig[i], *nbkg[i]) );
    if( sel_fun==15 || sel_fun==151 ) fit_result[i] = pdf[i]     ->fitTo( *data[i], Extended(), SumW2Error(kFALSE) );
  } 

  // ------------------------- Draw ----------------------------

  RooPlot**    frame       = new RooPlot*  [Nroohist];
  TPaveText**  box         = new TPaveText*[Nroohist];

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
      pdf[i] ->plotOn ( frame[i], ProjWData(*data[i]), Components(*gauss[i]), LineStyle(7), LineColor(2), LineWidth(1) );
      pdf[i] ->plotOn ( frame[i], ProjWData(*data[i]), Components(*modargus[i]), LineStyle(7), LineColor(4), LineWidth(1) );
      pdf[i] ->paramOn( frame[i], Format("NELU", AutoPrecision(2)), Layout(0.40, 0.99, 0.99) );
    }

    frame[i]->addObject( box[i] );
    c2->cd(i+1);
    frame[i]->Draw();
  }

  c1->Update();
  c2->Update();

  if( flag_save ){
    c1->Print( Form("pic/%s_fit_func%d_s0%s_set%s_c1.eps",  axis, sel_fun, stream, setname) );
    c2->Print( Form("pic/%s_fit_func%d_s0%s_set%s_c2.eps",  axis, sel_fun, stream, setname) );
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

