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


using namespace RooFit;
using namespace gmc;
//using namespace Mbc_bkg;
using namespace Mbc_bkg_wide;

const Bool_t flag_save = true; // outfile.eps and outfile.root
const Bool_t flag_fit  = true;
const Bool_t flag_roof = false; // flase(skip roofit) // not supported yet
const Bool_t flag_norm = true;

  const Bool_t   flag_k4pi     = true; // 1(veto  K4pi    modes)
  const Bool_t   flag_unflavor = true; // 1(veto unflavor modes)

const Bool_t   flag_fix_end = true; // true(fix endpoint parameter)
const Double_t endpoint     = 5.289;

void manip_func( TF1* func ){
  if( func->GetNpar()==3 || func->GetNpar()==4 ){ // (modified-)argus
    if( flag_fix_end ) func->FixParameter( 1, endpoint );
    else               func->SetParLimits( 1, 5.285, 5.290 );
  }else if( func->GetNpar()==6 || func->GetNpar()==7 ){ // gaussian + (modified-)argus
    func->FixParameter( 1, 5.2794  );
    func->FixParameter( 2, 0.00273 );
    if( flag_fix_end ) func->FixParameter( 4, endpoint );
    else               func->SetParLimits( 4, 5.285, 5.290 );
  }
}

Int_t main( Int_t argc, Char_t** argv ){
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==4 || argc==5) ) std::cerr << "wrong input" << std::endl
					<< " Usage : ./draws_1d_Mbc_fit (char*)stream (double)nstream (int)sel_fun [(int)fl_appRun]" << std::endl
					<< " [sel_fun] : 50(argus), 51(modified-argus), 15(gauss+argus), 151(gauss+modified-argus)"
					<< std::endl, abort();

  Char_t*  stream       = argv[1];
  Double_t used_nstream = atof(argv[2]);
  Int_t    sel_fun      = atoi(argv[3]);
  Int_t    fl_appRun    = 1;
  if( argc==5 ) fl_appRun = atoi( argv[4] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  sTmp << indir << "/gMC_*_s0[" << stream << "]";
  Char_t tmp_infile[255];
  strcpy( tmp_infile, (Char_t*)sTmp.str().c_str() );
  sTmp.str("");
  sTmp.clear();
  const Char_t* infile = (const Char_t*)tmp_infile;

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t Nctgry        = 2*6;      // [two coslp regions] x [six q^2 regions]
  const Int_t Nchain        = 2*Nctgry; // [ee, mm]
  const Int_t Nplot         = 2;        // [two q^2 regions]
  const Int_t Nhist         = 3*Nplot;  // [ee, mm, ee+mm]
  const Double_t scale_event_bkg = used_nstream; // gmc : N -> N/alpha
  const Int_t fl_mode_ll[Nchain] = {1,1,1,1,1,1,1,1,1,1,1,1,
				    0,0,0,0,0,0,0,0,0,0,0,0}; // 1(e), 0(mu)
  const Int_t add[Nhist][Nchain] ={
    /* [six q^2 regions]
    {1,0,0,0,0,0,1,0,0,0,0,0},   // ee
    {0,1,0,0,0,0,0,1,0,0,0,0},   // ee
    {0,0,1,0,0,0,0,0,1,0,0,0},   // ee
    {0,0,0,1,0,0,0,0,0,1,0,0},   // ee
    {0,0,0,0,1,0,0,0,0,0,1,0},   // ee
    {0,0,0,0,0,1,0,0,0,0,0,1},   // ee
    {0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0},   // mm
    {0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0},   // mm
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0},   // mm
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0},   // mm
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0},   // mm
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1},   // mm
    {1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0},   // ee+mm
    {0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0},   // ee+mm
    {0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0},   // ee+mm
    {0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0},   // ee+mm
    {0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0},   // ee+mm
    {0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1},   // ee+mm
    */
    /* [two q^2 regions]
    {1,1,1,0,0,0,1,1,1,0,0,0},   // ee
    {0,0,0,1,1,1,0,0,0,1,1,1},   // ee
    {0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,1,1,1,0,0,0},   // mm
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,1,1,1},   // mm
    {1,1,1,0,0,0,1,1,1,0,0,0,1,1,1,0,0,0,1,1,1,0,0,0},   // ee+mm
    {0,0,0,1,1,1,0,0,0,1,1,1,0,0,0,1,1,1,0,0,0,1,1,1},   // ee+mm
    */

    ///* [coslp]
    {1,1,1,1,1,1,0,0,0,0,0,0},   // ee
    {0,0,0,0,0,0,1,1,1,1,1,1},   // ee
    {0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0},   // mm
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1},   // mm
    {1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0},   // ee+mm
    {0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1},   // ee+mm
    //*/
  };
  
  Char_t* LRNB_cut = new Char_t[2048];
  strcpy( LRNB_cut, "1" ); // lrnb cut should be already applied.
  std::cout << "LRNB cut : " << LRNB_cut << std::endl;

  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    add_cut[i] = new Char_t[4096];
    ///*
    sTmp << "    !( genbfl!=0 &&  (lpgt==3 && lmgt==3 ) && lpself==1 && lmself==1 && rest_sw==0 && gm_bg1<10 ) "; // charmonium events
    sTmp << " && !( genbfl!=0 && !(lpgt==3 && lmgt==3 ) && lpself==0 && lmself==0                                                         && rest_sw==0 && gm_bg1<10 && lporg==lmorg )"; // double miss-id
    sTmp << " && !( genbfl!=0 && !(lpgt==3 && lmgt==3 ) && ( (lpself==1&&lmself==0) || (lpself==0&&lmself==1) ) && ( lpgt==3 || lmgt==3 ) && rest_sw==0 && gm_bg1<10              && ( abs(korg)==lporg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg ) && lporg==lmorg )"; // single miss-id (cc)
    sTmp << " && !( genbfl!=0 && !(lpgt==3 && lmgt==3 ) && ( (lpself==1&&lmself==0) || (lpself==0&&lmself==1) ) &&   lpgt!=3 && lmgt!=3   && rest_sw==0 && gm_bg1<10 && rm_xs>999 && ( abs(korg)==lporg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg ) && lporg==lmorg )"; // single miss-id with neutrino
    //*/
    /*
    sTmp << "    !( genbfl!=0 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && ( pi0org==1 || (pi0org<0.5&&rm_xs>999) || (pi0org<0.5&&rm_xs<999)) )" // charmonium events(right-pi0)
	 << " && !( genbfl!=0 && !(lpgt==3 && lmgt==3 ) && (lpself==0 && lmself==0) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3) "  // double miss-id
	 << " && !( genbfl!=0 && !(lpgt==3 && lmgt==3 ) && ( (lpself==1&&lmself==0) || (lpself==0&&lmself==1) ) && ( lpgt==3 || lmgt==3 ) && rest_sw==0 && gm_bg1<3              && ( korg==lporg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg ) && lporg==lmorg)" // single miss-id (cc)
	 << " && !( genbfl!=0 && !(lpgt==3 && lmgt==3 ) && ( (lpself==1&&lmself==0) || (lpself==0&&lmself==1) ) &&   lpgt!=3 && lmgt!=3   && rest_sw==0 && gm_bg1<3 && rm_xs>999 && ( korg==lporg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg ) && lporg==lmorg)"; // single miss-id with neutrino
    */

    if     ( i%Nctgry ==  0 ) sTmp << " && " << makeCut_q2( fl_mode_ll[i], 1, -1 );
    else if( i%Nctgry ==  1 ) sTmp << " && " << makeCut_q2( fl_mode_ll[i], 2, -1 );
    else if( i%Nctgry ==  2 ) sTmp << " && " << makeCut_q2( fl_mode_ll[i], 3, -1 );
    else if( i%Nctgry ==  3 ) sTmp << " && " << makeCut_q2( fl_mode_ll[i], 5, -1 );
    else if( i%Nctgry ==  4 ) sTmp << " && " << makeCut_q2( fl_mode_ll[i], 7, -1 );
    else if( i%Nctgry ==  5 ) sTmp << " && " << makeCut_q2( fl_mode_ll[i], 8, -1 );
    else if( i%Nctgry ==  6 ) sTmp << " && " << makeCut_q2( fl_mode_ll[i], 1,  1 );
    else if( i%Nctgry ==  7 ) sTmp << " && " << makeCut_q2( fl_mode_ll[i], 2,  1 );
    else if( i%Nctgry ==  8 ) sTmp << " && " << makeCut_q2( fl_mode_ll[i], 3,  1 );
    else if( i%Nctgry ==  9 ) sTmp << " && " << makeCut_q2( fl_mode_ll[i], 5,  1 );
    else if( i%Nctgry == 10 ) sTmp << " && " << makeCut_q2( fl_mode_ll[i], 7,  1 );
    else if( i%Nctgry == 11 ) sTmp << " && " << makeCut_q2( fl_mode_ll[i], 8,  1 );

    /*
    if( i==0 ) sTmp << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl
		    << " [1,-1] " << makeCut_q2( fl_mode_ll[i], 1, -1 ) << std::endl
		    << " [2,-1] " << makeCut_q2( fl_mode_ll[i], 2, -1 ) << std::endl
		    << " [3,-1] " << makeCut_q2( fl_mode_ll[i], 3, -1 ) << std::endl
		    << " [5,-1] " << makeCut_q2( fl_mode_ll[i], 5, -1 ) << std::endl
		    << " [7,-1] " << makeCut_q2( fl_mode_ll[i], 7, -1 ) << std::endl
		    << " [8,-1] " << makeCut_q2( fl_mode_ll[i], 8, -1 ) << std::endl
		    << " [1, 1] " << makeCut_q2( fl_mode_ll[i], 1,  1 ) << std::endl
		    << " [2, 1] " << makeCut_q2( fl_mode_ll[i], 2,  1 ) << std::endl
		    << " [3, 1] " << makeCut_q2( fl_mode_ll[i], 3,  1 ) << std::endl
		    << " [5, 1] " << makeCut_q2( fl_mode_ll[i], 5,  1 ) << std::endl
		    << " [7, 1] " << makeCut_q2( fl_mode_ll[i], 7,  1 ) << std::endl
		    << " [8, 1] " << makeCut_q2( fl_mode_ll[i], 8,  1 ) << std::endl
		    << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    */
    if( flag_k4pi     ) sTmp << " && " << makeCut_5body_veto().c_str();
    if( flag_unflavor ) sTmp << " && " << makeCut_unflavor_veto().c_str();
    strcpy( add_cut[i], (Char_t*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }   
 
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain   = new MChain*[2]; // [mm, ee]
  TH1D**    tmphist = new TH1D*  [Nchain];
  TH1D**    hist    = new TH1D*  [Nhist];
  TF1**     func    = new TF1*   [Nhist];
  Int_t*    N0bin   = new Int_t  [Nhist]; // # of zero-bins for chi2 calculation in RooFit
  TCanvas*  c1      = Canvas( "c1","c1",3, flag_fit ? Nplot+1 : 1 );

  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++
  std::cout << Form(" ************************ make chain-tree ( %s, %s ) *************************************",tname,axis) << std::endl;
  for( Int_t m=0; m<2; m++ ){
    std::cout << Form( "<infile %d > ", m );
    chain[m] = new MChain( infile, tname, branch_table(), 0, tail );
    nominal_cut_selection( chain[m], m )( chain[m]->GetCut(), tname );
  }

  // ++++++++++++++++++++++++
  // cut change
  for( Int_t m=0; m<2; m++ ){
    chain[m]->GetCut()->Set( "Mbc", 1,  5.20, 0.0, 5.30 );
  }
  //
  // ++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ changed cut *************************************" << std::endl;
  for( Int_t m=0; m<2; m++ ){
    chain[m]->GetCut()->Display(0);
    chain[m]->MakeTree();
  }
  // +++++++ make tmphist ++++++++++++++++++++++++++++++++++
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<tmphist %d > ", j ) << std::endl;
    std::cout << "add_cut : " << add_cut[j] << std::endl;
    tmphist[j] = new TH1D( Form("tmphist%d",j), Form("%s",chain[fl_mode_ll[j]]->GetChange()), xbin,offset+xmin,offset+xmax );
    chain[fl_mode_ll[j]]->GetTree()->Project( Form("tmphist%d",j), axis, add_cut[j] );
    tmphist[j]->Sumw2();
    tmphist[j]->Scale( 1/scale_event_bkg );
  }
  
  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make hist *************************************" << std::endl;
  for( Int_t i=0; i<Nhist; i++ ){
    hist[i] = new TH1D( Form("hist%d",i),Form("hist%d",i), xbin,offset+xmin,offset+xmax );
    Deco( hist[i], 1, i%Nplot+1, i%Nplot+1 );
    for( Int_t j=0; j<Nchain; j++ ){
      if( add[i][j] ) hist[i]->Add( tmphist[j] );
    }
    if( flag_norm ) hist[i]->Scale( 1/hist[i]->Integral() );
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
      if( sel_fun==50 ){
	func[i]->SetParNames  ("norm","Ebeam","shape");
	func[i]->SetParameters( 100*hist[i]->GetBinContent(1), endpoint, -10);
	if( flag_fix_end ) func[i]->FixParameter( 1, endpoint );
	else               func[i]->SetParLimits( 1, 5.285, 5.290 );
      }else if( sel_fun==51 ){
	func[i]->SetParNames  ("norm","Ebeam","shape","new");
	func[i]->SetParameters( 100*hist[i]->GetBinContent(1), endpoint, -10, 0.5);
	if( flag_fix_end ) func[i]->FixParameter( 1, endpoint );
	else               func[i]->SetParLimits( 1, 5.285, 5.290 );
      }else if( sel_fun==15 ){
	Double_t mu    = 5.279;
	Double_t sigma = 0.003;
	Double_t area  = 0.1;
	func[i]->SetParNames  ( "area","mean","sigma","norm",                          "Ebeam","shape" );
	func[i]->SetParameters( area,   mu,    sigma,  1.5*hist[i]->GetBinContent(1),  endpoint,   -20 );
	func[i]->FixParameter( 1, 5.2794  );
	func[i]->FixParameter( 2, 0.00273 );
	if( flag_fix_end ) func[i]->FixParameter( 4, endpoint );
	else               func[i]->SetParLimits( 4, 5.285, 5.290 );
      }else if( sel_fun==151 ){
	Double_t mu    = 5.279;
	Double_t sigma = 0.003;
	Double_t area  = 0.1;
	func[i]->SetParNames  ( "area","mean","sigma","norm",                          "Ebeam","shape", "new" );
	func[i]->SetParameters( area,   mu,    sigma,  1.5*hist[i]->GetBinContent(1),  endpoint,   -10,     0.5 );
	func[i]->FixParameter( 1, 5.2794 );
	func[i]->FixParameter( 2, 0.00273 );
	if( flag_fix_end ) func[i]->FixParameter( 4, endpoint );
	else               func[i]->SetParLimits( 4, 5.285, 5.290 );
      }else func_set_parameters(sel_fun, func[i], hist[i], xbin, offset+xmin, offset+xmax);
      func[i]->SetLineColor(2);
    }
  }
  

  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();
  TH2D** waku = new TH2D*[3];
  waku[0] = Waku( Nhist, hist, Form("%s (ee)",       xlabel) );
  waku[1] = Waku( Nhist, hist, Form("%s (#mu#mu)",   xlabel) );
  waku[2] = Waku( Nhist, hist, Form("%s (ee+#mu#mu)",xlabel) );
    
  for(Int_t i=0; i<Nhist; i++ ){
    c1->cd(i/Nplot+1);
    if     ( i/Nplot==0 ) hist[i]->SetTitle( Form("ee"       ) );
    else if( i/Nplot==1 ) hist[i]->SetTitle( Form("#mu#mu"   ) );
    else if( i/Nplot==2 ) hist[i]->SetTitle( Form("ee+#mu#mu") );
    hist[i]->SetXTitle(xlabel);
    hist[i]->SetYTitle(waku[i/Nplot]->GetYaxis()->GetTitle());
    hist[i]->GetXaxis()->CenterTitle();
    hist[i]->GetYaxis()->CenterTitle();
    if( i%Nplot==0 ) waku[i/Nplot]->Draw();
    hist[i]->DrawCopy( "same" );
    if( flag_fit ){
      c1->cd( 4+ i/Nplot + 3*(i%Nplot) );
      Deco( hist[i], 1, 1, 1 );
      std::cout << Form( "================================== FUNC%d =================================", i ) << std::endl;
      Double_t init_var[n_fitfunc_par(sel_fun)];
      for( Int_t m=0; m<n_fitfunc_par(sel_fun); m++ ) init_var[m] = func[i]->GetParameter(m);
      iterative_fit( hist[i], func[i], n_fitfunc_par(sel_fun), init_var, "PE0", manip_func );
    }
  }
  
  if( !flag_roof ){
    c1->Update();
    c1->Print( Form("pic/%s_shape_func%d_s0%s_c1.eps", axis, sel_fun, stream) );
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
  TCanvas* c2 = Canvas( "c2","c2",Nhist, 2 );
  // --- Observable ---
  RooRealVar** obs = new RooRealVar*[Nhist];
  RooRealVar** wgt = new RooRealVar*[Nhist];
  for( Int_t i=0; i<Nhist; i++ ) obs[i] = new RooRealVar( axis, xlabel, offset+xmin_fit, offset+xmax_fit );
  for( Int_t i=0; i<Nhist; i++ ) wgt[i] = new RooRealVar( Form("wgt%d",i), Form("wg%d",i), 0.00, 100.0 );
  Float_t    x_obs;
  Float_t    cc_m, cc_morg, rm_l;
  // --- Data Set ---
  RooDataSet** data   = new RooDataSet*[Nhist];
  TTree**      chain2 = new TTree*[Nhist]; // apply add_cut[]
  for( Int_t i=0; i<Nhist;  i++ ) data[i]   = new RooDataSet( Form("data%d",i), Form("data%d",i), RooArgSet(*obs[i],*wgt[i]), WeightVar(*wgt[i]) );
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
	  if( add[i][j] ){
	    data[i]->add( RooArgSet(*obs[i]), 1.0 );
	  }
	}
      }
    }
  }

  for( Int_t i=0; i<Nhist; i++ ) std::cout << Form("< data%d> nevt = ", i ) << data[i]->numEntries()
					   << "( " << data[i]->sumEntries() << " )"
					   << std::endl; // sum of weight
  // ------------------------- Fit Fucntion ----------------------------

  const Int_t nparam = 5; // # of fit parameters

  // Sig-PDF(gauss)
  RooRealVar**  sig_mean  = new RooRealVar*[Nhist];
  RooRealVar**  sig_sigma = new RooRealVar*[Nhist];
  RooGaussian** gauss     = new RooGaussian*[Nhist];
  // Fixed Parameters
  sig_mean [0] = new RooRealVar ( "#mu",    "#mu",    5.27927);
  sig_mean [1] = new RooRealVar ( "#mu",    "#mu",    5.27925);
  sig_mean [2] = new RooRealVar ( "#mu",    "#mu",    5.27926);
  sig_sigma[0] = new RooRealVar ( "#sigma", "#sigma", 0.00259);
  sig_sigma[1] = new RooRealVar ( "#sigma", "#sigma", 0.00254);
  sig_sigma[2] = new RooRealVar ( "#sigma", "#sigma", 0.00256);
  for( Int_t i=0; i<Nhist; i++ ){
    //sig_mean [i] = new RooRealVar ( "#mu",    "#mu",    5.280, 5.275, 5.285 );
    //sig_sigma[i] = new RooRealVar ( "#sigma", "#sigma", 0.001, 0.000, 0.005 );
    gauss    [i] = new RooGaussian( "gauss", "gauss(x,mean,sigma)", *obs[i], *sig_mean[i], *sig_sigma[i] );
  }

  // Bkg-PDF
  RooRealVar** arg_end   = new RooRealVar*[Nhist];
  RooRealVar** arg_shape = new RooRealVar*[Nhist];
  RooRealVar** arg_new   = new RooRealVar*[Nhist];
  RooAbsPdf**  modargus  = new RooAbsPdf* [Nhist];

  for( Int_t i=0; i<Nhist; i++ ){
    arg_end  [i] = new RooRealVar( "end",   "argus endpoint",        endpoint, 5.280, 5.295 );
    if( sel_fun== 50 || sel_fun== 15 ){
      arg_shape[0] = new RooRealVar( "shape", "argus shape parameter", -20.92,  -25.0,  -5.0  );
      arg_shape[1] = new RooRealVar( "shape", "argus shape parameter", -20.96,  -25.0,  -5.0  );
      arg_shape[2] = new RooRealVar( "shape", "argus shape parameter", -12.88,  -18.0,  -8.0  );
    }else if( sel_fun== 51 || sel_fun==151 ){
      arg_shape[0] = new RooRealVar( "shape", "argus shape parameter", -13.10,  -18.0,  -8.0  );
      arg_shape[1] = new RooRealVar( "shape", "argus shape parameter", -20.18,  -25.0,  -8.0  );
      arg_shape[2] = new RooRealVar( "shape", "argus shape parameter", -9.53,  -15.0,  -5.0  );
      arg_new  [0] = new RooRealVar( "new",   "argus sqrt parameter",  0.473, 0.35, 0.60 );
      arg_new  [1] = new RooRealVar( "new",   "argus sqrt parameter",  0.450, 0.35, 0.60 );
      arg_new  [2] = new RooRealVar( "new",   "argus sqrt parameter",  0.450, 0.35, 0.60 );
    }else if( sel_fun==81 || sel_fun==82 || sel_fun==181 || sel_fun==182 ){
      arg_shape[i] = new RooRealVar( "shape", "threshold shape parameter", -0.48,  -0.55,  -0.40 );
      arg_new  [i] = new RooRealVar( "new",   "argus sqrt parameter", 0.35, 0.3, 0.55 );
    }
    if( flag_fix_end ) arg_end[i]->setConstant(kTRUE);
    if     ( sel_fun==50 || sel_fun== 15 ) modargus[i] = bindPdf( Form("modargus%d",i), func_roo_argus,      *obs[i], *arg_end[i], *arg_shape[i]);
    else if( sel_fun==51 || sel_fun==151 ) modargus[i] = bindPdf( Form("modargus%d",i), func_roo_modargus,   *obs[i], *arg_end[i], *arg_shape[i], *arg_new[i] );
    else if( sel_fun==81 || sel_fun==181 ) modargus[i] = bindPdf( Form("modargus%d",i), func_roo_threshold1, *obs[i], *arg_end[i], *arg_shape[i]);
    else if( sel_fun==82 || sel_fun==182 ) modargus[i] = bindPdf( Form("modargus%d",i), func_roo_threshold2, *obs[i], *arg_end[i], *arg_shape[i], *arg_new[i] );
  }

  // Total-PDF
  RooAddPdf**    pdf        = new RooAddPdf*[Nhist];
  RooRealVar**   nsig       = new RooRealVar*[Nhist];
  RooRealVar**   nbkg       = new RooRealVar*[Nhist];
  RooFitResult** fit_result = new RooFitResult*[Nhist];
  for( Int_t i=0; i<Nhist; i++ ){
    nsig[i] = new RooRealVar ( "N_{sig}", "N_{sig}", 0.01*data[i]->sumEntries(), -0.02*data[i]->sumEntries(), 0.05*data[i]->sumEntries() );
    nbkg[i] = new RooRealVar ( "N_{bkg}", "N_{bkg}", 0.99*data[i]->sumEntries(),  0.95*data[i]->sumEntries(), 1.02*data[i]->sumEntries() );
    pdf[i]  = new RooAddPdf( Form("pdf%d",i), Form("pdf%d",i), RooArgList(*gauss[i],*modargus[i]), RooArgList(*nsig[i], *nbkg[i]) );
    if     ( sel_fun==50 || sel_fun== 51 || sel_fun== 81 || sel_fun== 82 ) fit_result[i] = modargus[i]->fitTo( *data[i]             );
    else if( sel_fun==15 || sel_fun==151 || sel_fun==181 || sel_fun==182 ) fit_result[i] = pdf[i]     ->fitTo( *data[i], Extended() );
  } 

  // ------------------------- Draw ----------------------------

  RooPlot**    frame       = new RooPlot*  [Nhist];
  TPaveText**  box         = new TPaveText*[Nhist];
  TH1D**       pullhist    = new TH1D*[Nhist];
  TH2D**       pullwaku    = new TH2D*[Nhist];

  for( Int_t i=0; i<Nhist; i++ ){
    frame[i] = obs[i]->frame();
    frame[i]->GetXaxis()->CenterTitle();
    frame[i]->GetYaxis()->CenterTitle();
    frame[i]->SetTitleOffset( 1.00,"x" );
    frame[i]->SetTitleOffset( 1.30,"y" );
    if     ( i==0 ) frame[i]->SetTitle( "ee"        );
    else if( i==1 ) frame[i]->SetTitle( "#mu#mu   " );
    else if( i==2 ) frame[i]->SetTitle( "ee+#mu#mu" );
    data[i]->plotOn ( frame[i], Binning(xbin), LineWidth(1) );
    if( sel_fun==50 || sel_fun==51 || sel_fun==81 || sel_fun==82 ){
      modargus[i]->plotOn( frame[i], ProjWData(*data[i]), LineWidth(1) );      
      box[i] = new TPaveText( 0.05, 0.87, 0.35, 0.93,"BRNDC" );
      box[i]->AddText( Form("#chi^{2}/ndf = %f", frame[i]->chiSquare(nparam+N0bin[i])) );
      modargus[i] ->paramOn( frame[i], Format("NELU", AutoPrecision(2)), Layout(0.40, 0.99, 0.99) );
    }else if( sel_fun==15 || sel_fun==151 || sel_fun==181 || sel_fun==182 ){
      pdf[i] ->plotOn ( frame[i], ProjWData(*data[i]), LineWidth(1) );
      box[i] = new TPaveText( 0.05, 0.87, 0.35, 0.93,"BRNDC" );
      box[i]->AddText( Form("#chi^{2}/ndf = %f", frame[i]->chiSquare(nparam+N0bin[i])) );
      pdf[i] ->plotOn ( frame[i], ProjWData(*data[i]), Components(*gauss[i]), LineStyle(7), LineColor(2), LineWidth(1) );
      pdf[i] ->plotOn ( frame[i], ProjWData(*data[i]), Components(*modargus[i]), LineStyle(7), LineColor(4), LineWidth(1) );
      pdf[i] ->paramOn( frame[i], Format("NELU", AutoPrecision(2)), Layout(0.40, 0.99, 0.99) );
    }
    
    
    // -------------- calculation of pull -------------
    pullwaku[i] = new TH2D( Form("pull%d",    i), "data-func",             2, xmin, xmax, 2, -3, 3);
    pullhist[i] = new TH1D( Form("pullhist%d",i), Form("pullhist%d",i), xbin, xmin, xmax );
    RooAbsReal* int_pdf_bin = modargus[i]->createIntegral( *obs[i], NormSet(*obs[i]), Range("bin_int") );
    for( Int_t k=0; k<xbin; k++ ){
      obs[i]->setRange("bin_int", hist[i]->GetBinLowEdge(k+1), hist[i]->GetBinLowEdge(k+1) + hist[i]->GetBinWidth(k+1) );
      if( hist[i]->GetBinContent(k+1)==0 ) continue;
      pullhist[i]->SetBinContent( k+1, (hist[i]->GetBinContent(k+1) - (nsig[i]->getVal()+nbkg[i]->getVal()) * int_pdf_bin->getVal())/hist[i]->GetBinError(k+1) );
    }
    // ------------------------------------------------

    frame[i]->addObject( box[i] );
    c2->cd(i+1);
    frame[i]->Draw();

    c2->cd(i+4);
    pullwaku[i]->Draw();
    pullhist[i]->Draw("Psame");
  }

  c1->Update();
  c2->Update();

  if( flag_save ){
    c1->Print( Form("pic/%s_shape_func%d_s0%s_c1.eps",  axis, sel_fun, stream) );
    c2->Print( Form("pic/%s_shape_func%d_s0%s_c2.eps",  axis, sel_fun, stream) );
    TFile outfile( Form("pic/%s_shape_func%d_s0%s.root",  axis, sel_fun, stream), "RECREATE" );
    for( Int_t i=0; i<Nhist; i++ ) pullhist[i]->Write();
    outfile.Close();
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

