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
using namespace gmc_emu;

using namespace Mbc_bkg;
//using namespace Mbc_bkg_wide;

const Bool_t flag_save = true; // outfile.eps and outfile.root
const Bool_t flag_fit  = true;
const Bool_t flag_roof = true; // flase(skip roofit)

const Bool_t   flag_fix_end = true; // true(fix endpoint parameter)
const Double_t endpoint     = 5.289;

const Bool_t   flag_k4pi     = true; // 1(veto  K4pi    modes)
const Bool_t   flag_unflavor = true; // 1(veto unflavor modes)

void manip_func( TF1* func ){
  if( flag_fix_end ) func->FixParameter( 1, endpoint );
  else               func->SetParLimits( 1, 5.280, 5.295 );
}

Int_t main( Int_t argc, Char_t** argv ){

  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==3 || argc==4) ) std::cerr << "wrong input" << std::endl
					<< " Usage : ./draws_1d_Mbc_fit (int)fl_nbin (int)fl_bin [(int)fl_appRun]"
					<< std::endl, abort();
  Char_t*  stream       = "0-5";
  Double_t used_nstream = 6;
  Int_t    sel_fun      = 50;
  Int_t    fl_appRun    = 1;
  Int_t    fl_nbin      = atoi( argv[1] ); // 8 or 14
  Int_t    fl_bin       = atoi( argv[2] ); // 8 -> [0-7,8], 14 -> [0-14], last is total
  if( argc==4 ) fl_appRun = atoi( argv[3] );
  
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t Nchain        = 3; // [gmc(ee,mm,emu)]
  const Int_t Ncategory     = 2; //[gmc(qq,non-peak)]
  const Int_t Nplot         = 2; //[bkg(qq,non-peak)]
  const Int_t Ntmp          = Ncategory*3; // x[ee,mm,emu]
  const Int_t Nhist         = Nplot    *3; // x[ee,mm,emu]
  const Int_t fl_sb[Ntmp]      = {0,0,1,1,2,2}; // 0(bkg,ee), 1(bkg,mm), 2(bkg,emu)
  const Int_t fl_mode_ll[Nchain] = {1,0,2};   // 1(e), 0(mu), 2(emu)

  const Int_t add[Nhist][Ntmp] ={
    {1,0}, // ee(qq)
    {1,1}, // ee(non-peak)
    {0,0,1,0}, // mm(qq)
    {0,0,1,1}, // mm(non-peak)
    {0,0,0,0,1,0}, // emu(qq)
    {0,0,0,0,1,1}, // emu(non-peak)
  };

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  Char_t** infile = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    infile[i] = new Char_t[1024];
    if     ( i==0 ) sTmp << indir[0] << "gMC_*_s0["    << stream << "]";  // bkg(ee)
    else if( i==1 ) sTmp << indir[0] << "gMC_*_s0["    << stream << "]";  // bkg(mm)
    else if( i==2 ) sTmp << indir[1] << "gMC_*_s0["    << stream << "]_*_emu0";  // bkg(emu)
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
      /*
      sTmp << " && genbfl!=0"
	   << " && !( (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && ( pi0org==1 || (pi0org<0.5&&rm_xs>999) || (pi0org<0.5&&rm_xs<999)) )" // charmonium events(right-pi0)
	   << " && !(!(lpgt==3 && lmgt==3 ) && (lpself==0 && lmself==0) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3) "  // double miss-id
           << " && !(!(lpgt==3 && lmgt==3 ) && ( (lpself==1&&lmself==0) || (lpself==0&&lmself==1) ) && ( lpgt==3 || lmgt==3 ) && rest_sw==0 && gm_bg1<3              && ( korg==lporg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg ) && lporg==lmorg)" // single miss-id (cc)
	   << " && !(!(lpgt==3 && lmgt==3 ) && ( (lpself==1&&lmself==0) || (lpself==0&&lmself==1) ) &&   lpgt!=3 && lmgt!=3   && rest_sw==0 && gm_bg1<3 && rm_xs>999 && ( korg==lporg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg ) && lporg==lmorg)"; // single miss-id with neutrino
      */
      sTmp << " && genbfl!=0 && !(lpgt==3 && lmgt==3) && lpself==1 && lmself==1"; // gmc(semi-leptonic b) // tmpppppp
      //sTmp << " && genbfl!=0"; // tmpppppp for emu
    }
    
    Int_t tmp_fl_mode_ll = 1;
    if( fl_sb[i]==1 ) tmp_fl_mode_ll = 0;
    if( fl_nbin == 14 ){
      if     ( fl_bin ==  0 ) sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 1,  1 ).c_str();
      else if( fl_bin ==  1 ) sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 1, -1 ).c_str();
      else if( fl_bin ==  2 ) sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 2,  1 ).c_str();
      else if( fl_bin ==  3 ) sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 2, -1 ).c_str();
      else if( fl_bin ==  4 ) sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 3,  1 ).c_str();
      else if( fl_bin ==  5 ) sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 3, -1 ).c_str();
      else if( fl_bin ==  6 ) sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 5,  1 ).c_str();
      else if( fl_bin ==  7 ) sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 5, -1 ).c_str();
      else if( fl_bin ==  8 ) sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 7,  1 ).c_str();
      else if( fl_bin ==  9 ) sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 7, -1 ).c_str();
      else if( fl_bin == 10 ) sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 8,  1 ).c_str();
      else if( fl_bin == 11 ) sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 8, -1 ).c_str();
      else if( fl_bin == 12 ) sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 9,  1 ).c_str();
      else if( fl_bin == 13 ) sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 9, -1 ).c_str();
    }else if( fl_nbin == 8 ){
      if     ( fl_bin ==  0 ){
	sTmp << " && ( " << makeCut_q2( tmp_fl_mode_ll, 1,  1 ).c_str()
	     << " || "
	     << makeCut_q2( tmp_fl_mode_ll, 2,  1 ).c_str()
	     << " )";
      }else if( fl_bin ==  1 ){
	sTmp << " && ( " << makeCut_q2( tmp_fl_mode_ll, 1, -1 ).c_str()
	     << " || "
	     << makeCut_q2( tmp_fl_mode_ll, 2, -1 ).c_str()
	     << " )";
      }else if( fl_bin ==  2 ){
	sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 3,  1 ).c_str();
      }else if( fl_bin ==  3 ){
	sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 3, -1 ).c_str();
      }else if( fl_bin ==  4 ){
	sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 5,  1 ).c_str();
      }else if( fl_bin ==  5 ){
	sTmp << " && " << makeCut_q2( tmp_fl_mode_ll, 5, -1 ).c_str();
      }else if( fl_bin ==  6 ){
	sTmp << " && ( " << makeCut_q2( tmp_fl_mode_ll, 7,  1 ).c_str()
	     << " || "   << makeCut_q2( tmp_fl_mode_ll, 8,  1 ).c_str()
	     << " || "   << makeCut_q2( tmp_fl_mode_ll, 9,  1 ).c_str()
	     << " )";
      }else if( fl_bin ==  7 ){
	sTmp << " && ( " << makeCut_q2( tmp_fl_mode_ll, 7, -1 ).c_str()
	     << " || "   << makeCut_q2( tmp_fl_mode_ll, 8, -1 ).c_str()
	     << " || "   << makeCut_q2( tmp_fl_mode_ll, 9, -1 ).c_str()
	     << " )";
      }
    }else std::cerr << "[Abort] Wrong fl_nbin : " << fl_nbin << std::endl, abort();
    

    if( flag_k4pi     ) sTmp << " && " << makeCut_5body_veto().c_str();
    if( flag_unflavor ) sTmp << " && " << makeCut_unflavor_veto().c_str();
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
    chain[j] = new MChain( infile[j], tname, branch_table(), 0, tail );
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
  }
  
  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make hist *************************************" << std::endl;
  for( Int_t i=0; i<Nhist; i++ ){
    hist[i] = new TH1D( Form("hist%d",i),Form("hist%d",i), xbin,offset+xmin,offset+xmax );
    Deco( hist[i], 3, col_fil[i%Nplot+1], col_fil[i%Nplot+1] );
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
      if( sel_fun==50 ){
	func[i]->SetParNames  ("norm","Ebeam","shape");
	func[i]->SetParameters( 100*hist[i]->GetBinContent(1), endpoint, -20);
	if( flag_fix_end ) func[i]->FixParameter( 1, endpoint );
	else               func[i]->SetParLimits( 1, 5.285, 5.290 );
      }else func_set_parameters(sel_fun, func[i], hist[i], xbin, offset+xmin, offset+xmax);
      func[i]->SetLineColor(2);
    }
  }
  

  // +++++++ draw ++++++++++++++++++++++++++++++++++
  c1->Draw();
  TH2D** waku = new TH2D*[Nhist];
  waku[0] = Waku( Nhist, hist, Form("%s (ee)",     xlabel) );
  waku[1] = Waku( Nhist, hist, Form("%s (#mu#mu)", xlabel) );
  waku[2] = Waku( Nhist, hist, Form("%s (e#mu)",   xlabel) );
  
  for(Int_t i=Nhist-1; i>=0; i-- ){
    c1->cd(i/Nplot+1);
    if     ( 0*Nplot <= i && i < 1*Nplot ) hist[i]->SetTitle( Form("ee"    ) );
    else if( 1*Nplot <= i && i < 2*Nplot ) hist[i]->SetTitle( Form("#mu#mu") );
    else if( 2*Nplot <= i && i < 3*Nplot ) hist[i]->SetTitle( Form("e#mu"  ) );

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
    c1->Print( Form("pic/%s_fit_shape_%dbin_%d_func%d_s0%s_c1.eps", axis, fl_nbin, fl_bin, sel_fun, stream) );
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
  const Int_t    Nroohist        = 3; // ee, mm , emu
  const Int_t    rooad[Nroohist] = {1, 3, 5};
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
	for( Int_t i=0; i<Nroohist; i++ ){
	  obs[i]->setVal( x_obs );
	  if( add[rooad[i]][j] ){
	    data[i]->add( RooArgSet(*obs[i]) );
	    if( i==Nroohist-1 ) evt_cnt[j]++;
	  }
	}
      }
    }
  }
  
  for( Int_t i=0; i<Nroohist; i++ ) std::cout << Form("< data%d> nevt = ", i ) << data[i]->numEntries()
					      << "( " << data[i]->sumEntries() << " )"
					      << std::endl; // sum of weight
  // ------------------------- Fit Fucntion ----------------------------

  const Int_t nparam = 2; // # of fit parameters

  // Bkg-PDF
  RooRealVar** arg_end   = new RooRealVar*[Nroohist];
  RooRealVar** arg_shape = new RooRealVar*[Nroohist];
  RooAbsPdf**  modargus  = new RooAbsPdf* [Nroohist];

  for( Int_t i=0; i<Nroohist; i++ ){
    arg_end  [i] = new RooRealVar( "end",   "argus endpoint",        endpoint, 5.285, 5.290 );
    if( flag_fix_end ) arg_end[i]->setConstant(kTRUE);
    //Double_t tmp_shape = func[rooad[i]]->GetParameter(2);
    Double_t tmp_shape = -20;
    arg_shape[i] = new RooRealVar( "shape", "argus shape parameter", tmp_shape, tmp_shape-80.0, tmp_shape+30.0 );
    modargus [i] = bindPdf( Form("modargus%d",i), func_roo_argus,      *obs[i],   *arg_end[i], *arg_shape[i] );
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Total-PDF
  RooFitResult** fit_result = new RooFitResult*[Nroohist];
  for( Int_t i=0; i<Nroohist; i++ ){
    fit_result[i] = modargus[i]->fitTo( *data[i], Range(offset+xmin, endpoint), Save(true) );
  }
  // ------------------------- Draw ----------------------------

  RooPlot**    frame = new RooPlot*  [Nroohist];
  TPaveText**  box   = new TPaveText*[Nroohist];

  for( Int_t i=0; i<Nroohist; i++ ){
    frame[i] = obs[i]->frame();
    frame[i]->GetXaxis()->CenterTitle();
    frame[i]->GetYaxis()->CenterTitle();
    frame[i]->SetTitleOffset( 1.00,"x" );
    frame[i]->SetTitleOffset( 1.30,"y" );
    if     ( i==0 ) frame[i]->SetTitle( "ee"     );
    else if( i==1 ) frame[i]->SetTitle( "#mu#mu" );
    else if( i==2 ) frame[i]->SetTitle( "e#mu"   );
    data[i]->plotOn ( frame[i], Binning(xbin), LineWidth(1) );
    modargus[i]->plotOn( frame[i], ProjWData(*data[i]), LineWidth(1) );      
    box[i] = new TPaveText( 0.05, 0.87, 0.35, 0.93,"BRNDC" );
    box[i]->AddText( Form("#chi^{2}/ndf = %f", frame[i]->chiSquare(nparam+N0bin[i])) );
    modargus[i] ->paramOn( frame[i], Format("NELU", AutoPrecision(2)), Layout(0.40, 0.99, 0.99) );
    
    frame[i]->addObject( box[i] );
    c2->cd(i+1);
    frame[i]->Draw();
  }
  
  c1->Update();
  c2->Update();
  
  if( flag_save ){
    c1->Print( Form("pic/%s_fit_shape_%dbin_%d_func%d_s0%s_c1.eps",  axis, fl_nbin, fl_bin, sel_fun, stream) );
    c2->Print( Form("pic/%s_fit_shape_%dbin_%d_func%d_s0%s_c2.eps",  axis, fl_nbin, fl_bin, sel_fun, stream) );
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
  delete[] arg_end;
  delete[] arg_shape;
  delete[] modargus;
  delete[] fit_result;
  delete[] frame;
  delete[] box;
  delete   c2;

  return 0;
}

