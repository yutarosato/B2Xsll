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
  //using namespace sig_gmc_rd_cut2_beforebgsup;
  using namespace sig_gmc_rd_cut2;

  using namespace Mbc_bkg;
  using namespace RooFit;
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==5 || argc==6) ) std::cerr << "wrong input" << std::endl
    					<< " Usage : ./draws_1d_Mbc (char*)stream (char*) setname (double)used_nstream (double)used_nset [(int)fl_appRun]" << std::endl
					<< std::endl, abort();
  Char_t*  stream        = argv[1];
  Char_t*  setname       = argv[2];
  Double_t used_nstream = atof(argv[3]); 
  Double_t used_nset    = atof(argv[4]); 
  Int_t   fl_appRun     = 1;
  if( argc==6 ) fl_appRun = atoi( argv[5] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Bool_t flag_k4pi      = true; // 1(veto  K4pi    modes)
  const Bool_t flag_unflavor  = true; // 1(veto unflavor modes)
  const Bool_t flag_ccpi0     = !true; // 1(veto ccpi0 peak    )
  const Bool_t flag_save      = true; // outfile.eps and outfile.root
  const Bool_t flag_fit       = !true;
  const Int_t  sel_fun        = 15;   // 15(gauss+argus)
  const Int_t  Nfun           = 6;    // [MC,RD] x [ee,mm,ee+mm]
  const Int_t  sel_hist[Nfun] = {4,10,16,5,11,17};
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t    Nchain             = 6*2; // [gmc(qq,peak-bb,non-peak-bb),sigmc(true,scf),rd] x [ee,mm]
  const Int_t    Nhist              = 6*3; // [gmc(qq,peak-bb,non-peak-bb),sigmc(true,scf),rd] x [ee,mm,ee+mm]
  const Int_t    nfile[Nchain]      = {0};
  const Int_t    fl_sb[Nchain]      = {0,0,0,1,1,2,
				       0,0,0,1,1,2}; // 0(bkg), 1(sig), 2(rd)
  const Int_t    fl_mode_ll[Nchain] = {1,1,1,1,1,1,  // 1(ee)
				       0,0,0,0,0,0}; // 0(mm)
  const Double_t scale_event_sig = sigmc_amount*(used_nset/(Double_t)nset); // sigmc : N -> N/alpha
  const Double_t scale_event_bkg = used_nstream;                            //   gmc : N -> N/alpha
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Char_t* LRNB_cut = new Char_t[4096];
  LRNB_cut = "1";
  //LRNB_cut = "((lr1_bb>0.93 && lr1_qq>0.86 && xs_m>1.1) || (lr1_bb>0.56 && lr1_qq>0.83 && xs_m<1.1))";
  //LRNB_cut = "((rm_l==1 && lr1_bb>0.94 && lr1_qq>0.89 && xs_m>1.1) || (rm_l==1 && lr1_bb>0.70 && lr1_qq>0.88 && xs_m<1.1) || (rm_l==0 && lr1_bb>0.94 && lr1_qq>0.89 && xs_m>1.1) || (rm_l==0 && lr1_bb>0.70 && lr1_qq>0.88 && xs_m<1.1))";
  std::cout << "LRNB cut : " << LRNB_cut << std::endl;
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
  const Int_t add[Nhist][Nchain] ={
    {1},           // gmc(qq),          ee
    {1,1},         // gmc(non-peak-bb), ee
    {1,1,1},       // gmc(peak-bb),     ee
    {1,1,1,1},     // sigmc(scf),       ee
    {1,1,1,1,1},   // sigmc(true),      ee
    //{1,1,1,0,0},   // no-sigmc,         emu
    //{1,1,1,0,0},   // no-sigmc,         emu
    {0,0,0,0,0,1}, // rd,      ee
    {0,0,0,0,0,0,1},           // gmc(qq),          mm
    {0,0,0,0,0,0,1,1},         // gmc(non-peak-bb), mm
    {0,0,0,0,0,0,1,1,1},       // gmc(peak-bb),     mm
    {0,0,0,0,0,0,1,1,1,1},     // sigmc(scf)        mm
    {0,0,0,0,0,0,1,1,1,1,1},   // sigmc(true)       mm
    //{0,0,0,0,0,0,1,1,1,0,0},   // no-sigmc          emu
    //{0,0,0,0,0,0,1,1,1,0,0},   // no-sigmc          emu
    {0,0,0,0,0,0,0,0,0,0,0,1}, // rd,               mm
    {1,0,0,0,0,0,1},             // gmc(qq),          ee+mm
    {1,1,0,0,0,0,1,1},           // gmc(non-peak-bb), ee+mm
    {1,1,1,0,0,0,1,1,1},         // gmc(peak-bb),     ee+mm
    {1,1,1,1,0,0,1,1,1,1},       // sigmc(scf),       ee+mm
    {1,1,1,1,1,0,1,1,1,1,1},     // sigmc(true),      ee+mm
    //{1,1,1,0,0,0,1,1,1,0,0},     // no-sigmc(true),      ee+mm
    //{1,1,1,0,0,0,1,1,1,0,0},     // no-sigmc(true),      ee+mm
    {0,0,0,0,0,1,0,0,0,0,0,1},   // rd,               ee+mm
  };

  Char_t** add_cut = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ){
    add_cut[i] = new Char_t[4096];
    sTmp << LRNB_cut;
    if( flag_k4pi     ) sTmp << " && " << makeCut_5body_veto().c_str();
    if( flag_unflavor ) sTmp << " && " << makeCut_unflavor_veto().c_str();
    if( flag_ccpi0    ) sTmp << " && " << cut_ccpi0;

    if     ( i==0 || i== 6 ) sTmp << " && genbfl==0"; // gmc(qq)
    else if( i==1 || i== 7 ){
      sTmp << " && genbfl!=0"; // gmc(non-peak-bb)
      sTmp << " && !( (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && ( pi0org==1 || (pi0org<0.5&&rm_xs>999) || (pi0org<0.5&&rm_xs<999)) )" // charmonium events(right-pi0)
	     << " && !(!(lpgt==3 && lmgt==3 ) && (lpself==0 && lmself==0) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3) "  // double miss-id
	     << " && !(!(lpgt==3 && lmgt==3 ) && ( (lpself==1&&lmself==0) || (lpself==0&&lmself==1) ) && ( lpgt==3 || lmgt==3 ) && rest_sw==0 && gm_bg1<3              && ( korg==lporg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg ) && lporg==lmorg)" // single miss-id (cc)
	     << " && !(!(lpgt==3 && lmgt==3 ) && ( (lpself==1&&lmself==0) || (lpself==0&&lmself==1) ) &&   lpgt!=3 && lmgt!=3   && rest_sw==0 && gm_bg1<3 && rm_xs>999 && ( korg==lporg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg ) && lporg==lmorg)"; // single miss-id with neutrino
    }else if( i==2 || i== 8 ){
      sTmp << " && genbfl!=0"; // gmc(peak-bb)
      sTmp << " && ("
	   << "( (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && ( pi0org==1 || (pi0org<0.5&&rm_xs>999) || (pi0org<0.5&&rm_xs<999)) )" // charmonium events(right-pi0)
	   << " || (!(lpgt==3 && lmgt==3 ) && (lpself==0 && lmself==0) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3) "  // double miss-id
	   << " || (!(lpgt==3 && lmgt==3 ) && ( (lpself==1&&lmself==0) || (lpself==0&&lmself==1) ) && ( lpgt==3 || lmgt==3 ) && rest_sw==0 && gm_bg1<3              && ( korg==lporg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg ) && lporg==lmorg)" // single miss-id (cc)
	   << " || (!(lpgt==3 && lmgt==3 ) && ( (lpself==1&&lmself==0) || (lpself==0&&lmself==1) ) &&   lpgt!=3 && lmgt!=3   && rest_sw==0 && gm_bg1<3 && rm_xs>999 && ( korg==lporg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg ) && lporg==lmorg)" // single miss-id with neutrino
	   << " ) ";
    }else if( i==3 || i== 9 ) sTmp << " && self!=1";   // sigmc(scf)
    else if( i==4 || i==10 ) sTmp << " && self==1";   // sigmc(true)

    
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
      sTmp << "&& ( " // Kll
      << Form( " (rm_xs==  1) || " ) // K+
      << Form( " (rm_xs== 10)    " ) // Ks
      << " )";
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
    if     ( fl_sb[j]==1 ) tmphist[j]->Scale( 1/scale_event_sig );
    else if( fl_sb[j]==0 ) tmphist[j]->Scale( 1/scale_event_bkg );
  }
  // +++++++ make hist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ make hist *************************************" << std::endl;
  for( Int_t i=0; i<Nhist; i++ ){
    hist[i] = new TH1D( Form("hist%d",i),Form("hist%d",i), xbin,offset+xmin,offset+xmax );
    //if( i!=3 && i!=7 && i!=11 ) Deco( hist[i], 3, col_fil[i%4], col_fil[i%4] );
    Deco( hist[i], 3, col_fil[i%6], col_fil[i%6] );
    if     ( i%6==0 ) Deco( hist[i], 3, 4, 4 ); // qq
    else if( i%6==1 ) Deco( hist[i], 3, 9, 9 ); // non-peak-bb
    else if( i%6==2 ) Deco( hist[i], 3, 6, 6 ); // peak-bb
    else if( i%6==3 ) Deco( hist[i], 3, 7, 7 ); // scf
    else if( i%6==4 ) Deco( hist[i], 3, 3, 3 ); // true
    else if( i%6==5 ) Deco( hist[i], 3, 2, 2 ); // rd

    for( Int_t j=0; j<Nchain; j++ ){
      if( add[i][j] ) hist[i]->Add( tmphist[j] );
    }
  }
  hist[ 0]->SetTitle("gMC(qq),          ee"       );
  hist[ 1]->SetTitle("gMC(non-peak-bb), ee"       );
  hist[ 2]->SetTitle("gMC(peak-bb),     ee"       );
  hist[ 3]->SetTitle("siMC(scf),        ee"       );
  hist[ 4]->SetTitle("siMC(true),       ee"       );
  hist[ 5]->SetTitle("RD,               ee"       );
  hist[ 6]->SetTitle("gMC(qq),          #mu#mu"   );
  hist[ 7]->SetTitle("gMC(non-peak-bb), #mu#mu"   );
  hist[ 8]->SetTitle("gMC(peak-bb),     #mu#mu"   );
  hist[ 9]->SetTitle("siMC(scf),        #mu#mu"   );
  hist[10]->SetTitle("siMC(true),       #mu#mu"   );
  hist[11]->SetTitle("RD,               #mu#mu"   );
  hist[12]->SetTitle("gMC(qq),          ee+#mu#mu");
  hist[13]->SetTitle("gMC(non-peak-bb), ee+#mu#mu");
  hist[14]->SetTitle("gMC(peak-bb),     ee+#mu#mu");
  hist[15]->SetTitle("siMC(scf),        ee+#mu#mu");
  hist[16]->SetTitle("siMC(true),       ee+#mu#mu");
  hist[17]->SetTitle("RD,               ee+#mu#mu");

  // +++++++ display ++++++++++++++++++++++++++++++++++
  std::cout << "=================================" << std::endl
	    << "[Scale Factor]"                    << std::endl
    	    << "sig : " << scale_event_sig
	    << "("      << used_nset << " set)"    << std::endl
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
  TH2D** waku = new TH2D*[3];
  c1->Draw();
  c1->cd(1); waku[0] = Waku( Nhist/3, &hist[ 0], xlabel, "ee", "ee"        ); waku[0]->Draw();
  c1->cd(2); waku[1] = Waku( Nhist/3, &hist[ 6], xlabel, "mm", "#mu#mu"    ); waku[1]->Draw();
  c1->cd(3); waku[2] = Waku( Nhist/3, &hist[12], xlabel, "ll", "ee+#mu#mu" ); waku[2]->Draw();
  
  for(Int_t i=Nhist/3-2; i>=0; i-- ){
    c1->cd(3);
    //if( i==Nhist/3-2 ) hist[i+2*Nhist/3]->Draw( "hist");
    //else               hist[i+2*Nhist/3]->Draw( "hist same" );
    hist[i+2*Nhist/3]->Draw( "hist same" );
    c1->cd(2);
    //if( i==Nhist/3-2 ) hist[i+1*Nhist/3]->Draw( "hist" );
    //else               hist[i+1*Nhist/3]->Draw( "hist same" );
    hist[i+1*Nhist/3]->Draw( "hist same" );
    c1->cd(1);
    //if( i==Nhist/3-2 ) hist[i+0*Nhist/3]->Draw( "hist" );
    //else               hist[i+0*Nhist/3]->Draw( "hist same" );
    hist[i+0*Nhist/3]->Draw( "hist same" );
  }

  // rd
  c1->cd(1); hist[ 5]->SetLineColor(2); hist[ 5]->Draw("PE0same");
  c1->cd(2); hist[11]->SetLineColor(2); hist[11]->Draw("PE0same");
  c1->cd(3); hist[17]->SetLineColor(2); hist[17]->Draw("PE0same");

  // +++++++ tlegend1 ++++++++++++++++++++++++++++++++++
  TLegend* legend1 = new TLegend( 0.75,0.75,0.99,0.99 );
  legend1->AddEntry( hist[5],"rd",                "P" );
  legend1->AddEntry( hist[4],"sigmc(true)",       "P" );
  legend1->AddEntry( hist[3],"sigmc(scf)",        "P" );
  legend1->AddEntry( hist[2],"gmc.(peak-bb)",     "P" );
  legend1->AddEntry( hist[1],"gmc.(non-peak-bb)", "P" );
  legend1->AddEntry( hist[0],"gmc.(qq)",          "P" );
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

  /*
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++ RooFit ++++++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // how to deal with scale factor ?????????
  
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
	  //if( add[sel_hist[i]][j] ){ // tmppppppppp
	    //if     ( fl_sb[j]==0 ) data[i]->add( RooArgSet(*obs[i]), 1/scale_event_bkg );
	  //else if( fl_sb[j]==1 ) data[i]->add( RooArgSet(*obs[i]), 1/scale_event_sig );
	  //else if( fl_sb[j]==2 ) data[i]->add( RooArgSet(*obs[i])                    );
	  //}
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

  RooPlot**    frame       = new RooPlot*  [Nfun];
  TPaveText**  box         = new TPaveText*[Nfun];

  for( Int_t i=0; i<Nfun; i++ ){
    c2->cd(i+1);
    frame[i] = obs[i]->frame();
    frame[i]->GetXaxis()->CenterTitle();
    frame[i]->GetYaxis()->CenterTitle();
    frame[i]->SetTitleOffset( 1.00,"x" );
    frame[i]->SetTitleOffset( 1.30,"y" );
    frame[i]->SetTitle( hist_fit[i]->GetTitle() );
    data[i]->plotOn ( frame[i], Binning(xbin),       LineWidth(1), Rescale( 1/used_nstream ) );
    pdf[i] ->plotOn ( frame[i], ProjWData(*data[i]), LineWidth(1), Rescale( 1/used_nstream ) );
    box[i] = new TPaveText( 0.05, 0.87, 0.35, 0.93,"BRNDC" );
    box[i]->AddText( Form("#chi^{2}/ndf = %f", frame[i]->chiSquare(nparam+N0bin[i])) );
    pdf[i] ->plotOn ( frame[i], ProjWData(*data[i]), Components(*gauss[i]), LineStyle(7), LineColor(2), LineWidth(1), Rescale( 1/used_nstream ) );
    pdf[i] ->plotOn ( frame[i], ProjWData(*data[i]), Components(*argus[i]), LineStyle(7), LineColor(4), LineWidth(1), Rescale( 1/used_nstream ) );
    pdf[i] ->paramOn( frame[i], Format("NELU", AutoPrecision(2)), Layout(0.40, 0.99, 0.99) );

    frame[i]->addObject( box[i] );
    frame[i]->Draw();
  }
  */  
  c1->Update();
  //c2->Update();
  if( flag_save ){
    c1->Print( Form("pic/%s_sep_peak_scf_comb_s0%s_set%s_c1.eps",  axis, stream, setname) );
    //c2->Print( Form("pic/%s_sep_peak_scf_comb_s0%s_set%s_c2.eps",  axis, stream, setname) );
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
  /*
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
  */
  return 0;
}

