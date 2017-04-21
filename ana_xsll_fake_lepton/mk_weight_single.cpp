#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

#include "../Util/Style.h"
#include "../Util/Canvas.h"
#include "../Util/const.h"
#include "../Util/MChain.h"
#include "../Util/MCut.h"
#include "../Util/MCut_array.h"
#include "../Util/FitFunc.h"
#include "../Set_fake_lepton/Nominal_cut_selection.h"
#include "../Set_fake_lepton/Branch.h"
#include "../Set_fake_lepton/makeCut.h"

#include "draws_.h"

#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>
#include <TCut.h>
#include <TString.h>
#include <TFile.h>
#include <TH2D.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TString.h>
#include <TArrow.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TPaletteAxis.h>

Int_t main( Int_t argc, Char_t** argv ){
  if( argc!=6 ){
    std::cerr << "wrong input" << std::endl
	      << " Usage : ./merge_cut2 (char*)indir (char*)outdir (char*)exp (int)stream (char*)type"
	      << std::endl;
    abort();
  }

  const Char_t* indir  = argv[1];
  const Char_t* outdir = argv[2];
  const Char_t* exp    = argv[3];
  const Int_t   stream = atoi( argv[4] );
  const Char_t* type   = argv[5];
  // --------------------------------------------------------------------
  const Bool_t  fl_message    = !true;
  const Double_t low_fakerate = 0.0;
  // --------------------------------------------------------------------
  const Char_t* tname_B    = "h511";
  TFile* file_fl  = new TFile( "fakerate_dstr/fake_rate.root" );
  TH2D** hist_fl  = new TH2D*[2];
  TH2D** hist_flE = new TH2D*[2];
  hist_fl [0] =  (TH2D*)file_fl->Get( "hist2_fake" ); // mm
  hist_fl [1] =  (TH2D*)file_fl->Get( "hist3_fake" ); // ee
  hist_flE[0] =  (TH2D*)file_fl->Get( "hist2_err"  ); // mm
  hist_flE[1] =  (TH2D*)file_fl->Get( "hist3_err"  ); // ee
  TFile* file_fp_ee_svd1  = new TFile( "fakerate_twophoton/eff_table_eid_svd1.root"  );
  TFile* file_fp_ee_svd2  = new TFile( "fakerate_twophoton/eff_table_eid_svd2.root"  );
  TFile* file_fp_mm_svd1  = new TFile( "fakerate_twophoton/eff_table_muid_svd1.root" );
  TFile* file_fp_mm_svd2  = new TFile( "fakerate_twophoton/eff_table_muid_svd2.root" );

  TH2D** hist_fp_ee  = new TH2D*[4];
  hist_fp_ee [0] =  (TH2D*)file_fp_ee_svd1->Get( "hist1" ); // mm[SVD1]
  hist_fp_ee [1] =  (TH2D*)file_fp_ee_svd2->Get( "hist2" ); // mm[SVD2]

  TH2D** hist_fp_mm  = new TH2D*[4];
  hist_fp_mm [0] =  (TH2D*)file_fp_mm_svd1->Get( "hist1" ); // mm[SVD1]
  hist_fp_mm [1] =  (TH2D*)file_fp_mm_svd2->Get( "hist2" ); // mm[exp31-39  and   exp45a(run1-220)]
  hist_fp_mm [2] =  (TH2D*)file_fp_mm_svd2->Get( "hist3" ); // mm[exp41-49 except exp45a(run1-220)]
  hist_fp_mm [3] =  (TH2D*)file_fp_mm_svd2->Get( "hist4" ); // mm[exp51-73]
  // --------------------------------------------------------------------

  TChain* chain_B = new TChain( tname_B );
  
  std::stringstream sTmp;
  sTmp << indir << "/gMC_" << type << "_e0" << exp << "*s0" << stream << "*.root";

  Int_t nfile_B = chain_B->Add( sTmp.str().c_str() );
  if( !nfile_B ) std::cerr << "[ABORT] NO file : " << sTmp.str().c_str() << std::endl, abort();

  // --------------------------------------------------------------------  
  TCut cut_lrnb = "1";
  // --------------------------------------------------------------------
  const Char_t* outfile = Form( "%s/gMC_%s_e0%s_s0%d_caseB_single_weight.root",outdir,type, exp, stream );
  TFile* rootf  = new TFile( outfile, "RECREATE" );
  TTree* tree_B = new TTree();
  tree_B = chain_B->CopyTree( cut_lrnb );
  // ------------------------------------------------------------------------
  Float_t exprun, rm_l;                           // Mode
  Float_t epp, emp, mpp, mmp;                     // Momentum
  Float_t pi1p, pi2p, pi3p, pi4p;                 // Momentum
  Float_t lpc, lmc;                               // Direction
  Float_t pi1c, pi2c, pi3c, pi4c;                 // Direction
  Float_t lp_eid, lm_eid, lp_muid, lm_muid;       // Lepton-PID
  Float_t pi1_muid, pi2_muid, pi3_muid, pi4_muid; // Lepton-PID
  Float_t pi1_eid,  pi2_eid,  pi3_eid,  pi4_eid;  // Lepton-PID

  Float_t fl_p, fl_c;                       // Momentum and Direction of fake-lepton
  Float_t weight, weightE;                  // Weight ( New Branch )
  Float_t w_pi,   wE_pi;
  
  tree_B->SetBranchAddress( "exprun",   &exprun  );
  tree_B->SetBranchAddress( "rm_l",     &rm_l    );
  tree_B->SetBranchAddress( "epp",      &epp     );
  tree_B->SetBranchAddress( "emp",      &emp     );
  tree_B->SetBranchAddress( "mpp",      &mpp     );
  tree_B->SetBranchAddress( "mmp",      &mmp     );
  tree_B->SetBranchAddress( "lpc",      &lpc     );
  tree_B->SetBranchAddress( "lmc",      &lmc     );
  tree_B->SetBranchAddress( "lp_eid",   &lp_eid  );
  tree_B->SetBranchAddress( "lm_eid",   &lm_eid  );
  tree_B->SetBranchAddress( "lp_muid",  &lp_muid );
  tree_B->SetBranchAddress( "lm_muid",  &lm_muid );
  tree_B->SetBranchAddress( "pi1p",     &pi1p    );
  tree_B->SetBranchAddress( "pi2p",     &pi2p    );
  tree_B->SetBranchAddress( "pi3p",     &pi3p    );
  tree_B->SetBranchAddress( "pi4p",     &pi4p    );
  tree_B->SetBranchAddress( "pi1c",     &pi1c    );
  tree_B->SetBranchAddress( "pi2c",     &pi2c    );
  tree_B->SetBranchAddress( "pi3c",     &pi3c    );
  tree_B->SetBranchAddress( "pi4c",     &pi4c    );
  tree_B->SetBranchAddress( "pi1_eid",  &pi1_eid  );
  tree_B->SetBranchAddress( "pi2_eid",  &pi2_eid  );
  tree_B->SetBranchAddress( "pi3_eid",  &pi3_eid  );
  tree_B->SetBranchAddress( "pi4_eid",  &pi4_eid  );
  tree_B->SetBranchAddress( "pi1_muid", &pi1_muid );
  tree_B->SetBranchAddress( "pi2_muid", &pi2_muid );
  tree_B->SetBranchAddress( "pi3_muid", &pi3_muid );
  tree_B->SetBranchAddress( "pi4_muid", &pi4_muid );

  tree_B->Branch( "weight",  &weight,  "weight/F"  );
  tree_B->Branch( "weightE", &weightE, "weightE/F" );
  tree_B->Branch( "w_pi",    &w_pi,    "w_pi/F"    );
  tree_B->Branch( "wE_pi",   &wE_pi,   "wE_pi/F"   );

  TTree* newtree_B = tree_B->CloneTree( 0 );

  Int_t nev = 0;
  while (tree_B->GetEntry(nev, 0)) {
    nev++;
    if( rm_l == 1 ){ // electron mdoe
      if( lp_eid>=0.80 && epp>0.40 ){
	fl_p = epp;
	fl_c = lpc;
      }else if( lm_eid>=0.80 && emp>0.40 ){
	fl_p = emp;
	fl_c = lmc;
      }else{
	std::cerr << "[ABORT] : No fake-electron : "
		  << lp_eid << ", " << epp << ", "
		  << lm_eid << ", " << emp << std::endl, abort();
      }
    }else if( rm_l == 0 ){ // muon mode
      if( lp_muid>=0.97 && mpp>0.80 ){
	fl_p = mpp;
	fl_c = lpc;
      }else if( lm_muid>=0.97 && mmp>0.80 ){
	fl_p = mmp;
	fl_c = lmc;
      }else{
	std::cerr << "[ABORT] : No fake-muon : "
		  << lp_muid << ", " << mpp << ", "
		  << lm_muid << ", " << mmp << std::endl, abort();
      }
    }else{
      std::cerr << "[ABORT] : Invalid modes(rm_l) : " << rm_l << std::endl, abort();
    }
    if( fl_message ) std::cout << " ******************* [nev] " << nev << " / " << tree_B->GetEntries()
			       << " ******************* " << std::endl; // tmppppp
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Int_t sel_hist = -1;
    if( rm_l==1 ){ // electron
      if     ( exprun < 30*10000                            ) sel_hist = 0;
      else                                                    sel_hist = 1;
    }else if( rm_l==0 ){ // muon
      if     ( exprun < 30*10000                            ) sel_hist = 0;
      else if( exprun < 40*10000                            ) sel_hist = 1;
      else if( exprun > 45*10000 && exprun < 45*10000 + 221 ) sel_hist = 1;
      else if( exprun < 50*10000                            ) sel_hist = 2;
      else                                                    sel_hist = 3;
    }
    if( fl_message ) std::cout << "rm_l = " << rm_l << ", exprun = " << exprun << ", sel_hist = " << sel_hist << std::endl;
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    w_pi  = 0;
    wE_pi = 0;

    if( rm_l==1 ){
      if( fl_message ) std::cout << "[1st pion(electron)] p = " << pi1p << ", c = " << pi1c << std::endl;
      if( pi1p<0.4 ){
	w_pi  += low_fakerate;
	wE_pi += 0.0*0.0;
	if( fl_message ) std::cout << "low momentum" << std::endl;
      }else if( pi1_eid>0.80 ){
	Int_t tmp_x =  hist_fp_ee[sel_hist]->FindBin((Double_t)pi1c) - (hist_fp_ee[sel_hist]->GetNbinsX()+2); // "FindBin() does not work correctly (due to varable bin width ? )
	Int_t tmp_y = (hist_fp_ee[sel_hist]->FindBin((Double_t)pi1c, (Double_t)pi1p))/(hist_fp_ee[sel_hist]->GetNbinsX()+2);
	if( sel_hist==0 ){ // FindBin() works correctly to SVD1 histograms
	  tmp_x =  hist_fp_ee[sel_hist]->FindBin((Double_t)pi1c);
	  tmp_y = (hist_fp_ee[sel_hist]->FindBin((Double_t)pi1c, (Double_t)pi1p))/(hist_fp_ee[sel_hist]->GetNbinsX()+2);
	}
	if( fl_message ) std::cout << "eid : " << pi1_eid << ",  -> " << "( " << tmp_x << ", " << tmp_y << " ) : " << hist_fp_ee[sel_hist]->GetBinContent( tmp_x, tmp_y ) << std::endl;
	if( hist_fp_ee[sel_hist]->GetBinContent(tmp_x, tmp_y)!=0 ){
	  w_pi  += 1/hist_fp_ee[sel_hist]->GetBinContent( tmp_x, tmp_y ) -1; // (1-eff)/eff
	  wE_pi += hist_fp_ee[sel_hist]->GetBinError  ( tmp_x, tmp_y ) * hist_fp_ee[sel_hist]->GetBinError  ( tmp_x, tmp_y ) / hist_fp_ee[sel_hist]->GetBinContent( tmp_x, tmp_y ) / hist_fp_ee[sel_hist]->GetBinContent( tmp_x, tmp_y ) / hist_fp_ee[sel_hist]->GetBinContent( tmp_x, tmp_y ) / hist_fp_ee[sel_hist]->GetBinContent( tmp_x, tmp_y );
	}
      }
      
      if( fl_message ) std::cout << "[2nd pion(electron)] p = " << pi2p << ", c = " << pi2c << std::endl;
      if( pi2p<0.4 ){
	w_pi  += low_fakerate;
	wE_pi += 0.0*0.0;
	if( fl_message ) std::cout << "low momentum" << std::endl;
      }else if( pi2_eid>0.80 ){
	Int_t tmp_x =  hist_fp_ee[sel_hist]->FindBin((Double_t)pi2c) - (hist_fp_ee[sel_hist]->GetNbinsX()+2); // "FindBin() does not work correctly (due to varable bin width ? )
	Int_t tmp_y = (hist_fp_ee[sel_hist]->FindBin((Double_t)pi2c, (Double_t)pi2p))/(hist_fp_ee[sel_hist]->GetNbinsX()+2);
	if( sel_hist==0 ){ // FindBin() works correctly to SVD1 histograms
	  tmp_x =  hist_fp_ee[sel_hist]->FindBin((Double_t)pi2c);
	  tmp_y = (hist_fp_ee[sel_hist]->FindBin((Double_t)pi2c, (Double_t)pi2p))/(hist_fp_ee[sel_hist]->GetNbinsX()+2);
	}
	if( fl_message ) std::cout << "eid : " << pi2_eid << " -> " << "( " << tmp_x << ", " << tmp_y << " ) : " << hist_fp_ee[sel_hist]->GetBinContent( tmp_x, tmp_y ) << std::endl;
	if( hist_fp_ee[sel_hist]->GetBinContent(tmp_x, tmp_y)!=0 ){
	  w_pi  += 1/hist_fp_ee[sel_hist]->GetBinContent( tmp_x, tmp_y ) -1;
	  wE_pi += hist_fp_ee[sel_hist]->GetBinError  ( tmp_x, tmp_y ) * hist_fp_ee[sel_hist]->GetBinError  ( tmp_x, tmp_y ) / hist_fp_ee[sel_hist]->GetBinContent( tmp_x, tmp_y ) / hist_fp_ee[sel_hist]->GetBinContent( tmp_x, tmp_y ) / hist_fp_ee[sel_hist]->GetBinContent( tmp_x, tmp_y ) / hist_fp_ee[sel_hist]->GetBinContent( tmp_x, tmp_y );
	}
      }
      
      if( fl_message ) std::cout << "[3rd pion(electron)] p = " << pi3p << ", c = " << pi3c << std::endl;
      if( pi3p<0.4 ){
	w_pi  += low_fakerate;
	wE_pi += 0.0*0.0;
	if( fl_message ) std::cout << "low momentum" << std::endl;
      }else if( pi3_eid>0.80 ){
	Int_t tmp_x =  hist_fp_ee[sel_hist]->FindBin((Double_t)pi3c) - (hist_fp_ee[sel_hist]->GetNbinsX()+2); // "FindBin() does not work correctly (due to varable bin width ? )
	Int_t tmp_y = (hist_fp_ee[sel_hist]->FindBin((Double_t)pi3c, (Double_t)pi3p))/(hist_fp_ee[sel_hist]->GetNbinsX()+2);
	if( sel_hist==0 ){ // FindBin() works correctly to SVD1 histograms
	  tmp_x =  hist_fp_ee[sel_hist]->FindBin((Double_t)pi3c);
	  tmp_y = (hist_fp_ee[sel_hist]->FindBin((Double_t)pi3c, (Double_t)pi3p))/(hist_fp_ee[sel_hist]->GetNbinsX()+2);
	}
	if( fl_message ) std::cout << "eid : " << pi3_eid << " -> " << "( " << tmp_x << ", " << tmp_y << " ) : " << hist_fp_ee[sel_hist]->GetBinContent( tmp_x, tmp_y ) << std::endl;
	if( hist_fp_ee[sel_hist]->GetBinContent(tmp_x, tmp_y)!=0 ){
	  w_pi  += 1/hist_fp_ee[sel_hist]->GetBinContent( tmp_x, tmp_y ) -1;
	  wE_pi += hist_fp_ee[sel_hist]->GetBinError  ( tmp_x, tmp_y ) * hist_fp_ee[sel_hist]->GetBinError  ( tmp_x, tmp_y ) / hist_fp_ee[sel_hist]->GetBinContent( tmp_x, tmp_y ) / hist_fp_ee[sel_hist]->GetBinContent( tmp_x, tmp_y ) / hist_fp_ee[sel_hist]->GetBinContent( tmp_x, tmp_y ) / hist_fp_ee[sel_hist]->GetBinContent( tmp_x, tmp_y );
	}
      }
      
      if( fl_message ) std::cout << "[4th pion(electron)] p = " << pi4p  << ", c = " << pi4c << std::endl;
      if( pi4p<0.4 ){
	w_pi  += low_fakerate;
	wE_pi += 0.0*0.0;
	if( fl_message ) std::cout << "low momentum" << std::endl;
      }else if( pi4_eid>0.80 ){
	Int_t tmp_x =  hist_fp_ee[sel_hist]->FindBin((Double_t)pi4c) - (hist_fp_ee[sel_hist]->GetNbinsX()+2); // "FindBin() does not work correctly (due to varable bin width ? )
	Int_t tmp_y = (hist_fp_ee[sel_hist]->FindBin((Double_t)pi4c, (Double_t)pi4p))/(hist_fp_ee[sel_hist]->GetNbinsX()+2);
	if( sel_hist==0 ){ // FindBin() works correctly to SVD1 histograms
	  tmp_x =  hist_fp_ee[sel_hist]->FindBin((Double_t)pi4c);
	  tmp_y = (hist_fp_ee[sel_hist]->FindBin((Double_t)pi4c, (Double_t)pi4p))/(hist_fp_ee[sel_hist]->GetNbinsX()+2);
	}
	if( fl_message ) std::cout << "eid : " << pi4_eid << " -> " << "( " << tmp_x << ", " << tmp_y << " ) : " << hist_fp_ee[sel_hist]->GetBinContent( tmp_x, tmp_y ) << std::endl;
	if( hist_fp_ee[sel_hist]->GetBinContent(tmp_x, tmp_y)!=0 ){
	  w_pi  += 1/hist_fp_ee[sel_hist]->GetBinContent( tmp_x, tmp_y ) -1;
	  wE_pi += hist_fp_ee[sel_hist]->GetBinError  ( tmp_x, tmp_y ) * hist_fp_ee[sel_hist]->GetBinError  ( tmp_x, tmp_y ) / hist_fp_ee[sel_hist]->GetBinContent( tmp_x, tmp_y ) / hist_fp_ee[sel_hist]->GetBinContent( tmp_x, tmp_y ) / hist_fp_ee[sel_hist]->GetBinContent( tmp_x, tmp_y ) / hist_fp_ee[sel_hist]->GetBinContent( tmp_x, tmp_y );
	}
      }
    }else if( rm_l==0 ){
      if( fl_message ) std::cout << "[1st pion(muon)] p = " << pi1p << ", c = " << pi1c << std::endl;
      if( pi1p<0.8 ){
	w_pi  += low_fakerate;
	wE_pi += 0.0*0.0;
	if( fl_message ) std::cout << "low momentum" << std::endl;
      }else if( pi1_muid>0.97 ){
	Int_t tmp_x =  hist_fp_mm[sel_hist]->FindBin((Double_t)pi1c) - (hist_fp_mm[sel_hist]->GetNbinsX()+2); // "FindBin() does not work correctly (due to varable bin width ? )
	Int_t tmp_y = (hist_fp_mm[sel_hist]->FindBin((Double_t)pi1c, (Double_t)pi1p))/(hist_fp_mm[sel_hist]->GetNbinsX()+2);
	if( sel_hist==0 ){ // FindBin() works correctly to SVD1 histograms
	  tmp_x =  hist_fp_mm[sel_hist]->FindBin((Double_t)pi1c);
	  tmp_y = (hist_fp_mm[sel_hist]->FindBin((Double_t)pi1c, (Double_t)pi1p))/(hist_fp_mm[sel_hist]->GetNbinsX()+2);
	}
	if( fl_message ) std::cout << "muid : " << pi1_muid << " -> " << "( " << tmp_x << ", " << tmp_y << " ) : " << hist_fp_mm[sel_hist]->GetBinContent( tmp_x, tmp_y ) << std::endl;
	if( hist_fp_mm[sel_hist]->GetBinContent(tmp_x, tmp_y)!=0 ){
	  w_pi  += 1/hist_fp_mm[sel_hist]->GetBinContent( tmp_x, tmp_y ) -1; // (1-eff)/eff
	  wE_pi += hist_fp_mm[sel_hist]->GetBinError  ( tmp_x, tmp_y ) * hist_fp_mm[sel_hist]->GetBinError  ( tmp_x, tmp_y ) / hist_fp_mm[sel_hist]->GetBinContent( tmp_x, tmp_y ) / hist_fp_mm[sel_hist]->GetBinContent( tmp_x, tmp_y ) / hist_fp_mm[sel_hist]->GetBinContent( tmp_x, tmp_y ) / hist_fp_mm[sel_hist]->GetBinContent( tmp_x, tmp_y );
	}
      }
      
      if( fl_message ) std::cout << "[2nd pion(muon)] p = " << pi2p << ", c = " << pi2c << std::endl;
      if( pi2p<0.8 ){
	w_pi  += low_fakerate;
	wE_pi += 0.0*0.0;
	if( fl_message ) std::cout << "low momentum" << std::endl;
      }else if( pi2_muid>0.97 ){
	Int_t tmp_x =  hist_fp_mm[sel_hist]->FindBin((Double_t)pi2c) - (hist_fp_mm[sel_hist]->GetNbinsX()+2); // "FindBin() does not work correctly (due to varable bin width ? )
	Int_t tmp_y = (hist_fp_mm[sel_hist]->FindBin((Double_t)pi2c, (Double_t)pi2p))/(hist_fp_mm[sel_hist]->GetNbinsX()+2);
	if( sel_hist==0 ){ // FindBin() works correctly to SVD1 histograms
	  tmp_x =  hist_fp_mm[sel_hist]->FindBin((Double_t)pi2c);
	  tmp_y = (hist_fp_mm[sel_hist]->FindBin((Double_t)pi2c, (Double_t)pi2p))/(hist_fp_mm[sel_hist]->GetNbinsX()+2);
	}
	if( fl_message ) std::cout << "muid : " << pi2_muid << " -> " << "( " << tmp_x << ", " << tmp_y << " ) : " << hist_fp_mm[sel_hist]->GetBinContent( tmp_x, tmp_y ) << std::endl;
	if( hist_fp_mm[sel_hist]->GetBinContent(tmp_x, tmp_y)!=0 ){
	  w_pi  += 1/hist_fp_mm[sel_hist]->GetBinContent( tmp_x, tmp_y ) -1;
	  wE_pi += hist_fp_mm[sel_hist]->GetBinError  ( tmp_x, tmp_y ) * hist_fp_mm[sel_hist]->GetBinError  ( tmp_x, tmp_y ) / hist_fp_mm[sel_hist]->GetBinContent( tmp_x, tmp_y ) / hist_fp_mm[sel_hist]->GetBinContent( tmp_x, tmp_y ) / hist_fp_mm[sel_hist]->GetBinContent( tmp_x, tmp_y ) / hist_fp_mm[sel_hist]->GetBinContent( tmp_x, tmp_y );
	}
      }
      
      if( fl_message ) std::cout << "[3rd pion(muon)] p = " << pi3p << ", c = " << pi3c << std::endl;
      if( pi3p<0.8 ){
	w_pi  += low_fakerate;
	wE_pi += 0.0*0.0;
	if( fl_message ) std::cout << "low momentum" << std::endl;
      }else if( pi3_muid>0.97 ){
	Int_t tmp_x =  hist_fp_mm[sel_hist]->FindBin((Double_t)pi3c) - (hist_fp_mm[sel_hist]->GetNbinsX()+2); // "FindBin() does not work correctly (due to varable bin width ? )
	Int_t tmp_y = (hist_fp_mm[sel_hist]->FindBin((Double_t)pi3c, (Double_t)pi3p))/(hist_fp_mm[sel_hist]->GetNbinsX()+2);
	if( sel_hist==0 ){ // FindBin() works correctly to SVD1 histograms
	  tmp_x =  hist_fp_mm[sel_hist]->FindBin((Double_t)pi3c);
	  tmp_y = (hist_fp_mm[sel_hist]->FindBin((Double_t)pi3c, (Double_t)pi3p))/(hist_fp_mm[sel_hist]->GetNbinsX()+2);
	}
	if( fl_message ) std::cout << "muid : " << pi3_muid << " -> " << "( " << tmp_x << ", " << tmp_y << " ) : " << hist_fp_mm[sel_hist]->GetBinContent( tmp_x, tmp_y ) << std::endl;
	if( hist_fp_mm[sel_hist]->GetBinContent(tmp_x, tmp_y)!=0 ){
	  w_pi  += 1/hist_fp_mm[sel_hist]->GetBinContent( tmp_x, tmp_y ) -1;
	  wE_pi += hist_fp_mm[sel_hist]->GetBinError  ( tmp_x, tmp_y ) * hist_fp_mm[sel_hist]->GetBinError  ( tmp_x, tmp_y ) / hist_fp_mm[sel_hist]->GetBinContent( tmp_x, tmp_y ) / hist_fp_mm[sel_hist]->GetBinContent( tmp_x, tmp_y ) / hist_fp_mm[sel_hist]->GetBinContent( tmp_x, tmp_y ) / hist_fp_mm[sel_hist]->GetBinContent( tmp_x, tmp_y );
	}
      }
      
      if( fl_message ) std::cout << "[4th pion(muon)] p = " << pi4p << ", c = " << pi4c << std::endl;
      if( pi4p<0.8 ){
	w_pi  += low_fakerate;
	wE_pi += 0.0*0.0;
	if( fl_message ) std::cout << "low momentum" << std::endl;
      }else if( pi4_muid>0.97 ){
	Int_t tmp_x =  hist_fp_mm[sel_hist]->FindBin((Double_t)pi4c) - (hist_fp_mm[sel_hist]->GetNbinsX()+2); // "FindBin() does not work correctly (due to varable bin width ? )
	Int_t tmp_y = (hist_fp_mm[sel_hist]->FindBin((Double_t)pi4c, (Double_t)pi4p))/(hist_fp_mm[sel_hist]->GetNbinsX()+2);
	if( sel_hist==0 ){ // FindBin() works correctly to SVD1 histograms
	  tmp_x =  hist_fp_mm[sel_hist]->FindBin((Double_t)pi4c);
	  tmp_y = (hist_fp_mm[sel_hist]->FindBin((Double_t)pi4c, (Double_t)pi4p))/(hist_fp_mm[sel_hist]->GetNbinsX()+2);
	}
	if( fl_message ) std::cout << "muid : " << pi4_muid << " -> " << "( " << tmp_x << ", " << tmp_y << " ) : " << hist_fp_mm[sel_hist]->GetBinContent( tmp_x, tmp_y ) << std::endl;
	if( hist_fp_mm[sel_hist]->GetBinContent(tmp_x, tmp_y)!=0 ){
	  w_pi  += 1/hist_fp_mm[sel_hist]->GetBinContent( tmp_x, tmp_y ) -1;
	  wE_pi += hist_fp_mm[sel_hist]->GetBinError  ( tmp_x, tmp_y ) * hist_fp_mm[sel_hist]->GetBinError  ( tmp_x, tmp_y ) / hist_fp_mm[sel_hist]->GetBinContent( tmp_x, tmp_y ) / hist_fp_mm[sel_hist]->GetBinContent( tmp_x, tmp_y ) / hist_fp_mm[sel_hist]->GetBinContent( tmp_x, tmp_y ) / hist_fp_mm[sel_hist]->GetBinContent( tmp_x, tmp_y );
	}
      }
    }
    
    wE_pi = sqrt(wE_pi);
    if( fl_message ) std::cout << "weight(pion  ) = " << w_pi << " +- " << wE_pi << std::endl;
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    Int_t bin_x =  hist_fl[(Int_t)rm_l]->FindBin((Double_t)fl_c);
    Int_t bin_y = (hist_fl[(Int_t)rm_l]->FindBin((Double_t)fl_c, (Double_t)fl_p))/(hist_fl[(Int_t)rm_l]->GetNbinsX()+2);
    Double_t w_l = 0;
    Double_t wE_l = 0;
    w_l  = hist_fl [(Int_t)rm_l]->GetBinContent(bin_x, bin_y);
    wE_l = hist_flE[(Int_t)rm_l]->GetBinContent(bin_x, bin_y);
    if( bin_x==8 && bin_y==1 ){
      w_l  = (hist_fl [(Int_t)rm_l]->GetBinContent(7,1) + hist_fl [(Int_t)rm_l]->GetBinContent(9,1) )/2.0;
      wE_l = sqrt( hist_flE[(Int_t)rm_l]->GetBinContent(7,1)*hist_flE[(Int_t)rm_l]->GetBinContent(7,1)
		     +hist_flE[(Int_t)rm_l]->GetBinContent(9,1)*hist_flE[(Int_t)rm_l]->GetBinContent(9,1) )/2.0;
    }

    if( fl_message ) std::cout << "[ lepton ] p = " << fl_p << ", c = " << fl_c << " -> ( " << bin_x << ", " << bin_y << " )" << std::endl;


    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    weight  = w_l * w_pi;
    if( weight != 0 ) weightE = weight * sqrt( wE_l*wE_l/w_l/w_l + wE_pi*wE_pi/w_pi/w_pi);
    else              weightE = 0;
    if( fl_message ){
      std::cout << "weight(lepton) = " << w_l    << " +- " << wE_l    << std::endl;
      std::cout << "weight(total ) = " << weight << " +- " << weightE << std::endl;
    }

    
    /*    
    if( fl_message || (weight==0 && weightE==0) ) std::cout << "[CHECK] "
							    << std::setw( 7) << std::right << nev     << " : "
							    << std::setw( 2) << std::right << rm_l    << " : ("
							    << std::setw(12) << std::right << fl_c    << ", "
							    << std::setw(12) << std::right << fl_p    << ") -> ("
							    << std::setw( 3) << std::right << bin_x   << ", "
							    << std::setw( 3) << std::right << bin_y   << ") => "
							    << std::setw(12) << std::right << weight  << " +- "
							    << std::setw(12) << std::right << weightE << std::endl;
    */
    newtree_B->Fill();
  }

  std::cout << std::setw(7) << std::right << type
	    << ", s0"   << stream
	    << ", exp0" << std::setw(3) << std::left  << exp
	    << ", "     << std::setw(3) << std::right << nfile_B << " files(" << tname_B    << ", " << std::setw(9) << std::right << chain_B->GetEntries()
	    << " -> "
	    << std::setw(9) << std::right << newtree_B  ->GetEntries()
	    << "[" << std::setw(2) << std::right << int(100*((double)newtree_B->GetEntries()/chain_B->GetEntries())) << "%]"
	    << ")"
	    << std::endl;
  // --------------------------------------------------------------------

  newtree_B->Write();
  //rootf    ->Close();
  
  delete chain_B;
  delete tree_B;
  delete newtree_B;
  delete rootf;
  delete file_fl;

  return 0;
}
