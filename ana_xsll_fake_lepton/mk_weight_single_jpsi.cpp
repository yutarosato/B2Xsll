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


// --------------------------------------------------------------------
const Bool_t  fl_message = !true;
// --------------------------------------------------------------------
const Char_t* tname_B    = "h511";
TFile* file_fl  = new TFile( "fakerate_dstr/fake_rate.root" );
TH2D** hist_fl  = new TH2D*[2];
TH2D** hist_flE = new TH2D*[2];
TFile* file_fp_ee_svd1  = new TFile( "fakerate_twophoton/eff_table_eid_svd1.root"  );
TFile* file_fp_ee_svd2  = new TFile( "fakerate_twophoton/eff_table_eid_svd2.root"  );
TFile* file_fp_mm_svd1  = new TFile( "fakerate_twophoton/eff_table_muid_svd1.root" );
TFile* file_fp_mm_svd2  = new TFile( "fakerate_twophoton/eff_table_muid_svd2.root" );
TH2D** hist_fp_ee  = new TH2D*[4];
TH2D** hist_fp_mm  = new TH2D*[4];
TFile* file_kpi_svd1  = new TFile( "eff_table_kpi/table_eff_fake_kid_svd1.root"  );
TFile* file_kpi_svd2  = new TFile( "eff_table_kpi/table_eff_fake_kid_svd2.root"  );
TH2D** hist_kpi   = new TH2D*[2];

void Get_eff_lepton( Double_t c, Double_t p, Int_t fl_mode_ll, Int_t sel_hist, Double_t& eff, Double_t& effE ){
  hist_fp_ee [0] =  (TH2D*)file_fp_ee_svd1->Get( "hist1" ); // mm[SVD1]
  hist_fp_ee [1] =  (TH2D*)file_fp_ee_svd2->Get( "hist2" ); // mm[SVD2]
  
  hist_fp_mm [0] =  (TH2D*)file_fp_mm_svd1->Get( "hist1" ); // mm[SVD1]
  hist_fp_mm [1] =  (TH2D*)file_fp_mm_svd2->Get( "hist2" ); // mm[exp31-39  and   exp45a(run1-220)]
  hist_fp_mm [2] =  (TH2D*)file_fp_mm_svd2->Get( "hist3" ); // mm[exp41-49 except exp45a(run1-220)]
  hist_fp_mm [3] =  (TH2D*)file_fp_mm_svd2->Get( "hist4" ); // mm[exp51-73]
  
  Int_t xbin, ybin;
  if( fl_mode_ll==1 ){ // ee
    // e+
    if( sel_hist ) sel_hist=1;
    if( sel_hist==0 ){ // FindBin() works correctly to SVD1 histograms
      xbin =  hist_fp_ee[sel_hist]->FindBin(c);
      ybin = (hist_fp_ee[sel_hist]->FindBin(c, p))/(hist_fp_ee[sel_hist]->GetNbinsX()+2);
    }else{
      xbin =  hist_fp_ee[sel_hist]->FindBin(c) - (hist_fp_ee[sel_hist]->GetNbinsX()+2); // "FindBin() does not work correctly (due to varable bin width ? )
      ybin = (hist_fp_ee[sel_hist]->FindBin(c, p))/(hist_fp_ee[sel_hist]->GetNbinsX()+2);
    }
    eff  = hist_fp_ee[sel_hist]->GetBinContent( xbin, ybin );
    effE = hist_fp_ee[sel_hist]->GetBinError  ( xbin, ybin );
  }else{
    // mu+
    if( sel_hist==0 ){ // FindBin() works correctly to SVD1 histograms
      xbin =  hist_fp_mm[sel_hist]->FindBin(c);
      ybin = (hist_fp_mm[sel_hist]->FindBin(c, p))/(hist_fp_mm[sel_hist]->GetNbinsX()+2);
    }else{
      xbin =  hist_fp_mm[sel_hist]->FindBin(c) - (hist_fp_mm[sel_hist]->GetNbinsX()+2); // "FindBin() does not work correctly (due to varable bin width ? )
      ybin = (hist_fp_mm[sel_hist]->FindBin(c, p))/(hist_fp_mm[sel_hist]->GetNbinsX()+2);
    }
    eff  = hist_fp_mm[sel_hist]->GetBinContent( xbin, ybin );
    effE = hist_fp_mm[sel_hist]->GetBinError  ( xbin, ybin );
  }
  if( fl_message ) std::cout << "     <Get_eff_lepton (rm_l=" << fl_mode_ll << ", sel_hist=" << sel_hist << ")> "
			     << "( " << c << ", " << p << " ) -> (" << xbin << ", " << ybin << " ) -> "
			     << eff << " +- " << effE << std::endl;
  return;
}

void Cal_weight_lepton( Double_t eff, Double_t effE, Double_t& weight, Double_t& weightE ){ // Calculate (1-eff)/eff
  weight  = (eff==0 ? 0 : 1/eff -1    ); // (1-eff)/eff
  weightE = (eff==0 ? 0 : effE/eff/eff);
}

void Get_eff_pi( Double_t c, Double_t p, Int_t sel_hist, Double_t& eff, Double_t& effE ){
  hist_kpi[0] =  (TH2D*)file_kpi_svd1->Get( "eff_pi_4" ); // [SVD1]
  hist_kpi[1] =  (TH2D*)file_kpi_svd2->Get( "eff_pi_4" ); // [SVD2]
  Int_t xbin =  hist_kpi[sel_hist]->FindBin(c) - (hist_kpi[sel_hist]->GetNbinsX()+2); // "FindBin() does not work correctly (due to varable bin width ? )
  Int_t ybin = (hist_kpi[sel_hist]->FindBin(c, p))/(hist_kpi[sel_hist]->GetNbinsX()+2);
  eff  = hist_kpi[sel_hist]->GetBinContent( xbin, ybin );
  effE = hist_kpi[sel_hist]->GetBinError  ( xbin, ybin );
  if( fl_message ) std::cout << "     <Get_eff_pi (sel_hist=" << sel_hist << ")> "
			     << "( " << c << ", " << p << " ) -> (" << xbin << ", " << ybin << " ) -> "
			     << eff << " +- " << effE << std::endl;
}

void Get_fake_pi( Double_t c, Double_t p, Int_t fl_mode_ll, Double_t& fake, Double_t& fakeE ){
  hist_fl [0] =  (TH2D*)file_fl->Get( "hist2_fake" ); // mm
  hist_fl [1] =  (TH2D*)file_fl->Get( "hist3_fake" ); // ee
  hist_flE[0] =  (TH2D*)file_fl->Get( "hist2_err"  ); // mm
  hist_flE[1] =  (TH2D*)file_fl->Get( "hist3_err"  ); // ee
 
  Int_t    xbin =  hist_fl[fl_mode_ll]->FindBin(c);
  Int_t    ybin = (hist_fl[fl_mode_ll]->FindBin(c, p))/(hist_fl[fl_mode_ll]->GetNbinsX()+2);
  fake  = hist_fl [fl_mode_ll]->GetBinContent( xbin, ybin );
  fakeE = hist_flE[fl_mode_ll]->GetBinContent( xbin, ybin );
  if( fl_mode_ll==1 && xbin==8 && ybin==1 ){
    fake  = (hist_fl[fl_mode_ll]->GetBinContent(7,1) + hist_fl[fl_mode_ll]->GetBinContent(9,1) )/2.0;
    fakeE = sqrt( hist_flE[fl_mode_ll]->GetBinContent(7,1)*hist_flE[fl_mode_ll]->GetBinContent(7,1)
		  +hist_flE[fl_mode_ll]->GetBinContent(9,1)*hist_flE[fl_mode_ll]->GetBinContent(9,1) )/2.0;
  }
  if( fl_message ) std::cout << "     <Get_fake_pi (rm_l=" << fl_mode_ll << ")> "
			     << "( " << c << ", " << p << " ) -> (" << xbin << ", " << ybin << " ) -> "
			     << fake << " +- " << fakeE << std::endl;
}

void Cal_weight_pion( Double_t eff, Double_t effE, Double_t fake, Double_t fakeE, Double_t& weight, Double_t& weightE ){ // Calculate fake/eff
  weight  = (eff==0 ? 0 : fake/eff ); // fake/eff 
  weightE = sqrt( eff*eff*fakeE*fakeE + effE*effE*fake*fake );
}

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

  TChain* chain_B = new TChain( tname_B );
  
  std::stringstream sTmp;
  sTmp << indir << "/gMC_" << type << "_e0" << exp << "*s0" << stream << "*.root";

  Int_t nfile_B = chain_B->Add( sTmp.str().c_str() );
  if( !nfile_B ) std::cerr << "[ABORT] NO file : " << sTmp.str().c_str() << std::endl, abort();

  // --------------------------------------------------------------------  
  TCut cut_lrnb = "1";
  // --------------------------------------------------------------------
  const Char_t* outfile = Form( "%s/gMC_%s_e0%s_s0%d_caseB_single_jpsi_weight.root",outdir,type, exp, stream );
  TFile* rootf  = new TFile( outfile, "RECREATE" );
  TTree* tree_B = new TTree();

  if( chain_B->GetEntries(cut_lrnb) ){
    tree_B = chain_B->CopyTree( cut_lrnb ); 
  }else{
    tree_B = chain_B->CloneTree( 0 );
    Float_t tmp_weight, tmp_weightE;
    tree_B->Branch( "weight",  &tmp_weight,  "weight/F"  );
    tree_B->Branch( "weightE", &tmp_weightE, "weightE/F" );
    tree_B->Write();
    return 0;
  }

  
  // ------------------------------------------------------------------------
  Float_t exprun, event, rm_l, rm_xs;             // Mode
  Float_t epp, emp, mpp, mmp;                     // Momentum
  Float_t pi1p, pi2p, pi3p, pi4p;                 // Momentum
  Float_t lpc, lmc;                               // Direction
  Float_t pi1c, pi2c, pi3c, pi4c;                 // Direction
  Float_t lp_eid, lm_eid, lp_muid, lm_muid;       // Lepton-PID
  Float_t pi1_muid, pi2_muid, pi3_muid, pi4_muid; // Lepton-PID
  Float_t pi1_eid,  pi2_eid,  pi3_eid,  pi4_eid;  // Lepton-PID
  Float_t pi1chg,   pi2chg,   pi3chg,   pi4chg;   // pion charge
  Float_t lpph, lmph, pi1ph, pi2ph, pi3ph, pi4ph; // phi
  Float_t cc_m, cc_morg;
  Float_t dzll3d, dzll3d1, dzll3d2, dzll3d3, dzll3d4;
  Float_t nb_lep1_qq,  nb_lep1_bb,  nb_lep0_qq,  nb_lep0_bb;
  Float_t nb1_lep1_qq, nb1_lep1_bb, nb1_lep0_qq, nb1_lep0_bb;
  Float_t nb2_lep1_qq, nb2_lep1_bb, nb2_lep0_qq, nb2_lep0_bb;
  Float_t nb3_lep1_qq, nb3_lep1_bb, nb3_lep0_qq, nb3_lep0_bb;
  Float_t nb4_lep1_qq, nb4_lep1_bb, nb4_lep0_qq, nb4_lep0_bb;
  Float_t xs_m, xs1_m, xs2_m, xs3_m, xs4_m;
  Float_t lpi1_m, lpi2_m, lpi3_m, lpi4_m;

  Float_t fl_p, fl_c;                       // Momentum and Direction of fake-lepton
  Float_t weight, weightE;                  // Weight ( New Branch )
  
  tree_B->SetBranchAddress( "exprun",   &exprun  );
  tree_B->SetBranchAddress( "event",    &event   );
  tree_B->SetBranchAddress( "rm_l",     &rm_l    );
  tree_B->SetBranchAddress( "rm_xs",    &rm_xs   );
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
  tree_B->SetBranchAddress( "lpph",     &lpph    );  
  tree_B->SetBranchAddress( "lmph",     &lmph    );  
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
  tree_B->SetBranchAddress( "pi1chg",   &pi1chg   );
  tree_B->SetBranchAddress( "pi2chg",   &pi2chg   );
  tree_B->SetBranchAddress( "pi3chg",   &pi3chg   );
  tree_B->SetBranchAddress( "pi4chg",   &pi4chg   );
  tree_B->SetBranchAddress( "pi1ph",    &pi1ph    );
  tree_B->SetBranchAddress( "pi2ph",    &pi2ph    );
  tree_B->SetBranchAddress( "pi3ph",    &pi3ph    );
  tree_B->SetBranchAddress( "pi4ph",    &pi4ph    );
  tree_B->SetBranchAddress( "cc_m",     &cc_m     );
  tree_B->SetBranchAddress( "cc_morg",  &cc_morg  );
  tree_B->SetBranchAddress( "dzll3d",   &dzll3d   );
  tree_B->SetBranchAddress( "dzll3d1",  &dzll3d1  );
  tree_B->SetBranchAddress( "dzll3d2",  &dzll3d2  );
  tree_B->SetBranchAddress( "dzll3d3",  &dzll3d3  );
  tree_B->SetBranchAddress( "dzll3d4",  &dzll3d4  );
  tree_B->SetBranchAddress( "xs_m",     &xs_m     );
  tree_B->SetBranchAddress( "xs1_m",    &xs1_m    );
  tree_B->SetBranchAddress( "xs2_m",    &xs2_m    );
  tree_B->SetBranchAddress( "xs3_m",    &xs3_m    );
  tree_B->SetBranchAddress( "xs4_m",    &xs4_m    );
  tree_B->SetBranchAddress( "lpi1_m",   &lpi1_m   );
  tree_B->SetBranchAddress( "lpi2_m",   &lpi2_m   );
  tree_B->SetBranchAddress( "lpi3_m",   &lpi3_m   );
  tree_B->SetBranchAddress( "lpi4_m",   &lpi4_m   );
  tree_B->SetBranchAddress( "nb_lep1_orgksfw_vtxcl_fmiss1_qq",    &nb_lep1_qq  );
  tree_B->SetBranchAddress( "nb_lep1_orgksfw_vtxcl_fmiss1_bb",    &nb_lep1_bb  );
  tree_B->SetBranchAddress( "nb_lep0_orgksfw_vtxcl_fmiss1_qq",    &nb_lep0_qq  );
  tree_B->SetBranchAddress( "nb_lep0_orgksfw_vtxcl_fmiss1_bb",    &nb_lep0_bb  );
  tree_B->SetBranchAddress( "nb1_lep1_orgksfw_vtxcl_fmiss1_qq",   &nb1_lep1_qq  );
  tree_B->SetBranchAddress( "nb1_lep1_orgksfw_vtxcl_fmiss1_bb",   &nb1_lep1_bb  );
  tree_B->SetBranchAddress( "nb1_lep0_orgksfw_vtxcl_fmiss1_qq",   &nb1_lep0_qq  );
  tree_B->SetBranchAddress( "nb1_lep0_orgksfw_vtxcl_fmiss1_bb",   &nb1_lep0_bb  );
  tree_B->SetBranchAddress( "nb2_lep1_orgksfw_vtxcl_fmiss1_qq",   &nb2_lep1_qq  );
  tree_B->SetBranchAddress( "nb2_lep1_orgksfw_vtxcl_fmiss1_bb",   &nb2_lep1_bb  );
  tree_B->SetBranchAddress( "nb2_lep0_orgksfw_vtxcl_fmiss1_qq",   &nb2_lep0_qq  );
  tree_B->SetBranchAddress( "nb2_lep0_orgksfw_vtxcl_fmiss1_bb",   &nb2_lep0_bb  );
  tree_B->SetBranchAddress( "nb3_lep1_orgksfw_vtxcl_fmiss1_qq",   &nb3_lep1_qq  );
  tree_B->SetBranchAddress( "nb3_lep1_orgksfw_vtxcl_fmiss1_bb",   &nb3_lep1_bb  );
  tree_B->SetBranchAddress( "nb3_lep0_orgksfw_vtxcl_fmiss1_qq",   &nb3_lep0_qq  );
  tree_B->SetBranchAddress( "nb3_lep0_orgksfw_vtxcl_fmiss1_bb",   &nb3_lep0_bb  );
  tree_B->SetBranchAddress( "nb4_lep1_orgksfw_vtxcl_fmiss1_qq",   &nb4_lep1_qq  );
  tree_B->SetBranchAddress( "nb4_lep1_orgksfw_vtxcl_fmiss1_bb",   &nb4_lep1_bb  );
  tree_B->SetBranchAddress( "nb4_lep0_orgksfw_vtxcl_fmiss1_qq",   &nb4_lep0_qq  );
  tree_B->SetBranchAddress( "nb4_lep0_orgksfw_vtxcl_fmiss1_bb",   &nb4_lep0_bb  );

  tree_B->Branch( "weight",  &weight,  "weight/F"  );
  tree_B->Branch( "weightE", &weightE, "weightE/F" );

  TTree* newtree_B = tree_B->CloneTree( 0 );

  Int_t nev = 0;
  while (tree_B->GetEntry(nev, 0)) {
    if( fl_message ) std::cout << "=====================================[" << nev << " rm_l = " << rm_l << ", " << "event = " << event << " ]=========================================" << std::endl
			       << "=====================================[" << rm_xs << " : " << pi1chg << ", " << pi2chg << ", " << pi3chg << ", " << pi4chg << " ]=========================================" << std::endl;
    nev++;

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
    Int_t flag_svd2 = sel_hist ? 1 : 0;
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // lepton cand
    Double_t w_lp    = 0;
    Double_t w_lm    = 0;
    Double_t eff_lp,    eff_lm,    eff_lpE,    eff_lmE;
    Double_t weight_lp, weight_lm, weight_lpE, weight_lmE;

    if( rm_l==1 ){
      Get_eff_lepton( lpc, epp, (Int_t)rm_l, sel_hist, eff_lp, eff_lpE );
      Get_eff_lepton( lmc, emp, (Int_t)rm_l, sel_hist, eff_lm, eff_lmE );
    }else{
      Get_eff_lepton( lpc, mpp, (Int_t)rm_l, sel_hist, eff_lp, eff_lpE );
      Get_eff_lepton( lmc, mmp, (Int_t)rm_l, sel_hist, eff_lm, eff_lmE );
    }
    Cal_weight_lepton( eff_lp, eff_lpE, weight_lp, weight_lpE );
    Cal_weight_lepton( eff_lm, eff_lmE, weight_lm, weight_lmE );
    if( fl_message ) std::cout << "[l+] eff = " << eff_lp << " +- " << eff_lpE << " -> weight = " << weight_lp << " +- " << weight_lpE << std::endl
			       << "[l-] eff = " << eff_lm << " +- " << eff_lmE << " -> weight = " << weight_lm << " +- " << weight_lmE << std::endl;
    
    // pi1
    Double_t weight_pi1=0;
    Double_t weight_pi1E=0;
    Double_t fake_pi1=0;
    Double_t fake_pi1E=0;
    Double_t eff_pi1=0;
    Double_t eff_pi1E=0;
    if( pi1chg!=0 && dzll3d1<0.0190
	&& (
	    (xs1_m<1.1 && rm_l==1 && nb1_lep1_qq>0.91 && nb1_lep1_bb>0.39) ||
	    (xs1_m<1.1 && rm_l==0 && nb1_lep0_qq>0.86 && nb1_lep0_bb>0.56) ||
	    (xs1_m>1.1 && rm_l==1 && nb1_lep1_qq>0.93 && nb1_lep1_bb>0.91) ||
	    (xs1_m>1.1 && rm_l==0 && nb1_lep0_qq>0.92 && nb1_lep0_bb>0.87)
	    )
	&& !(rm_l==1 && PDGmass::jpsi-0.40 < lpi1_m && lpi1_m < PDGmass::jpsi+0.15)
	&& !(rm_l==0 && PDGmass::jpsi-0.25 < lpi1_m && lpi1_m < PDGmass::jpsi+0.10)
	&& lpi1_m > 0.2 && xs1_m < 2.1
	){
      
      Get_fake_pi( pi1c, pi1p, (Int_t)rm_l,fake_pi1, fake_pi1E );
      Get_eff_pi ( pi1c, pi1p, flag_svd2,  eff_pi1,  eff_pi1E  );
      Cal_weight_pion( eff_pi1, eff_pi1E,  fake_pi1, fake_pi1E, weight_pi1, weight_pi1E );
      if     ( rm_l==1 && pi1p<0.4 ) weight_pi1 = 0, weight_pi1E = 0;
      else if( rm_l==0 && pi1p<0.8 ) weight_pi1 = 0, weight_pi1E = 0;
      if( fl_message ) std::cout << "[1st pion] eff = " << eff_pi1    << " +- " << eff_pi1E << ", fake = " << fake_pi1 << " +- " << fake_pi1E
				 << " -> weight = "     << weight_pi1 << " +- " << weight_pi1E << std::endl;
    }
    
    // pi2
    Double_t weight_pi2=0;
    Double_t weight_pi2E=0;
    Double_t fake_pi2=0;
    Double_t fake_pi2E=0;
    Double_t eff_pi2=0;
    Double_t eff_pi2E=0;

    if( pi2chg!=0 && dzll3d2<0.0190
	&& (
	    (xs2_m<1.1 && rm_l==1 && nb2_lep1_qq>0.91 && nb2_lep1_bb>0.39) ||
	    (xs2_m<1.1 && rm_l==0 && nb2_lep0_qq>0.86 && nb2_lep0_bb>0.56) ||
	    (xs2_m>1.1 && rm_l==1 && nb2_lep1_qq>0.93 && nb2_lep1_bb>0.91) ||
	    (xs2_m>1.1 && rm_l==0 && nb2_lep0_qq>0.92 && nb2_lep0_bb>0.87)
	    )
	&& !(rm_l==1 && PDGmass::jpsi-0.40 < lpi2_m && lpi2_m < PDGmass::jpsi+0.15)
	&& !(rm_l==0 && PDGmass::jpsi-0.25 < lpi2_m && lpi2_m < PDGmass::jpsi+0.10)
	&& lpi2_m > 0.2 && xs2_m < 2.1
	){
      
      Get_fake_pi( pi2c, pi2p, (Int_t)rm_l,fake_pi2, fake_pi2E );
      Get_eff_pi ( pi2c, pi2p, flag_svd2,  eff_pi2,  eff_pi2E  );
      Cal_weight_pion( eff_pi2, eff_pi2E,  fake_pi2, fake_pi2E, weight_pi2, weight_pi2E );
      if     ( rm_l==1 && pi2p<0.4 ) weight_pi2 = 0, weight_pi2E = 0;
      else if( rm_l==0 && pi2p<0.8 ) weight_pi2 = 0, weight_pi2E = 0;
      if( fl_message ) std::cout << "[2nd pion] eff = " << eff_pi2    << " +- " << eff_pi2E << ", fake = " << fake_pi2 << " +- " << fake_pi2E
				 << " -> weight = "     << weight_pi2 << " +- " << weight_pi2E << std::endl;
    }
    
    // pi3
    Double_t weight_pi3=0;
    Double_t weight_pi3E=0;
    Double_t fake_pi3=0;
    Double_t fake_pi3E=0;
    Double_t eff_pi3=0;
    Double_t eff_pi3E=0;
    
    if( pi3chg!=0 && dzll3d3<0.0190
	&& (
	    (xs3_m<1.1 && rm_l==1 && nb3_lep1_qq>0.91 && nb3_lep1_bb>0.39) ||
	    (xs3_m<1.1 && rm_l==0 && nb3_lep0_qq>0.86 && nb3_lep0_bb>0.56) ||
	    (xs3_m>1.1 && rm_l==1 && nb3_lep1_qq>0.93 && nb3_lep1_bb>0.91) ||
	    (xs3_m>1.1 && rm_l==0 && nb3_lep0_qq>0.92 && nb3_lep0_bb>0.87)
	    )
	&& !(rm_l==1 && PDGmass::jpsi-0.40 < lpi3_m && lpi3_m < PDGmass::jpsi+0.15)
	&& !(rm_l==0 && PDGmass::jpsi-0.25 < lpi3_m && lpi3_m < PDGmass::jpsi+0.10)
	&& lpi3_m > 0.2 && xs3_m < 2.1
	){
      
      Get_fake_pi( pi3c, pi3p, (Int_t)rm_l,fake_pi3, fake_pi3E );
      Get_eff_pi ( pi3c, pi3p, flag_svd2,  eff_pi3,  eff_pi3E  );
      Cal_weight_pion( eff_pi3, eff_pi3E,  fake_pi3, fake_pi3E, weight_pi3, weight_pi3E );
      if     ( rm_l==1 && pi3p<0.4 ) weight_pi3 = 0, weight_pi3E = 0;
      else if( rm_l==0 && pi3p<0.8 ) weight_pi3 = 0, weight_pi3E = 0;
      if( fl_message ) std::cout << "[3rd pion] eff = " << eff_pi3    << " +- " << eff_pi3E << ", fake = " << fake_pi3 << " +- " << fake_pi3E
				 << " -> weight = "     << weight_pi3 << " +- " << weight_pi3E << std::endl;
    }

    // pi4
    Double_t weight_pi4=0;
    Double_t weight_pi4E=0;
    Double_t fake_pi4=0;
    Double_t fake_pi4E=0;
    Double_t eff_pi4=0;
    Double_t eff_pi4E=0;
    
    if( pi4chg!=0 && dzll3d4<0.0190
	&& (
	    (xs4_m<1.1 && rm_l==1 && nb4_lep1_qq>0.91 && nb4_lep1_bb>0.39) ||
	    (xs4_m<1.1 && rm_l==0 && nb4_lep0_qq>0.86 && nb4_lep0_bb>0.56) ||
	    (xs4_m>1.1 && rm_l==1 && nb4_lep1_qq>0.93 && nb4_lep1_bb>0.91) ||
	    (xs4_m>1.1 && rm_l==0 && nb4_lep0_qq>0.92 && nb4_lep0_bb>0.87)
	    )
	&& !(rm_l==1 && PDGmass::jpsi-0.40 < lpi4_m && lpi4_m < PDGmass::jpsi+0.15)
	&& !(rm_l==0 && PDGmass::jpsi-0.25 < lpi4_m && lpi4_m < PDGmass::jpsi+0.10)
	&& lpi4_m > 0.2 && xs4_m < 2.1
	){
      
      Get_fake_pi( pi4c, pi4p, (Int_t)rm_l,fake_pi4, fake_pi4E );
      Get_eff_pi ( pi4c, pi4p, flag_svd2,  eff_pi4,  eff_pi4E  );
      Cal_weight_pion( eff_pi4, eff_pi4E,  fake_pi4, fake_pi4E, weight_pi4, weight_pi4E );
      if     ( rm_l==1 && pi4p<0.4 ) weight_pi4 = 0, weight_pi4E = 0;
      else if( rm_l==0 && pi4p<0.8 ) weight_pi4 = 0, weight_pi4E = 0;
      if( fl_message ) std::cout << "[4th pion] eff = " << eff_pi4    << " +- " << eff_pi4E << ", fake = " << fake_pi4 << " +- " << fake_pi4E
				 << " -> weight = "     << weight_pi4 << " +- " << weight_pi4E << std::endl;
    }

    weight  = 0;
    weightE = 0;

    if( pi1chg== 1 ){
      weight  += weight_lp * weight_pi1;
      weightE += weight_lpE*weight_lpE*weight_pi1*weight_pi1 + weight_lp*weight_lp*weight_pi1E*weight_pi1E;
    }else if( pi1chg==-1 ){
      weight  += weight_lm * weight_pi1;
      weightE += weight_lmE*weight_lmE*weight_pi1*weight_pi1 + weight_lm*weight_lm*weight_pi1E*weight_pi1E;
    }

    if( pi2chg== 1 ){
      weight  += weight_lp * weight_pi2;
      weightE += weight_lpE*weight_lpE*weight_pi2*weight_pi2 + weight_lp*weight_lp*weight_pi2E*weight_pi2E;
    }else if( pi2chg==-1 ){
      weight  += weight_lm * weight_pi2;
      weightE += weight_lmE*weight_lmE*weight_pi2*weight_pi2 + weight_lm*weight_lm*weight_pi2E*weight_pi2E;
    }

    if( pi3chg== 1 ){
      weight  += weight_lp * weight_pi3;
      weightE += weight_lpE*weight_lpE*weight_pi3*weight_pi3 + weight_lp*weight_lp*weight_pi3E*weight_pi3E;
    }else if( pi3chg==-1 ){
      weight  += weight_lm * weight_pi3;
      weightE += weight_lmE*weight_lmE*weight_pi3*weight_pi3 + weight_lm*weight_lm*weight_pi3E*weight_pi3E;
    }
    if( pi4chg== 1 ){
      weight  += weight_lp * weight_pi4;
      weightE += weight_lpE*weight_lpE*weight_pi4*weight_pi4 + weight_lp*weight_lp*weight_pi4E*weight_pi4E;
    }else if( pi4chg==-1 ){
      weight  += weight_lm * weight_pi4;
      weightE += weight_lmE*weight_lmE*weight_pi4*weight_pi4 + weight_lm*weight_lm*weight_pi4E*weight_pi4E;
    }

    weightE = sqrt(weightE);
    if( fl_message ) std::cout << "WEIGHT = " << weight << " +- " << weightE << std::endl;

    if( weight!=0 ){

      if     ( weight_pi1 !=0 ) cc_m = lpi1_m;
      else if( weight_pi2 !=0 ) cc_m = lpi2_m;
      else if( weight_pi3 !=0 ) cc_m = lpi3_m;
      else if( weight_pi4 !=0 ) cc_m = lpi4_m;
      newtree_B->Fill();
    }
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
