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



Int_t main( Int_t argc, Char_t** argv ){
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==3 || argc==4) ) std::cerr << "wrong input" << std::endl, abort();
					
  //const Int_t  fl_mode_ll = 1;
  //const Char_t*  axis     = "cc_m";
  Int_t    fl_mode_ll = atoi(argv[1]);
  Int_t    fl_cut     = atoi(argv[2]);
  Int_t    fl_appRun  = 1;
  if( argc==4 ) fl_appRun = atoi( argv[3] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t  Nchain     = 2;
  const Int_t  fl_norm    = 1; // 1(scale), 0(no)
  const Bool_t fl_scan    = true;
  const Double_t scale_val[Nchain] = {6, 100};
  Char_t** infile = new Char_t*[Nchain];
  //infile[0] = "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut5/bkg_522/gMC_mixed_e0*_s00_";
  //infile[1] = "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut5/cc_522/CC_*psi2s_";
  infile[0] = "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut5/bkg_522/gMC_";
  infile[1] = "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut5/cc_522/CC_";
  //infile[0] = "~/ewp/ana/ana_xsll_comb/tmp_gmc/gMC_";
  //infile[1] = "~/ewp/ana/ana_xsll_comb/tmp_cc_comb/CC_";

  Char_t* add_cut_delta = "(cc_morg!=cc_m)"; // applied for (cc_m - cc_morg)  
  Char_t** add_cut = new Char_t*[Nchain];

  //( pi0org==1 || (pi0org<0.5&&rm_xs>999) || (pi0org<0.5&&rm_xs<999))
  // pi0org : 1(true), 0(false&no-pi0 mode), -1(virtual)

  if( fl_cut==0 ){
    add_cut[0] = "(lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3";
    add_cut[1] = "(lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3";
  }else if( fl_cut==1 ){
    add_cut[0] = "rm_xs<1000 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3";
    add_cut[1] = "rm_xs<1000 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3";
  }else if( fl_cut==2 ){
    add_cut[0] = "rm_xs>999 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3";
    add_cut[1] = "rm_xs>999 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3";
    
  }else if( fl_cut==3 ){ // same with cut0
    add_cut[0] = "(lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && ( pi0org==1 || (pi0org<0.5&&rm_xs>999) || (pi0org<0.5&&rm_xs<999))";
    add_cut[1] = "(lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && ( pi0org==1 || (pi0org<0.5&&rm_xs>999) || (pi0org<0.5&&rm_xs<999))";
  }else if( fl_cut==4 ){ // same with cut1
    add_cut[0] = "rm_xs<1000 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && ( pi0org==1 || (pi0org<0.5&&rm_xs>999) || (pi0org<0.5&&rm_xs<999))";
    add_cut[1] = "rm_xs<1000 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && ( pi0org==1 || (pi0org<0.5&&rm_xs>999) || (pi0org<0.5&&rm_xs<999))";
  }else if( fl_cut==5 ){ // same with cut2
    add_cut[0] = "rm_xs>999 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3  && ( pi0org==1 || (pi0org<0.5&&rm_xs>999) || (pi0org<0.5&&rm_xs<999))";
    add_cut[1] = "rm_xs>999 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3  && ( pi0org==1 || (pi0org<0.5&&rm_xs>999) || (pi0org<0.5&&rm_xs<999))";
      
  }else if( fl_cut==6 ){
    add_cut[0] = "(lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==1";
    add_cut[1] = "(lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==1";
  }else if( fl_cut==7 ){
    add_cut[0] = "rm_xs<1000 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==1";
    add_cut[1] = "rm_xs<1000 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==1";
  }else if( fl_cut==8 ){
    add_cut[0] = "rm_xs>999 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==1";
    add_cut[1] = "rm_xs>999 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==1";

  }else if( fl_cut==9 ){
    add_cut[0] = "(lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==0";
    add_cut[1] = "(lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==0";
  }else if( fl_cut==10 ){
    add_cut[0] = "rm_xs<1000 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==0";
    add_cut[1] = "rm_xs<1000 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==0";
  }else if( fl_cut==11 ){
    add_cut[0] = "rm_xs>999 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==0";
    add_cut[1] = "rm_xs>999 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==0";

  }else if( fl_cut==12 ){
    add_cut[0] = "(lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==-1";
    add_cut[1] = "(lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==-1";
  }else if( fl_cut==13 ){
    add_cut[0] = "rm_xs<1000 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==-1";
    add_cut[1] = "rm_xs<1000 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==-1";
  }else if( fl_cut==14 ){
    add_cut[0] = "rm_xs>999 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==-1";
    add_cut[1] = "rm_xs>999 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==-1";
    
  }else if( fl_cut==15 ){
    add_cut[0] = "rm_xs==1 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3";
    add_cut[1] = "rm_xs==1 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3";
  }else if( fl_cut==16 ){
    add_cut[0] = "rm_xs==10 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3";
    add_cut[1] = "rm_xs==10 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3";
  }else if( fl_cut==17 ){
    add_cut[0] = "rm_xs==1001 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3";
    add_cut[1] = "rm_xs==1001 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3";
  }else if( fl_cut==18 ){
    add_cut[0] = "rm_xs==1010 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3";
    add_cut[1] = "rm_xs==1010 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3";

  }else if( fl_cut==19 ){
    add_cut[0] = "rm_xs==1 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==1";
    add_cut[1] = "rm_xs==1 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==1";
  }else if( fl_cut==20 ){
    add_cut[0] = "rm_xs==10 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==1";
    add_cut[1] = "rm_xs==10 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==1";
  }else if( fl_cut==21 ){
    add_cut[0] = "rm_xs==1001 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==1";
    add_cut[1] = "rm_xs==1001 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==1";
  }else if( fl_cut==22 ){
    add_cut[0] = "rm_xs==1010 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==1";
    add_cut[1] = "rm_xs==1010 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==1";

  }else if( fl_cut==23 ){
    add_cut[0] = "rm_xs==1 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==0";
    add_cut[1] = "rm_xs==1 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==0";
  }else if( fl_cut==24 ){
    add_cut[0] = "rm_xs==10 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==0";
    add_cut[1] = "rm_xs==10 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==0";
  }else if( fl_cut==25 ){
    add_cut[0] = "rm_xs==1001 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==0";
    add_cut[1] = "rm_xs==1001 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==0";
  }else if( fl_cut==26 ){
    add_cut[0] = "rm_xs==1010 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==0";
    add_cut[1] = "rm_xs==1010 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==0";

  }else if( fl_cut==27 ){
    add_cut[0] = "rm_xs==1 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==-1";
    add_cut[1] = "rm_xs==1 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==-1";
  }else if( fl_cut==28 ){
    add_cut[0] = "rm_xs==10 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==-1";
    add_cut[1] = "rm_xs==10 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==-1";
  }else if( fl_cut==29 ){
    add_cut[0] = "rm_xs==1001 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==-1";
    add_cut[1] = "rm_xs==1001 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==-1";
  }else if( fl_cut==30 ){
    add_cut[0] = "rm_xs==1010 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==-1";
    add_cut[1] = "rm_xs==1010 && (lpgt==3 && lmgt==3 ) && (lpself==1 && lmself==1) && korg==lporg && lporg==lmorg && abs(pi1org)==lporg && abs(pi2org)==lporg && abs(pi3org)==lporg && abs(pi4org)==lporg && rest_sw==0 && gm_bg1<3 && pi0org==-1";
  };
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  const Char_t*  tname    = "h511";

  const Int_t Naxis = 3;  
  const Double_t offset  [Naxis] = { 0.0, 0.0,  0.0 };
  const Int_t    xbin    [Naxis] = { 100, 100,  100 }; 
  const Double_t xmin    [Naxis] = { 2.0, 2.0, -0.2 };
  const Double_t xmax    [Naxis] = { 4.0, 4.0,  2.3 };
  const Char_t*  axis    [Naxis] = { "cc_m",         "cc_morg",            "cc_m-cc_morg"        };
  const Char_t*  fname   [Naxis] = { "cc_m",         "cc_morg",            "cc_m-cc_morg"        };
  const Char_t*  xlabel  [Naxis] = { "M_{ll} [GeV]", "M_{ll}^{org} [GeV]", "#Delta M_{ll} [GeV]" };
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MChain**  chain    = new MChain*[Nchain];
  TH1D***    tmphist = new TH1D** [Naxis];
  for( Int_t i=0; i<Naxis; i++ ) tmphist[i] = new TH1D*[Nchain];
  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++
  std::cout << Form(" ************************ chain-tree ( %s ) *************************************",tname) << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<infile %d > ", j );
    chain[j] = new MChain( infile[j], tname, branch_table(), 0, "*.root" );
    nominal_cut_selection( chain[j], fl_mode_ll )( chain[j]->GetCut(), tname );
  }

  // ++++++++++++++++++++++++
  // cut change
  for( Int_t j=0; j<Nchain; j++ ){
    chain[j]->GetCut()->Set( "Mbc", 1,  5.22, 0.0, 5.30 );
    //chain[j]->GetCut()->Set( "de",  1, -0.20, 0.0, 0.20 );
    //chain[j]->GetCut()->Set(    441, 0 );
    //chain[j]->GetCut()->Set(    443, 0 );
    //chain[j]->GetCut()->Set( 100443, 0 );
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
  for( Int_t i=0; i<Naxis; i++ ){
    for( Int_t j=0; j<Nchain; j++ ){
      tmphist[i][j] = new TH1D( Form("tmphist%d_%d",i,j), Form("%s",chain[j]->GetChange()), xbin[i],offset[i]+xmin[i],offset[i]+xmax[i] );
      if( i==2 ) chain[j]->GetTree()->Project( Form("tmphist%d_%d",i,j), axis[i], Form( "%s && %s", add_cut[j], add_cut_delta) );
      else       chain[j]->GetTree()->Project( Form("tmphist%d_%d",i,j), axis[i], add_cut[j] );
      std::cout << Form("add_cut %d : ", j) << add_cut[j] << std::endl;
      Deco( tmphist[i][j], 1, j+1, j+1 );
    }
  }
  

  // +++++++ scan for axis[0] ++++++++++++++++++++++++++++++++++
  if( fl_scan ){
    for( Int_t j=0; j<Nchain; j++ ){
      std::cout << "+++++++++++++++++++++++++++++++++++++++++++++< SCAN tmphist" << j << " : "
		<< chain[j]->GetTree()->GetEntries( add_cut[j] )
		<< " entries >++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		<< std::endl;
      chain[j]->GetTree()->SetScanField(0);
      chain[j]->GetTree()->Scan(
				"rm_xs:cc_m:cc_morg:cc_m-cc_morg:Mbc:de:lpgt:lmgt:lporg:lmorg:lpself:lmself:lpselfid:lmselfid:lpmoid:lmmoid:kmoid:korg:pi0self:gb1_semi:gb2_semi:gb1nd:gb2nd:gb1d1_se:gb1d2_se:gb2d1_se:gb2d2_se:gm_bg1:gm_bg2:gm_b1:gm_b2:gm_l1:gm_l2:gm_nu1:gm_nu2:rest:rest_sw:rest2:rest2_sw:dntrk:ntrk",
				Form("2.7>%s && %s", axis[0], add_cut[j]) );
      chain[j]->GetTree()->Scan(
				"rm_xs:cc_m:cc_morg:cc_m-cc_morg:Mbc:de:lpgt:lmgt:lporg:lmorg:lpself:lmself:lpselfid:lmselfid:lpmoid:lmmoid:kmoid:korg:pi0self:gb1_semi:gb2_semi:gb1nd:gb2nd:gb1d1_se:gb1d2_se:gb2d1_se:gb2d2_se:gm_bg1:gm_bg2:gm_b1:gm_b2:gm_l1:gm_l2:gm_nu1:gm_nu2:rest:rest_sw:rest2:rest2_sw:dntrk:ntrk",
				Form("3.2<%s && %s", axis[0], add_cut[j]) );
      std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    }
  }

  // +++++++ display ++++++++++++++++++++++++++++++++++
  Double_t entry_all[Naxis][Nchain] = {0};
  for( Int_t i=0; i<Naxis; i++ ){
    for( Int_t j=0; j<Nchain; j++ ){
      std::cout << Form("<tmphist %d %d> ",i, j);
      entry_all[i][j] = tmphist[i][j]->GetEntries();
      Double_t entry_canvas = tmphist[i][j]->Integral();
      Double_t entry_under  = tmphist[i][j]->GetBinContent(0);
      Double_t entry_over   = tmphist[i][j]->GetBinContent(xbin[i]+1);
      std::cout << entry_all[i][j]
		<< " events ( canvas : "
		<< entry_canvas
		<< " / under : "
		<< entry_under
		<< " / over  : "
		<< entry_over;
      if( i!=2 ) std::cout << " / 2.0-2.7  : "
			   << tmphist[i][j]->Integral( 1, tmphist[i][j]->FindBin(2.7-0.0000001) );
      std::cout << "]"
		<< std::endl;
      if( fl_norm==1 ){
	tmphist[i][j]->Sumw2();
	tmphist[i][j]->Scale( 1/scale_val[j] );
      }
    }
  }
  // +++++++ draw ++++++++++++++++++++++++++++++++++
  TCanvas**  c1   = new TCanvas*[Naxis];
  TH2D***   waku  = new TH2D**[Naxis];
  TH2D***   wakun = new TH2D**[Naxis];

  for( Int_t i=0; i<Naxis; i++ ){
    c1[i] = Canvas( Form("c%d",i),Form("c%d"), 4, 3 );
    c1[i]->Draw();
    waku [i] = new TH2D*[2];
    wakun[i] = new TH2D*[2];
    if( i==2 ) waku [i][0] = Waku( Nchain, tmphist[i], xlabel[i] );
    else       waku [i][0] = new TH2D( Form("w_%d", i), xlabel[i], 2, 2.7, 3.2, 2, 0.0, 1.3*tmphist[i][1]->GetMaximum() );

    Double_t tmp_ymax = 0;
    if( i==2 ) tmp_ymax = tmphist[i][1]->GetMaximum()/4.0;
    else       tmp_ymax = 10.0*tmphist[i][1]->Integral(1,40)/40.0;
    if( tmp_ymax == 0 ) tmp_ymax = 10;
    if( i==2 ){
      waku [i][1] = new TH2D( Form("waku_%d", i), Form( "%s (  Scaled  )", xlabel[i]), 2, xmin[i], xmax[i], 2, 0.0, tmp_ymax );
      wakun[i][1] = new TH2D( Form("wakun_%d",i), Form( "%s (Normalized)", xlabel[i]), 2, xmin[i], xmax[i], 2, 0.0, 0.10     );
    }else if( i==1 ){
      waku [i][1] = new TH2D( Form("waku_%d", i), Form( "%s (  Scaled  )", xlabel[i]), 2, 2.0, 3.2, 2, 0.0, tmp_ymax );
      wakun[i][1] = new TH2D( Form("wakun_%d",i), Form( "%s (Normalized)", xlabel[i]), 2, 2.0, 3.2, 2, 0.0, 0.10     );
    }else if( i==0 ){
      waku [i][1] = new TH2D( Form("waku_%d", i), Form( "%s (  Scaled  )", xlabel[i]), 2, 2.0, 3.2, 2, 0.0, tmp_ymax );
      wakun[i][1] = new TH2D( Form("wakun_%d",i), Form( "%s (Normalized)", xlabel[i]), 2, 2.0, 3.2, 2, 0.0, 0.0005   );
    }
    
    c1[i]->cd(1);
    waku[i][0]->Draw();
    for(Int_t j=0; j<Nchain; j++ ) tmphist[i][j]->DrawCopy( "hist same" );

    c1[i]->cd(2);
    waku[i][1]->Draw();
    for(Int_t j=0; j<Nchain; j++ ) tmphist[i][j]->DrawCopy( "hist same" );

    c1[i]->cd(3);
    for(Int_t j=0; j<Nchain; j++ ) tmphist[i][j]->DrawNormalized( j==0 ? "" : "hist same" );
    
    c1[i]->cd(4);
    wakun[i][1]->Draw();
    for(Int_t j=0; j<Nchain; j++ ) tmphist[i][j]->DrawNormalized( "hist same" );

    c1[i]->cd(5);
    waku[i][0]->Draw();
    for(Int_t j=0; j<Nchain; j++ ) tmphist[i][j]->DrawCopy( "same" );

    c1[i]->cd(6);
    waku[i][1]->Draw();
    for(Int_t j=0; j<Nchain; j++ ) tmphist[i][j]->DrawCopy( "same" );

    c1[i]->cd(7);
    for(Int_t j=0; j<Nchain; j++ ) tmphist[i][j]->DrawNormalized( j==0 ? "" : "same" );

    c1[i]->cd(8);
    wakun[i][1]->Draw();
    for(Int_t j=0; j<Nchain; j++ ) tmphist[i][j]->DrawNormalized( "same" );

    c1[i]->cd(9);
    gPad->SetLogy();
    for(Int_t j=0; j<Nchain; j++ ) tmphist[i][j]->DrawCopy( j==0 ? "hist" : "hist same" );

    c1[i]->cd(10);
    gPad->SetLogy();
    for(Int_t j=0; j<Nchain; j++ ) tmphist[i][j]->DrawCopy( j==0 ? "" : "same" );

    c1[i]->cd(11);
    gPad->SetLogy();
    for(Int_t j=0; j<Nchain; j++ ) tmphist[i][j]->DrawNormalized( j==0 ? "hist" : "hist same" );

    c1[i]->cd(12);
    gPad->SetLogy();
    for(Int_t j=0; j<Nchain; j++ ){
      tmphist[i][j]->DrawNormalized( j==0 ? "" : "same");
    }

    
    c1[i]->Update();
    c1[i]->Print( Form("pic/M_ll_%s_lep%d_cut%d_c%d.eps",  fname[i], fl_mode_ll,fl_cut,i) );
  }
  //tmphist[1]->DrawNormalized();
  //tmphist[0]->DrawNormalized("same");

  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();

  return 0;
}

