#ifndef MAKECUT_H
#define MAKECUT_H

#include <iostream>
#include <sstream>
#include <string>
#include <TROOT.h>
#include "const.h"

const Int_t    q2_bin = 9;
///* default setting
const Double_t q2_bins[2][q2_bin+1] = { // [mm,ee]
  {0.0, 2.0, 4.3,
   (PDGmass::jpsi -0.25-0.0*0.15)*(PDGmass::jpsi -0.25-0.0*0.15),
   (PDGmass::jpsi +0.10+0.0*0.05)*(PDGmass::jpsi +0.10+0.0*0.05),
   (PDGmass::psi2s-0.15-0.0*0.10)*(PDGmass::psi2s-0.15-0.0*0.10),
   (PDGmass::psi2s+0.10                )*(PDGmass::psi2s+0.10                ),
   16.0, 19.0, 25.0},
  {0.0, 2.0, 4.3,
   (PDGmass::jpsi -0.25-1.0*0.15)*(PDGmass::jpsi -0.25-1.0*0.15),
   (PDGmass::jpsi +0.10+1.0*0.05)*(PDGmass::jpsi +0.10+1.0*0.05),
   (PDGmass::psi2s-0.15-1.0*0.10)*(PDGmass::psi2s-0.15-1.0*0.10),
   (PDGmass::psi2s+0.10                )*(PDGmass::psi2s+0.10                ),
   16.0, 19.0, 25.0},
};
//*/
/* temporal setting for q2(1-6 GeV^2)
const Double_t q2_bins[2][q2_bin+1] = { // [mm,ee]
  {0.0, 0.5, 1.0,
   6.0,
   (PDGmass::jpsi +0.10+0.0*0.05)*(PDGmass::jpsi +0.10+0.0*0.05),
   (PDGmass::psi2s-0.15-0.0*0.10)*(PDGmass::psi2s-0.15-0.0*0.10),
   (PDGmass::psi2s+0.10                )*(PDGmass::psi2s+0.10                ),
   16.0, 19.0, 25.0},
  {0.0, 0.5, 1.0,
   6.0,
   (PDGmass::jpsi +0.10+1.0*0.05)*(PDGmass::jpsi +0.10+1.0*0.05),
   (PDGmass::psi2s-0.15-1.0*0.10)*(PDGmass::psi2s-0.15-1.0*0.10),
   (PDGmass::psi2s+0.10                )*(PDGmass::psi2s+0.10                ),
   16.0, 19.0, 25.0},
};
*/
std::string makeCut_Xs( Int_t fl_xsid=0 );
// fl_xsid : 0(all), 1(K), 2(K*), 3(Xs)


std::string makeCut_mode_category( Int_t fl_mode_xs, Int_t fl_ctgry=-10 );

std::string makeCut_category( Int_t fl_ctgry );
// fl_ctgry  : 0(false), 1(true), 2(pi0-conv.), 3(off-diagonal), 4(other-mode), -10(total)

std::string makeCut_mode_q2fl( Int_t fl_mode_xs, Int_t fl_q2=-10, Int_t fl_flavor=-10, Int_t fl_inv=0 );
///// * [fl_q2, fl_flavor, fl_comb] are vaild only for [false, off-diagonal, other-mode].
std::string makeCut_q2fl( Int_t fl_q2=-10, Int_t fl_flavor=-10, Int_t fl_inv=0 );
///// fl_q2     : 0(false), 1(true),             -10(no cut),
///// fl_flavor : 0(false), 1(true), 2(unknown), -10(no cut),
///// fl_inv   :  1-> !(cut), 0 -> (cut)

std::string makeCut_pi0( Int_t fl_mode_pi0, Int_t fl_self=0 );
// fl_mode_pi0 : 1->modes including pi0(10modes), 0->modes not including pi0(8modes)

std::string makeCut_track( Int_t fl_mode_track, Int_t fl_self=0 );
// fl_mode_track : # of charged tracks from Xs (0~5)

std::string makeCut_body( Int_t fl_mode_track, Int_t fl_self=0 );
// fl_mode_body : # of particles from Xs (1~5)

std::string makeCut_5body_veto();
// veto K4pi modes

std::string makeCut_unflavor_veto();
// veto Ks(10), Kspi0(1010), Kspi+pi-(210), Kspi+pi-pi0(1210), Kspi+pi+pi-pi-(410)


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
std::string makeCut_LRNB_internal       (                                                    Int_t fl_q2, Int_t fl_xs, Int_t fl_mode_ll, Double_t q2=6, Double_t xs=1.1 );
std::string makeCut_LRNB_1d_internal    ( Char_t*  brname, Double_t var,                     Int_t fl_q2, Int_t fl_xs, Int_t fl_mode_ll, Double_t q2=6, Double_t xs=1.1 ); // brname must include %s for bcs
std::string makeCut_LRNB_1d_lep_internal( Char_t*  brname, Double_t var,                     Int_t fl_q2, Int_t fl_xs, Int_t fl_mode_ll, Double_t q2=6, Double_t xs=1.1 ); // brname must include %d for lep and %s for bcs
std::string makeCut_LRNB_2d_internal    ( Char_t*  brname, Double_t var_qq, Double_t var_bb, Int_t fl_q2, Int_t fl_xs, Int_t fl_mode_ll, Double_t q2=6, Double_t xs=1.1 ); // brname must include %s for bcs
std::string makeCut_LRNB_2d_lep_internal( Char_t*  brname, Double_t var_qq, Double_t var_bb, Int_t fl_q2, Int_t fl_xs, Int_t fl_mode_ll, Double_t q2=6, Double_t xs=1.1 ); // brname must include %d for lep and %s for bcs
// fl_q2      : +(high-q2), -(low-q2),     0(no-cut)
// fl_xs      : +(high-xs), -(low-xs),     0(no-cut)
// fl_mode_ll : 1(e),       0(mu),     other(no-cut)
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

std::string makeCut_LRNB_1d    ( Char_t* brname, Int_t fl_q2, Double_t var_lowxs,                                 Double_t var_highxs,                                                                                                                                            Double_t q2=6, Double_t xs=1.1 );
std::string makeCut_LRNB_1d_lep( Char_t* brname, Int_t fl_q2, Double_t var_lowxs_ee,                              Double_t var_highxs_ee,                               Double_t var_lowxs_mm,                              Double_t var_highxs_mm,                               Double_t q2=6, Double_t xs=1.1 );
std::string makeCut_LRNB_2d    ( Char_t* brname, Int_t fl_q2, Double_t var_lowxs_qq,    Double_t var_lowxs_bb,    Double_t var_highxs_qq,    Double_t var_highxs_bb,                                                                                                              Double_t q2=6, Double_t xs=1.1 );
std::string makeCut_LRNB_2d_lep( Char_t* brname, Int_t fl_q2, Double_t var_lowxs_qq_ee, Double_t var_lowxs_bb_ee, Double_t var_highxs_qq_ee, Double_t var_highxs_bb_ee, Double_t var_lowxs_qq_mm, Double_t var_lowxs_bb_mm, Double_t var_highxs_qq_mm, Double_t var_highxs_bb_mm, Double_t q2=6, Double_t xs=1.1 );

std::string makeCut_q2( Int_t fl_mode_ll, Int_t fl_q2, Int_t fl_cos=0 );
// fl_q2[1,2,3,5,7,8,9]
//       4 (  J/psi  veto window )
//       6 ( psi(2s) veto window )
// fl_cos[0,+,-]
//       0( no cut)
//       +( positive cos region)
//       -( negative cos region)

Int_t tag_q2_cos_region(Int_t fl_mode_ll, Double_t q2, Double_t fl_cos);
//  1 ; 1st q2 bin, +cos
//  2 ; 1st q2 bin, -cos
//  3 ; 2nd q2 bin, +cos
//  4 ; 2nd q2 bin, -cos
//  5 ; 3rd q2 bin, +cos
//  6 ; 3rd q2 bin, -cos
//  7 ; 5th q2 bin, +cos
//  8 ; 5th q2 bin, -cos
//  9 ; 7th q2 bin, +cos
// 10 ; 7th q2 bin, -cos
// 11 ; 8th q2 bin, +cos
// 12 ; 8th q2 bin, -cos
// 13 ; 9th q2 bin, +cos
// 14 ; 9th q2 bin, -cos

// -1 : invalid region ( j/psi  veto)
// -2 : invalid region (psi(2s) veto)
// -3 : invalid region (underflow-q2)
// -4 : invalid region ( overflow-q2)

Int_t tag_q2_cos_region2(Int_t fl_mode_ll, Double_t q2, Double_t fl_cos);
//  1 ;    1st & 2nd    q2 bin, +cos
//  2 ;    1st & 2nd    q2 bin, -cos
//  3 ;       3rd       q2 bin, +cos
//  4 ;       3rd       q2 bin, -cos
//  5 ;       5th       q2 bin, +cos
//  6 ;       5th       q2 bin, -cos
//  7 ; 7th & 8th & 9th q2 bin, +cos
//  8 ; 7th & 8th & 9th q2 bin, -cos

Int_t tag_q2_cos_region_totFB(Int_t fl_mode_ll, Double_t q2, Double_t fl_cos);
//  1 ;    q2 bin, +cos
//  2 ;    q2 bin, -cos
#endif
