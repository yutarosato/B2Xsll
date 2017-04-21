#include "Branch.h"

MCut_array branch_table(){
  Int_t n = 53;
  MCut_array branch_list(n);
  Int_t i=0;
  //                    <branch-name,tag,type>
  // function-type(1: window, 2:low, 3:high, 4:esc, 5:bool(0or1), 6:bool(-1or+1), 7:window + minus value, (0:optional)
  // detail of function-type is written in "MCut.h"
  // [ Gen ]
  branch_list.Init( i++, "gm_fl_xs*gm_m_xs", 0, 4 );
  branch_list.Init( i++, "gm_l",             0, 5 );
  // [ B ]
  branch_list.Init( i++, "Mbc",           0, 1 );
  branch_list.Init( i++, "de",            0, 1 );
  branch_list.Init( i++, "deorg",         0, 1 );
  branch_list.Init( i++, "decalib",       0, 1 );
  // [ QQ ]
  branch_list.Init( i++, "thrust",        0, 3 );
  branch_list.Init( i++, "R2",            0, 3 );
  // [ Pseudo Charmonium ] <- set it using tag. (not using branch-name)
  branch_list.Init( i++, "cc_m",        441, 2 ); // for veto low region
  branch_list.Init( i++, "cc_m",        443, 4 ); // for veto J/psi
  branch_list.Init( i++, "cc_m",     100443, 4 ); // for veto psi(2s)
  branch_list.Init( i++, "cc_morg",     441, 2 ); // for veto low region
  branch_list.Init( i++, "cc_morg",     443, 4 ); // for veto J/psi
  branch_list.Init( i++, "cc_morg",  100443, 4 ); // for veto psi(2s)
  // [ Xs ]
  branch_list.Init( i++, "xs_m",          0, 3 );
  // [ Lepton ]
  branch_list.Init( i++, "rm_l",          0, 5 );
  branch_list.Init( i++, "lp_eid",       11, 2 ); branch_list.Init_second( i-1, "epp", 2 );
  branch_list.Init( i++, "lm_eid",       11, 2 ); branch_list.Init_second( i-1, "emp", 2 );
  branch_list.Init( i++, "lp_muid",      13, 2 ); branch_list.Init_second( i-1, "mpp", 2 );
  branch_list.Init( i++, "lm_muid",      13, 2 ); branch_list.Init_second( i-1, "mmp", 2 );
  branch_list.Init( i++, "lp_kid",     321, 3 );
  branch_list.Init( i++, "lm_kid",     321, 3 );
  // [ Charged K ]
  branch_list.Init( i++, "k_kid",       321, 2 );
  branch_list.Init( i++, "k_eid",        11, 3 );  branch_list.Init_second( i-1, "k_ep", 3 );
  branch_list.Init( i++, "k_muid",       13, 3 );  branch_list.Init_second( i-1, "k_mp", 3 );
  // [ Ks ]
  branch_list.Init( i++, "ks_m",          0, 7 );
  branch_list.Init( i++, "kspi1id",    -321, 3 );
  branch_list.Init( i++, "kspi2id",    -321, 3 );
  // [ neutral pi ]
  branch_list.Init( i++, "pi0_m",         0, 7 );
  branch_list.Init( i++, "pi0e",          0, 4 );
  branch_list.Init( i++, "pi0gam1e",      0, 4 );
  branch_list.Init( i++, "pi0gam2e",      0, 4 );
  // [ Charged pi ]
  branch_list.Init( i++, "pi1_kid",     321, 3 );
  branch_list.Init( i++, "pi2_kid",     321, 3 );
  branch_list.Init( i++, "pi3_kid",     321, 3 );
  branch_list.Init( i++, "pi4_kid",     321, 3 );
  branch_list.Init( i++, "pi1_eid",      11, 3 );  branch_list.Init_second( i-1, "pi1_ep", 3 );
  branch_list.Init( i++, "pi2_eid",      11, 3 );  branch_list.Init_second( i-1, "pi2_ep", 3 );
  branch_list.Init( i++, "pi3_eid",      11, 3 );  branch_list.Init_second( i-1, "pi3_ep", 3 );
  branch_list.Init( i++, "pi4_eid",      11, 3 );  branch_list.Init_second( i-1, "pi4_ep", 3 );
  branch_list.Init( i++, "pi1_muid",     13, 3 );  branch_list.Init_second( i-1, "pi1_mp", 3 );
  branch_list.Init( i++, "pi2_muid",     13, 3 );  branch_list.Init_second( i-1, "pi2_mp", 3 );
  branch_list.Init( i++, "pi3_muid",     13, 3 );  branch_list.Init_second( i-1, "pi3_mp", 3 );
  branch_list.Init( i++, "pi4_muid",     13, 3 );  branch_list.Init_second( i-1, "pi4_mp", 3 );
  // [ Vertex ]
  branch_list.Init( i++, "dzll3d",        0, 3 );
  branch_list.Init( i++, "dzll3dorg",     0, 3 );
  branch_list.Init( i++, "dzll3dcalib",   0, 3 );
  branch_list.Init( i++, "kfbcl",         0, 2 );  // branch_list.Init( i++, "kfbchi/kfbdgf", 0, 1 );
  branch_list.Init( i++, "kflchi/kfldgf", 0, 1 );
  branch_list.Init( i++, "drmax",         0, 3 );
  branch_list.Init( i++, "dzmax",         0, 3 );


  // [ h12(Gen) ]
  branch_list.Init( i++, "llg_m",         0, 2 );
  branch_list.Init( i++, "gm_fl_xs*Xs_m", 0, 4 );
  
  if( i != n ){
    std::cout << Form("Check # of branches(%d:%d)!!",n,i) << std::flush << std::endl;
    abort();
  }
  return branch_list;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// cut->Set       ( tag,  flag,          val1, offset, val2 );
// cut->Set       ( name, flag,          val1, offset, val2 );
// cut->Set_second( name, flag, combine, val1, offset, val2 );
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void make_cut_common( MCut_array* cut, const Char_t* h ){
  cut->AllOff();
  if( strcmp(h,"h511")==0 ){
    // [ Gen ]
    cut->Set( "gm_fl_xs*gm_m_xs", 1, 0.0, 0.0, 1.10 );
    // [ PID ] 
    cut->Set( 321, 1, 0.6  ); // k-id
    cut->Set(  11, 1, 0.8  ); // e-id
    cut->Set(  13, 1, 0.97 ); // mu-id

    cut->Set_second(   "k_ep", 1, 1, 0.40 );
    cut->Set_second( "pi1_ep", 1, 1, 0.40 );
    cut->Set_second( "pi2_ep", 1, 1, 0.40 );
    cut->Set_second( "pi3_ep", 1, 1, 0.40 );
    cut->Set_second( "pi4_ep", 1, 1, 0.40 );
    cut->Set_second(   "k_mp", 1, 1, 0.80 );
    cut->Set_second( "pi1_mp", 1, 1, 0.80 );
    cut->Set_second( "pi2_mp", 1, 1, 0.80 );
    cut->Set_second( "pi3_mp", 1, 1, 0.80 );
    cut->Set_second( "pi4_mp", 1, 1, 0.80 );
    // [ Momentum and Energy ]
    cut->Set( "pi0e",     1, -0.5, 0.0, 0.40 ); // neutral pi
    cut->Set( "pi0gam1e", 1, -0.5, 0.0, 0.05 ); // gamma from neutral pi
    cut->Set( "pi0gam2e", 1, -0.5, 0.0, 0.05 ); // gamma from neutral p
    // [ Mass ]
    cut->Set( "ks_m",  1, -0.015, PDGmass::ks,  0.015 );
    cut->Set( "pi0_m", 1, -0.010, PDGmass::pi0, 0.010 );
    cut->Set( "xs_m",  1, 2.0                         );
    // [ QQ ]
    make_cut_QQ  ( cut );
    // [ Vertex ]
    cut->Set( "kfbcl",         0 );// 1.0e-18 <- This selection does not work correctly (maybe this cut value is too small for "Double_t"), so this selection is applied in through_cut.cpp(sigmc) and merge_cut.cpp(gmc)
    //cut->Set( "kfbchi/kfbdgf", 1, 0.0, 0.0, 10.0 );
    cut->Set( "dzll3d",        1, 0.019   ); // modified at 20120710 (150[um]->190[um])
    cut->Set( "drmax",         1, 1.0     );
    cut->Set( "dzmax",         1, 5.0     );
  }else if( strcmp(h,"h12")==0 ){
    cut->AllOff();
    cut->Set( "llg_m",         1, 0.20 ); // dilepton invariant mass by generator information
    cut->Set( "gm_fl_xs*Xs_m", 1, 0.00, 0.0, 1.10 );
  }else{
    std::cout << Form("hist-name(%s)",h) << std::endl;
  }

  return;
}

void make_cut_ee( MCut_array* cut, const Char_t* h ){
  make_cut_common( cut, h );

  if( strcmp(h,"h511")==0 ){
    // [ Lepton Mode ]
    cut->Set( "rm_l", 1, 1 );
    // [ Momentum and Energy ]
    cut->Set_second( "epp", 1, 0, 0.40 ); // for e
    cut->Set_second( "emp", 1, 0, 0.40 ); // for e
    // [ B ]
    cut->Set( "Mbc", 1,  5.27, 0.0, 5.29 );
    cut->Set( "de",  1, -0.10, 0.0, 0.05 ); // for e
    // [ Pseudo Charmonium ]
    make_cut_cc_ee ( cut );
  }else if( strcmp(h,"h12")==0 ){
    cut->Set( "gm_l", 1, 1 );
  }
  cut->ClearChange();
  return;
}

void make_cut_mm( MCut_array* cut, const Char_t* h ){
  make_cut_common( cut, h );
  
  if( strcmp(h,"h511")==0 ){
    // [ Lepton Mode ]
    cut->Set( "rm_l", 1, 0 );
    // [ PID ] 
    cut->SetFunc( 11, 3 );  // e-id veto
    // [ Momentum and Energy ]
    cut->SetFunc_second( "epp", 3 );
    cut->SetFunc_second( "emp", 3 );
    cut->Set_second( "epp", 1, 1, 0.40 ); // for mu
    cut->Set_second( "emp", 1, 1, 0.40 ); // for mu
    cut->Set_second( "mpp", 1, 0, 0.80 ); // for mu
    cut->Set_second( "mmp", 1, 0, 0.80 ); // for mu
    // [ B ]
    cut->Set( "Mbc", 1,  5.27, 0.0, 5.29 );
    cut->Set( "de",  1, -0.05, 0.0, 0.05 ); // for mu
    // [ Pseudo Charmonium ]
    make_cut_cc_mm ( cut );
  }else if( strcmp(h,"h12")==0 ){
    cut->Set( "gm_l", 1, 0 );
  }
  cut->ClearChange();
  return;
}

void make_cut_jpsi_ee( MCut_array* cut, const Char_t* h ){
  make_cut_common( cut, h );

  if( strcmp(h,"h511")==0 ){
    cut->Set( "dzll3d",      0);
    cut->Set( "dzll3dcalib", 1, 0.019 );
    // [ Lepton Mode ]
    cut->Set( "rm_l", 1, 1 );
    // [ Momentum and Energy ]
    cut->Set_second( "epp", 1, 0, 0.40 ); // for e
    cut->Set_second( "emp", 1, 0, 0.40 ); // for e
    // [ B ]
    cut->Set( "Mbc",     1,  5.27, 0.0, 5.29 );
    cut->Set( "decalib", 1, -0.10, 0.0, 0.05 ); // for e
    // [ Pseudo Charmonium ]
    cut->Set(    441,    0 );
    cut->SetFunc(443,    1 );
    cut->Set(    443, 1, -0.40, PDGmass::jpsi,  0.10 );
    cut->Set( "cc_morg", 0 );
    cut->Set( 100443,    0 );
  }
  cut->ClearChange();
  return;
}

void make_cut_jpsi_mm( MCut_array* cut, const Char_t* h ){
  make_cut_common( cut, h );
  
  if( strcmp(h,"h511")==0 ){
    cut->Set( "dzll3d",      0);
    cut->Set( "dzll3dcalib", 1, 0.019 );
    // [ Lepton Mode ]
    cut->Set( "rm_l", 1, 0 );
    // [ PID ] 
    cut->SetFunc( 11, 3 );  // e-id veto
    // [ Momentum and Energy ]
    cut->SetFunc_second( "epp", 3 );
    cut->SetFunc_second( "emp", 3 );
    cut->Set_second( "epp", 1, 1, 0.40 ); // for mu
    cut->Set_second( "emp", 1, 1, 0.40 ); // for mu
    cut->Set_second( "mpp", 1, 0, 0.80 ); // for mu
    cut->Set_second( "mmp", 1, 0, 0.80 ); // for mu
    // [ B ]
    cut->Set( "Mbc",     1,  5.27, 0.0, 5.29 );
    cut->Set( "decalib", 1, -0.05, 0.0, 0.05 ); // for mu
    // [ Pseudo Charmonium ]
    cut->Set(    441,    0 );
    cut->SetFunc(443,    1 );
    cut->Set(    443, 1, -0.25, PDGmass::jpsi,  0.10 );
    cut->Set( "cc_morg", 0 );
    cut->Set( 100443,    0 );
  }
  cut->ClearChange();
  return;
}

void make_cut_emu( MCut_array* cut, const Char_t* h ){
  make_cut_common( cut, h );
  
  if( strcmp(h,"h511")==0 ){
    // [ Lepton Mode ]
    cut->Set( "rm_l", 0 ); // flag-off !!
    // [ PID ] 
    cut->SetFunc( "lm_eid", 3 );  // e-id veto // for mu(l-)

    // [ Momentum and Energy ]
    // branches for positive charge tracks are used for electron in emu mode.
    // branches for negative charge tracks are used for muon     in emu mode.
    cut->Set_second( "epp", 1, 0, 0.40 ); // for e(l+)
    
    cut->SetFunc_second( "emp", 3 );      // for mu(l-)
    cut->Set_second( "emp", 1, 1, 0.40 ); // for mu(l-)
    cut->Set_second( "mmp", 1, 0, 0.80 ); // for mu(l-)

    // [ B ] the wider ranges (for ee) in de and cc-veto are selected for emu mode.
    cut->Set( "Mbc", 1,  5.27, 0.0, 5.29 );
    cut->Set( "de",  1, -0.10, 0.0, 0.05 ); // for e
    // [ Pseudo Charmonium ]
    make_cut_cc_ee ( cut );
  }
  cut->ClearChange();
  return;
}

void make_cut_ee_dembc( MCut_array* cut, const Char_t* h ){
  cut->AllOff();
  if( strcmp(h,"h511")==0 ){
    cut->Set( "rm_l", 1, 1 );
    cut->Set( "Mbc", 1,  5.27, 0.0, 5.29 );
    cut->Set( "de",  1, -0.10, 0.0, 0.05 ); // for e
  }else if( strcmp(h,"h12")==0 ){
    cut->Set( "llg_m",         1, 0.20 ); // dilepton invariant mass by generator information
    cut->Set( "gm_fl_xs*Xs_m", 1, 0.00, 0.0, 1.10 );
    cut->Set( "gm_l", 1, 1 );
  }else{
    std::cout << Form("hist-name(%s)",h) << std::endl;
  }
  cut->ClearChange();
  return;
}

void make_cut_mm_dembc( MCut_array* cut, const Char_t* h ){
  cut->AllOff();
  if( strcmp(h,"h511")==0 ){
    cut->Set( "rm_l", 1, 0 );
    cut->Set( "Mbc", 1,  5.27, 0.0, 5.29 );
    cut->Set( "de",  1, -0.05, 0.0, 0.05 ); // for mu
  }else if( strcmp(h,"h12")==0 ){
    cut->Set( "llg_m",         1, 0.20 ); // dilepton invariant mass by generator information
    cut->Set( "gm_fl_xs*Xs_m", 1, 0.00, 0.0, 1.10 );
    cut->Set( "gm_l", 1, 0 );
  }else{
    std::cout << Form("hist-name(%s)",h) << std::endl;
  }
  cut->ClearChange();
  return;
}

void make_cut_jpsi_ee_dembc( MCut_array* cut, const Char_t* h ){
  cut->AllOff();
  if( strcmp(h,"h511")==0 ){
    cut->Set( "rm_l",     1, 1 );
    cut->Set( "Mbc",      1,  5.27, 0.0, 5.29 );
    cut->Set( "decalib",  1, -0.10, 0.0, 0.05 ); // for e
  }else if( strcmp(h,"h12")==0 ){
    cut->Set( "llg_m",         1, 0.20 ); // dilepton invariant mass by generator information
    cut->Set( "gm_fl_xs*Xs_m", 1, 0.00, 0.0, 1.10 );
  }else{
    std::cout << Form("hist-name(%s)",h) << std::endl;
  }
  cut->ClearChange();
  return;
}

void make_cut_jpsi_mm_dembc( MCut_array* cut, const Char_t* h ){
  cut->AllOff();
  if( strcmp(h,"h511")==0 ){
    cut->Set( "rm_l",     1, 0 );
    cut->Set( "Mbc",      1,  5.27, 0.0, 5.29 );
    cut->Set( "decalib",  1, -0.05, 0.0, 0.05 ); // for mu
  }else if( strcmp(h,"h12")==0 ){
    cut->Set( "llg_m",         1, 0.20 ); // dilepton invariant mass by generator information
    cut->Set( "gm_fl_xs*Xs_m", 1, 0.00, 0.0, 1.10 );
  }else{
    std::cout << Form("hist-name(%s)",h) << std::endl;
  }
  cut->ClearChange();
  return;
}

void make_cut_emu_dembc( MCut_array* cut, const Char_t* h ){
  cut->AllOff();
  if( strcmp(h,"h511")==0 ){
    cut->Set( "Mbc", 1,  5.27, 0.0, 5.29 );
    cut->Set( "de",  1, -0.10, 0.0, 0.05 ); // for e (emu)
  }else if( strcmp(h,"h12")==0 ){
    cut->Set( "llg_m",         1, 0.20 ); // dilepton invariant mass by generator information
    cut->Set( "gm_fl_xs*Xs_m", 1, 0.00, 0.0, 1.10 );
  }else{
    std::cout << Form("hist-name(%s)",h) << std::endl;
  }
  cut->ClearChange();
  return;
}

void make_cut_cc_ee( MCut_array* cut ){
  cut->Set(    441, 1,  0.2                        );
  cut->Set(    443, 1, -0.40, PDGmass::jpsi,  0.15 );
  cut->Set( 100443, 1, -0.25, PDGmass::psi2s, 0.10 );
  return; 
}
void make_cut_cc_mm( MCut_array* cut ){
  cut->Set(    441, 1,  0.2                        );
  cut->Set(    443, 1, -0.25, PDGmass::jpsi,  0.10 );
  cut->Set( 100443, 1, -0.15, PDGmass::psi2s, 0.10 );
  return;
}

void make_cut_QQ( MCut_array* cut ){
  //cut->Set( "R2",     1, 0.50 );
  //cut->Set( "thrust", 1, 0.90 );
}
