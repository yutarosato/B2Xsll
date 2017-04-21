#include "Nominal_cut_selection.h"

Nominal_Cut nominal_cut_selection( MChain* chain ){
  const Char_t* basename = gSystem->BaseName( chain->GetFileName() );
  
  chain->SetModeName( "MC(tentative muon)" );
  std::cout << Form("     <Mode> %s (%s)",chain->GetModeName(),basename) << std::endl;
  return make_cut_mm;

  //if( strstr(basename, "RECB1_JD1_Kst0_EMU0") != NULL || strstr(basename, "_JpsiKs_EMU0_") != NULL ){
  //chain->SetModeName( "Jpsi(mu)Ks" );
  //std::cout << Form("     <Mode> %s (%s)",chain->GetModeName(),basename) << std::endl;
  //return make_cut_jpsiks_mu;
  //}
  
}


Nominal_Cut nominal_cut_selection( MChain* chain, int fl_mode ){
  const Char_t* basename = gSystem->BaseName( chain->GetFileName() );
  

  if( fl_mode == 0 ){
    chain->SetModeName( "MC(muon mode(DE-MBC))" );
    std::cout << Form("     <Mode> %s (%s)",chain->GetModeName(),basename) << std::endl;
    return make_cut_mm_dembc;
  }else if( fl_mode == 1 ){
    chain->SetModeName( "MC(electron mode(DE-MBC))" );
    std::cout << Form("     <Mode> %s (%s)",chain->GetModeName(),basename) << std::endl;
    return make_cut_ee_dembc;
  }else if( fl_mode == 2 ){
    chain->SetModeName( "MC(emu mode(DE-MBC))" );
    std::cout << Form("     <Mode> %s (%s)",chain->GetModeName(),basename) << std::endl;
    return make_cut_emu_dembc;
  }else if( fl_mode == 10 ){
    chain->SetModeName( "MC(muon mode)" );
    std::cout << Form("     <Mode> %s (%s)",chain->GetModeName(),basename) << std::endl;
    return make_cut_mm;
  }else if( fl_mode == 11 ){
    chain->SetModeName( "MC(electron mode)" );
    std::cout << Form("     <Mode> %s (%s)",chain->GetModeName(),basename) << std::endl;
    return make_cut_ee;
  }else if( fl_mode == 12 ){
    chain->SetModeName( "MC(e-mu mode)" );
    std::cout << Form("     <Mode> %s (%s)",chain->GetModeName(),basename) << std::endl;
    return make_cut_emu;
  }else if( fl_mode == 100 ){
    chain->SetModeName( "MC(J/psi[mm] mode(DE-MBC))" );
    std::cout << Form("     <Mode> %s (%s)",chain->GetModeName(),basename) << std::endl;
    return make_cut_jpsi_mm_dembc;
  }else if( fl_mode == 101 ){
    chain->SetModeName( "MC(J/psi[ee] mode(DE-MBC))" );
    std::cout << Form("     <Mode> %s (%s)",chain->GetModeName(),basename) << std::endl;
    return make_cut_jpsi_ee_dembc;
  }else if( fl_mode == 110 ){
    chain->SetModeName( "MC(J/psi[mm] mode)" );
    std::cout << Form("     <Mode> %s (%s)",chain->GetModeName(),basename) << std::endl;
    return make_cut_jpsi_mm;
  }else if( fl_mode == 111 ){
    chain->SetModeName( "MC(J/psi[ee] mode)" );
    std::cout << Form("     <Mode> %s (%s)",chain->GetModeName(),basename) << std::endl;
    return make_cut_jpsi_ee;
  }else{
    std::cout << Form("invalid mode flag(%d, %s)",fl_mode,basename) << std::endl;
    chain->SetModeName( "MC(tentative muon)" );
    return make_cut_mm;
  }
}
