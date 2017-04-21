#include "Nominal_cut_selection.h"

Nominal_Cut nominal_cut_selection( MChain* chain, int fl_mode ){
  const Char_t* basename = gSystem->BaseName( chain->GetFileName() );


  if( fl_mode == 1 ){
    chain->SetModeName( "MC(electron) mode [de&Mbc]" );
    std::cout << Form("     <Mode> %s (%s)",chain->GetModeName(),basename) << std::endl;
    return make_cut_ee_dembc;
  }else if( fl_mode == 0 ){
    chain->SetModeName( "MC(muon)  mode[de&Mbc])" );
    std::cout << Form("     <Mode> %s (%s)",chain->GetModeName(),basename) << std::endl;
    return make_cut_mm_dembc;
  }else if( fl_mode == 110 ){
    chain->SetModeName( "MC(double lepton veto(e+e-) mode[fake lepton])" );
    std::cout << Form("     <Mode> %s (%s)",chain->GetModeName(),basename) << std::endl;
    return make_cut_double_ee;
  }else if( fl_mode == 111 ){
    chain->SetModeName( "MC(single lepton veto(e+)  mode[fake lepton])" );
    std::cout << Form("     <Mode> %s (%s)",chain->GetModeName(),basename) << std::endl;
    return make_cut_single_ep;
  }else if( fl_mode == -111 ){
    chain->SetModeName( "MC(single lepton veto(e-)  mode[fake lepton])" );
    std::cout << Form("     <Mode> %s (%s)",chain->GetModeName(),basename) << std::endl;
    return make_cut_single_em;
  }else if( fl_mode == 100 ){
    chain->SetModeName( "MC(double lepton veto(mu+mu-) mode[fake lepton])" );
    std::cout << Form("     <Mode> %s (%s)",chain->GetModeName(),basename) << std::endl;
    return make_cut_double_mm;
  }else if( fl_mode == 101 ){
    chain->SetModeName( "MC(single lepton veto(mu+)  mode[fake lepton])" );
    std::cout << Form("     <Mode> %s (%s)",chain->GetModeName(),basename) << std::endl;
    return make_cut_single_mp;
  }else if( fl_mode == -101 ){
    chain->SetModeName( "MC(single lepton veto(mu-)  mode[fake lepton])" );
    std::cout << Form("     <Mode> %s (%s)",chain->GetModeName(),basename) << std::endl;
    return make_cut_single_mm;
  }else{
    std::cout << Form("invalid mode flag(%d, %s)",fl_mode,basename) << std::endl;
    abort();
  }
}
