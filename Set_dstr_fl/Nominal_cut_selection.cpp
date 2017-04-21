#include "Nominal_cut_selection.h"

Nominal_Cut nominal_cut_selection( MChain* chain, int fl_mode ){
  const Char_t* basename = gSystem->BaseName( chain->GetFileName() );
  
  if( fl_mode == 0 ){
    chain->SetModeName( "dstr_fl" );
    std::cout << Form("     <Mode> %s (%s)",chain->GetModeName(),basename) << std::endl;
    return make_cut_common;
  }
}
