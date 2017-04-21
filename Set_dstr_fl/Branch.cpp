#include "Branch.h"

MCut_array branch_table(){
  Int_t n = 24;
  MCut_array branch_list(n);
  Int_t i=0;
  //                    <branch-name,tag,type>
  // function-type(1: window, 2:low, 3:high, 4:esc, 5:bool(0or1), 6:bool(-1or+1), 7:window + minus value, (0:optional)
  // detail of function-type is written in "MCut.h"
  branch_list.Init( i++,       "m_d",    0, 1 );
  branch_list.Init( i++,    "m_d_vf",    0, 1 );
  branch_list.Init( i++,     "m_dst",    0, 1 );
  branch_list.Init( i++,  "m_dst_vf",    0, 1 );
  branch_list.Init( i++,        "dm",    0, 1 );
  branch_list.Init( i++,     "dm_vf",    0, 1 );
  branch_list.Init( i++,        "xd",    0, 2 );
  branch_list.Init( i++,     "xd_vf",    0, 2 );
  branch_list.Init( i++,      "xdst",    0, 2 );
  branch_list.Init( i++,   "xdst_vf",    0, 2 );
  branch_list.Init( i++,   "m_drev",     0, 4 );
  branch_list.Init( i++,   "cosd",       0, 3 );
  branch_list.Init( i++,  "npatc1_0",    0, 2 ); // K
  branch_list.Init( i++,  "npatc1_1",    0, 2 ); // Pi
  branch_list.Init( i++,  "npatc1_2",    0, 2 ); // Slow Pi
  branch_list.Init( i++,   "mulmu_0",    0, 2 ); // K
  branch_list.Init( i++,   "mulmu_1",    0, 2 ); // Pi
  branch_list.Init( i++,   "mulmu_2",    0, 2 ); // Slow Pi
  branch_list.Init( i++,  "muchi2_0",    0, 2 ); // K
  branch_list.Init( i++,  "muchi2_1",    0, 2 ); // Pi
  branch_list.Init( i++,  "muchi2_2",    0, 2 ); // Slow Pi
  branch_list.Init( i++,   "eprob_0",    0, 2 ); // K
  branch_list.Init( i++,   "eprob_1",    0, 2 ); // Pi
  branch_list.Init( i++,   "eprob_2",    0, 2 ); // Slow Pi

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
  if( strcmp(h,"h1")==0 ){
    cut->Set(  "m_d_vf", 1, -0.030,  PDGmass::d0,                0.030  );
    cut->Set(  "m_drev", 1, -0.030,  PDGmass::d0,                0.030  ); 
    cut->Set(   "dm_vf", 1, -0.0015, PDGmass::dstrp-PDGmass::d0, 0.0015 );
    cut->Set( "xdst_vf", 1, 0.5 );
    cut->Set(    "cosd", 1, 0.8 );
  }else{
    std::cout << Form("hist-name(%s)",h) << std::endl;
  }
  cut->ClearChange();
  return;
}

