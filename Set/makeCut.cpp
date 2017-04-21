#include "makeCut.h"

const Int_t nmode_ctgry = 18;
const Int_t mode_ctgry[nmode_ctgry] = {  1, 10,
					 101,110,1001,1010,
					 201,210,1101,1110,
					 301,310,1201,1210,
					 401,410,1301,1310
};
const Char_t* mode_name_ctgry[nmode_ctgry] = {"K",    "K_{S}",
					      "K#pi", "K_{S}#pi", "K#pi^{0}",    "K_{S}#pi^{0}",
					      "K2#pi","K_{S}2#pi","K1#pi#pi^{0}","K_{S}1#pi#pi^{0}",
					      "K3#pi","K_{S}3#pi","K2#pi#pi^{0}","K_{S}2#pi#pi^{0}",
					      "K4#pi","K_{S}4#pi","K3#pi#pi^{0}","K_{S}3#pi#pi^{0}",
};


std::string makeCut_Xs( Int_t fl_xsid ){
  std::stringstream sTmp;
  if(      fl_xsid == 1 ) sTmp << "( abs(gm_xsid)==  311 || abs(gm_xsid)==  321 )"; // (K )
  else if( fl_xsid == 2 ) sTmp << "( abs(gm_xsid)==  313 || abs(gm_xsid)==  323 )"; // (K*)
  else if( fl_xsid == 3 ) sTmp << "( abs(gm_xsid)==30343 || abs(gm_xsid)==30353 )"; // (Xs)
  else                    sTmp << "(1)";                                            // (total)
  return sTmp.str();
}

std::string makeCut_mode_category( Int_t fl_mode_xs, Int_t fl_ctgry ){
  std::stringstream sTmp;
  sTmp << "( ";
  sTmp << "rm_xs==" << fl_mode_xs;
  if( fl_ctgry!=-10 ) sTmp << " && " << makeCut_category( fl_ctgry ).c_str();
  sTmp << " )";
  return sTmp.str();
}

std::string makeCut_category( Int_t fl_ctgry ){
  std::stringstream sTmp;
  sTmp << "( ";
  if(      fl_ctgry==0 ) sTmp << "self!=1 && gm_bg==0 && rm_xs==gm_xs"; // false combination
  else if( fl_ctgry==1 ) sTmp << "self==1 && gm_bg==0 && rm_xs==gm_xs"; // true
  else if( fl_ctgry==2 ) std::cerr << "not used -> abort" << std::endl, abort();
  else if( fl_ctgry==3 ){ // off-diagonal modes
    sTmp << "self!=1 && gm_bg==0 && rm_xs!=gm_xs && ( ";
    for( Int_t k=0; k<nmode_ctgry; k++ ){
    sTmp << " gm_xs==" << mode_ctgry[k];
    if( k!= nmode_ctgry-1 ) sTmp << " || ";
    }
    sTmp << " )";
  }else if( fl_ctgry==4 ){ // other-mode
    sTmp << "self!=1 && ( ";
    sTmp << "gm_bg!=0 ||( gm_bg==0 && ";
    for( Int_t k=0; k<nmode_ctgry; k++ ){
      sTmp << " gm_xs!=" << mode_ctgry[k];
      if( k!= nmode_ctgry-1 ) sTmp << " && ";
    }
    sTmp << ")";
    sTmp << " )";
  }else{
    sTmp << 1;
  }
  sTmp << " )";
  return sTmp.str();
}

std::string makeCut_mode_q2fl( Int_t fl_mode_xs, Int_t fl_q2, Int_t fl_flavor, Int_t fl_inv ){
  std::stringstream sTmp;
  sTmp << "( ";
  sTmp << "rm_xs==" << fl_mode_xs;
  if( fl_q2!=-10 || fl_flavor!=-10 ) sTmp << " && " << makeCut_q2fl( fl_q2, fl_flavor, fl_inv ).c_str();
  sTmp << " )";
  return sTmp.str();
}

std::string makeCut_q2fl( Int_t fl_q2, Int_t fl_flavor, Int_t fl_inv ){
  std::stringstream sTmp;
  if( fl_inv ) sTmp << "( !(";
  else         sTmp << "   (";
  
  if(      fl_q2==1 && fl_flavor==1 ) sTmp << "q2self==1 && recbfl==genbfl";
  else if( fl_q2==1 && fl_flavor==2 ) sTmp << "q2self==1 && recbfl!=genbfl && recbfl==2";
  else if( fl_q2==1 && fl_flavor==0 ) sTmp << "q2self==1 && recbfl!=genbfl && recbfl!=2";
  else if( fl_q2==0 && fl_flavor==1 ) sTmp << "q2self!=1 && recbfl==genbfl";
  else if( fl_q2==0 && fl_flavor==2 ) sTmp << "q2self!=1 && recbfl!=genbfl && recbfl==2";
  else if( fl_q2==0 && fl_flavor==0 ) sTmp << "q2self!=1 && recbfl!=genbfl && recbfl!=2";
  else if( fl_q2==1                 ) sTmp << "q2self==1";
  else if( fl_q2==0                 ) sTmp << "q2self!=1";
  else if(             fl_flavor==1 ) sTmp << "recbfl==genbfl";
  else if(             fl_flavor==2 ) sTmp << "recbfl!=genbfl && recbfl==2";
  else if(             fl_flavor==0 ) sTmp << "recbfl!=genbfl && recbfl!=2";
  else                                sTmp << 1;
  
  
  if( fl_inv ) sTmp << " ) )";
  else         sTmp << " )  ";
  return sTmp.str();
}

std::string makeCut_pi0( Int_t fl_mode_pi0, Int_t fl_self ){

  std::stringstream sTmp;
  sTmp << "( ";
  const Int_t add[2][nmode_ctgry] ={
    {1,1, 1,1,0,0, 1,1,0,0, 1,1,0,0, 1,1,0,0,}, // w/o pi0
    {0,0, 0,0,1,1, 0,0,1,1, 0,0,1,1, 0,0,1,1,}, // w/  pi0
  };
  for( Int_t i=0; i<nmode_ctgry; i++ ){
    if( add[fl_mode_pi0][i] ){
      sTmp << "rm_xs==" << mode_ctgry[i];
      if( !(i==15 || i==17) ) sTmp << " || ";
    }
  }
  sTmp << " )";
  if( fl_self ) sTmp << " && self==1 && gm_bg==0 && gm_xs==rm_xs"; // true
  return sTmp.str();
}

std::string makeCut_track( Int_t fl_mode_track, Int_t fl_self ){
  std::stringstream sTmp;
  sTmp << "( ";
  const Int_t add[6][nmode_ctgry] ={
    {0,1, 0,0,0,1, 0,0,0,0, 0,0,0,0, 0,0,0,0,}, // 0 track
    {1,0, 0,1,1,0, 0,0,0,1, 0,0,0,0, 0,0,0,0,}, // 1 track
    {0,0, 1,0,0,0, 0,1,1,0, 0,0,0,1, 0,0,0,0,}, // 2 track
    {0,0, 0,0,0,0, 1,0,0,0, 0,1,1,0, 0,0,0,1,}, // 3 track
    {0,0, 0,0,0,0, 0,0,0,0, 1,0,0,0, 0,1,1,0,}, // 4 track
    {0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,0,0,}, // 5 track
  };
  for( Int_t i=0; i<nmode_ctgry; i++ ){
    if( add[fl_mode_track][i] ){
      sTmp << "rm_xs==" << mode_ctgry[i];
      if( !(i==5 || i==9 || i==13 || i==17 || i==16 || i==14) ) sTmp << " || ";
    }
  }
  sTmp << " )";
  if( fl_self ) sTmp << " && self==1 && gm_bg==0 && gm_xs==rm_xs"; // true
  return sTmp.str();
}

std::string makeCut_body( Int_t fl_mode_body, Int_t fl_self ){
  std::stringstream sTmp;
  sTmp << "( ";
  const Int_t add[6][nmode_ctgry] ={
    {0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,}, // 0 body(none)
    {1,1, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,}, // 1 body
    {0,0, 1,1,1,1, 0,0,0,0, 0,0,0,0, 0,0,0,0,}, // 2 body
    {0,0, 0,0,0,0, 1,1,1,1, 0,0,0,0, 0,0,0,0,}, // 3 body
    {0,0, 0,0,0,0, 0,0,0,0, 1,1,1,1, 0,0,0,0,}, // 4 body
    {0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 1,1,1,1,}, // 5 body
  };
  if( fl_mode_body==0 ) std::cerr << "0 body mode does not exist -> abort()" << std::endl, abort();
  for( Int_t i=0; i<nmode_ctgry; i++ ){
    if( add[fl_mode_body][i] ){
      sTmp << "rm_xs==" << mode_ctgry[i];
      if( !(i==1 || i==5 || i==9 || i==13 || i==17) ) sTmp << " || ";
    }
  }
  sTmp << " )";
  if( fl_self ) sTmp << " && self==1 && gm_bg==0 && gm_xs==rm_xs"; // true
  return sTmp.str();
}


std::string makeCut_5body_veto(){
  std::stringstream sTmp;
  sTmp << "( !"
       << makeCut_body(5)
       << " )";

  return sTmp.str();
}

std::string makeCut_unflavor_veto(){
  std::stringstream sTmp;
  /*
  sTmp << "( "
       << "rm_xs!= 10 &&"
       << "rm_xs!= 1010 &&"
       << "rm_xs!= 210 &&"
       << "rm_xs!= 1210 &&"
       << "rm_xs!= 410 "
       << " )";
  */
  sTmp << "( recbfl!=2 )";
  return sTmp.str();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

std::string makeCut_LRNB_internal( Int_t    fl_q2,  Int_t    fl_xs, Int_t fl_mode_ll,
				   Double_t q2,     Double_t xs ){
  std::stringstream sTmp;
  sTmp << "( "; // BEGIN
  
  // [ q2 region cut ]
  if(      fl_q2 == 0 ) sTmp << " 1 " ;                     //  no cut
  else if( fl_q2 >  0 ) sTmp << "cc_m>sqrt(" << q2 << ")";  //  high q2 region
  else if( fl_q2 <  0 ) sTmp << "cc_m<sqrt(" << q2 << ")";  //   low q2 retion
  
  sTmp << " && ";
  
  // [ lepton type ]
  if(      fl_mode_ll==1 ) sTmp << "rm_l==1"; // e
  else if( fl_mode_ll==0 ) sTmp << "rm_l==0"; // mu
  else                     sTmp << " 1 ";     // no-cut
  
  sTmp << " && ";
  
  // [ xs region cut ]
  if(      fl_xs == 0 ) sTmp << " 1 ";         // no cut
  else if( fl_xs >  0 ) sTmp << "xs_m>" << xs; // high xs region
  else if( fl_xs <  0 ) sTmp << "xs_m<" << xs; //  low xs region

  sTmp << " )"; // END
  
  return sTmp.str();
}

std::string makeCut_LRNB_1d_internal( Char_t*  brname, Double_t var,
				      Int_t    fl_q2,  Int_t    fl_xs, Int_t fl_mode_ll,
				      Double_t q2,     Double_t xs ){
  std::stringstream sTmp;
  sTmp << "( "; // BEGIN
  

  sTmp << makeCut_LRNB_internal( fl_q2, fl_xs, fl_mode_ll, q2, xs ) << " && ";

  //[ LRNB cut ]
  sTmp << Form( brname, "tot" ) << ">" << var;
  
  sTmp << " )"; // END
  return sTmp.str();
}

std::string makeCut_LRNB_1d_lep_internal( Char_t*  brname, Double_t var,
					  Int_t    fl_q2,  Int_t    fl_xs, Int_t fl_mode_ll,
					  Double_t q2,     Double_t xs ){
  std::stringstream sTmp;
  sTmp << "( "; // BEGIN
  
  
  sTmp << makeCut_LRNB_internal( fl_q2, fl_xs, fl_mode_ll, q2, xs ) << " && ";
  
  //[ LRNB cut ]
  sTmp << Form( brname, fl_mode_ll, "tot" ) << ">" << var;
  
  sTmp << " )"; // END
  return sTmp.str();
}

std::string makeCut_LRNB_2d_internal( Char_t*  brname, Double_t var_qq, Double_t var_bb,
				      Int_t    fl_q2,  Int_t    fl_xs, Int_t fl_mode_ll,
				      Double_t q2,     Double_t xs ){
  std::stringstream sTmp;
  sTmp << "( "; // BEGIN
  sTmp << makeCut_LRNB_internal( fl_q2, fl_xs, fl_mode_ll, q2, xs ) << " && ";  
  
  //[ LRNB cut ]
  sTmp << Form( brname,"qq" ) << ">" << var_qq;
  sTmp << " && ";
  sTmp << Form( brname,"bb" ) << ">" << var_bb;
  
  sTmp << " )"; // END
  return sTmp.str();
}

std::string makeCut_LRNB_2d_lep_internal( Char_t*  brname, Double_t var_qq, Double_t var_bb,
					  Int_t    fl_q2,  Int_t    fl_xs,  Int_t fl_mode_ll,
					  Double_t q2,     Double_t xs ){
  std::stringstream sTmp;
  sTmp << "( "; // BEGIN
  sTmp << makeCut_LRNB_internal( fl_q2, fl_xs, fl_mode_ll, q2, xs ) << " && ";  
  
  //[ LRNB cut ]
  sTmp << Form( brname, fl_mode_ll, "qq" ) << ">" << var_qq;
  sTmp << " && ";
  sTmp << Form( brname, fl_mode_ll, "bb" ) << ">" << var_bb;
  
  sTmp << " )"; // END
  return sTmp.str();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

std::string makeCut_LRNB_1d( Char_t*   brname,   Int_t    fl_q2,
			     Double_t var_lowxs, Double_t var_highxs,
			     Double_t q2,        Double_t xs ){
  std::stringstream sTmp;
  sTmp << "( "; // BEGIN
  
  sTmp << makeCut_LRNB_1d_internal(brname, var_lowxs,  fl_q2, -1, 2, q2, xs ).c_str()  //  low-xs
       << " || "
       << makeCut_LRNB_1d_internal(brname, var_highxs, fl_q2,  1, 2, q2, xs ).c_str(); // high-xs
  
  sTmp << " )"; // END
  return sTmp.str();
}

std::string makeCut_LRNB_1d_lep( Char_t*  brname,       Int_t    fl_q2,
				 Double_t var_lowxs_ee, Double_t var_highxs_ee,
				 Double_t var_lowxs_mm, Double_t var_highxs_mm,
				 Double_t q2,           Double_t xs ){
  std::stringstream sTmp;
  sTmp << "( "; // BEGIN
  
  sTmp << makeCut_LRNB_1d_lep_internal(brname, var_lowxs_ee,  fl_q2, -1, 1, q2, xs ).c_str()  //  low-xs, ee
       << " || "
       << makeCut_LRNB_1d_lep_internal(brname, var_highxs_ee, fl_q2,  1, 1, q2, xs ).c_str()  // high-xs, ee
       << " || "
       << makeCut_LRNB_1d_lep_internal(brname, var_lowxs_mm,  fl_q2, -1, 0, q2, xs ).c_str()  //  low-xs, mm
       << " || "
       << makeCut_LRNB_1d_lep_internal(brname, var_highxs_mm, fl_q2,  1, 0, q2, xs ).c_str(); // high-xs, mm
  
  sTmp << " )"; // END
  return sTmp.str();
}

std::string makeCut_LRNB_2d( Char_t*  brname,        Int_t    fl_q2,
			     Double_t var_lowxs_qq,  Double_t var_lowxs_bb,
			     Double_t var_highxs_qq, Double_t var_highxs_bb,
			     Double_t q2,            Double_t xs ){
  std::stringstream sTmp;
  sTmp << "( "; // BEGIN
  
  sTmp << makeCut_LRNB_2d_internal(brname, var_lowxs_qq,  var_lowxs_bb,  fl_q2, -1, 2, q2, xs ).c_str()  //  low-xs
       << " || "
       << makeCut_LRNB_2d_internal(brname, var_highxs_qq, var_highxs_bb, fl_q2,  1, 2, q2, xs ).c_str(); // high-xs
  
  sTmp << " )"; // END
  return sTmp.str();
}

std::string makeCut_LRNB_2d_lep( Char_t* brname, Int_t fl_q2,
				 Double_t var_lowxs_qq_ee, Double_t var_lowxs_bb_ee, Double_t var_highxs_qq_ee, Double_t var_highxs_bb_ee,
				 Double_t var_lowxs_qq_mm, Double_t var_lowxs_bb_mm, Double_t var_highxs_qq_mm, Double_t var_highxs_bb_mm,
				 Double_t q2, Double_t xs ){
  std::stringstream sTmp;
  sTmp << "( "; // BEGIN
  
  sTmp << makeCut_LRNB_2d_lep_internal(brname, var_lowxs_qq_ee,  var_lowxs_bb_ee,  fl_q2, -1, 1, q2, xs ).c_str()  //  low-xs, ee
       << " || "
       << makeCut_LRNB_2d_lep_internal(brname, var_highxs_qq_ee, var_highxs_bb_ee, fl_q2,  1, 1, q2, xs ).c_str()  // high-xs, ee
       << " || "
       << makeCut_LRNB_2d_lep_internal(brname, var_lowxs_qq_mm,  var_lowxs_bb_mm,  fl_q2, -1, 0, q2, xs ).c_str()  //  low-xs, mm
       << " || "
       << makeCut_LRNB_2d_lep_internal(brname, var_highxs_qq_mm, var_highxs_bb_mm, fl_q2,  1, 0, q2, xs ).c_str(); // high-xs, mm
  
  sTmp << " )"; // END
  return sTmp.str();
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

std::string makeCut_q2( Int_t fl_mode_ll, Int_t fl_q2, Int_t fl_cos ){
  if( fl_q2<=0 || fl_q2>q2_bin || fl_q2==4 || fl_q2==6 ) std::cerr << "[ABORT] Invalid q2 region : " << fl_q2 << std::endl
								   << "        fl_q2 = {1,2,3,5,7,8,9}"       << std::endl, abort();
								     
  std::stringstream sTmp;
  sTmp << "( "; // BEGIN
  
  sTmp << q2_bins[fl_mode_ll][fl_q2-1] << " < cc_m*cc_m " << "&& "
       << q2_bins[fl_mode_ll][fl_q2  ] << " > cc_m*cc_m ";
  if     ( fl_cos>0 ) sTmp << " && coslp>0";
  else if( fl_cos<0 ) sTmp << " && coslp<0";
  
  sTmp << " )"; // END
  
  return sTmp.str();
}

Int_t tag_q2_cos_region(Int_t fl_mode_ll, Double_t q2, Double_t cos){
  Int_t tag = -10;
  if     ( q2 > q2_bins[fl_mode_ll][q2_bin] ) tag = -4; // too large q2
  else if( q2 < q2_bins[fl_mode_ll][0]      ) tag = -3; // negative  q2
  else if( q2 > q2_bins[fl_mode_ll][0] && q2 < q2_bins[fl_mode_ll][1] && cos > 0 ) tag =  1; // 1st q2 bin, +cos
  else if( q2 > q2_bins[fl_mode_ll][0] && q2 < q2_bins[fl_mode_ll][1] && cos < 0 ) tag =  2; // 1st q2 bin, -cos
  else if( q2 > q2_bins[fl_mode_ll][1] && q2 < q2_bins[fl_mode_ll][2] && cos > 0 ) tag =  3; // 2nd q2 bin, +cos
  else if( q2 > q2_bins[fl_mode_ll][1] && q2 < q2_bins[fl_mode_ll][2] && cos < 0 ) tag =  4; // 2nd q2 bin, -cos
  else if( q2 > q2_bins[fl_mode_ll][2] && q2 < q2_bins[fl_mode_ll][3] && cos > 0 ) tag =  5; // 3rd q2 bin, +cos
  else if( q2 > q2_bins[fl_mode_ll][2] && q2 < q2_bins[fl_mode_ll][3] && cos < 0 ) tag =  6; // 3rd q2 bin, -cos
  else if( q2 > q2_bins[fl_mode_ll][3] && q2 < q2_bins[fl_mode_ll][4]            ) tag = -1; //    jpsi veto 
  else if( q2 > q2_bins[fl_mode_ll][4] && q2 < q2_bins[fl_mode_ll][5] && cos > 0 ) tag =  7; // 4th q2 bin, +cos
  else if( q2 > q2_bins[fl_mode_ll][4] && q2 < q2_bins[fl_mode_ll][5] && cos < 0 ) tag =  8; // 4th q2 bin, -cos
  else if( q2 > q2_bins[fl_mode_ll][5] && q2 < q2_bins[fl_mode_ll][6]            ) tag = -2; //   psi(2s) veto
  else if( q2 > q2_bins[fl_mode_ll][6] && q2 < q2_bins[fl_mode_ll][7] && cos > 0 ) tag =  9; // 4th q2 bin, +cos
  else if( q2 > q2_bins[fl_mode_ll][6] && q2 < q2_bins[fl_mode_ll][7] && cos < 0 ) tag = 10; // 4th q2 bin, -cos
  else if( q2 > q2_bins[fl_mode_ll][7] && q2 < q2_bins[fl_mode_ll][8] && cos > 0 ) tag = 11; // 5th q2 bin, +cos
  else if( q2 > q2_bins[fl_mode_ll][7] && q2 < q2_bins[fl_mode_ll][8] && cos < 0 ) tag = 12; // 5th q2 bin, -cos
  else if( q2 > q2_bins[fl_mode_ll][8] && q2 < q2_bins[fl_mode_ll][9] && cos > 0 ) tag = 13; // 6th q2 bin, +cos
  else if( q2 > q2_bins[fl_mode_ll][8] && q2 < q2_bins[fl_mode_ll][9] && cos < 0 ) tag = 14; // 6th q2 bin, -cos

  return tag;
}

Int_t tag_q2_cos_region2(Int_t fl_mode_ll, Double_t q2, Double_t cos){
  Int_t tag = -10;
  if     ( q2 > q2_bins[fl_mode_ll][q2_bin] ) tag = -4; // too large q2
  else if( q2 < q2_bins[fl_mode_ll][0]      ) tag = -3; // negative  q2
  else if( q2 > q2_bins[fl_mode_ll][0] && q2 < q2_bins[fl_mode_ll][1] && cos > 0 ) tag =  1; // 1st q2 bin, +cos
  else if( q2 > q2_bins[fl_mode_ll][0] && q2 < q2_bins[fl_mode_ll][1] && cos < 0 ) tag =  2; // 1st q2 bin, -cos
  else if( q2 > q2_bins[fl_mode_ll][1] && q2 < q2_bins[fl_mode_ll][2] && cos > 0 ) tag =  1; // 2nd q2 bin, +cos
  else if( q2 > q2_bins[fl_mode_ll][1] && q2 < q2_bins[fl_mode_ll][2] && cos < 0 ) tag =  2; // 2nd q2 bin, -cos
  else if( q2 > q2_bins[fl_mode_ll][2] && q2 < q2_bins[fl_mode_ll][3] && cos > 0 ) tag =  3; // 3rd q2 bin, +cos
  else if( q2 > q2_bins[fl_mode_ll][2] && q2 < q2_bins[fl_mode_ll][3] && cos < 0 ) tag =  4; // 3rd q2 bin, -cos
  else if( q2 > q2_bins[fl_mode_ll][3] && q2 < q2_bins[fl_mode_ll][4]            ) tag = -1; //    jpsi veto 
  else if( q2 > q2_bins[fl_mode_ll][4] && q2 < q2_bins[fl_mode_ll][5] && cos > 0 ) tag =  5; // 4th q2 bin, +cos
  else if( q2 > q2_bins[fl_mode_ll][4] && q2 < q2_bins[fl_mode_ll][5] && cos < 0 ) tag =  6; // 4th q2 bin, -cos
  else if( q2 > q2_bins[fl_mode_ll][5] && q2 < q2_bins[fl_mode_ll][6]            ) tag = -2; //   psi(2s) veto
  else if( q2 > q2_bins[fl_mode_ll][6] && q2 < q2_bins[fl_mode_ll][7] && cos > 0 ) tag =  7; // 4th q2 bin, +cos
  else if( q2 > q2_bins[fl_mode_ll][6] && q2 < q2_bins[fl_mode_ll][7] && cos < 0 ) tag =  8; // 4th q2 bin, -cos
  else if( q2 > q2_bins[fl_mode_ll][7] && q2 < q2_bins[fl_mode_ll][8] && cos > 0 ) tag =  7; // 5th q2 bin, +cos
  else if( q2 > q2_bins[fl_mode_ll][7] && q2 < q2_bins[fl_mode_ll][8] && cos < 0 ) tag =  8; // 5th q2 bin, -cos
  else if( q2 > q2_bins[fl_mode_ll][8] && q2 < q2_bins[fl_mode_ll][9] && cos > 0 ) tag =  7; // 6th q2 bin, +cos
  else if( q2 > q2_bins[fl_mode_ll][8] && q2 < q2_bins[fl_mode_ll][9] && cos < 0 ) tag =  8; // 6th q2 bin, -cos

  return tag;
}

Int_t tag_q2_cos_region_totFB(Int_t fl_mode_ll, Double_t q2, Double_t cos){
  Int_t tag = -10;
  if     ( q2 > q2_bins[fl_mode_ll][q2_bin] ) tag = -4; // too large q2
  else if( q2 < q2_bins[fl_mode_ll][0]      ) tag = -3; // negative  q2
  else if( q2 > q2_bins[fl_mode_ll][0] && q2 < q2_bins[fl_mode_ll][1] && cos > 0 ) tag =  1; // 1st q2 bin, +cos
  else if( q2 > q2_bins[fl_mode_ll][0] && q2 < q2_bins[fl_mode_ll][1] && cos < 0 ) tag =  2; // 1st q2 bin, -cos
  else if( q2 > q2_bins[fl_mode_ll][1] && q2 < q2_bins[fl_mode_ll][2] && cos > 0 ) tag =  1; // 2nd q2 bin, +cos
  else if( q2 > q2_bins[fl_mode_ll][1] && q2 < q2_bins[fl_mode_ll][2] && cos < 0 ) tag =  2; // 2nd q2 bin, -cos
  else if( q2 > q2_bins[fl_mode_ll][2] && q2 < q2_bins[fl_mode_ll][3] && cos > 0 ) tag =  1; // 3rd q2 bin, +cos
  else if( q2 > q2_bins[fl_mode_ll][2] && q2 < q2_bins[fl_mode_ll][3] && cos < 0 ) tag =  2; // 3rd q2 bin, -cos
  else if( q2 > q2_bins[fl_mode_ll][3] && q2 < q2_bins[fl_mode_ll][4]            ) tag = -1; //    jpsi veto 
  else if( q2 > q2_bins[fl_mode_ll][4] && q2 < q2_bins[fl_mode_ll][5] && cos > 0 ) tag =  1; // 4th q2 bin, +cos
  else if( q2 > q2_bins[fl_mode_ll][4] && q2 < q2_bins[fl_mode_ll][5] && cos < 0 ) tag =  2; // 4th q2 bin, -cos
  else if( q2 > q2_bins[fl_mode_ll][5] && q2 < q2_bins[fl_mode_ll][6]            ) tag = -2; //   psi(2s) veto
  else if( q2 > q2_bins[fl_mode_ll][6] && q2 < q2_bins[fl_mode_ll][7] && cos > 0 ) tag =  1; // 4th q2 bin, +cos
  else if( q2 > q2_bins[fl_mode_ll][6] && q2 < q2_bins[fl_mode_ll][7] && cos < 0 ) tag =  2; // 4th q2 bin, -cos
  else if( q2 > q2_bins[fl_mode_ll][7] && q2 < q2_bins[fl_mode_ll][8] && cos > 0 ) tag =  1; // 5th q2 bin, +cos
  else if( q2 > q2_bins[fl_mode_ll][7] && q2 < q2_bins[fl_mode_ll][8] && cos < 0 ) tag =  2; // 5th q2 bin, -cos
  else if( q2 > q2_bins[fl_mode_ll][8] && q2 < q2_bins[fl_mode_ll][9] && cos > 0 ) tag =  1; // 6th q2 bin, +cos
  else if( q2 > q2_bins[fl_mode_ll][8] && q2 < q2_bins[fl_mode_ll][9] && cos < 0 ) tag =  2; // 6th q2 bin, -cos

  return tag;
}
