#ifndef DRAWS_H
#define DRAWS_H

#include <iostream>
#include <sstream>
#include <TROOT.h>

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// if this section is changed, "Set/makeCut.cpp" must be changed.
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
const Int_t nmode       = 18;
const Int_t nmode_afb   = 10;
const Int_t mode[nmode] = {  1, 10,
			   101,110,1001,1010,
			   201,210,1101,1110,
			   301,310,1201,1210,
			   401,410,1301,1310
};
const Char_t* mode_name[nmode] = {"K",    "K_{S}",
				  "K#pi", "K_{S}#pi", "K#pi^{0}",    "K_{S}#pi^{0}",
				  "K2#pi","K_{S}2#pi","K1#pi#pi^{0}","K_{S}1#pi#pi^{0}",
				  "K3#pi","K_{S}3#pi","K2#pi#pi^{0}","K_{S}2#pi#pi^{0}",
				  "K4#pi","K_{S}4#pi","K3#pi#pi^{0}","K_{S}3#pi#pi^{0}",
};

const Int_t n_k  [nmode] = {1,0,
			    1,0,1,0,
			    1,0,1,0,
			    1,0,1,0,
			    1,0,1,0,
};
const Int_t n_pi[nmode] = {0,0,
			   1,1,0,0,
			   2,2,1,1,
			   3,3,2,2,
			   4,4,3,3,
};
const Int_t n_pi0[nmode] = {0,0,
			    0,0,1,1,
			    0,0,1,1,
			    0,0,1,1,
			    0,0,1,1,
};

const Int_t mode_afb[nmode_afb] = {  1,
				     101,110,1001,
				     201,    1101,1110,
				     301,310,1201,
};

const Char_t* mode_name_afb[nmode_afb] = {"K",
					  "K#pi", "K_{S}#pi", "K#pi^{0}",
					  "K2#pi",            "K1#pi#pi^{0}","K_{S}1#pi#pi^{0}",
					  "K3#pi","K_{S}3#pi","K2#pi#pi^{0}",
};

const Int_t n_k_afb[nmode_afb] = {1,
				  1,0,1,
				  1,  1,0,
				  1,0,1,
};
const Int_t n_pi_afb[nmode_afb] = {0,
				   1,1,0,
				   2,  1,1,
				   3,3,2,
};
const Int_t n_pi0_afb[nmode_afb] = {0,
				    0,0,1,
				    0,  1,1,
				    0,0,1,
};

const Int_t   nbkgtype          = 4;
const Char_t* bkgtype[nbkgtype] = {"uds","charm","mixed","charged"};

const Int_t   nexpno = 30;
const Char_t* expno[nexpno]  = {"07","09","11","13","15","17","19a","19b","21","23","25","27",
				"31","33","35","37a","37b","39","41a","41b","43","45a","45b","47","49","51","55","61","63","65"};

const Int_t   nexpno_gmc = 26;
const Char_t* expno_gmc[nexpno_gmc]  = {"07","09","11","13","15","17","19","21","23","25","27",
					"31","33","35","37","39","41","43","45","47","49","51","55","61","63","65"};
const Int_t expno_gmc_int[nexpno_gmc]  = {7,9,11,13,15,17,19,21,23,25,27,
					  31,33,35,37,39,41,43,45,47,49,51,55,61,63,65};
const Double_t nbb_gmc[nexpno_gmc] = {6.4588, 4.7597, 8.8509, 11.6998, 13.5679, 12.4588, 16.3023+10.8682, 4.3371, 6.4755, 28.0008, 28.1814,
				      19.6601, 18.3460, 17.3783, 30.5480+36.6339, 47.0818, 34.0657+29.9477, 61.5614, 6.64662+7.70718,
				      41.2186, 29.7205, 41.8919, 80.2472,
				      37.4460, 35.6231, 41.7867
}; // from nBB.txt in mcprod(sigMC)

const Double_t nbb_ccmc[nexpno_gmc] = {6.4587, 4.7597, 8.8509, 11.6998, 13.5679, 12.4588, 16.0577655+11.1127345, 4.3371, 6.4755, 28.0008, 28.1814,
				       19.6587, 19.3022, 18.5262, 30.5677645+36.6141355, 47.0818, 34.0551288+29.9582712, 61.5614, 5.0525376+9.3012624,
				       41.2186 , 29.7271, 41.8919, 80.2472,
				       37.4460, 35.6231, 41.7867 // old value is used ?
};
// Difference( exp : 31, 33, 35, 49)


const Char_t* lrnb[2] = {"lr", "nb"};

const Double_t AFB_gen_6qbin[3][4] = {
  {-0.0839, 0.0901, 0.2522, 0.2539}, // ee
  {-0.0861, 0.1241, 0.2654, 0.2606}, // mm
  {-0.0850, 0.1098, 0.2609, 0.2573}, // ee+mm
};

const Double_t AFB_gen_6qbinE[3][4] = {
  {0.0039, 0.0066, 0.0111, 0.0080}, // ee
  {0.0041, 0.0056, 0.0080, 0.0079}, // mm
  {0.0028, 0.0043, 0.0065, 0.0056}, // ee+mm
};
/* calculated from Gen*Eff
const Double_t AFB_eff_6qbin[3][4] = {
  {-0.0614, 0.0791, 0.2432, 0.2291}, // ee
  {-0.0384, 0.0992, 0.2549, 0.2404}, // mm
  {-0.0516, 0.0904, 0.2515, 0.2354}, // ee+mm
};
*/
///* calculated from Rec
const Double_t AFB_eff_6qbin[3][4] = {
  {-0.0452, 0.0341, 0.2733, 0.2276}, // ee
  {-0.0216, 0.0879, 0.2677, 0.2325}, // mm
  {-0.0350, 0.0887, 0.2693, 0.2304}, // ee+mm
};
//*/
const Double_t AFB_eff_6qbinE[3][4] = {
  {0.0049, 0.0072, 0.0119, 0.0087}, // ee
  {0.0057, 0.0064, 0.0084, 0.0084}, // mm
  {0.0037, 0.0048, 0.0069, 0.0061}, // ee+mm
};
const Double_t AFB_gen_9qbin[3][7] = {
  {-0.0927, -0.0655, 0.0901, 0.2522, 0.2931, 0.2578, 0.0254}, // ee
  {-0.1031, -0.0571, 0.1241, 0.2654, 0.3031, 0.2548, 0.0306}, // mm
  {-0.0975, -0.0612, 0.1098, 0.2609, 0.2982, 0.2563, 0.0278}, // ee+mm
};

const Double_t AFB_gen_9qbinE[3][7] = {
  {0.0048, 0.0070, 0.0066, 0.0111, 0.0119, 0.0117, 0.0286}, // ee
  {0.0052, 0.0068, 0.0056, 0.0080, 0.0117, 0.0115, 0.0309}, // mm
  {0.0035, 0.0048, 0.0043, 0.0065, 0.0084, 0.0082, 0.0210}, // ee+mm
};
///* calculated from Gen*Eff
const Double_t AFB_eff_9qbin[3][7] = {
  {-0.0673, -0.0537, 0.0791, 0.2432, 0.2743, 0.2413, 0.0076}, // ee
  {-0.0503, -0.0247, 0.0992, 0.2549, 0.2888, 0.2401, 0.0340}, // mm
  {-0.0602, -0.0408, 0.0904, 0.2515, 0.2825, 0.2407, 0.0215}, // ee+mm
};
//*/
const Double_t AFB_eff_9qbinE[3][7] = {
  {0.0062, 0.0078, 0.0072, 0.0119, 0.0127, 0.0124, 0.0313}, // ee
  {0.0076, 0.0086, 0.0064, 0.0084, 0.0123, 0.0120, 0.0321}, // mm
  {0.0048, 0.0058, 0.0048, 0.0069, 0.0089, 0.0087, 0.0225}, // ee+mm
};

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
const Char_t* cut_ccpi0   = "!( (rm_xs==1001 || rm_xs==1010 || rm_xs==1101 || rm_xs==1110) && (rm_l==1 && -0.15<cc_mheg -3.096916 && cc_mheg -3.096916<0.05) )";
/*
const Char_t* cut_ccpi0_heg        = "!( (rm_xs==1001 || rm_xs==1010 || rm_xs==1101 || rm_xs==1110) && ((rm_l==1 && -0.40<cc_mheg -3.096916 && cc_mheg -3.096916<0.15) || (rm_l==1 && -0.25<cc_mheg -3.68609  && cc_mheg -3.68609 <0.10) || (rm_l==0 && -0.25<cc_mheg  -3.096916 && cc_mheg  -3.096916<0.10) || (rm_l==0 && -0.15<cc_mheg -3.68609  && cc_mheg -3.68609 <0.10) ) )";
const Char_t* cut_ccorgpi0_heg     = "!( (rm_xs==1001 || rm_xs==1010 || rm_xs==1101 || rm_xs==1110) && ((rm_l==1 && -0.40<cc_morgh-3.096916 && cc_morgh-3.096916<0.15) || (rm_l==1 && -0.25<cc_morgh-3.68609  && cc_morgh-3.68609 <0.10) || (rm_l==0 && -0.25<cc_morgh -3.096916 && cc_morgh -3.096916<0.10) || (rm_l==0 && -0.15<cc_morgh-3.68609  && cc_morgh-3.68609 <0.10) ) )";
const Char_t* cut_ccpi0_heg_jpsi   = "!( (rm_xs==1001 || rm_xs==1010 || rm_xs==1101 || rm_xs==1110) && ((rm_l==1 && -0.40<cc_mheg -3.096916 && cc_mheg -3.096916<0.15) || (rm_l==0 && -0.25<cc_mheg -3.096916 && cc_mheg -3.096916<0.10)) )";
const Char_t* cut_ccpi0_heg_jpsi2  = "!( (rm_xs==1001 || rm_xs==1010 || rm_xs==1101 || rm_xs==1110) && ((rm_l==1 && -0.15<cc_mheg -3.096916 && cc_mheg -3.096916<0.05) || (rm_l==0 && -0.15<cc_mheg -3.096916 && cc_mheg -3.096916<0.05)) )";
const Char_t* cut_ccpi0_heg_jpsi3  = "!( (rm_xs==1001 || rm_xs==1010 || rm_xs==1101 || rm_xs==1110) && ((rm_l==1 && -0.10<cc_mheg -3.096916 && cc_mheg -3.096916<0.05) || (rm_l==0 && -0.10<cc_mheg -3.096916 && cc_mheg -3.096916<0.05)) )";
const Char_t* cut_ccpi0_heg_jpsi4  = "!( (rm_xs==1001 || rm_xs==1010 || rm_xs==1101 || rm_xs==1110) && ((rm_l==1 && -0.15<cc_mheg -3.096916 && cc_mheg -3.096916<0.10) || (rm_l==0 && -0.15<cc_mheg -3.096916 && cc_mheg -3.096916<0.10)) )";
const Char_t* cut_ccpi0_heg_psi2s  = "!( (rm_xs==1001 || rm_xs==1010 || rm_xs==1101 || rm_xs==1110) && ((rm_l==1 && -0.25<cc_mheg -3.68609  && cc_mheg -3.68609 <0.10) || (rm_l==0 && -0.15<cc_mheg -3.68609  && cc_mheg -3.68609 <0.10)) )";
const Char_t* cut_ccpi0_heg_psi2s2 = "!( (rm_xs==1001 || rm_xs==1010 || rm_xs==1101 || rm_xs==1110) && ((rm_l==1 && -0.10<cc_mheg -3.68609  && cc_mheg -3.68609 <0.05) || (rm_l==0 && -0.10<cc_mheg -3.68609  && cc_mheg -3.68609 <0.05)) )";
const Char_t* cut_ccpi0_heg_psi2s3 = "!( (rm_xs==1001 || rm_xs==1010 || rm_xs==1101 || rm_xs==1110) && ((rm_l==1 && -0.05<cc_mheg -3.68609  && cc_mheg -3.68609 <0.05) || (rm_l==0 && -0.05<cc_mheg -3.68609  && cc_mheg -3.68609 <0.05)) )";
*/
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

const Int_t nset = 21;
const Char_t* set_name[nset] = {"A","B","C","D","E","F","G",
				"H","I","J","K","L","M","N",
				"O","P","Q","R","S","T","U"
};
const Double_t sigmc_amount = 540.8;

//const Char_t* indir_gen = "~/ewp/ana/data/sigmc/hbk5/right/hbk_cut2/";
//const Char_t* indir_gen = "~/ewp/ana/data/sigmc/hbk5/right/hbk_cut2_c10flip/"; // c10-flip


namespace sig_bkg{
  const Char_t* indir[2] = {
    "~/ewp/ana/data/gmc/hbk5/right/hbk_cut2/",
    "~/ewp/ana/data/sigmc/hbk5/right/hbk_cut2/",
    //"~/ewp/ana/data/gmc/hbk5/right/hbk_calib_cut2/",
    //"~/ewp/ana/data/sigmc/hbk5/right/hbk_calib_cut2/",
    //"~/ewp/ana/ana_xsll_comb/NB_lep_calib/hbk/hbk_orgksfw_vtxcl_fmiss1_bb_widebcs/",       // bkg
    //"~/ewp/ana/ana_xsll_comb/NB_lep_calib/hbk/hbk_orgksfw_vtxcl_fmiss1_bb_widebcs_merge/", // sig
  };
  const Char_t* tail = "*.root";
}

namespace sigmc{
  //const Char_t* indir     = "~/ewp/ana/data/sigmc/hbk6/right/hbk_calib_cut2/";
  //const Char_t* indir     = "~/ewp/ana/ana_xsll_comb/hbk_bcs/sig/sig1_tot_522_merge/";
  //const Char_t* indir     = "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut2/sig_522/";
  //const Char_t* indir     = "~/ewp/ana/ana_xsll_comb/hbk_afb/hbk_cut2/sig_522/";
  const Char_t* indir     = "~/ewp/ana/ana_xsll_comb/tmp_hbk_forbcs_recmode/hbk_cut2/sig_522_lrnb_bb_bcs_merge_lrnb/";
  const Char_t* tail      = "*.root";
  const Char_t* indir_gen = "~/ewp/ana/data/sigmc/hbk5/right/hbk_org/";
}

namespace sigmc_afterbgsup{
  const Char_t* indir     = "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut2/sig_522/";
  const Char_t* tail      = "*.root";
  const Char_t* indir_gen = "~/ewp/ana/data/sigmc/hbk5/right/hbk_org/";
}

namespace sigmc_afterbgsup_c10sym{
  const Char_t* indir     = "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut2/sig_c10sym_522/";
  //const Char_t* indir     = "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_ff6/sig_522_ff6/";
  const Char_t* tail      = "*.root"; // not be used, should not change
  const Char_t* indir_gen = "~/ewp/ana/data/sigmc/hbk5/right/hbk_org_c10sym/"; // c10-symmetry
  //const Char_t* indir_gen = "~/ewp/ana/data/sigmc_ff/hbk6/right/hbk_org_ff6_c10sym/"; // c10-symmetry
}

namespace sigmc_afterbgsup_c10sym_cut3{
  const Char_t* indir     = "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut3/sig_c10sym_522/";
  const Char_t* tail      = "*.root"; // not be used, should not change
  const Char_t* indir_gen = "~/ewp/ana/data/sigmc/hbk5/right/hbk_org_c10sym/"; // c10-symmetry
}

namespace sigmc_afterbgsup_c10sym_all{
  const Char_t* indir     = "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut2/sig_c10sym_522_all/";
  //const Char_t* indir     = "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_lambdaone/sig_522_lambdaone429/";
  const Char_t* tail      = "*.root"; // not be used, should not change
  const Char_t* indir_gen = "~/ewp/ana/data/sigmc/hbk5/right/hbk_org_c10sym_all_tmp/"; // c10-symmetry ( "Xs_m and m_ll cut" is applied to "the hbk_org_c10sym_all" for data suppression )
  //const Char_t* indir_gen = "~/ewp/ana/data/sigmc_lambdaone/hbk6/right/hbk_org_lambdaone429/"; // c10-symmetry ( "Xs_m and m_ll cut" is applied to "the hbk_org_c10sym_all" for data suppression )
}

namespace sigmc_afterbgsup_cut6_c10sym{
  const Char_t* indir     = "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut6/sig_522_c10sym/";
  const Char_t* tail      = "*.root"; // not be used, should not change
  const Char_t* indir_gen = "~/ewp/ana/data/sigmc/hbk5/right/hbk_org_c10sym/"; // c10-symmetry ( "Xs_m and m_ll cut" is applied to "the hbk_org_c10sym_all" for data suppression )
}

namespace sigmc_afterbgsup_fm200{
  const Char_t* indir     = "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_fm/sig_522_fm200/"; // c10-symmetry (including c10flip)
  const Char_t* tail      = "*.root"; // not be used, should not change
  const Char_t* indir_gen = "~/ewp/ana/data/sigmc_fm/hbk6/right/hbk_org_fm200/"; // c10-symmetry (including c10flip)
}

namespace sigmc_afterbgsup_fm480{
  const Char_t* indir     = "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_fm/sig_522_fm480/"; // c10-symmetry (including c10flip)
  const Char_t* tail      = "*.root"; // not be used, should not change
  const Char_t* indir_gen = "~/ewp/ana/data/sigmc_fm/hbk6/right/hbk_org_fm480/"; // c10-symmetry (including c10flip)
}

namespace sigmc_afterbgsup_mb465{
  const Char_t* indir     = "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_mb/sig_522_mb465/"; // c10-symmetry (including c10flip)
  const Char_t* tail      = "*.root"; // not be used, should not change
  const Char_t* indir_gen = "~/ewp/ana/data/sigmc_mb/hbk6/right/hbk_org_mb465/"; // c10-symmetry (including c10flip)
}

namespace sigmc_afterbgsup_mb495{
  const Char_t* indir     = "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_mb/sig_522_mb495/"; // c10-symmetry (including c10flip)
  const Char_t* tail      = "*.root"; // not be used, should not change
  const Char_t* indir_gen = "~/ewp/ana/data/sigmc_mb/hbk6/right/hbk_org_mb495/"; // c10-symmetry (including c10flip)
}

namespace sigmc_calib_cut2{
  const Char_t* indir     = "~/ewp/ana/data/sigmc/hbk5/right/hbk_calib_cut2/";
  const Char_t* tail      = "*.root"; // not be used, should not change
  const Char_t* indir_gen = "~/ewp/ana/data/sigmc/hbk5/right/hbk_org/";
}

namespace sigmc_org{
  const Char_t* indir     = "~/ewp/ana/data/sigmc/hbk5/right/hbk_org/";
  const Char_t* tail      = "*.root"; // not be used, should not change
  const Char_t* indir_gen = "~/ewp/ana/data/sigmc/hbk5/right/hbk_org/";
}

namespace sigmc_org_c10flip{
  const Char_t* indir     = "~/ewp/ana/data/sigmc/hbk5/right/hbk_org_c10flip/";
  const Char_t* tail      = "*.root"; // not be used, should not change
  const Char_t* indir_gen = "~/ewp/ana/data/sigmc/hbk5/right/hbk_org_c10flip/"; // c10-flip
}

namespace sigmc_org_c10sym{
  const Char_t* indir     = "~/ewp/ana/data/sigmc/hbk5/right/hbk_org_c10sym/";
  const Char_t* tail      = "*.root"; // not be used, should not change
  const Char_t* indir_gen = "~/ewp/ana/data/sigmc/hbk5/right/hbk_org_c10sym/";
}

namespace sigmc_cut2{
  const Char_t* indir     = "~/ewp/ana/data/sigmc/hbk5/right/hbk_cut2/";
  const Char_t* tail      = "*.root"; // not be used, should not change
  const Char_t* indir_gen = "~/ewp/ana/data/sigmc/hbk5/right/hbk_org/";
}

namespace sigmc_cut2_c10flip{
  const Char_t* indir     = "~/ewp/ana/data/sigmc/hbk5/right/hbk_cut2_c10flip/";
  const Char_t* tail      = "*.root"; // not be used, should not change
  const Char_t* indir_gen = "~/ewp/ana/data/sigmc/hbk5/right/hbk_org_c10flip/"; // c10-flip
}

namespace sigmc_cut2_c10sym{
  const Char_t* indir     = "~/ewp/ana/data/sigmc/hbk5/right/hbk_cut2_c10sym/";
  const Char_t* tail      = "*.root"; // not be used, should not change
  const Char_t* indir_gen = "~/ewp/ana/data/sigmc/hbk5/right/hbk_org_c10sym/"; // c10-symmetry
}

namespace sigmc_cut4{
  const Char_t* indir     = "~/ewp/ana/data/sigmc/hbk5/right/hbk_cut4/";
  const Char_t* tail      = "*.root"; // not be used, should not change
  const Char_t* indir_gen = "~/ewp/ana/data/sigmc/hbk5/right/hbk_org/";
}

namespace sigmc_cut4_c10flip{
  const Char_t* indir     = "~/ewp/ana/data/sigmc/hbk5/right/hbk_cut4_c10flip/";
  const Char_t* tail      = "*.root"; // not be used, should not change
  const Char_t* indir_gen = "~/ewp/ana/data/sigmc/hbk5/right/hbk_org_c10flip/"; // c10-flip

}

namespace sigmc_cut4_c10sym{
  const Char_t* indir     = "~/ewp/ana/data/sigmc/hbk5/right/hbk_cut4_c10sym/";
  const Char_t* tail      = "*.root"; // not be used, should not change
  const Char_t* indir_gen = "~/ewp/ana/data/sigmc/hbk5/right/hbk_org_c10sym/"; // c10-symmetry
}

namespace gmc{
  //const Char_t* indir = "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut2/bkg_emu_522/";
  const Char_t* indir = "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut2/bkg_522/";
  //const Char_t* indir = "~/ewp/ana/ana_xsll_comb/hbk_afb/hbk_cut2/bkg_522/";
  const Char_t* tail  = "*.root";
}

namespace gmc_emu{
  const Char_t* indir[2] = {
    "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut2/bkg_522/",     // ee or mm
    "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut2/bkg_emu_522/", // emu
  };
  const Char_t* tail  = "*.root";
}

namespace gmc_rd{
  const Char_t* indir[2] = {
    "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut2/bkg_522/", // gmc
    "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut2/rd_522/",  // rd
  };
  const Char_t* tail = "*.root";
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

namespace sig_gmc_rd{
  const Char_t* indir[3] = {
    //"~/ewp/ana/ana_xsll_comb/NB_lep_calib/hbk/hbk_orgksfw_vtxcl_fmiss1_bb_522_190_bcs/",     // gmc
    //"~/ewp/ana/ana_xsll_comb/NB_lep_calib/hbk/hbk_orgksfw_vtxcl_fmiss1_bb_522_190_bcs_merge/",     // sigmc
    //"~/ewp/ana/ana_xsll_comb/hbk_bcs/bkg_lrnb/bkg1_bb_522/",     // gmc
    //"~/ewp/ana/ana_xsll_comb/hbk_bcs/sig_lrnb/sig1_bb_522_merge/",     // sigmc
    "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut2/bkg_522/",     // gmc
    "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut2/sig_522/",     // sigmc
    "~/ewp/ana/ana_xsll_comb/NB_lep/hbk_rd/hbk_orgksfw_vtxcl_fmiss1/",  // rd
  };
  const Char_t* tail = "*.root";
}

namespace sig_gmc_rd_calib{
  const Char_t* indir[3] = {
    "~/ewp/ana/data/gmc/hbk5/right/hbk_calib/all/", // gmc
    "~/ewp/ana/data/sigmc/hbk5/right/hbk_calib/",   // sigmc 
    "~/ewp/ana/data/gmc/hbk5/right/hbk_calib/all/", // rd
  };
  const Char_t* tail = "*.root";
}

namespace sig_gmc_rd_cut2{
  const Char_t* indir[3] = {
    "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut2/bkg_522/", // gmc
    "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut2/sig_522/", // sigmc 
    //"~/ewp/ana/ana_xsll_comb/hbk_bcs/bkg_lrnb/bkg1_bb_522/",       // gmc   // LR
    //"~/ewp/ana/ana_xsll_comb/hbk_bcs/sig_lrnb/sig1_bb_522_merge/", // sigmc // LR
    //"~/ewp/ana/ana_xsll_comb/hbk_bcs/bkg/bkg1_bb_522/",       // gmc   // LR
    //"~/ewp/ana/ana_xsll_comb/hbk_bcs/sig/sig1_bb_522_merge/", // sigmc // LR
    "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut2/rd_522/",  // rd
  };
  const Char_t* tail = "*.root";
}

namespace sig_gmc_rd_cut2_all{
  const Char_t* indir[3] = {
    "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut2/bkg_522/",       // gmc
    "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut2/sig_522_all/", // sigmc 
    "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut2/rd_522/",        // rd
  };
  const Char_t* tail = "*.root";
}

namespace sig_gmc_rd_cut2_c10sym{
  const Char_t* indir[3] = {
    "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut2/bkg_522/", // gmc
    "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut2/sig_c10sym_522/", // sigmc 
    //"~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut2/sig_c10sym_522_all/", // sigmc 
    "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut2/rd_522/",  // rd
  };
  const Char_t* tail = "*.root";
}

namespace sig_gmc_rd_cut2_cc{
  const Char_t* indir[3] = {
    "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut2/cc_522/",  // cc mc
    "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut2/sig_522/", // sigmc 
    "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut2/rd_522/",  // rd
  };
  const Char_t* tail = "*.root";
}

namespace sig_gmc_rd_cut3{
  const Char_t* indir[3] = {
    "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut3/bkg_522/", // gmc
    "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut3/sig_522/", // sigmc 
    "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut3/rd_522/",  // rd
  };
  const Char_t* tail = "*.root";
}

namespace sig_gmc_rd_emu{
  const Char_t* indir[3] = {
    "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut2/bkg_emu/", // gmc
    "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut2/sig/",     // sigmc(dummy data)
    "~/ewp/ana/ana_xsll_comb/hbk_afb_calib/hbk_cut2/rd_emu/",  // rd
  };
  const Char_t* tail = "*.root";
}

namespace sig_gmc_rd_cut2_beforebgsup{
  const Char_t* indir[3] = {
    "~/ewp/ana/ana_xsll_comb/NB_lep_calib/hbk/hbk_orgksfw_vtxcl_fmiss1/",     // gmc
    "~/ewp/ana/ana_xsll_comb/NB_lep_calib/hbk/hbk_orgksfw_vtxcl_fmiss1/",     // sigmc 
    "~/ewp/ana/ana_xsll_comb/NB_lep_calib/hbk_rd/hbk_orgksfw_vtxcl_fmiss1/",  // rd
    //"~/ewp/ana/data/gmc/hbk5/right/hbk_calib_cut2/",     // gmc
    //"~/ewp/ana/data/sigmc/hbk5/right/hbk_calib_cut2/",   // sigmc
    //"~/ewp/ana/data/gmc/hbk5/right/hbk_calib_cut2/",     // rd
  };
  const Char_t* tail = "*.root";
}

namespace sig_gmc_rd_cut3_beforebgsup{
  const Char_t* indir[3] = {
    "~/ewp/ana/ana_xsll_comb/NB_lep_calib/hbk_cc/hbk_orgksfw_vtxcl_fmiss1/",     // gmc
    "~/ewp/ana/ana_xsll_comb/NB_lep_calib/hbk/hbk_orgksfw_vtxcl_fmiss1/",        // sigmc
    "~/ewp/ana/ana_xsll_comb/NB_lep_calib/hbk_rd_cc/hbk_orgksfw_vtxcl_fmiss1/",  // rd
    //"~/ewp/ana/data/gmc/hbk5/right/hbk_calib_cut3/",     // gmc
    //"~/ewp/ana/data/sigmc/hbk5/right/hbk_calib_cut3/",   // sigmc
    //"~/ewp/ana/data/gmc/hbk5/right/hbk_calib_cut3/",     // rd
  };
  const Char_t* tail = "*.root";
}

namespace sig_gmc_rd_emu_beforebgsup{
  const Char_t* indir[3] = {
    "~/ewp/ana/ana_xsll_comb/NB_lep_calib/hbk_emu/hbk_orgksfw_vtxcl_fmiss1/",    // gmc
    "~/ewp/ana/ana_xsll_comb/NB_lep_calib/hbk/hbk_orgksfw_vtxcl_fmiss1/",        // sigmc(dummy data)
    "~/ewp/ana/ana_xsll_comb/NB_lep_calib/hbk_rd_emu/hbk_orgksfw_vtxcl_fmiss1/", // rd
    //"~/ewp/ana/data/gmc/hbk5/emu/hbk_calib_cut2/",     // gmc
    //"~/ewp/ana/data/sigmc/hbk5/right/hbk_calib_cut2/",   // sigmc(dummy data)
    //"~/ewp/ana/data/gmc/hbk5/emu/hbk_calib_cut2/",     // rd
  };
  const Char_t* tail = "*.root";
}

namespace sig_gmc_rd_cut2_double{
  const Char_t* indir[3] = {
    "~/ewp/ana/data/gmc_fl/hbk5/right/hbk_double_nb/hbk_orgksfw_vtxcl_fmiss1_weight/", // gmc
    "~/ewp/ana/data/gmc_fl/hbk5/right/hbk_double_nb/hbk_orgksfw_vtxcl_fmiss1_weight/", // sigmc(dummy data)
    "~/ewp/ana/data/gmc_fl/hbk5/right/hbk_double_nb/hbk_orgksfw_vtxcl_fmiss1_weight/", // rd
  };
  const Char_t* tail = "*.root";
}

namespace sig_gmc_rd_cut4_swap{
  const Char_t* indir[3] = {
    "~/ewp/ana/data/gmc2/hbk5/right/hbk_cut4_nb/hbk_orgksfw_vtxcl_fmiss1_weight/", // gmc
    "~/ewp/ana/data/gmc2/hbk5/right/hbk_cut4_nb/hbk_orgksfw_vtxcl_fmiss1_weight/", // sigmc(dummy data)
    "~/ewp/ana/data/gmc2/hbk5/right/hbk_cut4_nb/hbk_orgksfw_vtxcl_fmiss1_weight/", // rd
  };
  const Char_t* tail = "*.root";
}

namespace sig_gmc_rd_comp{
  const Char_t* indir[6] = {
    "~/ewp/ana/ana_xsll_comb/NB_lep_calib/hbk/hbk_orgksfw_vtxcl_fmiss1/",        // gmc
    "~/ewp/ana/ana_xsll_comb/NB_lep_calib/hbk_emu/hbk_orgksfw_vtxcl_fmiss1/",    // gmc(emu)
    "~/ewp/ana/ana_xsll_comb/NB_lep_calib/hbk_cc/hbk_orgksfw_vtxcl_fmiss1/",     // gmc(jpsi[cut3])
    "~/ewp/ana/ana_xsll_comb/NB_lep_calib/hbk/hbk_orgksfw_vtxcl_fmiss1/",        // sigmc
    "~/ewp/ana/ana_xsll_comb/NB_lep_calib/hbk_rd_emu/hbk_orgksfw_vtxcl_fmiss1/", // rd(emu)
    "~/ewp/ana/ana_xsll_comb/NB_lep_calib/hbk_rd_cc/hbk_orgksfw_vtxcl_fmiss1/",  // rd(jpsi[cut3])
  };
  const Char_t* tail = "*.root";
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
const Int_t col_fil[30] = {4,9,3,2,15,12,1,16,5,13,6,11,7,8,10,14, 17,18,19,20,21,22,23,24,25,26,27,28,29,30};
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

namespace Mbc_deltaE{
  const Char_t*  tname    = "h511";
  const Char_t*  axis     = "de:Mbc";
  const Int_t    ybin     =   300;
  const Double_t ymin     = -0.20;
  const Double_t ymax     =  0.20;
  const Char_t*  ylabel   = "#DeltaE [GeV]";
  const Double_t offset   =   0.0;
  const Int_t    xbin     =   300;
  const Double_t xmin     =  5.22;
  const Double_t xmax     =  5.30;
  const Char_t*  xlabel   = "M_{bc} [GeV]";
}

namespace Mbc{
  const Char_t*  tname    = "h511";
  const Char_t*  axis     = "Mbc";
  const Double_t offset   =  0.0;
  const Int_t    xbin     =   80; // bin width = 1 [MeV]
  const Double_t xmin     = 5.22;
  const Double_t xmax     = 5.30;
  const Double_t xmin_fit = xmin;
  const Double_t xmax_fit = xmax;
  const Char_t*  xlabel   = "M_{bc} [GeV]";
}

namespace Mbc_comb{
  const Char_t*  tname    = "h511";
  const Char_t*  axis     = "Mbc";
  const Double_t offset   =  0.0;
  const Int_t    xbin     =   40;
  const Double_t xmin     = 5.22; // bin width = 2 [MeV]
  const Double_t xmax     = 5.30;
  const Double_t xmin_fit = xmin;
  const Double_t xmax_fit = xmax;
  const Char_t*  xlabel   = "M_{bc} [GeV]";
}

namespace Mbc_bkg{
  const Char_t*  tname    = "h511";
  const Char_t*  axis     = "Mbc";
  const Double_t offset   =  0.0;
  const Int_t    xbin     =   32;
  const Double_t xmin     = 5.22; // bin width = 2.5 [MeV]
  const Double_t xmax     = 5.30;
  const Double_t xmin_fit = xmin;
  const Double_t xmax_fit = xmax;
  const Char_t*  xlabel   = "M_{bc} [GeV]";
}

namespace Mbc_wide{
  const Char_t*  tname    = "h511";
  const Char_t*  axis     = "Mbc";
  const Double_t offset   =  0.0;
  const Int_t    xbin     =  100; // bin width = 1 [MeV]
  const Double_t xmin     = 5.20;
  const Double_t xmax     = 5.30;
  const Double_t xmin_fit = xmin;
  const Double_t xmax_fit = xmax;
  const Char_t*  xlabel   = "M_{bc} [GeV]";
}

namespace Mbc_comb_wide{
  const Char_t*  tname    = "h511";
  const Char_t*  axis     = "Mbc";
  const Double_t offset   =  0.0;
  const Int_t    xbin     =   50;
  const Double_t xmin     = 5.20; // bin width = 2 [MeV]
  const Double_t xmax     = 5.30;
  const Double_t xmin_fit = xmin;
  const Double_t xmax_fit = xmax;
  const Char_t*  xlabel   = "M_{bc} [GeV]";
}

namespace Mbc_bkg_wide{
  const Char_t*  tname    = "h511";
  const Char_t*  axis     = "Mbc";
  const Double_t offset   =  0.0;
  const Int_t    xbin     =   40;
  const Double_t xmin     = 5.20; // bin width = 2.5 [MeV]
  const Double_t xmax     = 5.30;
  const Double_t xmin_fit = xmin;
  const Double_t xmax_fit = xmax;
  const Char_t*  xlabel   = "M_{bc} [GeV]";
}

namespace deltaE{
  const Char_t*  tname    = "h511";
  const Char_t*  axis     = "de";
  const Double_t offset   = 0.0;
  const Int_t    xbin     = 120; 
  const Double_t xmin     = -0.12;
  const Double_t xmax     =  0.12;
  const Double_t xmin_fit = xmin;
  const Double_t xmax_fit = xmax;
  const Char_t*  xlabel   = "#DeltaE [GeV]";
}

namespace deltaE_wide{
  const Char_t*  tname    = "h511";
  const Char_t*  axis     = "de";
  const Double_t offset   = 0.0;
  const Int_t    xbin     =  80; 
  const Double_t xmin     = -0.20;
  const Double_t xmax     =  0.20;
  const Double_t xmin_fit = -0.10;
  const Double_t xmax_fit =  0.10;
  const Char_t*  xlabel   = "#DeltaE [GeV]";
}

namespace deltaE_calib{
  const Char_t*  tname    = "h511";
  const Char_t*  axis     = "de";
  const Double_t offset   = 0.0;
  const Int_t    xbin     =  200; 
  const Double_t xmin     = -0.20;
  const Double_t xmax     =  0.20;
  const Double_t xmin_fit = -0.10;
  const Double_t xmax_fit =  0.10;
  const Char_t*  xlabel   = "#DeltaE [GeV]";
}

namespace dzll3d{
  const Char_t*  tname    = "h511";
  const Char_t*  axis     = "10000*dzll3d";
  const Char_t*  axisrd   = "10000*dzll3d";
  const Char_t*  fname    = "dzll3d";
  const Double_t offset   = 0.0;
  const Int_t    xbin     =  85; 
  const Double_t xmin     =   0;
  const Double_t xmax     = 190;
  const Double_t xmin_fit =   0;
  const Double_t xmax_fit = 150;
  const Char_t*  xlabel   = "|#Deltaz(l^{+}l^{-})| [#mum]";
}

namespace dzll3dorg{
  const Char_t*  tname    = "h511";
  const Char_t*  axis     = "10000*dzll3dorg";
  const Char_t*  axisrd   = "10000*dzll3d";
  const Char_t*  fname    = "dzll3dorg";
  const Double_t offset   = 0.0;
  const Int_t    xbin     =  85; 
  const Double_t xmin     =   0;
  const Double_t xmax     = 190;
  const Double_t xmin_fit =   0;
  const Double_t xmax_fit = 190;
  const Char_t*  xlabel   = "|#Deltaz(l^{+}l^{-})| [#mum]";
}

namespace  bvtxchi2{
  const Char_t*  tname    = "h511";
  const Char_t*  axis     = "kfbchi/kfbdgf";
  const Char_t*  fname    = "kfbchi";
  const Double_t offset   =  0.0;
  const Int_t    xbin     =  200; 
  const Double_t xmin     =  0.0;
  const Double_t xmax     = 10.0;
  const Double_t xmin_fit =  0.0;
  const Double_t xmax_fit = 10.0;
  const Char_t*  xlabel   = "#chi^{2}/NDF";
}

namespace  bvtxchi2org{
  const Char_t*  tname    = "h511";
  const Char_t*  axis     = "kfbchiorg/kfbdgf";
  const Char_t*  fname    = "kfbchi";
  const Double_t offset   =  0.0;
  const Int_t    xbin     =  200; 
  const Double_t xmin     =  0.0;
  const Double_t xmax     = 10.0;
  const Double_t xmin_fit =  0.0;
  const Double_t xmax_fit = 10.0;
  const Char_t*  xlabel   = "#chi^{2}/NDF";
}

namespace M_Xs{
  const Char_t*  tname    = "h511";
  const Char_t*  axis     = "xs_m";
  const Double_t offset   = 0.0;
  const Int_t    xbin     = 100; 
  const Double_t xmin     = 0.0;
  const Double_t xmax     = 3.5;
  const Double_t xmin_fit = xmin;
  const Double_t xmax_fit = xmax;
  const Char_t*  xlabel   = "M_{Xs} [GeV]";
}
namespace M_ll{
  const Char_t*  tname    = "h511";
  const Char_t*  axis     = "cc_m";
  const Double_t offset   = 0.0;
  const Int_t    xbin     = 100; 
  const Double_t xmin     = 0.0;
  const Double_t xmax     = 5.0;
  const Double_t xmin_fit = xmin;
  const Double_t xmax_fit = xmax;
  const Char_t*  xlabel   = "M_{ll} [GeV]";
}

namespace bgsup{
  const Int_t    n_axis     = 16;
  const Char_t*  tname      = "h511";
  const Char_t*  axis[n_axis]     = {         "evis1",   "abs(mmiss1)",              "deroe1", "kfbchi/kfbdgf",                 "10000*dzll3d",         "bccm",      "fmiss1",       "ksfw_tot",       "ksfw_qq",       "ksfw_bb",      "ksfw_tot_lep1",      "ksfw_qq_lep1",      "ksfw_bb_lep1",     "ksfw_tot_lep0",      "ksfw_qq_lep0",      "ksfw_bb_lep0" };
  const Char_t*  fname[n_axis]    = {         "evis1",     "absmmiss1",              "deroe1",       "kfbchi2",                       "dzll3d",         "bccm",      "Fmiss1",       "ksfw_tot",       "ksfw_qq",       "ksfw_bb",      "ksfw_tot_lep1",      "ksfw_qq_lep1",      "ksfw_bb_lep1",     "ksfw_tot_lep0",      "ksfw_qq_lep0",      "ksfw_bb_lep0" }; // name for save-file
  const Double_t offset[n_axis]   = {            0.0,              0.0,                   0.0,             0.0,                            0.0,            0.0,           0.0,              0.0,             0.0,             0.0,                  0.0,                 0.0,                 0.0,                 0.0,                 0.0,                 0.0 };
  const Int_t    xbin[n_axis]     = {            100,              100,                   100,             200,                            100,            100,           100,              200,             200,             200,                  200,                 200,                 200,                 200,                 200,                 200 };
  const Double_t xmin[n_axis]     = {            0.0,              0.0,                  -5.0,             0.0,                            0.0,           -1.0,          -3.5,             -5.0,            -5.0,            -5.0,                 -5.0,                -5.0,                -5.0,                -5.0,                -5.0,                -5.0 };
  const Double_t xmax[n_axis]     = {           14.0,             10.0,                   2.0,            10.0,                            150,            1.0,           1.5,             15.0,            15.0,            15.0,                 15.0,                15.0,                15.0,                15.0,                15.0,                15.0 };
  const Double_t xmin_fit[n_axis] = {        xmin[0],          xmin[1],               xmin[2],         xmin[3],                        xmin[4],        xmin[5],       xmin[6],          xmin[7],         xmin[8],         xmin[9],             xmin[10],            xmin[11],            xmin[12],            xmin[13],            xmin[14],            xmin[15] };
  const Double_t xmax_fit[n_axis] = {        xmax[0],          xmax[1],               xmax[2],         xmax[3],                        xmax[4],        xmax[5],       xmax[6],          xmax[7],         xmax[8],         xmax[9],             xmax[10],            xmax[11],            xmax[12],            xmax[13],            xmax[14],            xmax[15] };
  const Char_t*  xlabel[n_axis]   = {"E_{vis} [GeV]", "M_{miss} [GeV]", "#DeltaE^{roe} [GeV]",  "#chi^{2}/NDF", "|#Deltaz(l^{+}l^{-})| [#mum]", "cos#theta_{B}",   "F_{miss}", "F_{ksfw}^{tot}", "F_{ksfw}^{qq}", "F_{ksfw}^{bb}", "F_{ksfw}^{tot}(ee)", "F_{ksfw}^{qq}(ee)", "F_{ksfw}^{bb}(ee)","F_{ksfw}^{tot}(mm)", "F_{ksfw}^{qq}(mm)", "F_{ksfw}^{bb}(mm)" };
}
namespace Evis{
  const Char_t*  tname    = "h511";
  const Char_t*  axis     = "evis"; // evis1ks
  const Double_t offset   = 0.0;
  const Int_t    xbin     = 100; 
  const Double_t xmin     =  0.0;
  const Double_t xmax     = 20.0;
  const Double_t xmin_fit = xmin;
  const Double_t xmax_fit = xmax;
  const Char_t*  xlabel   = "E_{vis} [GeV]";
}

namespace Evis_multi{
  const Char_t*  tname    = "h511";
  const Char_t*  axis[6]  = {"evis1","evis1ks","evis2","evis2ks","evis3","evis3ks"};
  const Double_t offset   = 0.0;
  const Int_t    xbin     = 100; 
  const Double_t xmin     =  0.0;
  const Double_t xmax     = 20.0;
  const Double_t xmin_fit = xmin;
  const Double_t xmax_fit = xmax;
  const Char_t*  xlabel   = "E_{vis} [GeV]";
}

namespace Mmiss{
  const Char_t*  tname    = "h511";
  const Char_t*  axis     = "abs(mmiss)";
  const Double_t offset   = 0.0;
  const Int_t    xbin     =  80; 
  const Double_t xmin     =  0.0;
  const Double_t xmax     =  8.0;
  const Double_t xmin_fit = xmin;
  const Double_t xmax_fit = xmax;
  const Char_t*  xlabel   = "M_{miss} [GeV]";
}

namespace Mmiss_multi{
  const Char_t*  tname    = "h511";
  const Char_t*  axis[6]  = {"abs(mmiss1)","abs(mmiss1ks)","abs(mmiss2)","abs(mmiss2ks)","abs(mmiss3)","abs(mmiss3ks)"};
  const Double_t offset   = 0.0;
  const Int_t    xbin     = 80; 
  const Double_t xmin     =  0.0;
  const Double_t xmax     = 8.0;
  const Double_t xmin_fit = xmin;
  const Double_t xmax_fit = xmax;
  const Char_t*  xlabel   = "M_{miss} [GeV]";
}

namespace Mmiss_multi_sign{
  const Char_t*  tname    = "h511";
  const Char_t*  axis[6]  = {"mmiss1","mmiss1ks","mmiss2","mmiss2ks","mmiss3","mmiss3ks"};
  const Double_t offset   = 0.0;
  const Int_t    xbin     =  80; 
  const Double_t xmin     = -8.0;
  const Double_t xmax     =  8.0;
  const Double_t xmin_fit = xmin;
  const Double_t xmax_fit = xmax;
  const Char_t*  xlabel   = "M_{miss} [GeV]";
}

namespace Deroe{
  const Char_t*  tname    = "h511";
  const Char_t*  axis     = "deroe";
  const Double_t offset   = 0.0;
  const Int_t    xbin     = 100; 
  const Double_t xmin     = -5.0;
  const Double_t xmax     =  5.0;
  const Double_t xmin_fit = xmin;
  const Double_t xmax_fit = xmax;
  const Char_t*  xlabel   = "#Delta E^{roe} [GeV]";
}

namespace Deroe_multi{
  const Char_t*  tname    = "h511";
  const Char_t*  axis[6]  = {"deroe1","deroe1ks","deroe2","deroe2ks","deroe3","deroe3ks"};
  const Double_t offset   = 0.0;
  const Int_t    xbin     = 100; 
  const Double_t xmin     = -5.0;
  const Double_t xmax     =  5.0;
  const Double_t xmin_fit = xmin;
  const Double_t xmax_fit = xmax;
  const Char_t*  xlabel   = "#Delta E^{roe} [GeV]";
}

namespace ksfw_lr_tot{
  const Char_t*  tname    = "h511";
  const Char_t*  axis     = "ksfwsig_tot/(ksfwsig_tot+ksfwbkg_tot)";
  const Double_t offset   = 0.0;
  const Int_t    xbin     = 100;
  const Double_t xmin     = 0.0;
  const Double_t xmax     = 1.0;
  const Double_t xmin_fit = xmin;
  const Double_t xmax_fit = xmax;
  const Char_t*  fname    = "ksfw_tot";
  const Char_t*  xlabel   = "LR of F_{ksfw}^{tot}";
}

namespace ksfw_var_tot{
  const Char_t*  tname    = "h511";
  const Char_t*  axis     = "ksfw_tot";
  const Double_t offset   = 0.0;
  const Int_t    xbin     = 100;
  const Double_t xmin     = -10.0;
  const Double_t xmax     =  10.0;
  const Double_t xmin_fit = xmin;
  const Double_t xmax_fit = xmax;
  const Char_t*  fname    = "ksfw_tot";
  const Char_t*  xlabel   = "ksfw^{tot}";
}

namespace ksfw{
  const Int_t    n_axis = 18; 
  const Char_t*  tname  = "h511";
  const Char_t*  axis[n_axis]   = { "k0mm2", "k0et", "k0hso00", "k0hso10", "k0hso20", "k0hso01", "k0hso02", "k0hso12", "k0hso22", "k0hso03", "k0hso04", "k0hso14", "k0hso24", "k0hoo0", "k0hoo1", "k0hoo2", "k0hoo3", "k0hoo4" };
  const Char_t*  fname[n_axis]  = { "k0mm2", "k0et", "k0hso00", "k0hso10", "k0hso20", "k0hso01", "k0hso02", "k0hso12", "k0hso22", "k0hso03", "k0hso04", "k0hso14", "k0hso24", "k0hoo0", "k0hoo1", "k0hoo2", "k0hoo3", "k0hoo4" };
  const Double_t offset[n_axis] = {     0.0,    0.0,       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,      0.0,      0.0,      0.0,      0.0,      0.0 };
  const Int_t    xbin[n_axis]   = {      35,     40,        30,        30,        40,        25,        25,        35,        35,        30,        35,        40,        50,       40,       30,       40,       30,       30 };
  const Double_t xmin[n_axis]   = {    -4.0,    3.0,       0.0,       0.0,       0.0,      -0.5,      -0.4,      -0.3,      -0.3,      -0.3,      -0.3,     -0.20,      -0.2,      0.0,   -0.015,    -0.02,   -0.015,    -0.02 };
  const Double_t xmax[n_axis]   = {    10.0,   11.0,       1.5,       1.2,       0.8,       0.5,       0.6,       0.4,       0.4,       0.3,       0.4,      0.20,       0.3,      0.2,    0.015,     0.06,    0.015,     0.04 };
  const Char_t*  xlabel[n_axis] = { "k0mm2", "k0et", "k0hso00", "k0hso10", "k0hso20", "k0hso01", "k0hso02", "k0hso12", "k0hso22", "k0hso03", "k0hso04", "k0hso14", "k0hso24", "k0hoo0", "k0hoo1", "k0hoo2", "k0hoo3", "k0hoo4" };
}

namespace ksfw_fsp{ // use_final_state_particle
  const Int_t    n_axis = 18; 
  const Char_t*  tname  = "h511";
  const Char_t*  axis[n_axis]   = { "k1mm2", "k1et", "k1hso00", "k1hso10", "k1hso20", "k1hso01", "k1hso02", "k1hso12", "k1hso22", "k1hso03", "k1hso04", "k1hso14", "k1hso24", "k1hoo0", "k1hoo1", "k1hoo2", "k1hoo3", "k1hoo4" };
  const Char_t*  fname[n_axis]  = { "k1mm2", "k1et", "k1hso00", "k1hso10", "k1hso20", "k1hso01", "k1hso02", "k1hso12", "k1hso22", "k1hso03", "k1hso04", "k1hso14", "k1hso24", "k1hoo0", "k1hoo1", "k1hoo2", "k1hoo3", "k1hoo4" };
  const Double_t offset[n_axis] = {     0.0,    0.0,       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,      0.0,      0.0,      0.0,      0.0,      0.0 };
  const Int_t    xbin[n_axis]   = {      35,     40,        30,        30,        40,        25,        25,        35,        35,        30,        35,        40,        50,       50,       30,       40,       30,       30 };
  const Double_t xmin[n_axis]   = {    -4.0,    3.0,       0.0,       0.0,       0.0,      -0.5,      -0.4,      -0.3,      -0.3,      -0.3,      -0.3,     -0.20,      -0.2,      0.0,   -0.015,    -0.02,   -0.015,    -0.02 };
  const Double_t xmax[n_axis]   = {    10.0,   11.0,       1.5,       1.2,       0.8,       0.5,       0.6,       0.4,       0.4,       0.3,       0.4,      0.20,       0.3,      0.2,    0.015,     0.06,    0.015,     0.04 };
  const Char_t*  xlabel[n_axis] = { "k1mm2", "k1et", "k1hso00", "k1hso10", "k1hso20", "k1hso01", "k1hso02", "k1hso12", "k1hso22", "k1hso03", "k1hso04", "k1hso14", "k1hso24", "k1hoo0", "k1hoo1", "k1hoo2", "k1hoo3", "k1hoo4" };
}

namespace ksfw_comb{
  const Int_t    n_axis = 18; 
  const Char_t*  tname  = "h511";
  const Char_t*  axis[n_axis]   = { "k%dmm2","k%det","k%dhso00","k%dhso10","k%dhso20","k%dhso01","k%dhso02","k%dhso12","k%dhso22","k%dhso03","k%dhso04","k%dhso14","k%dhso24","k%dhoo0","k%dhoo1","k%dhoo2","k%dhoo3","k%dhoo4" };
  const Char_t*  fname[n_axis]  = { "k_mm2", "k_et", "k_hso00", "k_hso10", "k_hso20", "k_hso01", "k_hso02", "k_hso12", "k_hso22", "k_hso03", "k_hso04", "k_hso14", "k_hso24", "k_hoo0", "k_hoo1", "k_hoo2", "k_hoo3", "k_hoo4"  };
  const Double_t offset[n_axis] = {     0.0,    0.0,       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,      0.0,      0.0,      0.0,      0.0,      0.0  };
  const Int_t    xbin[n_axis]   = {      35,     40,        30,        30,        40,        25,        25,        35,        35,        30,        35,        40,        50,       50,       30,       40,       30,       30  };
  const Double_t xmin[n_axis]   = {    -4.0,    3.0,       0.0,       0.0,       0.0,      -0.5,      -0.4,      -0.3,      -0.3,      -0.3,      -0.3,     -0.20,      -0.2,      0.0,   -0.015,    -0.02,   -0.015,    -0.02  };
  const Double_t xmax[n_axis]   = {    10.0,   11.0,       1.5,       1.2,       0.8,       0.5,       0.6,       0.4,       0.4,       0.3,       0.4,      0.20,       0.3,      0.2,    0.015,     0.06,    0.015,     0.04  };
  const Char_t*  xlabel[n_axis] = { "k_mm2", "k_et", "k_hso00", "k_hso10", "k_hso20", "k_hso01", "k_hso02", "k_hso12", "k_hso22", "k_hso03", "k_hso04", "k_hso14", "k_hso24", "k_hoo0", "k_hoo1", "k_hoo2", "k_hoo3", "k_hoo4"  };
}

namespace fmiss{
  const Int_t    n_axis = 6;
  const Char_t*  tname  = "h511";
  const Char_t*  axis[n_axis]   = { "evis1",    "abs(mmiss1)",          "deroe1",    "mmiss1", "mmiss1*mmiss1", "abs(mmiss1)*mmiss1*mmiss1/mmiss1" };
  const Char_t*  fname[n_axis]  = { "evis1",      "absmmiss1",          "deroe1",    "mmiss1", "mmiss1_square",                "mmiss1_signsquare" };
  const Double_t offset[n_axis] = {     0.0,              0.0,               0.0,         0.0,             0.0,                                0.0 };
  const Int_t    xbin[n_axis]   = {     280,              250,               280,         240,             240,                                300 };
  const Double_t xmin[n_axis]   = {     4.0,              0.0,              -5.0,        -3.0,             0.0,                               -3.0 };
  const Double_t xmax[n_axis]   = {    13.0,              5.0,               2.0,         5.0,            12.0,                               12.0 };
  const Char_t*  xlabel[n_axis] = { "E_{vis1}", "|M_{miss1}|", "#Delta E^{ROE1}", "M_{miss1}",  "M_{mis1}^{2}",                  "+-M_{miss1}^{2}" };
}

namespace M_Xs_gen{
  const Char_t*  tname    = "h12";
  const Char_t*  axis     = "Xs_m";
  const Double_t offset   = 0.0;
  const Int_t    xbin     = 100; 
  const Double_t xmin     = 0.0;
  const Double_t xmax     = 3.5;
  const Double_t xmin_fit = xmin;
  const Double_t xmax_fit = xmax;
  const Char_t*  xlabel   = "M_{Xs} [GeV]";
}
namespace M_ll_gen{
  const Char_t*  tname    = "h12";
  const Char_t*  axis     = "llg_m";
  const Double_t offset   = 0.0;
  const Int_t    xbin     = 100; 
  const Double_t xmin     = 0.0;
  const Double_t xmax     = 5.0;
  const Double_t xmin_fit = xmin;
  const Double_t xmax_fit = xmax;
  const Char_t*  xlabel   = "M_{ll(#gamma)} [GeV]";
}


namespace q2_theta_uniform{
  const Char_t*  tname  = "h511";
  const Char_t*  axis   = "coslp:cc_m*cc_m";
  const Int_t    ybin   =    10;
  const Double_t ymin   = -1.00;
  const Double_t ymax   =  1.00;
  const Char_t*  ylabel = "cos#theta";
  const Int_t    xbin   =    25;
  const Double_t xmin   =   0.0;
  const Double_t xmax   =  25.0;
  const Char_t*  xlabel = "q^{2} [GeV^{2}]";
}

namespace q2_theta_uniform_gen{
  const Char_t*  tname  = "h12";
  const Char_t*  axis    = "coslp:llg_m*llg_m";
  const Int_t    ybin   =    30;
  const Double_t ymin   = -1.00;
  const Double_t ymax   =  1.00;
  const Char_t*  ylabel = "cos#theta";
  const Int_t    xbin   =    10;
  const Double_t xmin   =   0.0;
  const Double_t xmax   =  25.0;
  const Char_t*  xlabel = "q^{2} [GeV^{2}]";
}


namespace q2_theta_nonuniform{
  
  const Char_t*  tname   = "h511";
  const Char_t*  axis    = "coslp:cc_m*cc_m";
  const Int_t    ybin    =    10;
  const Double_t ymin    = -1.00;
  const Double_t ymax    =  1.00;
  const Char_t*  ylabel  = "cos#theta";
  const Double_t xmax    =  25.0;
  const Char_t*  xlabel  = "q^{2} [GeV^{2}]";
  
  const Int_t    xbin    =    9;
  
  const Double_t xbins[2][xbin+1] = { // [mm,ee]
    {0.0, 2.0, 4.3,
     (PDGmass::jpsi -0.25-0.0*0.15)*(PDGmass::jpsi -0.25-0.0*0.15), // ~ 8.10
     (PDGmass::jpsi +0.10+0.0*0.05)*(PDGmass::jpsi +0.10+0.0*0.05), // ~10.22
     (PDGmass::psi2s-0.15-0.0*0.10)*(PDGmass::psi2s-0.15-0.0*0.10), // ~12.50
     (PDGmass::psi2s+0.10                )*(PDGmass::psi2s+0.10                ), // ~14.33
     16.0, 19.0, 25.0},
    {0.0, 2.0, 4.3,
     (PDGmass::jpsi -0.25-1.0*0.15)*(PDGmass::jpsi -0.25-1.0*0.15), // ~ 7.27
     (PDGmass::jpsi +0.10+1.0*0.05)*(PDGmass::jpsi +0.10+1.0*0.05), // ~10.54
     (PDGmass::psi2s-0.15-1.0*0.10)*(PDGmass::psi2s-0.15-1.0*0.10), // ~11.81
     (PDGmass::psi2s+0.10                )*(PDGmass::psi2s+0.10                ), // ~14.33
     16.0, 19.0, 25.0},
  };
}  

namespace q2_theta_nonuniform_gen{
  
  const Char_t*  tname   = "h12";
  const Char_t*  axis    = "coslp:llg_m*llg_m";
  const Int_t    ybin    =    30;
  const Double_t ymin    = -1.00;
  const Double_t ymax    =  1.00;
  const Char_t*  ylabel  = "cos#theta";
  const Double_t offset   =  0.0;
  const Double_t xmin     =  0.0;
  const Double_t xmax    =  25.0;
  const Char_t*  xlabel  = "q^{2} [GeV^{2}]";
  const Int_t    xbin    =    9;
  const Double_t xbins[2][xbin+1] = { // [mm,ee]
    {0.0, 2.0, 4.3,
     (PDGmass::jpsi -0.25-0.0*0.15)*(PDGmass::jpsi -0.25-0.0*0.15), // ~ 8.10
     (PDGmass::jpsi +0.10+0.0*0.05)*(PDGmass::jpsi +0.10+0.0*0.05), // ~10.22
     (PDGmass::psi2s-0.15-0.0*0.10)*(PDGmass::psi2s-0.15-0.0*0.10), // ~12.50
     (PDGmass::psi2s+0.10                )*(PDGmass::psi2s+0.10                ), // ~14.33
     16.0, 19.0, 25.0},
    {0.0, 2.0, 4.3,
     (PDGmass::jpsi -0.25-1.0*0.15)*(PDGmass::jpsi -0.25-1.0*0.15), // ~ 7.27
     (PDGmass::jpsi +0.10+1.0*0.05)*(PDGmass::jpsi +0.10+1.0*0.05), // ~10.54
     (PDGmass::psi2s-0.15-1.0*0.10)*(PDGmass::psi2s-0.15-1.0*0.10), // ~11.81
     (PDGmass::psi2s+0.10                )*(PDGmass::psi2s+0.10                ), // ~14.33
     16.0, 19.0, 25.0},
  };

}  



namespace q2_theta_nonuniform_eff{
  const Char_t*  tname[2] = { "h511", "h12" };
  const Char_t*  axis[2]  = { "coslp:cc_m*cc_m", "coslp:llg_m*llg_m"};
  const Int_t    ybin     =    30;
  const Double_t ymin     = -1.00;
  const Double_t ymax     =  1.00;
  const Char_t*  ylabel   = "cos#theta";
  const Double_t offset   =   0.0;
  const Double_t xmin     =   0.0;
  const Double_t xmax     =  25.0;
  const Char_t*  xlabel   = "q^{2} [GeV^{2}]";
  
  const Int_t    xbin     = 21;
  const Double_t xbins[2][xbin+1] = { // [mm,ee]
    { 0.0,  1.0,                                                           2.0,                                          3.0,  4.3,  5.0,  6.0,  7.0,  (PDGmass::jpsi -0.25-0.0*0.15)*(PDGmass::jpsi -0.25-0.0*0.15), (PDGmass::jpsi +0.10+0.0*0.05)*(PDGmass::jpsi +0.10+0.0*0.05),
      11.0, (PDGmass::psi2s-0.15-0.0*0.10)*(PDGmass::psi2s-0.15-0.0*0.10), (PDGmass::psi2s+0.10)*(PDGmass::psi2s+0.10), 15.0, 16.0, 17.0, 18.0, 19.0,  20.0,                                                          21.5,
      23.0, 25.0,
    },
    { 0.0,  1.0,                                                           2.0,                                          3.0,  4.3,  5.0,  6.0, (PDGmass::jpsi -0.25-1.0*0.15)*(PDGmass::jpsi -0.25-1.0*0.15), 8.10, (PDGmass::jpsi +0.10+1.0*0.05)*(PDGmass::jpsi +0.10+1.0*0.05),
      11.2, (PDGmass::psi2s-0.15-1.0*0.10)*(PDGmass::psi2s-0.15-1.0*0.10), (PDGmass::psi2s+0.10)*(PDGmass::psi2s+0.10), 15.0, 16.0, 17.0, 18.0, 19.0,                                                          20.0, 21.5,
      23.0, 25.0,
    },
  };

  /*
  const Int_t    xbin_afb    =    9;
  const Double_t xbins_afb[2][xbin_afb+1] = { // [mm,ee]
    {0.0, 2.0, 4.3,
     (PDGmass::jpsi -0.25-0.0*0.15)*(PDGmass::jpsi -0.25-0.0*0.15), // ~ 8.10
     (PDGmass::jpsi +0.10+0.0*0.05)*(PDGmass::jpsi +0.10+0.0*0.05), // ~10.22
     (PDGmass::psi2s-0.15-0.0*0.10)*(PDGmass::psi2s-0.15-0.0*0.10), // ~12.50
     (PDGmass::psi2s+0.10                )*(PDGmass::psi2s+0.10                ), // ~14.33
     16.0, 19.0, 25.0},
    {0.0, 2.0, 4.3,
     (PDGmass::jpsi -0.25-1.0*0.15)*(PDGmass::jpsi -0.25-1.0*0.15), // ~ 7.27
     (PDGmass::jpsi +0.10+1.0*0.05)*(PDGmass::jpsi +0.10+1.0*0.05), // ~10.54
     (PDGmass::psi2s-0.15-1.0*0.10)*(PDGmass::psi2s-0.15-1.0*0.10), // ~11.81
     (PDGmass::psi2s+0.10                )*(PDGmass::psi2s+0.10                ), // ~14.33
     16.0, 19.0, 25.0},
  };
  const Int_t xbins_convert[2][xbin_afb+1] = {
    {1,3,5,9,10,12,13,15,18,22}, // mm
    {1,3,5,8,10,12,13,15,18,22}, // ee
  };
  */
  ///* default setting
  const Int_t    xbin_afb    =    6;
  const Double_t xbins_afb[2][xbin_afb+1] = { // [mm,ee]
    {0.0, 4.3,
     (PDGmass::jpsi -0.25-0.0*0.15)*(PDGmass::jpsi -0.25-0.0*0.15), // ~ 8.10
     (PDGmass::jpsi +0.10+0.0*0.05)*(PDGmass::jpsi +0.10+0.0*0.05), // ~10.22
     (PDGmass::psi2s-0.15-0.0*0.10)*(PDGmass::psi2s-0.15-0.0*0.10), // ~12.50
     (PDGmass::psi2s+0.10                )*(PDGmass::psi2s+0.10                ), // ~14.33
     25.0},
    {0.0, 4.3,
     (PDGmass::jpsi -0.25-1.0*0.15)*(PDGmass::jpsi -0.25-1.0*0.15), // ~ 7.27
     (PDGmass::jpsi +0.10+1.0*0.05)*(PDGmass::jpsi +0.10+1.0*0.05), // ~10.54
     (PDGmass::psi2s-0.15-1.0*0.10)*(PDGmass::psi2s-0.15-1.0*0.10), // ~11.81
     (PDGmass::psi2s+0.10                )*(PDGmass::psi2s+0.10                ), // ~14.33
     25.0},
  };
  const Int_t xbins_convert[2][xbin_afb+1] = {
    {1,5,9,10,12,13,22}, // mm
    {1,5,8,10,12,13,22}, // ee
  };
  //*/
  /* temporal setting for q2(1-6 GeV^2)
  const Int_t    xbin_afb    =    6;
  const Double_t xbins_afb[2][xbin_afb+1] = { // [mm,ee]
    {1.0, 4.3,
     6.0,
     (PDGmass::jpsi +0.10+0.0*0.05)*(PDGmass::jpsi +0.10+0.0*0.05), // ~10.22
     (PDGmass::psi2s-0.15-0.0*0.10)*(PDGmass::psi2s-0.15-0.0*0.10), // ~12.50
     (PDGmass::psi2s+0.10                )*(PDGmass::psi2s+0.10                ), // ~14.33
     25.0},
    {1.0, 4.3,
     6.0,
     (PDGmass::jpsi +0.10+1.0*0.05)*(PDGmass::jpsi +0.10+1.0*0.05), // ~10.54
     (PDGmass::psi2s-0.15-1.0*0.10)*(PDGmass::psi2s-0.15-1.0*0.10), // ~11.81
     (PDGmass::psi2s+0.10                )*(PDGmass::psi2s+0.10                ), // ~14.33
     25.0},
  };
  const Int_t xbins_convert[2][xbin_afb+1] = {
    {1,2,7,10,12,13,22}, // mm
    {1,2,7,10,12,13,22}, // ee
  };
  */
}

namespace nb_lep{
  const Char_t*  tname  = "h511";
  const Char_t*  axis   = "nb_lep%d_orgksfw_vtxcl_fmiss1_qq:nb_lep%d_orgksfw_vtxcl_fmiss1_bb";
  const Int_t    ybin   =    50;
  const Double_t ymin   = -1.00;
  const Double_t ymax   =  1.00;
  const Char_t*  ylabel = "NB_{qq}";
  const Int_t    xbin   =    30;
  const Double_t xmin   = -1.00;
  const Double_t xmax   =  1.00;
  const Char_t*  xlabel = "NB_{bb}";
}

namespace bgsup_var{
  const Int_t    n_axis     = 7;
  const Char_t*  tname      = "h511";
  const Char_t*  axis[n_axis]     = {         "evis1",   "abs(mmiss1)", "kfbchi/kfbdgf",      "kfbcl",                  "10000*dzll3d",          "bccm",                "de" };
  const Char_t*  fname[n_axis]    = {         "evis1",     "absmmiss1",       "kfbchi2",      "kfbcl",                        "dzll3d",          "bccm",                "de" };
  const Double_t offset[n_axis]   = {            0.0,              0.0,             0.0,          0.0,                             0.0,             0.0,                 0.0 };
  //const Int_t    xbin[n_axis]     = {             40,               40,              50,           25,                              38,              50,                  50 };
  //const Double_t xmin[n_axis]     = {            5.0,              0.0,             0.0,          0.0,                             0.0,            -1.0,               -0.05 };
  //const Double_t xmax[n_axis]     = {           13.0,              6.0,            10.0,          1.0,                             190,             1.0,                0.05 };
  const Int_t    xbin[n_axis]     = {             40,               40,              50,           25,                              38,              50,                  75 };
  const Double_t xmin[n_axis]     = {            5.0,              0.0,             0.0,          0.0,                             0.0,            -1.0,               -0.10 };
  const Double_t xmax[n_axis]     = {           13.0,              6.0,            10.0,          1.0,                             190,             1.0,                0.05 };
  const Double_t xmin_fit[n_axis] = {        xmin[0],          xmin[1],         xmin[2],      xmin[3],                         xmin[4],         xmin[5],             xmin[6] };
  const Double_t xmax_fit[n_axis] = {        xmax[0],          xmax[1],         xmax[2],      xmax[3],                         xmax[4],         xmax[5],             xmax[6] };
  const Char_t*  xlabel[n_axis]   = {"E_{vis} [GeV]", "M_{miss} [GeV]",  "#chi^{2}/NDF",  "CL(B-vtx)",  "|#Deltaz(l^{+}l^{-})| [#mum]", "cos#theta_{B}",    "#Delta E [GeV]" };
}

namespace bgsup_var_org{
  const Int_t    n_axis     = 7;
  const Char_t*  tname      = "h511";
  const Char_t*  axis[n_axis]     = {         "evis1",   "abs(mmiss1)", "kfbchiorg/kfbdgf",      "kfbclorg",               "10000*dzll3dorg",          "bccm",          "deorg" };
  const Char_t*  axisrd[n_axis]   = {         "evis1",   "abs(mmiss1)",    "kfbchi/kfbdgf",         "kfbcl",                  "10000*dzll3d",          "bccm",             "de" };
  const Char_t*  fname[n_axis]    = {         "evis1",     "absmmiss1",          "kfbchi2",         "kfbcl",                        "dzll3d",          "bccm",             "de" };
  const Double_t offset[n_axis]   = {            0.0,              0.0,                0.0,             0.0,                             0.0,             0.0,              0.0 };
  const Int_t    xbin[n_axis]     = {             40,               40,                 50,              25,                              38,              50,               75 };
  const Double_t xmin[n_axis]     = {            5.0,              0.0,                0.0,             0.0,                             0.0,            -1.0,            -0.10 };
  const Double_t xmax[n_axis]     = {           13.0,              6.0,               10.0,             1.0,                             190,             1.0,             0.05 };
  const Double_t xmin_fit[n_axis] = {        xmin[0],          xmin[1],            xmin[2],         xmin[3],                         xmin[4],         xmin[5],          xmin[6] };
  const Double_t xmax_fit[n_axis] = {        xmax[0],          xmax[1],            xmax[2],         xmax[3],                         xmax[4],         xmax[5],          xmax[6] };
  const Char_t*  xlabel[n_axis]   = {"E_{vis} [GeV]", "M_{miss} [GeV]",     "#chi^{2}/NDF",     "CL(B-vtx)",  "|#Deltaz(l^{+}l^{-})| [#mum]", "cos#theta_{B}", "#Delta E [GeV]" };
}

namespace bgsup_var_calib{
  const Int_t    n_axis     = 7;
  const Char_t*  tname      = "h511";
  const Char_t*  axis[n_axis]     = {         "evis1",   "abs(mmiss1)", "kfbchicalib/kfbdgf", "kfbclcalib",              "10000*dzll3dcalib",          "bccm",        "decalib" };
  const Char_t*  axisrd[n_axis]   = {         "evis1",   "abs(mmiss1)",      "kfbchi/kfbdgf",      "kfbcl",                   "10000*dzll3d",          "bccm",             "de" };
  const Char_t*  fname[n_axis]    = {         "evis1",     "absmmiss1",            "kfbchi2",      "kfbcl",                         "dzll3d",          "bccm",             "de" };
  const Double_t offset[n_axis]   = {            0.0,              0.0,                  0.0,          0.0,                              0.0,             0.0,              0.0 };
  const Int_t    xbin[n_axis]     = {             40,               40,                   50,           25,                               38,              50,               75 };
  const Double_t xmin[n_axis]     = {            5.0,              0.0,                  0.0,          0.0,                              0.0,            -1.0,            -0.10 };
  const Double_t xmax[n_axis]     = {           13.0,              6.0,                 10.0,          1.0,                              190,             1.0,             0.05 };
  const Double_t xmin_fit[n_axis] = {        xmin[0],          xmin[1],              xmin[2],      xmin[3],                          xmin[4],         xmin[5],          xmin[6] };
  const Double_t xmax_fit[n_axis] = {        xmax[0],          xmax[1],              xmax[2],      xmax[3],                          xmax[4],         xmax[5],          xmax[6] };
  const Char_t*  xlabel[n_axis]   = {"E_{vis} [GeV]", "M_{miss} [GeV]",       "#chi^{2}/NDF",  "CL(B-vtx)",  "|#Deltaz(l^{+}l^{-})| [#mum]", "cos#theta_{B}", "#Delta E [GeV]" };
}

#endif
