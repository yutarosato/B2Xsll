#ifndef DRAWS_H
#define DRAWS_H

#include <iostream>
#include <sstream>
#include <TROOT.h>

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// if this section is changed, "Set/makeCut.cpp" must be changed.
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
const Int_t nmode       = 18;
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


const Int_t   nbkgtype          = 4;
const Char_t* bkgtype[nbkgtype] = {"uds","charm","mixed","charged"};

const Int_t   nexpno = 30;
const Char_t* expno[nexpno]  = {"07","09","11","13","15","17","19a","19b","21","23","25","27",
				"31","33","35","37a","37b","39","41a","41b","43","45a","45b","47","49","51","55","61","63","65"};

const Int_t   nexpno_gmc = 26;
const Char_t* expno_gmc[nexpno_gmc]  = {"07","09","11","13","15","17","19","21","23","25","27",
					"31","33","35","37","39","41","43","45","47","49","51","55","61","63","65"};
const Double_t nbb_gmc[nexpno_gmc] = {6.4588, 4.7597, 8.8509, 11.6998, 13.5679, 12.4588, 16.3023+10.8682, 4.3371, 6.4755, 28.0008,
				      28.1814, 19.6601, 18.3460, 17.3783, 30.5480+36.6339, 47.0818, 34.0657+29.9477, 61.5614, 6.64662+7.70718,
				      41.2186, 29.7205, 41.8919, 80.2472, 37.4460, 35.6231, 41.7867};
const Char_t* lrnb[2] = {"lr", "nb"};

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

const Int_t nset = 21;
const Char_t* set_name[nset] = {"A","B","C","D","E","F","G",
				"H","I","J","K","L","M","N",
				"O","P","Q","R","S","T","U"
};
const Double_t sigmc_amount = 540.8;

namespace gmc_fl{
  const Char_t* indir = "hbk_double_nb/hbk_orgksfw_vtxcl_fmiss1_weight/";
  //const Char_t* indir = "hbk_swap_nb/hbk_orgksfw_vtxcl_fmiss1_weight/";
  const Char_t* tail  = "*.root";
}

namespace gmc_fl_double{
  const Char_t* indir = "hbk_double_orgksfw_vtxcl_fmiss1_lrnb/";
  const Char_t* tail  = "*.root";
}

namespace gmc_fl_single{
  const Char_t* indir = "hbk_single_orgksfw_vtxcl_fmiss1_lrnb/";
  const Char_t* tail  = "*.root";
}

const Int_t col_fil[30] = {4,9,3,2,15,12,1,16,5,13,6,11,7,8,10,14, 17,18,19,20,21,22,23,24,25,26,27,28,29,30};

namespace Mbc_comb{
  const Char_t*  tname    = "h511";
  const Char_t*  axis     = "Mbc";
  const Double_t offset   =  0.0;
  const Int_t    xbin     =   40;
  const Double_t xmin     = 5.22;
  const Double_t xmax     = 5.30;
  const Double_t xmin_fit = xmin;
  const Double_t xmax_fit = xmax;
  const Char_t*  xlabel   = "M_{bc} [GeV]";
}

namespace q2_theta_uniform{
  const Char_t*  tname  = "h511";
  const Char_t*  axis   = "-coslp:cc_m*cc_m";
  const Int_t    ybin   =    10;
  const Double_t ymin   = -1.00;
  const Double_t ymax   =  1.00;
  const Char_t*  ylabel = "cos#theta";
  const Int_t    xbin   =    25;
  const Double_t xmin   =   0.0;
  const Double_t xmax   =  25.0;
  const Char_t*  xlabel = "q^{2} [GeV^{2}]";
}

namespace q2_theta_nonuniform{
  
  const Char_t*  tname   = "h511";
  const Char_t*  axis    = "-coslp:cc_m*cc_m";
  const Int_t    ybin    =    10;
  const Double_t ymin    = -1.00;
  const Double_t ymax    =  1.00;
  const Char_t*  ylabel  = "cos#theta";
  const Double_t xmax    =  25.0;
  const Char_t*  xlabel  = "q^{2} [GeV^{2}]";
  
  const Int_t    xbin    =    9;
  
  const Double_t xbins[2][xbin+1] = { // [mm,ee]
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
}  

#endif
