#ifndef DRAWS_H
#define DRAWS_H

#include <iostream>
#include <sstream>
#include <TROOT.h>


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


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

namespace gmc{
  const Char_t* indir = "~/ewp/ana/data/dstr_fl/hbk/";
  const Char_t* tail  = "*.root";
}

namespace rd{
  const Char_t* indir = "~/ewp/ana/data/dstr_fl/hbk_data/";
  const Char_t* tail  = "*.root";
}

const Int_t col_fil[30] = {4,9,3,2,15,12,1,16,5,13,6,11,7,8,10,14, 17,18,19,20,21,22,23,24,25,26,27,28,29,30};



namespace var{
  const Int_t    n_axis     = 11;
  const Char_t*  tname      = "h1";
  const Char_t*  axis[n_axis]     = {               "m_d",      "m_d_vf",              "m_dst",      "m_dst_vf",                   "dm",          "dm_vf",          "xd", "xd_vf",         "xdst", "xdst_vf",               "m_drev" };
  const Char_t*  fname[n_axis]    = {             "dmorg",          "dm",            "dstmorg",          "dstm",            "deltamorg",         "deltam",       "xdorg",    "xd",      "xdstorg",    "xdst",               "dm_rev" };
  const Double_t offset[n_axis]   = {         PDGmass::d0,   PDGmass::d0,       PDGmass::dstrp,  PDGmass::dstrp,                    0.0,              0.0,           0.0,     0.0,            0.0,       0.0,            PDGmass::d0 };
  const Int_t    xbin[n_axis]     = {                 120,           120,                  120,             120,                    100,              100,           100,     100,            100,       100,                    120 };
  const Double_t xmin[n_axis]     = {              -0.065,        -0.065,               -0.060,          -0.060,                  0.135,            0.135,           0.0,     0.0,            0.0,       0.0,                 -0.060 };
  const Double_t xmax[n_axis]     = {               0.055,         0.055,                0.060,           0.060,                  0.155,            0.155,           1.0,     1.0,            1.0,       1.0,                  0.060 };
  const Char_t*  xlabel[n_axis]   = { "M_{D}^{org} [Gev}", "M_{D} [GeV]", "M_{D*}^{org} [GeV]",  "M_{D*} [GeV]", "#Delta M^{org} [GeV]", "#Delta M [GeV]", "X_{D}^{org}", "X_{D}", "X_{D*}^{org}",  "X_{D*}",    "M_{D}^{rev} [GeV]" };
}

namespace pid{
  const Int_t    n_axis     = 3;
  const Char_t*  tname      = "h1";
  const Char_t*  axis[n_axis]     = {  "mulmu_%d", "eprob_%d", "npatc1_%d" };
  const Char_t*  fname[n_axis]    = {     "mu-id",     "e-id",       "kid" };
  const Double_t offset[n_axis]   = {         0.0,        0.0,         0.0 };
  const Int_t    xbin[n_axis]     = {         100,        100,         100 };
  const Double_t xmin[n_axis]     = {         0.0,        0.0,         0.0 };
  const Double_t xmax[n_axis]     = {         1.0,        1.0,         1.0 };
  const Char_t*  xlabel[n_axis]   = {    "#mu-ID",     "e-ID",      "K-ID" };
}

namespace d0mass{
  const Char_t*  tname    = "h1";
  const Char_t*  axis     = "m_d_vf";
  const Char_t*  fname    = "dm";
  const Double_t offset   = PDGmass::d0;
  const Int_t    xbin     = 120;
  const Double_t xmin     = -0.065;
  const Double_t xmin_fit = -0.065;
  const Double_t xmax     =  0.055;
  const Double_t xmax_fit =  0.055;
  const Char_t*  xlabel   = "M_{D} [GeV]";

  const Int_t    n_axis          =    2; // 0.4 GeV width, 0.6 GeV width
  const Int_t    n_mom[n_axis]   =  {  15,   10};
  const Double_t mom_min[n_axis] =  {0.40, 0.40};
  const Double_t mom_max[n_axis] =  {6.40, 6.40};
  const Double_t del_mom[n_axis] =  { (mom_max[0]-mom_min[0])/n_mom[0],
				      (mom_max[1]-mom_min[1])/n_mom[1] };
  const Int_t    n_cos           =   10;
  const Double_t cos_min         = -1.0;
  const Double_t cos_max         =  1.0;
  const Double_t del_cos         = (cos_max-cos_min)/n_cos;
}

namespace p_cos{
  const Int_t    n_axis       = 2; // 0.4 GeV width, 0.6 GeV width
  const Char_t*  tname        = "h1";
  const Char_t*  axis         = "pmag_%d:cos_%d";
  const Int_t    ybin[n_axis] =  {  15,   10};
  const Double_t ymin[n_axis] =  {0.40, 0.40};
  const Double_t ymax[n_axis] =  {6.40, 6.40};
  const Char_t*  ylabel       = "P [GeV]";
  const Double_t offset       =   0.0;
  const Int_t    xbin         =    10;
  const Double_t xmin         =  -1.0;
  const Double_t xmax         =   1.0;
  const Char_t*  xlabel       = "cos#theta";
}

#endif
