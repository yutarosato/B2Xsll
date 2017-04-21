#include <iostream>
#include <math.h>

#include "../Util/Style.h"
#include "../Util/Canvas.h"
#include "../Util/const.h"
#include "../Util/MChain.h"
#include "../Util/MCut.h"
#include "../Util/MCut_array.h"
#include "../Util/FitFunc.h"
#include "../Util/Manip.h"
#include "../Set/Nominal_cut_selection.h"
#include "../Set/Branch.h"
#include "../Set/makeCut.h"

#include "draws_.h"

#include <vector>
#include <stdlib.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TString.h>
#include <TArrow.h>

using namespace q2_theta_nonuniform_eff;
const Bool_t flag_save      = true;
const Bool_t fl_corr        = true; // correct PID systematic

Int_t main( Int_t argc, Char_t** argv ){
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==5 || argc==6) ) std::cerr << "wrong input" << std::endl
					<< " Usage : ./draws_2d (int)fl_mode_ll (char*)fl_A7 (char*)fl_A9 (char*)fl_A10 [(int)fl_appRun]" << std::endl
					<< "[fl_mode_ll] : 1(e), 0(mu)"
					<< std::endl, abort();
  
  Int_t    fl_mode_ll = atoi( argv[1] ); // 1(e), 0(mu)
  Char_t*  fl_A7      = argv[2];
  Char_t*  fl_A9      = argv[3];
  Char_t*  fl_A10     = argv[4];
  Char_t*  setname    = "A";
  Double_t used_nset  = 1;
  Int_t   fl_appRun  = 1;
  if( argc==6 ) fl_appRun = atoi( argv[5] );
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t  Nchain  = 1; // fixed
  const Int_t  Ntmp    = 19; // K, K*, Xs(16mode), missing
  const Int_t  nfile[Nchain] = {0};
  const Double_t scale_event_sig = sigmc_amount*(used_nset/(Double_t)nset); // sigmc : N -> N/alpha
  const Double_t add[Nchain][Ntmp] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
  
  Char_t** infile = new Char_t*[Nchain];
  for( Int_t i=0; i<Nchain; i++ ) infile[i] = new Char_t[1024]; 
  
  std::stringstream sTmp;
  for( Int_t i=0; i<Nchain; i++ ){
    sTmp << "~/Store/ewp/ana/data/gen_mb/hbk1/hbk/hbk_" << fl_A7 << "_" << fl_A9 << "_" << fl_A10 << "_set" << setname << "_mb495/Gen_";
    strcpy( infile[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Char_t** add_cut = new Char_t*[Ntmp];
  for( Int_t i=0; i<Ntmp; i++ ){
    add_cut[i] = new Char_t[4096];
    if     ( i==0 ) sTmp << "gm_fl_xs<0 && Xs_m<0.50"; // K ll
    else if( i==1 ) sTmp << "gm_fl_xs<0 && Xs_m>0.50"; // K*ll
    else if( i==18 ){
      sTmp << "( gm_bg!=0 || ( gm_bg==0 && "; // Xsll [missing modes]
      for( Int_t k=0; k<nmode; k++ ){
	sTmp << " gm_xs!=" << mode[k];
	if( k!= nmode-1 ) sTmp << " && ";
      }
      sTmp << ") ) && gm_fl_xs>0";
    }else sTmp << "gm_fl_xs>0 && gm_bg==0 && gm_xs==" << mode[i%nmode]; // Xsll [each mode]    
    strcpy( add_cut[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Double_t eff_corr_lep_9qbin [3][15] = {0}; // [mm,ee,ee+mm]
  Double_t eff_corr_lepE_9qbin[3][15] = {0}; // [mm,ee,ee+mm]
  eff_corr_lep_9qbin[0][ 0] =     0.938681; eff_corr_lepE_9qbin[0][ 0] =    0.0444138;
  eff_corr_lep_9qbin[0][ 1] =     0.938007; eff_corr_lepE_9qbin[0][ 1] =    0.0442265;
  eff_corr_lep_9qbin[0][ 2] =     0.919607; eff_corr_lepE_9qbin[0][ 2] =    0.0499837;
  eff_corr_lep_9qbin[0][ 3] =     0.918577; eff_corr_lepE_9qbin[0][ 3] =    0.0505896;
  eff_corr_lep_9qbin[0][ 4] =     0.913339; eff_corr_lepE_9qbin[0][ 4] =     0.049207;
  eff_corr_lep_9qbin[0][ 5] =     0.912478; eff_corr_lepE_9qbin[0][ 5] =    0.0494998;
  eff_corr_lep_9qbin[0][ 6] =     0.912588; eff_corr_lepE_9qbin[0][ 6] =    0.0453009;
  eff_corr_lep_9qbin[0][ 7] =      0.91212; eff_corr_lepE_9qbin[0][ 7] =     0.045706;
  eff_corr_lep_9qbin[0][ 8] =     0.901967; eff_corr_lepE_9qbin[0][ 8] =    0.0454063;
  eff_corr_lep_9qbin[0][ 9] =     0.901785; eff_corr_lepE_9qbin[0][ 9] =    0.0456503;
  eff_corr_lep_9qbin[0][10] =     0.895083; eff_corr_lepE_9qbin[0][10] =    0.0454866;
  eff_corr_lep_9qbin[0][11] =     0.895195; eff_corr_lepE_9qbin[0][11] =    0.0459444;
  eff_corr_lep_9qbin[0][12] =     0.885115; eff_corr_lepE_9qbin[0][12] =    0.0467886;
  eff_corr_lep_9qbin[0][13] =      0.88744; eff_corr_lepE_9qbin[0][13] =    0.0458132;
  eff_corr_lep_9qbin[0][14] =      0.91277; eff_corr_lepE_9qbin[0][14] =    0.0469858;
  eff_corr_lep_9qbin[1][ 0] =     0.952387; eff_corr_lepE_9qbin[1][ 0] =    0.0363613;
  eff_corr_lep_9qbin[1][ 1] =      0.95209; eff_corr_lepE_9qbin[1][ 1] =    0.0364467;
  eff_corr_lep_9qbin[1][ 2] =     0.958766; eff_corr_lepE_9qbin[1][ 2] =     0.035795;
  eff_corr_lep_9qbin[1][ 3] =     0.958868; eff_corr_lepE_9qbin[1][ 3] =    0.0361764;
  eff_corr_lep_9qbin[1][ 4] =     0.962943; eff_corr_lepE_9qbin[1][ 4] =    0.0343084;
  eff_corr_lep_9qbin[1][ 5] =     0.962774; eff_corr_lepE_9qbin[1][ 5] =    0.0343057;
  eff_corr_lep_9qbin[1][ 6] =     0.965919; eff_corr_lepE_9qbin[1][ 6] =    0.0298484;
  eff_corr_lep_9qbin[1][ 7] =      0.96589; eff_corr_lepE_9qbin[1][ 7] =    0.0299268;
  eff_corr_lep_9qbin[1][ 8] =     0.965587; eff_corr_lepE_9qbin[1][ 8] =    0.0294586;
  eff_corr_lep_9qbin[1][ 9] =     0.965681; eff_corr_lepE_9qbin[1][ 9] =    0.0295039;
  eff_corr_lep_9qbin[1][10] =     0.964128; eff_corr_lepE_9qbin[1][10] =    0.0296922;
  eff_corr_lep_9qbin[1][11] =     0.964012; eff_corr_lepE_9qbin[1][11] =    0.0298549;
  eff_corr_lep_9qbin[1][12] =     0.957566; eff_corr_lepE_9qbin[1][12] =     0.029983;
  eff_corr_lep_9qbin[1][13] =     0.958045; eff_corr_lepE_9qbin[1][13] =     0.030115;
  eff_corr_lep_9qbin[1][14] =     0.960233; eff_corr_lepE_9qbin[1][14] =     0.033481;
  eff_corr_lep_9qbin[2][ 0] =     0.946697; eff_corr_lepE_9qbin[2][ 0] =    0.0397114;
  eff_corr_lep_9qbin[2][ 1] =     0.946211; eff_corr_lepE_9qbin[2][ 1] =    0.0397021;
  eff_corr_lep_9qbin[2][ 2] =      0.94094; eff_corr_lepE_9qbin[2][ 2] =    0.0422952;
  eff_corr_lep_9qbin[2][ 3] =     0.940714; eff_corr_lepE_9qbin[2][ 3] =    0.0427139;
  eff_corr_lep_9qbin[2][ 4] =     0.935077; eff_corr_lepE_9qbin[2][ 4] =    0.0427339;
  eff_corr_lep_9qbin[2][ 5] =      0.93462; eff_corr_lepE_9qbin[2][ 5] =    0.0428686;
  eff_corr_lep_9qbin[2][ 6] =     0.928253; eff_corr_lepE_9qbin[2][ 6] =    0.0408143;
  eff_corr_lep_9qbin[2][ 7] =     0.927372; eff_corr_lepE_9qbin[2][ 7] =    0.0412825;
  eff_corr_lep_9qbin[2][ 8] =     0.928985; eff_corr_lepE_9qbin[2][ 8] =     0.038711;
  eff_corr_lep_9qbin[2][ 9] =     0.928843; eff_corr_lepE_9qbin[2][ 9] =    0.0388913;
  eff_corr_lep_9qbin[2][10] =     0.925384; eff_corr_lepE_9qbin[2][10] =    0.0386404;
  eff_corr_lep_9qbin[2][11] =     0.925187; eff_corr_lepE_9qbin[2][11] =    0.0390185;
  eff_corr_lep_9qbin[2][12] =     0.916564; eff_corr_lepE_9qbin[2][12] =    0.0395898;
  eff_corr_lep_9qbin[2][13] =     0.918008; eff_corr_lepE_9qbin[2][13] =    0.0391044;
  eff_corr_lep_9qbin[2][14] =     0.934461; eff_corr_lepE_9qbin[2][14] =    0.0408627;

  Double_t eff_corr_lep_6qbin [3][9] = {0}; // [mm,ee,ee+mm]
  Double_t eff_corr_lepE_6qbin[3][9] = {0}; // [mm,ee,ee+mm]
  eff_corr_lep_6qbin[0][ 0] =     0.929641; eff_corr_lepE_6qbin[0][ 0] =    0.0470587;
  eff_corr_lep_6qbin[0][ 1] =     0.928832; eff_corr_lepE_6qbin[0][ 1] =    0.0472332;
  eff_corr_lep_6qbin[0][ 2] =     0.913339; eff_corr_lepE_6qbin[0][ 2] =     0.049207;
  eff_corr_lep_6qbin[0][ 3] =     0.912478; eff_corr_lepE_6qbin[0][ 3] =    0.0494998;
  eff_corr_lep_6qbin[0][ 4] =     0.912588; eff_corr_lepE_6qbin[0][ 4] =    0.0453009;
  eff_corr_lep_6qbin[0][ 5] =      0.91212; eff_corr_lepE_6qbin[0][ 5] =     0.045706;
  eff_corr_lep_6qbin[0][ 6] =     0.896945; eff_corr_lepE_6qbin[0][ 6] =    0.0455876;
  eff_corr_lep_6qbin[0][ 7] =     0.897226; eff_corr_lepE_6qbin[0][ 7] =    0.0458032;
  eff_corr_lep_6qbin[0][ 8] =      0.91277; eff_corr_lepE_6qbin[0][ 8] =    0.0469858;
  eff_corr_lep_6qbin[1][ 0] =     0.955137; eff_corr_lepE_6qbin[1][ 0] =    0.0361064;
  eff_corr_lep_6qbin[1][ 1] =     0.955041; eff_corr_lepE_6qbin[1][ 1] =    0.0363215;
  eff_corr_lep_6qbin[1][ 2] =     0.962943; eff_corr_lepE_6qbin[1][ 2] =    0.0343084;
  eff_corr_lep_6qbin[1][ 3] =     0.962774; eff_corr_lepE_6qbin[1][ 3] =    0.0343057;
  eff_corr_lep_6qbin[1][ 4] =     0.965919; eff_corr_lepE_6qbin[1][ 4] =    0.0298484;
  eff_corr_lep_6qbin[1][ 5] =      0.96589; eff_corr_lepE_6qbin[1][ 5] =    0.0299268;
  eff_corr_lep_6qbin[1][ 6] =     0.964034; eff_corr_lepE_6qbin[1][ 6] =     0.029623;
  eff_corr_lep_6qbin[1][ 7] =     0.964094; eff_corr_lepE_6qbin[1][ 7] =    0.0297312;
  eff_corr_lep_6qbin[1][ 8] =     0.960233; eff_corr_lepE_6qbin[1][ 8] =     0.033481;
  eff_corr_lep_6qbin[2][ 0] =     0.944094; eff_corr_lepE_6qbin[2][ 0] =      0.04087;
  eff_corr_lep_6qbin[2][ 1] =     0.943708; eff_corr_lepE_6qbin[2][ 1] =    0.0410601;
  eff_corr_lep_6qbin[2][ 2] =     0.935077; eff_corr_lepE_6qbin[2][ 2] =    0.0427339;
  eff_corr_lep_6qbin[2][ 3] =      0.93462; eff_corr_lepE_6qbin[2][ 3] =    0.0428686;
  eff_corr_lep_6qbin[2][ 4] =     0.928253; eff_corr_lepE_6qbin[2][ 4] =    0.0408143;
  eff_corr_lep_6qbin[2][ 5] =     0.927372; eff_corr_lepE_6qbin[2][ 5] =    0.0412825;
  eff_corr_lep_6qbin[2][ 6] =     0.925959; eff_corr_lepE_6qbin[2][ 6] =    0.0387662;
  eff_corr_lep_6qbin[2][ 7] =     0.926001; eff_corr_lepE_6qbin[2][ 7] =    0.0389701;
  eff_corr_lep_6qbin[2][ 8] =     0.934461; eff_corr_lepE_6qbin[2][ 8] =    0.0408627;
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Double_t eff_corr_kpi_9qbin [3][15] = {0}; // [mm,ee,ee+mm]
  Double_t eff_corr_kpiE_9qbin[3][15] = {0}; // [mm,ee,ee+mm]
  eff_corr_kpi_9qbin[0][ 0] =     0.966477; eff_corr_kpiE_9qbin[0][ 0] =    0.0143941;
  eff_corr_kpi_9qbin[0][ 1] =     0.966155; eff_corr_kpiE_9qbin[0][ 1] =    0.0144256;
  eff_corr_kpi_9qbin[0][ 2] =     0.973781; eff_corr_kpiE_9qbin[0][ 2] =    0.0140941;
  eff_corr_kpi_9qbin[0][ 3] =     0.973953; eff_corr_kpiE_9qbin[0][ 3] =    0.0140505;
  eff_corr_kpi_9qbin[0][ 4] =      0.98174; eff_corr_kpiE_9qbin[0][ 4] =    0.0143237;
  eff_corr_kpi_9qbin[0][ 5] =     0.981359; eff_corr_kpiE_9qbin[0][ 5] =    0.0143095;
  eff_corr_kpi_9qbin[0][ 6] =     0.993737; eff_corr_kpiE_9qbin[0][ 6] =    0.0149821;
  eff_corr_kpi_9qbin[0][ 7] =      0.99401; eff_corr_kpiE_9qbin[0][ 7] =    0.0151417;
  eff_corr_kpi_9qbin[0][ 8] =      1.00549; eff_corr_kpiE_9qbin[0][ 8] =    0.0166249;
  eff_corr_kpi_9qbin[0][ 9] =      1.00533; eff_corr_kpiE_9qbin[0][ 9] =    0.0162339;
  eff_corr_kpi_9qbin[0][10] =      1.00951; eff_corr_kpiE_9qbin[0][10] =     0.016527;
  eff_corr_kpi_9qbin[0][11] =       1.0094; eff_corr_kpiE_9qbin[0][11] =    0.0162872;
  eff_corr_kpi_9qbin[0][12] =      1.01311; eff_corr_kpiE_9qbin[0][12] =    0.0123365;
  eff_corr_kpi_9qbin[0][13] =      1.01414; eff_corr_kpiE_9qbin[0][13] =    0.0119787;
  eff_corr_kpi_9qbin[0][14] =     0.986483; eff_corr_kpiE_9qbin[0][14] =    0.0142676;
  eff_corr_kpi_9qbin[1][ 0] =     0.965816; eff_corr_kpiE_9qbin[1][ 0] =    0.0142642;
  eff_corr_kpi_9qbin[1][ 1] =      0.96641; eff_corr_kpiE_9qbin[1][ 1] =    0.0143336;
  eff_corr_kpi_9qbin[1][ 2] =     0.974646; eff_corr_kpiE_9qbin[1][ 2] =    0.0133292;
  eff_corr_kpi_9qbin[1][ 3] =     0.975211; eff_corr_kpiE_9qbin[1][ 3] =    0.0136397;
  eff_corr_kpi_9qbin[1][ 4] =     0.981267; eff_corr_kpiE_9qbin[1][ 4] =    0.0137068;
  eff_corr_kpi_9qbin[1][ 5] =     0.980834; eff_corr_kpiE_9qbin[1][ 5] =    0.0137672;
  eff_corr_kpi_9qbin[1][ 6] =     0.994031; eff_corr_kpiE_9qbin[1][ 6] =    0.0145763;
  eff_corr_kpi_9qbin[1][ 7] =     0.993064; eff_corr_kpiE_9qbin[1][ 7] =    0.0141824;
  eff_corr_kpi_9qbin[1][ 8] =      1.00537; eff_corr_kpiE_9qbin[1][ 8] =    0.0160447;
  eff_corr_kpi_9qbin[1][ 9] =      1.00588; eff_corr_kpiE_9qbin[1][ 9] =    0.0159039;
  eff_corr_kpi_9qbin[1][10] =      1.00977; eff_corr_kpiE_9qbin[1][10] =    0.0162345;
  eff_corr_kpi_9qbin[1][11] =       1.0102; eff_corr_kpiE_9qbin[1][11] =    0.0163009;
  eff_corr_kpi_9qbin[1][12] =      1.01245; eff_corr_kpiE_9qbin[1][12] =    0.0120263;
  eff_corr_kpi_9qbin[1][13] =      1.01326; eff_corr_kpiE_9qbin[1][13] =    0.0124965;
  eff_corr_kpi_9qbin[1][14] =     0.982732; eff_corr_kpiE_9qbin[1][14] =     0.013756;
  eff_corr_kpi_9qbin[2][ 0] =      0.96609; eff_corr_kpiE_9qbin[2][ 0] =     0.014318;
  eff_corr_kpi_9qbin[2][ 1] =     0.966304; eff_corr_kpiE_9qbin[2][ 1] =    0.0143719;
  eff_corr_kpi_9qbin[2][ 2] =     0.974254; eff_corr_kpiE_9qbin[2][ 2] =    0.0136754;
  eff_corr_kpi_9qbin[2][ 3] =     0.974647; eff_corr_kpiE_9qbin[2][ 3] =    0.0138237;
  eff_corr_kpi_9qbin[2][ 4] =     0.981531; eff_corr_kpiE_9qbin[2][ 4] =    0.0140514;
  eff_corr_kpi_9qbin[2][ 5] =     0.981126; eff_corr_kpiE_9qbin[2][ 5] =     0.014069;
  eff_corr_kpi_9qbin[2][ 6] =     0.993824; eff_corr_kpiE_9qbin[2][ 6] =    0.0148617;
  eff_corr_kpi_9qbin[2][ 7] =     0.993739; eff_corr_kpiE_9qbin[2][ 7] =    0.0148668;
  eff_corr_kpi_9qbin[2][ 8] =      1.00544; eff_corr_kpiE_9qbin[2][ 8] =     0.016376;
  eff_corr_kpi_9qbin[2][ 9] =      1.00556; eff_corr_kpiE_9qbin[2][ 9] =    0.0160928;
  eff_corr_kpi_9qbin[2][10] =      1.00963; eff_corr_kpiE_9qbin[2][10] =    0.0163973;
  eff_corr_kpi_9qbin[2][11] =      1.00975; eff_corr_kpiE_9qbin[2][11] =    0.0162932;
  eff_corr_kpi_9qbin[2][12] =      1.01282; eff_corr_kpiE_9qbin[2][12] =    0.0122003;
  eff_corr_kpi_9qbin[2][13] =      1.01375; eff_corr_kpiE_9qbin[2][13] =    0.0122053;
  eff_corr_kpi_9qbin[2][14] =     0.984757; eff_corr_kpiE_9qbin[2][14] =    0.0140322;

  Double_t eff_corr_kpi_6qbin [3][9] = {0}; // [mm,ee,ee+mm]
  Double_t eff_corr_kpiE_6qbin[3][9] = {0}; // [mm,ee,ee+mm]
  eff_corr_kpi_6qbin[0][ 0] =     0.969916; eff_corr_kpiE_6qbin[0][ 0] =    0.0141573;
  eff_corr_kpi_6qbin[0][ 1] =     0.969791; eff_corr_kpiE_6qbin[0][ 1] =    0.0141594;
  eff_corr_kpi_6qbin[0][ 2] =      0.98174; eff_corr_kpiE_6qbin[0][ 2] =    0.0143237;
  eff_corr_kpi_6qbin[0][ 3] =     0.981359; eff_corr_kpiE_6qbin[0][ 3] =    0.0143095;
  eff_corr_kpi_6qbin[0][ 4] =     0.993737; eff_corr_kpiE_6qbin[0][ 4] =    0.0149821;
  eff_corr_kpi_6qbin[0][ 5] =      0.99401; eff_corr_kpiE_6qbin[0][ 5] =    0.0151417;
  eff_corr_kpi_6qbin[0][ 6] =      1.00804; eff_corr_kpiE_6qbin[0][ 6] =    0.0157863;
  eff_corr_kpi_6qbin[0][ 7] =       1.0079; eff_corr_kpiE_6qbin[0][ 7] =    0.0155454;
  eff_corr_kpi_6qbin[0][ 8] =     0.986483; eff_corr_kpiE_6qbin[0][ 8] =    0.0142676;
  eff_corr_kpi_6qbin[1][ 0] =     0.969619; eff_corr_kpiE_6qbin[1][ 0] =    0.0137762;
  eff_corr_kpi_6qbin[1][ 1] =     0.970237; eff_corr_kpiE_6qbin[1][ 1] =    0.0139455;
  eff_corr_kpi_6qbin[1][ 2] =     0.981267; eff_corr_kpiE_6qbin[1][ 2] =    0.0137068;
  eff_corr_kpi_6qbin[1][ 3] =     0.980834; eff_corr_kpiE_6qbin[1][ 3] =    0.0137672;
  eff_corr_kpi_6qbin[1][ 4] =     0.994031; eff_corr_kpiE_6qbin[1][ 4] =    0.0145763;
  eff_corr_kpi_6qbin[1][ 5] =     0.993064; eff_corr_kpiE_6qbin[1][ 5] =    0.0141824;
  eff_corr_kpi_6qbin[1][ 6] =      1.00816; eff_corr_kpiE_6qbin[1][ 6] =    0.0154065;
  eff_corr_kpi_6qbin[1][ 7] =      1.00849; eff_corr_kpiE_6qbin[1][ 7] =    0.0153822;
  eff_corr_kpi_6qbin[1][ 8] =     0.982732; eff_corr_kpiE_6qbin[1][ 8] =     0.013756;
  eff_corr_kpi_6qbin[2][ 0] =     0.969747; eff_corr_kpiE_6qbin[2][ 0] =    0.0139406;
  eff_corr_kpi_6qbin[2][ 1] =     0.970045; eff_corr_kpiE_6qbin[2][ 1] =    0.0140376;
  eff_corr_kpi_6qbin[2][ 2] =     0.981531; eff_corr_kpiE_6qbin[2][ 2] =    0.0140514;
  eff_corr_kpi_6qbin[2][ 3] =     0.981126; eff_corr_kpiE_6qbin[2][ 3] =     0.014069;
  eff_corr_kpi_6qbin[2][ 4] =     0.993824; eff_corr_kpiE_6qbin[2][ 4] =    0.0148617;
  eff_corr_kpi_6qbin[2][ 5] =     0.993739; eff_corr_kpiE_6qbin[2][ 5] =    0.0148668;
  eff_corr_kpi_6qbin[2][ 6] =       1.0081; eff_corr_kpiE_6qbin[2][ 6] =    0.0156204;
  eff_corr_kpi_6qbin[2][ 7] =      1.00816; eff_corr_kpiE_6qbin[2][ 7] =    0.0154745;
  eff_corr_kpi_6qbin[2][ 8] =     0.984757; eff_corr_kpiE_6qbin[2][ 8] =    0.0140322;

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Double_t eff_corr_ks_9qbin [3][15] = {0}; // [mm,ee,ee+mm]
  Double_t eff_corr_ksE_9qbin[3][15] = {0}; // [mm,ee,ee+mm]
  eff_corr_ks_9qbin[0][ 0] =     0.997432; eff_corr_ksE_9qbin[0][ 0] =   0.00603736;
  eff_corr_ks_9qbin[0][ 1] =     0.997421; eff_corr_ksE_9qbin[0][ 1] =   0.00603514;
  eff_corr_ks_9qbin[0][ 2] =     0.997821; eff_corr_ksE_9qbin[0][ 2] =   0.00602076;
  eff_corr_ks_9qbin[0][ 3] =     0.997837; eff_corr_ksE_9qbin[0][ 3] =   0.00602042;
  eff_corr_ks_9qbin[0][ 4] =     0.997604; eff_corr_ksE_9qbin[0][ 4] =   0.00601739;
  eff_corr_ks_9qbin[0][ 5] =     0.997532; eff_corr_ksE_9qbin[0][ 5] =   0.00601935;
  eff_corr_ks_9qbin[0][ 6] =     0.997256; eff_corr_ksE_9qbin[0][ 6] =   0.00602106;
  eff_corr_ks_9qbin[0][ 7] =     0.997295; eff_corr_ksE_9qbin[0][ 7] =   0.00602123;
  eff_corr_ks_9qbin[0][ 8] =     0.997588; eff_corr_ksE_9qbin[0][ 8] =   0.00602803;
  eff_corr_ks_9qbin[0][ 9] =     0.997371; eff_corr_ksE_9qbin[0][ 9] =   0.00602483;
  eff_corr_ks_9qbin[0][10] =      0.99798; eff_corr_ksE_9qbin[0][10] =   0.00603558;
  eff_corr_ks_9qbin[0][11] =     0.997952; eff_corr_ksE_9qbin[0][11] =   0.00603451;
  eff_corr_ks_9qbin[0][12] =     0.999791; eff_corr_ksE_9qbin[0][12] =   0.00600512;
  eff_corr_ks_9qbin[0][13] =     0.999822; eff_corr_ksE_9qbin[0][13] =   0.00601276;
  eff_corr_ks_9qbin[0][14] =     0.997627; eff_corr_ksE_9qbin[0][14] =   0.00601673;
  eff_corr_ks_9qbin[1][ 0] =     0.997509; eff_corr_ksE_9qbin[1][ 0] =   0.00604053;
  eff_corr_ks_9qbin[1][ 1] =     0.997577; eff_corr_ksE_9qbin[1][ 1] =   0.00603974;
  eff_corr_ks_9qbin[1][ 2] =     0.997922; eff_corr_ksE_9qbin[1][ 2] =   0.00602011;
  eff_corr_ks_9qbin[1][ 3] =      0.99792; eff_corr_ksE_9qbin[1][ 3] =   0.00602174;
  eff_corr_ks_9qbin[1][ 4] =     0.997701; eff_corr_ksE_9qbin[1][ 4] =   0.00601756;
  eff_corr_ks_9qbin[1][ 5] =     0.997689; eff_corr_ksE_9qbin[1][ 5] =   0.00601746;
  eff_corr_ks_9qbin[1][ 6] =     0.997221; eff_corr_ksE_9qbin[1][ 6] =   0.00601953;
  eff_corr_ks_9qbin[1][ 7] =     0.997281; eff_corr_ksE_9qbin[1][ 7] =   0.00602099;
  eff_corr_ks_9qbin[1][ 8] =      0.99737; eff_corr_ksE_9qbin[1][ 8] =   0.00602823;
  eff_corr_ks_9qbin[1][ 9] =     0.997481; eff_corr_ksE_9qbin[1][ 9] =   0.00602674;
  eff_corr_ks_9qbin[1][10] =     0.998083; eff_corr_ksE_9qbin[1][10] =   0.00603116;
  eff_corr_ks_9qbin[1][11] =     0.998035; eff_corr_ksE_9qbin[1][11] =   0.00603477;
  eff_corr_ks_9qbin[1][12] =     0.999822; eff_corr_ksE_9qbin[1][12] =   0.00600552;
  eff_corr_ks_9qbin[1][13] =     0.999792; eff_corr_ksE_9qbin[1][13] =   0.00601015;
  eff_corr_ks_9qbin[1][14] =     0.997733; eff_corr_ksE_9qbin[1][14] =   0.00601781;
  eff_corr_ks_9qbin[2][ 0] =     0.997477; eff_corr_ksE_9qbin[2][ 0] =    0.0060392;
  eff_corr_ks_9qbin[2][ 1] =     0.997512; eff_corr_ksE_9qbin[2][ 1] =   0.00603779;
  eff_corr_ks_9qbin[2][ 2] =     0.997876; eff_corr_ksE_9qbin[2][ 2] =    0.0060204;
  eff_corr_ks_9qbin[2][ 3] =     0.997883; eff_corr_ksE_9qbin[2][ 3] =   0.00602115;
  eff_corr_ks_9qbin[2][ 4] =     0.997647; eff_corr_ksE_9qbin[2][ 4] =   0.00601746;
  eff_corr_ks_9qbin[2][ 5] =     0.997601; eff_corr_ksE_9qbin[2][ 5] =    0.0060185;
  eff_corr_ks_9qbin[2][ 6] =     0.997245; eff_corr_ksE_9qbin[2][ 6] =    0.0060206;
  eff_corr_ks_9qbin[2][ 7] =     0.997291; eff_corr_ksE_9qbin[2][ 7] =   0.00602116;
  eff_corr_ks_9qbin[2][ 8] =     0.997494; eff_corr_ksE_9qbin[2][ 8] =   0.00602812;
  eff_corr_ks_9qbin[2][ 9] =     0.997418; eff_corr_ksE_9qbin[2][ 9] =   0.00602563;
  eff_corr_ks_9qbin[2][10] =     0.998026; eff_corr_ksE_9qbin[2][10] =   0.00603359;
  eff_corr_ks_9qbin[2][11] =     0.997988; eff_corr_ksE_9qbin[2][11] =   0.00603462;
  eff_corr_ks_9qbin[2][12] =     0.999805; eff_corr_ksE_9qbin[2][12] =    0.0060053;
  eff_corr_ks_9qbin[2][13] =     0.999809; eff_corr_ksE_9qbin[2][13] =   0.00601158;
  eff_corr_ks_9qbin[2][14] =     0.997676; eff_corr_ksE_9qbin[2][14] =   0.00601722;

  Double_t eff_corr_ks_6qbin [3][9] = {0}; // [mm,ee,ee+mm]
  Double_t eff_corr_ksE_6qbin[3][9] = {0}; // [mm,ee,ee+mm]
  eff_corr_ks_6qbin[0][ 0] =     0.997614; eff_corr_ksE_6qbin[0][ 0] =   0.00602834;
  eff_corr_ks_6qbin[0][ 1] =     0.997615; eff_corr_ksE_6qbin[0][ 1] =   0.00602724;
  eff_corr_ks_6qbin[0][ 2] =     0.997604; eff_corr_ksE_6qbin[0][ 2] =   0.00601739;
  eff_corr_ks_6qbin[0][ 3] =     0.997532; eff_corr_ksE_6qbin[0][ 3] =   0.00601935;
  eff_corr_ks_6qbin[0][ 4] =     0.997256; eff_corr_ksE_6qbin[0][ 4] =   0.00602106;
  eff_corr_ks_6qbin[0][ 5] =     0.997295; eff_corr_ksE_6qbin[0][ 5] =   0.00602123;
  eff_corr_ks_6qbin[0][ 6] =      0.99798; eff_corr_ksE_6qbin[0][ 6] =   0.00602403;
  eff_corr_ks_6qbin[0][ 7] =     0.997851; eff_corr_ksE_6qbin[0][ 7] =   0.00602245;
  eff_corr_ks_6qbin[0][ 8] =     0.997627; eff_corr_ksE_6qbin[0][ 8] =   0.00601673;
  eff_corr_ks_6qbin[1][ 0] =     0.997687; eff_corr_ksE_6qbin[1][ 0] =   0.00603029;
  eff_corr_ks_6qbin[1][ 1] =     0.997727; eff_corr_ksE_6qbin[1][ 1] =   0.00603062;
  eff_corr_ks_6qbin[1][ 2] =     0.997701; eff_corr_ksE_6qbin[1][ 2] =   0.00601756;
  eff_corr_ks_6qbin[1][ 3] =     0.997689; eff_corr_ksE_6qbin[1][ 3] =   0.00601746;
  eff_corr_ks_6qbin[1][ 4] =     0.997221; eff_corr_ksE_6qbin[1][ 4] =   0.00601953;
  eff_corr_ks_6qbin[1][ 5] =     0.997281; eff_corr_ksE_6qbin[1][ 5] =   0.00602099;
  eff_corr_ks_6qbin[1][ 6] =     0.997948; eff_corr_ksE_6qbin[1][ 6] =   0.00602284;
  eff_corr_ks_6qbin[1][ 7] =     0.997956; eff_corr_ksE_6qbin[1][ 7] =   0.00602347;
  eff_corr_ks_6qbin[1][ 8] =     0.997733; eff_corr_ksE_6qbin[1][ 8] =   0.00601781;
  eff_corr_ks_6qbin[2][ 0] =     0.997656; eff_corr_ksE_6qbin[2][ 0] =   0.00602944;
  eff_corr_ks_6qbin[2][ 1] =     0.997679; eff_corr_ksE_6qbin[2][ 1] =   0.00602914;
  eff_corr_ks_6qbin[2][ 2] =     0.997647; eff_corr_ksE_6qbin[2][ 2] =   0.00601746;
  eff_corr_ks_6qbin[2][ 3] =     0.997601; eff_corr_ksE_6qbin[2][ 3] =    0.0060185;
  eff_corr_ks_6qbin[2][ 4] =     0.997245; eff_corr_ksE_6qbin[2][ 4] =    0.0060206;
  eff_corr_ks_6qbin[2][ 5] =     0.997291; eff_corr_ksE_6qbin[2][ 5] =   0.00602116;
  eff_corr_ks_6qbin[2][ 6] =     0.997966; eff_corr_ksE_6qbin[2][ 6] =   0.00602351;
  eff_corr_ks_6qbin[2][ 7] =     0.997897; eff_corr_ksE_6qbin[2][ 7] =   0.00602289;
  eff_corr_ks_6qbin[2][ 8] =     0.997676; eff_corr_ksE_6qbin[2][ 8] =   0.00601722;

 //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Double_t eff_corr_pi0_9qbin [3][15] = {0}; // [mm,ee,ee+mm]
  Double_t eff_corr_pi0E_9qbin[3][15] = {0}; // [mm,ee,ee+mm]
  eff_corr_pi0_9qbin[0][ 0] =     0.991626; eff_corr_pi0E_9qbin[0][ 0] =   0.00305709;
  eff_corr_pi0_9qbin[0][ 1] =     0.991316; eff_corr_pi0E_9qbin[0][ 1] =   0.00310706;
  eff_corr_pi0_9qbin[0][ 2] =     0.993613; eff_corr_pi0E_9qbin[0][ 2] =   0.00235197;
  eff_corr_pi0_9qbin[0][ 3] =     0.993415; eff_corr_pi0E_9qbin[0][ 3] =    0.0024183;
  eff_corr_pi0_9qbin[0][ 4] =     0.993766; eff_corr_pi0E_9qbin[0][ 4] =    0.0022499;
  eff_corr_pi0_9qbin[0][ 5] =     0.994154; eff_corr_pi0E_9qbin[0][ 5] =   0.00216309;
  eff_corr_pi0_9qbin[0][ 6] =     0.994423; eff_corr_pi0E_9qbin[0][ 6] =   0.00204091;
  eff_corr_pi0_9qbin[0][ 7] =     0.994541; eff_corr_pi0E_9qbin[0][ 7] =   0.00201946;
  eff_corr_pi0_9qbin[0][ 8] =     0.996604; eff_corr_pi0E_9qbin[0][ 8] =    0.0014777;
  eff_corr_pi0_9qbin[0][ 9] =     0.996657; eff_corr_pi0E_9qbin[0][ 9] =   0.00146865;
  eff_corr_pi0_9qbin[0][10] =     0.998614; eff_corr_pi0E_9qbin[0][10] =   0.00090661;
  eff_corr_pi0_9qbin[0][11] =     0.998784; eff_corr_pi0E_9qbin[0][11] =  0.000831859;
  eff_corr_pi0_9qbin[0][12] =      1.00002; eff_corr_pi0E_9qbin[0][12] =  4.59449e-05;
  eff_corr_pi0_9qbin[0][13] =     0.999978; eff_corr_pi0E_9qbin[0][13] =  4.73056e-05;
  eff_corr_pi0_9qbin[0][14] =     0.994739; eff_corr_pi0E_9qbin[0][14] =   0.00199785;
  eff_corr_pi0_9qbin[1][ 0] =     0.992409; eff_corr_pi0E_9qbin[1][ 0] =   0.00275292;
  eff_corr_pi0_9qbin[1][ 1] =     0.992194; eff_corr_pi0E_9qbin[1][ 1] =   0.00280218;
  eff_corr_pi0_9qbin[1][ 2] =     0.994505; eff_corr_pi0E_9qbin[1][ 2] =   0.00201682;
  eff_corr_pi0_9qbin[1][ 3] =     0.994693; eff_corr_pi0E_9qbin[1][ 3] =   0.00195701;
  eff_corr_pi0_9qbin[1][ 4] =     0.994814; eff_corr_pi0E_9qbin[1][ 4] =   0.00185837;
  eff_corr_pi0_9qbin[1][ 5] =     0.994989; eff_corr_pi0E_9qbin[1][ 5] =   0.00178091;
  eff_corr_pi0_9qbin[1][ 6] =     0.994744; eff_corr_pi0E_9qbin[1][ 6] =   0.00194468;
  eff_corr_pi0_9qbin[1][ 7] =     0.995288; eff_corr_pi0E_9qbin[1][ 7] =   0.00178468;
  eff_corr_pi0_9qbin[1][ 8] =     0.996612; eff_corr_pi0E_9qbin[1][ 8] =   0.00144686;
  eff_corr_pi0_9qbin[1][ 9] =     0.996666; eff_corr_pi0E_9qbin[1][ 9] =   0.00144439;
  eff_corr_pi0_9qbin[1][10] =     0.998593; eff_corr_pi0E_9qbin[1][10] =   0.00085778;
  eff_corr_pi0_9qbin[1][11] =      0.99866; eff_corr_pi0E_9qbin[1][11] =  0.000854885;
  eff_corr_pi0_9qbin[1][12] =      1.00001; eff_corr_pi0E_9qbin[1][12] =    3.524e-05;
  eff_corr_pi0_9qbin[1][13] =      1.00002; eff_corr_pi0E_9qbin[1][13] =  4.86258e-05;
  eff_corr_pi0_9qbin[1][14] =     0.994999; eff_corr_pi0E_9qbin[1][14] =   0.00187765;
  eff_corr_pi0_9qbin[2][ 0] =     0.992085; eff_corr_pi0E_9qbin[2][ 0] =   0.00287891;
  eff_corr_pi0_9qbin[2][ 1] =     0.991828; eff_corr_pi0E_9qbin[2][ 1] =   0.00292918;
  eff_corr_pi0_9qbin[2][ 2] =     0.994101; eff_corr_pi0E_9qbin[2][ 2] =   0.00216852;
  eff_corr_pi0_9qbin[2][ 3] =      0.99412; eff_corr_pi0E_9qbin[2][ 3] =   0.00216364;
  eff_corr_pi0_9qbin[2][ 4] =     0.994229; eff_corr_pi0E_9qbin[2][ 4] =   0.00207705;
  eff_corr_pi0_9qbin[2][ 5] =     0.994525; eff_corr_pi0E_9qbin[2][ 5] =   0.00199358;
  eff_corr_pi0_9qbin[2][ 6] =     0.994519; eff_corr_pi0E_9qbin[2][ 6] =   0.00201236;
  eff_corr_pi0_9qbin[2][ 7] =     0.994755; eff_corr_pi0E_9qbin[2][ 7] =   0.00195217;
  eff_corr_pi0_9qbin[2][ 8] =     0.996607; eff_corr_pi0E_9qbin[2][ 8] =   0.00146448;
  eff_corr_pi0_9qbin[2][ 9] =     0.996661; eff_corr_pi0E_9qbin[2][ 9] =   0.00145828;
  eff_corr_pi0_9qbin[2][10] =     0.998605; eff_corr_pi0E_9qbin[2][10] =  0.000884957;
  eff_corr_pi0_9qbin[2][11] =     0.998729; eff_corr_pi0E_9qbin[2][11] =     0.000842;
  eff_corr_pi0_9qbin[2][12] =      1.00002; eff_corr_pi0E_9qbin[2][12] =  4.12464e-05;
  eff_corr_pi0_9qbin[2][13] =     0.999996; eff_corr_pi0E_9qbin[2][13] =  4.78834e-05;
  eff_corr_pi0_9qbin[2][14] =     0.994859; eff_corr_pi0E_9qbin[2][14] =   0.00194254;

  Double_t eff_corr_pi0_6qbin [3][9] = {0}; // [mm,ee,ee+mm]
  Double_t eff_corr_pi0E_6qbin[3][9] = {0}; // [mm,ee,ee+mm]
  eff_corr_pi0_6qbin[0][ 0] =     0.992561; eff_corr_pi0E_6qbin[0][ 0] =   0.00272515;
  eff_corr_pi0_6qbin[0][ 1] =     0.992298; eff_corr_pi0E_6qbin[0][ 1] =   0.00278473;
  eff_corr_pi0_6qbin[0][ 2] =     0.993766; eff_corr_pi0E_6qbin[0][ 2] =    0.0022499;
  eff_corr_pi0_6qbin[0][ 3] =     0.994154; eff_corr_pi0E_6qbin[0][ 3] =   0.00216309;
  eff_corr_pi0_6qbin[0][ 4] =     0.994423; eff_corr_pi0E_6qbin[0][ 4] =   0.00204091;
  eff_corr_pi0_6qbin[0][ 5] =     0.994541; eff_corr_pi0E_6qbin[0][ 5] =   0.00201946;
  eff_corr_pi0_6qbin[0][ 6] =     0.997903; eff_corr_pi0E_6qbin[0][ 6] =   0.00105846;
  eff_corr_pi0_6qbin[0][ 7] =     0.997986; eff_corr_pi0E_6qbin[0][ 7] =   0.00102514;
  eff_corr_pi0_6qbin[0][ 8] =     0.994739; eff_corr_pi0E_6qbin[0][ 8] =   0.00199785;
  eff_corr_pi0_6qbin[1][ 0] =     0.993315; eff_corr_pi0E_6qbin[1][ 0] =   0.00243491;
  eff_corr_pi0_6qbin[1][ 1] =     0.993284; eff_corr_pi0E_6qbin[1][ 1] =   0.00243342;
  eff_corr_pi0_6qbin[1][ 2] =     0.994814; eff_corr_pi0E_6qbin[1][ 2] =   0.00185837;
  eff_corr_pi0_6qbin[1][ 3] =     0.994989; eff_corr_pi0E_6qbin[1][ 3] =   0.00178091;
  eff_corr_pi0_6qbin[1][ 4] =     0.994744; eff_corr_pi0E_6qbin[1][ 4] =   0.00194468;
  eff_corr_pi0_6qbin[1][ 5] =     0.995288; eff_corr_pi0E_6qbin[1][ 5] =   0.00178468;
  eff_corr_pi0_6qbin[1][ 6] =     0.997926; eff_corr_pi0E_6qbin[1][ 6] =   0.00101242;
  eff_corr_pi0_6qbin[1][ 7] =     0.997962; eff_corr_pi0E_6qbin[1][ 7] =   0.00101763;
  eff_corr_pi0_6qbin[1][ 8] =     0.994999; eff_corr_pi0E_6qbin[1][ 8] =   0.00187765;
  eff_corr_pi0_6qbin[2][ 0] =      0.99299; eff_corr_pi0E_6qbin[2][ 0] =   0.00256013;
  eff_corr_pi0_6qbin[2][ 1] =     0.992859; eff_corr_pi0E_6qbin[2][ 1] =   0.00258472;
  eff_corr_pi0_6qbin[2][ 2] =     0.994229; eff_corr_pi0E_6qbin[2][ 2] =   0.00207705;
  eff_corr_pi0_6qbin[2][ 3] =     0.994525; eff_corr_pi0E_6qbin[2][ 3] =   0.00199358;
  eff_corr_pi0_6qbin[2][ 4] =     0.994519; eff_corr_pi0E_6qbin[2][ 4] =   0.00201236;
  eff_corr_pi0_6qbin[2][ 5] =     0.994755; eff_corr_pi0E_6qbin[2][ 5] =   0.00195217;
  eff_corr_pi0_6qbin[2][ 6] =     0.997913; eff_corr_pi0E_6qbin[2][ 6] =   0.00103835;
  eff_corr_pi0_6qbin[2][ 7] =     0.997976; eff_corr_pi0E_6qbin[2][ 7] =   0.00102188;
  eff_corr_pi0_6qbin[2][ 8] =     0.994859; eff_corr_pi0E_6qbin[2][ 8] =   0.00194254;
 
 //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Double_t eff_corr_tot_6qbin [3][9] = {0}; // [mm,ee,ee+mm]
  Double_t eff_corr_totE_6qbin[3][9] = {0}; // [mm,ee,ee+mm]
  Double_t eff_corr_tot_9qbin [3][15] = {0}; // [mm,ee,ee+mm]
  Double_t eff_corr_totE_9qbin[3][15] = {0}; // [mm,ee,ee+mm]
  for( Int_t i=0; i<3; i++ ){
    for( Int_t j=0; j<15; j++ ){
      eff_corr_tot_9qbin [i][j] = eff_corr_lep_9qbin[i][j]*eff_corr_kpi_9qbin[i][j]*eff_corr_ks_9qbin[i][j]*eff_corr_pi0_9qbin[i][j];
      eff_corr_totE_9qbin[i][j] = eff_corr_tot_9qbin[i][j]*sqrt(
								eff_corr_lepE_9qbin [i][j]*eff_corr_lepE_9qbin[i][j]/eff_corr_lep_9qbin[i][j]/eff_corr_lep_9qbin[i][j]
								+eff_corr_kpiE_9qbin[i][j]*eff_corr_kpiE_9qbin[i][j]/eff_corr_kpi_9qbin[i][j]/eff_corr_kpi_9qbin[i][j]
								+eff_corr_ksE_9qbin [i][j]*eff_corr_ksE_9qbin [i][j]/eff_corr_ks_9qbin [i][j]/eff_corr_ks_9qbin [i][j]
								+eff_corr_pi0E_9qbin[i][j]*eff_corr_pi0E_9qbin[i][j]/eff_corr_pi0_9qbin[i][j]/eff_corr_pi0_9qbin[i][j]
								);
    }
    for( Int_t j=0; j<9; j++ ){
      eff_corr_tot_6qbin [i][j] = eff_corr_lep_6qbin[i][j]*eff_corr_kpi_6qbin[i][j]*eff_corr_ks_6qbin[i][j]*eff_corr_pi0_6qbin[i][j];
      eff_corr_totE_6qbin[i][j] = eff_corr_tot_6qbin[i][j]*sqrt(
								eff_corr_lepE_6qbin [i][j]*eff_corr_lepE_6qbin[i][j]/eff_corr_lep_6qbin[i][j]/eff_corr_lep_6qbin[i][j]
								+eff_corr_kpiE_6qbin[i][j]*eff_corr_kpiE_6qbin[i][j]/eff_corr_kpi_6qbin[i][j]/eff_corr_kpi_6qbin[i][j]
								+eff_corr_ksE_6qbin [i][j]*eff_corr_ksE_6qbin [i][j]/eff_corr_ks_6qbin [i][j]/eff_corr_ks_6qbin [i][j]
								+eff_corr_pi0E_6qbin[i][j]*eff_corr_pi0E_6qbin[i][j]/eff_corr_pi0_6qbin[i][j]/eff_corr_pi0_6qbin[i][j]
								);
    }
  }

 //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  MChain**  chain            = new MChain*[Nchain];
  TH2D**    tmphist          = new TH2D*  [Ntmp  ];
  TH2D**    hist             = new TH2D*  [Nchain]; // q2-cos distribution before reconstruction efficiency correction
  TH2D**    hist_eff         = new TH2D*  [Nchain]; //                     after
  TH1D**    histx_gen        = new TH1D*  [Nchain]; //   q2   distribution before reconstruction efficiency correction
  TH1D**    histx_gen_eff    = new TH1D*  [Nchain]; //                     after
  TH1D**    histy_gen        = new TH1D*  [Nchain]; //   cos  distribution before reconstruction efficiency correction
  TH1D**    histy_gen_eff    = new TH1D*  [Nchain]; //                     after
  TH1D**    hist_afb         = new TH1D*  [Nchain]; //   afb  distribution before reconstruction efficiency correction
  TH1D**    hist_afb_eff     = new TH1D*  [Nchain]; //                     after
  TH1D**    hist_afb_bin     = new TH1D*  [Nchain]; //                     before                                      in afb binning
  TH1D**    hist_afb_eff_bin = new TH1D*  [Nchain]; //                     after                                       in afb binning
  TCanvas*  c1               = Canvas( "c1","c1",3, 2 );

  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++

  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form(" ************************ make tmphist ( %s, %s ) *************************************",tname[1],axis[1]) << std::endl;
    std::cout << Form( "<infile %d > ", j );
    chain[j] = new MChain( infile[j], tname[1], branch_table(), nfile[j], "*.root" );
    nominal_cut_selection( chain[j], fl_mode_ll )( chain[j]->GetCut(), tname[1] );
  }
  
  // ++++++++++++++++++++++++
  // cut change
  //
  //chain[0]->GetCut()->Set(    441, 0 );
  //chain[0]->GetCut()->Set(    443, 0 );
  //chain[0]->GetCut()->Set( 100443, 0 );
  //
  // ++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ changed cut *************************************" << std::endl;
  for( Int_t j=0; j<Nchain; j++ ){
    std::cout << Form( "<infile %d > ", j ) << std::endl;
    chain[j]->GetCut()->Display(0);
    chain[j]->MakeTree();
  }
  
  // +++++++ Open reconstruction-efficiency data ++++++++++++++++++++++
  TFile* file_eff_lep0  = new TFile( "eff_table/mb495/2d_q2_theta_lep0_setA-U.root" );
  TFile* file_eff_lep1  = new TFile( "eff_table/mb495/2d_q2_theta_lep1_setA-U.root" );
  if( file_eff_lep0->IsZombie() ) std::cerr << "[ABORT] can not find input-file : " << file_eff_lep0->GetName() << std::endl, abort();
  if( file_eff_lep1->IsZombie() ) std::cerr << "[ABORT] can not find input-file : " << file_eff_lep1->GetName() << std::endl, abort();
  TH2D** hist_rec_eff  = new TH2D*[2];
  hist_rec_eff[0] = (TH2D*)file_eff_lep0->Get( "h511_h12" ); // mm
  hist_rec_eff[1] = (TH2D*)file_eff_lep1->Get( "h511_h12" ); // ee
  if( hist_rec_eff[0] == NULL ) std::cerr << "[ABORT] can not find histgram(mm) : " << std::endl, abort();
  if( hist_rec_eff[1] == NULL ) std::cerr << "[ABORT] can not find histgram(ee) : " << std::endl, abort();
  
  // +++++++ make tmphist ++++++++++++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ tmphist *************************************" << std::endl;
  for( Int_t j=0; j<Ntmp; j++ ){
    tmphist[j] = new TH2D( Form("tmphist%d",j), Form("%s",chain[0]->GetChange()), xbin,xbins[fl_mode_ll], ybin,ymin,ymax );
    chain[0]->GetTree()->Project( Form("tmphist%d",j), axis[1], add_cut[j] );
    tmphist[j]->Sumw2();
    std::cout << "add_cut : " << add_cut[j] << std::endl;
  }
  
  // +++++++ make hist  ++++++++++++++++++++++++++++++++++
  for( Int_t i=0; i<Nchain; i++ ){
    hist[i] = new TH2D( Form("hist%d",i),Form("hist%d",i), xbin,xbins[fl_mode_ll], ybin,ymin,ymax );
    hist[i]->GetXaxis()->CenterTitle();
    hist[i]->GetYaxis()->CenterTitle();
    hist[i]->SetXTitle( xlabel );
    hist[i]->SetYTitle( ylabel );

    std::cout << "Added-Files( ";
    for( Int_t j=0; j<Ntmp; j++ ){
      if( add[i][j] ) std::cout << j <<  ",", hist[i]->Add( tmphist[j], add[i][j] );
    }
    std::cout << ")" << std::endl;
    // +++++++ display ++++++++++++++++++++++++++++++++++
    std::cout << std::setw(10) << Form("<hist%d>",i);
    std::cout << std::setw(9)  << hist[i]->GetEntries( ) << " events(GetEntries)"
	      << " ( "
	      << std::setw(9)  << hist[i]->Integral()    << " events(Integral)"
	      << " ) "
	      << std::endl; 
  }

  // +++++++ correct PID systematics  ++++++++++++++++++++
  if( fl_corr ){
    for( Int_t iy=1; iy<=ybin; iy++ ){
      for( Int_t ix=1; ix<=xbin; ix++ ){
	if( xbin_afb==9 ){
	  Int_t tmp_q2 = 0;
	  if( fl_mode_ll==0 ){
	    if     ( ix <  3 ) tmp_q2 =  0; // 0.0 ~ 2.0 GeV^2                // 1,2
	    else if( ix <  5 ) tmp_q2 =  1; // 2.0 ~ 4.3 GeV^2                // 3,4
	    else if( ix <  9 ) tmp_q2 =  2; // 4.3 ~ J/psi(low edge)          // 5,6,7,8
	    else if( ix < 10 ) tmp_q2 = -1; // J/psi veto                     // 9
	    else if( ix < 12 ) tmp_q2 =  3; // between J/psi and psi(2S)      // 10,11
	    else if( ix < 13 ) tmp_q2 = -1; // psi(2S) veto                   // 12
	    else if( ix < 15 ) tmp_q2 =  4; // psi(2S)(low edge) ~ 16.0 GeV^2 // 13,14
	    else if( ix < 18 ) tmp_q2 =  5; // 16.0 ~ 19.0 GeV^2              // 15,16,17
	    else if( ix < 22 ) tmp_q2 =  6; // 16.0 ~ 19.0 GeV^2              // 18,19,20,21
	  }else{
	    if     ( ix <  3 ) tmp_q2 =  0; // 0.0 ~ 2.0 GeV^2                // 1,2
	    else if( ix <  5 ) tmp_q2 =  1; // 2.0 ~ 4.3 GeV^2                // 3,4
	    else if( ix <  8 ) tmp_q2 =  2; // 4.3 ~ J/psi(low edge)          // 5,6,7,
	    else if( ix < 10 ) tmp_q2 = -1; // J/psi veto                     // 8,9
	    else if( ix < 12 ) tmp_q2 =  3; // between J/psi and psi(2S)      // 10,11
	    else if( ix < 13 ) tmp_q2 = -1; // psi(2S) veto                   // 12
	    else if( ix < 15 ) tmp_q2 =  4; // psi(2S)(low edge) ~ 16.0 GeV^2 // 13,14
	    else if( ix < 18 ) tmp_q2 =  5; // 16.0 ~ 19.0 GeV^2              // 15,16,17
	    else if( ix < 22 ) tmp_q2 =  6; // 16.0 ~ 19.0 GeV^2              // 18,19,20,21
	  }
	  if( tmp_q2>=0 ){
	    Int_t    tmp_cosq = iy > ybin/2 ? 1 : 0;
	    Int_t    tmp_bin  = 2*tmp_q2 + tmp_cosq;
	    Double_t tmp_eff  = hist_rec_eff[fl_mode_ll]->GetBinContent(ix,iy);
	    Double_t tmp_effE = hist_rec_eff[fl_mode_ll]->GetBinError  (ix,iy);
	    Double_t tmp_var  = tmp_eff*eff_corr_tot_9qbin[fl_mode_ll][tmp_bin];
	    Double_t tmp_varE = tmp_var*sqrt(
					     tmp_effE*tmp_effE/tmp_eff/tmp_eff
					     +eff_corr_totE_9qbin[fl_mode_ll][tmp_bin]*eff_corr_totE_9qbin[fl_mode_ll][tmp_bin]/eff_corr_totE_9qbin[fl_mode_ll][tmp_bin]/eff_corr_totE_9qbin[fl_mode_ll][tmp_bin]
					     );
	    if( tmp_effE==0 ) tmp_varE = 0;
	    //std::cout << Form("[ix=%d, iy=%d, fl_mode_ll=%d,q2=%d,cosq=%d,bin=%d] (%f +- %f) * (%f +- %f) = (%f +- %f)", ix,iy,fl_mode_ll,tmp_q2,tmp_cosq,tmp_bin,tmp_eff,tmp_effE,eff_corr_tot_9qbin[fl_mode_ll][tmp_bin],eff_corr_totE_9qbin[fl_mode_ll][tmp_bin],tmp_var, tmp_varE);
	    hist_rec_eff[fl_mode_ll]->SetBinContent( ix,iy, tmp_var  );
	    hist_rec_eff[fl_mode_ll]->SetBinError  ( ix,iy, tmp_varE );
	    //std::cout << Form( " = (%f +- %f)", hist_rec_eff[fl_mode_ll]->GetBinContent(ix,iy),hist_rec_eff[fl_mode_ll]->GetBinError(ix,iy) ) << std::endl;
	  }
	}else if( xbin_afb==6 ){
	  Int_t tmp_q2 = 0;
	  if( fl_mode_ll==0 ){
	    if     ( ix <  3 ) tmp_q2 =  0; // 0.0 ~ 2.0 GeV^2                // 1,2
	    else if( ix <  5 ) tmp_q2 =  0; // 2.0 ~ 4.3 GeV^2                // 3,4
	    else if( ix <  9 ) tmp_q2 =  1; // 4.3 ~ J/psi(low edge)          // 5,6,7,8
	    else if( ix < 10 ) tmp_q2 = -1; // J/psi veto                     // 9
	    else if( ix < 12 ) tmp_q2 =  2; // between J/psi and psi(2S)      // 10,11
	    else if( ix < 13 ) tmp_q2 = -1; // psi(2S) veto                   // 12
	    else if( ix < 15 ) tmp_q2 =  3; // psi(2S)(low edge) ~ 16.0 GeV^2 // 13,14
	    else if( ix < 18 ) tmp_q2 =  3; // 16.0 ~ 19.0 GeV^2              // 15,16,17
	    else if( ix < 22 ) tmp_q2 =  3; // 16.0 ~ 19.0 GeV^2              // 18,19,20,21
	  }else{
	    if     ( ix <  3 ) tmp_q2 =  0; // 0.0 ~ 2.0 GeV^2                // 1,2
	    else if( ix <  5 ) tmp_q2 =  0; // 2.0 ~ 4.3 GeV^2                // 3,4
	    else if( ix <  8 ) tmp_q2 =  1; // 4.3 ~ J/psi(low edge)          // 5,6,7,
	    else if( ix < 10 ) tmp_q2 = -1; // J/psi veto                     // 8,9
	    else if( ix < 12 ) tmp_q2 =  2; // between J/psi and psi(2S)      // 10,11
	    else if( ix < 13 ) tmp_q2 = -1; // psi(2S) veto                   // 12
	    else if( ix < 15 ) tmp_q2 =  3; // psi(2S)(low edge) ~ 16.0 GeV^2 // 13,14
	    else if( ix < 18 ) tmp_q2 =  3; // 16.0 ~ 19.0 GeV^2              // 15,16,17
	    else if( ix < 22 ) tmp_q2 =  3; // 16.0 ~ 19.0 GeV^2              // 18,19,20,21
	  }
	  if( tmp_q2>=0 ){
	    Int_t    tmp_cosq = iy > ybin/2 ? 1 : 0;
	    Int_t    tmp_bin  = 2*tmp_q2 + tmp_cosq;
	    Double_t tmp_eff  = hist_rec_eff[fl_mode_ll]->GetBinContent(ix,iy);
	    Double_t tmp_effE = hist_rec_eff[fl_mode_ll]->GetBinError  (ix,iy);
	    Double_t tmp_var  = tmp_eff*eff_corr_tot_6qbin[fl_mode_ll][tmp_bin];
	    Double_t tmp_varE = tmp_var*sqrt(
					     tmp_effE*tmp_effE/tmp_eff/tmp_eff
					     +eff_corr_totE_6qbin[fl_mode_ll][tmp_bin]*eff_corr_totE_6qbin[fl_mode_ll][tmp_bin]/eff_corr_totE_6qbin[fl_mode_ll][tmp_bin]/eff_corr_totE_6qbin[fl_mode_ll][tmp_bin]
					     );
	    if( tmp_effE==0 ) tmp_varE = 0;
	    //std::cout << Form("[ix=%d, iy=%d, fl_mode_ll=%d,q2=%d,cosq=%d,bin=%d] (%f +- %f) * (%f +- %f) = (%f +- %f)", ix,iy,fl_mode_ll,tmp_q2,tmp_cosq,tmp_bin,tmp_eff,tmp_effE,eff_corr_tot_6qbin[fl_mode_ll][tmp_bin],eff_corr_totE_6qbin[fl_mode_ll][tmp_bin],tmp_var, tmp_varE);
	    hist_rec_eff[fl_mode_ll]->SetBinContent( ix,iy, tmp_var  );
	    hist_rec_eff[fl_mode_ll]->SetBinError  ( ix,iy, tmp_varE );
	    //std::cout << Form( " = (%f +- %f)", hist_rec_eff[fl_mode_ll]->GetBinContent(ix,iy),hist_rec_eff[fl_mode_ll]->GetBinError(ix,iy) ) << std::endl;
	  }
	}else std::cerr << "[ABORT] Wrong xbin_afb : " << xbin_afb << std::endl, abort();
	
      }
    }
  }
  // +++++++ make hist after rec. eff. correction  ++++++++++++++++++++
  for( Int_t i=0; i<Nchain; i++ ) hist_eff[i] = multiply( hist[i], hist_rec_eff[fl_mode_ll], 0, Form("hist_corr%d",i), Form("hist_corr%d",i) );
  // +++++++ projection hist ++++++++++++++++++
  for( Int_t i=0; i<Nchain; i++ ){
    histx_gen    [i] = new TH1D(*(hist    [i]->ProjectionX( Form("_px%d",i), 0, -1, "e")) ); Deco( histx_gen[i],     0, 2, 2 );
    histy_gen    [i] = new TH1D(*(hist    [i]->ProjectionY( Form("_py%d",i), 0, -1, "e")) ); Deco( histy_gen[i],     0, 2, 2 );
    histx_gen_eff[i] = new TH1D(*(hist_eff[i]->ProjectionX( Form("_px%d",i), 0, -1, "e")) ); Deco( histx_gen_eff[i], 0, 1, 1 );
    histy_gen_eff[i] = new TH1D(*(hist_eff[i]->ProjectionY( Form("_py%d",i), 0, -1, "e")) ); Deco( histy_gen_eff[i], 0, 1, 1 );
  }

  // +++++++ make hist(afb) ++++++++++++++++++
  for( Int_t i=0; i<Nchain; i++ ){
    hist_afb        [i] = new TH1D( Form("afb_rec%d",        i), Form("afb_rec%d",        i), xbin,     xbins    [fl_mode_ll]);
    hist_afb_eff    [i] = new TH1D( Form("afb_rec_eff%d",    i), Form("afb_rec_eff%d",    i), xbin,     xbins    [fl_mode_ll]);
    hist_afb_bin    [i] = new TH1D( Form("afb_rec_bin%d",    i), Form("afb_rec_bin%d",    i), xbin_afb, xbins_afb[fl_mode_ll]);
    hist_afb_eff_bin[i] = new TH1D( Form("afb_rec_eff_bin%d",i), Form("afb_rec_eff_bin%d",i), xbin_afb, xbins_afb[fl_mode_ll]);
    hist_afb        [i]->SetXTitle( xlabel   ); hist_afb        [i]->SetYTitle( "A_{FB}" ); Deco( hist_afb        [i], 0, 2, 2 );
    hist_afb_bin    [i]->SetXTitle( xlabel   ); hist_afb_bin    [i]->SetYTitle( "A_{FB}" ); Deco( hist_afb_bin    [i], 0, 2, 2 );
    hist_afb_eff    [i]->SetXTitle( xlabel   ); hist_afb_eff    [i]->SetYTitle( "A_{FB}" ); Deco( hist_afb_eff    [i], 0, 1, 1 );
    hist_afb_eff_bin[i]->SetXTitle( xlabel   ); hist_afb_eff_bin[i]->SetYTitle( "A_{FB}" ); Deco( hist_afb_eff_bin[i], 0, 1, 1 );

    // for afb(gen)
    std::cout << "[Entry of AFB(gen)]" << std::endl;
    Double_t tmp_cnt=0;
    for( Int_t m=0; m<xbin; m++ ){
      Double_t tmp_theta_pE;
      Double_t tmp_theta_mE;
      Double_t tmp_theta_p = hist[i]->IntegralAndError( m+1,m+1,ybin/2+1, ybin,   tmp_theta_pE );
      Double_t tmp_theta_m = hist[i]->IntegralAndError( m+1,m+1,       1, ybin/2, tmp_theta_mE );
      std::cout << m << " : "
		<< " ( " << tmp_theta_p+tmp_theta_m << " +- " << sqrt(tmp_theta_pE*tmp_theta_pE+tmp_theta_mE*tmp_theta_mE) << " ) = " 
		<< " ( " << tmp_theta_p             << " +- " << tmp_theta_pE                                              << " ) +- "
		<< " ( " << tmp_theta_m             << " +- " << tmp_theta_mE                                              << " ) ";
      tmp_cnt += tmp_theta_p;
      tmp_cnt += tmp_theta_m;
      if( tmp_theta_p+tmp_theta_m==0 ){ std::cout << std::endl; continue;}
      Double_t afb     = (tmp_theta_p - tmp_theta_m ) / (tmp_theta_p + tmp_theta_m );
      Double_t afbE    = 2/(tmp_theta_p+tmp_theta_m)/(tmp_theta_p+tmp_theta_m)*sqrt(tmp_theta_p*tmp_theta_p*tmp_theta_mE*tmp_theta_mE+tmp_theta_m*tmp_theta_m*tmp_theta_pE*tmp_theta_pE);
      hist_afb[i]->SetBinContent( m+1, afb  );
      hist_afb[i]->SetBinError  ( m+1, afbE );
      std::cout << "AFB = " << afb << " +- " << afbE << std::endl;
    }
    std::cout << "total : " << tmp_cnt << " events" << std::endl;

    // for afb(gen*eff)
    std::cout << "[Entry of AFB(gen*eff)]" << std::endl;
    tmp_cnt=0;
    for( Int_t m=0; m<xbin; m++ ){
      Double_t tmp_theta_pE;
      Double_t tmp_theta_mE;
      Double_t tmp_theta_p = hist_eff[i]->IntegralAndError( m+1,m+1,ybin/2+1, ybin,   tmp_theta_pE );
      Double_t tmp_theta_m = hist_eff[i]->IntegralAndError( m+1,m+1,       1, ybin/2, tmp_theta_mE );
      std::cout << m << " : "
		<< " ( " << tmp_theta_p+tmp_theta_m << " +- " << sqrt(tmp_theta_pE*tmp_theta_pE+tmp_theta_mE*tmp_theta_mE) << " ) = " 
		<< " ( " << tmp_theta_p             << " +- " << tmp_theta_pE                                              << " ) +- "
		<< " ( " << tmp_theta_m             << " +- " << tmp_theta_mE                                              << " ) ";
      tmp_cnt += tmp_theta_p;
      tmp_cnt += tmp_theta_m;
      if( tmp_theta_p+tmp_theta_m==0 ){ std::cout << std::endl; continue;}
      Double_t afb     = (tmp_theta_p - tmp_theta_m ) / (tmp_theta_p + tmp_theta_m );
      Double_t afbE    = 2/(tmp_theta_p+tmp_theta_m)/(tmp_theta_p+tmp_theta_m)*sqrt(tmp_theta_p*tmp_theta_p*tmp_theta_mE*tmp_theta_mE+tmp_theta_m*tmp_theta_m*tmp_theta_pE*tmp_theta_pE);
      hist_afb_eff[i]->SetBinContent( m+1, afb  );
      hist_afb_eff[i]->SetBinError  ( m+1, afbE );
      std::cout << "AFB = " << afb << " +- " << afbE << std::endl;
    }
    std::cout << "total : " << tmp_cnt << " events" << std::endl;

    // for afb(gen) in afb binning
    std::cout << "[Entry of AFB(gen) in AFB bin]" << std::endl;
    tmp_cnt=0;
    for( Int_t m=0; m<xbin_afb; m++ ){
      Double_t tmp_theta_pE;
      Double_t tmp_theta_mE;
      Double_t tmp_theta_p = hist[i]->IntegralAndError( xbins_convert[fl_mode_ll][m], xbins_convert[fl_mode_ll][m+1]-1, ybin/2+1, ybin,   tmp_theta_pE );
      Double_t tmp_theta_m = hist[i]->IntegralAndError( xbins_convert[fl_mode_ll][m], xbins_convert[fl_mode_ll][m+1]-1,        1, ybin/2, tmp_theta_mE );
      std::cout << m << " : "
		<< " ( " << tmp_theta_p+tmp_theta_m << " +- " << sqrt(tmp_theta_pE*tmp_theta_pE+tmp_theta_mE*tmp_theta_mE) << " ) = " 
		<< " ( " << tmp_theta_p             << " +- " << tmp_theta_pE                                              << " ) +- "
		<< " ( " << tmp_theta_m             << " +- " << tmp_theta_mE                                              << " ) ";
      tmp_cnt += tmp_theta_p;
      tmp_cnt += tmp_theta_m;
      if( tmp_theta_p+tmp_theta_m==0 ){ std::cout << std::endl; continue;}
      Double_t afb     = (tmp_theta_p - tmp_theta_m ) / (tmp_theta_p + tmp_theta_m );
      Double_t afbE    = 2/(tmp_theta_p+tmp_theta_m)/(tmp_theta_p+tmp_theta_m)*sqrt(tmp_theta_p*tmp_theta_p*tmp_theta_mE*tmp_theta_mE+tmp_theta_m*tmp_theta_m*tmp_theta_pE*tmp_theta_pE);
      hist_afb_bin[i]->SetBinContent( m+1, afb  );
      hist_afb_bin[i]->SetBinError  ( m+1, afbE );
      std::cout << "AFB = " << afb << " +- " << afbE << std::endl;
    }
    std::cout << "total : " << tmp_cnt << " events" << std::endl;

    // for afb(gen*eff) in afb binning
    std::cout << "[Entry of AFB(gen*eff) in AFB bin]" << std::endl;
    tmp_cnt=0;
    for( Int_t m=0; m<xbin_afb; m++ ){
      Double_t tmp_theta_pE;
      Double_t tmp_theta_mE;
      Double_t tmp_theta_p = hist_eff[i]->IntegralAndError( xbins_convert[fl_mode_ll][m], xbins_convert[fl_mode_ll][m+1]-1, ybin/2+1, ybin,   tmp_theta_pE );
      Double_t tmp_theta_m = hist_eff[i]->IntegralAndError( xbins_convert[fl_mode_ll][m], xbins_convert[fl_mode_ll][m+1]-1,        1, ybin/2, tmp_theta_mE );
      std::cout << m << " : "
		<< " ( " << tmp_theta_p+tmp_theta_m << " +- " << sqrt(tmp_theta_pE*tmp_theta_pE+tmp_theta_mE*tmp_theta_mE) << " ) = " 
		<< " ( " << tmp_theta_p         << " +- " << tmp_theta_pE                                  << " ) +- "
		<< " ( " << tmp_theta_m         << " +- " << tmp_theta_mE                                  << " ) ";
      tmp_cnt += tmp_theta_p;
      tmp_cnt += tmp_theta_m;
      if( tmp_theta_p+tmp_theta_m==0 ){ std::cout << std::endl; continue;}
      Double_t afb     = (tmp_theta_p - tmp_theta_m ) / (tmp_theta_p + tmp_theta_m );
      Double_t afbE    = 2/(tmp_theta_p+tmp_theta_m)/(tmp_theta_p+tmp_theta_m)*sqrt(tmp_theta_p*tmp_theta_p*tmp_theta_mE*tmp_theta_mE+tmp_theta_m*tmp_theta_m*tmp_theta_pE*tmp_theta_pE);
      hist_afb_eff_bin[i]->SetBinContent( m+1, afb  );
      hist_afb_eff_bin[i]->SetBinError  ( m+1, afbE );
      std::cout << "AFB = " << afb << " +- " << afbE << std::endl;
    }
    std::cout << "total : " << tmp_cnt << " events" << std::endl;
  }

  // +++++++ check ++++++++++++++++++++++++++++++++++
  std::cout << "Effifciency = "
	    << hist_rec_eff[fl_mode_ll]->GetBinContent(10,10) << " +- "
	    << hist_rec_eff[fl_mode_ll]->GetBinError  (10,10) << std::endl;
  std::cout << "N(gen) = "
	    << hist[0]->GetBinContent(10,10) << " +- "
	    << hist[0]->GetBinError  (10,10) << std::endl;
  std::cout << "N(gen*eff) = "
	    << hist_eff[0]->GetBinContent(10,10) << " +- "
	    << hist_eff[0]->GetBinError  (10,10) << std::endl;
  // +++++++ draw ++++++++++++++++++++++++++++++++++
  TLegend* legend = new TLegend( 0.75,0.75,0.99,0.99 );
  legend->SetHeader( Form("lep%d, %s, %s, %s", fl_mode_ll, fl_A7, fl_A9, fl_A10) );
  legend->AddEntry( hist_afb    [0], "gen",     "P" );
  legend->AddEntry( hist_afb_eff[0], "gen*eff", "P" );
  
  c1->Draw();
  c1->cd(1);
  TH2D* waku = new TH2D("AFB","A_{FB}",2,0,25,2,-0.4,0.4);
  waku->SetXTitle( xlabel   );
  waku->SetYTitle( "A_{FB}" );
  waku->GetXaxis()->CenterTitle();
  waku->GetYaxis()->CenterTitle();
  waku->Draw();
  hist_afb        [0]->Draw("PE0same");
  hist_afb_eff    [0]->Draw("PE0same");
  
  c1->cd(2);
  waku->Draw();
  hist_afb_bin    [0]->Draw("PE0same");
  hist_afb_eff_bin[0]->Draw("PE0same");
  
  for( Int_t i=0; i<Nchain; i++ ){
    histx_gen    [i]->Scale(1/histx_gen    [i]->Integral());
    histy_gen    [i]->Scale(1/histy_gen    [i]->Integral());
    histx_gen_eff[i]->Scale(1/histx_gen_eff[i]->Integral());
    histy_gen_eff[i]->Scale(1/histy_gen_eff[i]->Integral());
  }

  c1->cd(3);
  TH2D* wakugenx = Waku( Nchain, histx_gen, xlabel, Form("%s (gen)", xlabel), Form("%s (gen)", xlabel) );
  wakugenx->Draw();
  legend->Draw();
  histx_gen    [0]->Draw("Psame");
  histx_gen_eff[0]->Draw("Psame");

  c1->cd(4);
  TH2D* wakugeny = Waku( Nchain, histy_gen_eff, ylabel, Form("%s (gen)", ylabel), Form("%s (gen)", ylabel) );
  wakugeny->Draw();
  histy_gen    [0]->Draw("Psame");
  histy_gen_eff[0]->Draw("Psame");

  c1->cd(5); hist    [0]->Draw("COLZ");
  c1->cd(6); hist_eff[0]->Draw("COLZ");

  c1->Update();
  if( flag_save ) c1->Print( Form("pic/2d_gen_q2_theta_correction_table_lep%d_set%s_%s_%s_%s_mb495.eps", fl_mode_ll, setname, fl_A7, fl_A9, fl_A10) );

  // LOG for correction function
    std::cout << std::setw( 5) << std::right << "rm_l"
	      << std::setw( 5) << std::right << "A7"
	      << std::setw( 5) << std::right << "A9"
	      << std::setw( 5) << std::right << "A10"
	      << std::setw( 5) << std::right << "q2"
	      << std::setw(15) << std::right << "AFB(meas.)"
	      << std::setw(15) << std::right << "error"
	      << std::setw(15) << std::right << "AFB(gen)"
	      << std::setw(15) << std::right << "error"

	      << std::setw(15) << std::right << "N+(meas)"
	      << std::setw(15) << std::right << "error"
	      << std::setw(15) << std::right << "N-(meas)"
	      << std::setw(15) << std::right << "error"      

	      << std::setw(15) << std::right << "N+(gen)"
	      << std::setw(15) << std::right << "error"
	      << std::setw(15) << std::right << "N-(gen)"
	      << std::setw(15) << std::right << "error"      

	      << std::endl;
  for( Int_t i=0; i<xbin_afb; i++ ){
    
    Double_t theta_pE_gen;
    Double_t theta_mE_gen;
    Double_t theta_p_gen = hist[0]->IntegralAndError( xbins_convert[fl_mode_ll][i], xbins_convert[fl_mode_ll][i+1]-1, ybin/2+1, ybin,   theta_pE_gen );
    Double_t theta_m_gen = hist[0]->IntegralAndError( xbins_convert[fl_mode_ll][i], xbins_convert[fl_mode_ll][i+1]-1,        1, ybin/2, theta_mE_gen );

    Double_t theta_pE_meas;
    Double_t theta_mE_meas;
    Double_t theta_p_meas = hist_eff[0]->IntegralAndError( xbins_convert[fl_mode_ll][i], xbins_convert[fl_mode_ll][i+1]-1, ybin/2+1, ybin,   theta_pE_meas );
    Double_t theta_m_meas = hist_eff[0]->IntegralAndError( xbins_convert[fl_mode_ll][i], xbins_convert[fl_mode_ll][i+1]-1,        1, ybin/2, theta_mE_meas );

    Int_t tmp_A7 = 0;
    if     ( !strcmp(fl_A7,"norm") ) tmp_A7 =  100;
    else if( !strcmp(fl_A7,"flip") ) tmp_A7 = -100;
    else                       std::cerr << "[ABORT] Wrong fl_A7 : " << fl_A7 << std::endl, abort();
    std::cout << std::setw( 5) << std::right << fl_mode_ll
	      << std::setw( 5) << std::right << tmp_A7
	      << std::setw( 5) << std::right << atoi(fl_A9)
	      << std::setw( 5) << std::right << atoi(fl_A10)
	      << std::setw( 5) << std::right << i
	      << std::setw(15) << std::right << hist_afb_eff_bin[0]->GetBinContent(i+1)
	      << std::setw(15) << std::right << hist_afb_eff_bin[0]->GetBinError  (i+1)
	      << std::setw(15) << std::right << hist_afb_bin    [0]->GetBinContent(i+1)
	      << std::setw(15) << std::right << hist_afb_bin    [0]->GetBinError  (i+1)
	      << std::setw(15) << std::right << theta_p_meas  / scale_event_sig
	      << std::setw(15) << std::right << theta_pE_meas / scale_event_sig
	      << std::setw(15) << std::right << theta_m_meas  / scale_event_sig
	      << std::setw(15) << std::right << theta_mE_meas / scale_event_sig
	      << std::setw(15) << std::right << theta_p_gen   / scale_event_sig
	      << std::setw(15) << std::right << theta_pE_gen  / scale_event_sig
	      << std::setw(15) << std::right << theta_m_gen   / scale_event_sig
	      << std::setw(15) << std::right << theta_mE_gen  / scale_event_sig
	      << std::setw( 5) << std::right << " HOGE fl_corr" << (Int_t)fl_corr
	      << std::endl;
  }
  
  std::cout << "finish" << std::flush;
  if( fl_appRun ) app.Run();

  delete[] chain;
  delete[] tmphist;
  delete[] hist;
  delete   c1;
    
  return 0;
}
