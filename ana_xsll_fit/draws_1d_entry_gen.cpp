#include <iostream>
#include <math.h>

#include "../Util/Style.h"
#include "../Util/Canvas.h"
#include "../Util/const.h"
#include "../Util/MChain.h"
#include "../Util/MCut.h"
#include "../Util/MCut_array.h"
#include "../Util/FitFunc.h"
#include "../Set/Nominal_cut_selection.h"
#include "../Set/Branch.h"
#include "../Set/makeCut.h"

#include "draws_.h"

#include <vector>
#include <stdlib.h>
#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TString.h>
#include <TF1.h>
#include <TArrow.h>

Int_t main( Int_t argc, Char_t** argv ){
  TApplication app( "app", &argc, argv );
  Style();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( !(argc==5 || argc==6) ) std::cerr << "wrong input" << std::endl
    					<< " Usage : ./draws_1d_Mbc (int)fl_mc (int)fl_exp (int)fl_linei (int)fl_nline [(int)fl_appRun]" << std::endl, abort();
  Int_t   fl_mc      = atoi(argv[1]); // 0(gmc), 1-4(ccmc)
  Int_t   fl_exp     = atoi(argv[2]); // remove initial 0 (07 -> 7)
  Char_t* fl_exp2    = argv[2];
  Int_t   fl_nline   = atoi(argv[3]);
  Int_t   fl_linei   = atoi(argv[4]);
  Int_t   fl_appRun  = 1;
  if( argc==6 ) fl_appRun = atoi( argv[5] );
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Char_t* tname         = "h12";
  const Int_t   fl_mode_ll    = 0; // no mean
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::stringstream sTmp;
  Char_t* infile = new Char_t[1024];
  if     ( fl_mc==0 ) sTmp << "~/ewp/ana/data/gmc_beforeskim/hbk5/right/hbk_org_gen/total/9999/gMC"; // gmc before skim
  else if( fl_mc==1 ) sTmp << "~/ewp/ana/data/gmc_cc/hbk5/right/hbk_org_gen/CC_mixedjpsi_f";           // cc-cm
  else if( fl_mc==2 ) sTmp << "~/ewp/ana/data/gmc_cc/hbk5/right/hbk_org_gen/CC_chargedjpsi_f";         // cc-mc
  else if( fl_mc==3 ) sTmp << "~/ewp/ana/data/gmc_cc/hbk5/right/hbk_org_gen/CC_mixedpsi2s";          // cc-mc
  else if( fl_mc==4 ) sTmp << "~/ewp/ana/data/gmc_cc/hbk5/right/hbk_org_gen/CC_chargedpsi2s";        // cc-mc
  else                std::cerr << "[ABORT] Wrong fl_mc : " << fl_mc << std::endl, abort();
    strcpy( infile, (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const Int_t Ncut = 6;
  Char_t** add_cut = new Char_t*[Ncut];
  for( Int_t i=0; i<Ncut; i++ ){
    add_cut[i] = new Char_t[4096];
    sTmp << Form("exprun>=%d0000 && exprun<%d0000", fl_exp, fl_exp+1);
    if     ( i==0 ) sTmp << " && ( (gb1_cc==443 && gm_b1==   1 && gm_bg1<3) || (gb2_cc==443 && gm_b2==   1 && gm_bg2<3) )";
    else if( i==1 ) sTmp << " && ( (gb1_cc==443 && gm_b1==  10 && gm_bg1<3) || (gb2_cc==443 && gm_b2==  10 && gm_bg2<3) )";
    else if( i==2 ) sTmp << " && ( (gb1_cc==443 && gm_b1== 101 && gm_bg1<3) || (gb2_cc==443 && gm_b2== 101 && gm_bg2<3) )";
    else if( i==3 ) sTmp << " && ( (gb1_cc==443 && gm_b1== 110 && gm_bg1<3) || (gb2_cc==443 && gm_b2== 110 && gm_bg2<3) )";
    else if( i==4 ) sTmp << " && ( (gb1_cc==443 && gm_b1==1001 && gm_bg1<3) || (gb2_cc==443 && gm_b2==1001 && gm_bg2<3) )";
    else if( i==5 ) sTmp << " && ( (gb1_cc==443 && gm_b1==1010 && gm_bg1<3) || (gb2_cc==443 && gm_b2==1010 && gm_bg2<3) )";
    strcpy( add_cut[i], (char*)sTmp.str().c_str() );
    sTmp.str("");
    sTmp.clear();
  }



  // +++++++ make chain-tree ++++++++++++++++++++++++++++++++++
  std::cout << Form(" ************************ make chain-tree ( %s ) *************************************",tname) << std::endl;
  std::cout << "<infile> ";
  MChain* chain = new MChain( infile, tname, branch_table(), fl_nline, "_*.root", fl_linei );
  nominal_cut_selection( chain, fl_mode_ll )( chain->GetCut(), tname );

  // ++++++++++++++++++++++++
  // cut change
  chain->GetCut()->Set( "gm_l",          0 );
  chain->GetCut()->Set( "llg_m",         0 );
  chain->GetCut()->Set( "gm_fl_xs*Xs_m", 0 );
  //
  // ++++++++++++++++++++++++
  std::cout << std::endl
	    << " ************************ changed cut *************************************" << std::endl;
  std::cout << "<chain-tree> " << std::endl;
  chain->GetCut()->Display(0);
  chain->MakeTree();

  std::cout << std::endl
	    << " ************************ add cut *************************************" << std::endl;
  for( Int_t i=0; i<Ncut; i++ ){
    std::cout << "add_cut" << std::setw(2) << std::right << i << " : " << add_cut[i] << std::endl;
  }

  // LOG
  std::cout << std::setw(3) << std::right << "mc"
	    << std::setw(3) << std::right << "mode"
	    << std::setw(3) << std::right << "exp"
	    << std::setw(8) << std::right << "nevt"
    
	    << std::endl;
  for( Int_t i=0; i<Ncut; i++ ){
    std::cout << std::setw(3) << std::right <<  fl_mc
	      << std::setw(3) << std::right <<  i
	      << std::setw(3) << std::right <<  fl_exp
	      << std::setw(8) << std::right << chain->GetTree()->GetEntries(add_cut[i])
	      << " HOGE"
	      << std::endl;
  }
  
  
  std::cout << "finish" << std::endl;
  if( fl_appRun ) app.Run();

  return 0;
}
